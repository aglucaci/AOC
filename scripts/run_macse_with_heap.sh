#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF' >&2
Usage: run_macse_with_heap.sh --xmx-mb N [--xms-mb N] [--macse-bin PATH] -- <macse args>
EOF
  exit 2
}

MACSE_BIN="macse"
XMS_MB="512"
XMX_MB=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --macse-bin)
      MACSE_BIN="$2"; shift 2;;
    --xms-mb)
      XMS_MB="$2"; shift 2;;
    --xmx-mb)
      XMX_MB="$2"; shift 2;;
    --)
      shift
      break;;
    *)
      usage;;
  esac
done

[[ -n "${XMX_MB}" ]] || usage
[[ $# -gt 0 ]] || usage

macse_path="$(command -v "${MACSE_BIN}")"
tmp_macse="$(mktemp)"
cleanup() {
  rm -f "${tmp_macse}"
}
trap cleanup EXIT

if ! grep -q '^default_jvm_mem_opts=' "${macse_path}"; then
  echo "[ERROR] Could not find default_jvm_mem_opts in ${macse_path}" >&2
  echo "[ERROR] Update scripts/run_macse_with_heap.sh for this MACSE launcher format." >&2
  exit 1
fi

sed "s/^default_jvm_mem_opts=.*/default_jvm_mem_opts=\"-Xms${XMS_MB}m -Xmx${XMX_MB}m\"/" "${macse_path}" > "${tmp_macse}"
chmod +x "${tmp_macse}"
"${tmp_macse}" "$@"
