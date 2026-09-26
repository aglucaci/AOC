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

error() {
  echo "[ERROR] $*" >&2
}

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

if ! macse_path="$(command -v "${MACSE_BIN}")"; then
  error "Could not resolve MACSE launcher '${MACSE_BIN}' on PATH."
  error "Activate the conda environment that provides 'macse' or set 'macse_launcher' to an executable path."
  exit 1
fi

if [[ ! -f "${macse_path}" ]]; then
  error "Resolved MACSE launcher is not a regular file: ${macse_path}"
  exit 1
fi

if [[ ! -r "${macse_path}" ]]; then
  error "Resolved MACSE launcher is not readable: ${macse_path}"
  exit 1
fi

macse_dir="$(cd "$(dirname "${macse_path}")" && pwd -P)"
macse_root="${macse_dir}"
if [[ "$(basename "${macse_dir}")" == "bin" ]]; then
  macse_root="$(cd "${macse_dir}/.." && pwd -P)"
fi

if ! tmp_dir="$(mktemp -d)"; then
  error "Failed to create temporary staging directory."
  exit 1
fi

launcher_rel="${macse_path#${macse_root}/}"
if [[ "${launcher_rel}" == "${macse_path}" ]]; then
  launcher_rel="$(basename "${macse_path}")"
fi

tmp_macse="${tmp_dir}/${launcher_rel}"
cleanup() {
  rm -rf "${tmp_dir}"
}
trap cleanup EXIT

if ! grep -q '^default_jvm_mem_opts=' "${macse_path}"; then
  error "Could not find default_jvm_mem_opts in ${macse_path}"
  error "Update scripts/run_macse_with_heap.sh for this MACSE launcher format."
  exit 1
fi

mkdir -p "$(dirname "${tmp_macse}")"
sed "s/^default_jvm_mem_opts=.*/default_jvm_mem_opts=\"-Xms${XMS_MB}m -Xmx${XMX_MB}m\"/" "${macse_path}" > "${tmp_macse}"
chmod +x "${tmp_macse}"

# MACSE launchers often resolve their resources relative to the script path ($0),
# often from within a conda env rooted above bin/. Preserve that layout in the
# temporary staging tree.
jar_count="$(find "${macse_root}" -maxdepth 6 -type f -iname '*macse*.jar' | wc -l | tr -d '[:space:]')"
if [[ "${jar_count}" == "0" ]]; then
  error "No MACSE jar files were found under: ${macse_root}"
  error "Expected to find files like macse_v*.jar somewhere under the active conda environment."
  exit 1
fi

while IFS= read -r jar_path; do
  jar_rel="${jar_path#${macse_root}/}"
  staged_jar="${tmp_dir}/${jar_rel}"
  mkdir -p "$(dirname "${staged_jar}")"
  cp "${jar_path}" "${staged_jar}"

  # Some MACSE launchers expect the jar to sit directly beside bin/macse.
  cp "${jar_path}" "$(dirname "${tmp_macse}")/$(basename "${jar_path}")"
done < <(find "${macse_root}" -maxdepth 6 -type f -iname '*macse*.jar')

staged_jar_count="$(find "${tmp_dir}" -type f -iname '*macse*.jar' | wc -l | tr -d '[:space:]')"
if [[ "${staged_jar_count}" == "0" ]]; then
  error "Failed to stage MACSE jar files into temporary directory: ${tmp_dir}"
  exit 1
fi

echo "[INFO] MACSE launcher: ${macse_path}" >&2
echo "[INFO] MACSE environment root: ${macse_root}" >&2
echo "[INFO] Staged launcher path: ${tmp_macse}" >&2
echo "[INFO] Staged ${staged_jar_count} jar file(s) with heap opts -Xms${XMS_MB}m -Xmx${XMX_MB}m" >&2

"${tmp_macse}" "$@"
