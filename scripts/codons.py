#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
codons.py
Match protein sequences to CDS by sliding translation windows.

Rules:
  - Organism-aware matching via FASTA headers
  - Only first transcript per organism (X1)
  - Each transcript matches at most one protein
  - FASTA header line must be <= 30 chars (ID only; no description in FASTA)
  - Full provenance written to a .map.tsv file
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import hashlib

# =============================================================================
# Snakemake params (unchanged)
# =============================================================================
PROT_FA = snakemake.params.Prot
CDS_FA  = snakemake.params.Nuc
OUT_FA  = snakemake.params.Out

LOG = OUT_FA + ".log"
MAP = OUT_FA + ".map.tsv"

# =============================================================================
# Options
# =============================================================================
HDR_MODE   = "auto"     # auto | tag | legacy | substring
USE_X1     = True      # only first transcript per organism (X1)
UNIQUE_TX  = True      # transcript maps to only one protein
ID_MAXLEN  = 30        # FASTA header line length limit (ID only)

# =============================================================================
# Global state
# =============================================================================
USED_TX   = set()      # transcript IDs already used
X1_BY_ORG = {}         # organism -> first transcript ID

results  = []
no_match = []

# =============================================================================
# Logging
# =============================================================================
def log(msg):
    with open(LOG, "a") as fh:
        print(msg, file=fh)

open(LOG, "w").close()

# =============================================================================
# FASTA ID shortening (<=30 chars)
# =============================================================================
def short_id(raw: str, max_len: int = 30) -> str:
    """
    Return an ID <= max_len.
    If truncation is needed, append a stable 6-char hash.
    """
    if len(raw) <= max_len:
        return raw
    h = hashlib.md5(raw.encode()).hexdigest()[:6]
    keep = max_len - 7  # "_" + 6 chars
    return f"{raw[:keep]}_{h}"

# =============================================================================
# Header parsing
# =============================================================================
BRACKET_RE = re.compile(r"\[([^\[\]]+)\]")

def parse_tags(desc):
    tags = {}
    for chunk in BRACKET_RE.findall(desc or ""):
        if "=" in chunk:
            k, v = chunk.split("=", 1)
            tags[k.lower().strip()] = v.strip()
    return tags

def get_org(desc, mode="auto"):
    desc = desc or ""
    tags = parse_tags(desc)

    if mode in ("auto", "tag"):
        for k in ("organism", "species", "taxname", "taxon"):
            if k in tags and tags[k]:
                return tags[k]
        if mode == "tag":
            return ""

    if mode in ("auto", "legacy"):
        for chunk in BRACKET_RE.findall(desc):
            if "=" in chunk:
                continue
            parts = chunk.split()
            if len(parts) >= 2 and parts[0][0].isupper():
                return " ".join(parts[:2])

    return ""

def org_match(prot_org, tx_desc, mode="auto"):
    if not prot_org:
        return True

    p = prot_org.strip().lower()
    t = get_org(tx_desc, mode).strip().lower()

    if t:
        return t == p
    if mode == "tag":
        return True
    return p in (tx_desc or "").lower()

# =============================================================================
# Core matcher
# =============================================================================
def match_protein(pid, pdesc, pseq, cds_fa, prot_org):
    """
    Returns:
      (SeqRecord, mapping_dict) or ("NO_MATCH", None)
    """
    nt_len = len(pseq) * 3

    with open(cds_fa) as fh:
        for tx in SeqIO.parse(fh, "fasta"):
            tx_id   = tx.id
            tx_desc = tx.description
            tx_seq  = tx.seq

            if not org_match(prot_org, tx_desc, HDR_MODE):
                continue

            tx_org = get_org(tx_desc, HDR_MODE) or prot_org
            org_k  = tx_org.lower() if tx_org else ""

            # X1 rule: only first transcript per organism
            if USE_X1 and org_k:
                if org_k not in X1_BY_ORG:
                    X1_BY_ORG[org_k] = tx_id
                if tx_id != X1_BY_ORG[org_k]:
                    continue

            # uniqueness rule: transcript can only be used once
            if UNIQUE_TX and tx_id in USED_TX:
                # if X1-only, there's no fallback transcript allowed
                if USE_X1:
                    return "NO_MATCH", None
                continue

            # slide CDS window
            for i in range(0, len(tx_seq) - nt_len + 1):
                try:
                    cds = tx_seq[i:i + nt_len]
                    if str(cds.translate()) == str(pseq):
                        if UNIQUE_TX:
                            USED_TX.add(tx_id)

                        # ID used in FASTA header (and ONLY this; no description)
                        rid = short_id(f"{pid}|{tx_id}", ID_MAXLEN)

                        # IMPORTANT: empty description so FASTA header line is ONLY rid
                        rec = SeqRecord(seq=cds, id=rid, description="")

                        mapping = {
                            "id30": rid,
                            "protein_id": pid,
                            "protein_desc": pdesc,
                            "transcript_id": tx_id,
                            "transcript_desc": tx_desc,
                            "organism": tx_org,
                            "cds_start": str(i),
                            "cds_end": str(i + nt_len),
                        }
                        return rec, mapping
                except Exception:
                    pass

            # If X1-only and we just tried X1 with no match, stop early for this protein
            if USE_X1 and org_k:
                return "NO_MATCH", None

    return "NO_MATCH", None

# =============================================================================
# Main
# =============================================================================
log(f"# HDR_MODE={HDR_MODE}")
log(f"# USE_X1={USE_X1}")
log(f"# UNIQUE_TX={UNIQUE_TX}")
log(f"# ID_MAXLEN={ID_MAXLEN}")
log(f"# PROT_FA={PROT_FA}")
log(f"# CDS_FA={CDS_FA}")
log(f"# OUT_FA={OUT_FA}")

# write mapping header
with open(MAP, "w") as mf:
    mf.write("\t".join([
        "id30", "organism",
        "protein_id", "transcript_id",
        "cds_start", "cds_end",
        "protein_desc", "transcript_desc"
    ]) + "\n")

with open(PROT_FA) as fh:
    for rec in SeqIO.parse(fh, "fasta"):
        pid   = rec.id
        pdesc = rec.description
        pseq  = rec.seq

        if "LOW QUALITY PROTEIN" in pdesc or "partial" in pdesc.lower():
            log(f"# skip {pid}")
            continue

        org = get_org(pdesc, HDR_MODE)
        log(f"\n# {pid} org={org}")

        out_rec, mapping = match_protein(pid, pdesc, pseq, CDS_FA, org)
        if isinstance(out_rec, SeqRecord):
            log(f"# match {out_rec.id}")
            results.append(out_rec)

            # write mapping line
            with open(MAP, "a") as mf:
                mf.write("\t".join([
                    mapping["id30"],
                    mapping["organism"],
                    mapping["protein_id"],
                    mapping["transcript_id"],
                    mapping["cds_start"],
                    mapping["cds_end"],
                    mapping["protein_desc"],
                    mapping["transcript_desc"],
                ]) + "\n")
        else:
            log("# no match")
            no_match.append(pdesc)

SeqIO.write(results, OUT_FA, "fasta")

log("\n# NO MATCH")
for x in no_match:
    log(x)

print(f"Done. {len(results)} matches written to {OUT_FA}")
print(f"Mapping written to {MAP}")
