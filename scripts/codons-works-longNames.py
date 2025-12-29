#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Feb  2 13:04:15 2020
@author: Alexander G. Lucaci

Given a protein FASTA and a transcript/CDS FASTA:
  - For each protein, scan transcripts and slide a window until translation matches protein
  - Output the matched coding DNA (codons)

Additions (keeping Snakemake compatibility):
  1) FASTA header/organism parsing:
       supports tags like [organism=Myotis davidii]
       also supports legacy headers with [Myotis davidii]
  2) If multiple transcripts for the same organism exist, keep ONLY transcript=X1 (first encountered)
  3) Each transcript may match ONLY ONE protein (one-to-one mapping)
"""

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

# =============================================================================
# Declares (Snakemake params)  -- keep compatibility
# =============================================================================
PROTEIN_FASTA     = snakemake.params.Prot
TRANSCRIPTS_FASTA = snakemake.params.Nuc
OUTPUT            = snakemake.params.Out

results  = []
no_match = []

logfile = OUTPUT + ".log"
with open(logfile, "w") as fh:
    fh.write("")

# =============================================================================
# Options
# =============================================================================
# Header parsing behavior:
#   "auto"      : try [organism=...]/[species=...]/[taxname=...]/[taxon=...], then legacy [Genus species], then substring fallback
#   "tag"       : only key/value tags
#   "legacy"    : only legacy [Genus species]
#   "substring" : substring match only (least strict)
HEADER_MODE = "auto"

# If there are multiple transcripts for the same organism in the transcript FASTA,
# keep only the first one (X1) for matching.
KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM = True

# Each transcript may match only ONE protein (global one-to-one constraint)
ENFORCE_TRANSCRIPT_UNIQUENESS = True

# =============================================================================
# Globals for constraints
# =============================================================================
# transcript_id -> already claimed by a protein?
USED_TRANSCRIPT_IDS = set()

# organism_key -> first transcript_id encountered for that organism (X1)
FIRST_TRANSCRIPT_ID_BY_ORG = {}

# =============================================================================
# Logging helper
# =============================================================================
def log(msg: str) -> None:
    with open(logfile, "a") as fh:
        print(msg, file=fh)

# =============================================================================
# Header parsing helpers
# =============================================================================
BRACKET_RE = re.compile(r"\[([^\[\]]+)\]")

def parse_bracket_keyvals(desc: str) -> dict:
    """
    Parse bracket tags like:
      "... [organism=Myotis davidii] [GeneID=102755104]"
    Returns dict {key.lower(): value}.
    """
    tags = {}
    for chunk in BRACKET_RE.findall(desc or ""):
        chunk = chunk.strip()
        if "=" in chunk:
            k, v = chunk.split("=", 1)
            tags[k.strip().lower()] = v.strip()
    return tags

def extract_organism(desc: str, mode: str = "auto") -> str:
    """
    Extract organism from a FASTA description using HEADER_MODE.
    Supports:
      - [organism=Genus species] (preferred)
      - [species=...], [taxname=...], [taxon=...]
      - legacy: [Genus species]
    """
    desc = desc or ""
    tags = parse_bracket_keyvals(desc)

    if mode in ("auto", "tag"):
        for k in ("organism", "species", "taxname", "taxon"):
            if k in tags and tags[k]:
                return tags[k]
        if mode == "tag":
            return ""

    if mode in ("auto", "legacy"):
        for chunk in BRACKET_RE.findall(desc):
            chunk = chunk.strip()
            if "=" in chunk:
                continue
            parts = chunk.split()
            if len(parts) >= 2 and parts[0][:1].isupper():
                return " ".join(parts[:2])
        if mode == "legacy":
            return ""

    return ""

def organism_matches(protein_org: str, transcript_desc: str, mode: str = "auto") -> bool:
    """
    Decide if a transcript header belongs to the same organism as the protein.
    - If protein_org is missing, don't filter.
    - If transcript has an organism tag, require exact equality.
    - Otherwise optionally fall back to substring matching in auto mode.
    """
    if not protein_org:
        return True

    p = protein_org.strip().lower()
    t_desc = transcript_desc or ""

    if mode in ("auto", "tag", "legacy"):
        t_org = extract_organism(t_desc, mode=mode).strip().lower()
        if t_org:
            return t_org == p
        if mode == "tag":
            # if transcript has no tag, don't exclude it (keeps backwards compatibility)
            return True
        # auto/legacy fallback: substring
        return p in t_desc.lower()

    # substring mode
    return p in t_desc.lower()

# =============================================================================
# Original duplicate guard (kept for compatibility)
# =============================================================================
def already_in_results(transcript_desc: str) -> bool:
    global results
    for record in results:
        if transcript_desc == record.description:
            return True
    return False

# =============================================================================
# Core matching logic
# =============================================================================
def Process(protein_id: str, protein_desc: str, protein_seq, transcripts_fasta: str, protein_org: str):
    """
    For a given protein:
      - consider transcripts that match organism
      - enforce transcript=X1 per organism if enabled
      - enforce each transcript matches only one protein if enabled
      - slide window until translation equals protein
    Returns SeqRecord or "NO_MATCH"
    """
    global USED_TRANSCRIPT_IDS, FIRST_TRANSCRIPT_ID_BY_ORG

    NT_SEQ_LENGTH = len(protein_seq) * 3

    with open(transcripts_fasta, "r") as transcript_handle:
        for transcript_record in SeqIO.parse(transcript_handle, "fasta"):
            transcript_id   = transcript_record.id
            transcript_desc = transcript_record.description
            transcript_seq  = transcript_record.seq

            # organism filter
            if not organism_matches(protein_org, transcript_desc, mode=HEADER_MODE):
                continue

            # Determine organism key for X1 bookkeeping
            tx_org = extract_organism(transcript_desc, mode=HEADER_MODE)
            if not tx_org and protein_org:
                tx_org = protein_org
            org_key = (tx_org or "").strip().lower()

            # Enforce transcript=X1 (first transcript encountered for that organism)
            if KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM and org_key:
                if org_key not in FIRST_TRANSCRIPT_ID_BY_ORG:
                    FIRST_TRANSCRIPT_ID_BY_ORG[org_key] = transcript_id  # claim X1
                # If this is not X1, skip
                if transcript_id != FIRST_TRANSCRIPT_ID_BY_ORG[org_key]:
                    continue

            # Enforce transcript uniqueness (one transcript -> one protein)
            if ENFORCE_TRANSCRIPT_UNIQUENESS and transcript_id in USED_TRANSCRIPT_IDS:
                # If we're forced to only use X1 and it's already used, there's no other transcript allowed.
                # That means this protein cannot be matched.
                if KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM and org_key:
                    return "NO_MATCH"
                continue

            # Slide window across transcript
            start = 0
            DONE = False
            coding_seq = None

            while start + NT_SEQ_LENGTH <= len(transcript_seq):
                try:
                    candidate = transcript_seq[start:start + NT_SEQ_LENGTH]
                    translated = candidate.translate()
                except Exception:
                    translated = None

                exists = already_in_results(transcript_desc)

                if translated is not None and str(translated) == str(protein_seq) and not exists:
                    DONE = True
                    coding_seq = candidate

                    # Claim this transcript for this protein (one-to-one)
                    if ENFORCE_TRANSCRIPT_UNIQUENESS:
                        USED_TRANSCRIPT_IDS.add(transcript_id)

                    break
                else:
                    start += 1

            if DONE:
                # Preserve provenance in header
                out_id = f"{protein_id}|{transcript_id}"
                out_desc = f"organism={tx_org} | protein={protein_desc} | transcript={transcript_desc}"
                return SeqRecord(seq=coding_seq, id=out_id, description=out_desc)

            # If we're enforcing X1, and we just attempted X1 for this organism but no match,
            # we must NOT consider X2/X3 for this organism; therefore we can stop early.
            if KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM and org_key and transcript_id == FIRST_TRANSCRIPT_ID_BY_ORG.get(org_key):
                return "NO_MATCH"

    return "NO_MATCH"

# =============================================================================
# Main sanity counts
# =============================================================================
def main(PROTEIN, TRANSCRIPTS):
    print("# TRANSCRIPT INPUT FILE:", TRANSCRIPTS)
    print("# PROTEIN INPUT FILE:", PROTEIN)

    trans_count = 0
    with open(TRANSCRIPTS, "r") as handle:
        for _ in SeqIO.parse(handle, "fasta"):
            trans_count += 1
    print("# Transcripts:", trans_count)

    prot_count = 0
    with open(PROTEIN, "r") as handle:
        for _ in SeqIO.parse(handle, "fasta"):
            prot_count += 1
    print("# Proteins:", prot_count)

    return trans_count, prot_count

# =============================================================================
# Main
# =============================================================================
print("# =============================================================================")
print("# Processing... ")
print("# =============================================================================")

trans_count, prot_count = main(PROTEIN_FASTA, TRANSCRIPTS_FASTA)

# Create empty output file
with open(OUTPUT, "w") as fh:
    fh.write("")

log(f"# HEADER_MODE={HEADER_MODE}")
log(f"# KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM={KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM}")
log(f"# ENFORCE_TRANSCRIPT_UNIQUENESS={ENFORCE_TRANSCRIPT_UNIQUENESS}")
log(f"# PROTEIN_FASTA={PROTEIN_FASTA}")
log(f"# TRANSCRIPTS_FASTA={TRANSCRIPTS_FASTA}")

# Iterate proteins
with open(PROTEIN_FASTA, "r") as prot_handle:
    for record in SeqIO.parse(prot_handle, "fasta"):
        protein_id   = record.id
        protein_desc = record.description
        protein_seq  = record.seq

        log("\n")
        log(f"# Processing: {protein_desc}")

        # keep original quality skip behavior, but make partial check case-insensitive
        if "LOW QUALITY PROTEIN" in str(protein_desc) or "partial" in str(protein_desc).lower():
            log(f"# Skipping due to quality issues: {protein_desc}")
            continue

        protein_org = extract_organism(protein_desc, mode=HEADER_MODE)
        if HEADER_MODE == "substring" and not protein_org:
            # best-effort: try auto extraction so substring mode still works
            protein_org = extract_organism(protein_desc, mode="auto")

        log(f"# Organism: {protein_org}")

        tx_record = Process(
            protein_id=protein_id,
            protein_desc=protein_desc,
            protein_seq=protein_seq,
            transcripts_fasta=TRANSCRIPTS_FASTA,
            protein_org=protein_org
        )

        if isinstance(tx_record, SeqRecord):
            log(f"# Match: {tx_record.id}")
            results.append(tx_record)
        else:
            log("# -- NO Match --")
            no_match.append(protein_desc)

# Write output
log(f"# Writing data to: {OUTPUT}")
SeqIO.write(results, OUTPUT, "fasta")

# Report no matches
log("--- The following had no matches")
for item in no_match:
    log(item)

print(f"# Done. Wrote {len(results)} matched CDS records to {OUTPUT}")
print(f"# No matches: {len(no_match)} (see {logfile})")
