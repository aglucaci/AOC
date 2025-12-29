#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:04:15 2020
@author: Alexander G. Lucaci

Original idea:
    Given a protein sequence and a transcript sequence
    find the codons by stepping over the transcript sequence until the translated sequence matches the protein sequence.
    Output the coding DNA segment (codons) from the transcript.

UPDATE (2025-12-27):
    - Keep Snakemake compatibility (snakemake.params.Prot / Nuc / Out).
    - Add an option to deal with FASTA headers in the CDS and protein file:
        * parse organism from headers like: [organism=Myotis davidii]
        * also supports headers like: [Myotis davidii] (legacy)
        * and can fall back to substring matching in the header text
    - If there is more than one transcript of the same organism in the transcript file,
      keep only transcript=X1 (the first encountered for that organism).
"""  # :contentReference[oaicite:0]{index=0}

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

# =============================================================================
# Declares (Snakemake params)
# =============================================================================
PROTEIN_FASTA     = snakemake.params.Prot
TRANSCRIPTS_FASTA = snakemake.params.Nuc
OUTPUT            = snakemake.params.Out

logfile = OUTPUT + ".log"

results  = []
no_match = []

with open(logfile, "w") as fh2:
    print("", file=fh2)

# =============================================================================
# Options (header handling + transcript X1 behavior)
# =============================================================================
# Header parsing behavior:
#   "auto"      : try bracket tags like [organism=...], then [species=...], then legacy [Genus species]
#   "tag"       : ONLY use bracket key/value tags (organism/species/taxname/taxon). If missing, no filtering.
#   "legacy"    : ONLY use legacy bracket species like [Myotis davidii]
#   "substring" : ignore tags; just use substring check of organism name in transcript header
HEADER_MODE = "auto"

# Keep only the FIRST transcript per organism encountered in TRANSCRIPTS_FASTA (X1)
KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM = True

# =============================================================================
# Regex + header parsing helpers
# =============================================================================
BRACKET_TAG_RE = re.compile(r"\[([^\[\]]+)\]")  # content inside [ ... ]

def log(msg: str) -> None:
    with open(logfile, "a") as fh2:
        print(msg, file=fh2)

def parse_bracket_keyvals(desc: str) -> dict:
    """
    Parse bracket tags like:
      "... [organism=Myotis davidii] [GeneID=102755104]"
    Returns dict of {key.lower(): value}.
    """
    tags = {}
    for chunk in BRACKET_TAG_RE.findall(desc or ""):
        chunk = chunk.strip()
        if "=" in chunk:
            k, v = chunk.split("=", 1)
            tags[k.strip().lower()] = v.strip()
    return tags

def extract_organism_from_desc(desc: str, mode: str = "auto") -> str:
    """
    Extract organism name from a FASTA description based on HEADER_MODE.

    Supported forms:
      - Key/val tags: [organism=Myotis davidii] (preferred)
      - Legacy: [Myotis davidii] (no key)
    """
    desc = desc or ""
    tags = parse_bracket_keyvals(desc)

    if mode in ("auto", "tag"):
        for key in ("organism", "species", "taxname", "taxon"):
            if key in tags and tags[key]:
                return tags[key]

        if mode == "tag":
            return ""  # strictly tags only

    if mode in ("auto", "legacy"):
        # legacy: find first bracket chunk that *looks like* "Genus species"
        for chunk in BRACKET_TAG_RE.findall(desc):
            chunk = chunk.strip()
            if "=" in chunk:
                continue
            # basic heuristic: "Genus species" (two tokens, first capitalized)
            parts = chunk.split()
            if len(parts) >= 2 and parts[0][:1].isupper():
                return " ".join(parts[:2])

        if mode == "legacy":
            return ""

    # substring mode doesn't extract; caller must supply organism separately
    return ""

def organism_matches(protein_org: str, transcript_desc: str, mode: str = "auto") -> bool:
    """
    Decide whether transcript belongs to same organism as protein.
    Case-insensitive matching.

    - If protein_org is empty -> no filtering (True)
    - If mode uses tags/legacy, try to extract transcript organism too; if present, require equality.
    - Otherwise fall back to substring check.
    """
    if not protein_org:
        return True

    p = protein_org.strip().lower()
    t_desc = transcript_desc or ""

    if mode in ("auto", "tag", "legacy"):
        t_org = extract_organism_from_desc(t_desc, mode=mode).strip().lower()
        if t_org:
            return t_org == p
        if mode == "tag":
            # strict tag mode: if transcript has no tag, don't filter it out
            return True
        # auto/legacy: if transcript org missing, fall back to substring
        return p in t_desc.lower()

    # substring mode
    return p in t_desc.lower()

# =============================================================================
# Helper functions from original script
# =============================================================================
def already_in_results(transcript_desc: str) -> bool:
    """
    Original behavior: prevent duplicates by transcript description identity.
    Keep this for backwards compatibility.
    """
    global results
    for record in results:
        if transcript_desc == record.description:
            return True
    return False

# =============================================================================
# Core matching logic
# =============================================================================
def Process(protein_desc: str, protein_seq, transcripts_fasta: str, protein_org: str):
    """
    Scan transcript FASTA. For transcripts that match organism, slide window until translation matches protein.
    If KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM is True and multiple transcripts exist for the same organism,
    only the FIRST transcript (X1) is eligible; the rest are ignored.
    """
    global logfile

    NT_SEQ_LENGTH = len(protein_seq) * 3

    # Track organism -> first transcript already seen
    seen_organisms = set()

    with open(transcripts_fasta, "r") as transcript_handle:
        for transcript_record in SeqIO.parse(transcript_handle, "fasta"):
            transcript_id   = transcript_record.id
            transcript_desc = transcript_record.description
            transcript_seq  = transcript_record.seq

            # Filter by organism (protein vs transcript)
            if not organism_matches(protein_org, transcript_desc, mode=HEADER_MODE):
                continue

            # Determine transcript organism key for "keep only X1"
            tx_org = extract_organism_from_desc(transcript_desc, mode=HEADER_MODE)
            if not tx_org and protein_org:
                tx_org = protein_org  # fallback to protein organism if transcript doesn't encode it
            tx_org_key = (tx_org or "").strip().lower()

            # Enforce transcript=X1 per organism
            if KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM and tx_org_key:
                if tx_org_key in seen_organisms:
                    continue
                seen_organisms.add(tx_org_key)

            # Slide window across transcript
            start = 0
            DONE = False
            coding_seq = None

            while start < len(str(transcript_seq)):
                try:
                    candidate = transcript_seq[start:start + NT_SEQ_LENGTH]
                    translated = candidate.translate()
                except Exception:
                    translated = None

                # original duplicate guard
                exists = already_in_results(transcript_desc)

                if translated is not None and str(translated) == str(protein_seq) and not exists:
                    DONE = True
                    coding_seq = candidate
                    break
                else:
                    start += 1

            # If enforcing X1, and we just tried X1 and it didn't match, do NOT try later transcripts
            if KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM and tx_org_key:
                # we have consumed X1 for that organism; stop searching other transcripts for this protein/org
                # because requirement says keep only transcript=X1
                if DONE:
                    break
                else:
                    # no match in X1; return NO_MATCH immediately
                    return "NO_MATCH"

            if DONE:
                break

    if DONE:
        # Preserve headers: create a NEW record so we don't lose provenance
        out_id = f"{transcript_id}"
        out_desc = (
            f"protein={protein_desc} | "
            f"transcript={transcript_desc} | "
            f"organism={protein_org}"
        )
        return SeqRecord(seq=coding_seq, id=out_id, description=out_desc)
    else:
        return "NO_MATCH"

# =============================================================================
# Main subroutine (sanity)
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

# Compute average protein length (AA) and corresponding NT length (kept from original)
prot_seq_lengths = []
with open(PROTEIN_FASTA, "r") as prot_handle:
    for record in SeqIO.parse(prot_handle, "fasta"):
        prot_seq_lengths.append(len(record.seq))
avg_sequence_length = sum(prot_seq_lengths) / len(prot_seq_lengths) if prot_seq_lengths else 0
avg_sequence_length_nt = avg_sequence_length * 3

log(f"# HEADER_MODE={HEADER_MODE}")
log(f"# KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM={KEEP_ONLY_FIRST_TRANSCRIPT_PER_ORGANISM}")
log(f"# Average sequence length is (PROTEIN AA ungapped): {avg_sequence_length}")
log(f"# Average sequence length is (NUCLEOTIDE NUC): {avg_sequence_length_nt}")

# Iterate proteins
with open(PROTEIN_FASTA, "r") as prot_handle:
    for record in SeqIO.parse(prot_handle, "fasta"):
        protein_id   = record.id
        protein_desc = record.description
        protein_seq  = record.seq

        log("\n")
        log(f"# Processing: {protein_desc}")

        if "LOW QUALITY PROTEIN" in str(protein_desc) or "partial" in str(protein_desc).lower():
            log(f"# Skipping this sequence due to quality issues: {protein_desc}")
            continue

        # Extract organism/species from protein header using selectable mode
        protein_org = ""
        if HEADER_MODE == "substring":
            # substring mode expects organism to be embedded; best effort: try auto anyway
            protein_org = extract_organism_from_desc(protein_desc, mode="auto")
        else:
            protein_org = extract_organism_from_desc(protein_desc, mode=HEADER_MODE)

        log(f"# Organism: {protein_org}")

        # Heavy lifting here (maintain original calling style/compatibility)
        tx_record = Process(protein_desc, protein_seq, TRANSCRIPTS_FASTA, protein_org)

        if type(tx_record) != str:
            log(f"# Match: {tx_record.description}\n")
            results.append(tx_record)
        else:
            log("# -- NO Match --\n")
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
