#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import time
import re
from io import StringIO
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError

import pandas as pd
from Bio import Entrez, Phylo, SeqIO
from Bio.SeqIO import SeqRecord

ACCESSION_COL_DEFAULT = "RefSeq Transcript accessions"
# Recognize accessions embedded in messy leaf labels.
# Handles:
#   NM_000546.6
#   NM_173527_3_Homo_sapiens_...   -> NM_173527.3
#   ref|NM_000546.6|TP53
ACCESSION_DOTTED = re.compile(r'([A-Z]{1,3}_[0-9]+\.[0-9]+)')
# Version encoded as underscore and followed by underscore or end (IMPORTANT: not a word boundary)
ACCESSION_UNDERSCORE_VER = re.compile(r'([A-Z]{1,3}_[0-9]+)_([0-9]+)(?=(_|$))')
ACCESSION_BASE = re.compile(r'([A-Z]{1,3}_[0-9]+)')

_last_entrez_call_t = 0.0

def log(msg: str, verbose: bool):
    if verbose:
        print(msg, file=sys.stderr)

def eprint(msg: str):
    print(msg, file=sys.stderr)

def throttle(min_interval_s: float, verbose: bool):
    """Global throttle to avoid hammering NCBI (esp. when no API key)."""
    global _last_entrez_call_t
    if min_interval_s <= 0:
        return
    now = time.time()
    dt = now - _last_entrez_call_t
    if dt < min_interval_s:
        sleep_for = min_interval_s - dt
        log(f"[throttle] sleeping {sleep_for:.3f}s", verbose)
        time.sleep(sleep_for)
    _last_entrez_call_t = time.time()

def normalize_accession(label: str, verbose: bool=False) -> str:
    """
    Extract a valid accession from a tree leaf label.

    Critical fix:
      The prior regex used a trailing \\b, which FAILS on strings like
      'NM_001731_3_Homo_sapiens_...' because '_' is a word character,
      so there's no word boundary after the version digit.

    This version correctly converts:
      NM_173527_3_Homo_sapiens_... -> NM_173527.3
    """
    s = str(label).strip()
    if not s or s.lower() == "nan":
        return ""

    m = ACCESSION_DOTTED.search(s)
    if m:
        acc = m.group(1)
        log(f"[normalize] '{label}' -> '{acc}' (dotted)", verbose)
        return acc

    m = ACCESSION_UNDERSCORE_VER.search(s)
    if m:
        acc = f"{m.group(1)}.{m.group(2)}"
        log(f"[normalize] '{label}' -> '{acc}' (underscore-version)", verbose)
        return acc

    m = ACCESSION_BASE.search(s)
    if m:
        acc = m.group(1)
        log(f"[normalize] '{label}' -> '{acc}' (base)", verbose)
        return acc

    # Last resort: sanitize first token
    tok = s.split()[0].strip().strip(",;|")
    tok = re.sub(r"[^A-Za-z0-9_.]", "", tok)
    log(f"[normalize] '{label}' -> '{tok}' (fallback)", verbose)
    return tok

def parse_tree_leaves(newick: str) -> List[str]:
    tree = Phylo.read(StringIO(newick), "newick")
    return [t.name for t in tree.get_terminals() if t.name]

def _efetch_genbank(accession: str, db: str, verbose: bool, min_interval_s: float) -> SeqRecord:
    throttle(min_interval_s, verbose)
    log(f"[efetch] db={db} id={accession}", verbose)
    handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
    try:
        rec = SeqIO.read(handle, "genbank")
    finally:
        handle.close()
    return rec

def fetch_taxid_organism_via_efetch(accession: str, db: str, verbose: bool, min_interval_s: float) -> Tuple[Optional[str], Optional[str]]:
    rec = _efetch_genbank(accession, db=db, verbose=verbose, min_interval_s=min_interval_s)

    organism = None
    if hasattr(rec, "annotations") and isinstance(rec.annotations, dict):
        organism = rec.annotations.get("organism")

    taxid = None
    for f in getattr(rec, "features", []):
        if f.type != "source":
            continue
        for x in f.qualifiers.get("db_xref", []):
            if isinstance(x, str) and x.startswith("taxon:"):
                taxid = x.split("taxon:", 1)[1].strip()
                break
        if taxid:
            break

    return organism, taxid

def fetch_taxid_organism_via_esummary(accession: str, db: str, verbose: bool, min_interval_s: float) -> Tuple[Optional[str], Optional[str]]:
    throttle(min_interval_s, verbose)
    log(f"[esummary] db={db} id={accession}", verbose)
    handle = Entrez.esummary(db=db, id=accession, retmode="xml")
    try:
        rec = Entrez.read(handle)
    finally:
        handle.close()

    if not rec:
        return None, None

    doc = rec[0]
    organism = doc.get("Organism")
    taxid = doc.get("TaxId")
    return organism, taxid

def fetch_lineage_from_taxid(taxid: str, verbose: bool, min_interval_s: float) -> Optional[List[str]]:
    if not taxid:
        return None
    throttle(min_interval_s, verbose)
    log(f"[taxonomy] taxid={taxid}", verbose)
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    try:
        rec = Entrez.read(handle)
    finally:
        handle.close()

    if not rec:
        return None
    r0 = rec[0]
    lineage = []
    lin = r0.get("Lineage", "")
    if lin:
        lineage.extend([x.strip() for x in lin.split(";") if x.strip()])
    sci = r0.get("ScientificName", "")
    if sci:
        lineage.append(sci)
    return lineage or None

def resolve_accession(accession: str, db: str, retries: int, sleep_s: float, verbose: bool, min_interval_s: float) -> Tuple[Optional[str], Optional[str], Optional[List[str]]]:
    """
    Resolve accession -> (organism, taxid, lineage)

    Candidates tried:
      - as given (e.g., NM_173527.3)
      - without version suffix (NM_173527)
    """
    candidates = [accession]
    if "." in accession:
        candidates.append(accession.split(".", 1)[0])

    # De-dup preserve order
    seen = set()
    candidates = [c for c in candidates if c and not (c in seen or seen.add(c))]

    last_err: Optional[Exception] = None

    for cand in candidates:
        for attempt in range(1, retries + 1):
            try:
                log(f"[resolve] cand={cand} attempt={attempt}/{retries} db={db}", verbose)

                organism, taxid = fetch_taxid_organism_via_efetch(
                    cand, db=db, verbose=verbose, min_interval_s=min_interval_s
                )

                if not taxid:
                    org2, tax2 = fetch_taxid_organism_via_esummary(
                        cand, db=db, verbose=verbose, min_interval_s=min_interval_s
                    )
                    organism = organism or org2
                    taxid = taxid or tax2

                lineage = fetch_lineage_from_taxid(
                    taxid, verbose=verbose, min_interval_s=min_interval_s
                ) if taxid else None

                log(f"[resolved] cand={cand} taxid={taxid} organism={organism}", verbose)
                return organism, taxid, lineage

            except HTTPError as e:
                last_err = e
                # 429: too many requests -> back off harder
                if e.code == 429:
                    backoff = max(5.0, sleep_s * attempt * 5)
                    eprint(f"[WARN] resolve failed for {cand} (attempt {attempt}/{retries}) HTTP 429: Too Many Requests; backing off {backoff:.1f}s")
                    time.sleep(backoff)
                    continue
                eprint(f"[WARN] resolve failed for {cand} (attempt {attempt}/{retries}) HTTP {e.code}: {e.reason}")
            except Exception as e:
                last_err = e
                eprint(f"[WARN] resolve failed for {cand} (attempt {attempt}/{retries}): {e}")

            time.sleep(sleep_s * attempt)

    eprint(f"[ERROR] could not resolve {accession}. last error: {last_err}")
    return None, None, None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthologs", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--accession_col", default=ACCESSION_COL_DEFAULT)
    ap.add_argument("--db", default="nuccore")
    ap.add_argument("--email", default="")
    ap.add_argument("--api_key", default="")
    ap.add_argument("--retries", type=int, default=8)
    ap.add_argument("--sleep", type=float, default=0.6)
    ap.add_argument("--min_interval", type=float, default=0.4,
                    help="Minimum seconds between Entrez requests (throttle). "
                         "Use ~0.4 without API key; with API key you can lower (e.g., 0.12).")
    ap.add_argument("--allow_unknown", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    Entrez.tool = "annotate_lineages"
    if not args.email:
        eprint("[WARN] --email not provided. NCBI may reject or throttle anonymous requests; provide --email (and ideally --api_key).")
    else:
        Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    df = pd.read_csv(args.orthologs)
    if args.accession_col not in df.columns:
        raise ValueError(f"Accession column '{args.accession_col}' not found. Columns={list(df.columns)}")

    with open(args.tree, "r") as f:
        newick = f.read().strip()
    leaves = parse_tree_leaves(newick)

    cache: Dict[str, Tuple[Optional[str], Optional[str], Optional[List[str]]]] = {}
    rows = []

    for leaf in leaves:
        acc = normalize_accession(leaf, verbose=args.verbose)
        if not acc:
            msg = f"Empty accession after normalization for leaf '{leaf}'"
            if args.allow_unknown:
                rows.append({"leaf": leaf, "accession": "", "organism": "Unknown", "taxid": "", "lineage": "Unknown"})
                continue
            raise RuntimeError(msg)

        if acc not in cache:
            cache[acc] = resolve_accession(
                acc, db=args.db, retries=args.retries, sleep_s=args.sleep,
                verbose=args.verbose, min_interval_s=args.min_interval
            )

        organism, taxid, lineage = cache[acc]
        if not taxid or not lineage:
            msg = f"Failed to resolve lineage for accession {acc}. Provide --email/--api_key or fix IDs."
            if args.allow_unknown:
                rows.append({"leaf": leaf, "accession": acc, "organism": organism or "Unknown", "taxid": taxid or "", "lineage": "Unknown"})
                continue
            raise RuntimeError(msg)

        rows.append({
            "leaf": leaf,
            "accession": acc,
            "organism": organism or "",
            "taxid": str(taxid),
            "lineage": ";".join(lineage),
        })

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
