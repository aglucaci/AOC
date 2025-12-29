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

ACCESSION_DOTTED = re.compile(r'([A-Z]{1,3}_[0-9]+\.[0-9]+)')
ACCESSION_UNDERSCORE_VER = re.compile(r'([A-Z]{1,3}_[0-9]+)_([0-9]+)(?=(_|$))')
ACCESSION_BASE = re.compile(r'([A-Z]{1,3}_[0-9]+)')

_last_entrez_call_t = 0.0

def log(msg: str, verbose: bool):
    if verbose:
        print(msg, file=sys.stderr)

def eprint(msg: str):
    print(msg, file=sys.stderr)

def throttle(min_interval_s: float, verbose: bool):
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
    candidates = [accession]
    if "." in accession:
        candidates.append(accession.split(".", 1)[0])

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
                if e.code == 429:
                    backoff = max(5.0, sleep_s * attempt * 5)
                    eprint(f"[WARN] resolve failed for {cand} (attempt {attempt}/{retries}) HTTP 429; backing off {backoff:.1f}s")
                    time.sleep(backoff)
                    continue
                eprint(f"[WARN] resolve failed for {cand} (attempt {attempt}/{retries}) HTTP {e.code}: {e.reason}")
            except Exception as e:
                last_err = e
                eprint(f"[WARN] resolve failed for {cand} (attempt {attempt}/{retries}): {e}")

            time.sleep(sleep_s * attempt)

    eprint(f"[ERROR] could not resolve {accession}. last error: {last_err}")
    return None, None, None

def load_map_id_to_transcript(map_tsv: str) -> Dict[str, str]:
    """
    map.tsv must contain columns:
      - id30
      - transcript_id

    Tree software often sanitizes leaf names (replacing '.' and '|' with '_').
    So we build a lookup that accepts BOTH:
      - id30 as-is
      - sanitized(id30): '.'->'_' and '|'->'_'
    """
    df = pd.read_csv(map_tsv, sep="\t", dtype=str).fillna("")
    if "id30" not in df.columns:
        raise ValueError(f"--map_tsv missing 'id30'. Columns={list(df.columns)}")
    if "transcript_id" not in df.columns:
        raise ValueError(f"--map_tsv missing 'transcript_id'. Columns={list(df.columns)}")

    m: Dict[str, str] = {}

    for _, r in df.iterrows():
        id30 = str(r["id30"]).strip()
        txid = str(r["transcript_id"]).strip()
        if not id30 or not txid:
            continue

        # 1) exact key
        if id30 not in m:
            m[id30] = txid

        # 2) sanitized key (matches IQ-TREE/others)
        safe = id30.replace("|", "_").replace(".", "_")
        if safe not in m:
            m[safe] = txid

    return m



def read_orthologs_robust(path: str, verbose: bool) -> Optional[pd.DataFrame]:
    """
    Orthologs is optional when --map_tsv is used.
    When provided, read it robustly to avoid crashing on bad lines.
    """
    if not path:
        return None

    try:
        df = pd.read_csv(
            path,
            dtype=str,
            engine="python",         # more tolerant than C engine
            on_bad_lines="warn"      # or "skip" if you prefer silence
        )
        return df
    except Exception as e:
        eprint(f"[WARN] Failed to read --orthologs='{path}': {e}")
        eprint("[WARN] Continuing without orthologs table (it is optional when --map_tsv is provided).")
        if verbose:
            eprint("[HINT] If orthologs is TSV, try converting or pass a clean delimiter.")
        return None

def main():
    ap = argparse.ArgumentParser()

    # orthologs becomes optional
    ap.add_argument("--orthologs", default="",
                    help="Optional. A CSV/TSV ortholog table. Not required if --map_tsv is provided.")
    ap.add_argument("--tree", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--accession_col", default=ACCESSION_COL_DEFAULT)

    ap.add_argument("--map_tsv", default="",
                    help="Optional: codons output .map.tsv with columns 'id30' and 'transcript_id'. "
                         "If provided, tree leaf labels are interpreted as id30 and mapped to transcript_id.")

    ap.add_argument("--db", default="nuccore")
    ap.add_argument("--email", default="")
    ap.add_argument("--api_key", default="")
    ap.add_argument("--retries", type=int, default=8)
    ap.add_argument("--sleep", type=float, default=0.6)
    ap.add_argument("--min_interval", type=float, default=0.4)
    ap.add_argument("--allow_unknown", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    Entrez.tool = "annotate_lineages"
    if not args.email:
        eprint("[WARN] --email not provided. NCBI may throttle/reject; provide --email (and ideally --api_key).")
    else:
        Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # orthologs is optional now
    df_orth = read_orthologs_robust(args.orthologs, args.verbose)
    if df_orth is not None and args.accession_col not in df_orth.columns:
        eprint(f"[WARN] --accession_col '{args.accession_col}' not found in orthologs columns={list(df_orth.columns)}")
        eprint("[WARN] Ortholog table will not be used for accession extraction (tree/map will drive accessions).")

    with open(args.tree, "r") as f:
        newick = f.read().strip()
    leaves = parse_tree_leaves(newick)

    id_to_tx: Dict[str, str] = {}
    if args.map_tsv:
        id_to_tx = load_map_id_to_transcript(args.map_tsv)
        log(f"[map] loaded {len(id_to_tx)} id30->transcript_id mappings", args.verbose)

    if not args.map_tsv and not leaves:
        raise RuntimeError("No leaves found in tree and no --map_tsv provided.")

    cache: Dict[str, Tuple[Optional[str], Optional[str], Optional[List[str]]]] = {}
    rows = []

    for leaf in leaves:
        if args.map_tsv:
            acc = id_to_tx.get(leaf, "")
            if not acc:
                msg = f"No transcript_id mapping for leaf '{leaf}' in map file {args.map_tsv}"
                if args.allow_unknown:
                    rows.append({"leaf": leaf, "accession": "", "organism": "Unknown", "taxid": "", "lineage": "Unknown"})
                    continue
                raise RuntimeError(msg)
            acc = normalize_accession(acc, verbose=args.verbose)
        else:
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
            msg = f"Failed to resolve lineage for accession {acc}."
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
