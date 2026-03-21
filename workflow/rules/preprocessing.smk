rule macse:
    input:
        fasta=lambda wc: codon_fasta_path(wc.sample)
    output:
        codons=os.path.join(str(OUTROOT), "{sample}", "{sample}.codons.fa"),
        aa=os.path.join(str(OUTROOT), "{sample}", "{sample}.aa.fa"),
    params:
        launcher=MACSE_LAUNCHER,
        runner=MACSE_RUNNER,
        xms_mb=MACSE_XMS_MB,
        xmx_mb=_macse_heap_mb(),
    resources:
        mem_mb=_macse_requested_mem_mb(),
    log:
        os.path.join(LOGDIR, "macse.{sample}.log")
    shell:
        r"""
        mkdir -p {OUTROOT}/{wildcards.sample}
        bash "{params.runner}" \
          --macse-bin "{params.launcher}" \
          --xms-mb {params.xms_mb} \
          --xmx-mb {params.xmx_mb} \
          -- \
          -prog alignSequences \
          -seq "{input.fasta}" \
          -out_NT "{output.codons}" \
          -out_AA "{output.aa}" \
          -max_refine_iter 3 \
          -local_realign_init 0.3 \
          -local_realign_dec 0.2 \
          > {log} 2>&1
        """


rule cln:
    input:
        codons=rules.macse.output.codons
    output:
        cln=os.path.join(str(OUTROOT), "{sample}", "{sample}.codons.cln.fa")
    log:
        os.path.join(LOGDIR, "cln.{sample}.log")
    shell:
        r"""
        {HYPHY} CLN Universal {input.codons} 'No/No' {output.cln} > {log} 2>&1
        """


rule strike_ambigs_msa:
    input:
        cln=rules.cln.output.cln
    output:
        sa=os.path.join(str(OUTROOT), "{sample}", "{sample}.SA.codons.cln.fa")
    log:
        os.path.join(LOGDIR, "strike_ambigs_msa.{sample}.log")
    shell:
        r"""
        {HYPHY} {STRIKE_AMBIGS_BF} --alignment {input.cln} --output {output.sa} > {log} 2>&1
        """


rule remove_duplicates_msa:
    input:
        sa=rules.strike_ambigs_msa.output.sa
    output:
        rd=os.path.join(str(OUTROOT), "{sample}", "{sample}.RD.SA.codons.cln.fa")
    log:
        os.path.join(LOGDIR, "remove_duplicates_msa.{sample}.log")
    shell:
        r"""
        {HYPHY} {REMOVE_DUPS_BF} --msa {input.sa} --output {output.rd} ENV='DATA_FILE_PRINT_FORMAT=9' > {log} 2>&1
        """


rule tn93_cluster:
    input:
        rd=rules.remove_duplicates_msa.output.rd
    output:
        cluster_json=os.path.join(str(OUTROOT), "{sample}", "{sample}.RD.SA.codons.cln.fa.cluster.json"),
        cluster_fasta=os.path.join(str(OUTROOT), "{sample}", "{sample}.RD.SA.codons.cln.fa.cluster.fasta"),
    params:
        threshold=lambda wc: config.get("tn93_threshold", 0.01),
        max_seqs=lambda wc: config.get("tn93_max_seqs", 20),
    log:
        os.path.join(LOGDIR, "tn93_cluster.{sample}.log")
    shell:
        r"""
        python {TN93_CLUSTER_PY} \
          --input {input.rd} \
          --output_fasta {output.cluster_fasta} \
          --output_json {output.cluster_json} \
          --threshold {params.threshold} \
          --max_retain {params.max_seqs} \
          > {log} 2>&1
        """


rule recombination:
    input:
        msa=rules.tn93_cluster.output.cluster_fasta
    output:
        gard_json=os.path.join(str(OUTROOT), "{sample}", "{sample}.RD.SA.codons.cln.fa.cluster.fasta.GARD.json"),
        bestgard=os.path.join(str(OUTROOT), "{sample}", "{sample}.RD.SA.codons.cln.fa.cluster.fasta.best-gard"),
    threads:
        _gard_mpi_procs()
    log:
        os.path.join(LOGDIR, "recombination.{sample}.log")
    run:
        procs = int(threads)
        if str(GARD_LAUNCHER).lower() == "hyphympi":
            cmd = (
                f"{GARD_MPI_RUNNER} -np {procs} {GARD_LAUNCHER} "
                f"GARD --alignment {input.msa} --rv GDD --output {output.gard_json} "
                f"--mode Faster ENV=TOLERATE_NUMERICAL_ERRORS=1"
            )
        else:
            cmd = (
                f"{GARD_LAUNCHER} GARD --alignment {input.msa} --rv GDD "
                f"--output {output.gard_json} --mode Faster "
                f"ENV=TOLERATE_NUMERICAL_ERRORS=1"
            )
        shell(cmd + f" > {log} 2>&1")


rule build_header_map_and_test_list:
    input:
        labels=lambda wc: labels_csv_path(wc.sample)
    output:
        map_csv=os.path.join(str(OUTROOT), "{sample}", "{sample}.sequence_header_map.csv"),
        test_txt=os.path.join(str(OUTROOT), "{sample}", "{sample}.test_sequences.txt"),
    log:
        os.path.join(LOGDIR, "build_header_map_and_test_list.{sample}.log")
    run:
        import os, re
        import pandas as pd
        from collections import Counter

        os.makedirs(os.path.join(str(OUTROOT), wildcards.sample), exist_ok=True)
        logp = log[0] if hasattr(log, "__iter__") else str(log)
        os.makedirs(os.path.dirname(logp), exist_ok=True)

        def sanitize_header(h):
            h = (h or "").strip().lstrip(">")
            h = re.sub(r"[^A-Za-z0-9]+", "_", h)
            h = re.sub(r"_+", "_", h)
            return h.strip("_")

        with open(logp, "w") as handle:
            handle.write(f"[INFO] sample={wildcards.sample}\n")
            handle.write(f"[INFO] labels={input.labels}\n")

        if (not input.labels) or (str(input.labels).strip() == ""):
            raise ValueError(f"[build_header_map_and_test_list] labels path is empty for sample={wildcards.sample}")
        if not os.path.exists(input.labels):
            raise FileNotFoundError(f"[build_header_map_and_test_list] labels CSV not found: {input.labels}")

        df = pd.read_csv(input.labels)
        df.columns = [c.strip().lower() for c in df.columns]

        label_candidates = ["label", "group", "set", "class", "foreground", "branchset"]
        header_candidates = ["header", "sequence_header", "seq_header", "id", "name", "taxon", "tip", "sequence"]

        label_col = next((c for c in label_candidates if c in df.columns), None)
        header_col = next((c for c in header_candidates if c in df.columns), None)

        if label_col is None or header_col is None:
            if df.shape[1] == 2:
                label_col = df.columns[0]
                header_col = df.columns[1]
            else:
                raise ValueError(
                    f"[build_header_map_and_test_list] Could not infer label/header columns.\n"
                    f"Columns found: {list(df.columns)}\n"
                    f"Expected something like LABEL + HEADER (or a 2-column CSV)."
                )

        labels_raw = df[label_col].astype(str).str.strip()
        headers_raw = df[header_col].astype(str).str.strip().str.lstrip(">")
        mapped = headers_raw.map(sanitize_header)

        counts = Counter()
        mapped_unique = []
        for header in mapped.tolist():
            counts[header] += 1
            mapped_unique.append(header if counts[header] == 1 else f"{header}_dup{counts[header]}")

        mapped_unique = pd.Series(mapped_unique)
        labels_norm = labels_raw.str.lower()
        is_test = labels_norm.isin({"test", "foreground", "fg", "case"}) | labels_norm.str.contains(r"\btest\b", regex=True)
        test_headers_mapped = mapped_unique[is_test].dropna().tolist()

        map_df = pd.DataFrame({"original_header": headers_raw, "mapped_header": mapped_unique})
        map_df.to_csv(output.map_csv, index=False)

        with open(output.test_txt, "w") as out:
            out.write("\n".join(test_headers_mapped) + ("\n" if test_headers_mapped else ""))

        with open(logp, "a") as handle:
            handle.write(f"[INFO] label_col={label_col}, header_col={header_col}\n")
            handle.write(f"[INFO] total_rows={len(df)}\n")
            handle.write(f"[INFO] test_rows={len(test_headers_mapped)}\n")

        if os.path.getsize(output.map_csv) == 0:
            raise RuntimeError(f"[build_header_map_and_test_list] map_csv empty: {output.map_csv}")
        if config.get("require_test_sequences", False) and len(test_headers_mapped) == 0:
            raise RuntimeError(
                f"[build_header_map_and_test_list] No Test sequences found for sample={wildcards.sample} "
                f"but require_test_sequences=True. See log: {logp}"
            )


checkpoint parse_gard:
    input:
        bestgard=lambda wc: rules.recombination.output.bestgard if _recombination_enabled(wc.sample) else [],
        msa=rules.remove_duplicates_msa.output.rd
    output:
        manifest=os.path.join(str(OUTROOT), "{sample}", "gard", "segments.json"),
        segdir=directory(os.path.join(str(OUTROOT), "{sample}", "gard", "segments")),
        log=os.path.join(str(OUTROOT), "{sample}", "{sample}.gard.log"),
    run:
        os.makedirs(os.path.dirname(output.log), exist_ok=True)

        records = list(SeqIO.parse(input.msa, "fasta"))
        if not records:
            raise ValueError(f"[ERROR] Empty MSA: {input.msa}")

        first_seq = str(records[0].seq)

        if _recombination_enabled(wildcards.sample):
            if not input.bestgard or not os.path.exists(input.bestgard):
                raise FileNotFoundError(f"[ERROR] GARD results file does not exist: {input.bestgard}")

            data = [line.strip() for line in open(input.bestgard) if "CHARSET" in line]
            if not data:
                raise ValueError("[ERROR] GARD results file is empty or improperly formatted (no CHARSET lines).")

            data_parsed = [[int(x) for x in line.split("=")[1].strip(" ;").split("-")] for line in data]
            data_adjusted = [[start - 1, end] for start, end in data_parsed]

            data_adjusted_coords = []
            carry_over = 0
            for start, end in data_adjusted:
                start = start - carry_over
                seg_len = len(first_seq[int(start): int(end)])
                remainder = seg_len % 3
                data_adjusted_coords.append([start, end - remainder])
                carry_over = remainder
        else:
            trimmed_len = len(first_seq) - (len(first_seq) % 3)
            data_adjusted_coords = [[0, trimmed_len]]

        os.makedirs(gard_segments_dir(wildcards.sample), exist_ok=True)
        manifest = {
            "sample": wildcards.sample,
            "msa": str(input.msa),
            "bestgard": str(input.bestgard) if input.bestgard else None,
            "recombination_enabled": _recombination_enabled(wildcards.sample),
            "n_parts": len(data_adjusted_coords),
            "parts": [],
        }

        with open(output.log, "w") as log_file:
            if not _recombination_enabled(wildcards.sample):
                log_file.write("Recombination disabled; using full clustered alignment as a single partition.\n")
            for i, (start, stop) in enumerate(data_adjusted_coords, start=1):
                out_fa = os.path.join(gard_segments_dir(wildcards.sample), f"{wildcards.sample}.part{i}.codon.fas")
                with open(out_fa, "w") as out_f:
                    for record in records:
                        partition = record[int(start):int(stop)]
                        out_f.write(f">{record.id}\n{partition.seq}\n")

                adjusted_len = len(records[0][int(start):int(stop)].seq)
                log_file.write(
                    f"Partition {i}: start={start}, stop={stop}, length={adjusted_len}, remainder={adjusted_len % 3}\n"
                )

                manifest["parts"].append({
                    "part": i,
                    "start0": int(start),
                    "end_excl": int(stop),
                    "fasta": out_fa,
                    "length": int(adjusted_len),
                })

        with open(output.manifest, "w") as fh:
            json.dump(manifest, fh, indent=2)
