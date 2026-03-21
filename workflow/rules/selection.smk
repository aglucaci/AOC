def _parts_for_sample(wc):
    ck = checkpoints.parse_gard.get(sample=wc.sample)
    with open(ck.output.manifest) as fh:
        manifest = json.load(fh)
    return [str(p["part"]) for p in manifest["parts"]]


def _partition_fasta(wc):
    return os.path.join(str(OUTROOT), wc.sample, "gard", "segments", f"{wc.sample}.part{wc.part}.codon.fas")


rule sanitize_partition_alignment:
    input:
        aln=_partition_fasta
    output:
        aln=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "partition.filtered.codon.fas")
    log:
        os.path.join(LOGDIR, "sanitize_partition_alignment.{sample}.part{part}.log")
    run:
        records = list(SeqIO.parse(input.aln, "fasta"))
        if not records:
            raise ValueError(f"[sanitize_partition_alignment] Empty partition alignment: {input.aln}")

        allowed = set("ACGTURYKMSWBDHVN")
        kept = []
        dropped = []
        for record in records:
            seq = str(record.seq).upper()
            informative = [ch for ch in seq if ch not in {"-", "?", "."}]
            if any(ch in allowed for ch in informative):
                kept.append(record)
            else:
                dropped.append(record.id)

        if len(kept) < 3:
            raise ValueError(
                f"[sanitize_partition_alignment] Fewer than 3 taxa remain after filtering all-gap/all-ambiguity sequences "
                f"for {wildcards.sample} part {wildcards.part}."
            )

        os.makedirs(os.path.dirname(output.aln), exist_ok=True)
        with open(output.aln, "w") as handle:
            SeqIO.write(kept, handle, "fasta")

        with open(log[0] if hasattr(log, "__iter__") else str(log), "w") as handle:
            handle.write(f"input_records={len(records)}\n")
            handle.write(f"kept_records={len(kept)}\n")
            handle.write(f"dropped_records={len(dropped)}\n")
            for taxon in dropped:
                handle.write(f"dropped\t{taxon}\n")


rule iqtree_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln
    output:
        tree=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "IQ-TREE.treefile")
    params:
        prefix=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "IQ-TREE")
    log:
        os.path.join(LOGDIR, "iqtree_partition.{sample}.part{part}.log")
    shell:
        r"""
        mkdir -p {OUTROOT}/{wildcards.sample}/selection/part{wildcards.part}
        {IQTREE} \
          -s {input.aln} \
          -st DNA \
          -m GTR+G \
          -nt AUTO \
          --prefix {params.prefix} \
          --redo \
          > {log} 2>&1
        """


rule label_tree_TestForeground_partition:
    input:
        tree=rules.iqtree_partition.output.tree,
        test_list=_test_list_for_sample
    output:
        tree=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "IQ-TREE.labelled.treefile")
    log:
        os.path.join(LOGDIR, "label_tree_TestForeground_partition.{sample}.part{part}.log")
    run:
        os.makedirs(os.path.dirname(output.tree), exist_ok=True)
        logp = log[0] if hasattr(log, "__iter__") else str(log)
        with open(logp, "a") as handle:
            handle.write("START\n")
        if (not os.path.exists(input.test_list)) or (os.path.getsize(input.test_list) == 0):
            shell("cp {input.tree} {output.tree}")
            with open(logp, "a") as handle:
                handle.write("No test list; copied tree\nEND\n")
        else:
            cmd = (
                f"{HYPHY} {LABEL_TREE_BF} "
                f"--tree {input.tree} "
                f"--list {input.test_list} "
                f"--output {output.tree} "
                f"--label Test"
            )
            shell(cmd + " > " + logp + " 2>&1")
            with open(logp, "a") as handle:
                handle.write("END\n")


rule FEL_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln,
        tree=rules.label_tree_TestForeground_partition.output.tree
    output:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "FEL.json")
    params:
        branch_arg=lambda wc: _branch_arg(wc.sample)
    log:
        os.path.join(LOGDIR, "FEL_partition.{sample}.part{part}.log")
    shell:
        r"""
        mkdir -p {OUTROOT}/{wildcards.sample}/selection/part{wildcards.part}
        {HYPHY} FEL --alignment {input.aln} --tree {input.tree} --output {output.json} --ci Yes {params.branch_arg} > {log} 2>&1
        """


rule MEME_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln,
        tree=rules.label_tree_TestForeground_partition.output.tree
    output:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "MEME.json")
    params:
        branch_arg=lambda wc: _branch_arg(wc.sample)
    log:
        os.path.join(LOGDIR, "MEME_partition.{sample}.part{part}.log")
    shell:
        r"""
        mkdir -p {OUTROOT}/{wildcards.sample}/selection/part{wildcards.part}
        {HYPHY} MEME --alignment {input.aln} --tree {input.tree} --output {output.json} \
          {params.branch_arg} ENV=TOLERATE_NUMERICAL_ERRORS=1 > {log} 2>&1
        """


rule ABSREL_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln,
        tree=rules.label_tree_TestForeground_partition.output.tree
    output:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "ABSREL.json")
    params:
        branch_arg=lambda wc: _branch_arg(wc.sample)
    log:
        os.path.join(LOGDIR, "ABSREL_partition.{sample}.part{part}.log")
    shell:
        r"""
        mkdir -p {OUTROOT}/{wildcards.sample}/selection/part{wildcards.part}
        {HYPHY} ABSREL --alignment {input.aln} --tree {input.tree} --output {output.json} {params.branch_arg} > {log} 2>&1
        """


rule BUSTEDS_MH_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln,
        tree=rules.label_tree_TestForeground_partition.output.tree
    output:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "BUSTEDS-MH.json")
    params:
        branch_arg=lambda wc: _branch_arg(wc.sample)
    log:
        os.path.join(LOGDIR, "BUSTEDS_MH_partition.{sample}.part{part}.log")
    shell:
        r"""
        mkdir -p {OUTROOT}/{wildcards.sample}/selection/part{wildcards.part}
        {HYPHY} BUSTED \
          --alignment {input.aln} \
          --tree {input.tree} \
          --output {output.json} \
          --srv Yes \
          --multiple-hits Double+Triple \
          {params.branch_arg} \
          > {log} 2>&1
        """


rule RELAX_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln,
        tree=rules.label_tree_TestForeground_partition.output.tree,
        test_list=_test_list_for_sample
    output:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "RELAX.json")
    log:
        os.path.join(LOGDIR, "RELAX_partition.{sample}.part{part}.log")
    run:
        os.makedirs(os.path.dirname(output.json), exist_ok=True)
        logp = log[0] if hasattr(log, "__iter__") else str(log)
        with open(logp, "a") as handle:
            handle.write("START\n")
        if (not os.path.exists(input.test_list)) or (os.path.getsize(input.test_list) == 0):
            with open(output.json, "w") as fh:
                json.dump({"skipped": True, "reason": "empty_or_missing_test_list"}, fh, indent=2)
            with open(logp, "a") as handle:
                handle.write("Skipped (no test list)\nEND\n")
        else:
            shell(
                f"{HYPHY} RELAX --alignment {input.aln} --tree {input.tree} "
                f"--output {output.json} --test Test --reference Background "
                f"> {logp} 2>&1"
            )
            with open(logp, "a") as handle:
                handle.write("END\n")


rule CFEL_partition:
    input:
        aln=rules.sanitize_partition_alignment.output.aln,
        tree=rules.label_tree_TestForeground_partition.output.tree,
        test_list=_test_list_for_sample
    output:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "CFEL.json")
    log:
        os.path.join(LOGDIR, "CFEL_partition.{sample}.part{part}.log")
    run:
        os.makedirs(os.path.dirname(output.json), exist_ok=True)
        logp = log[0] if hasattr(log, "__iter__") else str(log)
        with open(logp, "a") as handle:
            handle.write("START\n")
        if (not os.path.exists(input.test_list)) or (os.path.getsize(input.test_list) == 0):
            with open(output.json, "w") as fh:
                json.dump({"skipped": True, "reason": "empty_or_missing_test_list"}, fh, indent=2)
            with open(logp, "a") as handle:
                handle.write("Skipped (no test list)\nEND\n")
        else:
            shell(
                f"{HYPHY} contrast-fel --alignment {input.aln} --tree {input.tree} "
                f"--output {output.json} --branch-set Test "
                f"> {logp} 2>&1"
            )
            with open(logp, "a") as handle:
                handle.write("END\n")


def _all_partition_jsons_for_sample(wc):
    parts = _parts_for_sample(wc)
    outs = []
    for part in parts:
        outs.extend([
            os.path.join(selection_part_dir(wc.sample, part), "FEL.json"),
            os.path.join(selection_part_dir(wc.sample, part), "MEME.json"),
            os.path.join(selection_part_dir(wc.sample, part), "ABSREL.json"),
            os.path.join(selection_part_dir(wc.sample, part), "BUSTEDS-MH.json"),
            os.path.join(selection_part_dir(wc.sample, part), "RELAX.json"),
            os.path.join(selection_part_dir(wc.sample, part), "CFEL.json"),
        ])
    return outs


rule run_selection_all_partitions:
    input:
        _all_partition_jsons_for_sample
    output:
        done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done")
    log:
        os.path.join(LOGDIR, "run_selection_all_partitions.{sample}.log")
    run:
        os.makedirs(os.path.dirname(output.done), exist_ok=True)
        logp = log[0] if hasattr(log, "__iter__") else str(log)
        with open(logp, "a") as handle:
            handle.write("START\n")
        with open(output.done, "w") as fh:
            fh.write("ok\n")
        with open(logp, "a") as handle:
            handle.write("END\n")
