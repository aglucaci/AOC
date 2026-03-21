rule PlotFEL_partition:
    input:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "FEL.json")
    output:
        output_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "part{part}", "{sample}.part{part}.AOC.FEL_Results.csv"),
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "part{part}", "{sample}.part{part}.FEL.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "part{part}", "{sample}.part{part}.FEL.svg"),
        output_figureLegend=os.path.join(str(OUTROOT), "{sample}", "visualizations", "part{part}", "{sample}.part{part}.FEL.figure_legend.txt")
    log:
        os.path.join(LOGDIR, "PlotFEL_partition.{sample}.part{part}.log")
    run:
        data = _read_json(input.json)
        headers = _hyphy_headers(data)
        df = pd.DataFrame(_hyphy_content_table(data), columns=headers, dtype=float)
        df.index += 1
        df["CodonSite"] = df.index
        df["adjusted_p-value"] = _safe_fdrcorrection(df["p-value"].tolist(), alpha=0.10)
        os.makedirs(os.path.dirname(output.output_csv), exist_ok=True)
        df.to_csv(output.output_csv, index=False)

        source = df.copy().dropna().rename(columns={"p-value": "p_value", "adjusted_p-value": "adjusted_p_value"})
        line = alt.Chart(source).mark_circle(clip=True, opacity=0.9, size=70).encode(
            x=alt.X("CodonSite:Q", title="Codon site"),
            y=alt.Y("dN/dS MLE:Q", title="dN/dS estimate", scale=alt.Scale(domain=(0, 5), clamp=True, nice=False, type="sqrt")),
            color=alt.condition("datum.adjusted_p_value <= 0.1", alt.value("red"), alt.value("lightgray")),
            tooltip=["CodonSite", "alpha", "beta", "dN/dS MLE", "p_value", "adjusted_p_value"]
        ).properties(width=900, height=450, title=f"FEL: {wildcards.sample} part {wildcards.part}")
        band = alt.Chart(source).mark_area(opacity=0.25).encode(x="CodonSite:Q", y="dN/dS LB:Q", y2="dN/dS UB:Q")
        save_altair_with_vl_convert(band + line, svg_path=output.output_svg, png_path=output.output_png, scale=2)

        sig = df[df["adjusted_p-value"] <= 0.10]
        pos = sig[sig["dN/dS MLE"] > 1.0]
        neg = sig[sig["dN/dS MLE"] < 1.0]
        with open(output.output_figureLegend, "w") as fh:
            fh.write(
                f"FEL identified {len(sig)} significant codon sites in {wildcards.sample} partition {wildcards.part} at FDR <= 0.10, including {len(pos)} sites with dN/dS > 1 and {len(neg)} sites with dN/dS < 1.\n"
            )


rule PlotMEME_partition:
    input:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "MEME.json")
    output:
        output_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "part{part}", "{sample}.part{part}.AOC.MEME_Results.csv"),
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "part{part}", "{sample}.part{part}.MEME.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "part{part}", "{sample}.part{part}.MEME.svg"),
        output_figureLegend=os.path.join(str(OUTROOT), "{sample}", "visualizations", "part{part}", "{sample}.part{part}.MEME.figure_legend.txt")
    log:
        os.path.join(LOGDIR, "PlotMEME_partition.{sample}.part{part}.log")
    run:
        data = _read_json(input.json)
        headers = _hyphy_headers(data)
        df = pd.DataFrame(_hyphy_content_table(data), columns=headers, dtype=float)
        df.index += 1
        df["CodonSite"] = df.index
        df["adjusted_p-value"] = _safe_fdrcorrection(df["p-value"].tolist(), alpha=0.10)
        os.makedirs(os.path.dirname(output.output_csv), exist_ok=True)
        df.to_csv(output.output_csv, index=False)

        source = df.copy().rename(columns={"p-value": "p_value", "adjusted_p-value": "adjusted_p_value", "beta+": "beta_plus"})
        tooltip_fields = [alt.Tooltip("CodonSite:Q", title="Codon site")]
        if "alpha" in source.columns:
            tooltip_fields.append(alt.Tooltip("alpha:Q", title="alpha"))
        if "beta_plus" in source.columns:
            tooltip_fields.append(alt.Tooltip("beta_plus:Q", title="beta+"))
        if "p_value" in source.columns:
            tooltip_fields.append(alt.Tooltip("p_value:Q", title="p-value"))
        tooltip_fields.append(alt.Tooltip("adjusted_p_value:Q", title="adjusted p-value"))

        chart = alt.Chart(source).mark_point(size=55).encode(
            x=alt.X("CodonSite:Q", title="Codon site"),
            y=alt.Y("adjusted_p_value:Q", title="Adjusted p-value"),
            color=alt.Color("adjusted_p_value:Q", scale=alt.Scale(scheme="reds", reverse=True)),
            tooltip=tooltip_fields
        ).properties(width=900, height=450, title=f"MEME: {wildcards.sample} part {wildcards.part}")
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)

        sig = df[df["adjusted_p-value"] <= 0.10]
        with open(output.output_figureLegend, "w") as fh:
            fh.write(
                f"MEME identified {len(sig)} significant codon sites in {wildcards.sample} partition {wildcards.part} at FDR <= 0.10, consistent with episodic selection along a subset of branches.\n"
            )


rule TableABSREL_partition:
    input:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "ABSREL.json")
    output:
        csv=os.path.join(str(OUTROOT), "{sample}", "tables", "part{part}", "{sample}.part{part}.AOC.ABSREL_Results.csv")
    log:
        os.path.join(LOGDIR, "TableABSREL_partition.{sample}.part{part}.log")
    run:
        data = _read_json(input.json)
        rows = []
        if data.get("skipped"):
            rows.append({"Sample": wildcards.sample, "Partition": int(wildcards.part), "Branch": None, "Corrected P-value": None, "Raw P-value": None, "Rate classes": None, "omega_max": None, "significant_branch_0.10": None, "status": "skipped"})
        else:
            branch_attrs = data.get("branch attributes", {}).get("0", {})
            for branch, vals in branch_attrs.items():
                if not isinstance(vals, dict):
                    continue
                cp = vals.get("Corrected P-value", None)
                rp = vals.get("P-value", None)
                try:
                    cp_num = float(cp)
                except (TypeError, ValueError):
                    cp_num = None
                omega_vals = [value for key, value in vals.items() if isinstance(key, str) and key.startswith("omega") and isinstance(value, (int, float))]
                rows.append({
                    "Sample": wildcards.sample,
                    "Partition": int(wildcards.part),
                    "Branch": branch,
                    "Corrected P-value": cp,
                    "Raw P-value": rp,
                    "Rate classes": len([key for key in vals.keys() if isinstance(key, str) and key.startswith("omega")]),
                    "omega_max": max(omega_vals) if omega_vals else None,
                    "significant_branch_0.10": int(cp_num <= 0.10) if cp_num is not None else 0,
                    "status": "ok",
                })
        pd.DataFrame(rows).to_csv(output.csv, index=False)


rule TableBUSTEDS_MH_partition:
    input:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "BUSTEDS-MH.json")
    output:
        csv=os.path.join(str(OUTROOT), "{sample}", "tables", "part{part}", "{sample}.part{part}.AOC.BUSTEDS-MH_Results.csv")
    log:
        os.path.join(LOGDIR, "TableBUSTEDS_MH_partition.{sample}.part{part}.log")
    run:
        data = _read_json(input.json)
        if data.get("skipped"):
            rows = [{"Sample": wildcards.sample, "Partition": int(wildcards.part), "test": "BUSTEDS-MH", "p_value": None, "LRT": None, "tested_branches": None, "background_branches": None, "status": "skipped"}]
        else:
            tr = data.get("test results", {})
            tested = data.get("tested", None)
            background = data.get("background", None)

            def count_branches(branches):
                if isinstance(branches, dict):
                    return sum(len(v) if isinstance(v, list) else 0 for v in branches.values())
                if isinstance(branches, list):
                    return len(branches)
                if isinstance(branches, int):
                    return branches
                return None

            rows = [{"Sample": wildcards.sample, "Partition": int(wildcards.part), "test": "BUSTEDS-MH", "p_value": tr.get("p-value", None), "LRT": tr.get("LRT", None), "tested_branches": count_branches(tested), "background_branches": count_branches(background), "status": "ok"}]
        pd.DataFrame(rows).to_csv(output.csv, index=False)


rule TableRELAX_partition:
    input:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "RELAX.json")
    output:
        csv=os.path.join(str(OUTROOT), "{sample}", "tables", "part{part}", "{sample}.part{part}.AOC.RELAX_Results.csv")
    log:
        os.path.join(LOGDIR, "TableRELAX_partition.{sample}.part{part}.log")
    run:
        data = _read_json(input.json)
        if data.get("skipped"):
            rows = [{"Sample": wildcards.sample, "Partition": int(wildcards.part), "test": "RELAX", "k": None, "p_value": None, "LRT": None, "selection_shift": None, "status": "skipped"}]
        else:
            tr = data.get("test results", {})
            k = tr.get("relaxation or intensification parameter", None)
            try:
                k_num = float(k)
            except (TypeError, ValueError):
                k_num = None
            shift = None if k_num is None else "intensified" if k_num > 1 else "relaxed" if k_num < 1 else "no_change"
            rows = [{"Sample": wildcards.sample, "Partition": int(wildcards.part), "test": "RELAX", "k": k, "p_value": tr.get("p-value", None), "LRT": tr.get("LRT", None), "selection_shift": shift, "status": "ok"}]
        pd.DataFrame(rows).to_csv(output.csv, index=False)


rule TableCFEL_partition:
    input:
        json=os.path.join(str(OUTROOT), "{sample}", "selection", "part{part}", "CFEL.json")
    output:
        csv=os.path.join(str(OUTROOT), "{sample}", "tables", "part{part}", "{sample}.part{part}.AOC.CFEL_Results.csv")
    log:
        os.path.join(LOGDIR, "TableCFEL_partition.{sample}.part{part}.log")
    run:
        data = _read_json(input.json)
        if data.get("skipped"):
            df = pd.DataFrame([{"Sample": wildcards.sample, "Partition": int(wildcards.part), "CodonSite": None, "status": "skipped"}])
        else:
            headers = _hyphy_headers(data)
            df = pd.DataFrame(_hyphy_content_table(data), columns=headers)
            df.index += 1
            df["CodonSite"] = df.index
            df.insert(0, "Partition", int(wildcards.part))
            df.insert(0, "Sample", wildcards.sample)
            df["status"] = "ok"
            pcol = next((c for c in df.columns if isinstance(c, str) and "p-value" in c.lower()), None)
            df["significant_site_0.10"] = (pd.to_numeric(df[pcol], errors="coerce") <= 0.10).astype(int) if pcol is not None else 0
        df.to_csv(output.csv, index=False)


rule MergeFELCSVs:
    input:
        all_done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"),
        fel_csvs=lambda wc: [fel_table(wc.sample, part) for part in _parts_for_sample(wc)]
    output:
        merged_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_FEL_Results.csv")
    log:
        os.path.join(LOGDIR, "MergeFELCSVs.{sample}.log")
    run:
        all_dfs = []
        for file in input.fel_csvs:
            df = pd.read_csv(file)
            part_match = re.search(r"part(\d+)", os.path.basename(file).replace(".AOC.FEL_Results.csv", ""))
            part = int(part_match.group(1)) if part_match else None
            df.insert(0, "Partition", part)
            df.insert(0, "Sample", wildcards.sample)
            all_dfs.append(df)
        merged_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame(columns=["Sample", "Partition"])
        merged_df.to_csv(output.merged_csv, index=False)


rule MergeMEMECSVs:
    input:
        all_done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"),
        meme_csvs=lambda wc: [meme_table(wc.sample, part) for part in _parts_for_sample(wc)]
    output:
        merged_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_MEME_Results.csv")
    log:
        os.path.join(LOGDIR, "MergeMEMECSVs.{sample}.log")
    run:
        all_dfs = []
        for file in input.meme_csvs:
            df = pd.read_csv(file)
            part_match = re.search(r"part(\d+)", os.path.basename(file).replace(".AOC.MEME_Results.csv", ""))
            part = int(part_match.group(1)) if part_match else None
            df.insert(0, "Partition", part)
            df.insert(0, "Sample", wildcards.sample)
            all_dfs.append(df)
        merged_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame(columns=["Sample", "Partition"])
        merged_df.to_csv(output.merged_csv, index=False)


rule MergeABSRELCSVs:
    input:
        all_done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"),
        csvs=lambda wc: [absrel_table(wc.sample, part) for part in _parts_for_sample(wc)]
    output:
        merged_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_ABSREL_Results.csv")
    log:
        os.path.join(LOGDIR, "MergeABSRELCSVs.{sample}.log")
    run:
        pd.concat([pd.read_csv(f) for f in input.csvs], ignore_index=True).to_csv(output.merged_csv, index=False)


rule MergeBUSTEDS_MH_CSVs:
    input:
        all_done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"),
        csvs=lambda wc: [busted_table(wc.sample, part) for part in _parts_for_sample(wc)]
    output:
        merged_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_BUSTEDS-MH_Results.csv")
    log:
        os.path.join(LOGDIR, "MergeBUSTEDS_MH_CSVs.{sample}.log")
    run:
        pd.concat([pd.read_csv(f) for f in input.csvs], ignore_index=True).to_csv(output.merged_csv, index=False)


rule MergeRELAXCSVs:
    input:
        all_done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"),
        csvs=lambda wc: [relax_table(wc.sample, part) for part in _parts_for_sample(wc)]
    output:
        merged_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_RELAX_Results.csv")
    log:
        os.path.join(LOGDIR, "MergeRELAXCSVs.{sample}.log")
    run:
        pd.concat([pd.read_csv(f) for f in input.csvs], ignore_index=True).to_csv(output.merged_csv, index=False)


rule MergeCFELCSVs:
    input:
        all_done=os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"),
        csvs=lambda wc: [cfel_table(wc.sample, part) for part in _parts_for_sample(wc)]
    output:
        merged_csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_CFEL_Results.csv")
    log:
        os.path.join(LOGDIR, "MergeCFELCSVs.{sample}.log")
    run:
        pd.concat([pd.read_csv(f) for f in input.csvs], ignore_index=True).to_csv(output.merged_csv, index=False)


rule PlotFELMerge:
    input:
        merged_csv=rules.MergeFELCSVs.output.merged_csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.FEL.merged.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.FEL.merged.svg"),
        marginal_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.FEL.marginals.png"),
        marginal_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.FEL.marginals.svg")
    log:
        os.path.join(LOGDIR, "PlotFELMerge.{sample}.log")
    run:
        df = pd.read_csv(input.merged_csv)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No FEL results"]})).mark_text(size=18).encode(text="message:N")
            marginal_chart = chart
        else:
            source = df.dropna(subset=["CodonSite", "dN/dS MLE"]).copy().rename(columns={"adjusted_p-value": "adjusted_p_value", "p-value": "p_value"})
            source["Partition"] = source["Partition"].astype(str)
            source["is_significant"] = source["adjusted_p_value"].le(0.10).fillna(False)
            source["selection_class"] = source["dN/dS MLE"].apply(
                lambda x: (
                    "Positive" if pd.notna(x) and x > 1.05
                    else "Purifying" if pd.notna(x) and x < 0.95
                    else "Neutral"
                )
            )
            source["site_density"] = source.groupby("Partition")["CodonSite"].transform("count")

            base = alt.Chart().encode(
                x=alt.X("CodonSite:Q", title="Codon site"),
                tooltip=["Partition", "CodonSite", "dN/dS MLE", "alpha", "beta", "p_value", "adjusted_p_value"]
            )
            band = base.mark_area(opacity=0.18, color="#8ecae6").encode(
                y=alt.Y("dN/dS LB:Q", title="dN/dS estimate", scale=alt.Scale(domain=(0, 6), clamp=True, nice=False, type="sqrt")),
                y2="dN/dS UB:Q"
            )
            dots = base.mark_circle(size=90, opacity=0.95, stroke="#ffffff", strokeWidth=0.6).encode(
                y=alt.Y("dN/dS MLE:Q", title="dN/dS estimate", scale=alt.Scale(domain=(0, 6), clamp=True, nice=False, type="sqrt")),
                color=alt.Color(
                    "selection_class:N",
                    title="Signal class",
                    scale=alt.Scale(
                        domain=["Positive", "Neutral", "Purifying"],
                        range=["#bc4749", "#8d99ae", "#3d5a80"]
                    )
                ),
                fillOpacity=alt.condition("datum.is_significant", alt.value(1.0), alt.value(0.28))
            )
            chart = (
                alt.layer(band, dots, data=source)
                .properties(width=680, height=420)
                .facet(column=alt.Column("Partition:O", title="Partition", header=alt.Header(labelFontSize=12, titleFontSize=14)))
                .properties(
                    title=alt.TitleParams(
                        text=f"Merged FEL landscape: {wildcards.sample}",
                        subtitle=[
                            "Point size and color emphasize site-level behavior across partitions.",
                            "Filled points pass FDR <= 0.10; dN/dS > 1.05 is shown as positive, dN/dS < 0.95 as purifying, and intermediate values as neutral."
                        ]
                    )
                )
            )

            partition_counts = source.groupby("Partition", as_index=False).agg(
                total_sites=("CodonSite", "count"),
                significant_sites=("is_significant", "sum")
            )
            partition_counts["background_sites"] = partition_counts["total_sites"] - partition_counts["significant_sites"]
            part_long = partition_counts.melt(
                id_vars=["Partition"],
                value_vars=["significant_sites", "background_sites"],
                var_name="category",
                value_name="count"
            )
            top_hist = alt.Chart(source).mark_bar(opacity=0.9).encode(
                x=alt.X("CodonSite:Q", bin=alt.Bin(maxbins=40), title="Codon site"),
                y=alt.Y("count():Q", title="Site density"),
                color=alt.Color("is_significant:N", title="FEL significance", scale=alt.Scale(domain=[True, False], range=["#bc4749", "#d9d9d9"]))
            ).properties(
                width=980,
                height=170,
                title="Marginal site density with significant sites highlighted"
            )
            right_bar = alt.Chart(part_long).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8).encode(
                x=alt.X("Partition:O", title="Partition"),
                y=alt.Y("count:Q", title="Sites"),
                color=alt.Color("category:N", title="", scale=alt.Scale(domain=["significant_sites", "background_sites"], range=["#bc4749", "#d9d9d9"])),
                tooltip=["Partition", "category", "count"]
            ).properties(
                width=980,
                height=220,
                title="Partition-level hit burden"
            )
            marginal_chart = alt.vconcat(top_hist, right_bar, spacing=24).properties(
                title=alt.TitleParams(
                    text=f"FEL marginal summaries: {wildcards.sample}",
                    subtitle=["Upper panel: where hits concentrate along the coding sequence.", "Lower panel: how strongly each partition contributes to the total FEL signal."]
                )
            )
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)
        save_altair_with_vl_convert(marginal_chart, svg_path=output.marginal_svg, png_path=output.marginal_png, scale=2)


rule PlotMEMEMerge:
    input:
        merged_csv=rules.MergeMEMECSVs.output.merged_csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.MEME.merged.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.MEME.merged.svg"),
        marginal_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.MEME.marginals.png"),
        marginal_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.MEME.marginals.svg")
    log:
        os.path.join(LOGDIR, "PlotMEMEMerge.{sample}.log")
    run:
        df = pd.read_csv(input.merged_csv)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No MEME results"]})).mark_text(size=18).encode(text="message:N")
            marginal_chart = chart
        else:
            source = df.dropna(subset=["CodonSite"]).copy().rename(columns={"adjusted_p-value": "adjusted_p_value", "p-value": "p_value", "beta+": "beta_plus"})
            source["Partition"] = source["Partition"].astype(str)
            source["episodic_hit"] = source["adjusted_p_value"].le(0.10).fillna(False)
            source["evidence"] = source["adjusted_p_value"].apply(lambda x: -math.log10(max(x, 1e-12)) if pd.notna(x) and x > 0 else 0.0)

            tooltip_fields = ["Partition", "CodonSite", "p_value", "adjusted_p_value"]
            if "beta_plus" in source.columns:
                tooltip_fields.append("beta_plus")
            if "alpha" in source.columns:
                tooltip_fields.append("alpha")

            base = alt.Chart().encode(
                x=alt.X("CodonSite:Q", title="Codon site"),
                tooltip=tooltip_fields
            )
            points = base.mark_circle(size=100, opacity=0.95, stroke="#fffaf0", strokeWidth=0.7).encode(
                y=alt.Y("evidence:Q", title="-log10(adjusted p-value)"),
                color=alt.Color("evidence:Q", title="Evidence", scale=alt.Scale(scheme="oranges")),
                fillOpacity=alt.condition("datum.episodic_hit", alt.value(1.0), alt.value(0.32))
            )
            chart = (
                alt.layer(points, data=source)
                .properties(width=680, height=420)
                .facet(column=alt.Column("Partition:O", title="Partition", header=alt.Header(labelFontSize=12, titleFontSize=14)))
                .properties(
                    title=alt.TitleParams(
                        text=f"Merged MEME landscape: {wildcards.sample}",
                        subtitle=[
                            "Higher points carry stronger episodic-selection evidence.",
                            "Filled points pass FDR <= 0.10 across the partitioned alignment."
                        ]
                    )
                )
            )

            part_hits = source.groupby("Partition", as_index=False).agg(
                episodic_hits=("episodic_hit", "sum"),
                mean_evidence=("evidence", "mean")
            )
            site_dist = alt.Chart(source).mark_bar(opacity=0.9).encode(
                x=alt.X("CodonSite:Q", bin=alt.Bin(maxbins=40), title="Codon site"),
                y=alt.Y("count():Q", title="Site density"),
                color=alt.Color("episodic_hit:N", title="MEME significance", scale=alt.Scale(domain=[True, False], range=["#d8572a", "#e9ecef"]))
            ).properties(
                width=980,
                height=170,
                title="Marginal site density with episodic hits highlighted"
            )
            burden = alt.Chart(part_hits).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8).encode(
                x=alt.X("Partition:O", title="Partition"),
                y=alt.Y("episodic_hits:Q", title="Significant sites"),
                color=alt.Color("mean_evidence:Q", title="Mean evidence", scale=alt.Scale(scheme="goldorange")),
                tooltip=["Partition", "episodic_hits", "mean_evidence"]
            ).properties(
                width=980,
                height=220,
                title="Partition-level episodic-selection burden"
            )
            marginal_chart = alt.vconcat(site_dist, burden, spacing=24).properties(
                title=alt.TitleParams(
                    text=f"MEME marginal summaries: {wildcards.sample}",
                    subtitle=["Upper panel: spatial clustering of candidate sites.", "Lower panel: partitions ranked by the number and strength of episodic hits."]
                )
            )
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)
        save_altair_with_vl_convert(marginal_chart, svg_path=output.marginal_svg, png_path=output.marginal_png, scale=2)


rule PlotABSRELMerge:
    input:
        merged_csv=rules.MergeABSRELCSVs.output.merged_csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.ABSREL.merged.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.ABSREL.merged.svg")
    log:
        os.path.join(LOGDIR, "PlotABSRELMerge.{sample}.log")
    run:
        df = pd.read_csv(input.merged_csv)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No ABSREL results"]})).mark_text(size=18).encode(text="message:N")
        else:
            df["significant_branch_0.10"] = pd.to_numeric(df.get("significant_branch_0.10", 0), errors="coerce").fillna(0)
            df["omega_max"] = pd.to_numeric(df.get("omega_max", 0), errors="coerce").fillna(0)
            summary = (
                df.groupby("Partition", as_index=False)
                .agg(significant_branches=("significant_branch_0.10", "sum"), max_omega=("omega_max", "max"))
            )
            summary["Partition"] = summary["Partition"].astype(str)
            bars = alt.Chart(summary).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8).encode(
                x=alt.X("Partition:O", title="Partition"),
                y=alt.Y("significant_branches:Q", title="Significant branches"),
                color=alt.Color("max_omega:Q", title="Max omega", scale=alt.Scale(scheme="goldred")),
                tooltip=["Partition", "significant_branches", "max_omega"]
            ).properties(width=900, height=420, title=f"aBSREL branch signals: {wildcards.sample}")
            points = alt.Chart(summary).mark_point(size=180, shape="diamond", color="#1c1c1c").encode(
                x="Partition:O",
                y="max_omega:Q",
                tooltip=["Partition", "significant_branches", "max_omega"]
            )
            chart = alt.layer(bars, points).resolve_scale(y="independent")
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


rule PlotBUSTEDSMHMerge:
    input:
        merged_csv=rules.MergeBUSTEDS_MH_CSVs.output.merged_csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.BUSTEDS-MH.merged.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.BUSTEDS-MH.merged.svg")
    log:
        os.path.join(LOGDIR, "PlotBUSTEDSMHMerge.{sample}.log")
    run:
        df = pd.read_csv(input.merged_csv)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No BUSTEDS-MH results"]})).mark_text(size=18).encode(text="message:N")
        else:
            df["p_value"] = pd.to_numeric(df.get("p_value"), errors="coerce")
            df["LRT"] = pd.to_numeric(df.get("LRT"), errors="coerce")
            df["score"] = df["p_value"].apply(lambda x: -math.log10(max(x, 1e-12)) if pd.notna(x) and x > 0 else 0.0)
            df["significant"] = df["p_value"].le(0.10).fillna(False)
            df["Partition"] = df["Partition"].astype(str)
            chart = alt.Chart(df).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8).encode(
                x=alt.X("Partition:O", title="Partition"),
                y=alt.Y("score:Q", title="-log10(p-value)"),
                color=alt.condition("datum.significant", alt.value("#b22222"), alt.value("#c9ced6")),
                tooltip=["Partition", "p_value", "LRT", "tested_branches", "background_branches"]
            ).properties(width=900, height=420, title=f"BUSTED-S-MH gene-wide evidence: {wildcards.sample}")
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


rule PlotRELAXMerge:
    input:
        merged_csv=rules.MergeRELAXCSVs.output.merged_csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.RELAX.merged.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.RELAX.merged.svg")
    log:
        os.path.join(LOGDIR, "PlotRELAXMerge.{sample}.log")
    run:
        df = pd.read_csv(input.merged_csv)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No RELAX results"]})).mark_text(size=18).encode(text="message:N")
        else:
            df["k"] = pd.to_numeric(df.get("k"), errors="coerce")
            df["p_value"] = pd.to_numeric(df.get("p_value"), errors="coerce")
            df["log2_k"] = df["k"].apply(lambda x: math.log(x, 2) if pd.notna(x) and x > 0 else 0.0)
            df["Partition"] = df["Partition"].astype(str)
            chart = alt.Chart(df).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8).encode(
                x=alt.X("Partition:O", title="Partition"),
                y=alt.Y("log2_k:Q", title="log2(k)"),
                color=alt.Color("selection_shift:N", title="Shift", scale=alt.Scale(domain=["intensified", "relaxed", "no_change"], range=["#b22222", "#1d4e89", "#9aa0a8"])),
                opacity=alt.condition("datum.p_value <= 0.1", alt.value(1.0), alt.value(0.5)),
                tooltip=["Partition", "k", "p_value", "LRT", "selection_shift"]
            ).properties(width=900, height=420, title=f"RELAX intensity shifts: {wildcards.sample}")
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


rule PlotCFELMerge:
    input:
        merged_csv=rules.MergeCFELCSVs.output.merged_csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.CFEL.merged.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.CFEL.merged.svg")
    log:
        os.path.join(LOGDIR, "PlotCFELMerge.{sample}.log")
    run:
        df = pd.read_csv(input.merged_csv)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No CFEL results"]})).mark_text(size=18).encode(text="message:N")
        else:
            df["Partition"] = df["Partition"].astype(str)
            if "significant_site_0.10" in df.columns:
                df["significant_site_0.10"] = pd.to_numeric(df["significant_site_0.10"], errors="coerce").fillna(0)
            else:
                pcol = next((c for c in df.columns if isinstance(c, str) and "p-value" in c.lower()), None)
                if pcol is not None:
                    df["significant_site_0.10"] = (
                        pd.to_numeric(df[pcol], errors="coerce").le(0.10).fillna(False).astype(int)
                    )
                else:
                    df["significant_site_0.10"] = 0

            if "CodonSite" in df.columns:
                df["CodonSite"] = pd.to_numeric(df["CodonSite"], errors="coerce")
                total_sites_col = "CodonSite"
            else:
                df["row_count"] = 1
                total_sites_col = "row_count"

            summary = (
                df.groupby("Partition", as_index=False)
                .agg(significant_sites=("significant_site_0.10", "sum"), total_sites=(total_sites_col, "count" if total_sites_col == "CodonSite" else "sum"))
            )
            summary["background_sites"] = summary["total_sites"] - summary["significant_sites"]
            source = summary.melt(id_vars=["Partition"], value_vars=["significant_sites", "background_sites"], var_name="category", value_name="count")
            chart = alt.Chart(source).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8).encode(
                x=alt.X("Partition:O", title="Partition"),
                y=alt.Y("count:Q", title="Sites"),
                color=alt.Color("category:N", title="", scale=alt.Scale(domain=["significant_sites", "background_sites"], range=["#0b6e4f", "#d7dee5"])),
                tooltip=["Partition", "category", "count"]
            ).properties(width=900, height=420, title=f"CFEL site contrast summary: {wildcards.sample}")
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


def _overview_inputs_for_sample(wc):
    outs = [
        os.path.join(str(OUTROOT), wc.sample, "tables", f"{wc.sample}.AOC.merged_FEL_Results.csv"),
        os.path.join(str(OUTROOT), wc.sample, "tables", f"{wc.sample}.AOC.merged_MEME_Results.csv"),
    ]
    for part in _parts_for_sample(wc):
        for test in ["ABSREL", "BUSTEDS-MH", "RELAX", "CFEL"]:
            outs.append(os.path.join(selection_part_dir(wc.sample, part), f"{test}.json"))
    return outs


rule SelectionOverviewTable:
    input:
        _overview_inputs_for_sample
    output:
        csv=os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.selection_overview.csv")
    log:
        os.path.join(LOGDIR, "SelectionOverviewTable.{sample}.log")
    run:
        rows = []
        fel = pd.read_csv(input[0]) if os.path.exists(input[0]) else pd.DataFrame()
        meme = pd.read_csv(input[1]) if os.path.exists(input[1]) else pd.DataFrame()
        if not fel.empty and "adjusted_p-value" in fel.columns and "Partition" in fel.columns:
            fel["adjusted_p-value"] = pd.to_numeric(fel["adjusted_p-value"], errors="coerce")
            for part, subdf in fel.groupby("Partition"):
                rows.append({"sample": wildcards.sample, "partition": part, "method": "FEL", "metric": "significant_sites_FDR_0.10", "value": int((subdf["adjusted_p-value"] <= 0.10).sum())})
        if not meme.empty and "adjusted_p-value" in meme.columns and "Partition" in meme.columns:
            meme["adjusted_p-value"] = pd.to_numeric(meme["adjusted_p-value"], errors="coerce")
            for part, subdf in meme.groupby("Partition"):
                rows.append({"sample": wildcards.sample, "partition": part, "method": "MEME", "metric": "significant_sites_FDR_0.10", "value": int((subdf["adjusted_p-value"] <= 0.10).sum())})
        for path in input[2:]:
            test = Path(path).stem
            data = _read_json(path)
            part_match = re.search(r"/part(\d+)/", str(path))
            partition = int(part_match.group(1)) if part_match else None
            if data.get("skipped"):
                rows.append({"sample": wildcards.sample, "partition": partition, "method": test, "metric": "status", "value": "skipped"})
                continue
            if test == "BUSTEDS-MH":
                rows.append({"sample": wildcards.sample, "partition": partition, "method": test, "metric": "gene_wide_pvalue", "value": data.get("test results", {}).get("p-value", None)})
            elif test == "RELAX":
                tr = data.get("test results", {})
                rows.append({"sample": wildcards.sample, "partition": partition, "method": test, "metric": "k", "value": tr.get("relaxation or intensification parameter", None)})
                rows.append({"sample": wildcards.sample, "partition": partition, "method": test, "metric": "pvalue", "value": tr.get("p-value", None)})
            elif test == "ABSREL":
                sig = 0
                for _, vals in data.get("branch attributes", {}).get("0", {}).items():
                    if not isinstance(vals, dict):
                        continue
                    try:
                        cp = float(vals.get("Corrected P-value", 1.0))
                    except (TypeError, ValueError):
                        cp = 1.0
                    if cp <= 0.10:
                        sig += 1
                rows.append({"sample": wildcards.sample, "partition": partition, "method": test, "metric": "significant_branches", "value": sig})
            elif test == "CFEL":
                try:
                    headers = [x[0] for x in data["MLE"]["headers"]]
                    df = pd.DataFrame(data["MLE"]["content"]["0"], columns=headers)
                    pcol = next((c for c in headers if "p-value" in c.lower()), None)
                    sig = int((pd.to_numeric(df[pcol], errors="coerce") <= 0.10).sum()) if pcol else 0
                except Exception:
                    sig = 0
                rows.append({"sample": wildcards.sample, "partition": partition, "method": test, "metric": "significant_sites_p_0.10", "value": sig})
        pd.DataFrame(rows, columns=["sample", "partition", "method", "metric", "value"]).to_csv(output.csv, index=False)


rule PlotSelectionMethodHeatmap:
    input:
        overview=rules.SelectionOverviewTable.output.csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.method_heatmap.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.method_heatmap.svg")
    log:
        os.path.join(LOGDIR, "PlotSelectionMethodHeatmap.{sample}.log")
    run:
        df = pd.read_csv(input.overview)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No integrated selection summary"]})).mark_text(size=18).encode(text="message:N")
        else:
            metric_map = {
                "significant_sites_FDR_0.10": "Site hits",
                "significant_sites_p_0.10": "Site hits",
                "significant_branches": "Branch hits",
                "gene_wide_pvalue": "Gene p-value",
                "pvalue": "Gene p-value",
                "k": "RELAX k",
                "status": "Status",
            }
            df["display_metric"] = df["metric"].map(metric_map).fillna(df["metric"])
            numeric = pd.to_numeric(df["value"], errors="coerce")
            df["score"] = numeric.fillna(0)
            mask = df["metric"].isin(["gene_wide_pvalue", "pvalue"])
            df.loc[mask, "score"] = numeric[mask].apply(lambda x: -math.log10(max(x, 1e-12)) if pd.notna(x) and x > 0 else 0.0)
            df["interpretation"] = df["metric"].map({
                "significant_sites_FDR_0.10": "Count of significant codon sites after FDR correction",
                "significant_sites_p_0.10": "Count of significant codon sites at p <= 0.10",
                "significant_branches": "Number of branches with evidence for episodic selection",
                "gene_wide_pvalue": "Gene-wide signal converted to -log10(p)",
                "pvalue": "RELAX significance converted to -log10(p)",
                "k": "RELAX intensity parameter; >1 intensified, <1 relaxed",
                "status": "Skipped analysis"
            }).fillna("Method-specific summary metric")
            df["value_label"] = df["value"].astype(str)
            df["partition"] = df["partition"].astype(str)
            heat = alt.Chart(df).mark_rect(cornerRadius=6).encode(
                x=alt.X("partition:O", title="Partition"),
                y=alt.Y("method:N", title="Selection test"),
                color=alt.Color("score:Q", title="Evidence score", scale=alt.Scale(scheme="teals")),
                tooltip=[
                    alt.Tooltip("partition:N", title="Partition"),
                    alt.Tooltip("method:N", title="Method"),
                    alt.Tooltip("metric:N", title="Metric"),
                    alt.Tooltip("value_label:N", title="Raw value"),
                    alt.Tooltip("score:Q", title="Heatmap score"),
                    alt.Tooltip("interpretation:N", title="Interpretation")
                ]
            )
            labels = alt.Chart(df).mark_text(fontSize=12, fontWeight="bold", color="#0f172a").encode(
                x="partition:O",
                y="method:N",
                text=alt.Text("value:N")
            )
            explainer = alt.Chart(pd.DataFrame([
                {"line": "How to read this heatmap", "order": 1},
                {"line": "Columns are partitions; rows are HyPhy methods.", "order": 2},
                {"line": "Text inside each cell is the raw value reported by that method.", "order": 3},
                {"line": "Cell color is a standardized evidence score so site counts, branch counts, RELAX k, and p-values can be scanned together.", "order": 4},
                {"line": "For p-values the color score is -log10(p), so darker cells mean smaller p-values.", "order": 5},
                {"line": "For FEL, MEME, aBSREL, and CFEL, darker cells mean more significant sites or branches in that partition.", "order": 6},
                {"line": "For RELAX k, values above 1 suggest intensified selection and values below 1 suggest relaxation; color reflects distance from 0 on the shared score axis, not direction.", "order": 7},
                {"line": "Use repeated dark cells in the same column as the strongest clue that a partition is consistently interesting across tests.", "order": 8}
            ])).mark_text(align="left", baseline="top", fontSize=14).encode(
                y=alt.Y("order:O", axis=None),
                text="line:N"
            ).properties(width=980, height=188)
            chart = alt.vconcat(
                explainer,
                alt.layer(heat, labels).properties(width=980, height=360)
            ).properties(
                title=alt.TitleParams(
                    text=f"Integrated selection evidence matrix: {wildcards.sample}",
                    subtitle=["Read across rows to compare methods within a gene partition.", "Read down columns to find partitions repeatedly highlighted by different tests."]
                )
            )
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


rule PlotSelectionConcordance:
    input:
        overview=rules.SelectionOverviewTable.output.csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.concordance.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.concordance.svg")
    log:
        os.path.join(LOGDIR, "PlotSelectionConcordance.{sample}.log")
    run:
        df = pd.read_csv(input.overview)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No concordance data"]})).mark_text(size=18).encode(text="message:N")
        else:
            df["partition"] = df["partition"].astype(str)
            numeric = pd.to_numeric(df["value"], errors="coerce")
            df["score"] = numeric.fillna(0)
            pval_mask = df["metric"].isin(["gene_wide_pvalue", "pvalue"])
            df.loc[pval_mask, "score"] = numeric[pval_mask].apply(lambda x: -math.log10(max(x, 1e-12)) if pd.notna(x) and x > 0 else 0.0)
            evidence_mask = (
                ((df["metric"].isin(["significant_sites_FDR_0.10", "significant_sites_p_0.10", "significant_branches"])) & (numeric.fillna(0) > 0)) |
                ((df["metric"].isin(["gene_wide_pvalue", "pvalue"])) & (numeric.fillna(1) <= 0.10)) |
                ((df["metric"] == "k") & (numeric.fillna(1).sub(1.0).abs() >= 0.1))
            )
            df["supports_partition"] = evidence_mask.astype(int)
            summary = df.groupby("partition", as_index=False).agg(
                methods_supporting=("supports_partition", "sum"),
                cumulative_score=("score", "sum")
            )
            summary["support_band"] = summary["methods_supporting"].apply(
                lambda x: "broad support" if x >= 4 else "moderate support" if x >= 2 else "limited support"
            )
            chart = alt.Chart(summary).mark_circle(opacity=0.95, stroke="#ffffff", strokeWidth=1.0).encode(
                x=alt.X("partition:O", title="Partition"),
                y=alt.Y("methods_supporting:Q", title="Methods with nonzero / significant signal"),
                size=alt.Size("cumulative_score:Q", title="Cumulative evidence", scale=alt.Scale(range=[150, 2200])),
                color=alt.Color("support_band:N", title="Concordance", scale=alt.Scale(domain=["limited support", "moderate support", "broad support"], range=["#c9ced6", "#ee9b00", "#9b2226"])),
                tooltip=["partition", "methods_supporting", "cumulative_score", "support_band"]
            ).properties(
                width=980,
                height=380,
                title=alt.TitleParams(
                    text=f"Selection concordance by partition: {wildcards.sample}",
                    subtitle=[
                        "Higher points indicate that more tests independently flagged the partition.",
                        "Larger circles indicate stronger total evidence accumulated across all methods."
                    ]
                )
            )
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


rule PlotSelectionDashboard:
    input:
        overview=rules.SelectionOverviewTable.output.csv
    output:
        output_png=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.dashboard.png"),
        output_svg=os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.dashboard.svg")
    log:
        os.path.join(LOGDIR, "PlotSelectionDashboard.{sample}.log")
    run:
        df = pd.read_csv(input.overview)
        if df.empty:
            chart = alt.Chart(pd.DataFrame({"message": ["No dashboard data"]})).mark_text(size=18).encode(text="message:N")
        else:
            numeric = pd.to_numeric(df["value"], errors="coerce").fillna(0)
            df["score"] = numeric
            pval_mask = df["metric"].isin(["gene_wide_pvalue", "pvalue"])
            df.loc[pval_mask, "score"] = numeric[pval_mask].apply(lambda x: -math.log10(max(x, 1e-12)) if pd.notna(x) and x > 0 else 0.0)

            method_summary = df.groupby("method", as_index=False)["score"].sum()
            partition_summary = df.groupby("partition", as_index=False)["score"].sum()
            partition_summary["partition"] = partition_summary["partition"].astype(str)
            top_partition = partition_summary.sort_values("score", ascending=False).head(1)
            top_method = method_summary.sort_values("score", ascending=False).head(1)

            left = alt.Chart(method_summary).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8, color="#9b2226").encode(
                x=alt.X("method:N", title="Method"),
                y=alt.Y("score:Q", title="Cumulative evidence"),
                tooltip=["method", "score"]
            ).properties(width=470, height=360, title="Method burden")

            right = alt.Chart(partition_summary).mark_bar(cornerRadiusTopLeft=8, cornerRadiusTopRight=8, color="#005f73").encode(
                x=alt.X("partition:O", title="Partition"),
                y=alt.Y("score:Q", title="Cumulative evidence"),
                tooltip=["partition", "score"]
            ).properties(width=470, height=360, title="Partition burden")

            summary_lines = []
            if not top_method.empty:
                summary_lines.append(f"Top method by cumulative evidence: {top_method.iloc[0]['method']} ({top_method.iloc[0]['score']:.2f})")
            if not top_partition.empty:
                summary_lines.append(f"Top partition by cumulative evidence: {top_partition.iloc[0]['partition']} ({top_partition.iloc[0]['score']:.2f})")
            summary_lines.append("Cumulative evidence pools raw counts, RELAX k values, and transformed p-values into one comparative score.")
            summary_lines.append("Use this dashboard to find where multiple tests point to the same partition-level story.")
            explain_df = pd.DataFrame({"line": summary_lines, "order": list(range(1, len(summary_lines) + 1))})
            explainer = alt.Chart(explain_df).mark_text(align="left", baseline="top", fontSize=14).encode(
                y=alt.Y("order:O", axis=None),
                text="line:N"
            ).properties(width=980, height=24 * len(summary_lines) + 8)

            chart = alt.vconcat(
                explainer,
                alt.hconcat(left, right, spacing=28).resolve_scale(y="shared")
            ).properties(
                title=alt.TitleParams(
                    text=f"Selection dashboard: {wildcards.sample}",
                    subtitle=["Left: which tests contribute the strongest aggregate signal.", "Right: which partitions accumulate the strongest cross-method evidence."]
                )
            )
        save_altair_with_vl_convert(chart, svg_path=output.output_svg, png_path=output.output_png, scale=2)


rule executive_summary:
    input:
        fel=rules.MergeFELCSVs.output.merged_csv,
        meme=rules.MergeMEMECSVs.output.merged_csv,
        fel_plot=rules.PlotFELMerge.output.output_png,
        meme_plot=rules.PlotMEMEMerge.output.output_png,
        absrel_plot=rules.PlotABSRELMerge.output.output_png,
        busted_plot=rules.PlotBUSTEDSMHMerge.output.output_png,
        relax_plot=rules.PlotRELAXMerge.output.output_png,
        cfel_plot=rules.PlotCFELMerge.output.output_png,
        heatmap=rules.PlotSelectionMethodHeatmap.output.output_png,
        concordance=rules.PlotSelectionConcordance.output.output_png,
        dashboard=rules.PlotSelectionDashboard.output.output_png,
        overview=rules.SelectionOverviewTable.output.csv
    output:
        html=os.path.join(str(OUTROOT), "{sample}", "summary", "executive_summary.html")
    log:
        os.path.join(LOGDIR, "executive_summary.{sample}.log")
    run:
        os.makedirs(os.path.dirname(output.html), exist_ok=True)
        logp = log[0] if hasattr(log, "__iter__") else str(log)
        if os.path.exists(EXEC_SUMMARY_PY):
            shell(f"python {EXEC_SUMMARY_PY} --fel {input.fel} --meme {input.meme} --out {output.html} > {logp} 2>&1")
        else:
            fel = pd.read_csv(input.fel) if os.path.exists(input.fel) else pd.DataFrame()
            meme = pd.read_csv(input.meme) if os.path.exists(input.meme) else pd.DataFrame()
            overview = pd.read_csv(input.overview) if os.path.exists(input.overview) else pd.DataFrame()
            fel_sig = int((fel.get("adjusted_p-value", pd.Series(dtype=float)) <= 0.10).sum()) if not fel.empty else 0
            meme_sig = int((meme.get("adjusted_p-value", pd.Series(dtype=float)) <= 0.10).sum()) if not meme.empty else 0
            overview_html = overview.to_html(index=False) if not overview.empty else "<p>No overview metrics available.</p>"
            with open(output.html, "w") as fh:
                fh.write(f"""<html><body>
                <h1>Executive Summary: {wildcards.sample}</h1>
                <p>FEL significant sites (FDR <= 0.10): {fel_sig}</p>
                <p>MEME significant sites (FDR <= 0.10): {meme_sig}</p>
                <h2>Selection overview</h2>
                {overview_html}
                <h2>Plots</h2>
                <p><img src="../visualizations/{wildcards.sample}.FEL.merged.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.FEL.marginals.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.MEME.merged.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.MEME.marginals.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.ABSREL.merged.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.BUSTEDS-MH.merged.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.RELAX.merged.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.CFEL.merged.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.selection.method_heatmap.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.selection.concordance.png" width="1200"></p>
                <p><img src="../visualizations/{wildcards.sample}.selection.dashboard.png" width="1200"></p>
                </body></html>""")
