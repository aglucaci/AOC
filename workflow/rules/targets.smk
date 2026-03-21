def _default_recombination_targets():
    return [
        os.path.join(str(OUTROOT), sample, f"{sample}.RD.SA.codons.cln.fa.cluster.fasta.best-gard")
        for sample in SAMPLES
        if _recombination_enabled(sample)
    ]


rule all:
    default_target: True
    input:
        expand(os.path.join(str(OUTROOT), "{sample}", "{sample}.RD.SA.codons.cln.fa"), sample=SAMPLES),
        _default_recombination_targets(),
        expand(os.path.join(str(OUTROOT), "{sample}", "selection", "ALL_PARTITIONS.done"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "{sample}.test_sequences.txt"), sample=SAMPLES_WITH_LABELS),
        expand(os.path.join(str(OUTROOT), "{sample}", "{sample}.sequence_header_map.csv"), sample=SAMPLES_WITH_LABELS),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_FEL_Results.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_MEME_Results.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_ABSREL_Results.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_BUSTEDS-MH_Results.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_RELAX_Results.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.AOC.merged_CFEL_Results.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "tables", "{sample}.selection_overview.csv"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.FEL.merged.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.FEL.marginals.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.MEME.merged.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.MEME.marginals.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.ABSREL.merged.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.BUSTEDS-MH.merged.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.RELAX.merged.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.CFEL.merged.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.method_heatmap.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.concordance.png"), sample=SAMPLES),
        expand(os.path.join(str(OUTROOT), "{sample}", "visualizations", "{sample}.selection.dashboard.png"), sample=SAMPLES)


ruleorder:
    macse >
    cln >
    strike_ambigs_msa >
    remove_duplicates_msa >
    tn93_cluster >
    recombination >
    build_header_map_and_test_list >
    parse_gard >
    sanitize_partition_alignment >
    iqtree_partition >
    label_tree_TestForeground_partition >
    FEL_partition >
    MEME_partition >
    ABSREL_partition >
    BUSTEDS_MH_partition >
    RELAX_partition >
    CFEL_partition >
    run_selection_all_partitions >
    PlotFEL_partition >
    PlotMEME_partition >
    TableABSREL_partition >
    TableBUSTEDS_MH_partition >
    TableRELAX_partition >
    TableCFEL_partition >
    MergeFELCSVs >
    MergeMEMECSVs >
    MergeABSRELCSVs >
    MergeBUSTEDS_MH_CSVs >
    MergeRELAXCSVs >
    MergeCFELCSVs >
    PlotFELMerge >
    PlotMEMEMerge >
    PlotABSRELMerge >
    PlotBUSTEDSMHMerge >
    PlotRELAXMerge >
    PlotCFELMerge >
    SelectionOverviewTable >
    PlotSelectionMethodHeatmap >
    PlotSelectionConcordance >
    PlotSelectionDashboard >
    executive_summary
