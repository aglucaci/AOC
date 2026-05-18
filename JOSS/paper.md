---
title: 'AOC: A Snakemake workflow for the characterization of natural selection in protein-coding genes'
tags:
  - Python
  - Snakemake
  - Molecular Evolution
  - Bioinformatics
authors:
  - name: Alexander G. Lucaci
    orcid: 0000-0002-4896-6088
    corresponding: true
    affiliation: 1
  - name: Sergei Pond
    orcid: 0000-0003-4817-4029
    corresponding: true
    affiliation: 2
affiliations:
 - name: Department of Systems and Computational Biomedicine, Weill Cornell Medicine, Cornell University, New York, NY 10021, United States of America
   index: 1
 - name: Institute for Genomics and Evolutionary Medicine, Temple University, Philadelphia, PA, United States of America
   index: 2
date: 22 April 2026
bibliography: paper.bib

---

# Summary

Modern molecular sequence analysis increasingly relies on automated and robust software tools for interpretation, annotation, and biological insight. The Analysis of Orthologous Collections (AOC) application automates the identification of genomic sites and species/lineages influenced by natural selection in coding sequence analysis. AOC quantifies different types of selection: negative, diversifying or directional positive, or differential selection between groups of branches. We include all steps necessary to go from unaligned homologous sequences to complete results and interactive visualizations that are designed to aid in the interpretation and contextualization of inferred selection signals (e.g., site-level dN/dS estimates, lineage-specific selection patterns, and statistical support), enabling users to relate these results to functional, evolutionary, and biological hypotheses. We are motivated by a desire to make evolutionary analyses as simple as possible, and to close the disparity in the literature between genes which draw a significant amount of interest and those that are largely overlooked and underexplored. We believe that such underappreciated and understudied genetic datasets can hold rich biological information and offer substantial insights into the diverse patterns and processes of evolution, especially if domain experts are able to perform the analyses themselves. 

# 1 Introduction

Genomic research is inevitably biased towards certain organisms (humans, model organisms, agriculturally important species, pathogens), and genes (biomedically important, functionally understood) [@Stoeger2018ignored]. For example, GeneRif -- a database of the reference set of articles describing the function of a gene [@GeneRIF2023, last accessed July 6, 2023], is dominated by 5 species: Humans, Mouse, Rat, Arabidopsis, Drosophila corresponding to about 92% of total coverage; Humans alone represent 63% of all GeneRifs. A highly skewed coverage of protein functional information concentrated in a largely anthropocentric fashion fails to benefit from the potential knowledge gained from studying the diversity of the natural world.
The Analysis of Orthologous Collections (AOC) application is designed to be a one-stop shop for molecular sequence evaluation using state of the art methods and techniques. The pipeline is fully automated and incorporates recombination detection, a powerful force in shaping gene evolution which can produce spurious results if not considered. The application is simple to install and use, requiring few dependencies and few input files or configuration. We differentiate ourselves from other approaches in the field [@Picard2020dgin] through both comprehensive data preparation (Figure 1) and the breadth of selection analyses performed. Specifically, AOC integrates lineage-specific and site-level inference to detect pervasive and episodic selection, including negative and positive (diversifying and directional) selection, shifts in amino acid preferences, differential selection between predefined groups of branches, and changes in selection intensity (relaxation or intensification) [@Lucaci2022bdnf].

# 2 Methods 
## 2.1 Implementation
The AOC application is designed for comprehensive protein-coding molecular sequence analysis. AOC (a Snakemake workflow), allows for the inclusion of recombination detection, which is a powerful force in shaping gene evolution and critically important to correctly interpreting analytic results which are vulnerable to changing recombinant topologies. Lineage assignment allows for between-group comparisons of selective pressures using selection analysis. The application accepts CDS FASTA files in which each file corresponds to a single gene (i.e., a set of homologous sequences), typically retrieved from public databases such as NCBI Gene or curated by the user. AOC supports both single-gene analyses (one file) and multi-gene analyses (multiple files) within the same workflow.

## 2.2 Pre-processing
To generate multiple sequence alignments, we use MACSEv2 [@Ranwez2018macse] due to its ability to create codon-aware multiple sequence alignment. We also measure the Tamura-Nei 1993 (TN93) genetic distance of alignments using the HyPhy implmentation of [TN93](https://github.com/veg/tn93). Recombination detection is automatically performed using Genetic Algorithm for Recombination Detection (GARD) [@KosakovskyPond2006]. A recombination-free set of alignment fragments is placed in the results folder where phylogenetic tree inference and downstream selection analysis are performed. For datasets where recombination is not detected this results in a single file for analysis. In datasets where recombination is detected, we parse out recombinant partitions into multiple files correcting for recombinant breakpoints which occur within a codon. Next, phylogenetic tree inference is done for all the recombination-free FASTA files, we perform phylogenetic inference via FastTree2 [@price_fasttree_2010].
 We perform tree labeling via the hyphy-analyses script Label-Trees method and results in one annotated tree with a designation for all lineages:[(HyPhy-analyses): Label Trees](https://github.com/veg/hyphy-analyses/tree/master/LabelTrees).


## 2.3 Selection analysis
All recombination-free alignments and unrooted phylogenetic trees are evaluated using a suite of molecular evolutionary methods, each designed to address specific biological and statistical questions (described in Table 1; [@Spielman2019evolution; @KosakovskyPond2020]).

**Table 1. Summary of selection analysis methods**

| Method              | Description |
|---------------------|-------------|
| **FEL**             | Locates codon sites with evidence of pervasive positive diversifying or negative selection. Answers: Which site(s) in a gene are subject to pervasive diversifying selection? [@KosakovskyPond2005] |
| **BUSTED[+S+MH]**   | Tests for gene-wide episodic selection while accounting for synonymous rate variation and multiple instantaneous substitutions. [@Wisotsky2020synonymous; @Lucaci2023shortcuts] |
| **MEME**            | Detects codon sites under episodic positive diversifying selection. Answers: Which site(s) are subject to episodic or pervasive diversifying selection? [@Murrell2012meme] |
| **aBSREL**          | Tests if positive selection has occurred on a proportion of branches. [@Smith2015absrel] |
| **SLAC**            | Performs substitution mapping to detect pervasive diversifying selection. [@KosakovskyPond2005] |
| **BGM**             | Identifies groups of sites that are co-evolving. [@Poon2008spidermonkey] |
| **RELAX**           | Compares gene-wide selection pressure between a query clade and background lineages to detect relaxation/intensification. [@Wertheim2015relax] |
| **Contrast-FEL**    | Compares site-by-site selection pressure between query and background sequences. [@ContrastFEL2021] |
| **FitMultiModel**   | Tests model fit by allowing multiple instantaneous substitutions. [@Lucaci2021extrabase] |
| **FUBAR**           | Identifies sites under pervasive selection using a fast Bayesian approach. [@Murrell2013fubar] |

## 2.4 Visualizations and Tables
We provide a high-level executive summary and multiple-test correction of the selection analyses and on input files where available for information such as sequence divergence. In addition, we generate figures from all selection analyses along with accompanying summary result tables and figure legends which describe the results. Individual results, specifically output JSON files from HyPhy analyses may also be visualized using [Hyphy-Vision](http://vision.hyphy.org) or interactive ObservableHQ [@Perkel2021notebooks] notebooks [HyPhy: Interactive Observable Notebooks](https://observablehq.com/@hyphy).

![Flowchart diagram of the AOC (Snakemake) workflow and an example using Primate ACE2 data. The workflow consists of three parts, the first of which does quality control, and converts input transcript and protein files from the NCBI ortholog database into codon-aware alignments and checks for phylogenetic evidence of genetic recombination. The second part performs full maximum-likelihood phylogenetic inference and lineage annotation based on NCBI Taxonomy and runs a full suite of selection detection methods using HyPhy. The last part consists of summarizing results into useful tables and visualizations that can be used for post-hoc interpretation and interactions.](figures/AOC-Fig1.png)


## 2.5 Testing and benchmarking
As an example, using an application of AOC, we were able to report on novel sites of adaptive evolution, broad relationships of coevolution, and independently verify previously reported results on the signatures of purifying selection in the mammalian BDNF [@Lucaci2022bdnf] gene, which plays a critical role in brain development. 
We also explored the evolutionary history of the primate ACE2 protein. Data was accessed from NCBI via the Ortholog database. We downloaded FASTA files from 32 species, with RefSeq Transcripts and RefSeq Proteins (one sequence per species) and metadata in tabular form (CSV). Additional details of our analysis, including all intermediate and HyPhy JSON files are available in our GitHub repository.
 For more information on how selection analysis scales along with dataset complexity and size, we refer the reader to HyPhy benchmarking results available at [HyPhy: Benchmarks and Profiling](https://observablehq.com/@stevenweaver/hyphy-benchmarks-and-profiling.).

# 3 State of the field
Codon-based models of molecular evolution are widely implemented in established tools such as PAML[@yang_paml_2007] and HyPhy[@KosakovskyPond2020], which provide statistically rigorous frameworks for detecting selection through site, branch, and branch-site models. While these platforms offer powerful inference engines, they are primarily designed for model execution rather than standardized, large-scale, and reproducible workflow orchestration. Existing wrappers and graphical interfaces, such as the HyPhy graphical user interface (HYPHY Vision), Datamonkey, PAML front-ends, and general workflow platforms like Galaxy, lower the barrier to running individual analyses but do not integrate alignment quality control, phylogenetic reconstruction, coordinated execution of multiple complementary selection tests, structured result aggregation, and publication-ready visualization within a unified, reproducible framework. AOC addresses this gap by building on established inference engines (particularly HyPhy), and embedding them within a modular Snakemake-based system that formalizes best practices for codon-aware evolutionary analysis. Rather than introducing new statistical methodology, AOC contributes a reproducible infrastructure that enables scalable, consistent, and comparative application of state-of-the-art evolutionary models across genes and ortholog collections, addressing a critical workflow-level bottleneck in the field.

# 4 Statement of need
Comparative evolutionary analyses of protein-coding genes often require coordinating multiple complex steps: alignment quality control, phylogenetic reconstruction, branch labeling, and complementary selection tests across many genes and datasets. While powerful tools such as HyPhy exist, running these analyses reproducibly at scale typically requires substantial scripting, manual intervention, and fragmented result handling. AOC addresses this gap by providing a unified, reproducible workflow that integrates codon-aware alignment processing, tree inference, structured branch labeling, multiple selection analyses, and standardized result aggregation within a single framework. The target audience includes evolutionary biologists, comparative genomicists, molecular evolution researchers, and computational biologists who require scalable, transparent, and publication-ready evolutionary inference across large gene collections or multi-species datasets.

# 5 Software design
AOC is implemented using Snakemake to ensure reproducibility, scalability, and transparent dependency management across multi-step evolutionary analyses. The workflow adopts a modular, per-gene architecture in which each CDS FASTA file represents a single gene (a set of homologous sequences), enabling straightforward parallelization and fault isolation while supporting both single- and multi-gene analyses. Rather than reimplementing statistical models, AOC integrates established tools in the HyPhy suite (Table 1), prioritizing methodological robustness and community validation over less reliable custom implementations. These design choices allow users to move from raw sequence data to interpretable selection inferences in a reproducible and extensible manner, particularly for large or underexplored gene sets.

# 6 Research impact statement
AOC is designed to support reproducible and scalable molecular evolutionary analyses and has been validated through application to diverse coding sequence datasets, with fully reproducible workflows and example datasets provided to demonstrate end-to-end functionality. The software integrates established methods in HyPhy within a unified pipeline, lowering the barrier to performing complex selection analyses across genes and lineages. Early indicators of community readiness include its use in internal and collaborative research projects, and interest from external users seeking reproducible evolutionary analysis workflows [@Lucaci2022bdnf; @Zehr2023feline; @Martin2022omicron; @Silva2023utrigen]. By enabling standardized, multi-gene selection analyses with minimal setup, AOC addresses a clear need for accessible, extensible pipelines in comparative genomics and evolutionary biology.

# 7 Conclusion
Modern molecular sequence analysis pipelines enable the detection of natural selection and generation of testable biological hypotheses [@Martin2022omicron; @Viana2022omicron; @Tegally2022ba45; @Martin2021n501y; @Silva2023utrigen; @Zehr2023feline]. The AOC workflow is designed to play a role in scientific and medical discovery by providing a simple-to-use software application for molecular sequence analysis, especially for insights into unexplored genetic datasets.

# Acknowledgements
We would like to thank members of the [HyPhy](http://lab.hyphy.org/) and [Datamonkey](https://www.datamonkey.org/) teams for their contributions to this project, method development, and the maintenance of state-of-the-art molecular sequence analysis software. This work was supported by a NIH grant (GM151683) to SLKP.

# AI usage disclosure
Generative AI tools were used to assist with aspects of code development and manuscript preparation. All outputs were critically evaluated, tested, and validated by the authors to ensure correctness and reproducibility. The authors retain full responsibility for the software design, implementation, and scientific conclusions.

# References
