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
 - name: Department of Physiology and Biophysics, Weill Cornell Medicine, Cornell University, New York, NY 10021, USA
   index: 1
 - name: Institute for Genomics and Evolutionary Medicine, Temple University, Philadelphia, PA, USA
   index: 2
date: 01 July 2025
bibliography: paper.bib

---

# Summary

Modern molecular sequence analysis increasingly relies on automated and robust software tools for interpretation, annotation, and biological insight. The Analysis of Orthologous Collections (AOC) application automates the identification of genomic sites and species/lineages influenced by natural selection in coding sequence analysis. AOC quantifies different types of selection: negative, diversifying or directional positive, or differential selection between groups of branches. We include all steps necessary to go from unaligned homologous sequences to complete results and interactive visualizations that are designed to aid in the useful interpretation and contextualization. We are motivated by a desire to make evolutionary analyses as simple as possible, and to close the disparity in the literature between genes which draw a significant amount of interest and those that are largely overlooked and underexplored. We believe that such underappreciated and understudied genetic datasets can hold rich biological information and offer substantial insights into the diverse patterns and processes of evolution, especially if domain experts are able to perform the analyses themselves. 

# 1 Introduction

Genomic research is inevitably biased towards certain organisms (humans, model organisms, agriculturally important species, pathogens), and genes (biomedically important, functionally understood) [@Stoeger2018ignored]. For example, GeneRif -- a database of the reference set of articles describing the function of a gene [@GeneRIF2023, last accessed July 6, 2023], is dominated by 5 species: Humans, Mouse, Rat, Arabidopsis, Drosophila corresponding to about 92% of total coverage; Humans alone represent 63% of all GeneRifs. A highly skewed coverage of protein functional information concentrated in a largely anthropocentric fashion fails to benefit from the potential knowledge gained from studying the diversity of the natural world.
The AOC application is designed to be a one-stop shop for molecular sequence evaluation using state of the art methods and techniques. The pipeline is fully automated and incorporates recombination detection, a powerful force in shaping gene evolution which can produce spurious results if not considered. The application is simple to install and use, requiring few dependencies and few input files or configuration. We differentiate ourselves from other approaches in the field [@Picard2020dgin] by data preparation steps we take (see Figure 1), and the selection analysis modalities we take advantage of which include lineage-specific and site-level information, and search for pervasive or episodic selective patterns with consideration of positive, negative, directional, biochemical, between-group comparison, and relaxed evolutionary forces [@Lucaci2022bdnf].
 We are also motivated by the so-called “day science” and “night science” [@Yanai2019nightscience] scientific duality. Here, “day science” is the application and evaluation of a priori hypotheses which are validated or falsified by the available data. We apply this kind of evaluation because each of the selection analysis methods we use are designed to ask and answer biological and statistical questions (we highlight these in the Implementation section). However, we also focus on “night science” where a user can explore the “unstructured realm of possible hypotheses, of ideas not yet fully fleshed out” [@Yanai2019nightscience] which may not have occurred to the user when they first set out to evaluate their gene of interest. Therefore, AOC is designed as a blend between the two philosophical lines of inquiry, where a user can approach the application with a particular hypothesis in mind, but also allows for data exploration to serve as a guide on a scientific adventure not previously considered. In addition, as the AOC application use grows, the results of experiments can become part of a kind of genetic profile, allowing for placement in a repository and subsequent meta-analysis.

# 2 Methods 
## 2.1 Implementation
The application is designed for use with the NCBI Gene database via [www.ncbi.nlm.nih.gov/gene](www.ncbi.nlm.nih.gov/gene) and retrieve gene orthologs. This can be done based on a single sequence per species, which is recommended if multiple transcripts are available, to limit data bias. Depending on study design we may also limit our search to only include species specific taxonomic groups (birds, turtles, lizards, mammals, etc). These queries return full gene transcript (RefSeq transcript) and protein sequence (RefSeq protein) files with tabular data (CSV-format) containing useful metadata (including NCBI accession numbers). Other sources of genomic information can also be used. We use protein sequences and full gene transcripts to derive coding sequences (CDS) via a custom script: scripts/codons.py. We also recommend using only high-quality protein sequences, as “PREDICTED” or “PARTIAL” sequence files may contain errors and are not appropriate for downstream selection analysis. Our application removes low-quality protein sequences from downstream analysis, as they may inflate rates of nonsynonymous change or otherwise bias the analyses. 
The AOC application is designed for comprehensive protein-coding molecular sequence analysis. AOC allows for the inclusion of recombination detection, which is a powerful force in shaping gene evolution and critically important to correctly interpreting analytic results which are vulnerable to changing recombinant topologies. We also include an automated method for lineage assignment and annotation which relies on input tabular data (e.g. from NCBI Gene) and NCBI Taxonomy information. Lineage assignment allows for between-group comparisons of selective pressures using selection analysis. 
The application accepts two input data files from the NCBI Orthologs database: a protein sequence unaligned FASTA file, and a transcript sequence unaligned FASTA file for the same gene. Typically, this can be retrieved from public databases such as NCBI Gene (described above). Although this is the recommended route, other methods of data compilation are also acceptable. If protein sequence and transcript sequence files are provided, a custom script scripts/codons.py is executed and returns a CDS FASTA file. Note that the application is easily modifiable to accept a single CDS input, if such data are available to the user. This script is currently set to assume the standard genetic code, this can be modified for alternate codon tables. This script also removes low-quality sequences including those where no match is found.

## 2.2 Pre-processing
To generate multiple sequence alignments, we use MACSEv2 [@Ranwez2018macse] due to its ability to create codon-aware multiple sequence alignment. We also measure the Tamura-Nei 1993 (TN93) genetic distance of alignments using the HyPhy implmentation of [TN93](https://github.com/veg/tn93). Recombination detection is automatically performed using Genetic Algorithm for Recombination Detection (GARD) [@KosakovskyPond2006]. A recombination-free set of alignment fragments is placed in the results folder where phylogenetic tree inference and downstream selection analysis are performed. For datasets where recombination is not detected this results in a single file for analysis. In datasets where recombination is detected, we parse out recombinant partitions into multiple files correcting for recombinant breakpoints which occur within a codon. Next, phylogenetic tree inference is done for all the recombination-free FASTA files, we perform maximum-likelihood (ML) phylogenetic inference via IQ-TREE [@Minh2020iqtree2].
For all the unrooted phylogenetic trees via an automated lineage annotation script that uses the NCBI and the python package ete3 toolkit [@HuertaCepas2016]. Lineages are binned into taxonomic groups. Here, the aim is to have a broad representation of taxonomic groups, rather than the species being heavily clustered into a single group. We perform tree labeling via the hyphy-analyses script Label-Trees method and results in one annotated tree with a designation for all lineages [(HyPhy-analyses): Label Trees](https://github.com/veg/hyphy-analyses/tree/master/LabelTrees).


## 2.3 Selection analysis
All recombination-free alignment and unrooted phylogenetic tree is evaluated through a suite of molecular evolutionary methods designed to ask and answer specific biological and statistical questions including (Table 1) [@Spielman2019evolution, @KosakovskyPond2020].  

**Table 1. Summary of selection analysis methods**

| Method              | Description |
|---------------------|-------------|
| **FEL**             | Locates codon sites with evidence of pervasive positive diversifying or negative selection. Answers: Which site(s) in a gene are subject to pervasive diversifying selection? [@KosakovskyPond2005] |
| **BUSTED[+S+MH]**   | Tests for gene-wide episodic selection while accounting for synonymous rate variation and multiple instantaneous substitutions. [@Wisotsky2020synonymous; @Lucaci2023shortcuts] |
| **MEME**            | Detects codon sites under episodic positive diversifying selection. Answers: Which site(s) are subject to episodic or pervasive diversifying selection? [@Murrell2012meme] |
| **aBSREL**          | Tests if positive selection has occurred on a proportion of branches. [Smith2015absrel] |
| **SLAC**            | Performs substitution mapping to detect pervasive diversifying selection. [@KosakovskyPond2005] |
| **BGM**             | Identifies groups of sites that are co-evolving. [@Poon2008spidermonkey] |
| **RELAX**           | Compares gene-wide selection pressure between a query clade and background lineages to detect relaxation/intensification. [@Wertheim2015relax] |
| **Contrast-FEL**    | Compares site-by-site selection pressure between query and background sequences. [@ContrastFEL2021] |
| **FitMultiModel**   | Tests model fit by allowing multiple instantaneous substitutions. [@Lucaci2021extrabase] |
| **FUBAR**           | Identifies sites under pervasive selection using a fast Bayesian approach. [@Murrell2013fubar] |

## 2.4 Visualizations and Tables
We provide a high-level executive summary and multiple-test correction of the selection analyses and on input files where available for information such as sequence divergence. In addition, we generate figures from all selection analyses along with accompanying summary result tables and figure legends which describe the results. Individual results, specifically output JSON files from HyPhy analyses may also be visualized using [Hyphy-Vision](http://vision.hyphy.org) or interactive ObservableHQ [@Perkel2021notebooks] notebooks [HyPhy: Interactive Observable Notebooks](https://observablehq.com/@hyphy).

![Flowchart diagram of the AOC workflow and an example using Primate ACE2 data. The workflow consists of three parts, the first of which does quality control, and converts input transcript and protein files from the NCBI ortholog database into codon-aware alignments and checks for phylogenetic evidence of genetic recombination. The second part performs full maximum-likelihood phylogenetic inference and lineage annotation based on NCBI Taxonomy and runs a full suite of selection detection methods using HyPhy. The last part consists of summarizing results into useful tables and visualizations that can be used for post-hoc interpretation and interactions.](figures/AOC-Fig1.png)


## 2.5 Testing and benchmarking
As an example, using an application of AOC, we were able to report on novel sites of adaptive evolution, broad relationships of coevolution, and independently verify previously reported results on the signatures of purifying selection in the mammalian BDNF [@Lucaci2022bdnf] gene, which plays a critical role in brain development. 
We also explored the evolutionary history of the primate ACE2 protein. Data was accessed from NCBI via the Ortholog database. We downloaded FASTA files from 32 species, with RefSeq Transcripts and RefSeq Proteins (one sequence per species) and metadata in tabular form (CSV). Additional details of our analysis, including all intermediate and HyPhy JSON files are available in our GitHub repository.
 For more information on how selection analysis scales along with dataset complexity and size, we refer the reader to HyPhy benchmarking results available at [HyPhy: Benchmarks and Profiling](https://observablehq.com/@stevenweaver/hyphy-benchmarks-and-profiling.).

# 3 Conclusion
The application of modern pipelines for molecular sequence evaluation is of critical importance. These methods have proven to be powerful [@Martin2022omicron; @Viana2022omicron; @Tegally2022ba45; @Martin2021n501y; @Silva2023utrigen; @Zehr2023feline] to detect the role of natural selection in shaping proteins and offer the ability to further interrogate their results with carefully designed experimental approaches. The combination of computational and experimental biology has the potential to drive significant innovation and discovery in both the basic and translational sciences. AOC is designed to play a role in scientific and medical discovery by providing a simple-to-use software application for molecular sequence analysis especially for insights into unexplored genetic datasets.

# Acknowledgements
We would like to thank members of the [HyPhy](http://lab.hyphy.org/) and [Datamonkey](https://www.datamonkey.org/) teams for their contributions to this project, method development, and the maintenance of state-of-the-art molecular sequence analysis software. This work was supported by a NIH grant (GM151683) to SLKP.

# References