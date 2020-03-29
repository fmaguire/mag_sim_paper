---
author-meta:
- Finlay Maguire*
- Baofeng Jia*
- Kristen Gray
- Wing Yin Venus Lau
- Robert G. Beiko
- Fiona S.L. Brinkman
bibliography:
- content/manual-references.json
date-meta: '2020-03-29'
header-includes: '<!--

  Manubot generated metadata rendered from header-includes-template.html.

  Suggest improvements at https://github.com/manubot/manubot/blob/master/manubot/process/header-includes-template.html

  -->

  <meta name="dc.format" content="text/html" />

  <meta name="dc.title" content="Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids and Genomic Islands" />

  <meta name="citation_title" content="Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids and Genomic Islands" />

  <meta property="og:title" content="Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids and Genomic Islands" />

  <meta property="twitter:title" content="Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids and Genomic Islands" />

  <meta name="dc.date" content="2020-03-29" />

  <meta name="citation_publication_date" content="2020-03-29" />

  <meta name="dc.language" content="en-CA" />

  <meta name="citation_language" content="en-CA" />

  <meta name="dc.relation.ispartof" content="Manubot" />

  <meta name="dc.publisher" content="Manubot" />

  <meta name="citation_journal_title" content="Manubot" />

  <meta name="citation_technical_report_institution" content="Manubot" />

  <meta name="citation_author" content="Finlay Maguire*" />

  <meta name="citation_author_institution" content="Faculty of Computer Science, Dalhousie University" />

  <meta name="citation_author_orcid" content="0000-0002-1203-9514" />

  <meta name="twitter:creator" content="@fmaguire" />

  <meta name="citation_author" content="Baofeng Jia*" />

  <meta name="citation_author_institution" content="Department of Molecular Biology and Biochemistry, Simon Fraser University" />

  <meta name="citation_author" content="Kristen Gray" />

  <meta name="citation_author_institution" content="Department of Molecular Biology and Biochemistry, Simon Fraser University" />

  <meta name="citation_author_orcid" content="0000-0002-1962-189X" />

  <meta name="citation_author" content="Wing Yin Venus Lau" />

  <meta name="citation_author_institution" content="Department of Molecular Biology and Biochemistry, Simon Fraser University" />

  <meta name="citation_author" content="Robert G. Beiko" />

  <meta name="citation_author_institution" content="Faculty of Computer Science, Dalhousie University" />

  <meta name="citation_author" content="Fiona S.L. Brinkman" />

  <meta name="citation_author_institution" content="Department of Molecular Biology and Biochemistry, Simon Fraser University" />

  <link rel="canonical" href="https://fmaguire.github.io/mag_sim_paper/" />

  <meta property="og:url" content="https://fmaguire.github.io/mag_sim_paper/" />

  <meta property="twitter:url" content="https://fmaguire.github.io/mag_sim_paper/" />

  <meta name="citation_fulltext_html_url" content="https://fmaguire.github.io/mag_sim_paper/" />

  <meta name="citation_pdf_url" content="https://fmaguire.github.io/mag_sim_paper/manuscript.pdf" />

  <link rel="alternate" type="application/pdf" href="https://fmaguire.github.io/mag_sim_paper/manuscript.pdf" />

  <link rel="alternate" type="text/html" href="https://fmaguire.github.io/mag_sim_paper/v/e760c7d5b40b64f1d18715c7dcc93c4a1e142a95/" />

  <meta name="manubot_html_url_versioned" content="https://fmaguire.github.io/mag_sim_paper/v/e760c7d5b40b64f1d18715c7dcc93c4a1e142a95/" />

  <meta name="manubot_pdf_url_versioned" content="https://fmaguire.github.io/mag_sim_paper/v/e760c7d5b40b64f1d18715c7dcc93c4a1e142a95/manuscript.pdf" />

  <meta property="og:type" content="article" />

  <meta property="twitter:card" content="summary_large_image" />

  <link rel="icon" type="image/png" sizes="192x192" href="https://manubot.org/favicon-192x192.png" />

  <link rel="mask-icon" href="https://manubot.org/safari-pinned-tab.svg" color="#ad1457" />

  <meta name="theme-color" content="#ad1457" />

  <!-- end Manubot generated metadata -->'
keywords:
- markdown
- publishing
- manubot
lang: en-CA
manubot-clear-requests-cache: false
manubot-output-bibliography: output/references.json
manubot-output-citekeys: output/citations.tsv
manubot-requests-cache-path: ci/cache/requests-cache
title: Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids and Genomic Islands
...






<small><em>
This manuscript
([permalink](https://fmaguire.github.io/mag_sim_paper/v/e760c7d5b40b64f1d18715c7dcc93c4a1e142a95/))
was automatically generated
from [fmaguire/mag_sim_paper@e760c7d](https://github.com/fmaguire/mag_sim_paper/tree/e760c7d5b40b64f1d18715c7dcc93c4a1e142a95)
on March 29, 2020.
</em></small>

## Authors



+ **Finlay Maguire***<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0002-1203-9514](https://orcid.org/0000-0002-1203-9514)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [fmaguire](https://github.com/fmaguire)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [fmaguire](https://twitter.com/fmaguire)<br>
  <small>
     Faculty of Computer Science, Dalhousie University
     · Funded by ['Genome Canada', 'Donald Hill Family Fellowship']
  </small>

+ **Baofeng Jia***<br>
    · ![GitHub icon](images/github.svg){.inline_icon}
    [imasianxd](https://github.com/imasianxd)<br>
  <small>
     Department of Molecular Biology and Biochemistry, Simon Fraser University
  </small>

+ **Kristen Gray**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0002-1962-189X](https://orcid.org/0000-0002-1962-189X)<br>
  <small>
     Department of Molecular Biology and Biochemistry, Simon Fraser University
  </small>

+ **Wing Yin Venus Lau**<br><br>
  <small>
     Department of Molecular Biology and Biochemistry, Simon Fraser University
  </small>

+ **Robert G. Beiko**<br><br>
  <small>
     Faculty of Computer Science, Dalhousie University
  </small>

+ **Fiona S.L. Brinkman**<br><br>
  <small>
     Department of Molecular Biology and Biochemistry, Simon Fraser University
  </small>



## Abstract {#abstract}

# Motivation

Metagenomic methods emerges as a key tool in public-health microbiology for surveillance of virulence and antimicrobial resistance (AMR) genes. 
Metagenomic data, even when assembled, results in a complex, mixed set of DNA fragments rather than nicely resolved individual genomes. 
Recently, metagenome-assembled genomes (MAGs) emerged as approaches that attempt to group these DNA fragments into bins that are likely derived from the same underlying genome. 
The most important sequences for surveillance purposes are those associated with mobile genetic elements (MGEs) such as plasmids and genomic islands (GIs). 
we hypothesise that due to the different copy number and sequence composition of plasmids and GIs compared to their corresponding chromosomes, these types of sequences will be highly under-represented in MAG-based approaches.

# Results

To evaluate the impact of MAG recovery methods on recovery of AMR genes and MGEs, we generated a simulated metagenomic dataset comprised of 30 genomes with up to 16.65% of the chrosomomal DNA consisting of GIs and 65 associated plasmids. 
MAGs were then recovered from this data using 12 different MAG pipelines and evaluated for recovery accuracies. 
Across all pipelines, 81.9-94.3% of chromosomes were recovered and binned. However, only 37.8-44.1% of GIs and 1.5-29.2% of plasmids were recovered and correctly binned at >50% coverage. 
In terms of AMR and virulence factor (VF) genes associated with MGEs, 0-45% of GI-associated AMR genes and 0-16% of GI-associated VF genes were correctly assigned. 
More strikingly, 0% of plasmid-borne VF or AMR genes were recovered.
This work shows that regardless of the MAG recovery approach used, plasmid and GI dominated sequences will disproportionately be left unbinned or incorrectly binned. 
From a public-health perspective, this means MAG approaches are unsuited for analysis of mobile genes, especially vital groups such as AMR and VF genes. 
This underlines the utility of read-based and long-read approaches to thoroughly evaluate the resistome in metagenomic data.


## Introduction {#intro}

Metagenomics, the sequencing of fragments of DNA from within an environmental sample, is widely used for characterising viral and microbial communities [@doi:10.1073/pnas.202488399; @doi:10.1038/nbt.3935]. 
By randomly sampling from the total genomic content these methods allow researchers to simultaneously profile the functional potential and the taxonomic identity of a large proportion of the organisms in a sample. 
Metagenomic techniques are now bein used to profile antimicrobial resistance (AMR) and pathogen virulence. 
These approaches have been instrumental in developing our understanding of the distribution and evolutionary history of AMR genes [@doi:10.1016/j.cell.2014.08.032; @doi:10.1016/j.mib.2007.08.009; @doi:10.1038/nature10388], as well as tracking pathogen outbreaks [@doi:10.1001/jama.2013.3231].

While long-read DNA sequencing technology (e.g., Oxford Nanopore [@doi:10.1111/1755-0998.12324], PacBio [@doi:10.1126/science.1162986]) is now being used for metagenomic sequencing [@doi:10.1093/gigascience/giz043; @doi:10.1186/s12866-019-1500-0],  high-throughput sequencing of relatively short reads (150-250bp) in platforms such as the Illumina MiSeq currently dominate metagenomic analyses.
Inference of taxonomic and functional diversity can be assessed directly from sequenced reads using reference databases and BLAST-based sequence similarity search tools (e.g. DIAMOND [@doi:10.1038/nmeth.3176]), read mapping (e.g. Bowtie 2 [@doi:10.1038/nmeth.1923]), Hidden Markov Models (e.g. HMMER3 [@doi:10.1093/bioinformatics/btt403]) or k-mer hashing (e.g. CLARK [@doi:10.1186/s12864-015-1419-2]). 
These read-based approaches allow analysis of all reads with detectable similarity to genes of interest even if the organism has relatively low abundance in the sample.
Since these reads are shorter than most genes, however, read-based methods provide very little information about the genomic organisation of genes.
This lack of contextual information is particularly problematic in the study of AMR genes and virulence factors as the genomic context plays a role in function [@doi:10.1128/AAC.01710-09], selective pressures [@doi:10.1016/j.tim.2006.02.006], and then likelihood of lateral gene transfer (LGT) [@doi:10.1111/j.1574-6976.2011.00273.x].

Sequence assembly is often used to genereate information about genomic context [@doi:10.1038/nature02340]. 
de Bruijn graph-based assemblers have been developed to handle the particular challenges of this type of data including metaSPAdes [@doi:10.1101/gr.213959.116] , IDBA-UD [@doi:10.1093/bioinformatics/bts174], and megahit [@doi:10.1093/bioinformatics/btv033].
A crucial challenge in metagenomic analysis is that reads from different origanisms must be distentangled to avoid hybrid assemblies.
A common way to deal with this challenge is to assign all contigs from a given source genomes to a cluster or "bin" based on similarities in the relative abudance and sequence composition.
These resulting bins are often known as metagenome-assembled genomes (MAGs). 
This binning is typically performed by grouping all the contigs with similar abundance and similar sequence composition into the same bin. 
A range of tools have been released to perform this binning including CONCOCT [@doi:10.1093/bioinformatics/btw290], MetaBAT 2 [@doi:10.7287/peerj.preprints.27522v1], and MaxBin 2 [@doi:10.1093/bioinformatics/btv638]. 
There is also the meta-binning tool DAS Tool [@doi:10.1038/s41564-018-0171-1] which combines predictions from multiple binning tools together. 
These MAG approaches have been used to great effect in unveiling huge amounts of previously uncharacterised genomic diversity [@doi:10.1038/nature14486; @doi:10.1038/s41564-017-0012-7; @doi:10.1101/489443].

There is loss of information at both the metagenomic assembly and binning steps.
This compounded data loss means that only a relatively small proportion of reads are successfully assembled and binned in large complex metagenome datasets, for example, 24.2-36.4% of reads from permafrost [@doi:10.1038/s41586-018-0338-1] and soil metagenomes [@doi:10.1038/s41564-019-0449-y]. 
Additionally, a large number of detected genomes are not reconstructed at all with ~23% of all detected genomes recovered in some examples [@doi:10.1038/s41564-019-0449-y]. 
The Critical Assessment of Metagenome Interpretation (CAMI) challenge's (https://data.cami-challenge.org/) Assessment of Metagenome BinnERs (AMBER) [@doi:10.1093/gigascience/giy069] assesses the global completeness and purity of recovered MAGs across methods.
However, to our best knowledge, there hasn't been a specific assessment of the impact of metagenomic assembly and binning on the loss of specific genomic elements.
In particular, the impact on mobile genetic elements (MGEs), such as genomic islands (GIs) and plasmids, which can be of great health and research importance, has not been evaluated.

Genomic islands (GIs) are clusters of genes that are known or predicted to have been acquired through LGT events. 
GIs can arise following the integation of MGEs, such as integrons, transposons, integrative and conjugative elements (ICEs) and prophages (integrated phages) [@doi:10.1038/nrmicro2350; @doi:10.1038/nrg3962]; they disproportionately encode virulence factors [@doi:10.1371/journal.pone.0008094] and are a major mechanism of LGT of AMR genes [@doi:10.3389/fmicb.2016.00173; @doi:10.1016/j.plasmid.2015.01.001]. 
GIs often have different nucleotide composition compared to the rest of the genome [@doi:10.1038/nrmicro2350].
This compositional difference is exploited by tools designed to detect GIs such as SIGI-HMM [@doi:10.1186/1471-2105-5-22] and IslandPath-DIMOB [@doi:10.1093/bioinformatics/bty095]. 
GIs may exist as multiple copies within a genome [@doi:10.1093/bib/bby042] leading to potential difficulties in correctly assembling these regions in metagenome assemblies as well as likely biases in the calculation of coverage statistics. 

Plasmids are circular or linear extrachromosomal self-replicating pieces of DNA. 
Similar to GIs, plasmids's sequence composition are markley different compared to the genome they are associated with [@doi:10.1093/bioinformatics/btq299; @doi:10.1093/molbev/msp281]. 
This is largely attributable to their repetitive sequences, variable copy number, and different selection pressures [@doi:10.1038/s41559-016-0010; @doi:10.1128/AAC.00235-15].
Plasmids are of great research importance, these elements are a major source of the lateral dissemination of AMR genes throughout microbial ecosystems [@pmid:26603922; @doi:10.3389/fmicb.2016.00173]. 
Due to these reasons, the correct assembly of DNA sequence of plasmid origin has proven to be difficult from short read data [@doi:10.1099/mgen.0.000128,@doi:10.1099/mgen.0.000206].

GIs and plasmids pose significant challenges in MAG recovery due to their unusual sequence compositon and relative abundance; as these MGEs are key to the function and spread of pathogenic traits such as AMR and virulence, it is vital that we assess the impact of metagenome assembly and binning on the representation of these specific elements. 
This is particularly important with the increasing popularity of MAG approaches within microbial and public-health research. 
Therefore, to address this issue we performed an analysis of GI and plasmid recovery accuracy across a set of state-of-the-art short-read metagenome assembly and binning approaches using a simulated metagenome comprised of GI- and plasmid-rich taxa.


## Materials and Methods {#methods}

All analyses presented in this paper can be reproduced and inspected with the associated github repository [github.com/fmaguire/MAG_gi_plasmid_analysis](github.com/fmaguire/MAG_gi_plasmid_analysis) and data repository [osf.io/nrejs/](osf.io/nrejs/).

### Metagenome Simulation

All genomes were selected from the set of completed RefSeq genomes as of April 2019.
Genomic islands for these genomes were previously predicted using IslandPath-DIMOB [@doi:10.1093/bioinformatics/bty095] and collated into the IslandViewer database [www.pathogenomics.sfu.ca/islandviewer](www.pathogenomics.sfu.ca/islandviewer) [@doi:10.1093/nar/gkv401].
Plasmid sequences were recovered for each genome using the linked GenBank Project IDs.
Thirty genomes were manually selected to satisfy the following criteria: 

	1) 10 genomes with 1-10 plasmids.

	2) 10 genomes with >10% of chromosomal DNA predicted to reside in GIs.

	3) 10 genomes with <1% of chromosomal DNA predicted to reside in GIs.

The data used to select the taxa is listed in Supplemental Table 1 and the details of the selected subset taxa are listed in Supplemental Table 2 with their NCBI accessions.
The sequences themselves are available in the data repository [osf.io/nrejs/](osf.io/nrejs/) under "data/sequences".

In accordance with the recommendation in the CAMI challenge [@doi:10.1038/nmeth.4458] the genomes were randomly assigned a relative abundance following a log-normal distribution (μ = 1, σ = 2).
Plasmid copy number estimates could not be accurately found for all organisms, therefore, plasmids were randomly assigned a copy number regime: low (1-20), medium (20-100), or high (500-1000) at a 2:1:1 rate.
Within each regime, the exact copy number was selected using an appropriately scaled gamma distribution (α = 4, β = 1) truncated to the regime range.

Finally, the effective plasmid relative abundance was determined by multiplying the plasmid copy number with the genome relative abundance.
The full set of randomly assigned relative abundances and copy numbers can be found in Supplemental Table 3.
Sequences were then concatenated into a single FASTA file with the appropriate relative abundance.
MiSeq v3 250bp paired-end reads with a mean fragment length of 1000bp (standard deviation of 50bp) were then simulated using art_illumina (v2016.06.05) [@doi:10.1093/bioinformatics/btr708] resulting in a simulated metagenome of 31,174,411 read pairs.
The selection of relative abundance and metagenome simulation itself was performed using the "data_simluation/simulate_metagenome.py" script.

### Metagenome Assembled Genome Recovery

Reads were trimmed using sickle (v1.33) [@url:sickle] resulting in 25,682,644 surviving read pairs.
The trimmed reads were then assembled using 3 different metagenomic assemblers: metaSPAdes (v3.13.0) [@doi:10.1101/gr.213959.116], IDBA-UD (v1.1.3) [@doi:10.1093/bioinformatics/bts174], and megahit (v1.1.3) [@doi:10.1093/bioinformatics/btv033]).
The resulting assemblies were summarised using metaQUAST (v5.0.2) [@doi:10.1093/bioinformatics/btv697].
The assemblies were then indexed and reads mapped back using Bowtie 2 (v2.3.4.3) [@doi:10.1038/nmeth.1923].

Samtools (v1.9) was used to sort the read mappings and the read coverage calculated using the MetaBAT2 accessory script (jgi\_summarize\_bam\_contig\_depths).
The three metagenome assemblies were then separately binned using MetaBAT2 (v2.13) [@doi:10.7287/peerj.preprints.27522v1], and MaxBin 2 (v2.2.6) [@doi:10.1093/bioinformatics/btv638]. 
MAGs were also recovered using CONCOCT (v0.4.2) [@doi:10.1093/bioinformatics/btw290] following the recommended protocol in the user manual.
Briefly, the supplied CONCOCT accessory scripts were used to cut contigs into 10 kilobase fragments (cut\_up\_fasta.py) and read coverage calculated for the fragments (CONCOCT\_coverage\_table.py).
These fragment coverages were then used to bin the 10kb fragments before the clustered fragments were merged (merge\_cutup\_clustering.py) to create the final CONCOCT MAG bins (extra\_fasta\_bins.py).
Finally, for each metagenome assembly the predicted bins from these three binners (Maxbin2, MetaBAT 2, and CONCOCT) were combined using the DAS Tool (v1.1.1) meta-binner [@doi:10.1038/s41564-018-0171-1].
This resulted in 12 separate sets of MAGs (one set for each assembler and binner pair).

### MAG assessment

#### Chromosomal Coverage

The MAG assessment for chromosomal coverage was performed by creating a BLASTN 2.9.0+ [@doi:10.1186/1471-2105-10-421] database consisting of all the chromosomes of the input reference genomes.
Each MAG contig was then used as a query against this database and the coverage of the underlying chromosomes tallied by merging the overlapping aligning regions and summing the total length of aligned MAG contigs.
The most represented genome in each MAG was assigned as the “identity” of that MAG for further analyses.
Coverages less than 5% were filtered out and the number of different genomes that contigs from a given MAG aligned to were tallied.
Finally, the overall proportion of chromosomes that were not present in any MAG was tallied for each binner and assembler.

In order to investigate the impact of close relatives in the metagenome on ability to bin chromosomes we generated a phylogenetic tree for all the input genomes.
Specifically, single copy universal bacterial proteins were identified in the reference genomes using BUSCO v4.0.2 with the Bacteria Odb10 data [@doi:10.1093/bioinformatics/btv351].
The 86 of these proteins that were found in every reference genome were concatenated and aligned using MAFFT v7.427 [@doi:10.1093/bioinformatics/bty121] and masked with trimal v1.4.1-3 [@doi:10.1093/bioinformatics/btp348].
A maximum-likelihood phylogeny was then inferred with IQ-Tree v1.6.12 [@doi:10.1093/molbev/msu300] with the in-built ModelFinder determined partitioning [@doi:10.1093/molbev/mss020].
Pairwise branch distances were then extracted from the resulting tree using ETE3 v3.1.1 [@doi:10.1093/molbev/msw046] and regressed using a linear model against coverage and contamination in seaborn v0.10.0 [@doi:10.5281/zenodo.3629446].

#### Plasmid and GI Coverage

Plasmid and GI coverage were assessed in the same way.
Firstly, a BLASTN database was generated for each set of MAG contigs.
Then each MAG database was searched for plasmid and GI sequences with greater than 50% coverage.
All plasmids or GIs which could be found in the unbinned contigs or MAGs were recorded as having been successfully assembled.
The subset of these that were found in the binned MAGs was then separately tallied.
Finally, we evaluated the proportion of plasmids or GIs that were correctly assigned to the bin that was maximally composed of chromosomes from the same source genome.

### Antimicrobial Resistance and Virulence Factors Assessment

#### Detection of AMR/VF Genes

For the reference genomes, as well as 12 sets of MAGs prodigal [@doi:10.1186/1471-2105-11-119] was used to predict open reading frames (ORFs) using the default parameters. 
AMR genes were predicted using Resistance Gene Identifier (RGI v5.0.0; default parameters) and the Comprehensive Antibiotic Resistance Database (CARD v3.0.3) [@doi:10.1093/nar/gkw1004].
Virulence factors were predicted using the predicted ORFs and BLASTX 2.9.0+ [@doi:10.1186/1471-2105-10-421] against the Virulence Factor Database (VFDB; obtained on Aug 26, 2019) with an e-value cut-off of 0.001 and a minimum identity of 90% [@doi:10.1093/nar/gky1080].
Each MAG was then assigned to a reference chromosome using the above mentioned mapping criteria for downstream analysis.

#### AMR/VF Gene Recovery

For each MAG set, we counted the total number of AMR/VF genes recovered in each metagenomic assembly and each MAG and compared this to the number predicted in their assigned reference chromosome and plasmids.
We then assessed the ability for MAGs to correctly bin AMR/VF genes of chromosomal, plasmid and GI origin by mapping the location of the reference replicon’s predicted genes to the location of the same genes in the MAGs.

#### Protein subcellular localization predictions

We then sought to assess what the impact of a protein's predicted subcellular localization was on its recovery and binning in MAGs.
The MAG bins from megahit-DAS Tool assembler-binner combination was selected (as generally best performing) and ORFs predicted using prodigal [@doi:10.1186/1471-2105-11-119] as above.
Subcellular localisation of these proteins were then predicted using PSORTb v3.0 with default parameters and the appropriate Gram setting for that bin's assigned taxa [@doi:10.1093/bioinformatics/btq249]. 




## Results {#results}

### Recovery of Genomic Elements
#### Chromosomes 

The overall ability of MAG methods to recapitulate the original chromosomal source genome results varied widely.
We considered the "identity" of a given MAG bin to be that of the genome that composes the largest proportion of sequence within that bin.
In other words if a bin is identifiably 70% species A and 30% species B we consider that to be a bin of species A.
Ideally, we wish to generate a single bin for each source genome comprised of the entire genome and no contigs from other genomes.
Some genomes are cleanly and accurately binned regardless of the assembler and binning method used (see Fig. @fig:speciescov).
Specifically, greater than 90% of _Streptomyces parvulus_ (minimum 91.8%) and _Clostridium baratii_ (minimum 96.4%) chromosomes are represented in individual bins across all methods.
However, no other genomes were consistently recovered at >30% chromosomal coverage across methods.
The three _Streptococcus_ genomes were particularly problematic with the best recovery for each ranging from 1.7% to 47.49%.
Contrary to what might be expected the number of closely relatives to a given genome in the metagenome did not clearly affect the MAG coverage (Fig. @fig:coverphylo).


![Top genome coverage for input genomes across MAG binners. Each dot represents the coverage of a specified genome when it comprised the plurality of the sequences in a bin. The binning tool is indicated by the colour of the dot as per the legend. Genomes such as _Clostridium baratti_ were accurately recovered across all binner-assembler combinations whereas genomes such as _Streptococcus macedonicus_ were systematically poorly recovered.](images/top_hits_per_bin.png){#fig:speciescov}

In terms of the impact of different metagenome assemblers, megahit resulted in the highest median chromosomal coverage across all binners (81.9%) with metaSPAdes performing worst (76.8%) (Fig. @fig:chromcover).
In terms of binning tool, CONCOCT performed very poorly with a median 26% coverage for top hit per bin, followed by maxbin2 (83.1%), and MetaBAT2 (88.5%).
It is perhaps unsurprising that the best-performing binner in terms of bin top hit coverage was the metabinner DASTool that combines predictions from the other 3 binners (94.3% median top hit chromosome coverage per bin; (Fig. @fig:chromcover)).

![Chromosomal coverages of most prevalent genome in each bin across binners and metagenome assemblies. Of the 3 assemblers (y-axis), megahit resulted in the highest median chromosomal coverage (x-axis) across all binners (colored bars) at 81.9% with metaSPAdes performing the worst (76.8%). Of the 4 binners, CONCOCT (blue) performed poorly with a median coverage, followed by maxbin2 (yellow), MetaBAT2 (red) and DASTool (green) performing the best. Diamonds in the figure represents outliers (greater or lower than the interquartile range marked by the error bars) and box represents the lower quartile, median, and upper quartile respectively.](images/bin_coverage.png){#fig:chromcover}

Bin purity, i.e. the number of genomes present in a bin at >5% coverage, was largely equivalent across assemblers, with a very marginally higher purity for IDBA.
In terms of binning tools, however, maxbin2 proved an exception with nearly twice as many bins containing multiple species as the next binner (Fig. @fig:purity).
The remaining binning tools were largely equivalent, producing chimeric bins at approximately the same rates.
Unlike coverage, purity was strongly affected by the number of close relatives in the metagenome to a given input genome. 
Specifically, the closer the nearest relative the less pure the bin (Fig. @fig:purityphylo).

![Distribution of bin purity across assemblers and binners. The total number of genomes present in a bin at >5% coverage (y-axis) was largely equivalent across assemblers (x-axis). In term of binning tools, maxbin2 (orange) produced nearly twice as many bins containing multiple species compared to CONCOCT (blue), MetaBAT2 (red) and DASTool (green), which all produced chimeric bins at roughly the same rate. Similar to above, outliers outside the interquartile range marked by the error bars are shown as diamonds.](images/bin_purity.png){#fig:purity}

#### Plasmids

Regardless of method, a very small proportion of plasmids were correctly grouped in the bin that was principally comprised of chromosomal contigs from the same source genome.
Specifically, between 1.5% (IDBA-UD assembly with DASTool bins) and 29.2% (metaSPAdes with CONCOCT bins) were correctly binned at over 50% coverage.
In terms of metagenome assembly, metaSPAdes was by far the most successful assembler at assembling plasmids with 66.2% of plasmids identifiable at greater than 50% coverage.
IDBA-UD performed worst with 17.1% of plasmids recovered, and megahit recovered 36.9%.
If the plasmid was successfully assembled, it was, with one exception, placed in a MAG bin by maxbin2 and CONCOCT, although a much smaller fraction were correctly binned (typically less than 1/3rd).
Interestingly, the MetaBAT2 and DASTool binners were more conservative in assigning plasmid contigs to bins; however, of those assigned to bins nearly all were correctly binned (Fig. @fig:plasmids).

![The performance of metagenomic assembly and binning to recover plasmid sequences. Each plot represents a different metagenome assembler, with the groups of bars along the x-axes showing the plasmid recovery performance of each binning tool when applied to the assemblies produced by that tool.  For each of these 12 assembler-binner-pair-produced MAGs the grouped bars from left to right show the percentage of plasmids assembled, assigned to any bin, and binned with the correct chromosomes.  These stages of the evaluation are indicated by the bar colours as per the legend.  Across all tools the assembly process resulted in the largest loss of plasmid sequences and only a small proportion of the assembled plasmids were correctly binned.](images/plasmid_recovery.png){#fig:plasmids}

#### Genomic Islands

GIs displayed a similar pattern of assembly and correct binning performance as plasmids (Fig. @fig:gis).
Assembly of GIs with >50% coverage was consistently poor (37.8-44.1%) with metaSPAdes outperforming the other two assembly approaches.
For the CONCOCT and maxbin2 binning tools, all GIs that were assembled were assigned to a bin, although the proportion of binned GIs that were correctly binned was lower than for DASTool and MetaBAT2.
DASTool, MetaBAT2 and CONCOCT did not display the same precipitous drop between those assembled and those correctly binned as was observed for plasmids.
In terms of overall correct binning with the chromosomes from the same genome the metaSPAdes assembly with CONCOCT (44.1%) and maxbin2 (43.3%) binners performed best.

![Impact of metagenomic assembly and MAG binning on recovery of genomic islands. GIs were recovered in a similarly poor fashion to plasmids. Generally, \<40% were correctly assigned to the same bin majorly comprised of chromosomal contigs from the same source genome regardless of binning (x-axis) and assembly (facet) methods at >50% coverage. metaSPAdes performed the best at assembling GIs (blue). Maxbin2 and CONCOCT placed GIs in a bin majority of the time (orange) however a very small fraction was correctly binned (green). Generally, GIs were correctly binned better than plasmids with DASTool, MetaBAT2 and CONCOCT.](images/GI_recovery.png){#fig:gis}

#### AMR Genes

The recovery of AMR genes in MAGs was poor with only ~49-55% of all AMR genes predicted in our reference genomes regardless of the assembly tool used, and metaSPAdes performing marginally better than other assemblers (Fig. @fig:AMRGenePercentRecoveryStage).
Binning the contigs resulted in a ~1-15% loss in AMR gene recovery with the CONCOCT-metaSPAdes pair performing best at only 1% loss and DASTool-megahit performing the worst at 15% reduction of AMR genes recovered.
Overall, only 24% - 40% of all AMR genes were correctly binned.
This was lowest with the maxbin2-IDBA-UDA pair (24%) and highest in the CONCOCT-metaSPAdes pipe (40%).

![Recovery of AMR genes across assemblers and binners. The proportion of reference AMR genes recovered (y-axis) was largely similar across assembly tools (blue), at roughly 50% with metaSPAdes performing marginally better. Binning tools resulted in a small reduction in AMR genes recovered (orange), however only 24-40% of all AMR genes were correctly binned (green). metaSPAdes-CONCOCT was the best performing MAG binning pipeline. ](images/amr_recovery.png){#fig:AMRGenePercentRecoveryStage}

Moreover, focusing on only the AMR genes that were correctly binned (Fig. @fig:AMRGenePercentRecoveryCorrectlyBinned) we can evaluate the impact of different genomic contexts (i.e. chromosomal, plasmid, GI).
Across all methods only 30%-53% of all chromosomally located AMR genes (n=120), 0-45% of genomic island located AMR genes (n=11) and none of the plasmid-localised AMR genes (n=20) were correctly binned.

![Percent of correctly binned AMR genes recovered by genomic context. MAG methods were best at recovering chromosomally located AMR genes (light blue) regardless of metagenomic assembler or binning tool used. Recovery of AMR genes in GIs showed a bigger variation between tools (light green). None of the 12 evaluated MAG recovery methods were able to recover plasmid located AMR genes (orange).](images/amr_localization_recovery.png){#fig:AMRGenePercentRecoveryCorrectlyBinned}

#### Virulence Factor Genes

We also examined the impact of MAG approaches on recovery of virulence factor (VF) genes as identified using the Virulence Factor Database (VFDB).
We saw a similar trend as AMR genes (Fig. @fig:VFGenePercentRecoveryStage).
Between 56% and 64% of VFs were identifiable in the metagenomic assemblies (with megahit recovering the greatest proportion).
The binning process further reduced the number of recovered VFs by 4-26% with DASTool-megahit performing the worst (26%) and CONCOCT-metaSPAdes performing the best (4%).
Unlike AMR genes, the majority of VF genes assigned to a bin were assigned to the correct bin (i.e. that bin largely made up of contigs from the same input genome).
Overall, CONCOCT-metaSPAdes again performed best with 43% of all VFs correctly assigned.

![Percent of reference virulence factor (VF) genes recovered across assemblers and binners. The proportion of reference VF genes recovered (y-axis) exhibited a similar trend as AMR genes. Recovery was greatest after the assembling stage (blue), with megahit performing best. Binning tools resulted in a larger reduction in VF genes recovered (orange) compared to AMR genes. However, in majority of cases, VF genes that are binned are correctly binned (green). metaSPAdes-CONCOCT was again the best performing pair.](images/vf_recovery.png){#fig:VFGenePercentRecoveryStage}

As with AMR genes, the genomic context (chromosome, plasmid, GI) of a given VF largely determined how well it was binned (Fig. @fig:VFGenePercentRecoveryStage).
The majority (73%-98%) of all chromosomally located VF genes (n=757) were correctly binned.
However, 0-16% of GI-localised VF genes (n=809) and again none of the plasmid-associated VF genes (n=3) were recovered across all 12 MAG pipelines.

![Percent of correctly binned VF genes recovered in each genomic region. Metagenome assembled genomes (MAGs) were again best at recovering chromosomally located VF genes (light blue), able to correctly bin majority of chromosomally located VFs. GIs recovered again performed very poorly (light green) and again none of the plasmid located AMR genes (orange) was correctly binned.](images/vf_localization_recovery.png){#fig:VFGenePercentRecoveryStage}

### Comparisons of Rates of Loss

We combined the performance metrics for Figs. @fig:plasmids, @fig:gis, @fig:AMRGenePercentRecoveryStage, and @fig:VFGenePercentRecoveryStage to compare the rates of loss of different components (see Fig. @fig:rateofloss).
This highlighted that genomic components (GIs and plasmids) and plasmids in particular are lost at a higher rate than individual gene types during MAG recovery.



## Discussion {#discussion}

In this paper, we evaluated the ability and accuracy of metagenome-assembled genome (MAGs) binning methods to correctly recover mobile genetic elements (i.e. genomic islands and plasmids) from metagenomic samples across different tools used to assemble and bin MAGs.

Overall, the best assembler-binner pair was megahit-DASTOOL in terms of both chromosomal coverage (94.3%) and bin purity (1).
Looking at genomes with the lowest coverage, the three Streptococcus genomes that were recovered poorly are likely due to their similarity (Fig. @fig:coverphylo, @fig:purityphylo).
This supports the intuition that MAG recovery approaches struggle to distinguish closely related species.
While CONCOCT performed significantly worse than other binners in terms of chromosomal coverage and bin purity, we did notice that CONCOCT was prone to generating many small partial bins.
Potentially, CONCOCT binning could be used to distinguish closely related species but at a cost of more fragmented genomes.

While the overall recovery and binning of chromosomes was likely sufficient for some use-cases, we were specifically interested in the ability of MAG methods to appropriately recover MGEs. This was due to the importance of MGEs in the function and spread of pathogen traits such as AMR and virulence, as well as our hypothesis that these sequences may prove difficult to bin. Regardless of the the metagenomic assembly approach or MAG binning method used, both plasmids and GIs were disproportionately lost compared to chromosomes in general. At best (with metaSPAdes and CONCOCT) 29.2% of plasmids and 44.1% of GIs were identifiable at >50% coverage in the correct bin (i.e. grouped with a bin that was mostly made up of contigs from the same genome). The >50% coverage requirement set a high bar and more GIs and plasmids were likely recovered in more incomplete forms. 
Partial MGEs may be useful for some research, but for researchers interested in selective pressures and lateral gene transfer this may lead to inaccurate inferences. 

This poor result is not unexpected as genomic islands and plasmids have known divergent compositional features and are often repetitive with variable copy numbers relative to the chromosome.
Furthermore, the difference between the percentages suggests that binning plasmids is harder than binning GIs.
This difference might be attributed to the known difficulties in assembly of plasmids from short-read data [@pmid:29177087].
Therefore, binning efficiency might improve if we use DNA sequencing and assembly methods optimised for recovering plasmids [@doi:10.1099/mgen.0.000128] (such as SCAPP [@doi:10.1101/2020.01.12.903252]).

Due to the importance of MGEs in the dissemination of clinically relevant AMR genes and VFs, we explored whether or not MAG approaches can be used to provide useful insight into the LGT of these genes.
With respect to AMR genes, MAG methods were able to recover roughly 40% of all AMR genes present in our reference genomes.
We noted a sharp drop in the number of AMR genes detected between assemblies and MAGs, suggesting that many of these genes were left in the unbinned portion.
Overall, the CONCOCT-metaSPAdes combination, while it did not recover the highest amount of AMR genes at the assembly stage, performed the best in correctly binning an AMR gene to the right species.
Regardless of tools, chromosomally located AMR genes were most frequently correctly binned (as expected from the relative performance of MAGs at recovering chromosomes).
While there was variability in performance, AMR genes located on GIs were correctly binned slightly less well than chromosomally located AMR genes.
This variability might be explained by the fact that there were only 11 AMR genes located on GIs in our reference genomes.
All 20 of the plasmid-borne AMR genes were assembled, but none were placed into MAG bins.
We included high-threat MGEs-associated AMR genes such as the KPC and OXA carbapenemases. 
We intended on a systematic review of which AMR genes are more or less likely end up correctly binned, however, MAGs was not able to correctly bin enough AMR genes on plasmids or GIs to allow this. 

Virulence factors showed a similar trend to the AMR genes, with a recovery of ~63% of virulence factors present in the reference genomes.
There still is a sharp decline in the number of VF detected between assemblies and MAGs and CONCOCT-metaSPAdes again produced the highest binning accuracy.
A majority (73%-98%) of chromosomally located VF genes were also able to be correctly binned to the right species for the MAGs.
However, the MAG approach performed much worse in correctly recovering GI located and plasmid located VFs, with <16% of GI VFs (n=809) correctly recovered and none of the plasmid VFs (n=3).
This drastic reduction in recovery accuracy of mobile elements, especially GIs, is expected.
Previous studies have found that VFs are disproportionally present on GIs[@doi:10.1371/journal.pone.0008094], which might be the reason why the recovery accuracy was worse compared to AMR genes.
Together, this and the AMR gene results suggests that MAG-based methods might be of limited utility in public health research focused on the transmission and dissemination of AMR genes and VFs. 

One potential caveat is that some AMR genes and VFs successfully assembled in the MAGs may no longer be annotated as such due to issues with ORF prediction (see suppl. discussion & Fig. @fig:geneContent). 
Previous studies have observed that ORF predictions in draft genomes are more fragmented, which can lead to downstream over- or under-annotation with functional labels depending on the approach used [@doi:10.1186/1471-2164-13-14]. 
Similarly, if the ORFs predicted in the MAGs differ in sequence or degree of fragmentation from the corresponding ORFs predicted in the original reference genomes (or are no longer predicted at all), this could impact recovery of AMR/VF predictions, even though the sequences themselves may be partially or fully present in the assembly.

It should also be noted that while CONCOCT performed the best in terms of recovery of both chromosomes and MGEs, it created many relatively clean but fragmentary partial MAGs.
While this might be ideal for some users, caution should be taken in using CONCOCT when assuming a bin represents a whole genome.

With the recovery of plasmids, GIs, VFs, and AMR genes the same pattern was observed, a progressive loss of data in each analytical step.
The process of metagenomic assembly itself generally resulted in the loss of most of these elements/genes regardless of the assembly method used.
With repetitive DNA sequence particularly difficult to correctly assemble from short reads [@doi:10.1093/jac/dkw184,@doi:10.1099/mgen.0.000128].
Across binning tools, the binning process resulted in further loss with a large proportion of MGEs and genes left unbinned.
Finally, only a very small proportion of these elements/genes were generally correctly binned with the appropriate host chromosomes.
This follows the well known, but rarely explicitly stated, idea that the more analysis you perform the more of the original data gets lost.
Indeed, this is one of the reasons why the huge amount of redundancy in metagenomic sequencing is necessary (i.e. many more base-pairs of DNA must be sequenced than are in the underlying sample).


## Conclusions {#conclusions}

Using a simulated medium-complexity metagenome, this study has shown that MAG-based approaches provide a useful tool to study a bacterial species’ chromosomal elements, but have severe limitations in the recovery of MGEs.
The majority of these MGEs will either both fail to assemble or be incorrectly binned.
The consequence of this is the disproportionate loss of key public health priority genes like VF and AMR genes.
This is particularly acute as the VF and AMR genes found on these poorly recovered MGEs are generally considered the most important due to their propensity for lateral gene transfer between unrelated bacteria.
Therefore, it is vital that we utilize a combination of MAGs and other methods (e.g. read-based methods) in public health metagenomic research when short-read sequencing is used.
For example, targeted AMR [@doi:10.1099/mgen.0.000131], plasmid specialised assembly approaches [@doi:10.1101/2020.01.12.903252], and read-based sequence homology search [@doi:10.1038/nmeth.3176].
Without this, MAG-based methods are insufficient to thoroughly profile the resistome and provide vital epidemiological data for metagenomic data.


## Supplementals {#supplementals}

### Recovery of Specific Gene Content

We explored the ability of different approaches to find open reading frames (ORFs) within MAGs.
Overall, the total number of predicted ORFs in MAGs followed a similar trend (Fig. @fig:geneContent) as the chromosomal coverage (Fig. @fig:chromcover) and purity (Fig. @fig:purity).
Of the four binning tools, CONCOCT performed the worst, finding <30% of the number of ORFs in our reference genomes used to construct the synthetic data.
MetaBAT2 performed second worst at ~80%.
DASTool recovered a similar number to our reference and Maxbin2 detected 7-46% more genes.
The Assembler method did not significantly impact the number of genes predicted with the exception of Maxbin2, in which IDBA_UD was the closest to reference and metaSPAdes predicted 46% more ORFs.
Given that there is reason to suspect that there are some issues with the ORF calling in the MAGs. i.e. some tools produced more predicted ORFs than reference, it could be the case that some of these sequences are present in the assemblies (with errors/gaps), but are not being identified as ORFs, or are broken into mulpitle ORFs, leading to issues downstream labeling them correctly as AMR/VF genes. 
Regardless of different tools producing a different number of ORFs, the recovery of AMR/VF is pretty consistent regardless of how many ORFs are predicted.

![Predicted Gene Content. The total number of open reading frames (ORF) predicted followed the same trend as chromosomal coverage and purity. The assemblers (colored bars) did not contribute to variability in the number of ORFs detected. Of the 4 binners, CONCOCT recovered \<30\% of our reference genome ORFs. DASTool and MetaBAT2 predicted a similar number as our reference genomes.](images/number_of_predicted_genes.png){#fig:geneContent}

### Impact of Related Genomes on MAG

By generating a phylogeny of universal single copy genes in our input genomes we analysed the relationship between the presence of closely related genomes and the ability of the different MAG-recovery methods to bin chromosomal sequences.
Specifically, we regressed phylogenetic distance on this phylogeny with per-bin chromosomal coverage (Fig. @fig:coverphylo) and bin purity (Fig. @fig:purityphylo).
This identified no clear relationship between chromosomal coverage and the phylogenetic distance to the nearest relative in the metagenome (Fig. @fig:coveragephylo), however, there did seem to be a negative correlation between phylogenetic distance to closest relative and the purity of a MAG bin (Fig. @fig:purityphylo).
In other words, across all methods, a MAG bin was more likely to have multiple genomes present if there were close relatives.

![Relationship between phylogenetic distance to closest neighbour input genome on genomic coverage in MAG majority comprised of that taxa. Each dot represents the genomic coverage of a particular taxa and the branch distance on an 86-protein concatenated phylogeny between that taxa and its nearest neighbour. Rows indicate the binning software and columns the metagenomic assembler. Regression line is a simple linear model fitted in seaborn. ](images/coverage_phylo_dist.png){#fig:coverphylo}

![Relationship between phylogenetic distance to closest neighbour input genome on bin purity.  Each dot shows the number of other input genomes detectable in a given MAG bin in relation to the branch distance on an 86-protein concatenated phylogeny between the majority taxa in that bin and its nearest neighbour.](images/purity_phylo_dist.png){#fig:purityphylo}

### Comparisons of Rates of Loss

Combining the performance metrics for Figs. @fig:plasmids, @fig:gis, @fig:AMRGenePercentRecoveryStage, and @fig:VFGenePercentRecoveryStage to compare the rates of loss of different components emphasises some of the observed patterns (see Fig. @fig:rateofloss).
This highlights that genomic components (GIs and plasmids) and plasmids in particular are lost at a higher rate than individual gene types during MAG recovery.

![Comparison of rates of loss for different genomic components and gene types across assemblers and binning tools. Each line represents a different component as indicated by the legend with assemblers indicated by row and binning tool by column. This shows that regardless of approach genomic components (GIs and plasmids) are lost at a higher rate than individual VF or AMR genes.](images/rate_of_loss.png){#fig:rateofloss}



[@tag:deep-review]: doi:10.1098/rsif.2017.0387