---
author-meta:
- Finlay Maguire*
- Baofeng Jia*
- Kristen Gray
- Venus Lau
- Robert G. Beiko
- Fiona S.L. Brinkman
date-meta: '2020-01-24'
keywords:
- markdown
- publishing
- manubot
lang: en-CA
title: Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids
  and Genomic Islands
...






<small><em>
This manuscript
([permalink](https://fmaguire.github.io/mag_sim_paper/v/a9f44e07ce1a0f52820132149a6ecace9264a152/))
was automatically generated
from [fmaguire/mag_sim_paper@a9f44e0](https://github.com/fmaguire/mag_sim_paper/tree/a9f44e07ce1a0f52820132149a6ecace9264a152)
on January 24, 2020.
</em></small>

## Authors
Please note the current author order is chronological and does not reflect the final order.



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
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [imasianxd](https://github.com/imasianxd)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>

+ **Kristen Gray**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>

+ **Venus Lau**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>

+ **Robert G. Beiko**<br><br>
  <small>
     Faculty of Computer Science, Dalhousie University
  </small>

+ **Fiona S.L. Brinkman**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>




## Abstract {#abstract}

Metagenomic methods, in which all the DNA in sample is simultaneously sequenced, is an increasingly popular method in the life sciences. 
They have a major advantage over genomic or phenotypic methods, as they do not require time-intensive and bias-inducing culturing steps. 
This means a much greater diversity can be achieved with minimal _a priori_ assumptions. 
Due to this strength, metagenomics is emerging as a key tool in public health microbiology for surveillance of virulence and antimicrobial resistance (AMR) genes. 
The most important sequences for surveillance purposes are those associated with mobile genetic elements such as plasmids and genomic islands (GIs). 
Unfortunately, metagenomic data, even when assembled, results in complex mixed set of DNA fragments rather than nicely resolved individual genomes. 
Recently, methods have been developed that attempt to group these DNA fragments into bins that likely to have been derived from the same underlying genome. 
These bins are referred to as metagenome-assembled genomes (MAGs). 
MAG based approaches have been used to great effect in revealing huge amounts of previously uncharacterised microbial diversity.
These methods perform this grouping using aspects of the sequence composition and the relative abundance of that sequence in the dataset. 
Unfortunately, plasmids are often represented at different copy numbers than the corresponding chromosomes. 
Additionally, both plasmids and genomic islands often feature significantly different sequence composition than the rest of the source genome as a whole. 
Due to this we hypothesise, these types of sequences will be highly under-represented in MAG based approaches.

To evaluate this we generated a simulated metagenomic dataset comprised of 30 genomes with up to 16.65% of chrosomomal DNA consisting of GIs and 65 associated plasmids. 
MAGs were then recovered from this data using 12 different MAG pipelines. 
The recovery and correct binning of mobile genetic elements was then evaluated for each pipeline. 
Across all pipelines, 81.9-94.3% of chromosomes were recovered and binned. However, only 37.8-44.1% of GIs and 1.5-29.2% of plasmids were recovered and correctly binned at >50% coverage. 
In terms of AMR and VF genes associated with MGEs, 0-45% of GI-associated AMR genes and 0-16% of GI-associated VF genes were correctly assigned. 
More strikingly, 0% of plasmid-borne VF or AMR genes were recovered.

This work shows that regardless of the MAG recovery approach used, plasmid and GI dominated sequences will disproportionately be left unbinned or incorrectly binned. 
From a public health perspective, this means MAG approaches are unsuited for analysis of mobile genes, especially vital groups such as AMR and VF genes. 
This underlines the utility of read-based and long-read approaches to thorougly evaluate the resistome in metagenomic data.

## Introduction {#intro}

Metagenomics, the untargeted sequencing of all DNA within a sample, has become the dominant approach for characterising viral and microbial communities over the last 17 years [@8PLOeAH6; @7RV1Ygsv]. 
By sampling from the total genomic content these methods allow researchers to simultaneously profile the functional potential and the taxonomic identity of all organisms ina sample. 
This is in contrast to bar-coding based approaches such as 16S or 18S rRNA sequencing which only provide taxonomic information [@Trbw5ZnC] (although you can attempt to predict functional potential from taxonomic data [@93ojQnSg; @L09Dmaq3]). 
One of many areas where metagenomics has been very useful is in the analysis of antimicrobial resistance (AMR) and pathogen virulence. 
These approaches have been instrumental in developing our understanding of the distribution and evolutionary history of AMR genes [@rwhLEYRY; @QSe5BqFk; @2xaXclNM]. 
It has also formed a key tool for pathogen tracking in public health outbreak analyses [@khJQfjDf].

While 3rd generation long-read technology has begun to be adopted in metagenomics analyses [@U4vhNZoB; @1759XyDVi] the majority of analyses still involve high-throughput 2nd generation sequencing. 
These 2nd generation platforms such as Illumina's MiSeq provide high numbers (10s-100s of millions) of relatively short reads (150-250bp) randomly sampled from the underlying DNA in the sample. 
This sampling is, therefore, in proportion to the relative abundance of different organisms (i.e. more abundant organisms will be more represented in the reads). 
There are two main approaches for the analysis of 2nd generation metagenomic data: read homology and metagenome assembly. 
Read-based approaches involve using reference databases and BLAST-based sequence similarity search tools (e.g. DIAMOND [@4R96QRcV]), read mapping (e.g. Bowtie 2 [@PiS0h6Mu]), Hidden Markov Models (e.g. HMMER3 [@77xWEk9S]) or k-mer hashing (e.g. CLARK [@OoKZ0WcH]). 
These read-based approaches allow analysis of all reads with detectable similarity to genes of interest even if the organism has relatively low abundance in the sample.
However, read-based methods are reliant on quality of the reference database (i.e. you don't detect things you don't already know about) and does not provide any information about the genomic organisation of the genes.
This lack of contextual information is particularly problematic in the study of AMR genes and virulence factors as the genomic context plays a role in function [@17Dww6tOF], selective pressures [@1BcFmOfCH], and how liable the sequence is to lateral gene transfer (LGT) [@SACLvb9k].

In order to get more data about the relative genomic context and organisation of your genes of interest it is possible (although computationally demanding) to assemble the short reads into longer fragments of DNA (contigs). 
This approach has been used successfully since early metagenomic analysis papers [@F7RexqdF]. 
There are a variety of specialised _de Bruijn_ graph assemblers developed to handle the particular challenges of this type of assembly (such as metaSPAdes [@KP5SjPXN] , IDBA-UD [@a4mT7fuU], and megahit [@1EUV0Ejkr]) each with a range of different strengths and weaknesses [@f0X2CPKM]. 
While this approach does result in longer contigs it still leaves you with a large collection of fragmentary data derived from many different organisms.

An increasingly common way to deal with this is to attempt to group these assembled contigs into bins all derived from the same underlying genome in the sample. 
These resulting bins are known as metagenome assembled genomes (MAGs). 
This binning is typically performed by grouping all the contigs with similar abundance and similar sequence composition into the same bin. 
A range of tools have been released to perform this binning including CONCOCT [@WRoCf6pg], MetaBAT 2 [@b2WO18xh], and MaxBin 2 [@sG4CX8Uj]. 
There is also the meta-binning tool DAS Tool [@DfIRBmdF] which combines predictions from multiple binning tools together. 
These MAG approaches have been used to great e ect in unveiling huge amounts of previously uncharacterised genomic diversity [@4rsFboY4; @wrBRBdFb; @Rk2NATlI].

Unfortunately, there is loss of information at both the metagenomic assembly step (e.g. repetitive DNA sequences that are difficult to correctly assemble with short reads) [@UWOTvAMZ,@12zFifp5x] and in binning.
This compounded data loss means that only a relatively small proportion of reads are successfully assembled and binned in large complex metagenome datasets e.g. 24.2-36.4% of reads from permafrost [@buqrbdBh] and soil metagenomes [@d5Hh0941]. 
Additionally, a large number of detected genomes are not reconstructed at all with ~23% of all detected genomes recovered in some examples [@d5Hh0941]. 
There have been attempts to benchmark and compare the assembly and binning tools such as the Critical Assessment of Metagenome Interpretation (CAMI) challenge's (https://data.cami-challenge.org/) Assessment of Metagenome BinnERs (AMBER) [@Y8sHlHi] however these largely investigate the overall completeness and purity of recovered MAGs relative to the known genomes in the evaluation samples.
To our best knowledge, there hasn't been a specific assessment of the impact of metagenomic assembly and binning on the loss of specific genomic elements.
Two such genomic elements of great health and research importance are mobile genetic elements (MGEs) such as genomic islands (GIs) and plasmids.

Genomic islands (GIs) are clusters of genes that are known or predicted to have been acquired through LGT events. 
These include integrons, transposons, integrative and conjugative elements (ICEs) and prophages (integrated phages) [@DET3tBYj; @1Af4oXwEX]. 
They have been shown to disproportionately encode virulence factors [@LxGqo7iq] and are a major mechanism of LGT of AMR genes [@x7HhCKyS; @17U91060Y]. 
However, these GIs often have different nucleotide composition compared to the rest of the genome [@DET3tBYj].
This compositional difference is exploited by tools designed to detect GIs such as SIGI-HMM [@19UeQywMr] and IslandPath-DIMOB [@M1pdcdMy]. 
GIs may exist as multiple copies within a genome [@5g9Xc4ot] leading to potential difficulties in correctly assembling these regions in metagenome assemblies as well as likely biases in the calculation of coverage statistics. 
Similarly, plasmids, circular or linear extrachromosomal self-replicating pieces of DNA, are a major source of the dissemination of AMR genes throughout microbial ecosystems [@X9j9vETu; @x7HhCKyS]. 

Due to their research importance, lots of work has identified the difficulty of assembling these sequences correctly from short-read data [@12zFifp5x,@bkwNETh8].
This is largely attributable to their repetitive sequences, variable copy number [@qtpTcNWp; @Z1irb7eF] and often markedly different sequence composition to the genome they are associated with [@QK9dmRUA; @ps1aOiRU]. 
As MAG binning is performed on the basis of sequence composition and relative abundance this suggests that these types of sequences are liable to being incorrectly binned or lost in the process of recovering MAGs. 
As these MGEs are key to the function and spread of pathogenic traits such as AMR and virulence it is vital that we assess the impact of metagenome assembly and binning on the representation of these specific elements. 
This is particularly important with the increasing popularity of MAG approaches within microbial and public health research. 
Therefore, to address this issue we performed an analysis of GI and plasmid recovery accuracy across a broad-set of current state-of-the-art short-read metagenome assembly and binning approaches using a simulated medium complexity metagenome comprised of GI- and plasmid-rich taxa.

## Materials and Methods {#methods}

All analyses presented in this paper can be reproduced and inspected with the associated github repository [github.com/fmaguire/MAG_gi_plasmid_analysis](github.com/fmaguire/MAG_gi_plasmid_analysis) and data repository [osf.io/nrejs/](osf.io/nrejs/).

### Metagenome Simulation

All genomes were selected from the set of completed RefSeq genomes as of April 2019.
Genomic islands for these genomes were previously predicted using IslandPath-DIMOB [@M1pdcdMy] and collated into the IslandViewer database [www.pathogenomics.sfu.ca/islandviewer](www.pathogenomics.sfu.ca/islandviewer) [@4eEyIkDg].
Plasmid sequences and numbers were recovered for each genome using the linked GenBank Project IDs.
Thirty genomes were manually selected to exemplify the following criteria: 

	1) 10 genomes with high numbers of plasmids.

	2) 10 genomes with a very high proportion (>10%) of chromosomes corresponding to GIs detected by compositional features.

	3) 10 genomes with a very low proportion (<1%) of chromosomes corresponding to GIs detected by compositional features.

The data used to select the taxa is listed in Supplemental Table 1 and the details of the selected subset taxa are listed in Supplemental Table 2 with their NCBI accessions.
The sequences themselves are available in the data repository [osf.io/nrejs/](osf.io/nrejs/) under "data/sequences".

In accordance to the recommendation in the CAMI challenge [@lsbnKJf8] the genomes were randomly assigned a relative abundance following a log-normal distribution (μ = 1, σ = 2).
Plasmid copy number estimates could not be accurately found for all organisms, therefore, plasmids were randomly assigned a copy number regime: low (1-20), medium (20-100), or high (500-1000) at a 2:1:1 rate.
Within each regime, the exact copy number was selected using an appropriately scaled gamma distribution (α = 4, β = 1) or the minimum edge of the regime.

Finally, the effective plasmid relative abundance was determined by multiplying the plasmid copy number with the genome relative abundance.
The full set of randomly assigned relative abundances and copy numbers can be found in Supplemental Table 3.
Sequences were then concatenated into a single FASTA file with the appropriate relative abundance.
MiSeq v3 250bp paired-end reads with a mean fragment length of 1000bp (standard deviation of 50bp) were then simulated using art_illumina (v2016.06.05) [@znONJtTo] at a fold coverage of 2.9 resulting in a simulate metagenome of 31,174,411 read pairs.
The selection of relative abundance and metagenome simulation itself was performed using the "data_simluation/simulate_metagenome.py" script.

### Metagenome Assembled Genome Recovery

Reads were trimmed using sickle (v1.33) [@1CBlSILo4] resulting in 25,682,644 surviving read pairs.
The trimmed reads were then assembled using 3 different metagenomic assemblers: metaSPAdes (v3.13.0)[@KP5SjPXN], IDBA-UD (v1.1.3) [@a4mT7fuU], and megahit (v1.1.3) [@1EUV0Ejkr]).
The resulting assemblies were summarised using metaQUAST (v5.0.2) [@TeRvtMCl].
The assemblies were then indexed and reads mapped back using Bowtie 2 (v2.3.4.3) [@PiS0h6Mu].

Samtools (v1.9) was used to sort the read mappings and the read coverage calculated using the MetaBAT2 accessory script (jgi\_summarize\_bam\_contig\_depths).
The three metagenome assemblies were then separately binned using MetaBAT2 (v2.13) [@b2WO18xh], and MaxBin 2 (v2.2.6) [@sG4CX8Uj]. 
MAGs were also recovered using CONCOCT (v0.4.2) [@WRoCf6pg] following the recommended protocol in the user manual.
Briefly, the supplied CONCOCT accessory scripts were used to cut contigs into 10 kilobase fragments (cut\_up\_fasta.py) and read coverage calculated for the fragments (CONCOCT\_coverage\_table.py).
These fragment coverages were then used to bin the 10kb fragments before the clustered fragments were merged (merge\_cutup\_clustering.py) to create the final CONCOCT MAG bins (extra\_fasta\_bins.py).
Finally, for each metagenome assembly the predicted bins from these three binners (Maxbin2, MetaBAT 2, and CONCOCT) were combined using the DAS Tool (v1.1.1) meta-binner [@DfIRBmdF].
This resulted in 12 separate sets of MAGs (one set for each assembler and binner pair).

### MAG assessment

#### Chromosomal Coverage

The MAG assessment for chromosomal coverage was performed by creating a BLASTN 2.9.0+ [@nEsJGUWa] database consisting of all the chromosomes of the input reference genomes.
Each MAG contig was then used as a query against this database and the coverage of the underlying chromosomes tallied by merging the overlapping aligning regions and summing the total length of aligned MAG contigs.
The most represented genome in each MAG was assigned as the “identity” of that MAG for further analyses.
Coverages less than 5% were filtered out and the number of different genomes that contigs from a given MAG aligned to were tallied.
Finally, the overall proportion of chromosomes that were not present in any MAG were tallied for each binner and assembler.

#### Plasmid and GI Coverage

Plasmid and GI coverage were assessed in the same way.
Firstly, a BLASTN database was generated for each set of MAG contigs.
Then each MAG database was searched for plasmid and GI sequences.
Any plasmid or GI with greater than 50% coverage in a MAG was retained.
All plasmids or GIs which could be found in the unbinned contigs or the MAGs was recorded as having been successfully assembled.
The subset of these which were found in the binned MAGs was then separately tallied.
Finally, we evaluated the proportion of plasmids or GIs that were binned correctly in the bin that was maximally composed of chromosomes from the same source genome.
This was determined using the bin “identities” from the chromosomal coverage analysis.

### Antimicrobial Resistance and Virulence Factors Assessment

#### Detection of AMR/VF Genes

For each of the 12 MAGs sets, and the reference chromosome and plasmids, prodigal [@lX665mdh] was used to predict open reading frames (ORFs) using the default parameters. 
AMR genes were predicted using Resistance Gene Identifier (RGI v5.0.0; default parameters) and the Comprehensive Antibiotic Resistance Database (CARD v3.0.3) [@1ByMfX8Y1].
Virulence factors were predicted using the predicted ORFs and BLASTX 2.9.0+ [@nEsJGUWa] against the Virulence Factors Database (VFDB; obtained on Aug 26, 2019) with an e-value cut-off of 0.001 and a minimum identity of 90% [@pYB1SP5].
Each MAG was then assigned to a reference chromosome using the above mentioned mapping criteria for downstream analysis.

#### AMR/VF Gene Recovery

For each MAG set, we counted the total number of AMR/VF genes recovered in each assembly and each MAG and compared this number to the number predicted in their assigned reference chromosome and plasmids.
We then assessed the ability for MAGs to correctly bin AMR/VF genes of chromosomal, plasmid and GI origin by mapping the location of the reference replicon’s predicted genes to the location of the same genes in the MAGs.

#### Protein subcellular localization predictions

We then sought to assess what the impact of a proteins predicted subcellular localization was on its recovery and binning in MAGs.
The MAG bins from megahit-DAS Tool assembler-binner combination was selected (as generally best performing) and ORFs predicted using prodigal [@lX665mdh] as above.
Subcellular localisation of these proteins were then predicted using PSORTb v3.0 with default parameters and the appropriate Gram setting for that bin's assigned taxa [@19bHWbO47]. 




## Results {#results}

### Recovery of Genomic Elements
#### Chromosomes 

The overall ability of MAG methods to recapitulate the original chromosomal source genome results varied widely.
We considered the "identity" of a given MAG bin to be that of the genome that composes the largest proportion of sequence within that bin.
In other words if a bin is identifiably 70% species A and 30% species B we consider that to be a bin of species A.
Ideally, we wish to generate a single bin for each source genome comprised of the entire genome and no contigs from other genomes.
Some genomes are cleanly and accurately binned regardless of the assembler and binning method used (see Fig. @fig:speciescov).
Specifically, greater than 90% of Streptomyces parvulus (minimum 91.8%) and Clostridium baratii (minimum 96.4%) chromosomes are represented in individual bins across all methods.
However, no other genomes were consistently recovered by all methods for more than a 1/3rd of the chromosomes.
The three Streptococcus genomes were particularly problematic with the best recovery for each ranging from 1.7% to 47.49%.

![Top genome coverage for input genomes across MAG binners. Each dot represents the coverage of a specified genome when it comprised the majority of the sequences in a bin. The binning tool is indicated by the colour of the dot as per the legend. Genomes such as _Clostridium baratti_ were accurately recovered across all binner-assembler combinations whereas genomes such as _Streptococcus macedonicus_ were systematically poorly recovered.](images/top_hits_per_bin.png){#fig:speciescov}

In terms of the impact of different metagenome assemblers, megahit resulted in the highest median chromosomal coverage across all binners (81.9%) with metaSPAdes performing worst (76.8%) (Fig. @fig:chromcover).
In terms of binning tool, CONCOCT performed very poorly with a median 26% coverage for top hit per bin, followed by maxbin2 (83.1%), and MetaBAT2 (88.5%).
It is perhaps unsurprising that the best performing binner in terms of bin top hit coverage was the metabinner DASTool that combines predictions from the other 3 binners (94.3% median top hit chromosome coverage per bin; (Fig. @fig:chromcover)).

![Chromosomal coverages of most prevalent genome in each bin across binners and metagenome assemblies. Of the 3 assemblers (y-axis), megahit resulted in the highest median chromosomal coverage (x-axis) across all binners (colored bars) at 81.9% with metaSPAdes performing the worst (76.8%). Of the 4 binners, CONCOCT (blue) performed poorly with a median coverage, followed by maxbin2 (yellow), MetaBAT2 (red) and DASTool (green) performing the best. Diamonds in the figure represents outliers (greater or lower than the interquartile range marked by the error bars) and box represents the lower quartile, median, and upper quartile respectively.](images/bin_coverage.png){#fig:chromcover width="5in"}


Bin purity, i.e. the number of genomes present in a bin at >5% coverage, was largely equivalent across assemblers (Fig. @fig:purity), with a very marginally higher purity for IDBA.
In terms of binning tool, however, maxbin2 proved an outlier with nearly twice as many bins containing multiple species as the next binner.
The remaining binning tools were largely equivalent, producing chimeric bins at approximately the same rates.

![Distribution of bin purities across assemblers and binners. The total number of genomes present in a bin at >5% coverage (y-axis) was largely equivalent across assemblers (x-axis). In term of binning tools, maxbin2 (orange) produced nearly twice as many bins containing multiple species compared to CONCOCT (blue), MetaBAT2 (red) and DASTool (green), which all produced chimeric bins at roughly the same rate. Similar to above, outliers outwith of the interquartile range marked by the error bars are shown as diamonds.](images/bin_purity.png){#fig:purity width="5in"}

#### Plasmids

Regardless of method, a very small proportion of plasmids were correctly grouped in the bin that was principally comprised of chromosomal contigs from the same source genome.
Specifically, between 1.5% (IDBA-UD assembly with DASTool bins) and 29.2% (metaSPAdes with CONCOCT bins) were correctly binned at over 50% coverage.
In terms of metagenome assembly, metaSPAdes was by far the most successful assembler at assembling plasmids with 66.2% of plasmids identifiable at greater than 50% coverage.
IDBA-UD performed worst with 17.1% of plasmids recovered, and megahit recovered 36.9%.
If the plasmid was successfully assembled, it was fairly consistently placed in a MAG bin by maxbin2 and CONCOCT, although a much smaller fraction were correctly binned (typically less than 1/3rd).
Interestingly, MetaBAT2 and DASTool binners were a lot more conservative in assigning plasmid contigs to bins; however, of those assigned to bins nearly all were correctly binned (Fig. @fig:plasmids)

![The performance of metagenomic assembly and binning to recover plasmid sequences. Each plot represents a different metagenome assembler, with the groups of bars along the x-axes showing the plasmid recovery performance of each binning tool when applied to the assemblies produced by that tool.  For each of these 12 assembler-binner pair produced MAGs the grouped bars from left to right show the percentage of plasmids assembled, binned in any bin, and binned with the correct chromosomes.  These stages of the evaluation are indicated by the bar colours as per the legend.  Across all tools the assembly process resulted in the largest loss of plasmid sequences and only a small proportion of the assembled plasmids were correctly binned.](images/plasmid_recovery.png){#fig:plasmids width="5in"}

#### Genomic Islands

GIs displayed a similar pattern of assembly and correct binning performance as plasmids (Fig. @fig:gis).
These sequences were assembled uniformly badly (37.8-44.1%) with metaSPAdes outperforming the other two assembly approaches.
For CONCOCT and maxbin2 binning tools all GIs that were assembled were assigned to a bin although the proportion of binned GIs that were correctly binned was lower than for DASTool and MetaBAT2.
DASTool, MetaBAT2 and CONCOCT did not display the same precipitous drop between those assembled and those correctly binned as was observed for plasmids.
In terms of overall correct binning with the chromosomes from the same genome the metaSPAdes assembly with CONCOCT (44.1%) and maxbin2 (43.3%) binners performed best.

![Impact of metagenomic assembly and MAG binning on recovery of genomic islands. GIs were recovered in a similarly poor fashion to plasmids. Generally, \<40% were correctly binned to the same bin majority commprised of chromosomal contigs from the same source genome regardless of binning (x-axis) and assembly (facet) methods at >50% coverage. metaSPAdes performed the best at assembling GIs (blue). Maxbin2 and CONCOCT placed GIs in a bin majority of the time (orange) however a very small fraction was correctly binned (green). Generally, GIs were correctly binned better than plasmids with DASTool, MetaBAT2 and CONCOCT.](images/GI_recovery.png){#fig:gis width="5in"}

### Recovery of Specific Gene Content

In term of gene content, we first explored the ability to find open reading frames (ORFs) within MAGs.
Overall, the total number of predicted ORFs in MAGs followed a similar trend (Fig @fig:geneContent) as the chromosomal coverage (Fig. @fig:chromcover) and purity (Fig. @fig:purity).
Of the four binning tools, CONCOCT performed the worst, finding <30% of the number of ORFs in our reference genomes.
MetaBAT2 performed second worst at ~80%.
DASTool recovered a similar number to our reference and Maxbin2 seemed to predicted 7-46% more genes.
The Assembler method did not significantly impact the number of genes predicted with the exception of Maxbin2, in which IDBA_UD was the closest to reference and metaSPAdes predicted 46% more ORFs.

![Predicted Gene Content. The total number of open reading frames (ORF) predicted followed the same trend as chromosomal coverage and purity. The assemblers (colored bars) did not contribute to a big variance in the number of ORFs. Of the 4 binners, CONCOCT recovered \<30\% of our reference genome ORFs. DASTool and MetaBAT2 predicted a similar number as our reference genomes.](images/number_of_predicted_genes.png){#fig:geneContent width="8in"}

#### AMR Genes

First, we focused on the ability of MAGs to recover clinically relevant AMR genes (Fig. @fig:AMRGenePercentRecoveryStage).
After the assembly stage, we were only able to recover between ~49-55% of the AMR genes predicted in our reference genomes regardless of the assembly tool used, with metaSPAdes performing marginally better than other assemblers.
Binning the contigs resulted in a ~1-15% loss in AMR gene recovery with the CONCOCT-metaSPAdes pair performing the best at only 1% loss and DASTool-megahit performing the worst at 15% reduction of AMR genes recovered.
Overall, only 24% - 40% of all AMR genes were correctly binned.
This was lowest with the maxbin2-IDBA-UDA pair (24%) and highest in the CONCOCT-metaSPAdes pipe (40%).

![Percent of reference antimicrobial resistance genes (AMR) recovered across assemblers and binners. The proportion of reference AMR genes recovered (y-axis) was largely similar across assembly tools (blue), at roughly 50% with metaSPAdes performing marginally better. Binning tools resulted in a small reduction in AMR genes recovered (orange), however only 24-40% of all AMR genes were correctly binned (green). metaSPAdes-CONCOCT was the best performing MAG binning pipeline. ](images/amr_recovery.png){#fig:AMRGenePercentRecoveryStage width="15in"}

Moreover, focusing on only the AMR genes that were correctly binned (Fig. @fig:AMRGenePercentRecoveryCorrectlyBinned) we can evaluate the impact of different genomic contexts (i.e. chromosomal, plasmid, GI).
Across all methods only 30%-53% of all chromosomally located AMR genes (n=120), 0-45% of genomic island located AMR genes (n=11) and none of the plasmid located AMR genes (n=20) were correctly binned.

![Percent of correctly binned AMR genes recovered by genomic context. MAG methods were best at recovering chromosomally located AMR genes (light blue) regardless of metagenomic assembler or binning tool used. Recovery of AMR genes in GIs showed a bigger variation between tools (light green). None of the 12 evaluated MAG recovery methods were able to recover plasmid located AMR genes (orange).](images/amr_localization_recovery.png){#fig:AMRGenePercentRecoveryCorrectlyBinned width="15in"}

#### VF Genes

Aside from AMR genes, we also examined the impact of MAG approaches on recovery of virulence factor (VF) genes as identified using the Virulence Factor Database (VFDB).
We saw a similar trend as AMR genes (Fig. @fig:VFGenePercentRecoveryStage).
Between 56% and 64% of VFs were identifiable in the metagenomic assemblies (with megahit recovering the greatest proportion).
The binning process further reduced the number of recovered VFs by 4-26% with DASTool-megahit performing the worst (26%) and CONCOCT-metaSPAdes performing the best (4%).
Unlike AMR genes, the majority of VF genes assigned to a bin were assigned to the correct bin (i.e. that bin largely made up of contigs from the same input genome).
Overall, CONCOCT-metaSPAdes again performed best with 43% of all VFs correctly assigned.

![Percent of reference virulence factor (VF) genes recovered across assemblers and binners. The proportion of reference VF genes recovered (y-axis) exhibited a similar trend as AMR genes. Recovery was greatest after the assembling stage (blue), with megahit performing best. Binning tools resulted in a larger reduction in VF genes recovered (orange) compared to AMR genes. However, in majority of cases, VF genes that are binned are correctly binned (green). metaSPAdes-CONCOCT was again the best performing pair.](images/vf_recovery.png){#fig:VFGenePercentRecoveryStage width="15in"}

Again, the genomic context (chromosome, plasmid, GI) of a given VFs largely determined how well it was binned (Fig. @fig:VFGenePercentRecoveryStage).
The majority (73%-98%) of all chromosomally located VF genes (n=757) were correctly binned.
However, 0-16% of GI located VF genes (n=809) and again none of the plasmid located VF genes (n=3) were recovered across all 12 MAG pipelines.

![Percent of correctly binned VF genes Recovered in each genomic region. Metagenome assembled genomes (MAGs) were again best at recovering chromosomally located VF genes (light blue), able to correctly bin majority of chromosomally located VFs. GIs recovered again performed very poorly (light green) and again none of the plasmid located AMR genes (orange) was correctly binned.](images/vf_localization_recovery.png){#fig:VFGenePercentRecoveryStage width="15in"}


## Discussion {#discussion}

In this paper, we evaluated the ability and accuracy of metagenome-assembled genomes (MAGs) to correctly recover mobile genetic elements (i.e. genomic islands and plasmids) from metagenomic samples across different tools used to assemble and bin MAGs.

Overall, the best assembler-binner pair was megahit-DASTOOL in term of both chromosomal coverage (94.3%) and bin purity (1).
Looking at genomes with the lowest coverage, the three Streptococcus genomes that was recovered poorly are likely due to their similarity.
This supports the intuition that MAG recovery approaches struggle to distinguish closely related species.
While CONCOCT performed significantly worse than other binners in terms of chromosomal coverage and bin purity, we did notice that CONCOCT seems to display a trend of generating many small partial bins.
Potentially, CONCOCT binning could be used to distinguish closely related species but at a cost of more fragmented genomes.
 
While the overall recovery and binning of chromosomes were acceptable, we were specifically interested in the ability of MAG methods to appropriately recover MGEs. 
This was due to the importance of MGEs in the function and spread of pathogen traits such as AMR and virulence, as well our hypothesis that their sequence characteristics (composition and copy number) would prove difficult to bin.
Unfortunately, our hypothesis was confirmed, despite the metagenomic assembly approach or MAG binning method used both plasmids and GIs were disproportionately lost compared to chromosomes in general.
At best (via metaSPAdes and CONCOCT) 29.2% of plasmids and 44.1% GIs were identifiable at >50% coverage in the correct bin (i.e. grouped with a bin that was mostly made up of contigs from the same genome).
The >50% coverage requirement might have been a high-bar and there is a possibility that more GIs and plasmids were recovered but in very incomplete forms.

This poor result is not unexpected as genomic islands and plasmids have known divergent compositional features and are often repetitive with variable copy numbers relative to the chromosomes.
Furthermore, the difference between the percentages suggest that binning plasmids are harder than binning GIs.
This difference is partially attributal to the known difficulties in assembly of plasmids from short-read data [@1EFAqjrRj].
Therefore, binning efficiency might improve if we use DNA sequencing and assembly methods optimised for recovering plasmids [@12zFifp5x].

Due to the importance of mobile genetic elements to disseminate clinically relevant AMR genes and VFs, we explored whether or not MAG approaches can be used to provide useful insight LGT of these genes.
With respect to AMR genes, MAG methods were able to recover roughly 40% of all AMR genes present in our reference genomes.
We noted a sharp drop in the number of AMR genes detected between assemblies and MAGs, suggesting that many of these genes were left in the unbinned portion.
Overall, CONCOCT-metaSPAdes combination, while it did not recover the highest amount of AMR genes at the assembly stage, performed the best in correctly binning an AMR gene to the right species.
Regardless of tools, chromosomally-located AMR genes were most frequently correctly binned (as expected from the relative performance of MAGs at recovering chromosomes).
While there was a lot of variability in performance, AMR genes located on GIs were correctly binned slightly less well than chromosomally located AMR genes.
This variability might be explained by the fact there were only 11 GI located AMR genes in our reference genomes.
All 20 of the plasmid-borne AMR genes were assembled but none were placed into MAG bins.
This included high-threat MGEs-based AMR genes such as the KPC and OXA carbapenemases.

Virulence factors had shown a similar trend as AMR genes, recovering ~63% of virulence factors present in the reference genome.
There still is a sharp decline in the number of VF detected between assemblies and MAGs and CONCOCT-metaSPAdes again produced the highest binning accuracy.
MAGs were also able to correctly bin majority (73%-98%) of chromosomally located VF genes to the right species.
However, MAGs performed much worse in correctly recovering GI located and plasmid located VFs, with <16% of GI VFs (n=809) correctly recovered and none of the plasmid VFs (n=3).
This drastic reduction in recovery accuracy of mobile elements, especially GI, is expected.
Previous studies has found that VFs are disproportionally present on GIs[@LxGqo7iq], which might be the reason to why the recovery accuracy was worse compared to AMR genes.
Together, this and the AMR gene results suggests that MAG-based methods might be of limited utility in public health research focused on the transmission and dissemination of AMR genes and VFs. 

It should also be noted that while CONCOCT performed the best in terms of recovery of both chromosomes and MGEs, it created lots of relatively clean but fragmentary partial MAGs.
While this might be ideal for some users, caution should be taken in using CONCOCT when assuming a bin represents a whole genome.

With the recovery of plasmids, GIs, VFs, and AMR genes the same pattern was observed, a progressive loss of data in each analytical step.
The act of metagenomic assembly itself generally resulted in the loss of the majority of these elements/genes regardless of assembly method used.
Across binning tools, the binning process resulted in further loss with a large proportion of MGEs and genes left unbinned.
Finally, only a very small proportion of these elements/genes were generally correctly binned with the appropriate host chromosomes.
While the concept that analysis is "lossy" and that the more analysis you do the more of the input data you are likely to lose is fairly well known it is rarely explicitly stated.
Indeed, this is one of the reasons why the huge amount of redundancy in metagenomic sequencing is necessary i.e. many more base-pairs of DNA than are in the underlying sample.


## Conclusions {#conclusions}

Using a simulated medium complexity metagenome, this study had shown that MAG based approaches provides a useful tool to study a bacterial species’ chromosomal elements but have severe limitations in the recovery of MGEs.
The majority of these mobile genetic elements will both fail to assemble or be correctly binned.
The consequence of this is the disproportionate loss of key public health priority genes like VF and AMR genes.
This is particularly acute as the VF and AMR genes found on these poorly recovered MGEs are generally considered as the most important due to their propensity for lateral gene transfer between unrelated bacteria.
Therefore, it is vital that we utilize a combination of MAGs and other methods (e.g. read-based methods) in public health metagenomic research.
Without this, MAG-based methods are insufficient to thoroughly profile the resistome and provide vital epidemiological data for metagenomic data.


## Supplementals {#supplementals}


![Top Species Coverage](images/s1_species_top_coverage.png){#fig:supspeciescov}


We looked at the ability for MAGs to predict subcellular localization of proteins using PSORTb. Overall, the localization distribution of predicted proteins were very similar in MAGs compared to the reference genome (Fig. (@fig:12subcellularLocalization)). Previous works have shown that AMR genes that are on mobile genetic elements disproportionally encode secrete proteins. Given that the recovery of plasmid-borne genes were not great, we asked if MAGs would affect the ability to predict the subcellular localization of proteins. We found that the proportion of predicted localizations were very similar between MAGs and our reference genomes, suggesting that there is not a significant penalty to use MAGs as input for protein localization predictions.

![Distribution of Predicted Protein Subcellular Localization](images/12subcellularLocalization.png){#fig:12subcellularLocalization width="15in"} 
