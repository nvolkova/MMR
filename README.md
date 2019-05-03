# MMR
Mismatch repair signature analysis in *C. elegans* and ICGC colorectal and stomach cancer data for [Meier, Volkova et al. 2017](https://www.biorxiv.org/content/biorxiv/early/2017/06/13/149153.full.pdf) (subm.).

***C. elegans* data**: contains .RData files with count and design matrices for MMR experiment on C. elegans. Raw sequencing data are available through the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) ENA Study Accession Numbers ERP000975 and ERP004086 as detailed in Suppl. Table 1 of the paper. Requires *C. elegans* WB235cel genome, instructions in corresponding README file.

**ICGC data**:
- **COAD**: data from ICGC DCC portal on colorectal adenocarcinoma, contains donor and sample codes, MSI status, VCF with all variants per sample; requires downloading data from ICGC portal, instructions in corresponding README file.
- **STAD**: data from ICGC DCC portal on stomach adenocarcinoma, contains donor and sample codes, MSI status, VCF with all variants per sample; requires downloading data from ICGC portal, instructions in corresponding README file.
- **profiles_and_decomposition.RData** - contains mutation count matrices, signatures, and decomposition over these signatures.
All this data is publicly available from [ICGC](http://dcc.icgc.org) and [TCGA CE](http://genomeportal.stanford.edu/pan-tcga).
- requires downloading the meta-VCF with precise description of all the mutations in ICGC data, available [here](https://dcc.icgc.org/releases/current/Summary) under **simple_somatic_mutation.aggregated.vcf.gz**.

**plotting_functions**: contains some codes for signature/profile plotting

**exome.RData** - contains pre-calculated positions of well-covered part of human exome (according to Agilent SureSelect V5 Human All Exon, freely available [here](https://earray.chem.agilent.com/earray/) after registration). 

**Cosine_similarities_explained.R** - simulations for explanation of signature threshold; stability analysis of similarity score between *C. elegans* and COAD/STAD signatures by random drawing from 95% confidence intervals and jack-knife bootstrapping.

**homopolymer_analysis.R** - analysis of homopolymers and 1-bp indels in homopolymers for *C. elegans* and human cancer samples.

**nmSolve.R** - contains function for non-negative matrix factorization.

**Signatures_plots.R** - calculation of adjustment coefficients between human exome and *C. elegans* genome, and plotting of respective contexts and adjusted signatures.

**SNP_signature_analysis.R** - analysis of COAD/STAD signature 8 and the ratio of synonimous to non-synonimous changes in the samples with a prevalence of this signature.

**Signatures.xlsx** - table with *C. elegans* mutation patterns, original and exome-adjusted, with and without 1-bp indel fractions; COAD/STAD signatures, with and without 1-bp indel fractions.

Session info with package versions for running the scripts:

```
R version 3.4.4 (2018-03-15)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
 [1] splines   stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SNPlocs.Hsapiens.dbSNP144.GRCh37_0.99.20 bayesm_3.1-0.1                           MASS_7.3-49                              pROC_1.11.0                             
 [5] BSgenome.Hsapiens.UCSC.hg19_1.4.0        BSgenome_1.44.2                          rtracklayer_1.36.6                       gam_1.15                                
 [9] foreach_1.4.4                            dplyr_0.7.4                              deconstructSigs_1.8.0                    VariantAnnotation_1.22.3                
[13] Rsamtools_1.28.0                         Biostrings_2.44.2                        XVector_0.16.0                           SummarizedExperiment_1.6.5              
[17] DelayedArray_0.3.0                       matrixStats_0.53.1                       GenomicRanges_1.28.6                     GenomeInfoDb_1.12.3                     
[21] IRanges_2.10.5                           S4Vectors_0.14.7                         reshape2_1.4.3                           ggplot2_2.2.1                           
[25] mg14_0.0.5                               devtools_1.13.5                          NMF_0.21.0                               Biobase_2.36.2                          
[29] BiocGenerics_0.22.1                      cluster_2.0.6                            rngtools_1.2.4                           pkgmaker_0.22                           
[33] registry_0.5                             tsne_0.1-3                              

loaded via a namespace (and not attached):
 [1] bit64_0.9-7              assertthat_0.2.0         blob_1.1.1               GenomeInfoDbData_0.99.0  yaml_2.1.18              pillar_1.2.1             RSQLite_2.0              lattice_0.20-35         
 [9] glue_1.2.0               digest_0.6.15            RColorBrewer_1.1-2       colorspace_1.3-2         Matrix_1.2-12            plyr_1.8.4               XML_3.98-1.10            pkgconfig_2.0.1         
[17] biomaRt_2.32.1           zlibbioc_1.22.0          xtable_1.8-2             scales_0.5.0             BiocParallel_1.10.1      tibble_1.4.2             withr_2.1.2              GenomicFeatures_1.28.5  
[25] lazyeval_0.2.1           magrittr_1.5             memoise_1.1.0            doParallel_1.0.11        tools_3.4.4              gridBase_0.4-7           stringr_1.3.0            munsell_0.4.3           
[33] AnnotationDbi_1.38.2     bindrcpp_0.2             compiler_3.4.4           rlang_0.2.0              grid_3.4.4               RCurl_1.95-4.10          iterators_1.0.9          bitops_1.0-6            
[41] gtable_0.2.0             codetools_0.2-15         DBI_0.8                  R6_2.2.2                 GenomicAlignments_1.12.2 knitr_1.20               bit_1.1-12               bindr_0.1.1             
[49] stringi_1.1.7            Rcpp_0.12.16         
```

In case of any questions or enquiries, please contact Nadezda Volkova (nvolkova@ebi.ac.uk).

N. Volkova, B. Meier, M. Gerstung. EMBL-EBI, 2017.
