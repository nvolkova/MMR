# MMR
Mismatch repair signature analysis in *C. elegans* and ICGC colorectal and stomach cancer data for [Meier, Volkova et al. 2017](https://www.biorxiv.org/content/biorxiv/early/2017/06/13/149153.full.pdf) (subm.).

*C. elegans* data: contains .RData files with count and design matrices for MMR experiment on C. elegans, raw data available at the NCBI Sequence Read Archive [SRA](http://www.ncbi.nlm.nih.gov/sra) under accession number SRP020555. Requires *C. elegans* WB235cel genome, instructions in corresponding README file.

**ICGC data**:
- **COAD**: data from ICGC DCC portal on colorectal adenocarcinoma, contains donor and sample codes, MSI status, processed VCF with all variants per sample, methylation info for MLH1 gene; requires downloading data from ICGC portal, instructions in corresponding README file.
- **STAD**: data from ICGC DCC portal on stomach adenocarcinoma, contains donor and sample codes, MSI status, processed VCF with all variants per sample, methylation info for MLH1 gene; requires downloading data from ICGC portal, instructions in corresponding README file.
- **profiles_and_decomposition.RData** - contains mutation count matrices, signatures, and decomposition over these signatures.
All this data is publicly available from [ICGC](http://dcc.icgc.org) and [TCGA CE](http://genomeportal.stanford.edu/pan-tcga).

**plotting_functions**: contains some codes for signature/profile plotting

**exome.RData** - contains pre-calculated positions of well-covered part of human exome (according to Agilent SureSelect V5 Human All Exon, freely available [here](https://earray.chem.agilent.com/earray/) after registration). 

**Cosine_similarities_explained.R** - simulations for explanation of signature threshold; stability analysis of similarity score between *C. elegans* and COAD/STAD signatures by random drawing from 95% confidence intervals and jack-knife bootstrapping.

**homopolymer_analysis.R** - analysis of homopolymers and 1-bp indels in homopolymers for *C. elegans* and human cancer samples.

**nmSolve.R** - contains function for non-negative matrix factorization.

**Signatures_plots.R** - calculation of adjustment coefficients between human exome and *C. elegans* genome, and plotting of respective contexts and adjusted signatures.

**SNP_signature_analysis.R** - analysis of COAD/STAD signature 8 and the ratio of synonimous to non-synonimous changes in the samples with a prevalence of this signature.

**Signatures.xlsx** - table with *C. elegans* mutation patterns, original and exome-adjusted, with and without 1-bp indel fractions; COAD/STAD signatures, with and without 1-bp indel fractions.

In case of any questions or enquiries, please contact Nadezda Volkova (nvolkova@ebi.ac.uk).

N. Volkova, B. Meier, M. Gerstung. EMBL-EBI, 2017.
