###################################################################################################################
#####  Checking the presence of human population SNPs identified in 1000 genome project in COAD/STAD samples: #####
#####  providing evidence for de novo signature 8 representing SNP contamination.                             #####
##### N. Volkova, EMBL-EBI, 2017                                                                              #####
###################################################################################################################

load('ICGC data/profiles_and_decomposition.RData')
load('ICGC data/COAD/vcf_list_COAD.RData')
load('ICGC data/STAD/vcf_list_STAD.RData')

## SNP signature
# Doublecheck the SNPs in the samples

# Get the SNPs
library('SNPlocs.Hsapiens.dbSNP144.GRCh37')
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# Find substitutions in COAD variants
regs <- sapply(vcf_list_COAD, function(vcf) {
  vcf[width(vcf$REF)==1 & width(unlist(vcf$ALT))==1,]
})
# Adjust chromosome names
for (sample in names(regs)) {
  seqlevels(regs[[sample]]) <- substr(seqlevels(regs[[sample]]),4,nchar(seqlevels(regs[[sample]])))
}
# Intersect SNPs and substitutions
snps_in_regs <- list()
for (i in seq_along(vcf_list_COAD)) {
  genome(regs[[i]]) <- "GRCh37.p13"
  snps_in_regs[[i]] <- snpsByOverlaps(snps,regs[[i]],minoverlap=1)
  inds <- match(pos(snps_in_regs[[i]]),start(regs[[i]]))
  alts <- as.character(unlist(regs[[i]][inds]$ALT))
  snpalts <- as.character(IUPAC_CODE_MAP[snps_in_regs[[i]]$alleles_as_ambig])
  if (length(alts)==0) next
  snps_in_regs[[i]] <- snps_in_regs[[i]][sapply(1:length(inds), function(j) {
    grepl(alts[j],snpalts[j])
  })]
}
names(snps_in_regs) <- names(regs)
# Find substitutions in STAD variants
STregs <- sapply(vcf_list_STAD, function(vcf) {
  vcf[width(vcf$REF)==1 & width(unlist(vcf$ALT))==1,]
})
# Intersect SNPs and substitutions
STsnps_in_regs <- list()
for (i in seq_along(vcf_list_STAD)) {
  seqlevels(STregs[[i]]) <- substr(seqlevels(STregs[[i]]),4,nchar(seqlevels(STregs[[i]])))
  genome(STregs[[i]]) <-  "GRCh37.p13"
  STsnps_in_regs[[i]] <- snpsByOverlaps(snps,STregs[[i]],minoverlap=1)
  print(i)
  inds <- match(pos(STsnps_in_regs[[i]]),start(STregs[[i]]))
  alts <- as.character(unlist(STregs[[i]][inds]$ALT))
  snpalts <- as.character(IUPAC_CODE_MAP[STsnps_in_regs[[i]]$alleles_as_ambig])
  if (length(alts)==0) next
  STsnps_in_regs[[i]] <- STsnps_in_regs[[i]][sapply(1:length(inds), function(j) {
    grepl(alts[j],snpalts[j])
  })]
}
names(STsnps_in_regs) <- names(STregs)
# plot the fraction of potential SNPs to all variants per sample
hist( sapply(STsnps_in_regs,length) / sapply(STregs,length), breaks=20)

# Get coding and non-coding SNPs from 1000 Genomes project (http://www.internationalgenome.org/category/vcf/)
coding_vcf <- readVcf("ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf")
non_coding_vcf <- readVcf("ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz")

# Count the coding SNPs
regs <- c(regs, STregs)
cod_snps_by_overlaps <- sapply(regs, function(reg) {
  tmp <- findOverlaps(granges(coding_vcf), reg)
  alts <- as.character(unlist(reg[subjectHits(tmp)]$ALT))
  snpalts <- as.character(unlist(granges(coding_vcf)$ALT))[queryHits(tmp)]
  tmp[alts==snpalts]
})
#non_cod_snps_by_overlaps <- sapply(regs, function(reg) findOverlaps(non_coding_vcf, reg))
counts <- sapply(cod_snps_by_overlaps, function(x) 
  return(c(length(grep(":synonymous:",unlist(info(coding_vcf[queryHits(x)])$VA))),
           length(grep(":nonsynonymous:",unlist(info(coding_vcf[queryHits(x)])$VA))))))
# Check nonsysnonymous to synonimous
length(grep(":nonsynonymous:",unlist(info(coding_vcf)$VA))) / length(grep(":synonymous:",unlist(info(coding_vcf)$VA)))

# Take the ones with signature 8
suspects <- rownames(decomposition)[decomposition[,8]>0.3]
barplot(counts[2,suspects[counts[1,suspects]!=0]] / (counts[1,suspects[counts[1,suspects]!=0]]))
# Looks like it is SNP contamination!