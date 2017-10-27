##########################################################################################
######### Homopolymer analysis for C. elegans and for human exome data          ##########
######### M. Gerstung, B. Meier, N. Volkova, EMBL-EBI, Univ. of Dundee, 2016-17 ##########
##########################################################################################

# Libraries
library(Biostrings)
library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(ggplot2)

# Human genomes
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome.Hsapiens.UCSC.hg19)

# Human exome
#exome <- read.table("S04380110_Covered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, skip=2,quote="")
a <- getSeq(get(ref_genome))
a <- a[1:23] # get rid of Y!!!
gr <- as(seqinfo(a), "GRanges") # turn into GRanges object
genome(gr) <- "hg19"
#exactexomelist <- list()
#for (j in 1:23) {
#  tmp <- exome[exome$V1==seqlevels(a)[j],]
#  exactexomelist[[j]] <- lapply(1:nrow(tmp), function(i)
#    a[[j]][tmp$V2[i]:tmp$V3[i]])
#}
load('exome.RData') # contains pre-calculated coordinates of exome

############################################################################################################
# C. elegans
# define homopolymers
motif_6bp_HP <- DNAStringSet(c("BAAAAAAB", "VTTTTTTV", "HGGGGGGH", "DCCCCCCD"))
homopolymer_pool_length_4to35 <- lapply(seq(4,35),function(y) lapply(c("A","C","G","T"),function(x)
  if (x=="A") {paste0("B",paste(rep(x,y),collapse=""),"B")
  }
  else
    if (x=="T") {paste0("V",paste(rep(x,y),collapse=""),"V")
    }
  else
    if (x=="G") {paste0("H",paste(rep(x,y),collapse=""),"H")
    }
  else
    if (x=="C") {paste0("D",paste(rep(x,y),collapse=""),"D")
    }
))
whomopolymers <- unlist(homopolymer_pool_length_4to35)
list(whomopolymers)
# turn into DNA stringset #
pattern_all <- DNAStringSet(whomopolymers)

# Find homopolymers in each chromosome
whits_all <- lapply(WBcel235, function(chr) sapply(pattern_all, matchPattern, subject=chr, fixed=FALSE))
wlengths <- lapply(whits_all, function(chr) sapply(chr, length)) # numbers of homopolymers on each chromosome
wnonzero_hits_all <- lapply(1:7, function(i) whits_all[[i]][which(wlengths[[i]]>0)]) # 89 82 85 87 13 92 89 - numbers of classes for each chromosome
wnonzero_lengths <- lapply(wnonzero_hits_all, function(chr) sapply(chr, length)) # numbers of homopolymers on each chromosome for non empty classes only
names(wnonzero_hits_all) <- c("I","II","III","IV","MtDNA","V","X")

# Generate GRanges object out of list of all motif hits
sites.gr.worm <- lapply(1:7, function(i) do.call("c",lapply(lapply(wnonzero_hits_all[[i]], as, "IRanges"),
                                                            GRanges,seqnames=names(wnonzero_hits_all)[i])))

# Add corresponding motif sequences and genomic homopolymer length to the GRanges object
for (j in 1:7) {
  sites.gr.worm[[j]]$pattern.searched <- rep(as.character(pattern_all), sapply(whits_all[[j]], length))
  sites.gr.worm[[j]]$motif.found <- unlist(lapply(wnonzero_hits_all[[j]], as.character))
  sites.gr.worm[[j]]$pattern.length <- unlist(lapply(wnonzero_hits_all[[j]], width)) # includes the two flanking bases
  sites.gr.worm[[j]]$homopolymer.length <- (sites.gr.worm[[j]]$pattern.length)-2 # remove 2 bases, 5' and 3' are not part of homopolymer
  #sites.gr.worm[[j]] # granges object
  #elementMetadata(sites.gr[[j]]) # data frame
}

# combine information from all Chromosomes into one GRanges object showing all hits across the genome
sites.gr_all.worm <- do.call("c",sites.gr.worm[-5]) # exclude MtDNA
genome(sites.gr_all.worm) <- "WBcel235" # add genome info
##########################################################################################
# Human exome
homopolymer_pool_length_4to55 <- lapply(seq(4,55),function(y) lapply(c("A","C","G","T"),function(x)
  if (x=="A") {paste0("B",paste(rep(x,y),collapse=""),"B")
  }
  else
    if (x=="T") {paste0("V",paste(rep(x,y),collapse=""),"V")
    }
  else
    if (x=="G") {paste0("H",paste(rep(x,y),collapse=""),"H")
    }
  else
    if (x=="C") {paste0("D",paste(rep(x,y),collapse=""),"D")
    }
))
homopolymers <- unlist(homopolymer_pool_length_4to55)
pattern_all <- DNAStringSet(homopolymers)

# Find all HPs in the genome
hits_all <- list()
for (i in 1:23) {
  maskMotif(a[[i]], "N") -> masked
  hits_all[[i]] <- sapply(pattern_all, 
                          matchPattern, 
                          subject=masked,
                          fixed=F)
  print(i)
}
lengths <- lapply(hits_all, function(l) sapply(l,length))
hits_all <- lapply(1:length(hits_all), function(i) hits_all[[i]][which(lengths[[i]]>0)])
nonzero_lengths <- lapply(hits_all, function(chr) sapply(chr, length)) # numbers of homopolymers on each chromosome for non empty classes only
names(hits_all) <- seqlevels(gr)
# Get the context and coordinates of all Hps
sites.gr <- lapply(1:23, function(i) do.call("c",lapply(lapply(hits_all[[i]], as, "IRanges"),
                                                        GRanges,seqnames=names(hits_all)[i])))
# Add corresponding motif sequences and genomic homopolymer length to the GRanges object
for (j in 1:23) {
  sites.gr[[j]]$pattern.searched <- rep(as.character(pattern_all), lengths[[j]])
  sites.gr[[j]]$motif.found <- unlist(lapply(hits_all[[j]], as.character))
  sites.gr[[j]]$pattern.length <- unlist(lapply(hits_all[[j]], width)) # includes the two flanking bases
  sites.gr[[j]]$homopolymer.length <- (sites.gr[[j]]$pattern.length)-2 # remove 2 bases, 5' and 3' are not part of homopolymer
}
sites.gr_all <- do.call("c",sites.gr)
genome(sites.gr_all) <- "hg19" # add genome info
# Intersect with exonic regions
exomranges <- GRanges(seqnames=exome$V1,ranges=IRanges(start = exome$V2,end=exome$V3))
hits <- findOverlaps(exomranges,sites.gr_all,minoverlap = 2)
sites.gr_all <- sites.gr_all[sort(unique(subjectHits(hits)))]
all_sites <- as.data.frame(sites.gr_all) # make a data frame

# Create a dataframe with homopolymer information
homopolymers.worm <- as.data.frame(table(all_sites_worm$pattern.searched))
homopolymers.worm$length <- width(as.character(homopolymers.worm$Var1))-2 # include length of motif, -2 removes flanking 5'and 3' bases
homopolymers.worm$base <- ifelse(grepl("AAA",homopolymers.worm$Var1),"A",ifelse(grepl("TTT", homopolymers.worm$Var1), "T", ifelse(grepl("GGG", homopolymers.worm$Var1), "G", "C")))
homopolymer_frequency_by_length_worm <- homopolymers.worm %>%
  group_by(length) %>%
  summarise(total_frequency =sum(Freq))
homopolymer_frequency_by_length_worm$total_number_of_bases <- (homopolymer_frequency_by_length_worm$length)*as.numeric(homopolymer_frequency_by_length_worm$total_frequency) # calculate the number of bases for each HP
##########################################################################################
# C. elegans indels in homopolymers
# Upload vcf from C. elegans with indels in a list of vcfs called MMRvcf

indel1bp <- sapply(MMRvcf, function(vcf) vcf[abs(width(granges(vcf)$REF)-width(unlist(granges(vcf)$ALT)))==1 &
                                               (width(granges(vcf)$REF)==1 | width(unlist(granges(vcf)$ALT))==1)])
barplot(sapply(indel1bp, length))

indel1bp <- sapply(indel1bp, rowRanges)
for (k in 1:length(indel1bp)) {
  indel1bp[[k]]$insertion <- width(unlist(indel1bp[[k]]$ALT))-1
  indel1bp[[k]]$deletion <- width(indel1bp[[k]]$REF)-1
}

# Now check the intersections of homoplymers and indels;
# watch out for cases when an indel is between 2 homopolymers and can be counted twice
worm_hit_list <- sapply(indel1bp, function(x) findOverlaps(x,sites.gr_all.worm))
indels.to.check.list <- sapply(worm_hit_list, function(x) unique(as.matrix(x)[,1])[which(table(as.matrix(x)[,1])>1)])
for (t in which(sapply(indels.to.check.list,length)>0)) {
  to.delete.final <- NULL
  indels.to.check <- indels.to.check.list[[t]]
  for (j in 1:length(indels.to.check)) {
    homopolymers.to.check <- as.matrix(worm_hit_list[[t]])[which(as.matrix(worm_hit_list[[t]])[,1] == (indels.to.check[j])),2]
    
    homopolymers.types <- sapply(sites.gr_all.worm[homopolymers.to.check]$pattern.searched,function(x) substr(x,2,2))
    
    ref <- unlist(strsplit(as.data.frame(indel1bp[[t]][indels.to.check[j]])$REF,split = ""))
    alt <- unlist(strsplit(as.character(as.data.frame(indel1bp[[t]][indels.to.check[j]])$ALT[[1]]),split = ""))
    
    if (indel1bp[[t]][indels.to.check[j]]$insertion==1) {
      indel.type <- setdiff(alt,ref)
      to.delete <- which(as.matrix(worm_hit_list[[t]])[,1]==indels.to.check[j]) [homopolymers.types != indel.type]
      if (length(to.delete)==0) {
        if (indel.type==alt[1]) to.delete <- which(as.matrix(worm_hit_list[[t]])[,1]==indels.to.check[j])[2]
        if (indel.type==alt[2]) to.delete <- which(as.matrix(worm_hit_list[[t]])[,1]==indels.to.check[j])[1]
      }
    }
    if (indel1bp[[t]][indels.to.check[j]]$deletion==1) {
      indel.type <- setdiff(ref,alt)
      to.delete <- which(as.matrix(worm_hit_list[[t]])[,1]==indels.to.check[j]) [homopolymers.types != indel.type]
      if (length(to.delete)==0) {
        if (indel.type==ref[1]) to.delete <- which(as.matrix(worm_hit_list[[t]])[,1]==indels.to.check[j])[2]
        if (indel.type==ref[2]) to.delete <- which(as.matrix(worm_hit_list[[t]])[,1]==indels.to.check[j])[1]
      }
    }
    to.delete.final <- c(to.delete.final,to.delete)
  }
  if (length(to.delete.final)>0) worm_hit_list[[t]] <- worm_hit_list[[t]][-to.delete.final]
}
worm_hit_list <- worm_hit_list[sapply(worm_hit_list,length)>0]
# Turn the hist into a list of dataframes
worm_hits_data_frame <- list()
for (i in 1:length(worm_hit_list)) {
  queries <- unique(queryHits(worm_hit_list[[i]]))
  idx <- sapply(queries, function(x) which(queryHits(worm_hit_list[[i]])==x)[1])
  idxs <- subjectHits(worm_hit_list[[i]])[idx]
  motifs <- DataFrame(sites.gr_all.worm$motif.found[idxs])
  length <- unlist(lapply(sites.gr_all.worm$motif.found[idxs], width))-2
  worm_hits_data_frame[[i]] <- as.data.frame(indel1bp[[names(worm_hit_list)[i]]][queries,])
  worm_hits_data_frame[[i]]$motif.found <- sites.gr_all.worm$motif.found[idxs]
  worm_hits_data_frame[[i]]$motif.length <- width(worm_hits_data_frame[[i]]$motif.found)
  worm_hits_data_frame[[i]]$homopolymer.length <- worm_hits_data_frame[[i]]$motif.length-2
  print(i)
}
names(worm_hits_data_frame) <- names(worm_hit_list)
# Add homopolymer type
for (i in 1:length(worm_hits_data_frame)) {
  worm_hits_data_frame[[i]]$base <- ifelse(grepl("AAA",worm_hits_data_frame[[i]]$motif.found),"A",
                                           ifelse(grepl("TTT", worm_hits_data_frame[[i]]$motif.found), "T", 
                                                  ifelse(grepl("GGG", worm_hits_data_frame[[i]]$motif.found), "G", "C")))
  worm_hits_data_frame[[i]]$hpn <- paste0(worm_hits_data_frame[[i]]$homopolymer.length,worm_hits_data_frame[[i]]$base)
}

# Group the data
mlh1_20 = do.call("rbind",worm_hits_data_frame[c("CD0134a","CD0134c","CD0134d")])
pms2_20 = do.call("rbind",worm_hits_data_frame[c("CD0135a","CD0135c","CD0135d")])
pms2_10 = do.call("rbind",worm_hits_data_frame[c("CD0244a","CD0244c","CD0244d")])
pole4pms2_10 = do.call("rbind",worm_hits_data_frame[c("CD0246d","CD0246e")])
# Variant calling in long homopymers is infeasible (too few reads covering the whole homopolymer)
mlh1_20 = mlh1_20[mlh1_20$homopolymer.length<19,]
pms2_20 = pms2_20[pms2_20$homopolymer.length<19,]
pms2_10 = pms2_10[pms2_10$homopolymer.length<19,]
pole4pms2_10 = pole4pms2_10[pole4pms2_10$homopolymer.length<19,]

# Plot raw frequencies and ratios

total_frequency = homopolymer_frequency_by_length_worm$total_frequency[1:15]
to.show <- data.frame(homopolymer.length = 4:(length(total_frequency)+3), 
                      pms2.20 = vector("numeric",length(total_frequency)),
                      mlh1.20 = vector("numeric",length(total_frequency)),
                      pms2.10 = vector("numeric",length(total_frequency)),
                      pole4pms2.10 = vector("numeric",length(total_frequency)))
to.show$mlh1.20[as.numeric(names(table(mlh1_20$homopolymer.length)))-3] <- (table(mlh1_20$homopolymer.length)/3)
to.show$pms2.20[as.numeric(names(table(pms2_20$homopolymer.length)))-3] <- (table(pms2_20$homopolymer.length)/3) 
to.show$pms2.10[as.numeric(names(table(pms2_10$homopolymer.length)))-3] <- (table(pms2_10$homopolymer.length)/3)
to.show$pole4pms2.10[as.numeric(names(table(pole4pms2_10$homopolymer.length)))-3] <- round(table(pole4pms2_10$homopolymer.length)/2)


library(ggplot2)
library(reshape2)
df <- melt(to.show,id.vars = "homopolymer.length")
#df$value <- log10(df$value)
#df = df[-which(df$value==-Inf),]
pdf("~/Documents/Frequency of 1 bp indels in homopolymers FILT.pdf",10,3)
ggplot(data=df,aes(x=homopolymer.length,y=value)) + geom_bar(stat='identity', width = 0.75) +
  facet_grid(. ~ factor(variable,levels=c("mlh1.20","pms2.20","pms2.10","pole4pms2.10"))) +theme_bw() +
  ylab("Number of 1 bp indels") + xlab("Homopolymer length") + 
  theme(axis.text = element_text(size=12, family='ArialMT'), 
        axis.title = element_text(size=16, family='ArialMT'), 
        strip.text = element_text(size=16, family='ArialMT'),
        panel.grid = element_blank(), strip.background = element_blank())
dev.off()

for (i in 2:ncol(to.show)) to.show[,i] = to.show[,i] / total_frequency
df <- melt(to.show,id.vars = "homopolymer.length")
ggplot(data=df,aes(x=homopolymer.length,y=value)) + 
  facet_grid(. ~ factor(variable,levels=c("mlh1.20","pms2.20","pms2.10","pole4pms2.10"))) +
  geom_point() + geom_line() + theme_bw() + 
  ylab("Ratio \n Mutations/HP") + xlab("Homopolymer length") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), strip.text = element_text(size=16))


# Apply generalized additive models with splines

total_frequency = homopolymer_frequency_by_length_worm$total_frequency[1:15]
to.show <- data.frame(homopolymer.length = 4:(length(total_frequency)+3), 
                      pms2 = vector("numeric",length(total_frequency)),
                      mlh1= vector("numeric",length(total_frequency)),
                      pole4pms2= vector("numeric",length(total_frequency)),
                      pms210 = vector("numeric",length(total_frequency)),
                      total_frequency = total_frequency,
                      total_number_of_bases = homopolymer_frequency_by_length_worm$total_number_of_bases[1:15])
to.show$mlh1[as.numeric(names(table(mlh1_20$homopolymer.length)))-3] <- round(table(mlh1_20$homopolymer.length)/3)
to.show$pms2[as.numeric(names(table(pms2_20$homopolymer.length)))-3] <- round(table(pms2_20$homopolymer.length)/3)
to.show$pms210[as.numeric(names(table(pms2_10$homopolymer.length)))-3] <- round(table(pms2_10$homopolymer.length)/3)
to.show$pole4pms2[as.numeric(names(table(pole4pms2_10$homopolymer.length)))-3] <- round(table(pole4pms2_10$homopolymer.length)/2)


library(gam)
par(mfrow = c(1,4))
fit <- gam(cbind(mlh1, total_frequency-mlh1) ~ s(homopolymer.length, 6), data=to.show, family=binomial)
plot(to.show$homopolymer.length, to.show$mlh1/to.show$total_frequency, xlim=c(4,18), ylim=c(0,0.015), ylab="Ratio \n Mutations/HP", pch=16, col='grey', main='mlh-1',font.main=3, xlab="Homopolymer length")
segments(to.show$homopolymer.length, qbeta(0.025, .5+to.show$mlh1, .5+to.show$total_frequency- to.show$mlh1),
         to.show$homopolymer.length, qbeta(0.975, .5+to.show$mlh1, .5+to.show$total_frequency- to.show$mlh1), col='grey')
p <-  predict(fit, type='link', se.fit=TRUE)
lines(to.show$homopolymer.length, exp(p$fit), col='red', lwd=2)
lines(to.show$homopolymer.length, exp( p$fit +  2*p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit -  2* p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit), col='#FF00000d', lwd=2)
for(q in seq(0.025, 0.475,0.025))
  polygon(c(to.show$homopolymer.length, rev(to.show$homopolymer.length)), c(exp(qnorm(q, p$fit, p$se.fit)), exp(rev(qnorm(1-q, p$fit, p$se.fit)))), border=NA, col='#FF00000d')

fit <- gam(cbind(pms2, total_frequency-pms2) ~ s(homopolymer.length, 6), data=to.show, family=binomial)
plot(to.show$homopolymer.length, to.show$pms2/to.show$total_frequency, yaxt="n",
     xlim=c(4,18), ylim=c(0,0.015), ylab=NA,
     pch=16, col='grey', main='pms-2',font.main=3, xlab="Homopolymer length")
segments(to.show$homopolymer.length, qbeta(0.025, .5+to.show$pms2, .5+to.show$total_frequency- to.show$pms2),
         to.show$homopolymer.length, qbeta(0.975, .5+to.show$pms2, .5+to.show$total_frequency- to.show$pms2), col='grey')
p <-  predict(fit, type='link', se.fit=TRUE)
lines(to.show$homopolymer.length, exp(p$fit), col='red', lwd=2)
lines(to.show$homopolymer.length, exp( p$fit +  2*p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit -  2* p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit), col='#FF00000d', lwd=2)
for(q in seq(0.025, 0.475,0.025))
  polygon(c(to.show$homopolymer.length, rev(to.show$homopolymer.length)), c(exp(qnorm(q, p$fit, p$se.fit)), exp(rev(qnorm(1-q, p$fit, p$se.fit)))), border=NA, col='#FF00000d')

fit <- gam(cbind(pms210, total_frequency-pms210) ~ s(homopolymer.length, 6), data=to.show, family=binomial)
plot(to.show$homopolymer.length, to.show$pms210/to.show$total_frequency, yaxt="n", xlim=c(4,18), ylim=c(0,0.015), ylab=NA, pch=16, col='grey', main='pms-2 10',font.main=3, xlab="Homopolymer length")
segments(to.show$homopolymer.length, qbeta(0.025, .5+to.show$pms210, .5+to.show$total_frequency- to.show$pms210),
         to.show$homopolymer.length, qbeta(0.975, .5+to.show$pms210, .5+to.show$total_frequency- to.show$pms210), col='grey')
p <-  predict(fit, type='link', se.fit=TRUE)
lines(to.show$homopolymer.length, exp(p$fit), col='red', lwd=2)
lines(to.show$homopolymer.length, exp( p$fit +  2*p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit -  2* p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit), col='#FF00000d', lwd=2)
for(q in seq(0.025, 0.475,0.025))
  polygon(c(to.show$homopolymer.length, rev(to.show$homopolymer.length)), c(exp(qnorm(q, p$fit, p$se.fit)), exp(rev(qnorm(1-q, p$fit, p$se.fit)))), border=NA, col='#FF00000d')

fit <- gam(cbind(pole4pms2, total_frequency-pole4pms2) ~ s(homopolymer.length, 6), data=to.show, family=binomial)
plot(to.show$homopolymer.length, to.show$pole4pms2/to.show$total_frequency, yaxt="n", xlim=c(4,18), ylim=c(0,0.015), ylab=NA, pch=16, col='grey', main='pole-4; pms-2',font.main=3, xlab="Homopolymer length")
segments(to.show$homopolymer.length, qbeta(0.025, .5+to.show$pole4pms2, .5+to.show$total_frequency- to.show$pole4pms2),
         to.show$homopolymer.length, qbeta(0.975, .5+to.show$pole4pms2, .5+to.show$total_frequency- to.show$pole4pms2), col='grey')
p <-  predict(fit, type='link', se.fit=TRUE)
lines(to.show$homopolymer.length, exp(p$fit), col='red', lwd=2)
lines(to.show$homopolymer.length, exp( p$fit +  2*p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit -  2* p$se.fit), col='red', lty=2)
lines(to.show$homopolymer.length, exp(p$fit), col='#FF00000d', lwd=2)
for(q in seq(0.025, 0.475,0.025))
  polygon(c(to.show$homopolymer.length, rev(to.show$homopolymer.length)), c(exp(qnorm(q, p$fit, p$se.fit)), exp(rev(qnorm(1-q, p$fit, p$se.fit)))), border=NA, col='#FF00000d')
##########################################################################################
# Human indels in homopolymers
# Get vcf files for STAD and COAD cohorts as described in MMR_paper_analysis.Rmd

indel_list_STAD <- sapply(vcf_list_STAD, function(vcf) vcf[abs(width(vcf$REF)-width(unlist(vcf$ALT)))==1 & (width(vcf$REF)==1 | width(unlist(vcf$ALT))==1),])
indel_list_COAD <- sapply(vcf_list_COAD, function(vcf) vcf[abs(width(vcf$REF)-width(unlist(vcf$ALT)))==1 & (width(vcf$REF)==1 | width(unlist(vcf$ALT))==1),])
for (k in 1:length(indel_list_STAD)) {
  genome(indel_list_STAD[[k]]) <- "hg19"
  indel_list_STAD[[k]]$insertion <- width(unlist(indel_list_STAD[[k]]$ALT))-1
  indel_list_STAD[[k]]$deletion <- width(indel_list_STAD[[k]]$REF)-1
}
for (k in 1:length(indel_list_COAD)) {
  genome(indel_list_COAD[[k]]) <- "hg19"
  indel_list_COAD[[k]]$insertion <- width(unlist(indel_list_COAD[[k]]$ALT))-1
  indel_list_COAD[[k]]$deletion <- width(indel_list_COAD[[k]]$REF)-1
}

# Check the intersections with homopolymers in exome
indels.in.hp.stad <- vector("numeric",length(indel_list_STAD))
for (i in 1:length(indel_list_STAD)) {
  x <- subsetByOverlaps(indel_list_STAD[[i]],sites.gr_all)
  indels.in.hp.stad[i] <- nrow(as.data.frame(x))
}
names(indels.in.hp.stad) <- names(indel_list_STAD)
indels.in.hp.coad <- vector("numeric",length(indel_list_COAD))
for (i in 1:length(indel_list_COAD)) {
  x <- subsetByOverlaps(indel_list_COAD[[i]],sites.gr_all)
  indels.in.hp.coad[i] <- nrow(as.data.frame(x))
}
names(indels.in.hp.coad) <- names(indel_list_COAD)

# Plot everything
df <- data.frame(sample=rownames(STAD.mutation.counts), Non.HP.indels = rowSums(STAD.mutation.counts[,97:104])-indels.in.hp.stad[rownames(STAD.mutation.counts)], HP.indels = indels.in.hp.stad[rownames(STAD.mutation.counts)])
df = melt(df,id.vars="sample")
p1 <- ggplot(data=df[df$sample %in% mmr.ucsc,], aes(x=factor(sample, levels = df$sample[order(rowSums(STAD.mutation.counts[,97:104]), decreasing = F)]),y=value,fill=variable)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust=0, size=6, family='ArialMT'),
        legend.title=element_blank(),
        axis.title = element_text(size=8, family='ArialMT'),
        legend.text = element_text(size=8,family='ArialMT'),
        title = element_text(size=10, family='ArialMT', face='bold')) +
  scale_fill_discrete(labels=c("Indels not in HP","Indels in HP")) + xlab("MMR deficient sample") + 
  ylab("1bp indel counts") + ggtitle("1 bp indels in homopolymers in STAD samples with MMR deficiency")
p1

df <- data.frame(sample=rownames(COAD.mutation.counts), Non.HP.indels = rowSums(COAD.mutation.counts[,97:104])-indels.in.hp.coad[rownames(COAD.mutation.counts)], HP.indels = indels.in.hp.coad[rownames(COAD.mutation.counts)])
df = melt(df,id.vars="sample")
p2 <- ggplot(data=df[df$sample %in% mmr.ucsc,], aes(x=factor(sample, levels = df$sample[order(rowSums(COAD.mutation.counts[,97:104]), decreasing = F)]),y=value,fill=variable)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust=0, size=7, family='ArialMT'),
        legend.title=element_blank(),
        axis.title = element_text(size=8, family='ArialMT'),
        legend.text = element_text(size=8,family='ArialMT'),
        title = element_text(size=10, family='ArialMT', face='bold')) + 
  scale_fill_discrete(labels=c("Indels not in HP","Indels in HP")) + xlab("MMR deficient sample") + 
  ylab("1bp indel counts") + ggtitle("1 bp indels in homopolymers in COAD samples with MMR deficiency")
p2