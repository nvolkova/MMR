# Plotting signatures and trinucleotide contexts

# Making comparative plot, worms vs human
library(deconstructSigs)
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
types <- as.character(cancer_signatures$Trinucleotide) # trinucleotide classes
types.full <- as.character(cancer_signatures$Somatic.Mutation.Type) # substitution types
row.names(cancer_signatures) <- types.full
cancer_signatures = as.matrix(cancer_signatures[,4:33])

worm.trinucleotides.32 <- sapply(unique(types), function(x) {
  return(worm.trinucleotides[x] + worm.trinucleotides[as.character(reverseComplement(DNAString(x)))])
})
names(worm.trinucleotides.32) <- unique(types)
types.order <- c(rep(types[1:16],3),rep(types[49:64],3))
# human.trinucleotides in deconstructSigs

WBcel235 <- readDNAStringSet("Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa")
worm.trinucleotides <- colSums(trinucleotideFrequency(WBcel235)[-5,])
#worm.bases <- colSums(oligonucleotideFrequency(WBcel235,width=1)[-5,])
human.trinucleotides <- as.vector(t(tri.counts.genome)) # / sum(tri.counts.genome))) # counts from "deconstructSigs" package
names(human.trinucleotides) <- row.names(tri.counts.genome)
trinucleotide.freq.factor <- sapply(unique(types), function(x) {
  freq.worm <- worm.trinucleotides[x] + worm.trinucleotides[as.character(reverseComplement(DNAString(x)))]
  return(freq.worm /  human.trinucleotides[x]) # tri.counts.genome is already classified w.r.t. pyrimidine reference
})
human.trinucleotides <- as.vector(t(tri.counts.exome)) # / sum(tri.counts.genome))) # counts from "deconstructSigs" package
names(human.trinucleotides) <- row.names(tri.counts.exome)
trinucleotide.freq.factor.ex <- sapply(unique(types), function(x) {
  freq.worm <- worm.trinucleotides[x] + worm.trinucleotides[as.character(reverseComplement(DNAString(x)))]
  return(freq.worm / human.trinucleotides[x]) # tri.counts.genome is already classified w.r.t. pyrimidine reference
})
names(trinucleotide.freq.factor.ex) = names(trinucleotide.freq.factor) <- unique(types)



df = data.frame(worm.counts=c(worm.trinucleotides.32[types.order], rep(0.15*sum(worm.trinucleotides.32[types.order][33:64]),16)),
                human.counts=c(t(tri.counts.exome)[1,][types.order], rep(0.15*sum(t(tri.counts.exome)[1,][types.order][33:64]),16)),
                Type = rep(c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'VOID'), each = 16),
                Trinucleotide = c(types.order,types.order[1:16]))
df$worm.counts <- df$worm.counts / sum(df$worm.counts[33:64])
df$human.counts = df$human.counts / sum(df$human.counts[33:64])
df$human.counts = -df$human.counts

df <- melt(df, id.vars = c('Type','Trinucleotide'))
colnames(df) <- c('Type', 'Trinucleotide', 'Species', 'Fraction')
levels(df$Species) <- c('C. elegans', 'H. sapiens exome')
rownames(df) = NULL

p <- ggplot(data = df,aes(x = Trinucleotide,y=Fraction,fill=Species)) + 
  geom_bar(stat="identity",colour="black",position = "dodge",size=0.1,width = 0.5) + 
  scale_fill_manual(values = c("darkred","lightblue")) +
  facet_grid(Species ~ Type, scales = "free") +
  scale_x_discrete(labels=types.order) +
  theme_bw() + coord_cartesian() +
  theme(text = element_text(family='ArialMT'),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        strip.text = element_text(size = 20),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=24,face="bold"),
        panel.grid = element_blank(),
        strip.background = element_rect(colour='white', fill='white'),
        panel.border = element_rect(colour='black', size=0.1),
        panel.spacing = unit(0.01,'lines'))
p

# Plotting signatures
source('plotting_functions/plot_sigs.R')
library(xlsx)
learned.sigs <- read.xlsx('Signatures.xlsx', sheetIndex = 1)
learned.sigs.exome <- read.xlsx('Signatures.xlsx', sheetIndex = 3)
defsigs <- read.xlsx('Signatures.xlsx', sheetIndex = 5)
pdf('~/Documents/worm_signatures_humanized_with_indels1.pdf', 14,8)
to.show <- cbind(learned.sigs[2,1:96],learned.sigs[2,1:96],learned.sigs[3,1:96],learned.sigs[3,1:96],learned.sigs[5,1:96],learned.sigs[5,1:96])
colnames(to.show) <- c('mlh-1', 'mlh-1 hum', 'pms-2', 'pms-2 hum', 'pole-4;pms-2', 'pole-4;pms-2 hum')
plot_sig_wb(to.show, ymax = 0.1, colors = c("grey29", "#DE1C14", "green3", "blue", "#2EBAED", "magenta"))
to.show <- cbind(learned.sigs.exome[2,1:96],learned.sigs.exome[2,1:96],learned.sigs.exome[3,1:96],learned.sigs.exome[3,1:96],learned.sigs.exome[5,1:96],learned.sigs.exome[5,1:96])
colnames(to.show) <- c('mlh-1', 'mlh-1 hum', 'pms-2', 'pms-2 hum', 'pole-4;pms-2', 'pole-4;pms-2 hum')
plot_sig_wb(to.show,colors = c("grey29", "#DE1C14", "green3", "blue", "#2EBAED", "magenta"), ymax = 0.1)
plot_sig_wb(to.show,colors = c("grey29", "#DE1C14", "green3", "blue", "#2EBAED", "magenta"), flip = T, ymax=-0.1)
to.show <- cbind(defsigs[1:96,5],learned.sigs.exome[3,1:96],learned.sigs.exome[2,1:96],defsigs[1:96,5])
colnames(to.show) <- c('MMR-1', 'pms-2 hum', 'mlh-1 hum','pms-2 hum1')
plot_sig_wb(to.show, ymax = 0.15)
plot_sig_wb(to.show, ymax = 0.15, colors = c("grey29", "#DE1C14", "green3", "blue", "#2EBAED", "magenta"))
dev.off()