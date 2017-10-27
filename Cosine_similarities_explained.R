############################################################
## Cosine similarity simulations for signature comparison ##
## N. Volkova, EMBL-EBI, 2017                             ##
############################################################

# Functions
cosine <- function(x,y) {
  return(sum(x * y) / sqrt(sum(x**2)) / sqrt(sum(y**2)))
}
cosineM <- function(X) {
  return(sapply(1:ncol(X), function(i) sapply(1:ncol(X), function(j) cosine(X[,i], X[,j]))))
}
poisI <- function(X, theta, Y){
  lambda <- as.numeric(X %*% theta)
  t(Y/lambda^2 * X) %*% X
} # Fisher information matrix for Poisson model

# Plotting
source('plotting functions/plot_sigs.R')

# Generate uniformly distributed profiles in positive cone
maxs <- sample(1:5000, 1000, replace=T)
X <- sapply(1:1000, function(x) runif(n = 104, min = 0, max=maxs[x]))
X <- apply(X,2,function(y) y/sum(y))
XM <- cosineM(X)
#hist(XM[upper.tri(XM)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
#abline(v = quantile(as.vector(XM[upper.tri(XM)]), 0.95), col='red', lty=2) # (0.80)

df = data.frame(x = 1:length(XM[upper.tri(XM)]), val = XM[upper.tri(XM)])
p <- ggplot(data = df, aes(val)) + 
  geom_histogram(binwidth = 0.005, col='black', fill='white') + 
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(family='ArialMT'), panel.border = element_rect(colour='white')) +
  geom_vline(xintercept = quantile(as.vector(XM[upper.tri(XM)]), 0.95), col='red', linetype='dashed') +
  ggtitle('Distribution of similarities between uniform random vectors from positive cone') +
  xlab('Cosine similarity score') + ylab('Frequency') + xlim(c(0,1))
ggsave(plot=p, height=4,width=7, filename = '~/Documents/uniform_sim.pdf', useDingbats=FALSE)

# Distribution of angles in COSMIC signatures
# Download signatures
library(deconstructSigs)
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
types <- as.character(cancer_signatures$Trinucleotide) # trinucleotide classes
types.full <- as.character(cancer_signatures$Somatic.Mutation.Type) # substitution types
row.names(cancer_signatures) <- types.full
cancer_signatures = as.matrix(cancer_signatures[,4:33])
# Get the histogram and heatmap
cosmic <- cosineM(cancer_signatures)
hist(cosmic[upper.tri(cosmic)], prob=T, breaks=20, main='Similarities between COSMIC signatures', xlab = 'Cosine similarity')
lines(density(cosmic[upper.tri(cosmic)], adjust=2), lty="dotted") 
# Heatmap
image.plot(cosmic)
for (x in 1:ncol(cosmic))
  for (y in 1:ncol(cosmic))
    text((x-1)/(ncol(cosmic)-1), (y-1)/(ncol(cosmic)-1), sprintf("%0.2f", cosmic[x,y]))
# ggplot histogram
df = data.frame(x = 1:length(cosmic[upper.tri(cosmic)]), val=cosmic[upper.tri(cosmic)])
p <- ggplot(data = df, aes(val)) + 
  geom_histogram(binwidth = 0.05, col='black', fill='white') + 
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(family='ArialMT'), panel.border = element_rect(colour='white')) +
  geom_vline(xintercept = 0.8, col='red', linetype='dashed') +
  ggtitle('Distribution of similarities between COSMIC cancer signatures') +
  xlab('Cosine similarity score') + ylab('Frequency') + xlim(c(0,1))
ggsave(plot=p, height=4,width=7, filename = '~/Documents/cosmic_sim.pdf', useDingbats=FALSE)

# Another approach: simulation of signatures using Dirichlet(1) distribution,
# it generates signatures uniformly in simplex => cosines are centered in 0.5,
# the higher the dimension - the higher the peak in 0.5.
# Already at 96 the distribution is close to normal.
library(bayesm)
# 96D case - already roughly normal, but centre not exactly in 0.5 (0.51)
Dr1 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(1,96)))
Dr1M <- cosineM(Dr1)
hist(Dr1M[upper.tri(Dr1M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr1M[upper.tri(Dr1M)]), 0.95), col='red', lty=2) # (0.61)
# 1000D case - normal, centre very close to 0.5
Dr1000 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(1,1000)))
Dr1000M <- cosineM(Dr1000)
hist(Dr1000M[upper.tri(Dr1000M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr1000M[upper.tri(Dr1000M)]), 0.95), col='red', lty=2)

####################################################################################################################

# Now draw a signature from bootstrapped COADSTAD signatures and compare to randomly drawn within CI worm signature
load('~/profiles_and_decomposition.RData')
mm <- rbind(COAD.mutation.counts,STAD.mutation.counts)
mut_mat = t(mm) + 0.0001
res <- nmf(mut_mat,rank=8,seed=123456,method='brunet')
defsigs <- NMF::basis(res)
defsigs <- defsigs[,c(8,1,2,5,4,7,3,6)]
colnames(defsigs) <- c("Clock-1", "Clock-2", "POLE", "17-like", "MMR-1", "MMR-2", "MMR-3", "SNP")

###################################################################################################################

# Bootstrapping
sigs = list()
for (i in 1:505) {
  
  mm = rbind(COAD.mutation.counts,STAD.mutation.counts)
  mut_mat <- t(mm[-i,]) + 0.0001
  res <- nmf(mut_mat, rank = 8, method='brunet')
  sigs[[i]] <- NMF::basis(res)
  
}

# Check out stability of signatures in those bootstraps
similarity_hist <- function(defsigs,set,k) {
  hist(sapply(1:ncol(set[[k]]), function(s) cosine(defsigs[,k],set[[k]][,s])), breaks=100)
}

jkset = list()
for (k in 1:8) {
  jkset[[k]] <- sapply(1:505, function(i) 
    sigs[[i]][,which.max(sapply(1:8, function(j) cosine(sigs[[i]][,j],defsigs[,k])))] / sum(sigs[[i]][,which.max(sapply(1:8, function(j) cosine(sigs[[i]][,j],defsigs[,k])))]))
}
similarity_hist(defsigs,jkset,5)
summary(sapply(1:ncol(jkset[[5]]), function(s) cosine(defsigs[,5],jkset[[5]][,s])))

###################################################################################################################

# Now take the worm sigs randomly from their CIs
load('~/Learned_signatures_indel.RData')
learned.sigs <- nmSolve(worm.donor.mut.mat[,1:96],small.X,maxIter=10000, tol = 1e-06, div.err = 1e-10)
Y = worm.donor.mut.mat
mut_matrix_lower <- matrix(NA,nrow=96,ncol=nrow(learned.sigs),dimnames=list(colnames(learned.sigs)[1:96],row.names(learned.sigs)))
mut_matrix_upper <- matrix(NA,nrow=96,ncol=nrow(learned.sigs),dimnames=list(colnames(learned.sigs)[1:96],row.names(learned.sigs)))
library(MASS)
for (i in 1:96) {
  cov.mat <- poisI(as.matrix(small.X),learned.sigs[,i],Y[,i])
  to.keep <- colnames(cov.mat)[which(diag(cov.mat)!=0)]
  if (length(to.keep)!=0) {
    cov.mat <- cov.mat[to.keep,to.keep]
    SVD <- svd(cov.mat)$d
    if (length(which(SVD<1e-6))>0)
      cov.mat <- ginv(cov.mat)
    if (length(which(SVD<1e-6))==0)
      cov.mat <- solve(cov.mat)
    row.names(cov.mat) <- to.keep
    colnames(cov.mat) <- to.keep
    sterr <- sqrt(diag(cov.mat))
    mut_matrix_lower[i,to.keep] <- learned.sigs[to.keep,i] - sterr * qt(0.975,30)
    mut_matrix_upper[i,to.keep] <- learned.sigs[to.keep,i] + sterr * qt(0.975,30)
    mut_matrix_lower[i,to.keep][mut_matrix_lower[i,to.keep]<0] <- 0
  }
} 
# Adjust for trinucleotide difference
load('exome.RData')
learned.sigs.exome <- learned.sigs[,1:96]
mut_matrix_lower.exome <- t(mut_matrix_lower)
mut_matrix_upper.exome <- t(mut_matrix_upper)
for (i in 1:nrow(learned.sigs.exome)) {
  learned.sigs.exome[i,] <- learned.sigs.exome[i,] / trinucleotide.freq.factor.ex[types]
  mut_matrix_lower.exome[i,1:96] <- mut_matrix_lower.exome[i,1:96] / trinucleotide.freq.factor.ex[types] # + learned.sigs.exome[i,1:96] - learned.sigs[i,1:96]
  mut_matrix_upper.exome[i,1:96] <- mut_matrix_upper.exome[i,1:96] / trinucleotide.freq.factor.ex[types] # + learned.sigs.exome[i,1:96] - learned.sigs[i,1:96]
}
mut_matrix_lower.exome[mut_matrix_lower.exome<0] <- 0
# plot signatures with CIs
plot_96_profile_CI(mut_matrix = t(learned.sigs[c(2:3,5),1:96]), mut_matrix_lower = mut_matrix_lower[,c(2:3,5)], mut_matrix_upper = mut_matrix_upper[,c(2:3,5)], CI=TRUE)
plot_96_profile_CI(mut_matrix = t(learned.sigs.exome[c(2:3,5),1:96]), mut_matrix_lower = t(mut_matrix_lower.exome[c(2:3,5),]), mut_matrix_upper = t(mut_matrix_upper.exome[c(2:3,5),]), CI=TRUE)

###################################################################################################################

# Now draw worm and signatures randomly from their confidence intervals

# Similarities for fixed signatures
cosine(learned.sigs.exome[3,1:96],defsigs[1:96,5]) # # 0.85
cosine(learned.sigs.exome[2,1:96],defsigs[1:96,5]) # # 0.81
cosine(learned.sigs.exome[3,-c(35,39,43,47,97:104)],defsigs[-c(35,39,43,47,97:104),5]) # 0.92
cosine(learned.sigs.exome[2,-c(35,39,43,47,97:104)],defsigs[-c(35,39,43,47,97:104),5]) # 0.90
# Drawing random signatures
random_signature <- function(sig, sig_low, sig_up) {
  res <- sapply(1:length(sig), function(i) rnorm(1,mean=sig[i],sd=(sig_up[i]-sig[i])/2))
  res[res<0] <- 0
  res[is.nan(res)] <- sig[is.nan(res)]
  return(res)
}

###################################################################################################################

# pole4;pms-2 double mutant versus single signatures

a1 <- sapply(1:1000, function(s) cosine(random_signature(sig=learned.sigs[2,1:96],sig_low = mut_matrix_lower[,2],sig_up = mut_matrix_upper[,2]),
                                        random_signature(sig=learned.sigs[5,1:96],sig_low = mut_matrix_lower[,5],sig_up = mut_matrix_upper[,5])))
a2 <- sapply(1:1000, function(s) cosine(random_signature(sig=learned.sigs[3,1:96],sig_low = mut_matrix_lower[,3],sig_up = mut_matrix_upper[,3]),
                                        random_signature(sig=learned.sigs[5,1:96],sig_low = mut_matrix_lower[,5],sig_up = mut_matrix_upper[,5])))
a3 <- sapply(1:1000, function(s) cosine(random_signature(sig=learned.sigs[2,1:96],sig_low = mut_matrix_lower[1:96,2],sig_up = mut_matrix_upper[1:96,2]),
                                        random_signature(sig=learned.sigs[3,1:96],sig_low = mut_matrix_lower[1:96,3],sig_up = mut_matrix_upper[1:96,3])))

df = data.frame(x = 1:length(a1), a1, a2, a3)
df = melt(df, id.vars = "x")
levels(df$variable) <- c('mlh-1 vs pole-4;pms-2', 'pms-2 vs pole-4;pms-2', 'mlh-1 vs pms-2')
df2 <- data.frame(cos = c(cosine(learned.sigs[2,1:96],learned.sigs[5,1:96]),
                           cosine(learned.sigs[3,1:96],learned.sigs[5,1:96]),
                           cosine(learned.sigs[2,1:96],learned.sigs[3,1:96])),
                  variable = c('mlh-1 vs pole-4;pms-2', 'pms-2 vs pole-4;pms-2', 'mlh-1 vs pms-2'))
p <- ggplot(data = df, aes(value)) + 
  geom_histogram(binwidth = 0.01, col='black', fill='white') + 
  facet_grid(variable ~ ., scales='fixed') +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(family='ArialMT'), strip.background = element_rect(colour='white', fill='white'),
                     panel.border = element_rect(colour = 'white'), strip.text = element_text(size=12)) +
  geom_vline(data=df2, aes(xintercept=cos), col='red', linetype='dashed') +
  ggtitle('Distribution of similarities between C. elegans signatures') +
  xlab('Cosine similarity score') + ylab('Frequency') + xlim(c(0,1))
p
ggsave(plot=p, height=7,width=7, filename = '~/Documents/worm-sig-sim.pdf', useDingbats=FALSE)

###################################################################################################################

# mlh-1 and pms-2 vs MMR-1

a1 <- sapply(1:505, function(s) cosine(random_signature(sig=learned.sigs.exome[3,1:96],sig_low = mut_matrix_lower.exome[3,],sig_up = mut_matrix_upper.exome[3,]),jkset[[5]][-c(97:104),s]))
a2 <- sapply(1:505, function(s) cosine(random_signature(sig=learned.sigs.exome[3,c(1:96)[-c(35,39,43,47)]],sig_low = mut_matrix_lower.exome[3,-c(35,39,43,47)],
                                                                   sig_up = mut_matrix_upper.exome[3,-c(35,39,43,47)]),jkset[[5]][-c(35,39,43,47,97:104),s]))
b1 <- sapply(1:505, function(s) cosine(random_signature(sig=learned.sigs.exome[2,1:96],sig_low = mut_matrix_lower.exome[3,],
                                                                   sig_up = mut_matrix_upper.exome[3,]),jkset[[5]][-c(97:104),s]))
b2 <- sapply(1:505, function(s) cosine(random_signature(sig=learned.sigs.exome[2,c(1:96)[-c(35,39,43,47)]],sig_low = mut_matrix_lower.exome[3,-c(35,39,43,47)],
                                                                   sig_up = mut_matrix_upper.exome[3,-c(35,39,43,47)]),jkset[[5]][-c(35,39,43,47,97:104),s]))


a1 <- sapply(1:505, function(s) cosine(learned.sigs.exome[3,1:96],jkset[[5]][-c(97:104),s]))
a2 <- sapply(1:505, function(s) cosine(learned.sigs.exome[3,c(1:96)[-c(35,39,43,47)]],jkset[[5]][-c(35,39,43,47,97:104),s]))
b1 <- sapply(1:505, function(s) cosine(learned.sigs.exome[2,1:96],jkset[[5]][-c(97:104),s]))
b2 <- sapply(1:505, function(s) cosine(learned.sigs.exome[2,c(1:96)[-c(35,39,43,47)]],jkset[[5]][-c(35,39,43,47,97:104),s]))


df = data.frame(x = 1:length(a1), a1, a2)
df = melt(df, id.vars = "x")
levels(df$variable) <- c('Full signatures', 'Without C>T at CpG sites')
df2 <- data.frame(cos = c(cosine(learned.sigs.exome[3,1:96],defsigs[1:96,5]),
                          cosine(learned.sigs.exome[3,-c(35,39,43,47,97:104)],defsigs[-c(35,39,43,47,97:104),5])),
                  variable =c('Full signatures', 'Without C>T at CpG sites'))
p <- ggplot(data = df, aes(value)) + 
  geom_histogram(binwidth = 0.01, col='black', fill='white') + 
  facet_grid(. ~ variable, scales='fixed') +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(family='ArialMT'), strip.background = element_rect(colour='white', fill='white'),
                     panel.border = element_rect(colour = 'white'), strip.text = element_text(size=12)) +
  geom_vline(data=df2, aes(xintercept = cos), col='red', linetype='dashed') +
  ggtitle('Distribution of similarities between pms-2 and MMR-1 signatures') +
  xlab('Cosine similarity score') + ylab('Frequency') + xlim(c(0,1))
p
ggsave(plot=p, height=4,width=8, filename = '~/Documents/pms2jk.pdf', useDingbats=FALSE)

df = data.frame(x = 1:length(a1), b1, b2)
df = melt(df, id.vars = "x")
levels(df$variable) <- c('Full signatures', 'Without C>T at CpG sites')
df2 <- data.frame(cos = c(cosine(learned.sigs.exome[2,1:96],defsigs[1:96,5]),
                          cosine(learned.sigs.exome[2,-c(35,39,43,47,97:104)],defsigs[-c(35,39,43,47,97:104),5])),
                  variable =c('Full signatures', 'Without C>T at CpG sites'))
p <- ggplot(data = df, aes(value)) + 
  geom_histogram(binwidth = 0.01, col='black', fill='white') + 
  facet_grid(. ~ variable, scales='fixed') +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(family='ArialMT'), strip.background = element_rect(colour='white', fill='white'),
                     panel.border = element_rect(colour = 'white'), strip.text = element_text(size=12)) +
  geom_vline(data=df2, aes(xintercept = cos), col='red', linetype='dashed') +
  ggtitle('Distribution of similarities between mlh-1 and MMR-1 signatures') +
  xlab('Cosine similarity score') + ylab('Frequency') + xlim(c(0,1))
p
ggsave(plot=p, height=4,width=8, filename = '~/Documents/mlh1jk.pdf', useDingbats=FALSE)