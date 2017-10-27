# Function for plotting 96-long signatures of relative contributions of different mutation types
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), samples x signatures
# colors - colors to use for 6 mutation classes
# ymin, ymax = max and min value on the plot
# flip - FALSE or TRUE, flips the plot vertically
# font - size of mutation class titles
plot_sig <- function (mut_matrix, colors, ymin=0,ymax = 0.15,flip=F,font=11) # plotting 96-signatures
{
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  IND_TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  TRIPLETS_112 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3),IND_TRIPLETS)
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  if (missing(colors)) {
    colors = c("#2EBAED", "#000000", "#DE1C14",
                    "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  context = TRIPLETS_96
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  if (flip==T) norm_mut_matrix = -norm_mut_matrix
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (ymax>0) breaks = c(0,0.1) else breaks = c(-0.1,0)
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
          geom_bar(stat = "identity", colour = "black",size = 0.2) +
          scale_fill_manual(values = colors) + facet_grid(variable ~ substitution) +
          ylab("Relative contribution") + coord_cartesian(ylim = c(ymin,ymax)) + scale_y_continuous(breaks = breaks) +
          guides(fill = FALSE) + theme_bw() +
          theme(axis.title.y = element_text(size = 14,vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 14), 
                  axis.text.x = element_text(size = 7, angle = 90, vjust = 0.4), 
                  strip.text.x = element_text(size = 11), strip.text.y = element_text(size = font), 
                  panel.grid.major.x = element_blank())
  return(plot)
}
# Function for plotting 96-long signatures of relative contributions of different mutation types: different design
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), samples x signatures
# colors - colors to use for 6 mutation classes
# ymin, ymax = max and min value on the plot
# flip - FALSE or TRUE, flips the plot vertically
# size - size of mutation class titles
plot_sig_wb <- function (mut_matrix, colors, ymin=0,ymax = 0.15,flip=F,size=12) # plotting 96-signatures with no lines
{
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  IND_TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  TRIPLETS_112 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3),IND_TRIPLETS)
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  if (missing(colors)) {
    colors = c("#2EBAED", "#000000", "#DE1C14",
               "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  context = TRIPLETS_96
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  if (flip==T) norm_mut_matrix = -norm_mut_matrix
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (ymax>0) breaks = c(0,0.1) else breaks = c(-0.1,0)
  if (ymax>0) lbls = c('0','0.1') else lbls = c('0.1','0')
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution) +
    ylab("Relative contribution") + coord_cartesian(ylim = c(ymin,ymax)) + 
    guides(fill = FALSE) + theme_bw() + scale_y_continuous(breaks = breaks, labels = lbls) +
    theme(text = element_text(family='ArialMT'),
          axis.title.y = element_text(size = 14,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_text(size = 14), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  return(plot)
}

# Function for plotting 96-long mutational profiles (in absolute numbers)
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), samples x signatures
# colors - colors to use for 6 mutation classes
plot_profiles <- function (mut_matrix, colors) # plot mutational profiles without normalizing
{
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  if (missing(colors)) {
    colors = c(
      "#2EBAED", "#000000", "#DE1C14",
      "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6")
  }
  context = TRIPLETS_96
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution) +
    ylab("Relative contribution") + coord_cartesian(ylim = c(0,max(mut_matrix))) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 14,vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 14), 
          axis.text.x = element_text(size = 7, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 13), 
          panel.grid.major.x = element_blank())
  return(plot)
}
# Function for plotting 104-long signatures of relative contributions of different mutation types
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), mutation types x signatures
# colors - colors to use for 8 mutation classes
# size - size of mutation class titles
# ymax - maximal value on the plot
plot_sig_104 <- function(mut_matrix,colors=NA,size=8,ymax=0.2) { # plot 104-signatures
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  if (is.na(colors)) colors=c("#2EBAED", "#000000", "#DE1C14","orange","purple","#D4D2D2", "#ADCC54", "#F0D0CE")
  
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  df = as.data.frame(norm_mut_matrix)
  df$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df$Base = c(types,rep(c("A","C","G","T"),2))
  substring(df$Base[nchar(df$Base)>1], 2, 2) = "*"
  #colnames(df)[1:ncol(mut_matrix)] = paste0("Signature.",c(1:ncol(mut_matrix)),sep="")
  df2 = melt(df,id.vars = c("Type","Base"))
  df2$Type_f = factor(df2$Type, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DEL","INS"))
  plot = ggplot(data = df2, aes(x = Base, y = value, fill = Type, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ Type_f,scales = "free_x") +
    ylab("Relative contribution") + coord_cartesian(ylim = c(0,ymax)) + 
    guides(fill = FALSE) + theme_bw() + scale_y_continuous(breaks = c(0,0.15)) +
    theme(text = element_text(family='ArialMT'),
          axis.title.y = element_text(size = 14,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_text(size = 14), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  return(plot)
}
# Function for plotting 104-long signatures of relative contributions of different mutation types with confidence intervals
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), mutation types x signatures
# mut_matrix_lower - matrix of lower CI values for signatures, mutation types x signatures
# mut_matrix_upper - matrix of upper CI values for signatures, mutation types x signatures
# colors - colors to use for 8 mutation classes
# size - size of mutation class titles
# ymax - maximal value on the plot
plot_sig_104_CI <- function(mut_matrix,mut_matrix_lower,mut_matrix_upper,size=8,ymax=0.2,colors) { # plot 104-signatures with confidence intervals
  
  colors = c("#2EBAED", "#000000", "#DE1C14","orange","purple","#D4D2D2", "#ADCC54", "#F0D0CE")
  
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  df = as.data.frame(mut_matrix)
  df$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df$Base = c(types,rep(c("A","C","G","T"),2))
  substring(df$Base[nchar(df$Base)>1], 2, 2) = "*"
  df2 = melt(df,id.vars = c("Type","Base"))
  df2$Type_f = factor(df2$Type, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DEL","INS"))
  plot = ggplot(data = df2, aes(x = Base, y = value, fill = Type, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ Type_f,scales = "free_x") +
    ylab("Relative contribution") + coord_cartesian(ylim = c(0,ymax)) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 14,vjust = 1), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 14), 
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 10), strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank())
  
  mut_matrix_lower <- mut_matrix_lower[row.names(mut_matrix),]
  mut_matrix_upper <- mut_matrix_upper[row.names(mut_matrix),]
  rownames(mut_matrix_lower) = rownames(mut_matrix_upper) = NULL
  df_lower = as.data.frame(mut_matrix_lower)
  df_upper = as.data.frame(mut_matrix_upper)
  df_lower$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df_lower$Base = c(types,rep(c("A","C","G","T"),2))
  df_upper$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df_upper$Base = c(types,rep(c("A","C","G","T"),2))
  df_lower = melt(df_lower, id.vars = c("Type","Base"))
  df_upper = melt(df_upper, id.vars = c("Type","Base"))
  df3 <- cbind(df2, value_min = df_lower$value, value_max = df_upper$value)
  df3$Type_f = factor(df3$Type, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DEL","INS"))
  plot2 = plot + geom_pointrange(data=df3, aes(ymin=value,ymax=value_max,colour=factor(Type)), fatten = 0.01, size=0.5, show.legend = F, colour="white") +
    scale_color_manual(values=colors) +
    geom_pointrange(data=df3, aes(ymin=value_min,ymax=value), size=0.5, fatten = 0.01,show.legend = F, colour="white")
  
  return(plot2)
}
# Function for plotting 104-long mutational profiles (normalized, i.e. in relative contribution)
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), samples x signatures
# colors - colors to use for 8 mutation classes
# boxplot - FALSE or TRUE, merges all the values across profiles into a boxplot
# normalize - FALSE or TRUE, turns absolute numbers into relative contributions
# size - regulates font size of the text on the plot
plot_profiles_104 <- function (mut_matrix, colors=NA, boxplot=F, normalize=F, size=6) # plot 104-long mutational profiles
{
  
  if (is.na(colors)) colors = c("#2EBAED", "#000000", "#DE1C14","orange","purple","#D4D2D2", "#ADCC54", "#F0D0CE")
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  if (normalize) {
    mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  }
  df = as.data.frame(mut_matrix)
  df$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df$Base = c(types,rep(c("A","C","G","T"),2))
  substring(df$Base[nchar(df$Base)>1], 2, 2) = "*"
  df2 = melt(df,id.vars = c("Type","Base"))
  df2$Type_f = factor(df2$Type, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DEL","INS"))
  if (boxplot) {
    plot = ggplot(data = df2, aes(x = Base, y = value, fill=Type)) +
      geom_boxplot() +
      scale_fill_manual(values = colors) + 
      facet_grid(. ~ Type_f, scales="free") +
      ylab("Mutation counts") + 
      guides(fill = FALSE) + theme_bw() +
      theme(axis.title.y = element_text(size = 20,vjust = 1), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 20), 
            axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
            strip.text.x = element_text(size = size), strip.text.y = element_text(size = size), 
            panel.grid.major.x = element_blank())
  } else {
    plot = ggplot(data = df2, aes(x = Base, y = value, fill = Type, width = 0.6)) +
      geom_bar(stat = "identity", colour = "black",size = 0.2) +
      scale_fill_manual(values = c("#2EBAED", "#000000", "#DE1C14","orange","purple","#D4D2D2", "#ADCC54", "#F0D0CE")) + 
      facet_grid(variable ~ Type_f,scales = "free_x") +
      ylab("Mutation counts") + coord_cartesian(ylim = c(0,max(mut_matrix))) + 
      guides(fill = FALSE) + theme_bw() +
      theme(axis.title.y = element_text(size = 14,vjust = 1), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 14), 
            axis.text.x = element_text(size = 6, angle = 90, vjust = 0.4), 
            strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6), 
            panel.grid.major.x = element_blank())
  }
  return(plot)
}
# Function for plotting 104-long mutational profiles (in absolute numbers) with confidence intervals
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), samples x signatures
# size - regulates font size of the text on the plot
plot_104_profile_CI <- function (mut_matrix, mut_matrix_lower = NULL, 
                                mut_matrix_upper = NULL,  
                                colors=c("deepskyblue2","black","red3","orange","purple","grey","olivedrab3","pink"),size=6,ymax=0.25) # plot 104-long profiles with confidence intervals
{
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  
  norm_mut_matrix = mut_matrix
  mult <- colSums(norm_mut_matrix)
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  df = as.data.frame(norm_mut_matrix)
  df$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df$Base = c(types,rep(c("A","C","G","T"),2))
  df2 = melt(df,id.vars = c("Type","Base"))
  df2$Type_f = factor(df2$Type, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DEL","INS"))
  plot = ggplot(data = df2, aes(x = Base, y = value, fill = Type, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ Type_f,scales = "free_x") +
    ylab("Mutation counts") + coord_cartesian(ylim = c(0,ymax)) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 16,vjust = 1), axis.text.y = element_text(size = size), axis.title.x = element_text(size = 16), 
          axis.text.x = element_text(size = size, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size), strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank())
  
    mut_matrix_lower <- mut_matrix_lower[row.names(mut_matrix),]
    mut_matrix_upper <- mut_matrix_upper[row.names(mut_matrix),]
    norm_mut_matrix_lower <- mut_matrix_lower
    norm_mut_matrix_upper <- mut_matrix_upper
    for (i in 1:ncol(norm_mut_matrix_lower)) {
      norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
      norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
    }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = as.data.frame(norm_mut_matrix_lower)
    df_upper = as.data.frame(norm_mut_matrix_upper)
    df_lower$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
    df_lower$Base = c(types,rep(c("A","C","G","T"),2))
    df_upper$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
    df_upper$Base = c(types,rep(c("A","C","G","T"),2))
    df_lower = melt(df_lower, id.vars = c("Type","Base"))
    df_upper = melt(df_upper, id.vars = c("Type","Base"))
    df3 <- cbind(df2, value_min = df_lower$value, value_max = df_upper$value)
    plot2 = plot + geom_pointrange(data=df3, aes(ymin=value,ymax=value_max,colour = Type_f), fatten = 0.01, size=1, show.legend = F) +
      scale_color_manual(values=c(colors,"white")) +
      geom_pointrange(data=df3, aes(ymin=value_min,ymax=value,colour="white"), size=1, fatten = 0.01,show.legend = F)
  
  return(plot2)
}
# Function for plotting 96-long mutational profiles of different mutation types with confidence intervals
# requires ggplot2 and reshape2 libraries
# mut_matrix - matrix of signatures (or profiles), mutation types x signatures
# mut_matrix_lower - matrix of lower CI values for signatures, mutation types x signatures
# mut_matrix_upper - matrix of upper CI values for signatures, mutation types x signatures
# colors - colors to use for 8 mutation classes
# size - size of mutation class titles
# ymax - maximal value on the plot
# CI - FALSE or TRUE, do not / do plot CIs
plot_96_profile_CI <- function (mut_matrix, mut_matrix_lower = NULL, 
                                mut_matrix_upper = NULL, CI=FALSE, 
                                colors=c("deepskyblue2","black","red3","grey","olivedrab3","pink"), ymax=NA, size=9) 
{
  if(is.na(ymax)) ymax = max(mut_matrix_upper)
  if (ymax>1) step = 10
  if (ymax<=1) step = 0.1
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  context = types
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ substitution) + 
    coord_cartesian(ylim = c(0,ymax)) + 
    scale_y_continuous(breaks = seq(0, ymax, step)) +
    guides(fill = FALSE) + theme_bw() + 
    theme(text = element_text(family='ArialMT'),
          axis.title.y = element_text(size = 14,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_text(size = 14), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines')) +
    labs(y="Number of mutations", x="Context")
  
  if (CI) { # add confidence intervals
    rownames(mut_matrix_upper) = rownames(mut_matrix_lower) = NULL
    df_lower = cbind(df, as.data.frame(mut_matrix_lower))
    df_upper = cbind(df, as.data.frame(mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df2 = cbind(df, as.data.frame(mut_matrix[,colnames(mut_matrix_lower)]))
    df3 = melt(df2, id.vars = c("substitution", "context"))
    df3 <- cbind(df3, value_min = df_lower$value, value_max = df_upper$value)
    plot = plot + geom_pointrange(data=df3, aes(ymin=value,ymax=value_max,colour = substitution), size=1.5, fatten = 0.001,show.legend = F) +
      scale_color_manual(values=c(colors,"white")) +
      geom_pointrange(data=df3, aes(ymin=value_min,ymax=value,colour="white"), size=1.5, fatten = 0.001,show.legend = F)
  }
  return(plot)
}
