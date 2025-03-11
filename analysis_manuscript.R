library(ggplot2)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(ggpubr)
library(viridis)
library(scales)
library(gmodels)
library(DescTools)
library(MetBrewer)
library(svglite)



theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}



scaleFUN <- function(x) sprintf("%.2f", x)

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")

palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[7], palette[3])
pallete3 <- c(palette[1], palette[3], palette[5],palette[7])

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(pzfx)


setwd("Z:/kacottre/Manuscripts/Off-target_manuscript/Source Data/")



#read western (immunoblot) data for knockdowns
immunoblot_skbr3_adar_ddx54 <- fread("immunoblot_skbr3_adar_ddx54.txt")

immunoblot_skbr3_adar_ddx54 <- na.omit(immunoblot_skbr3_adar_ddx54)

immunoblot_skbr3_adar_ddx54 <- subset(immunoblot_skbr3_adar_ddx54, !`fold_change` == "pending")


immunoblot_skbr3_adar_ddx54$sample <- paste(immunoblot_skbr3_adar_ddx54$shRNA_1, immunoblot_skbr3_adar_ddx54$shRNA_2, sep = "/")

immunoblot_skbr3_adar_ddx54$sample <- gsub("-", ".", immunoblot_skbr3_adar_ddx54$sample)

immunoblot_skbr3_adar_ddx54$fold_change <- as.numeric(immunoblot_skbr3_adar_ddx54$fold_change)

#group and summarise data
immunoblot_skbr3_adar_ddx54_sum <- dplyr::group_by(immunoblot_skbr3_adar_ddx54, sample, Gene)

immunoblot_skbr3_adar_ddx54_sum <- dplyr::summarise(immunoblot_skbr3_adar_ddx54_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_skbr3_adar_ddx54_sum$ci <- qnorm(0.975)*(immunoblot_skbr3_adar_ddx54_sum$sd_expression/sqrt(4))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx54, Gene == "p-PKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_skbr3_adar_ddx54$fold_change[immunoblot_skbr3_adar_ddx54$Gene == "p-PKR"]*1.1)

tppkr <- tppkr[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tppkr)), ]

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shSCR/shDDX54.5", "shADAR/shSCR", "shADAR/shDDX54.4", "shADAR/shDDX54.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shSCR\nshDDX54-5", "shADAR1\nshSCR", "shADAR1\nshDDX54-4", "shADAR1\nshDDX54-5")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2)


tppkr <- tppkr[!grepl("shDDX54.5", rownames(tppkr)), ]


ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_ddx54_skbr3-4.tiff", height = 4.2, width = 4.2)



aov_peif2a <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx54, Gene == "p-eIF2a"))
tpeif2a <- TukeyHSD(aov_peif2a)
tpeif2a <- as.data.frame(tpeif2a$sample[,1:4])

out <- strsplit(as.character(row.names(tpeif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpeif2a$stars <- makeStars(tpeif2a$`p adj`)

tpeif2a$xmin <- out$V1
tpeif2a$xmax <- out$V2
tpeif2a$y <- max(immunoblot_skbr3_adar_ddx54$fold_change[immunoblot_skbr3_adar_ddx54$Gene == "p-eIF2a"]*1.1)

tpeif2a <- tpeif2a[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tpeif2a)), ]

tpeif2a <- subset(tpeif2a, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shSCR/shDDX54.5", "shADAR/shSCR", "shADAR/shDDX54.4", "shADAR/shDDX54.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shSCR\nshDDX54-5", "shADAR1\nshSCR", "shADAR1\nshDDX54-4", "shADAR1\nshDDX54-5")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/peif2a_adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2)



aov_cparp <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx54, Gene == "c-PARP"))
tcparp <- TukeyHSD(aov_cparp)
tcparp <- as.data.frame(tcparp$sample[,1:4])

out <- strsplit(as.character(row.names(tcparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tcparp$stars <- makeStars(tcparp$`p adj`)

tcparp$xmin <- out$V1
tcparp$xmax <- out$V2
tcparp$y <- max(immunoblot_skbr3_adar_ddx54$fold_change[immunoblot_skbr3_adar_ddx54$Gene == "c-PARP"]*1.1)

tcparp <- tcparp[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tcparp)), ]

tcparp <- subset(tcparp, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shSCR/shDDX54.5", "shADAR/shSCR", "shADAR/shDDX54.4", "shADAR/shDDX54.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shSCR\nshDDX54-5", "shADAR1\nshSCR", "shADAR1\nshDDX54-4", "shADAR1\nshDDX54-5")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tcparp$xmin, 
               xmax = tcparp$xmax, 
               y.position = tcparp$y , label = tcparp$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2)


aov_ddx54 <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx54, Gene == "DDX54"))
tddx54 <- TukeyHSD(aov_ddx54)
tddx54 <- as.data.frame(tddx54$sample[,1:4])

out <- strsplit(as.character(row.names(tddx54)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tddx54$stars <- makeStars(tddx54$`p adj`)

tddx54$xmin <- out$V1
tddx54$xmax <- out$V2
tddx54$y <- max(immunoblot_skbr3_adar_ddx54$fold_change[immunoblot_skbr3_adar_ddx54$Gene == "DDX54"]*1.1)

tddx54 <- tddx54[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tddx54)), ]

tddx54 <- subset(tddx54, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "DDX54"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shSCR/shDDX54.5", "shADAR/shSCR", "shADAR/shDDX54.4", "shADAR/shDDX54.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shSCR\nshDDX54-5", "shADAR1\nshSCR", "shADAR1\nshDDX54-4", "shADAR1\nshDDX54-5")) + 
  labs(x = "", y = "Relative DDX54 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "DDX54"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "DDX54"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tddx54$xmin, 
               xmax = tddx54$xmax, 
               y.position = tddx54$y , label = tddx54$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ddx54_adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2)


aov_adar <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx54, Gene == "ADARp110"))
tadar <- TukeyHSD(aov_adar)
tadar <- as.data.frame(tadar$sample[,1:4])

out <- strsplit(as.character(row.names(tadar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tadar$stars <- makeStars(tadar$`p adj`)

tadar$xmin <- out$V1
tadar$xmax <- out$V2
tadar$y <- max(immunoblot_skbr3_adar_ddx54$fold_change[immunoblot_skbr3_adar_ddx54$Gene == "ADARp110"]*1.1)

tadar <- tadar[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tadar)), ]

tadar <- subset(tadar, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shSCR/shDDX54.5", "shADAR/shSCR", "shADAR/shDDX54.4", "shADAR/shDDX54.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shSCR\nshDDX54-5", "shADAR1\nshSCR", "shADAR1\nshDDX54-4", "shADAR1\nshDDX54-5")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tadar$xmin, 
               xmax = tadar$xmax, 
               y.position = tadar$y , label = tadar$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2)


aov_pkr <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx54, Gene == "total PKR"))
tpkr <- TukeyHSD(aov_pkr)
tpkr <- as.data.frame(tpkr$sample[,1:4])

out <- strsplit(as.character(row.names(tpkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpkr$stars <- makeStars(tpkr$`p adj`)

tpkr$xmin <- out$V1
tpkr$xmax <- out$V2
tpkr$y <- max(immunoblot_skbr3_adar_ddx54$fold_change[immunoblot_skbr3_adar_ddx54$Gene == "total PKR"]*1.1)

tpkr <- tpkr[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tpkr)), ]

tpkr <- subset(tpkr, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx54, Gene == "Total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shSCR/shDDX54.5", "shADAR/shSCR", "shADAR/shDDX54.4", "shADAR/shDDX54.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shSCR\nshDDX54-5", "shADAR1\nshSCR", "shADAR1\nshDDX54-4", "shADAR1\nshDDX54-5")) + 
  labs(x = "", y = "Relative total PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx54_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/pkr_adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2)

aov_pkr <- NULL
tpkr <- NULL
aov_peif2a <- NULL
t_peif2a <- NULL
aov_ppkr <- NULL
tppkr <- NULL
aov_cparp <- NULL
tcparp <- NULL
aov_adar <- NULL
tadar <- NULL
aov_pkr <- NULL
tpkr <- NULL
aov_ddx54 <- NULL
tddx54 <- NULL




######MCF-7 knockdowns

#read western (immunoblot) data for knockdowns
immunoblot_mcf7_adar_ddx54 <- fread("immunoblot_mcf7_adar_ddx54.txt")

immunoblot_mcf7_adar_ddx54 <- na.omit(immunoblot_mcf7_adar_ddx54)

immunoblot_mcf7_adar_ddx54 <- subset(immunoblot_mcf7_adar_ddx54, !`fold_change` == "pending")


immunoblot_mcf7_adar_ddx54$sample <- paste(immunoblot_mcf7_adar_ddx54$shRNA_1, immunoblot_mcf7_adar_ddx54$shRNA_2, sep = "/")

immunoblot_mcf7_adar_ddx54$sample <- gsub("-", ".", immunoblot_mcf7_adar_ddx54$sample)

immunoblot_mcf7_adar_ddx54$fold_change <- as.numeric(immunoblot_mcf7_adar_ddx54$fold_change)

#group and summarise data
immunoblot_mcf7_adar_ddx54_sum <- dplyr::group_by(immunoblot_mcf7_adar_ddx54, sample, Gene)

immunoblot_mcf7_adar_ddx54_sum <- dplyr::summarise(immunoblot_mcf7_adar_ddx54_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_mcf7_adar_ddx54_sum$ci <- qnorm(0.975)*(immunoblot_mcf7_adar_ddx54_sum$sd_expression/sqrt(4))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_ddx54, Gene == "p-PKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_mcf7_adar_ddx54$fold_change[immunoblot_mcf7_adar_ddx54$Gene == "p-PKR"]*1.1)

tppkr <- tppkr[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tppkr)), ]

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_ddx54, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2)




aov_peif2a <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_ddx54, Gene == "p-eIF2a"))
tpeif2a <- TukeyHSD(aov_peif2a)
tpeif2a <- as.data.frame(tpeif2a$sample[,1:4])

out <- strsplit(as.character(row.names(tpeif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpeif2a$stars <- makeStars(tpeif2a$`p adj`)

tpeif2a$xmin <- out$V1
tpeif2a$xmax <- out$V2
tpeif2a$y <- max(immunoblot_mcf7_adar_ddx54$fold_change[immunoblot_mcf7_adar_ddx54$Gene == "p-eIF2a"]*1.1)

tpeif2a <- tpeif2a[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tpeif2a)), ]

tpeif2a <- subset(tpeif2a, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_ddx54, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) +
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/peif2a_adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2)



aov_cparp <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_ddx54, Gene == "c-PARP"))
tcparp <- TukeyHSD(aov_cparp)
tcparp <- as.data.frame(tcparp$sample[,1:4])

out <- strsplit(as.character(row.names(tcparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tcparp$stars <- makeStars(tcparp$`p adj`)

tcparp$xmin <- out$V1
tcparp$xmax <- out$V2
tcparp$y <- max(immunoblot_mcf7_adar_ddx54$fold_change[immunoblot_mcf7_adar_ddx54$Gene == "c-PARP"]*1.1)

tcparp <- tcparp[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tcparp)), ]

tcparp <- subset(tcparp, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_ddx54, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) +
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tcparp$xmin, 
               xmax = tcparp$xmax, 
               y.position = tcparp$y , label = tcparp$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2)


aov_ddx54 <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_ddx54, Gene == "DDX54"))
tddx54 <- TukeyHSD(aov_ddx54)
tddx54 <- as.data.frame(tddx54$sample[,1:4])

out <- strsplit(as.character(row.names(tddx54)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tddx54$stars <- makeStars(tddx54$`p adj`)

tddx54$xmin <- out$V1
tddx54$xmax <- out$V2
tddx54$y <- max(immunoblot_mcf7_adar_ddx54$fold_change[immunoblot_mcf7_adar_ddx54$Gene == "DDX54"]*1.1)

tddx54 <- tddx54[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tddx54)), ]

tddx54 <- subset(tddx54, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_ddx54, Gene == "DDX54"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) +
  labs(x = "", y = "Relative DDX54 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "DDX54"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "DDX54"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tddx54$xmin, 
               xmax = tddx54$xmax, 
               y.position = tddx54$y , label = tddx54$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ddx54_adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2)


aov_adar <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_ddx54, Gene == "ADARp110"))
tadar <- TukeyHSD(aov_adar)
tadar <- as.data.frame(tadar$sample[,1:4])

out <- strsplit(as.character(row.names(tadar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tadar$stars <- makeStars(tadar$`p adj`)

tadar$xmin <- out$V1
tadar$xmax <- out$V2
tadar$y <- max(immunoblot_mcf7_adar_ddx54$fold_change[immunoblot_mcf7_adar_ddx54$Gene == "ADARp110"]*1.1)

tadar <- tadar[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tadar)), ]

tadar <- subset(tadar, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_ddx54, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) +
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tadar$xmin, 
               xmax = tadar$xmax, 
               y.position = tadar$y , label = tadar$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2)


aov_pkr <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_ddx54, Gene == "total PKR"))
tpkr <- TukeyHSD(aov_pkr)
tpkr <- as.data.frame(tpkr$sample[,1:4])

out <- strsplit(as.character(row.names(tpkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpkr$stars <- makeStars(tpkr$`p adj`)

tpkr$xmin <- out$V1
tpkr$xmax <- out$V2
tpkr$y <- max(immunoblot_mcf7_adar_ddx54$fold_change[immunoblot_mcf7_adar_ddx54$Gene == "total PKR"]*1.1)

tpkr <- tpkr[!grepl("shDDX54.4.*shDDX54.5|shDDX54.5.*shDDX54.4", rownames(tpkr)), ]

tpkr <- subset(tpkr, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_ddx54, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) +  labs(x = "", y = "Relative total PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_ddx54_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/pkr_adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2)

aov_pkr <- NULL
tpkr <- NULL
aov_peif2a <- NULL
t_peif2a <- NULL
aov_ppkr <- NULL
tppkr <- NULL
aov_cparp <- NULL
tcparp <- NULL
aov_adar <- NULL
tadar <- NULL
aov_pkr <- NULL
tpkr <- NULL
aov_ddx54 <- NULL
tddx54 <- NULL






####Immunoblot TNBC lines
immunoblot_tnbc_ddx54 <- fread("immunoblot_bt549_mb231_ddx54.txt")

immunoblot_tnbc_ddx54 <- na.omit(immunoblot_tnbc_ddx54)

immunoblot_tnbc_ddx54$shRNA <- gsub("-", ".", immunoblot_tnbc_ddx54$shRNA)

immunoblot_tnbc_ddx54$Cell_line <- factor(immunoblot_tnbc_ddx54$Cell_line, levels = c("MB231", "BT549"))

#group and summarise data
immunoblot_tnbc_ddx54_sum <- dplyr::group_by(immunoblot_tnbc_ddx54, Cell_line, Gene, shRNA)

immunoblot_tnbc_ddx54_sum<- dplyr::summarise(immunoblot_tnbc_ddx54_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

tnbc_ddx54_ppkr_mb231 <-subset(immunoblot_tnbc_ddx54, Gene == "p-PKR" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_ppkr_mb231)
summary(mb231)

tnbc_ddx54_ppkr_mb231$shRNA <- factor(as.factor(tnbc_ddx54_ppkr_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231ppkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_ppkr_mb231)
mb231ppkr <- as.data.frame(mb231ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231ppkr$group1 <- out$V1
mb231ppkr$group2 <- out$V2
mb231ppkr$pval <- makeStars(mb231ppkr$pval)
mb231ppkr$y.position <- max(tnbc_ddx54_ppkr_mb231$fold_change*1.1)
mb231ppkr$Cell_line = "MB231"


tnbc_ddx54_ppkr_bt549 <-subset(immunoblot_tnbc_ddx54, Gene == "p-PKR" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_ppkr_bt549)
summary(bt549)

tnbc_ddx54_ppkr_bt549$shRNA <- factor(as.factor(tnbc_ddx54_ppkr_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549ppkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_ppkr_bt549)
bt549ppkr <- as.data.frame(bt549ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ppkr$group1 <- out$V1
bt549ppkr$group2 <- out$V2
bt549ppkr$pval <- makeStars(bt549ppkr$pval)
bt549ppkr$y.position <- max(tnbc_ddx54_ppkr_bt549$fold_change*1.1)
bt549ppkr$Cell_line = "BT549"


tnbcppkr <- rbind(mb231ppkr, bt549ppkr)
tnbcppkr$Cell_line <- factor(tnbcppkr$Cell_line, levels = c("MB231", "BT549"))
tnbcppkr <- subset(tnbcppkr, !pval == "ns")

ggplot(subset(immunoblot_tnbc_ddx54, Gene == "p-PKR"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") +  facet_grid(cols = vars(Cell_line)) +
  scale_colour_manual(values = pallete3, 
                      labels = c(MB231 = "MDA-MB-231", BT549 = "BT-549")) +
  geom_col(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "p-PKR", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "p-PKR", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcppkr, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line") 
ggsave("Western_plots/ppkr_tnbc_ddx54_facet.tiff", height =4.5, width = 4, units = "in")


tnbc_ddx54_cparp_mb231 <-subset(immunoblot_tnbc_ddx54, Gene == "c-PARP" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_cparp_mb231)
summary(mb231)

tnbc_ddx54_cparp_mb231$shRNA <- factor(as.factor(tnbc_ddx54_cparp_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231cparp <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_cparp_mb231)
mb231cparp <- as.data.frame(mb231cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231cparp$group1 <- out$V1
mb231cparp$group2 <- out$V2
mb231cparp$pval <- makeStars(mb231cparp$pval)
mb231cparp$y.position <- max(tnbc_ddx54_cparp_mb231$fold_change*1.1)
mb231cparp$Cell_line = "MB231"


tnbc_ddx54_cparp_bt549 <-subset(immunoblot_tnbc_ddx54, Gene == "c-PARP" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_cparp_bt549)
summary(bt549)

tnbc_ddx54_cparp_bt549$shRNA <- factor(as.factor(tnbc_ddx54_cparp_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549cparp <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_cparp_bt549)
bt549cparp <- as.data.frame(bt549cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549cparp$group1 <- out$V1
bt549cparp$group2 <- out$V2
bt549cparp$pval <- makeStars(bt549cparp$pval)
bt549cparp$y.position <- max(tnbc_ddx54_cparp_bt549$fold_change*1.1)
bt549cparp$Cell_line = "BT549"




tnbccparp <- rbind(mb231cparp, bt549cparp)
tnbccparp $Cell_line <- factor(tnbccparp$Cell_line, levels = c("MB231", "BT549"))
tnbccparp <- subset(tnbccparp, !pval == "ns")

ggplot(subset(immunoblot_tnbc_ddx54, Gene == "c-PARP"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + 
  scale_colour_manual(values = pallete3, 
                      labels = c(BT549 = "BT-549", MB231 = "MDA-MB-231")) +
  geom_col(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "c-PARP", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "c-PARP", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbccparp, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/cparp_tnbc_ddx54_facet.tiff", height =5, width = 5, units = "in")




tnbc_ddx54_peif2a_mb231 <-subset(immunoblot_tnbc_ddx54, Gene == "p-eIF2a" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_peif2a_mb231)
summary(mb231)

tnbc_ddx54_peif2a_mb231$shRNA <- factor(as.factor(tnbc_ddx54_peif2a_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231peif2a <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_peif2a_mb231)
mb231peif2a <- as.data.frame(mb231peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231peif2a$group1 <- out$V1
mb231peif2a$group2 <- out$V2
mb231peif2a$pval <- makeStars(mb231peif2a$pval)
mb231peif2a$y.position <- max(tnbc_ddx54_peif2a_mb231$fold_change*1.1)
mb231peif2a$Cell_line = "MB231"


tnbc_ddx54_peif2a_bt549 <-subset(immunoblot_tnbc_ddx54, Gene == "p-eIF2a" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_peif2a_bt549)
summary(bt549)

tnbc_ddx54_peif2a_bt549$shRNA <- factor(as.factor(tnbc_ddx54_peif2a_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549peif2a <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_peif2a_bt549)
bt549peif2a <- as.data.frame(bt549peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549peif2a$group1 <- out$V1
bt549peif2a$group2 <- out$V2
bt549peif2a$pval <- makeStars(bt549peif2a$pval)
bt549peif2a$y.position <- max(tnbc_ddx54_peif2a_bt549$fold_change*1.1)
bt549peif2a$Cell_line = "BT549"




tnbcpeif2a <- rbind(mb231peif2a, bt549peif2a)
tnbcpeif2a $Cell_line <- factor(tnbcpeif2a$Cell_line, levels = c("MB231", "BT549"))
tnbcpeif2a <- subset(tnbcpeif2a, !pval == "ns")

ggplot(subset(immunoblot_tnbc_ddx54, Gene == "p-eIF2a"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + 
  labs(x = "", y = "Relative p-eIF2A/eIF2A Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + 
  scale_colour_manual(values = pallete3, 
                      labels = c(BT549 = "BT-549", MB231 = "MDA-MB-231")) +
  geom_col(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "p-eIF2a", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "p-eIF2a", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcpeif2a, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/peif2a_tnbc_ddx54_facet.tiff", height =5, width = 5, units = "in")


tnbc_ddx54_adar_mb231 <-subset(immunoblot_tnbc_ddx54, Gene == "ADAR" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_adar_mb231)
summary(mb231)

tnbc_ddx54_adar_mb231$shRNA <- factor(as.factor(tnbc_ddx54_adar_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231adar <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_adar_mb231)
mb231adar <- as.data.frame(mb231adar$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231adar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231adar$group1 <- out$V1
mb231adar$group2 <- out$V2
mb231adar$pval <- makeStars(mb231adar$pval)
mb231adar$y.position <- max(tnbc_ddx54_adar_mb231$fold_change*1.1)
mb231adar$Cell_line = "MB231"


tnbc_ddx54_adar_bt549 <-subset(immunoblot_tnbc_ddx54, Gene == "ADAR" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_adar_bt549)
summary(bt549)

tnbc_ddx54_adar_bt549$shRNA <- factor(as.factor(tnbc_ddx54_adar_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549adar <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_adar_bt549)
bt549adar <- as.data.frame(bt549adar$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549adar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549adar$group1 <- out$V1
bt549adar$group2 <- out$V2
bt549adar$pval <- makeStars(bt549adar$pval)
bt549adar$y.position <- max(tnbc_ddx54_adar_bt549$fold_change*1.1)
bt549adar$Cell_line = "BT549"




tnbcadar <- rbind(mb231adar, bt549adar)
tnbcadar$Cell_line <- factor(tnbcadar$Cell_line, levels = c("MB231", "BT549"))
tnbcadar <- subset(tnbcadar, !pval == "ns")

ggplot(subset(immunoblot_tnbc_ddx54, Gene == "ADAR"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + 
  labs(x = "", y = "Relative ADAR Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + 
  scale_colour_manual(values = pallete3, 
                      labels = c(BT549 = "BT-549", MB231 = "MDA-MB-231")) +
  geom_col(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "ADAR", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "ADAR", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcadar, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/adar_tnbc_ddx54_facet.tiff", height =5, width = 5, units = "in")




tnbc_ddx54_pkr_mb231 <-subset(immunoblot_tnbc_ddx54, Gene == "total PKR" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_pkr_mb231)
summary(mb231)

tnbc_ddx54_pkr_mb231$shRNA <- factor(as.factor(tnbc_ddx54_pkr_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231pkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_pkr_mb231)
mb231pkr <- as.data.frame(mb231pkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231pkr$group1 <- out$V1
mb231pkr$group2 <- out$V2
mb231pkr$pval <- makeStars(mb231pkr$pval)
mb231pkr$y.position <- max(tnbc_ddx54_pkr_mb231$fold_change*1.1)
mb231pkr$Cell_line = "MB231"


tnbc_ddx54_pkr_bt549 <-subset(immunoblot_tnbc_ddx54, Gene == "total PKR" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_pkr_bt549)
summary(bt549)

tnbc_ddx54_pkr_bt549$shRNA <- factor(as.factor(tnbc_ddx54_pkr_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549pkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_pkr_bt549)
bt549pkr <- as.data.frame(bt549pkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549pkr$group1 <- out$V1
bt549pkr$group2 <- out$V2
bt549pkr$pval <- makeStars(bt549pkr$pval)
bt549pkr$y.position <- max(tnbc_ddx54_pkr_bt549$fold_change*1.1)
bt549pkr$Cell_line = "BT549"




tnbcpkr <- rbind(mb231pkr, bt549pkr)
tnbcpkr$Cell_line <- factor(tnbcpkr$Cell_line, levels = c("MB231", "BT549"))
tnbcpkr <- subset(tnbcpkr, !pval == "ns")

ggplot(subset(immunoblot_tnbc_ddx54, Gene == "total PKR"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + 
  labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + 
  scale_colour_manual(values = pallete3, 
                      labels = c(BT549 = "BT-549", MB231 = "MDA-MB-231")) +
  geom_col(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "total PKR", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "total PKR", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcpkr, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/pkr_tnbc_ddx54_facet.tiff", height =5, width = 5, units = "in")





tnbc_ddx54_ddx54_mb231 <-subset(immunoblot_tnbc_ddx54, Gene == "DDX54" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_ddx54_mb231)
summary(mb231)

tnbc_ddx54_ddx54_mb231$shRNA <- factor(as.factor(tnbc_ddx54_ddx54_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231ddx54 <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_ddx54_mb231)
mb231ddx54 <- as.data.frame(mb231ddx54$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231ddx54)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231ddx54$group1 <- out$V1
mb231ddx54$group2 <- out$V2
mb231ddx54$pval <- makeStars(mb231ddx54$pval)
mb231ddx54$y.position <- max(tnbc_ddx54_ddx54_mb231$fold_change*1.1)
mb231ddx54$Cell_line = "MB231"


tnbc_ddx54_ddx54_bt549 <-subset(immunoblot_tnbc_ddx54, Gene == "DDX54" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_ddx54_ddx54_bt549)
summary(bt549)

tnbc_ddx54_ddx54_bt549$shRNA <- factor(as.factor(tnbc_ddx54_ddx54_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549ddx54 <- DunnettTest(fold_change ~ shRNA, data = tnbc_ddx54_ddx54_bt549)
bt549ddx54 <- as.data.frame(bt549ddx54$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549ddx54)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ddx54$group1 <- out$V1
bt549ddx54$group2 <- out$V2
bt549ddx54$pval <- makeStars(bt549ddx54$pval)
bt549ddx54$y.position <- max(tnbc_ddx54_ddx54_bt549$fold_change*1.1)
bt549ddx54$Cell_line = "BT549"




tnbcddx54 <- rbind(mb231ddx54, bt549ddx54)
tnbcddx54$Cell_line <- factor(tnbcddx54$Cell_line, levels = c("MB231", "BT549"))
tnbcddx54 <- subset(tnbcddx54, !pval == "ns")

ggplot(subset(immunoblot_tnbc_ddx54, Gene == "DDX54"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + 
  labs(x = "", y = "Relative DDX54 Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + 
  scale_colour_manual(values = pallete3, 
                      labels = c(BT549 = "BT-549", MB231 = "MDA-MB-231")) +
  geom_col(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "DDX54", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_ddx54_sum, Gene == "DDX54", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcddx54, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/ddx54_tnbc_ddx54_facet.tiff", height =5, width = 5, units = "in")







#read ff data for knockdowns
ff_mcf7_adar_ddx54 <- fread("DDX54 Foci Formation/ff_mcf7_adar_ddx54.txt")

ff_mcf7_adar_ddx54 <- na.omit(ff_mcf7_adar_ddx54)

ff_mcf7_adar_ddx54$sample <- paste(ff_mcf7_adar_ddx54$shRNA_1, ff_mcf7_adar_ddx54$shRNA_2, sep = "/")

ff_mcf7_adar_ddx54$sample <- gsub("-", ".", ff_mcf7_adar_ddx54$sample)

#group and summarise data
ff_mcf7_adar_ddx54_sum <- dplyr::group_by(ff_mcf7_adar_ddx54, sample)

ff_mcf7_adar_ddx54_sum <- dplyr::summarise(ff_mcf7_adar_ddx54_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_mcf7_adar_ddx54_sum$ci <- qnorm(0.975)*(ff_mcf7_adar_ddx54_sum$sd_area/sqrt(4))

aov_ff <- aov(Relative_area ~ sample, ff_mcf7_adar_ddx54)
tff <- TukeyHSD(aov_ff)
tff <- as.data.frame(tff$sample[,1:4])

out <- strsplit(as.character(row.names(tff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tff$stars <- makeStars(tff$`p adj`)

tff$xmin <- out$V1
tff$xmax <- out$V2
tff$y <- max(ff_mcf7_adar_ddx54$Relative_area*1.1)

tff <- subset(tff, !stars == "ns")

ggplot(ff_mcf7_adar_ddx54, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR1/shSCR", "shADAR1/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_mcf7_adar_ddx54_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_mcf7_adar_ddx54_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) +
  geom_bracket(xmin = tff$xmin, 
               xmax = tff$xmax, 
               y.position = tff$y , label = tff$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("FF_plots/adar_ddx54_mcf7.tiff", height = 4.2, width = 4.2, units = "in")


aov_ff <- NULL
tff <- NULL



#read ff data for knockdowns
ff_skbr3_adar_ddx54 <- fread("DDX54 Foci Formation/ff_skbr3_adar_ddx54.txt")

ff_skbr3_adar_ddx54 <- na.omit(ff_skbr3_adar_ddx54)

ff_skbr3_adar_ddx54$sample <- paste(ff_skbr3_adar_ddx54$shRNA_1, ff_skbr3_adar_ddx54$shRNA_2, sep = "/")

ff_skbr3_adar_ddx54$sample <- gsub("-", ".", ff_skbr3_adar_ddx54$sample)

#group and summarise data
ff_skbr3_adar_ddx54_sum <- dplyr::group_by(ff_skbr3_adar_ddx54, sample)

ff_skbr3_adar_ddx54_sum <- dplyr::summarise(ff_skbr3_adar_ddx54_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_skbr3_adar_ddx54_sum$ci <- qnorm(0.975)*(ff_skbr3_adar_ddx54_sum$sd_area/sqrt(4))

aov_ff <- aov(Relative_area ~ sample, ff_skbr3_adar_ddx54)
tff <- TukeyHSD(aov_ff)
tff <- as.data.frame(tff$sample[,1:4])

out <- strsplit(as.character(row.names(tff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tff$stars <- makeStars(tff$`p adj`)

tff$xmin <- out$V1
tff$xmax <- out$V2
tff$y <- max(ff_skbr3_adar_ddx54$Relative_area*1.1)

tff <- tff[!grepl("shDDX54.5", rownames(tff)), ]

tff <- subset(tff, !stars == "ns")

ggplot(ff_skbr3_adar_ddx54, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX54.4", "shADAR/shSCR", "shADAR/shDDX54.4"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX54-4", "shADAR1\nshSCR", "shADAR1\nshDDX54-4")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_skbr3_adar_ddx54_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_skbr3_adar_ddx54_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) +
  geom_bracket(xmin = tff$xmin, 
               xmax = tff$xmax, 
               y.position = tff$y , label = tff$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("FF_plots/adar_ddx54_skbr3.tiff", height = 4.2, width = 4.2, units = "in")


aov_ff <- NULL
tff <- NULL



#read ff data for knockdowns

ff_tnbc_ddx54 <- fread("tnbc_ddx54_ff.txt")


ff_tnbc_ddx54$shRNA <- gsub("-", ".", ff_tnbc_ddx54$shRNA)

#group and summarise data
ff_tnbc_ddx54_sum <- dplyr::group_by(ff_tnbc_ddx54, Cell_line, shRNA)

ff_tnbc_ddx54_sum <- dplyr::summarise(ff_tnbc_ddx54_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_tnbc_ddx54_sum$ci_area <- qnorm(0.975)*(ff_tnbc_ddx54_sum$sd_area/sqrt(3))

tnbc_ddx54_ff_mb231 <-subset(ff_tnbc_ddx54, Cell_line == "MB231")
mb231 <- aov(Relative_area ~ shRNA, data = tnbc_ddx54_ff_mb231)
summary(mb231)

tnbc_ddx54_ff_mb231$shRNA <- factor(as.factor(tnbc_ddx54_ff_mb231$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

mb231ff <- DunnettTest(Relative_area ~ shRNA, data = tnbc_ddx54_ff_mb231)
mb231ff <- as.data.frame(mb231ff$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231ff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231ff$group1 <- out$V1
mb231ff$group2 <- out$V2
mb231ff$pval <- makeStars(mb231ff$pval)
mb231ff$y.position <- max(tnbc_ddx54_ff_mb231$Relative_area*1.1)
mb231ff$Cell_line = "MB231"


tnbc_ddx54_ff_bt549 <-subset(ff_tnbc_ddx54, Cell_line == "BT549")
bt549 <- aov(Relative_area ~ shRNA, data = tnbc_ddx54_ff_bt549)
summary(bt549)

tnbc_ddx54_ff_bt549$shRNA <- factor(as.factor(tnbc_ddx54_ff_bt549$shRNA), levels = c("shSCR", "shDDX54.4", "shDDX54.5"))

bt549ff <- DunnettTest(Relative_area ~ shRNA, data = tnbc_ddx54_ff_bt549)
bt549ff <- as.data.frame(bt549ff$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549ff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ff$group1 <- out$V1
bt549ff$group2 <- out$V2
bt549ff$pval <- makeStars(bt549ff$pval)
bt549ff$y.position <- max(tnbc_ddx54_ff_bt549$Relative_area*1.1)
bt549ff$Cell_line = "BT549"




tnbcff <- rbind(mb231ff, bt549ff)
tnbcff <- subset(tnbcff, !pval == "ns")

ggplot(ff_tnbc_ddx54, aes(shRNA, Relative_area, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDDX54.4", "shDDX54.5"), labels = c("shSCR", "shDDX54-4", "shDDX54-5")) + labs(x = "", y = "Relative Foci Area", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = ff_tnbc_ddx54_sum, aes(shRNA, mean_area, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = ff_tnbc_ddx54_sum, aes(x = shRNA, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area, colour = Cell_line), width=0.2, position = position_dodge(width = 0.7)) + 
  stat_pvalue_manual(data = tnbcff, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("FF_plots/tnbc_ddx54_facet.tiff", height =4.5, width = 4, units = "in")





######BT549 Rescue

#read western (immunoblot) data for knockdowns
immunoblot_BT549_rescue <- fread("immunoblot_BT549_rescue.txt")

immunoblot_BT549_rescue$V6 <- NULL

immunoblot_BT549_rescue <- na.omit(immunoblot_BT549_rescue)


immunoblot_BT549_rescue$sample <- paste(immunoblot_BT549_rescue$Overexpression, immunoblot_BT549_rescue$shRNA, sep = "/")

immunoblot_BT549_rescue$sample <- gsub("-", ".", immunoblot_BT549_rescue$sample)

immunoblot_BT549_rescue$fold_change <- as.numeric(immunoblot_BT549_rescue$fold_change)

#group and summarise data
immunoblot_BT549_rescue_sum <- dplyr::group_by(immunoblot_BT549_rescue, sample, Gene)

immunoblot_BT549_rescue_sum <- dplyr::summarise(immunoblot_BT549_rescue_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_BT549_rescue_sum$ci <- qnorm(0.975)*(immunoblot_BT549_rescue_sum$sd_expression/sqrt(3))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_BT549_rescue, Gene == "pPKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_BT549_rescue$fold_change[immunoblot_BT549_rescue$Gene == "pPKR"]*1.1)

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_BT549_rescue, Gene == "pPKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("EV/shSCR", "EV/shDDX54", "DDX54/shSCR", "DDX54/shDDX54"),
                   labels = c("EV\nshSCR", "EV\nshDDX54", "DDX54\nshSCR", "DDX54\nshDDX54")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_BT549_rescue_sum, Gene == "pPKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_BT549_rescue_sum, Gene == "pPKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_ddx54_rescue.tiff", height = 4.2, width = 4.2)


pkr_kinase_df <- read.csv("pkr_kinase_df.csv")

pkr_kinase_df_sum <- dplyr::group_by(pkr_kinase_df, shRNA, RNA)

pkr_kinase_df_sum <- dplyr::summarise(pkr_kinase_df_sum, mean_signal = mean(Signal, na.rm=T), sd_signal = sd(Signal, na.rm=T))

ggplot(pkr_kinase_df, aes(RNA, Signal, colour = shRNA)) + geom_point(aes(shape = as.factor(Replicate)), alpha = 0.5) +
  geom_line(data = pkr_kinase_df_sum, aes(RNA, mean_signal, colour = shRNA)) + 
  geom_pointrange(data = pkr_kinase_df_sum, aes(RNA, y = mean_signal, ymin = mean_signal-sd_signal, 
                                                ymax = mean_signal + sd_signal,
                                                colour = shRNA)) +
  theme_science() + scale_colour_manual(values = pallete3) + guides(shape = "none") + 
  labs(x = "RNA (nM)", y = "% PKR Activation") + theme(legend.position = c(0.5,0.7), legend.title=element_blank())
ggsave("pkr_activation.tiff", height = 3, width = 5)


aov_ppkr <- NULL
tppkr <- NULL


#read cell_proliferation data for knockdowns
cell_proliferation_mb231 <- fread("Cell_counts.txt")

cell_proliferation_mb231 <- na.omit(cell_proliferation_mb231)

#group and summarise data
cell_proliferation_mb231_sum <- dplyr::group_by(cell_proliferation_mb231, Sample, Timepoint)

cell_proliferation_mb231_sum <- dplyr::summarise(cell_proliferation_mb231_sum, mean_count = mean(Cell_count), sd_count = sd(Cell_count))

ggplot(cell_proliferation_mb231_sum, aes(Timepoint, mean_count, colour = Sample)) +
  geom_line() + theme_science(base_size = 13)  +
  labs(x = "Days after plating", y = "Cell Count", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), 
        legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(shape = "none") + scale_colour_manual(limits = c("EV-Scr", "EV-sh54", "54-Scr", "54-sh54"),
                                            labels = c("EV\nshSCR", "EV\nshDDX54", "DDX54\nshSCR", "DDX54\nshDDX54"), values = pallete3) +
  geom_errorbar(aes(x = Timepoint, y = mean_count, ymin = mean_count - sd_count, ymax = mean_count + sd_count), width=0.2) +
  geom_point(data =  cell_proliferation_mb231, aes(Timepoint, Cell_count, colour = Sample, shape = as.factor(Replicate)), alpha = 0.5)
ggsave("cell_proliferation_mb231.tiff", height = 3, width = 6, units = "in")



######MB231 Rescue

#read western (immunoblot) data for knockdowns
immunoblot_MB231_rescue <- fread("revisionexp.csv")

immunoblot_MB231_rescue <- na.omit(immunoblot_MB231_rescue)


immunoblot_MB231_rescue$sample <- paste(immunoblot_MB231_rescue$Overexpression, immunoblot_MB231_rescue$shRNA, sep = "/")

immunoblot_MB231_rescue$sample <- gsub("-", ".", immunoblot_MB231_rescue$sample)

immunoblot_MB231_rescue$fold_change <- as.numeric(immunoblot_MB231_rescue$fold_change)

#group and summarise data
immunoblot_MB231_rescue_sum <- dplyr::group_by(immunoblot_MB231_rescue, sample, Gene)

immunoblot_MB231_rescue_sum <- dplyr::summarise(immunoblot_MB231_rescue_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_MB231_rescue_sum$ci <- qnorm(0.975)*(immunoblot_MB231_rescue_sum$sd_expression/sqrt(3))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_MB231_rescue, Gene == "pPKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_MB231_rescue$fold_change[immunoblot_MB231_rescue$Gene == "pPKR"]*1.1)

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_MB231_rescue, Gene == "pPKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("EV/shSCR", "EV/shDDX54", "DDX54/shSCR", "DDX54/shDDX54"),
                   labels = c("EV\nshSCR", "EV\nshDDX54", "DDX54\nshSCR", "DDX54\nshDDX54")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_MB231_rescue_sum, Gene == "pPKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_MB231_rescue_sum, Gene == "pPKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_ddx54_rescue_mb231.tiff", height = 4.2, width = 4.2)
