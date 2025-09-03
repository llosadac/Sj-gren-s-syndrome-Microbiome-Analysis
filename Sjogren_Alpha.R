#Sjogren-s-syndrome-Microbiome-Analysis

#Import libraries
library(MicrobiotaProcess)
library(phyloseq)
library(ggpubr)
library(tibble)
library(dplyr)
library(microViz)
library(vegan)
library(pairwiseAdonis)
library(rstatix)
library(data.table)
library(tidyr)
library(microbiome)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(scales)
library(vcd)
library(stats)
library(car)
library(Hmisc)
library(magrittr)
library(microbiomeMarker)
library(microbiomeutilities)
library(jeevanuDB)
library(ANCOMBC)
library(tidyverse)
library(nlme)
library(tidyverse)
library(compositions)
library(ANCOMBC)
library(mia)
library(dplyr)
library(tidyr)
library(microbiomeMarker)
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(kableExtra)
library(xtable)
library(ggpubr)
library(phyloseq)
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggpubr)
library(BiocManager)
library(DESeq2)
library(microbiome)
library(philr)
library(btools)
library(fantaxtic)
library(ampvis2)
library(tsnemicrobiota)
library(cowplot)
library(ggpubr)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library("phyloseq")
library("ggplot2")
library("readxl")
library("dplyr")
library(devtools)
library(SpiecEasi)
library(seqtime)
library("MicEco")
library("vegan")
library(phyloseq)
library(metagenomeSeq)
library(iNEXT)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
library(decontam); packageVersion("decontam")
library(paletteer)

setwd("C:/Users/Laura L/Documents/TG Microbiota Sjogren/DADA2_Outputs")

# Load data necessary to phyloseq object
seqtab_all_no_chimeras <- readRDS("dd2.seqtab_all_no_chimeras.rds")
ASV_seq_rdp_set18_all <- readRDS("ASV_seq_rdp_set18.rds")
filereport_read_run_PRJNA525566_tsv <- read.csv("Metadata.csv")

# Set row names
rownames(filereport_read_run_PRJNA525566_tsv) <- filereport_read_run_PRJNA525566_tsv$X

# Create Phyloseq object
Phy_obj_RDP_all <- merge_phyloseq(sample_data(filereport_read_run_PRJNA525566_tsv), 
                                  otu_table(seqtab_all_no_chimeras, taxa_are_rows = FALSE),
                                  tax_table(ASV_seq_rdp_set18_all))

ascaris <- Phy_obj_RDP_all

# Verificate zero abundance
any(taxa_sums(ascaris) == 0)
sum(taxa_sums(ascaris) == 0)

PS <- ascaris
summarize_phyloseq(PS)

# HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL)

# Transpose the ASVs that are in columns
otu_mat <- as(otu_table(PS), "matrix")
if (!taxa_are_rows(PS)) {
  otu_mat <- t(otu_mat)
}

# How many reads for off-target eukaryotes and archaea
by(rowSums(otu_mat), tax_table(PS)[, "Kingdom"], sum)

# Filtering: Eliminate "empty" samples
PS <- prune_samples(sample_sums(PS)>0, PS)

# Sample filtering: Filtering samples with low counts
PS <- prune_samples(sample_sums(PS)>=2000, PS)
nsamples(PS)

# Taxa filtering
PS<- subset_taxa(PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "NA"))
PS <- subset_taxa(PS, !Genus %in% c("Stenotrophomonas", "Stenotrophomonas 1", "Stenotrophomonas 2"))

# General check
plot_richness(PS, x= "group" , color = "group", measures = c("Observed","Simpson", "Shannon")) +
  geom_jitter(alpha = 0.005) +
  scale_color_manual(values = c("Control"= "#A6CEE3","pSS" = "#B2DF8A")) +
  xlab(" ") +
  labs(color="Group")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

# Calculate observed richness
div_df <- estimate_richness(PS, measures = c("Observed", "Shannon", "Simpson"))
meta_df <- data.frame(sample_data(PS))
div_df$group <- meta_df$group

# Mann-Whitney-Wilcoxon test
wilcox_results <- list(
  wilcox.test(Observed ~ group, data = div_df),
  wilcox.test(Shannon ~ group, data = div_df),
  wilcox.test(Simpson ~ group, data = div_df)
)
wilcox_results

# Remove low prevalent taxa
# Create a prevalence dataframe
phyla2Filter<- c("Euryarchaeota") 
PS<- subset_taxa(PS, !Phylum %in% phyla2Filter)
Prevdf<- apply(X = otu_table(PS),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
# Verificate orientation
otu_mat <- otu_table(PS)
if (!taxa_are_rows(PS)) {
  otu_mat <- t(otu_mat)
}

# Calculate prevalence by taxon
Prevdf <- apply(X = otu_mat, MARGIN = 1, FUN = function(x){sum(x > 0)})

# Data.frame with prevalence, total abundance, and taxonomy
Prevdf <- data.frame(
  Prevalence = Prevdf,
  TotalAbundance = taxa_sums(PS),
  tax_table(PS)
)

plyr::ddply(Prevdf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

# Prevalence plot by phylum
ggplot(Prevdf, aes(TotalAbundance, Prevalence / nsamples(PS),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Log 10 Total Reads") + ylab("Prevalence [Prop. of Samples]") +
  theme_bw()+
  facet_wrap(~Phylum) + theme(legend.position="none")

# Normalization of proportions
# Normalization transformation to an even sample size
PS.Norm<- transform_sample_counts(PS, function(x) 1E6 * x/sum(x))

# Normalization transformation to an even sample size
PS.Norm<- transform_sample_counts(PS, function(x) 1E6 * x/sum(x))

# Check how many samples ended after filtering
table(PS@sam_data$group)

# Merge ASVs that have the same taxonomy at a certain taxonomic rank
PS.Fam<-  tax_glom(PS, "Family", NArm = T)
PS.Gen<-  tax_glom(PS, "Genus", NArm = T)
PS.Phy<-  tax_glom(PS, "Phylum", NArm = T)
PS.Spe<-  tax_glom(PS, "Species", NArm = T)

# Transform to relative abundance
PS.Spe.rel <- transform_sample_counts(PS.Spe, function(x) x / sum(x))
PS.Gen.rel <- transform_sample_counts(PS.Gen, function(x) x / sum(x))
PS.Fam.rel <- transform_sample_counts(PS.Fam, function(x) x / sum(x))
PS.Phy.rel <- transform_sample_counts(PS.Phy, function(x) x / sum(x))

# Plot relative abundance
plot_bar(PS.Phy.rel, fill = "Phylum") +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(legend.position = "right") +
  ylab("Abundancia relativa")

# Top 10 relative abundance by Phylum
abundance_phylum <- tax_glom(PS.Phy.rel, taxrank = "Phylum")
abundance_totalp <- taxa_sums(abundance_phylum)
top10_phylum <- names(sort(abundance_totalp, decreasing = TRUE))[1:10]
PS.top_phylum <- prune_taxa(top10_phylum, abundance_phylum)
PS.top_phylum_plot <- transform_sample_counts(PS.top_phylum, function(x) 100* x / sum(x))

plot_bar(PS.top_phylum_plot, fill = "Phylum")+
  facet_wrap(~group, scales = "free_x", ncol=1) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(12)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90, size=5)) +
  ylab("Relative Abundace (%)")

# Top 10 relative abundance by Genus
abundance_genus <- tax_glom(PS.Gen.rel, taxrank = "Genus")
abundance_totalg <- taxa_sums(abundance_genus)
top10_genus <- names(sort(abundance_totalg, decreasing = TRUE))[1:10]
PS.top_genus <- prune_taxa(top10_genus, abundance_genus)
PS.top_genus_plot <- transform_sample_counts(PS.top_genus, function(x) 100* x / sum(x))
plot_bar(PS.top_genus_plot, fill = "Genus")+
  facet_wrap(~group, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(12))+
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90, size=6)) +
  ylab("Relative Abundace (%)")

# Estimate alpha diversity for individual samples
 alphadiv<- estimate_richness(PS)
 
# Add sample data into a single data frame
as.data.frame(PS@sam_data)->tmp
alphadiv<-cbind(alphadiv, tmp)
table(alphadiv$group)

saveRDS(PS, "C:/Users/Laura L/Documents/TG Microbiota Sjogren/DADA2_Outputs/phyloseq_outputs/PS.PA.Rds") ##Total
saveRDS(PS.Norm, "C:/Users/Laura L/Documents/TG Microbiota Sjogren/DADA2_Outputs/phyloseq_outputs/PS.Norm.Rds") ##Normalized
saveRDS(alphadiv, "C:/Users/Laura L/Documents/TG Microbiota Sjogren/DADA2_Outputs/phyloseq_outputs/alphadiv.rds") ##Alpha diverisity tables


## ALPHA AND BETA DIVERSITY

# Check packages 
library(phyloseq)
library(microbiome)
library(tidyverse)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(ggsci)
library(microbiome)
library(lme4)
library(circlize)

setwd("C:/Users/Laura L/Documents/TG Microbiota Sjogren/DADA2_Outputs/phyloseq_outputs")

# Load data
PS.PA<- readRDS("PS.PA.Rds")  #PS original
PS.PA.Norm<- readRDS("PS.Norm.Rds") #PS Normalizado
alphadiv.PA.rare<- readRDS("alphadiv.rds") # Alpha diversity tables with sample information

# Alpha diversity by groups

#Statistical evaluation of diferences between groups, regarging Shannon Index
alphadiv.PA.rare%>%
  wilcox_test(Shannon ~ group)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "group") %>% print(n = 30)

# Plot Observed
alphadiv.PA.rare%>%
  mutate(group = factor(group, levels = c("pSS", "Control")))%>%
  ggplot(aes(x= group, y= Observed, color= group, fill= group))+
  geom_boxplot(outlier.shape=NA, width = 0.6)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("Control" = "#030303", "pSS" = "#030303"))+
  scale_fill_manual(values = c("Control"= "#A6CEE3","pSS" = "#B2DF8A"))+
  xlab("")+
  ylab("ASV Diversity (Observed Index)")+
  labs(fill= "Group")+
  theme_bw() + 
  theme(legend.position="none") +

# Plot Shannon
alphadiv.PA.rare%>%
  mutate(group = factor(group, levels = c("pSS", "Control")))%>%
  ggplot(aes(x= group, y= Shannon, color= group, fill= group))+
  geom_boxplot(outlier.shape=NA, width = 0.6)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("Control" = "#030303", "pSS" = "#030303"))+
  scale_fill_manual(values = c("Control"= "#A6CEE3","pSS" = "#B2DF8A"))+
  xlab("")+
  ylab("ASV Diversity (Shannon Index)")+
  labs(fill= "Group")+
  theme_bw() + 
  theme(legend.position="none") +

# Plot simpson
alphadiv.PA.rare%>%
  mutate(group = factor(group, levels = c("pSS", "Control")))%>%
  ggplot(aes(x= group, y= Simpson, color= group, fill= group))+
  geom_boxplot(outlier.shape=NA, width = 0.6)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("Control" = "#030303", "pSS" = "#030303"))+
  scale_fill_manual(values = c("Control"= "#A6CEE3","pSS" = "#B2DF8A"))+
  xlab("")+
  ylab("ASV Diversity (Simpson Index)")+
  labs(fill= "Group")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  theme_bw()

# Beta diversity by groups

set.seed(1024)

# Dissimilarity metric
dist_bray = phyloseq::distance(PS, method="bray")

# Distance matrix
ordination = ordinate(PS, method="NMDS", distance=dist_bray, trymax=100)

# Plot NMDS by Bray-Curtis
plot_ordination(PS, ordination, color = "group") +
  geom_point(size = 3, shape = 21, stroke = 1, aes(fill = group, color = group)) +
  scale_fill_manual(values = c("pSS" = "#B2DF8A", "Control" = "#A6CEE3")) +
  scale_color_manual(values = c("pSS" = "#66A94D", "Control" = "#5A9ACF")) +
  stat_ellipse(aes(group = group), type = "norm", linetype = 2) +
  labs(title = "Beta Diversity NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ordination$stress, 3)),
       fill = "Group", color = "Group") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.1) )

# Permanova test
perm_braynmds <- adonis2(dist_bray ~ group, data = data.frame(sample_data(PS)), permutations = 999)
perm_braynmds

# Dispersion test
disp <- betadisper(dist_bray, group = sample_data(PS)$group)
permutest(disp)

# Differential abundance analysis Lefse

lefse_model <- run_lefse(
  PS.PA.Norm,
  taxa_rank = "all",
  wilcoxon_cutoff = 0.01,
  group = "group",
  multigrp_strat = FALSE,
  lda_cutoff = 4
)
lefse_model

marker_list <- lefse_model@marker_table@.Data
marker_names <- lefse_model@marker_table@names
lefse_df <- setNames(as.data.frame(marker_list), marker_names)
write.csv(lefse_df, file = "lefse_resultados.csv", row.names = FALSE)

# Plot lefse at the ASV level
# Barplot
p_bar <- plot_ef_bar(lefse_model)
p_bar + scale_fill_manual(values = c("pSS" = "#B2DF8A", "Control" = "#A6CEE3")) +
  theme_bw()

# Additional graphics
plot_ef_dot(lefse_model)
plot_abundance(lefse_model, group = "group")

# Heatmap
group_colors <- c("pSS" = "#B2DF8A", "Control" = "#A6CEE3")

plot_heatmap(lefse_model,
             transform = "log10p",
             group = "group",
             cluster_marker = T,
             cluster_sample = F,
             col = brewer.pal(9,"RdPu"),
             border = "black",
             annotation_col = group_colors,
             row_names_gp = gpar(fontsize = 9),
             height = unit(12, "cm"))
             
