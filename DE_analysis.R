# Load libraries----
library(tidyverse)
library(tximport)
library(gt)
library(edgeR)
library(matrixStats)
library(cowplot)
library(DT)
library(plotly)
library(wesanderson)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(pheatmap)

# Importing Kallisto mapping data into R----

# read in study design file
targets <- read_tsv("studydesign.txt")
# set file paths to your Kallisto output folders that contain quantification data
files <- file.path("../read_mapping", targets$long_name, "abundance.tsv")
all(file.exists(files))

# use TxImport package to read Kallisto data into R
txi_kallisto <- tximport(files, 
                         type = "kallisto",
                         txOut = TRUE, # represent your data at transcript level
                         countsFromAbundance = "lengthScaledTPM")

# save the resulting R data object for later use
save(txi_kallisto, file = "txi_kallisto")

# capture variables of interest from the study design
condition <- as.factor(targets$condition)
condition <- factor(condition, levels = c("WT_Lac", "Mut_Lac", "WT_LNnT", "Mut_LNnT"))
batch <- as.factor(targets$batch)
strain <- as.factor(targets$strain)
carb <- as.factor(targets$carb)
# capture sample labels for later use
sampleLabels <- targets$sample

# use gt package to produce table of study design
gt(targets) %>% 
  cols_align(
    align = "left",
    columns = everything()
  )

# Create the DGElist (count matrix)----

myDGEList <- DGEList(txi_kallisto$counts)
# use the 'cpm' function from EdgeR to get counts per million
log2.cpm <- cpm(myDGEList, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = WT_Lac1:Mut_LNnT3,
                                  # column names to be stored as a SINGLE variable
                                  names_to = "samples", 
                                  # name of that new variable (column)
                                  values_to = "expression") 
                                  # name of new variable (column)
                                  # storing all the values (data)

# plot unfiltered, non-normalized CPM
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 124, 
               size = 5.5, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Unfiltered, non-normalized") +
  coord_flip() +
  theme_bw()

# Filter counts and plot filtered, non-normalized CPM----

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=3 # only keep genes that have cpm>1 (== not zeroes) in more than 3 samples (minimal group size)
myDGEList.filtered <- myDGEList[keepers,]
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, 
                                           # dataframe to be pivoted
                                           cols = WT_Lac1:Mut_LNnT3,
                                           # column names to be stored as a SINGLE variable
                                           names_to = "samples",
                                           # name of that new variable (column)
                                           values_to = "expression")
                                           # name of new variable (column) storing all the values (data)

# plot filtered, non-normalized CPM
p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 124, 
               size = 5.5, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized") +
  coord_flip() +
  theme_bw()

# Normalize counts and plot filtered, normalized CPM----

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                # dataframe to be pivoted
                                                cols = WT_Lac1:Mut_LNnT3,
                                                # column names to be stored as a SINGLE variable
                                                names_to = "samples",
                                                # name of that new variable (column)
                                                values_to = "expression") 
                                                # name of new variable (column) 
                                                # storing all the values (data)

# plot of signal distribution again to see effect of normalization
p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 124, 
               size = 5.5, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Filtered, TMM normalized") +
  coord_flip() +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

# Perform PCA and plot the results----

# running PCA
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1)
# converting PCA result into a tibble for plotting
pca.res.df <- as_tibble(pca.res$x)

# plotting PCA
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = carb, shape = strain) +
  geom_point(size=4) +
  scale_color_manual(name = "Carbon source",
                     breaks=c("Lac","LNnT"), 
                     values=c("slateblue", "sienna2"), 
                     labels=c("Lactose", "LNnT")) +
  scale_shape_manual(name = "Strain",
                     breaks=c("JCM1222","M3"),
                     values=c(16, 17),
                     labels=c("WT", "ΔnagR")) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title= "PCA of B. infantis JCM1222: WT vs ΔnagR",
       subtitle = "Principal component analysis (PCA) showing clear separation \nbetween growth on Lac and LNnT and between WT and ΔNagR",
       color = "carb", shape="strain") +
  coord_fixed(ratio=1.2) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"))

# save the figure in the eps format
#ggsave("pca.eps", width = 5, height = 5)

# create small multiples plot
pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels,
             condition = condition)
             #you can insert other groupings here based on your study design file
pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=condition) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  scale_fill_manual(name = "Condition",
                    breaks=c("WT_LNnT", "WT_Lac", "Mut_LNnT","Mut_Lac"), 
                    #values=c("slateblue", "sienna2", "red", "blue"),
                    values=c(wes_palette(n=4, name="Moonrise2")),
                    labels=c("WT LNnT", "WT Lac", "ΔnagR LNnT", "ΔnagR_Lac")) +
  labs(title="PCA small multiples plot of B. infantis JCM1222: WT vs ΔNagR") +
  theme_bw() +
  coord_flip()

# Perform DE analysis----

# setting up model matrix without an intercept
design <- model.matrix(~0 + condition)
colnames(design) <- levels(condition)
# using VOOM function from Limma package to apply precision weights to each gene
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
fit <- lmFit(v.DEGList.filtered.norm, design)
# setting up contrast matrix for two main pairwise comparisons
contrast.matrix <- makeContrasts(Mut_WT_Lac = Mut_Lac - WT_Lac,
                                 Mut_WT_LNnT = Mut_LNnT - WT_LNnT,
                                 LNnT_WT = WT_LNnT - WT_Lac,
                                 LNnT_Mut = Mut_LNnT - Mut_Lac,
                                 levels=design)
fits <- contrasts.fit(fit, contrast.matrix)
# extracting stats 
ebFit <- eBayes(fits)

# listing stats for all genes in the dataset to be used for making volcano plot
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=2600, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
# select only genes with significant logFC and adj.P.Val
myTopHits.df.de <- subset(myTopHits.df, (logFC > 1 | logFC < -1) & adj.P.Val < 0.01)
# save the DE data
#write_tsv(myTopHits.df, "DE_Mut_WT_Lac_all_genes.txt")
#write_tsv(myTopHits.df.de, "DE_Mut_WT_Lac.txt")

# Create a volcano plot----

# create a vector containing locus_tags of genes potentially controlled by NagR
targets.nagR <- c("Blon_0879", "Blon_0881", "Blon_0882", "Blon_0883", "Blon_0884", "Blon_0885", "Blon_2171", "Blon_2172","Blon_2173", "Blon_2174", "Blon_2175", "Blon_2176", "Blon_2177", "Blon_2331", 
                  "Blon_2332", "Blon_2341", "Blon_2342", "Blon_2343", "Blon_2344", "Blon_2345", "Blon_2346", "Blon_2347", "Blon_2348", "Blon_2349", "Blon_2350", "Blon_2351", "Blon_2352", "Blon_2354")
# create a vector containing locus_tags of genes potentially controlled by CscR
targets.cscR <- c("Blon_0789", "Blon_0788", "Blon_0787", "Blon_0786")
# subset volcano plot data based targets.nagR and targets.cscR
myTopHits.nagR <- subset(myTopHits.df, geneID %in% targets.nagR)
myTopHits.cscR <- subset(myTopHits.df, geneID %in% targets.cscR)
# subset volcano plot data labels (NagR-controlled genes)
myTopHits.df$nagR <- myTopHits.df$geneID
myTopHits.nagR_selected <- myTopHits.df$nagR %in% myTopHits.nagR$geneID
myTopHits.df$nagR[!myTopHits.nagR_selected] <- NA
# subset volcano plot data labels (CscR-controlled genes)
myTopHits.df$cscR <- myTopHits.df$geneID
myTopHits.cscR_selected <- myTopHits.df$cscR %in% myTopHits.cscR$geneID
myTopHits.df$cscR[!myTopHits.cscR_selected] <- NA

# Create a volcano plot----
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=3, alpha=.3) +
  geom_point(mapping=NULL, myTopHits.nagR, size = 3, shape = 21, fill= "sienna2", inherit.aes = TRUE) +
  geom_point(mapping=NULL, myTopHits.cscR , size = 3, shape = 21, fill= "slateblue", inherit.aes = TRUE) +
  geom_text_repel(aes(label = nagR), size = 3, fontface=2, color="black") +
  geom_text_repel(aes(label = cscR), size = 3, fontface=2, color="black") +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=0.6) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=0.6) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=0.6) +
  annotate("text", x=-6, y=-log10(0.01)+0.3,
           label=paste("Padj<0.01"), size=5, fontface="bold") +
  scale_x_continuous(limits=c(-7,5), breaks = -7:5) +
  labs(title="Volcano plot",
       subtitle = "B. infantis JCM1222 grown on lactose: ΔnagR vs WT") +
  theme(plot.title = element_text(face="bold")) +
  theme_cowplot() + 
  coord_fixed(0.55)
  
ggplotly(vplot)


# listing stats for all genes in the dataset to be used for making volcano plot
myTopHits3 <- topTable(ebFit, adjust ="BH", coef=3, number=2600, sort.by="logFC")
myTopHits.df3 <- myTopHits3 %>%
  as_tibble(rownames = "geneID")
# select only genes with significant logFC and adj.P.Val
myTopHits.df3.de <- subset(myTopHits.df3, (logFC > 1 | logFC < -1) & adj.P.Val < 0.01)

# save the DE data
#write_tsv(myTopHits.df3, "DE_WT_LNnT_vs_Lac_all_genes.txt")
write_tsv(myTopHits.df3.de, "DE_WT_LNnT_vs_Lac.txt")

# subset volcano plot data based targets.nagR and targets.cscR
myTopHits.nagR3 <- subset(myTopHits.df3, geneID %in% targets.nagR)
myTopHits.cscR3 <- subset(myTopHits.df3, geneID %in% targets.cscR)
# subset volcano plot data labels (NagR-controlled genes)
myTopHits.df3$nagR <- myTopHits.df3$geneID
myTopHits.nagR_selected3 <- myTopHits.df3$nagR %in% myTopHits.nagR3$geneID
myTopHits.df3$nagR[!myTopHits.nagR_selected3] <- NA
# subset volcano plot data labels (CscR-controlled genes)
myTopHits.df3$cscR <- myTopHits.df3$geneID
myTopHits.cscR_selected3 <- myTopHits.df3$cscR %in% myTopHits.cscR3$geneID
myTopHits.df3$cscR[!myTopHits.cscR_selected3] <- NA

# create a volcano plot
vplot2 <- ggplot(myTopHits.df3) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=3, shape = 16, color = "black", alpha =.3) +
  geom_point(mapping=NULL, myTopHits.nagR3, size = 3, shape = 16, color = "sienna2", inherit.aes = TRUE) +
  geom_point(mapping=NULL, myTopHits.cscR3 , size = 3, shape = 16, color = "slateblue", inherit.aes = TRUE) +
  geom_text_repel(aes(label = nagR), size = 3, fontface=2, color="black", max.overlaps = 100) +
  geom_text_repel(aes(label = cscR), size = 3, fontface=2, color="black") +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=0.6) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=0.6) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=0.6) +
  annotate("text", x=-6, y=-log10(0.01)+0.3,
           label=paste("Padj<0.01"), size=5, fontface="bold") +
  scale_x_continuous(limits=c(-5.5,5), breaks = -6:5) +
  labs(title="Volcano plot",
       subtitle = "B. infantis JCM1222 WT: LNnT vs Lac") +
  theme(plot.title = element_text(face="bold")) +
  theme_cowplot()

ggplotly(vplot2)

# save the figure in the eps format
#ggsave("DE_Mut_WT_Lac.eps", width = 5, height = 5)

# Create a heatmap----

colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
diffGenes.df.nagR <- subset(diffGenes.df, geneID %in% targets.nagR | geneID %in% targets.cscR) %>%
  relocate(geneID, WT_Lac1, WT_Lac2, WT_Lac3, WT_LNnT1, WT_LNnT2, WT_LNnT3, Mut_Lac1, Mut_Lac2, Mut_Lac3, Mut_LNnT1, Mut_LNnT2, Mut_LNnT3)
diffGenes.nagR <- as.matrix(diffGenes.df.nagR[,-1])
rownames(diffGenes.nagR) <- diffGenes.df.nagR$geneID

myheatcolors <- colorRampPalette(c('slateblue', 'white', 'orange'))(100)

heatmap.2(diffGenes.nagR, 
          Rowv=FALSE,
          Colv=FALSE,
          col=myheatcolors, scale='row',
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

myTopHits.df.de # first
myTopHits.df3.de  # second
diffGenes.df.two_conditions <- subset(diffGenes.df, geneID %in% myTopHits.df.de$geneID | geneID %in% myTopHits.df3.de$geneID) %>%
    relocate(geneID, WT_Lac1, WT_Lac2, WT_Lac3, Mut_Lac1, Mut_Lac2, Mut_Lac3, WT_LNnT1, WT_LNnT2, WT_LNnT3, Mut_LNnT1, Mut_LNnT2, Mut_LNnT3)
diffGenes.two_conditions <- as.matrix(diffGenes.df.two_conditions[,-1])
rownames(diffGenes.two_conditions) <- diffGenes.df.two_conditions$geneID

breaklist_nagR <- seq(0, 20, by = 1)
colorz_nagR <- colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(length(breaklist_nagR))

pheatmap(diffGenes.two_conditions,
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(3, 6, 9),
         cellwidth = 10,
         cellheight = 10)

Ploto <- pheatmap(diffGenes.two_conditions,
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         #gaps_row = c(4, 10, 17, 19),
         gaps_col = c(3, 6, 9),
         cellwidth = 20,
         cellheight = 20,
         border_color = "white",
         na_col = "#000000"
)


ggsave("heatmap.eps", width = 5, height = 8)

# Annotate DE results with mcSEED annotations----
seed.ann <- read_tsv('SEED_annotations.tsv')
corr <- read_tsv('corr.txt')
seed.ann.corr <- right_join(seed.ann, corr, by = c('seed_id' = 'seed_id'))

DE_Mut_WT_Lac <- read_tsv('DE_Mut_WT_Lac.txt')
DE_Mut_WT_Lac.out <- right_join(seed.ann.corr, DE_Mut_WT_Lac, by = c('prokka_id' = 'geneID'))
write_tsv(DE_Mut_WT_Lac.out , "DE_Mut_WT_Lac_ann.tsv")

DE_Mut_WT_LNnT <- read_tsv('DE_Mut_WT_LNnT.txt')
DE_Mut_WT_LNnT.out <- right_join(seed.ann.corr, DE_Mut_WT_LNnT, by = c('prokka_id' = 'geneID'))
write_tsv(DE_Mut_WT_LNnT.out , "DE_Mut_WT_LNnT_ann.tsv")

DE_WT_LNnT_vs_Lac<- read_tsv('DE_WT_LNnT_vs_Lac.txt')
DE_WT_LNnT_vs_Lac.out <- right_join(seed.ann.corr, DE_WT_LNnT_vs_Lac, by = c('prokka_id' = 'geneID'))
write_tsv(DE_WT_LNnT_vs_Lac.out , "DE_WT_LNnT_vs_Lac_ann.tsv")

DE_Mut_LNnT_vs_Lac<- read_tsv('DE_Mut_LNnT_vs_Lac.txt')
DE_Mut_LNnT_vs_Lac.out <- right_join(seed.ann.corr, DE_Mut_LNnT_vs_Lac, by = c('prokka_id' = 'geneID'))
write_tsv(DE_Mut_LNnT_vs_Lac.out , "DE_Mut_LNnT_vs_Lac_ann.tsv")
