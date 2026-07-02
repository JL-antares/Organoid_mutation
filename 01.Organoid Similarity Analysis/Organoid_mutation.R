rm(list = ls())
library(BSgenome)
library(MutationalPatterns)
library(patchwork)
library(ggplot2)
library(ggsci)
library(scales)
library(gridExtra)
library(VennDiagram) 
library(venn)
library(ggsci)
library(ggplot2)
library(ggpolypath)
library(scales)
library(reshape2)
library(data.table)
library(ggupset)
library(ggplotify) 
library(ggimage) 
library(maftools)
library(ComplexHeatmap)
library(dplyr)
library(tidyverse)
library(circlize)
library(ggupset)
library(ggplotify) 
library(ggimage) 
library(UpSetR)
library(venn)
library(scales)
library(openxlsx)




mycolor = pal_npg("nrc")(10)

mycolor <- pal_aaas(alpha = 0.8)(10)


path <- 


if (F) {
  cnv_files <- list.files(path=path, pattern = ".cns", full.names = TRUE)
  sample_names <- c()
  names(cnv_files) <- sample_names
  
  lst <- list()
  for (i in names(cnv_files)) {
    df <- read.table(cnv_files[i], sep = "\t", header = T)
    df$Sample <- i
    lst[[i]] <- df
  }
  df <- do.call(rbind, lst)
  

  geneList <- readxl::read_xlsx("./GENE.xlsx", sheet = 1)
  

  select_gene <- c()
  
  genes <- c(geneList$Symbol[1:30], select_gene)
  
  df <- read.table("./data/", sep = "\t", header = T)
  dfT <- df[df$gene %in% genes, c(4, 6)]
  colnames(dfT) <- c("gene", "A210_T")
  dfT <- dfT[!duplicated(dfT$gene), ]
  rownames(dfT) <- dfT$gene
  
  
  df <- read.table("./data/", sep = "\t", header = T)
  dfP <- df[df$gene %in% genes, c(4, 6)]
  colnames(dfP) <- c("gene", "A210_P")
  dfP <- dfP[!duplicated(dfP$gene), ]
  rownames(dfP) <- dfP$gene
  
  
  df <- left_join(dfT, dfP, by = "gene") %>% column_to_rownames("gene")
  
  mycol <-  colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
  
  
  col <- colorRampPalette(c("#3D2E91", "white", "#EA3131"))(30)
  
  png("3.cnv_stat.png",res = 300)
  Heatmap(as.matrix(df),
          # row_km = 2,
          # column_km = 3,
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          row_names_side = "left",
          show_row_dend = F,
          name = "Copy Number",
          rect_gp = gpar(col = "white", lwd = 2),
          # row_names_gp = gpar(fontsize = 12,
          #                     # fontface='bold',
          #                     col = c(rep("red", 48), rep("#7D51FD", 2))),
          row_title = "Gene",
          row_dend_side = "right",
          
          column_dend_side = "top", 
          col = col
          # col = rev(rainbow(10))
  )
  
  dev.off()
  
}


file <- "./data/allsamples.vep.maf"
maf <- read.maf(file)

if (F) {
  samples <- c()
  
  data <- maf@data[maf@data$Tumor_Sample_Barcode %in% samples,]
  maf <- read.maf(data)
}

mycolor <- pal_aaas(alpha = 0.7)(8)
show_col(mycolor)
vc_cols <- mycolor[1:7]
vc_cols <- c("darkgreen", "#875CFC", "#FE6E4A", "#3D2E91", "purple", "orange", "#8F88B9", "#BB0021B2")
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)








geneList <- c()


geneList <- intersect(maf@data$Hugo_Symbol, geneList)
sampleOrder <- sort(maf@variants.per.sample$Tumor_Sample_Barcode)
last_two  <- substring(sampleOrder,
                         nchar(as.character(sampleOrder))-1,
                         nchar(as.character(sampleOrder)))  
sorted_indices <- order(last_two)  
sampleOrder <- sampleOrder[sorted_indices]  

if (T) {
  pdf("RESULT/B.cancer_genes_heatmap.pdf", height = 5, width = 4.5)
  oncoplot(maf = maf,
           

           top = 30, # 展示top基因
           
           colors = vc_cols,
          
           # showPct = TRUE,
           drawRowBar = F,
           drawColBar = F,
           rightBarData = F,
           rightBarLims = F,
           numericAnnoCol = F,
           includeColBarCN = F,
           topBarData = F,
           # drawBox = TRUE,
           # annoBorderCol = "white",
           # removeNonMutated = TRUE,
           #genesToIgnore = genes,
           
           # altered = TRUE,
           # colbar_pathway = TRUE,
          
           sampleOrder = sampleOrder,


           keepGeneOrder=TRUE,

           titleFontSize = 0.8,  #图题头的字体大小
           fontSize = 0.4, # 基因名的字体大小
          
           legendFontSize = 0.6,
           
           legend_height = 2.0,
           
           showTumorSampleBarcodes = T, # 是否展示样本名/样本ID
           SampleNamefontSize = 0.6, # 样本名字体大小
           sepwd_samples = 2,# 样本间距
           showPct = F
  )
  dev.off()
}

sampleInfo <- read.xlsx("../Groups.xlsx", sheet = 2)
maf@data$Group <- maf@data %>% 
  pull(Tumor_Sample_Barcode) %>% 
  plyr::mapvalues(., from = sampleInfo$Samples, to = sampleInfo$Group2)

vaf_df <- maf@data %>% mutate(t_vaf=t_alt_count/t_ref_count) 
df <- vaf_df %>% transmute(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification, t_vaf, Group)
# df <- df[df$Hugo_Symbol %in% geneList,]
wide_df <- dcast(df, value.var = "t_vaf")

ggplot(data = df, aes(x = Hugo_Symbol, y = t_vaf)) +
  geom_point(aes(color = Variant_Classification)) +
  labs(title = "Scatter Plot of Variant Allele Frequency",
       x = "Hugo Symbol",
       y = "t_vaf") +
  theme_minimal()

# unique(grep("Tissue", df$Tumor_Sample_Barcode, value = T, invert = F))

path <- "data"
files <- list.files(path, pattern =)

mut.list <- list()
mut.df <- list()
for (file in files){ 
  fileName <- sub(, "", file, perl = T)
  myfile <- paste0(path, "/", file)
  df <- readxl::read_excel(myfile, sheet = 1)
  df$mut <- paste(df$Chr, df$Start, df$End, df$Ref, df$Alt, sep = "_" )
  df$Samples <- fileName
  mut.list[[fileName]] <- df$mut
  mut.df[[fileName]] <- df
  # mut.list[[fileName]] <- unique(df$Gene)
}

mut.df <- do.call(rbind, mut.df)

p1 <- upset(fromList(mut.list), nsets = length(mut.list),
            order.by = "freq", 
            number.angles = 0, 
            point.size = 3, 
            line.size = 1, 
            mainbar.y.label = "Count of Intersection",  # y 标题
            sets.x.label = "Datasets Size",  # x 标题
            text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
            query.legend = "top",
          
)
p1
g1 <- as.ggplot(p1) 

ggsave("RESULT/D.samples.upset.mutSites.png", plot = g1, width = 12, height = 14)
ggsave("RESULT/D.samples.upset.mutSites.pdf", plot = g1, width = 12, height = 14)



p2 <- venn(x = mut.list, 

           zcolor=mycolor,

           opacity = 0.5,
 
           box=F, 
  
           sncs=0.8, 

           ilcs=0.8, ggplot = T
)
p2
g2 <- as.ggplot(p2)



p <- g1 + geom_subview(subview = g2 , x=.75, y=.75, w=.55, h=.55)
p

ggsave("RESULT/D.somatic_intersect.upset.venn.mutSites.png", plot = p, width = 12, height = 8)
ggsave("RESULT/D.somatic_intersect.upset.venn.mutSites.pdf", plot = p, width = 12, height = 8)

df_inter <- get.venn.partitions(mut.list)
for (i in 1:nrow(df_inter)) {
  df_inter[i,'values'] <- paste(df_inter[[i,'..values..']], collapse = ',')
}
df_inter <- df_inter[order(df_inter$..count.., decreasing = T),]
nsample <- length(mut.list)

diff <- mut.df
for (i in 1:nrow(df_inter)) {
  fileName <- paste(colnames(df_inter)[1:nsample][t(df_inter[i, 1:nsample])], collapse = "_")
  df_inter[i, 'fileName'] <- paste(fileName, df_inter[i, "..count.."], sep = ".") 
  filter_id <- unlist(strsplit(df_inter[i, "values"], split = ",", fixed = T))
  df2 <- diff[diff$mut %in% filter_id, ] %>% distinct()
  openxlsx::write.xlsx(df2, paste(, df_inter[i, "fileName"],  sep = ""), overwrite = T)
}



head(available.genomes())

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

library(ref_genome, character.only = TRUE)

vcf_files <- list.files(path=, full.names = F)
names(vcf_files) <- sub(, "", vcf_files, perl = T)
sample_names <- names(vcf_files)
vcfs <- read_vcfs_as_granges(paste(, vcf_files, sep = "/"), sample_names, ref_genome)


type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
type_occurrences$samples <- rownames(type_occurrences)
writexl::write_xlsx(type_occurrences,)

group <- sample_names
p <- plot_spectrum(type_occurrences, by = group, CT = TRUE)
ggsave("C.Mutational_Pattern.samples.pdf", p, width = 12, height = 8)
ggsave("C.Mutational_Pattern.samples.png", p, width = 12, height = 8)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
df <- as.data.frame(mut_mat)
df$Type <- rownames(df)
writexl::write_xlsx(df, "Mutational_spectrum.xlsx")

p <- plot_96_profile(mut_matrix = mut_mat, condensed = TRUE)
ggsave("Mutational_spectrum.pdf", p, width = 18, height = 24)
ggsave("Mutational_spectrum.png", p, width = 18, height = 24)
plot_compare_profiles(mut_mat[,1], mut_mat[,2], condensed = TRUE)
library(NMF)
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank=2:8, method="brunet", nrun=10, seed=123456)

plot(estimate)

nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 10)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")

df <- as.data.frame(nmf_res$signatures)
df$Type <- rownames(df)

writexl::write_xlsx(df, "Mutational_signatures.xlsx")

p <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
p <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative")

ggsave("Mutational_signatures.pdf", p, width = 12, height = 8)
ggsave("Mutational_signatures.png", p, width = 12, height = 8)
p <- plot_compare_profiles(nmf_res$signatures[,1], nmf_res$signatures[,2], condensed = TRUE)
ggsave("Mutational_signatures.compare.pdf", p, width = 12, height = 8)
ggsave("Mutational_signatures.compare.png", p, width = 12, height = 8)

pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "absolute")

p <- grid.arrange(pc1, pc2)

ggsave("Mutational_signatures.contribution.pdf", p, width = 12, height = 8)
ggsave("Mutational_signatures.contribution.png", p, width = 12, height = 8)
p <- plot_contribution_heatmap(nmf_res$contribution,sig_order = c("Signature B", "Signature A"))
ggsave("Mutational_signatures.contribution.pheatmap.pdf", p, width = 12, height = 8)
ggsave("Mutational_signatures.contribution.pheatmap.png", p, width = 12, height = 8)
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)

cancer_signatures = cancer_signatures[as.vector(new_order),]

row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type

cancer_signatures = as.matrix(cancer_signatures[,4:33])
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")

cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)

plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)

fit_res <-fit_to_signatures(mut_mat, cancer_signatures)

select <- which(rowSums(fit_res$contribution) > 10)
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")
plot_contribution_heatmap(fit_res$contribution,cluster_samples = TRUE,method = "complete")

library(maftools)
library(ggplot2)
file <- ""
maf <- read.maf(file)

laml.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)# plot titv summary

plotTiTv(res = laml.titv, plotNotch = T)

titv <- laml.titv$TiTv.fractions
titv$Group <- ifelse(grepl("^LCO", titv$Tumor_Sample_Barcode), "LCO", "T_LCO")

df <- reshape2::melt(titv, value.name = "Value")

ggplot(df, aes(x = variable, y = Value)) +
  geom_boxplot(fill = "grey",
               color = "black",
               width = 0.5,
               notch = FALSE,
               size = 0.5
  ) +
  labs(title = "", x = "", y = "Proportion(%)") +
  theme_classic() +
  facet_wrap(~Group) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5, size = 12),  
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12),  
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12)  
  )
ggsave("plot_Ti-Tv.G1.pdf", width = 6, height = 6)
