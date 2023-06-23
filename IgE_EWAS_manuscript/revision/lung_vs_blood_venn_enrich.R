###########################
# miss-methyl e-forge CpGs, Lung vs. Blood
###########################

rm(list=ls())

system("hostname")
print(Sys.Date())
print(Sys.time())

libs <- c("limma", "readxl", "ggVennDiagram", "ggpubr", "gProfileR", "gplots", "QQperm", "DOSE", "enrichplot", "clusterProfiler", "org.Hs.eg.db", "R.utils", "knitr", "gprofiler2","ggpubr","gridExtra", "grid", "ggplot2", "lattice", "ggnewscale", "magrittr", "msigdbr", "data.table")

for (l in libs) {
  if (require(l, character.only = T)) {
    print(paste0(l, " loaded successfully"))
  } else {
    install.packages(l)
    require(l, character.only = T)
    print(paste0(l, " installed and loaded successfully"))
  }
}

# Save result files with timeStamp
timeStamp <- as.character(round(unclass(Sys.time())))
print(timeStamp)

sig_digits <- 2
sum_sd <- function(data, varname) {
  eval(parse(text = str_c("data[, round(summary(", varname, "), digits=2)] %>% print()")))
  eval(parse(text = str_c("print(str_c('SD: ', data[, sd(", varname, ", na.rm = T) %>% 
                                round(sig_digits)]))")))
}

# downloaded result from e-forge2 and classified into Lung vs. Blood CpGs
lb <- read_xlsx(path =  "/proj/regeps/regep00/studies/CAMP/analyses/reprk/methylation/results/IgE_paper/193_CpGs.850k.erc2-H3-all.chart_sig.xlsx", sheet = 3)

library(ggVennDiagram)
library(data.table)


x <- list(FL_H3K4me1=lb$Fetal_Lung_H3K4me1,
          LungH3K4me1=lb$Lung_H3K4me1,
          LungH3K36me3=lb$Lung_H3K36me3,
          Bl_H3K4me1=lb$Blood_H3K4me1,
          Bl_H3K36me3=lb$Blood_H3K36me3)

#x,category.names = c("1","2","3", "4")

pdf(file=file.path(plots.dir,
                   "venn_lung_vs_blood.pdf"),
    width=10, height=6)
ggVennDiagram(x, label_size=3) + scale_fill_gradient(low="lightblue",high = "pink")
dev.off()

lung <- data.frame(union(lb$Lung_H3K4me1, lb$Lung_H3K36me3))
colnames(lung) <- c("lung_cgs")
blood <- data.frame(union(lb$Blood_H3K36me3, lb$Blood_H3K4me1))
colnames(blood) <- c("blood_cgs")

cra.dir="/proj/regeps/regep00/studies/CRA"
results.dir = file.path(cra.dir, "analysis/reprk/methylation/results/IgE_paper")

ann850k <- fread(file=file.path(results.dir, "ann850k.txt"))
ann850k<- data.frame(ann850k)

ann850k$Gene <- sub(";.*", "", ann850k$UCSC_RefGene_Name);
uniq.ge.plt <- unique(ann850k$Gene)

uniq.ge.plt <- uniq.ge.plt[uniq.ge.plt!="NA"]
uniq.ge.plt <- uniq.ge.plt[uniq.ge.plt!=""]
length(uniq.ge.plt)
head(uniq.ge.plt)

lung <- merge(lung, ann850k, by.x="lung_cgs", by.y="Name") # 132 CpGs
blood <- merge(blood, ann850k, by.x="blood_cgs", by.y="Name") # 161 CpGs

# LUNG
lung_ge <- unique(lung$Gene)
lung_ge <- lung_ge[lung_ge!="NA"]
lung_ge <- lung_ge[lung_ge!=""]
length(lung_ge) # 111

lung_enrich <- gprofiler2::gost(query=lung_ge, organism = "hsapiens", custom_bg=uniq.ge.plt, multi_query = FALSE, evcodes=TRUE, correction_method="fdr", sources = c("GO:BP","GO:CC","GO:MF","KEGG","REAC","MIRNA", "TRANSFAC", "CORUM", "HP", "HPA", "WP"))
names(lung_enrich)

gem <- lung_enrich$result[,c("term_id", "term_name", "term_size", "intersection_size", "source", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "term_size", "intersection_size", "Domain", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "term_size", "intersection_size", "Domain", "FDR", "Phenotype", "Genes")]
gem <- gem[order(gem$FDR),]
head(gem); dim(gem)
write.table(gem, file=file.path(results.dir, 
                                paste0("gProfiler2_lung_cgs_", timeStamp, ".txt")),
            sep="\t", quote=F, row.names=F)

table(gem$Domain)
lung.path <- gem[gem$Domain=="KEGG" | gem$Domain=="REAC" | gem$Domain=="WP",]
dim(lung.path)

# BLOOD
blood_ge <- unique(blood$Gene)
blood_ge <- blood_ge[blood_ge!="NA"]
blood_ge <- blood_ge[blood_ge!=""]
length(blood_ge) # 131

blood_enrich <- gprofiler2::gost(query=blood_ge, organism = "hsapiens", custom_bg=uniq.ge.plt, multi_query = FALSE, evcodes=TRUE, correction_method="fdr", sources = c("GO:BP","GO:CC","GO:MF","KEGG","REAC","MIRNA", "TRANSFAC", "CORUM", "HP", "HPA", "WP"))
names(blood_enrich)

gem <- blood_enrich$result[,c("term_id", "term_name", "term_size", "intersection_size", "source", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "term_size", "intersection_size", "Domain", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "term_size", "intersection_size", "Domain", "FDR", "Phenotype", "Genes")]
gem <- gem[order(gem$FDR),]
head(gem); dim(gem)
write.table(gem, file=file.path(results.dir, 
                                paste0("gProfiler2_blood_cgs_", timeStamp, ".txt")),
            sep="\t", quote=F, row.names=F)

table(gem$Domain)
blood.path <- gem[gem$Domain=="KEGG" | gem$Domain=="REAC" | gem$Domain=="WP",]
dim(blood.path)

