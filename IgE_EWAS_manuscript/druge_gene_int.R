# Code for drug-gene interactions

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("rDGIdb")

library(rDGIdb)
camp.dir="/proj/regeps/regep00/studies/CAMP"
results.camp.dir = file.path(camp.dir, "analyses/reprk/methylation/results/IgE_paper")
plots.camp.dir = file.path(results.camp.dir, "plots")
ges <- read.table(file.path(results.camp.dir, "CAMP_CRA_FHS_DMPs_fdr_logIgE_hg19_overlap_1630601015.txt"),sep="\t", header=T,stringsAsFactors=FALSE, quote="")
ges <- unique(ges$Gene.x)
length(ges)
ges <- ges[ges!="NA"]
ges <- ges[ges!=""]
result <- queryDGIdb(ges)

timeStamp <- as.character(round(unclass(Sys.time())))
print(timeStamp)

# all results
save(result, file=file.path(results.camp.dir,
                            paste0("rDGIdb_int_IgE_ge_overlap_", timeStamp, ".RData")))

## Result summary
result.summary <- resultSummary(result)
## Detailed results
detailed.results <- detailedResults(result)
## By gene
by.gene <- byGene(result)
## Search term summary
termsum <- searchTermSummary(result)

write.table(result.summary,file=file.path(results.camp.dir, 
                                      paste0("result.summary_rdgi_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

write.table(detailed.results,file=file.path(results.camp.dir, 
                                          paste0("detailed.results_rdgi_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

write.table(by.gene,file=file.path(results.camp.dir, 
                                          paste0("by.gene_rdgi_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

write.table(termsum,file=file.path(results.camp.dir, 
                                   paste0("termsum_rdgi_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

pdf(file = file.path(plots.camp.dir, "barplot_interactions_by_source_rdgi.pdf"),width = 6, height = 5)
plotInteractionsBySource(result, main = "Number of interactions by source")
dev.off()

# by specific attributes
resultFilter <- queryDGIdb(ges,
                           sourceDatabases = c("DrugBank","ChemblInteractions"),
                           geneCategories = "CLINICALLY ACTIONABLE",
                           interactionTypes = c("suppressor","activator","blocker"))

# no interactions to plot
#plotInteractionsBySource(resultFilter, main = "Number of interactions by source")

# by specific attributes
resultFilter <- queryDGIdb(ges,
                           geneCategories = "CLINICALLY ACTIONABLE",
                           interactionTypes = c("suppressor","activator","blocker"))

# no interactions to plot
#plotInteractionsBySource(resultFilter, main = "Number of interactions by source")

resultFilter <- queryDGIdb(ges,
                           geneCategories = "CLINICALLY ACTIONABLE")

## Result summary
result.summary <- resultSummary(resultFilter)
## Detailed results
detailed.results <- detailedResults(resultFilter)
## By gene
by.gene <- byGene(resultFilter)
## Search term summary
termsum <- searchTermSummary(resultFilter)

write.table(result.summary,file=file.path(results.camp.dir, 
                                          paste0("result.summary_rdgi_ca_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

write.table(detailed.results,file=file.path(results.camp.dir, 
                                            paste0("detailed.results_rdgi_ca_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

write.table(by.gene,file=file.path(results.camp.dir, 
                                   paste0("by.gene_rdgi_ca_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

write.table(termsum,file=file.path(results.camp.dir, 
                                   paste0("termsum_rdgi_ca_",timeStamp,".txt"))
            ,sep="\t",row.names=F,quote=F)

