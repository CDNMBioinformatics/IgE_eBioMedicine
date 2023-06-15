#---
#title: "Cell type estimation CRA and CAMP - based on Varun's reference"
#author: "Priyadarshini Kachroo"
#date: "12/10/2022"
#---

# restart R session
#.rs.restartR()
rm(list=ls())

system("hostname")
print(Sys.Date())
print(Sys.time())

## load libraries
libs <- c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylationEPICmanifest",
          "limma", "wateRmelon", "minfi", "gplots", "ggplot2", "knitr", "R.utils", "impute", 
          "stats", "tidyverse", "data.table", "here", "e1071", "GGally", "ggrepel", "ENmix",
          "geneplotter", "EpiDISH")

for (l in libs) {
  if (require(l, character.only = T)) {
    print(paste0(l, " loaded successfully"))
  } else {
    install.packages(l)
    require(l, character.only = T)
    print(paste0(l, " installed and loaded successfully"))
  }
}
setwd("/udd/reprk/projects/TOPMed/scripts")

##############################
# Estimate cell counts CRA
##############################

cra.res.dir = "/proj/regeps/regep00/studies/CRA/analysis/reprk/methylation/results"
plots.cra.dir = file.path(cra.res.dir,"../plots")

#setwd(results.dir)
# Save result files with timeStamp
timeStamp <- as.character(round(unclass(Sys.time())))
print(timeStamp)

# Load QCed betas
load(file=file.path(cra.res.dir, "norm.betas.cra_hg19_clean_NOsexchr_probands_1620830664.RData"))
dim(norm.betas.rcp.auto.prob)

load(file=file.path(cra.res.dir, "../TruD_clocks_10k_FINALclocks/cent12CT.Rd"))
ImmunFrac.m.cra <- epidish(beta.m = norm.betas.rcp.auto.prob, ref.m = cent12CT.m, method = "RPC", maxit=500)$estF

pdf(file = file.path(plots.cra.dir, "ImmunFrac_epidish_cra_clean.pdf"), width = 10, height = 5)
boxplot(ImmunFrac.m.cra, las=2)
dev.off()

save(ImmunFrac.m.cra, file=file.path(cra.res.dir,
                        paste0("CRA_clean_EPIC_estimatecellcounts_ImmunFrac_epidish_", timeStamp,".RData")))

immune_deconvolution = t(ImmunFrac.m.cra)
immune_deconvolution = immune_deconvolution[,order(colnames(immune_deconvolution))]

save(immune_deconvolution, file=file.path(cra.res.dir,
                                     paste0("CRA_clean_EPIC_estimatecellcounts_Immundeconv_epidish_", timeStamp,".RData")))

rm(norm.betas.rcp.auto.prob); rm(ImmunFrac.m.cra)

##############################
# Estimate cell counts CAMP
##############################

camp.res.dir = "/proj/regeps/regep00/studies/CAMP/analyses/reprk/methylation/results"
plots.camp.dir = file.path(camp.res.dir,"../plots")

load(file=file.path(camp.res.dir, "norm.betas.camp_hg19_clean_NOsexchr_probands_1620830489.RData"))
dim(norm.betas.rcp.auto.prob)

ImmunFrac.m.camp <- epidish(beta.m = norm.betas.rcp.auto.prob, ref.m = cent12CT.m, method = "RPC", maxit=500)$estF

pdf(file = file.path(plots.camp.dir, "ImmunFrac_epidish_camp_clean.pdf"), width = 10, height = 5)
boxplot(ImmunFrac.m.camp, las=2)
dev.off()

save(ImmunFrac.m.camp, file=file.path(camp.res.dir,
                                     paste0("CAMP_clean_EPIC_estimatecellcounts_ImmunFrac_epidish_", timeStamp,".RData")))

immune_deconvolution = t(ImmunFrac.m.camp)
immune_deconvolution = immune_deconvolution[,order(colnames(immune_deconvolution))]

save(immune_deconvolution, file=file.path(camp.res.dir,
                                          paste0("CAMP_clean_EPIC_estimatecellcounts_Immundeconv_epidish_", timeStamp,".RData")))

rm(norm.betas.rcp.auto.prob); rm(ImmunFrac.m.camp)

print(Sys.Date())
print(Sys.time())
sessionInfo()

