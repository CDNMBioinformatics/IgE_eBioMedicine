#---
#title: "Cell type estimation VDAART - based on Varun's reference"
#author: "Priyadarshini Kachroo"
#date: "12/13/2022"
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
# Estimate cell counts VDAART
##############################

results.dir <- "/proj/regeps/regep00/studies/VDAART/data/methylation/VDAART_850K/data/freezes/20210823/analysis/celltype_est_immdeconv"
plots.dir = file.path(results.dir, "plots")

#setwd(results.dir)
# Save result files with timeStamp
timeStamp <- as.character(round(unclass(Sys.time())))
print(timeStamp)

# Load QCed betas
load(file=file.path(results.dir, "../../850K_QC/minfi/norm.betas.rcp.VDAART_hg19_clean_NOsexchr_1639171641.RData"))
dim(norm.betas.rcp.auto)
#[1] 865859    574 # total
#786598    574 --> after removing bad probes and sex chr

qc.dir = "/proj/regeps/regep00/studies"
cra.res.dir = file.path(qc.dir,"CRA/analysis/reprk/methylation/results")

load(file=file.path(cra.res.dir, "../TruD_clocks_10k_FINALclocks/cent12CT.Rd"))
ImmunFrac.m.vda <- epidish(beta.m = norm.betas.rcp.auto, ref.m = cent12CT.m, method = "RPC", maxit=500)$estF

pdf(file = file.path(plots.dir, "ImmunFrac_epidish_vda.pdf"), width = 10, height = 5)
boxplot(ImmunFrac.m.vda, las=2)
dev.off()

save(ImmunFrac.m.vda, file=file.path(results.dir,
                        paste0("VDAART_EPIC_estimatecellcounts_ImmunFrac_epidish_", timeStamp,".RData")))

immune_deconvolution = t(ImmunFrac.m.vda)
immune_deconvolution = immune_deconvolution[,order(colnames(immune_deconvolution))]

save(immune_deconvolution, file=file.path(results.dir,
                                     paste0("VDAART_EPIC_estimatecellcounts_Immundeconv_epidish_", timeStamp,".RData")))

print(Sys.Date())
print(Sys.time())
sessionInfo()

