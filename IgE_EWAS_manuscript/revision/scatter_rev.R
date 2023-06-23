#################################
# IgE scatter plots revision
# datasets: CRA and CAMP
# Author: PKachroo, June 09
#################################

#########################
# SETUP, package loading
#########################

# load the betas from the timeStamped TOPMed CRA and CAMP IgE EWAS

#######################################
# Load files and generate plot for CRA
#######################################

library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

cra.dir="/proj/regeps/regep00/studies/CRA"
results.dir = file.path(cra.dir, "analysis/reprk/methylation/results/IgE_paper")
plots.dir = file.path(results.dir, "plots")

# load phenotype and methylation data
load(file=file.path(results.dir, 
                    "CRA_betas_pheno_forIgE.EWAS_1630531041.RData"))

sel.cgs <- beta.ewas[c("cg25087851","cg10159529","cg02427831", "cg25087851", "cg18927901", "cg20315954"),]
sel.cgs <- t(sel.cgs)

selcgs.plot <- merge(pData.pheno.cra, sel.cgs, by.x="toe_ids", by.y="row.names", sort=F)

pdf(file=file.path(plots.dir, "selcgs_IgE_ggscatter_CRA.pdf"),
    width = 14, height = 10)
a <- ggscatter(selcgs.plot, x = "log10Ige", y = "cg25087851",
               add = "reg.line", ylab="cg25087851:GPR44",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
b <- ggscatter(selcgs.plot, x = "log10Ige", y = "cg10159529",
               add = "reg.line", ylab="cg10159529:IL5RA",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
c <- ggscatter(selcgs.plot, x = "log10Ige", y = "cg02427831",
               add = "reg.line", ylab="cg02427831:SIGLEC8",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
d <- ggscatter(selcgs.plot, x = "log10Ige", y = "cg25087851",
               add = "reg.line", ylab="cg25087851:GPR44",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
e <- ggscatter(selcgs.plot, x = "log10Ige", y = "cg18927901",
               add = "reg.line", ylab="cg18927901:PNPLA1",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")

f <- ggscatter(selcgs.plot, x = "log10Ige", y = "cg20315954",
               add = "reg.line", ylab="cg20315954:PMP22",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)"); 
grid.arrange(a,b,c,d,e,f, nrow = 2, ncol = 3)
dev.off()

# clearing mem
rm(beta.ewas)

#######################################
# Load files and generate plot for CAMP
#######################################

camp.dir="/proj/regeps/regep00/studies/CAMP"
results.dir = file.path(camp.dir, "analyses/reprk/methylation/results/IgE_paper")
plots.dir = file.path(results.dir, "plots")
load(file=file.path(results.dir, 
                    "CAMP_betas_pheno_forIgE.EWAS_1630531235.RData"))

sel.cgs <- beta.ewas[c("cg25087851","cg10159529","cg02427831", "cg25087851", "cg18927901", "cg20315954"),]
sel.cgs <- t(sel.cgs)

selcgs.plot <- merge(pData.pheno.meth.camp, sel.cgs, by.x="toe_ids", by.y="row.names", sort=F)

pdf(file=file.path(plots.dir, "selcgs_IgE_ggscatter_CAMP.pdf"),
    width = 14, height = 10)
a <- ggscatter(selcgs.plot, x = "LOG10IGE_iuml_F48", y = "cg25087851",
               add = "reg.line", ylab="cg25087851:GPR44",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
b <- ggscatter(selcgs.plot, x = "LOG10IGE_iuml_F48", y = "cg10159529",
               add = "reg.line", ylab="cg10159529:IL5RA",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
c <- ggscatter(selcgs.plot, x = "LOG10IGE_iuml_F48", y = "cg02427831",
               add = "reg.line", ylab="cg02427831:SIGLEC8",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
d <- ggscatter(selcgs.plot, x = "LOG10IGE_iuml_F48", y = "cg25087851",
               add = "reg.line", ylab="cg25087851:GPR44",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")
e <- ggscatter(selcgs.plot, x = "LOG10IGE_iuml_F48", y = "cg18927901",
               add = "reg.line", ylab="cg18927901:PNPLA1",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)")

f <- ggscatter(selcgs.plot, x = "LOG10IGE_iuml_F48", y = "cg20315954",
               add = "reg.line", ylab="cg20315954:PMP22",                        # Add regression line
               conf.int = TRUE, point=T, xlab="log10IgE (IU/mL)"); 
grid.arrange(a,b,c,d,e,f, nrow = 2, ncol = 3)
dev.off()
