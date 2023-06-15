#.libPaths("/app/R-4.0.3@i86-rhel7.0/lib64/R/library")
#sessionInfo()
#assign(".lib.loc", "/app/R-4.0.3@i86-rhel7.0/lib64/R/library", envir = environment(.libPaths))
#.libPaths()
#system("hostname")
#sessionInfo()
#BiocManager::install("genefilter")
#Library.dynam.unload("Cario")
library(meffil)
library(base)
options(mc.cores=5)
system("hostname")
# Save result files with timeStamp
timeStamp <- as.character(round(unclass(Sys.time())))
print(timeStamp)

# relevant directories
setwd("/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K")
results.dir = "/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K"

qc.dir = "/proj/regeps/regep00/studies"
camp.dir = file.path(qc.dir,"CAMP/data/epigenetic/methylation/TopMed/data/freezes/20200117")

samplesheet.camp <- meffil.read.samplesheet(base="/proj/regeps/regep00/studies/CAMP/data/epigenetic/methylation/TopMed/data/freezes/20200117/LEVEL1",pattern="SampleSheet.csv$")

samplesheet.camp$Slide <- as.factor(samplesheet.camp$SentrixID)
samplesheet.camp$Sample_Plate <- as.factor(samplesheet.camp$BATCH)
samplesheet.camp$Array <- as.factor(samplesheet.camp$SentrixPosition)
samplesheet.camp$SentrixID <- NULL; samplesheet.camp$BATCH <- NULL; samplesheet.camp$SentrixPosition <- NULL

samplesheet.camp$Sex1 = samplesheet.camp$Gender
samplesheet.camp$Sex <- NA
samplesheet.camp$Sex=ifelse(samplesheet.camp$Sex1 == "M", yes = "M", no = samplesheet.camp$Sex)
samplesheet.camp$Sex=ifelse(samplesheet.camp$Sex1 == "F", yes = "F", no = samplesheet.camp$Sex)
table(samplesheet.camp$Sex)
samplesheet.camp$Sex <- as.factor(samplesheet.camp$Sex)
samplesheet.camp$Gender <- NULL; samplesheet.camp$Sex1 <- NULL
levels(samplesheet.camp$Sex)
batch_var<-c("Slide", "Array", "Sample_Plate", "Sex")

# sex mismatches identified using minfi by Siqin
sex.mismatch <- read.table(file=file.path(camp.dir, "LEVEL2/sex_mismatch.txt"), 
                           sep="\t", header=F,stringsAsFactors=FALSE)
sex.mismatch <- sex.mismatch$V1
length(sex.mismatch)
sex.mismatch

samplesheet.camp=samplesheet.camp[!samplesheet.camp$Sample_Name %in% sex.mismatch,]

# remove genotype concordance issues identified by John Ziniti
# report and stats in pairwise_concordance.txt and .html report in LEVEL2 (I removed 6 sex mismatches that would have been removed before)
rem <- c("TOE413480-BIS-v01_R08C01", "TOE926173-BIS-v01_R02C01", "TOE912812-BIS-v01_R01C01",
"TOE481991-BIS-v01_R05C01", "TOE512382-BIS-v01_R05C01", "TOE828804-BIS-v01_R05C01",
"TOE207810-BIS-v01_R08C01", "TOE780216-BIS-v01_R07C01", "TOE774703-BIS-v01_R06C01",
"TOE489508-BIS-v01_R05C01", "TOE125555-BIS-v01_R04C01", "TOE490542-BIS-v01_R01C01",
"TOE854851-BIS-v01_R03C01", "TOE982044-BIS-v01_R01C01", "TOE394568-BIS-v01_R03C01", 
"TOE471600-BIS-v01_R04C01", "TOE907865-BIS-v01_R06C01", "TOE951213-BIS-v01_R04C01",
"TOE912170-BIS-v01_R05C01", "TOE487557-BIS-v01_R03C01", "TOE903778-BIS-v01_R07C01",
"TOE409877-BIS-v01_R08C01", "TOE412218-BIS-v01_R08C01", "TOE228302-BIS-v01_R07C01",
"TOE758817-BIS-v01_R05C01", "TOE669181-BIS-v01_R01C01")
length(rem)
samplesheet.camp=samplesheet.camp[!samplesheet.camp$Sample_Name %in% rem,]

# check list of available cell type references
meffil.list.cell.type.references()

# Background and dye bias correction, sexprediction, cell counts estimates, I have chosen the one for whole blood
qc.objects <- meffil.qc(samplesheet.camp, cell.type.reference="blood gse35069 complete", verbose=TRUE)
save(qc.objects,file=file.path(results.dir, paste0("qc/meffil.qc.objects.CAMP.850K.hg19_beforeSamp_filtering_", timeStamp, ".RData")))

# Generate QC report
# you can pass genotypes generated using plink and parameters also to meffil.qc.summary, default parameters are:
#qc.parameters<-meffil.qc.parameters(
#  colour.code = NULL,
#  control.categories = NULL,
#  sex.outlier.sd = 3,
#  meth.unmeth.outlier.sd = 3,
#  control.means.outlier.sd = 5,
#  detectionp.samples.threshold = 0.2,
#  beadnum.samples.threshold = 0.2,
#  detectionp.cpgs.threshold = 0.2,
#  beadnum.cpgs.threshold = 0.2,
#  snp.concordance.threshold = 0.9,
#  sample.genotype.concordance.threshold = 0.9
#)

# update parameters
qc.parameters <- meffil.qc.parameters(sex.outlier.sd = 3, meth.unmeth.outlier.sd = 3, control.means.outlier.sd = 5, beadnum.samples.threshold = 0.2, beadnum.cpgs.threshold = 0.2, snp.concordance.threshold = 0.9, sample.genotype.concordance.threshold = 0.9, detectionp.samples.threshold = 0.25, detectionp.cpgs.threshold = 0.25)

# pass using: parameters = qc.parameters, generate QC report
qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
meffil.qc.report(qc.summary, output.file="/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K/qc/meffil.qc.CAMP.850Kreport.hg19_beforeSamp_filtering.html")
save(qc.summary, file=file.path(results.dir, paste0("qc/meffil.qcsummary.CAMP.850K.hg19_beforeSamp_filtering_", timeStamp, ".RData")))

# save removed samples and probes also as per meffil
cg.rem <- qc.summary$bad.cpgs$name
length(cg.rem)
cg.rem <- unique(cg.rem)
samp.rem <- qc.summary$bad.samples$sample.name
dim(qc.summary$bad.samples)

failed_samples1 <- as.data.frame(qc.summary$bad.samples)
failed_samples2 <- unique(failed_samples1$sample.name)
# unique failed samples
length(failed_samples2)

table(qc.summary$bad.samples$issue)
cat(paste(shQuote(failed_samples2, type="cmd"), collapse=", "))

write.table(cg.rem,file=file.path(results.dir, paste0("qc/camp_failed_cgs_hg19_",timeStamp,".txt")),sep="\t",row.names=F,quote=F)
write.table(failed_samples1,file=file.path(results.dir, paste0("qc/camp_failed_samples_metrics_hg19_",timeStamp,".txt")),sep="\t",row.names=F,quote=F)
save(cg.rem, failed_samples2, file=file.path(results.dir, paste0("qc/unique_failed_samps_CGs.CAMP.850K.hg19_",timeStamp,".RData",sep="")))

# Plot residuals remaining after fitting control matrix to decide on the number PCs
# to include in the normalization below.
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K/qc/plots/meffil.pc.CAMP.850K.hg19_beforefil.pdf",height=6,width=6)

# Remove all outlier samples if necessary
#qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

# Note: From previous reports, the sex outliers don't look that bad so we keep those samples
# Get the IDs for each quality metric. One sample may be classified as failed in more than 1 metric. Have just checked few overlaps. I have also saved the table of all failed samples in each metric so perhaps a table of unique failed samples would make sense, where mention every metric for which that sample failed, but not necessary though. I have picked metrics: failed bisulfite conversion, meth vs. unmeth, dye bias, specificity and hybridization
failed.xy <- failed_samples1[(failed_samples1$issue == "X-Y ratio outlier"),]
failed.xy <- unique(failed.xy$sample.name)
length(failed.xy)
cat(paste(shQuote(failed.xy, type="cmd"), collapse=", "))

failed.bis1 <- failed_samples1[(failed_samples1$issue == "Control probe (bisulfite1)"),]
failed.bis1 <- unique(failed.bis1$sample.name)
length(failed.bis1)
cat(paste(shQuote(failed.bis1, type="cmd"), collapse=", "))

unique(intersect(failed.xy, failed.bis1)) # 1 sample common

failed.bis2 <- failed_samples1[(failed_samples1$issue == "Control probe (bisulfite2)"),]
failed.bis2 <- unique(failed.bis2$sample.name)
length(failed.bis2)
cat(paste(shQuote(failed.bis2, type="cmd"), collapse=", "))

unique(intersect(failed.bis1, failed.bis2))
failed.mu <- failed_samples1[(failed_samples1$issue == "Methylated vs Unmethylated"),]
failed.mu <- unique(failed.mu$sample.name)
length(failed.mu)
cat(paste(shQuote(failed.mu, type="cmd"), collapse=", "))

unique(intersect(failed.mu, failed.bis1))
unique(intersect(failed.mu, failed.xy))

# selected samples to exclude based on QC report
index <- failed_samples1$issue %in% c("Control probe (dye.bias)", 
                                      "Methylated vs Unmethylated",
                                      "Control probe (bisulfite1)",
                                      "Control probe (bisulfite2)",
                                      "Control probe (hybe.21771417)",
                                      "Control probe (hybe.28684356)",
                                      "Control probe (hybe.39782321)")

outlier <- failed_samples1[index,]
dim(outlier); failed.ids <- unique(outlier$sample.name)

#failed.sel <- failed_samples1[(failed_samples1$issue == "Methylated vs Unmethylated" | failed_samples1$issue == "Control probe (bisulfite1)" | failed_samples1$issue == "Control probe (bisulfite2)" | failed_samples1$issue == "Control probe (dye.bias)" | failed_samples1$issue == "Control probe (hybe.21771417)" | failed_samples1$issue == "Control probe (hybe.28684356)" | failed_samples1$issue == "Control probe (hybe.39782321)" | failed_samples1$issue == "Control probe (spec2.G.34730329)" | failed_samples1$issue == "Control probe (spec1.ratio)" | failed_samples1$issue == "Control probe (spec1.ratio1)" | failed_samples1$issue == "Control probe (spec1.ratio2)"),]

length(failed.ids); failed.ids # finally samples removed

length(qc.objects)
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects)

save(qc.objects,file=file.path(results.dir, paste0("qc/meffil.qc.objects.CAMP.850K.hg19_afterSamp_filtering_", timeStamp, ".RData")))

# redo summary calculations and residual plot
qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
meffil.qc.report(qc.summary, output.file="/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K/qc/meffil.qc.CAMP.850Kreport.hg19_afterSamp_filtering.html")
save(qc.summary, file=file.path(results.dir, paste0("qc/meffil.qcsummary.CAMP.850K.hg19_afterSamp_filtering_", timeStamp, ".RData")))

# Further bad probes and samples and probes after filtering extreme outliers
# would not recommend any filtering at this stage
cg.rem <- qc.summary$bad.cpgs$name
samp.rem <- qc.summary$bad.samples$sample.name
dim(qc.summary$bad.samples)
failed_samples <- as.data.frame(qc.summary$bad.samples)
failed_samples <- unique(failed_samples$sample.name)
length(failed_samples)

table(qc.summary$bad.samples$issue)

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K/qc/plots/meffil.pc.CAMP.850K.hg19_afterfil.pdf",height=6,width=6)

#############################################
# Perform quantile normalization with 5 PCs
#############################################
pc <- 5
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)

str(norm.objects[[1]]$samplesheet)
for (i in 1:length(norm.objects)){
  norm.objects[[i]]$samplesheet$Slide<-as.factor(norm.objects[[i]]$samplesheet$Slide)
  norm.objects[[i]]$samplesheet$Sex<-as.factor(norm.objects[[i]]$samplesheet$Sex)
  norm.objects[[i]]$samplesheet$Sample_Plate<-as.factor(norm.objects[[i]]$samplesheet$Sample_Plate)
  norm.objects[[i]]$samplesheet$Array<-as.factor(norm.objects[[i]]$samplesheet$Array)
}

# seems meffil does quantile as a replacement for noob
# This dataset could be used for Horvath's epigenetic age estimations as well since all CGs should be in this file apart from B2 probes
save(norm.objects,file=file.path(results.dir, paste0("qc/meffil.norm.quant.CAMP.850K.hg19_afterSamp_filtering_5PCs_",timeStamp,".RData",sep="")))

# Generate normalized probe values, we just save these objects for timeStamping
# Final betas will be created using minfi pipeline though
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta,file=file.path(results.dir, paste0("qc/meffil.funnorm.betas.CAMP.850K.hg19_afterSamp_CG_filtering_5PCs_",timeStamp,".RData",sep="")))

norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables=batch_var,
  control.pcs=1:5,
  batch.pcs=1:5,
  batch.threshold=1e-30
)

# PCs calculated from autosomes
pcs <- meffil.methylation.pcs(norm.beta, probe.range=800000)
save(pcs,file=file.path(results.dir, paste0("qc/pcs5.norm.beta850K.hg19.CAMP_", timeStamp,".Robj",sep="")))

# Generate normalization report
norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs, parameters=norm.parameters)
meffil.normalization.report(norm.summary, output.file="/proj/regeps/regep00/studies/CAMP/analyses/reprk/meffil_850K/qc/meffil.norm.CAMP.850Kreport.hg19_afterSamp_CG_filtering_5PCs.html")

#########################################
# CpG site based QC - setting to missing
#########################################

# Load manifest: B4 hg19
festV1 <- read.csv("/proj/rerefs/reref00/Illumina/MethylationEPIC-v1-0-B4/lib/MethylationEPIC_v-1-0_B4.csv",skip=7,as.is=TRUE, sep=",", stringsAsFactors=FALSE)
rownames(festV1) <- festV1$IlmnID
length(grep("^ch.", rownames(festV1), value=TRUE))
## [1] 2932
#The number of "rs" probes
length(grep("^rs", rownames(festV1), value=TRUE))
## [1] 59

# Set ch and rs probes to NA
ch <- grep("^ch.", rownames(festV1), value=TRUE)
rs <- grep("^rs", rownames(festV1), value=TRUE)

length(grep("^ch.", rownames(norm.beta), value=TRUE)) # 2932
length(grep("^rs.", rownames(norm.beta), value=TRUE)) # 0

norm.beta[(rownames(norm.beta) %in% ch),] <- NA
norm.beta[(rownames(norm.beta) %in% rs),] <- NA

# Possible SNPs under probe to NA
sup.sel <- festV1[festV1$SNP_DISTANCE <5 & festV1$SNP_MinorAlleleFrequency>0.05,]
dim(sup.sel) # 50,478
sup.ids <- sup.sel$IlmnID
norm.beta["cg10136773",]
norm.beta[(rownames(norm.beta) %in% sup.ids),] <- NA

# check one of the CpGs of those
length(sup.ids);head(sup.ids)

norm.beta["cg10136773",] # should be all NAs now

# cross reactive probes to NA
cross_probes_file = paste(camp.dir, "/LEVEL2/cross_reactive_probes.txt",
                          sep = "")
if (!file.size(cross_probes_file) == 0){
  cross_probes = read.table(cross_probes_file, sep = "\t",
                            header = F, quote = "\"", fill = T)
  colnames(cross_probes) = c("sample")
  n_cross_probes = nrow(cross_probes)
  n_cross_probes
} else {
  n_cross_probes = 0
}
n_cross_probes
crossprobes <- cross_probes$sample
norm.beta[(rownames(norm.beta) %in% crossprobes),] <- NA

B2_probes_file = file.path(camp.dir, "LEVEL2/B2_probes.csv")
if (!file.size(B2_probes_file) == 0){
  B2_probes = read.csv(B2_probes_file)
  n_B2_probes = length(unique(B2_probes$x))
} else {
  n_B2_probes = 0
}
#number of B2 probes
n_B2_probes
## [1] 232
B2_probes <- B2_probes$x
norm.beta[(rownames(norm.beta) %in% B2_probes),] <- NA

dim(norm.beta) # kept sex chromosomes in

save(norm.beta,file=file.path(results.dir, paste0("qc/meffil.funnorm.betas.CAMP.850K.hg19_clean_",timeStamp,".RData",sep="")))

# if we want to set sex chromosomes to NA
meffil.list.featuresets()

# annotation hg19, b37
y<-meffil.get.features("epic")

featureset<-"epic"
autosomal.sites <- meffil.get.autosomal.sites(featureset)
autosomal.sites <- intersect(autosomal.sites, rownames(norm.beta))
#norm.beta <- norm.beta[autosomal.sites,]
norm.beta[!(rownames(norm.beta) %in% autosomal.sites),] <- NA

save(norm.beta,file=file.path(results.dir, paste0("qc/meffil.funnorm.betas.CAMP.850K.hg19_clean_without_sexchr_",timeStamp,".RData",sep="")))

sessionInfo()
gc()

timeStamp <- as.character(round(unclass(Sys.time())))
print(timeStamp)

