setwd("/spin1/users/huant/IgE/methy/EWAS/results/")

require(R.utils)

args <- cmdArgs()
cat("User command-line arguments used when invoking R:\n")
str(args)

# Retrieve command line argument 'n', e.g. '-n 13' or '--n=13'
start <- cmdArg("s", 1L)
printf("Argument start=%d\n", start)

end <- cmdArg("e", 20L)
printf("Argument end=%d\n", end)

##
load("/spin1/users/huant/mQTLs/data/methylation-offgen3-mat-4170_fhsids.Rdata")
ls()

#sv
load("/spin1/users/huant/IgE/methy/sv/sv_IgE_offGen3_sample3471.RData")

#pheno
#pheno=read.csv("/spin1/users/huant/IgE/pheno/pheno_for_COPD_IgE.csv")
pheno=read.table("/spin1/users/huant/Huant_projects/data/MasterFile_Gene_OffGenIII_5626_PredictedCBCbyPLS.txt", sep="\t",header=T)

pheno_1=pheno[,c("SabreID","EO_PER_Pred")]
dim(sv)
sv=merge(sv,pheno_1,by.x="sabreid",by.y="SabreID")
dim(sv)
## sv p01

sv_n=280
SV_p_sum=array(0,c(sv_n,2))
colnames(SV_p_sum)=c("SV","IGE_pval")

for (i in 1:sv_n) {

SVi=paste("SV", i, sep="")
SV_tmp=sv[,paste(SVi)]

lmSV_ige=lm(SV_tmp~LOGIGE, data=sv)
ige_p_sv=summary(lmSV_ige)$coefficients[2,4]


SV_p_sum[i,1]=SVi
SV_p_sum[i,2]=ige_p_sv


}

SV_p_sum_1=as.data.frame(SV_p_sum)
SV_p_sum_sig=SV_p_sum_1[as.numeric(paste(SV_p_sum_1$IGE_pval))< 0.1,]
dim(SV_p_sum_sig)## 42, 2

sv_p01=SV_p_sum_sig[1,1]
for(i in 2:dim(SV_p_sum_sig)[1]){
    j=i+1
	sv_p01=paste(sv_p01,"+",SV_p_sum_sig[i,1],sep="")
	}



sv$"Smoke8"[is.na(sv$"Smoke8")]=0
sv$"PKYRS8"[is.na(sv$"PKYRS8")]=0

pheno_meth=merge(sv[,c(2,3,32,33,42:323)],pheno_meth,by.x="FRAMID",by.y="framid")
dim(pheno_meth)

cor(pheno_meth$EO_PER_Pred,pheno_meth[,"LOGIGE"])##0.21

SEX=as.factor(pheno_meth[,"Sex"])
AGE=pheno_meth[,"Age"]

SMK=as.factor(pheno_meth[,"Smoke8"])
PKYRS=pheno_meth[,"PKYRS8"]
Gen=as.factor(pheno_meth$idtype)
IGE=pheno_meth[,"LOGIGE"]



library(nlme)


fitmodel=function(y){

    ##y=pheno_meth[,355]
	index=!is.na(y) & !is.na(SEX) & !is.na(AGE) & !is.na(IGE) & !is.na(PKYRS) & !is.na(SMK)& !is.na(Gen)
	y=y[index]
	y=qnorm(rank(y)/(length(y)+1),mean=0,sd=1)
	
	col_temp=c("FRAMID",paste("SV",1:sv_n,sep=""), "Age", "Smoke8",
	"Sex",  "LOGIGE", "PKYRS8","Gen","pedno","CD4T","CD8T","NK","Bcell","Mono", "Gran", "EO_PER_Pred")
	
	data_temp=pheno_meth[index,col_temp]
	data_temp_y=cbind(data_temp,y)
	
	data_temp_y$"Smoke8"=as.factor(data_temp_y$"Smoke8")
	data_temp_y$"Sex"=as.factor(data_temp_y$"Sex")
	data_temp_y$"Gen"=as.factor(data_temp_y$"Gen")
	data_temp_y$"pedno"=as.factor(data_temp_y$"pedno")
	
	
	ige_result=rep(NA,5)

	
	tryCatch({
	     fm_ige=paste("y~LOGIGE+Age+Sex+Smoke8+PKYRS8+CD4T+CD8T+NK+Bcell+Mono+Gran+EO_PER_Pred+", sv_p01, sep="")
	     
		 lm2_ige=lme(as.formula(fm_ige),random=~1|pedno,data=data_temp_y)
		 
		 ige_result=summary(lm2_ige)$tTable[2,]
		 


	}, error = function(e) {print(paste("error",sep=""))})
	c(ige_result,sum(index))
}

s_n=start+346
s_e=end+346

result=as.data.frame(t(apply(pheno_meth[,c(s_n:s_e)],2,fitmodel)))

colnames(result)=c("IgE.Estimate","IgE.Std.Error","IgE.DF","IgE.t.value","IgE.P","Sample_size")


write.table(result,paste("/spin1/users/huant/IgE/methy/EWAS/results_adjEOS/IgE_offGen3_sample3471_adjSVp01_adjEOS_",start,"_",end,".csv",sep=""), sep=",", row.names=T,quote=F)


