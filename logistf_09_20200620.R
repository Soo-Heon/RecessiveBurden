library("optparse", lib.loc="/humgen/diabetes2/users/skwak/lib/")

print(paste0("Analysis Started: ", Sys.time()))

option_list=list(
	make_option(c("-i", "--input"), action="store", default=NA, type='character', help="input file storing weighted burden"),
	make_option(c("-o", "--output"), action="store", default=NA, type='character', help="output file name prefix"),
	make_option(c("-p", "--pheno"), action="store", default=NA, type='character', help="phenotype file"),

	make_option(c("-y", "--yvar"), action="store", default=NA, type='character', help="outcome variable or dependent variable"),

	make_option(c("-x", "--xvar_cov"), action="store", default=NA, type='character', help="independent categorical variable")
	)

parser <-OptionParser(usage="Rscript --vanilla %prog -i chr1.weight.txt -o output.txt -p pheno.txt -y t2d -x age,sex,bmi", option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
cov <- as.list(strsplit(opt$xvar_cov, ","))

library("logistf", lib.loc="/humgen/diabetes2/users/skwak/lib/")
library("svMisc", lib.loc="/humgen/diabetes2/users/skwak/lib/")

#Read in input file, phenotype file
data1<-read.table(opt$input, header=T, check.names=F)
data2<-read.table(opt$pheno, header=T)
#Remove count summary statistics
data3<-data1[7:length(data1$"GENE"),-2]
#Rename column name from "GENE" to "SampleID"
names(data3)[1] <-"SampleID"
select_col<-c(names(data2)[1],opt$yvar, unlist(cov))
data4<-data2[select_col]

###THIS PART SHOULD BE ADJUSTED###
#Merge input file and phenotype file
data<-merge(data3, data4, by.x="SampleID", by.y=names(data4)[1], all.y=T)
#Remove rows with any NA
ndata<-data[complete.cases(data),]
mdata<-ndata
genes<-names(ndata)[2:(length(names(ndata))-length(select_col)+1)]

#Set empty vectors to store results of logistf
coef<-rep(NA, length(genes))
or<-rep(NA, length(genes))
se<-rep(NA, length(genes))
cil<-rep(NA, length(genes))
ciu<-rep(NA, length(genes))
pval<-rep(NA, length(genes))
gene<-rep(NA, length(genes))
ccase<-rep(NA, length(genes))
ccont<-rep(NA, length(genes))
ncase<-rep(NA, length(genes))
ncont<-rep(NA, length(genes))
wcase<-rep(NA, length(genes))
wcont<-rep(NA, length(genes))
tcase<-rep(NA, length(genes))
tcont<-rep(NA, length(genes))

#Run Firth's bias corrected logistic regression by columns
for (i in 2:(length(names(ndata))-length(select_col)+1)){
progress(i-1,(length(names(ndata))-length(select_col)), progress.bar=T)
ndata[,i]<-mdata[,i]/2
ncase[i-1]<-table(ndata[,i], ndata[,names(ndata)==opt$yvar])[1,2]
ncont[i-1]<-table(ndata[,i], ndata[,names(ndata)==opt$yvar])[1,1]
gene[i-1]<-genes[i-1]
wcase[i-1]<-tapply(ndata[,i], ndata[,names(ndata)==opt$yvar], FUN=sum, na.rm=T)[2]
wcont[i-1]<-tapply(ndata[,i], ndata[,names(ndata)==opt$yvar], FUN=sum, na.rm=T)[1]
ccase[i-1]<-sum(ndata[ndata[,opt$yvar]==1,][,i]!=0, na.rm=T)
ccont[i-1]<-sum(ndata[ndata[,opt$yvar]==0,][,i]!=0, na.rm=T)
tcase[i-1]<-ncase[i-1] + ccase[i-1]
tcont[i-1]<-ncont[i-1] + ccont[i-1]
#In case where there is at least one CompHet/Hom mutations
formula<-paste(opt$yvar, "~", "ndata[,i]+", paste(unlist(cov), sep="", collapse="+"), sep="", collapse="")
if (sum(ndata[!is.na(ndata[,opt$yvar]),][,i],na.rm=T) > 0){
model<-logistf(formula=as.formula(formula), data=ndata)
coef[i-1]<-round(model$coefficient[2], digits=4)
or[i-1]<-round(exp(coef[i-1]), digits=4)
se[i-1]<-round(sqrt(diag(vcov(model)))[2], digits=4)
cil[i-1]<-round(model$ci.lower[2], digits=4)
ciu[i-1]<-round(model$ci.upper[2], digits=4)
pval[i-1]<-model$prob[2]
}
cat("\n")
if (i == (length(names(ndata))-length(select_col)+1)) print(paste("A total of ", i, " Iteration Done!", sep=""))
}

result<-data.frame(cbind(gene, or, coef, se, cil, ciu, pval, ccase, wcase, ccont, wcont, tcase, tcont))
names(result)<-c("GENE", "OR", "BETA", "SE.BETA", "95%CIL", "95%CIU", "P", "Case.Count", "Case.Weight", "Cont.Count", "Cont.Weight", "Total.Case", "Total.Cont")
write.table(result, opt$output, col.names=T, row.names=F, sep="\t", quot=F)

print(paste0("Analysis Finished: ", Sys.time()))
