library("optparse", lib.loc="/humgen/diabetes2/users/skwak/lib/")

print(paste0("Analysis Started: ", Sys.time()))

option_list=list(
	make_option(c("-i", "--input"), action="store", default=NA, type='character', help="input file storing weighted burden"),
	make_option(c("-o", "--output"), action="store", default=NA, type='character', help="output file name prefix"),
	make_option(c("-p", "--pheno"), action="store", default=NA, type='character', help="phenotype file"),

	make_option(c("-y", "--yvar"), action="store", default=NA, type='character', help="outcome variable or dependent variable"),

	make_option(c("-x", "--xvar_cov"), action="store", default=NA, type='character', help="independent categorical variable")
	)

parser <-OptionParser(usage="Rscript --vanilla %prog -i chr1.weight.txt -o output.txt -p pheno.txt -y glucose -x age,sex,bmi", option_list=option_list)
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

#Merge input file and phenotype file
data<-merge(data3, data4, by.x="SampleID", by.y=names(data4)[1], all.y=T)
#Remove rows with any NA
ndata<-data[complete.cases(data),]
mdata<-ndata
genes<-names(ndata)[2:(length(names(ndata))-length(select_col)+1)]

#Set empty vectors to store results of logistf
coef<-rep(NA, length(genes))
se<-rep(NA, length(genes))
cil<-rep(NA, length(genes))
ciu<-rep(NA, length(genes))
pval<-rep(NA, length(genes))
gene<-rep(NA, length(genes))
tcount<-rep(NA, length(genes))
wcount<-rep(NA, length(genes))
wsum<-rep(NA, length(genes))
meancomphet<-rep(NA, length(genes))
meanwild<-rep(NA, length(genes))

#Run lm for outcome variable by genotype columns adjusted for covariates
for (i in 2:(length(names(ndata))-length(select_col)+1)){
ndata[,i]<-mdata[,i]/2
gene[i-1]<-genes[i-1]
tcount[i-1]<-sum(!is.na(ndata[,i]))
wcount[i-1]<-sum(ndata[,i]!=0, na.rm=T)
wsum[i-1]<-sum(ndata[,i], na.rm=T)
meanwild[i-1]<-round(mean(ndata[ndata[,i]==0,names(ndata)==opt$yvar], na.rm=T), digits=2)
if (as.character(NA) %in% cov){
formula<-paste(opt$yvar, "~", "ndata[,i]", sep="", collapse="")
} else{
formula<-paste(opt$yvar, "~", "ndata[,i]+", paste(unlist(cov), sep="", collapse="+"), sep="", collapse="")
#In case where there is at least one CompHet/Hom mutations
if (wsum[i-1] > 0){
meancomphet[i-1]<-round(mean(ndata[ndata[,i]!=0,names(ndata)==opt$yvar], na.rm=T), digits=2)
model<-lm(formula=as.formula(formula), data=ndata)
coef[i-1]<-round(model$coefficient[2], digits=4)
se[i-1]<-round(summary(model)$coefficients[2,2], digits=4)
cil[i-1]<-round(confint(model)[2,1], digits=4)
ciu[i-1]<-round(confint(model)[2,2], digits=4)
pval[i-1]<-summary(model)$coefficients[2,4]
}
if (i == (length(names(ndata))-length(select_col)+1)) print(paste("A total of ", i, " Iteration Done!", sep=""))
}

result<-data.frame(cbind(gene, coef, se, cil, ciu, pval, tcount, wcount, wsum, meancomphet, meanwild))
names(result)<-c("GENE", "BETA", "SE.BETA", "95%CIL", "95%CIU", "P", "Total.N", "Total.Count", "Total.Weight", "Mean.CompHet", "Mean.Wild")

write.table(result, opt$output, col.names=T, row.names=F, sep="\t", quot=F)

print(paste0("Analysis Finished: ", Sys.time()))
