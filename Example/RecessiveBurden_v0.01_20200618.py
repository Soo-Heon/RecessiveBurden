##### VarCount using PyVCF, version 0.15 April 16th 2020, by Sean Soo-Heon Kwak #####
"""
Created on Thur April 16th 15:47:26 2020

@author: Sean Soo-Heon Kwak, M.D., Ph.D., shkwak@snu.ac.kr
modified from codes by Allison Cox (allison.cox@yale.edu)
"""
"""
###Basic usage
usage   = "usage: python VarCountSHK_v0.15_20200225.py vcf outname min_AC max_AF phenoname weightname"

###Input parameters
1) vcf     = Input vcf filename. Can be either vcf, vcf.gz.
          The input file should be annotated with Gene Symbol in 'INFO/ANN' field.
2) outname = Prefix of output file name.
          Output files will be appended with '_counts.txt', '_weights.txt', '_variants.txt', '_sorted_genotypes.txt.gz'.
3) min_AC  = Minimum allele count for a variant to be included in the analysis.
          'INFO/AC' field will be used to filter variants based on the throshold.
          Variants with less than this threshold will be excluded.
          Should be an integer.
4) max_AF  = Maximum alternate allele frequency for a variant to be included in the analysis.
          'INFO/AF' field will be used to filter variants based on the threshold.
          Variants with higher allele frequency than this threshold will be excluded.
          Should be a float.
5) phenoname = Phenotype file. cases coded as '1', controls coded as '0', missing coded as 'NA'.
          The file should be 'TAB' delimited with two columns. 
          First column being sampleIDs and second column being phenotype codes.
6) weightname = File for variant weighting based on 7 masks defined as in Flannick J et al. Nature 2019 JUNE 6 (PMID:31118516).
          The file should be 'TAB' delimited with two columns. 
          First column being variantID and second column being weight for each variant.
          
###Output files
1) outname_counts.txt: 
        Main output that can be used for logistic regression and down stream analysis.
        Each individual is marked as compound heterozygous ('c') or homozygous ('h') or non-recessive ('0') for each gene.
        If an individual has homozygous mutation for a specific gene, it will be marked as homozygous ('h'). 
        If the individual does not have homozygous mutation, but have compound heterozygous mutation, it will be marked as compound heterozygous ('c').
        Otherwise it will be marked as non-recessive or normal ('0').
        First row has column names.
        Column names are GENE (gene symbol), TOTAL (total count of individuals with CompHet and/or Homozygous mutation),
        COMPHET (count of individuals with CompHet mutation), HOMO (count of individuals with Homozygous mutation), 
        PHENO=NA (count of individuals with CompHet and/or Homozygous mutation in missing phenotype),
        PHENO=0 (count of individuals with CompHet and/or Homozygous mutation in controls),
        PHENO=1 (count of individuals with CompHet and/or Homozygous mutation in cases), and
        sampleIDs. 
        Second row has phenotype information of the samples.   
2) outname_weights.txt: 
        Alternative main output that can be used for logistic regression and down stream analysis.
        The difference with 'outname_counts.txt' is that it uses variant weights to calculate burden of recessive mutations.
        Each individual has weight for recessive burden for each gene.
        For individuals with CompHet or Homozygous mutation, the weight is calculated by 
        adding the maximum weight from the maternally inherited alleles and the maximum weight from the paternally inherited alleles.
        Individuals who do not have CompHet or Homozygous mutation for the gene will have '0' value as their weight.
        First row has column names.
        Column names are GENE (gene symbol), TOTAL (total weight of individuals with CompHet and/or Homozygous mutation),
        COMPHET (weight of individuals with CompHet mutation), HOMO (weight of individuals with Homozygous mutation), 
        PHENO=NA (weight of individuals with CompHet and/or Homozygous mutation in missing phenotype),
        PHENO=0 (weight of individuals with CompHet and/or Homozygous mutation in controls),
        PHENO=1 (weight of individuals with CompHet and/or Homozygous mutation in cases), and
        sampleIDs.           
3) outname_variants.txt: 
        Information of variants that comprise the compound heterozygous or homozygous mutation 
        and list of individuals with that specific variant.
        (het) column: list of individuals heterozygous for the variant
        (hot) column: list of individuals homozygous for the variant
4) outname_sorted_genotypes.txt.gz: 
        Genotypes extracted from input vcf file. The file is sorted based on genename.
        Column names are CHROM POS ID REF ALT AC AF GENE IMPACT SampleIDs
"""

import vcf
import sys
import gzip
import os

vcf_filename=sys.argv[1]
output_filename=sys.argv[2]
min_AC=sys.argv[3]
max_AF=sys.argv[4]
phenoname=sys.argv[5]
weight_filename=sys.argv[6]

##Read in vcf file (either in vcf or vcf.gz)
vcf_reader=vcf.Reader(filename=vcf_filename)

##Designate output 'countfile' name and open with write
countfilename=output_filename+'_counts.txt'
countfile=open(countfilename, 'w')

##Designate output 'individual weight file' name and open with write
idvweightfilename=output_filename+'_weights.txt'
idvweightfile=open(idvweightfilename, 'w')

##Make list of sampleIDs
sampleIDs=vcf_reader.samples


##Designate intermediate file 'intermfile' name
intermfilename=output_filename+'_interm.txt'
#Open intermediate file with open
intermfile=open(intermfilename, 'w')
#Make header with sampleIDs for intermfile
#CHROM POS ID REF ALT AC AF GENE IMPACT SAMPLEIDs
intermfile.write('#CHROM'+'\t'+'POS'+'\t'+'ID:ALT'+'\t'+'REF'+'\t'+'ALT'+'\t'+'AC'+'\t'+'AF'+'\t'+'GENE'+'\t'+'IMPACT'+'\t')
for i in range(0,len(sampleIDs)-1):
    intermfile.write(sampleIDs[i]+'\t')
intermfile.write(sampleIDs[len(sampleIDs)-1] + '\n')

#Append individual genotype data to intermfile per variant
#CHROM POS ID REF ALT AC AF GENE IMPACT GENOTYPEs
for record in vcf_reader:
#Filter by AC>sys.argv[3] and AF<sys.argv[4]
    if ((record.INFO['AC'][0] > int(min_AC)) and (record.INFO['AF'][0] < float(max_AF))):
        intermfile.write(str(record.CHROM) + '\t' + str(record.POS) + '\t' + record.ID + ':' + str(record.ALT[0]) + '\t' + str(record.REF) + '\t' + str(record.ALT[0]) + '\t' + str(record.INFO['AC'][0]) + '\t' + str(record.INFO['AF'][0]) + '\t' + record.INFO['ANN'][0].split('|')[3] + '\t' + record.INFO['ANN'][0].split('|')[2] + '\t')
        for sample in record.samples[0:len(record.samples)-1]:
            intermfile.write(sample['GT']+'\t')
        intermfile.write(record.samples[len(record.samples)-1]['GT'] + '\n')
#Close intermfile
intermfile.close()

##Sort intermfile by "GENE"
#This is required for overlapping genes and variants will be searched for specific gene without interruption. 
#Read in the output_filename_interm.txt
intermfile=open(intermfilename, 'r')
#Write out as output_filename_sorted_genotypes.txt.gz
sortfilename=output_filename+'_sorted_genotypes.txt.gz'
#Read in and skip header
header=intermfile.readline()
#Read in lines
lines=intermfile.readlines()
#Write out as output_filename_sorted_genotypes.txt.gz
with gzip.open(sortfilename, 'wt') as sortfile:
#Paste header to sortfile
    sortfile.write(header)
#Sort by GENE
    for line in sorted(lines, key=lambda line: line.split('\t')[7]):
        sortfile.write(line)
sortfile.close()
intermfile.close()

##Remove intermfile
os.remove(intermfilename)


##Make header with sampleIDs for countfile
countfile.write('GENE'+'\t'+ 'TOTAL'+'\t'+'COMPHET'+'\t'+'HOMO'+'\t'+'PHENO=NA'+'\t'+'PHENO=0'+'\t'+'PHENO=1'+'\t')
for i in range(0,len(sampleIDs)-1):
    countfile.write(sampleIDs[i]+'\t')
countfile.write(sampleIDs[len(sampleIDs)-1]+'\n')
countfile.write('T2D'+'\t'+'TOTAL'+'\t'+'COMPHET'+'\t'+'HOMO'+'\t'+'MISSING'+'\t'+'CONTROL'+'\t'+'CASE'+'\t')

##Make header with sampleIDs for idvweightfile
idvweightfile.write('GENE'+'\t'+ 'TOTAL'+'\t'+'COMPHET'+'\t'+'HOMO'+'\t'+'PHENO=NA'+'\t'+'PHENO=0'+'\t'+'PHENO=1'+'\t')
for i in range(0,len(sampleIDs)-1):
    idvweightfile.write(sampleIDs[i]+'\t')
idvweightfile.write(sampleIDs[len(sampleIDs)-1]+'\n')
idvweightfile.write('T2D'+'\t'+'TOTALW'+'\t'+'COMPHETW'+'\t'+'HOMOW'+'\t'+'MISSINGW'+'\t'+'CONTROLW'+'\t'+'CASEW'+'\t')


##Read in sorted intermediate gzip file
with gzip.open(sortfilename, 'rt') as infile:
#Skip header
    line=infile.readline()

##Make empty dictionary to store 'sampleID':'pheno' pairs
#Initially, set all phenos to 'NA'
    d={}
    for i in range(0,len(sampleIDs)):
        d[sampleIDs[i]]='NA'
#Read in phenotype file
    phenos=open(phenoname, 'r')
    for line in phenos:
        (key,val)=line.rstrip('\n').split('\t')
        d[key]=val
    for i in range(0, len(sampleIDs)-1):
        countfile.write(d[sampleIDs[i]]+'\t')
        idvweightfile.write(d[sampleIDs[i]]+'\t')
    countfile.write(d[sampleIDs[len(sampleIDs)-1]]+'\n')
    idvweightfile.write(d[sampleIDs[len(sampleIDs)-1]]+'\n')
    phenos.close()

##Make empty dictionary to store 'variantID':'weight' pairs
#Initially, set all phenos to '0'
#Read in all lines of sorted intermediate gzip file
    w={}
    lines=infile.readlines()
    for vline in lines:
        words=vline.rstrip('\n').split('\t')
        w[words[2]]='0' #words[2] is variantID
#Read in variant weight file and generate dictionary
    weights=open(weight_filename, 'r')
    for line in weights:
        (key,val)=line.rstrip('\n').split('\t')
        w[key]=val	

#Re-read all lines of sorted intermediate gzip file from the beginning
    infile.seek(0)
    line=infile.readline()
    lines=infile.readlines()
#Set samplecount, samplecountnew, sampleweight related lists to 0s
    samplecount=[]
    samplecountnew=[]
    sampleweight_L=[]
    sampleweight_R=[]
    sampleweight_max=[]
    sampleweightnew_L=[]
    sampleweightnew_R=[]
    sampleweightnew_max=[]
    for i in range(0,len(sampleIDs)):
        samplecount.append('0')
        samplecountnew.append('0')
        sampleweight_L.append(['0']) #generate list of lists
        sampleweight_R.append(['0']) #generate list of lists
        sampleweight_max.append('0')
        sampleweightnew_L.append(['0']) #generate list of lists
        sampleweightnew_R.append(['0']) #generate list of lists
        sampleweightnew_max.append('0')
#Set initial genename as first variant GENE
    genename=lines[0].split('\t')[7]

##For each variant in same gene, test if each individual has compound heterozygous or homozygous mutation
    for line in lines:
#Split line into words[CHROM, POS, ID, REF, ALT, AC, AF, GENE, IMPACT, GTs...]
        words=line.rstrip('\n').split('\t')
#Set genename in line as genenamenew
        genenamenew=words[7]
        numwords=len(words)
    
#####THIS IS THE MAIN ITERATION#####
#Iterate each sample GTs to find CompHet/Homozygous individuals
        if genenamenew == genename:
            for i in range(9,numwords):
                if (words[i]=='0|1'):
                    sampleweight_R[i-9].append(w[words[2]]) #append variant weight to R list
                    if samplecount[i-9]=='left':
                        samplecount[i-9]='comp'
                    elif samplecount[i-9]=='right':
                        samplecount[i-9]='right'
                    elif samplecount[i-9]=='0':
                        samplecount[i-9]='right'
                    elif samplecount[i-9]=='hom':
                        samplecount[i-9]='hom'
                    elif samplecount[i-9]=='comp':
                        samplecount[i-9]='comp'
                elif (words[i]=='1|0'):
                    sampleweight_L[i-9].append(w[words[2]]) #append variant weight to L list
                    if samplecount[i-9]=='left':
                        samplecount[i-9]='left'
                    elif samplecount[i-9]=='right':
                        samplecount[i-9]='comp'
                    elif samplecount[i-9]=='0':
                        samplecount[i-9]='left'
                    elif samplecount[i-9]=='hom':
                        samplecount[i-9]='hom'
                    elif samplecount[i-9]=='comp':
                        samplecount[i-9]='comp'
                elif (words[i]=='0|0'):
                    samplecount[i-9]=samplecount[i-9]
                elif (words[i]=='1|1'):
                    sampleweight_R[i-9].append(w[words[2]]) #append variant weight to R list
                    sampleweight_L[i-9].append(w[words[2]]) #append variant weight to L list
                    samplecount[i-9]='hom'
#If new GENE starts, calculate total counts of either CompHet/Homozygous individuals
#and recode 'comp' as 'c', 'homo' as 'h', and others as '0'
#also find the maximal weight from L and R and add them to generate individual weight for the specific gene
        else:
            total=0 #total count
            total_h = 0 #homozygous count
            total_c = 0 #compound heterozygous count
            total_w = 0 #sum of total weights
            total_hw = 0 #sum of homozyogus weights
            total_cw = 0 #sum of compound heterozygous weights
            for i in range(0, len(sampleIDs)):
                if ((samplecount[i] == 'comp') or (samplecount[i]=='hom')):
#Count total number of CompHet/Homo individuals
                    total = total+1
#Compute individual weight by adding maximum of R and L
                    sampleweight_max[i]=round(float(max(sampleweight_L[i]))+float(max(sampleweight_R[i])),3)
                    total_w = round(total_w + float(sampleweight_max[i]) ,3)
                    if samplecount[i] == 'hom':
                        samplecount[i] = 'h'
#Count total number of Homo individuals
                        total_h = total_h+1
                        total_hw = round(total_hw + float(sampleweight_max[i]), 3) 
                    elif samplecount[i] == 'comp':
                        samplecount[i] = 'c'
#Count total number of CompHet individuals
                        total_c = total_c+1
                        total_cw =  round(total_cw + float(sampleweight_max[i]), 3)
                else:
                    samplecount[i]='0'
            total_na=0 #sum of counts in missing phenotype individuals
            total_control=0 #sum of counts in controls
            total_case=0 #sum of counts in cases
            total_naw=0 #sum of weights in missing phenotype individuals
            total_controlw=0 #sum of weights in controls
            total_casew=0 #sum of weights in cases
            for i in range(0, len(sampleIDs)):
                if d[sampleIDs[i]]=='NA':
                    if ((samplecount[i] == 'c') or (samplecount[i]=='h')):
#Count number of CompHet/Homo individuals in Missingness
                        total_na = total_na+1
                        total_naw = round(total_naw + float(sampleweight_max[i]), 3)
                elif d[sampleIDs[i]]=='0':
                    if ((samplecount[i] == 'c') or (samplecount[i]=='h')):
#Count number of CompHet/Homo individuals in Controls
                        total_control=total_control+1
                        total_controlw = round(total_controlw + float(sampleweight_max[i]), 3)
                elif d[sampleIDs[i]]=='1':
                    if ((samplecount[i] == 'c') or (samplecount[i]=='h')):
#Count number of CompHet/Homo individuals in Cases
                        total_case=total_case+1
                        total_casew = round(total_casew + float(sampleweight_max[i]), 3)
#Write out samplecount to countfile
            countfile.write(genename+'\t'+str(total)+'\t'+str(total_c)+'\t'+str(total_h)+'\t'+str(total_na)+'\t'+str(total_control)+'\t'+str(total_case)+'\t')
            for i in range(0, len(sampleIDs)-1):
                countfile.write(str(samplecount[i])+'\t')
            countfile.write(str(samplecount[len(sampleIDs)-1])+'\n')

#Write out individual weight to idvweightfile
            idvweightfile.write(genename+'\t'+str(total_w)+'\t'+str(total_cw)+'\t'+str(total_hw)+'\t'+str(total_naw)+'\t'+str(total_controlw)+'\t'+str(total_casew)+'\t')
            for i in range(0, len(sampleIDs)-1):
                idvweightfile.write(str(sampleweight_max[i])+'\t')
            idvweightfile.write(str(sampleweight_max[len(sampleIDs)-1])+'\n')

#Reset 'genename' and samplecounter
            genename=genenamenew
            total=0
            samplecount=[]
            samplecountnew=[]
            sampleweight_L=[]
            sampleweight_R=[]
            sampleweight_max=[]
            sampleweightnew_L=[]
            sampleweightnew_R=[]
            sampleweightnew_max=[]

            for i in range(0, len(sampleIDs)):
                samplecount.append('0')
                samplecountnew.append('0')
                sampleweight_L.append(['0'])
                sampleweight_R.append(['0'])
                sampleweight_max.append('0')
                sampleweightnew_L.append(['0'])
                sampleweightnew_R.append(['0'])
                sampleweightnew_max.append('0')

#Main Iteration for First Variant of GENENAMENEW
#Iterate each sample GTs to find CompHet/Homozygous individuals
            for i in range(9,numwords):
                if (words[i]=='0|1'):
                    sampleweightnew_R[i-9].append(w[words[2]])
                    if samplecountnew[i-9]=='left':
                        samplecountnew[i-9]='comp'
                    elif samplecountnew[i-9]=='right':
                        samplecountnew[i-9]='right'
                    elif samplecountnew[i-9]=='0':
                        samplecountnew[i-9]='right'
                    elif samplecountnew[i-9]=='hom':
                        samplecountnew[i-9]='hom'
                    elif samplecount[i-9]=='comp':
                        samplecount[i-9]='comp'
                elif (words[i]=='1|0'):
                    sampleweightnew_L[i-9].append(w[words[2]])
                    if samplecountnew[i-9]=='left':
                        samplecountnew[i-9]='left'
                    elif samplecountnew[i-9]=='right':
                        samplecountnew[i-9]='comp'
                    elif samplecountnew[i-9]=='0':
                        samplecountnew[i-9]='left'
                    elif samplecountnew[i-9]=='hom':
                        samplecountnew[i-9]='hom'
                    elif samplecount[i-9]=='comp':
                        samplecount[i-9]='comp'
                elif (words[i]=='0|0'):
                    samplecountnew[i-9]=samplecountnew[i-9]
                elif (words[i]=='1|1'):
                    sampleweightnew_L[i-9].append(w[words[2]])
                    sampleweightnew_R[i-9].append(w[words[2]])
                    samplecountnew[i-9]='hom'
            samplecount=samplecountnew
            sampleweight_L=sampleweightnew_L
            sampleweight_R=sampleweightnew_R
#This part works for the last gene that does not go to 'else' for genenamenew==genename    
    total = 0
    total_h = 0
    total_c = 0
    total_w = 0
    total_hw = 0
    total_cw = 0
    for i in range(0, len(sampleIDs)):
        if ((samplecount[i] == 'comp') or (samplecount[i]=='hom')):
#Count total number of CompHet/Homo individuals
            total=total+1
            sampleweight_max[i]=round(float(max(sampleweight_L[i]))+float(max(sampleweight_R[i])), 3)
            total_w = round(total_w + float(sampleweight_max[i]), 3)
            if samplecount[i] == 'hom':
                samplecount[i] = 'h'
#Count total number of Homo individuals
                total_h=total_h+1
                total_hw = round(total_hw + float(sampleweight_max[i]), 3)
            elif samplecount[i] == 'comp':
                samplecount[i] = 'c'
#Count total number of CompHet individuals
                total_c=total_c+1
                total_cw =  round(total_cw + float(sampleweight_max[i]), 3)
        else:
            samplecount[i] = '0'
    total_na=0
    total_control=0
    total_case=0
    total_naw=0
    total_controlw=0
    total_casew=0
    for i in range(0, len(sampleIDs)):
        if d[sampleIDs[i]]=='NA':
            if ((samplecount[i] == 'c') or (samplecount[i]=='h')):
#Count number of CompHet/Homo individuals in Missingness
                total_na=total_na+1
                total_naw = round(total_naw + float(sampleweight_max[i]), 3)
        elif d[sampleIDs[i]]=='0':
            if ((samplecount[i] == 'c') or (samplecount[i]=='h')):
#Count number of CompHet/Homo individuals in Controls
                total_control=total_control+1
                total_controlw = round(total_controlw + float(sampleweight_max[i]), 3)
        elif d[sampleIDs[i]]=='1':
            if ((samplecount[i] == 'c') or (samplecount[i]=='h')):
#Count number of CompHet/Homo individuals in Cases
                total_case=total_case+1
                total_casew = round(total_casew + float(sampleweight_max[i]), 3)
    countfile.write(genename+'\t'+str(total)+'\t'+str(total_c)+'\t'+str(total_h)+'\t'+str(total_na)+'\t'+str(total_control)+'\t'+str(total_case)+'\t')
    for i in range(0, len(sampleIDs)-1):
        countfile.write(str(samplecount[i])+'\t')
    countfile.write(str(samplecount[len(sampleIDs)-1])+'\n')

#Write out individual weight to idvweightfile
    idvweightfile.write(genename+'\t'+str(total_w)+'\t'+str(total_cw)+'\t'+str(total_hw)+'\t'+str(total_naw)+'\t'+str(total_controlw)+'\t'+str(total_casew)+'\t')
    for i in range(0, len(sampleIDs)-1):
        idvweightfile.write(str(sampleweight_max[i])+'\t')
    idvweightfile.write(str(sampleweight_max[len(sampleIDs)-1])+'\n')

    weights.close()

infile.close()

countfile.close()

idvweightfile.close()


#This part is to generate sampleIDs for each variant that comprise compound heterozygous or homozygous mutations
#Read in sorted intermediate gzip file
with gzip.open(sortfilename, 'rt') as infile:
#Skip header
    line=infile.readline() 
#Designate output 'variantfile' name
    variantfilename=output_filename+'_variants.txt'
    variantfile=open(variantfilename, 'w')
    variantfile.write('CHOM'+'\t'+'POS'+'\t'+'ID:ALT'+'\t'+'REF'+'\t'+'ALT'+'\t'+'GENE'+'\t'+'HET_SAMPLEIDs'+'\t'+'HOM_SAMPLEIDs'+'\n')

#Read in all lines of sorted intermediate gzip file
    lines=infile.readlines()
#For each variant, identify individuals with 'c' or 'h' in the same gene, and if their genotype is either '0|1', or '1|0', or'1|1'
#list their sampleIDs
    for line in lines:    
#Split line into words[CHROM, POS, ID, REF, ALT, AC, AF, GENE, IMPACT, GTs...]
        words=line.rstrip('\n').split('\t')
#Set genename in line as gene
        gene=words[7]
        numwords=len(words)
        
#Set variantsampleID list with variantfile columns: CHROM POS ID GENE
        variantsampleID=[words[0], words[1], words[2], words[3], words[4], words[7]]
#Open countfile
        countfile=open(countfilename, 'r')
#Skip header of countfile
#GENE TOTAL COMPHET HOM MISSING CONTROL CASE        
        countfile_header=countfile.readline()
        countfile_header=countfile.readline()
#Read in all lines of countfile
        countfile_lines=countfile.readlines()
#Set empty lists for Hets and Homs
        variantsampleID_het=[]
        variantsampleID_hom=[]
#For genes in countfile, test if the gene name is identical with the variant
#If it is identical, get sampleIDs who are marked as 'c', or 'h'
#Then, if that sampleID has '0|1' or '1|0' genotype for the variant, list it in the variantsampleID_het
#If that sampleID has '1|1' genotype for the variant, list it in the variantsampleID_hom        
        for clines in countfile_lines:
            cwords=clines.rstrip('\n').split('\t')
            if (gene == cwords[0]):
                for j in range(9, numwords):
                    if (((cwords[j-2]=='c') or (cwords[j-2]=='h')) and ((words[j]=='0|1') or (words[j]=='1|0'))):
                        variantsampleID_het.append(sampleIDs[j-9])
                    elif (((cwords[j-2]=='c') or (cwords[j-2]=='h')) and ((words[j]=='1|1'))):
                        variantsampleID_hom.append(sampleIDs[j-9])
        if ((len(variantsampleID_het) + len(variantsampleID_hom))>0):
            variantfile.write('\t'.join(variantsampleID)+'\t'+'(het)'+'\t'+','.join(variantsampleID_het)+'\t'+'(hom)'+'\t'+','.join(variantsampleID_hom))
            variantfile.write('\n')
        countfile.close()
    variantfile.close()
infile.close()
            
            
