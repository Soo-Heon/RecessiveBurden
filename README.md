# 1. Calculating RecessiveBurden From VCF
Python codes to calculate recessive burden (compound heterozygous + homozygous) in annotated exome sequenced data

## 1.1 Requirement
Python version 3.6 or above

## 1.2 Basic Usage
"python RecessiveBurden_v0.02_20200624.py -i input_file -o output -c minAC -f maxAF -w weight_file -p pheno_file"

## 1.3 Input Parameters
1) vcf_file [-i]: Input vcf filename. Can be either vcf, vcf.gz.

          The input file should be annotated with 
          Gene Symbol and IMPACT (moderate|high) in 'INFO/ANN' field preferably by SnpEff.
          Recessive burden is calculated in gene level by looping through variants. 
          Therefore, variants in a gene should be sorted in a consecutive manner.
2) output [-o]: Prefix of output file name.

          Output files will be appended with '_counts.txt', '_weights.txt', 
          '_variants.txt', and '_sorted_genotypes.txt.gz'.
3) minAC [-c]: Minimum allele count for a variant to be included in the analysis.

          'INFO/AC' field will be used to filter variants based on the throshold.
          Variants with minor allele count "GREATER OR EQUAL" to this threshold will be included.
          Should be an integer.
4) maxAF [-f]: Maximum alternate allele frequency for a variant to be included in the analysis.

          'INFO/AF' field will be used to filter variants based on the threshold.
          Variants with minor allele frequency "LESS" than this threshold will be included.
          Should be a float.
5) pheno_file [-p]: Phenotype file. cases coded as '1', controls coded as '0', missing coded as 'NA'.

          The file should be 'TAB' delimited with two columns. 
          First column being sampleIDs and second column being phenotype codes.
6) weight_file [-w]: File for variant weighting based on 7 masks defined as in Flannick J et al. Nature 2019 JUNE 6 (PMID:31118516).

          The file should be 'TAB' delimited with two columns. 
          First column being 'variantID:AlternativeAllele' (ex rs2233580:T or var10001:A) and 
          second column being weight for each variant.
          
## 1.4 Output Files
1) outname_counts.txt: 

        Main output that can be used for logistic regression and down stream analysis.
        Each individual is marked as compound heterozygous ('c') or homozygous ('h') or non-recessive ('0') for each gene.
        If an individual has homozygous mutation for a specific gene, it will be marked as homozygous ('h'). 
        If the individual does not have homozygous mutation, but have compound heterozygous mutation, 
        it will be marked as compound heterozygous ('c').
        Otherwise it will be marked as non-recessive or wild-type ('0').
        First row has the following column names:
        GENE (gene symbol),
        TOTAL (total count of individuals with CompHet and/or Homozygous mutation),
        COMPHET (total count of individuals with CompHet mutation), HOMO (total count of individuals with Homozygous mutation),
        PHENO=NA (total count of individuals with CompHet and/or Homozygous mutation in missing phenotype),
        PHENO=0 (total count of individuals with CompHet and/or Homozygous mutation in controls),
        PHENO=1 (total count of individuals with CompHet and/or Homozygous mutation in cases),
        sampleIDs. 
        Second row has phenotype information of the samples with cases coded as '1' and controls coded as '0'.
        Individual level count ('0' or '1') for each gene will be displayed from third row.
2) outname_weights.txt: 

        Alternative main output that can be used for logistic regression and down stream analysis.
        The difference with 'outname_counts.txt' is that it uses variant weights to calculate burden of recessive mutations.
        Each individual has weight for recessive burden for each gene.
        For individuals with compound heterozygous or homozygous mutation, the weight is calculated by 
        adding the maximum weight from the maternally inherited alleles and the maximum weight from the paternally inherited alleles.
        Individuals who do not have compound heterozygous or homozygous mutation for the gene will have '0' value as their weight.
        First row has the following column names:
        GENE (gene symbol),
        TOTALW (total weight of individuals with CompHet and/or Homozygous mutation),
        COMPHETW (total weight of individuals with CompHet mutation), HOMOW (total weight of individuals with Homozygous mutation),
        PHENOW=NA (total weight of individuals with CompHet and/or Homozygous mutation in missing phenotype),
        PHENOW=0 (total weight of individuals with CompHet and/or Homozygous mutation in controls),
        PHENOW=1 (total weight of individuals with CompHet and/or Homozygous mutation in cases),
        sampleIDs. 
        Second row has phenotype information of the samples with cases coded as '1' and controls coded as '0'.  
        Individual level weight (0.0 - 2.0) for each gene will be displayed from third row.
3) outname_variants.txt: 

        Information of variants that comprise the compound heterozygous or homozygous mutation 
        and list of individuals with that specific variant.
        (het) column: list of individuals heterozygous for the variant
        (hom) column: list of individuals homozygous for the variant
4) outname_sorted_genotypes.txt.gz: 

        Genotypes extracted from input vcf file. The file is sorted based on genename.
        Column names are CHROM POS ID:ALT REF ALT AC AF GENE IMPACT SampleIDs

# 2. Association Testing Using Weighted RecessiveBurden
- logistf_09_20200620.R:

          R code to perform Firth's penalized logistic regression analysis for binary trait using output from RecessiveBurden.
- linear_09_20200624.R:

          R code to perform linear regression analysis for continuous trait using output from RecessiveBurden.

## 2.1 Requirement
R version 3.6 or above
Packages logistf, svMisc
## 2.2 Basic Usage
"Rscript --vanilla logistf_09_20200620.R -i RecessiveBurden_weight.txt -o output_file -p pheno_file -y t2d -x age,sex,bmi"

## 2.3 Input Parameters
1) input_file [-i]:

          Input file should be transposed form of RecessiveBurden output either from 'outname_weights.txt' or 'outname_counts.txt'.
     Transpose can be done with "transpose.awk" script included in this repository.

          Command should be: "awk -f transpose.awk outname_weights.txt > RecessiveBurden_weight.txt".
          Once transposition is done, 'RecessiveBurden_weight.txt' can be used as input_file.

2) output_file [-o]:

          Output file name.
3) pheno_file [-p]:

          Phenotype file. 
          Cases coded as '1', controls coded as '0', missing coded as 'NA'.
          To prevent column name being identical with gene symbol, columns should be named in lower case.
4) outcome variable [-y]:

          Outcome variable or dependent variable name.
5) independent variable [-x]:

          Independent variable or covariate name.
          Multiple variables should be separated by "," without space.

## 2.4 Output Files
outname_file:

          Gene level Firth's penalized logistic regression result.
          First row has the following column names:
          GENE (gene symbol, CHROM (chromosome), POS (position), 
          OR (odds ratio), BETA (coefficient of logistic regression), SE.BETA (standard error of beta),
          95%CIL (lower 95% confidence interval of beta), 95%CIU (upper 95% confidence interval of beta), 
          P (significance P value), 
          Case.Count (number of cases being CompHet and/or Homozygous), 
          Case.Weight (sum of weights for cases being CompHet and/or Homozygous),
          Cont.Count (number of controls being CompHet and/or Homozygous), 
          Cont.Weight (sum of weights for controls being CompHet and/or Homozygous),
          Total.Case (total number of cases), Total.Cont (total number of controls)
