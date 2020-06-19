# RecessiveBurden
Python codes to calculate recessive burden (compound heterozygous + homozygous) in annotated exome sequenced data

###Basic usage
"python VarCountSHK_v0.15_20200225.py vcf outname min_AC max_AF phenoname weightname"

###Input parameters
1) vcf     = Input vcf filename. Can be either vcf, vcf.gz.
          The input file should be annotated with Gene Symbol in 'INFO/ANN' field preferably by SnpEff.
2) outname = Prefix of output file name.
          Output files will be appended with '_counts.txt', '_weights.txt', '_variants.txt', '_sorted_genotypes.txt.gz'.
3) min_AC  = Minimum allele count for a variant to be included in the analysis.
          'INFO/AC' field will be used to filter variants based on the throshold.
          Variants with minor allele count GREATER than this threshold will be included.
          Should be an integer.
4) max_AF  = Maximum alternate allele frequency for a variant to be included in the analysis.
          'INFO/AF' field will be used to filter variants based on the threshold.
          Variants with minor allele frequency LESS than this threshold will be included.
          Should be a float.
5) phenoname = Phenotype file. cases coded as '1', controls coded as '0', missing coded as 'NA'.
          The file should be 'TAB' delimited with two columns. 
          First column being sampleIDs and second column being phenotype codes.
6) weightname = File for variant weighting based on 7 masks defined as in Flannick J et al. Nature 2019 JUNE 6 (PMID:31118516).
          The file should be 'TAB' delimited with two columns. 
          First column being variantID:AlternativeAllele (ex rs2233580:T or var10001:A) and second column being weight for each variant.
          
###Output files
1) outname_counts.txt: 
        Main output that can be used for logistic regression and down stream analysis.
        Each individual is marked as compound heterozygous ('c') or homozygous ('h') or non-recessive ('0') for each gene.
        If an individual has homozygous mutation for a specific gene, it will be marked as homozygous ('h'). 
        If the individual does not have homozygous mutation, but have compound heterozygous mutation, it will be marked as compound heterozygous ('c').
        Otherwise it will be marked as non-recessive or wild-type ('0').
        First row has the following column names:
        GENE (gene symbol), TOTAL (total count of individuals with CompHet and/or Homozygous mutation),
        COMPHET (total count of individuals with CompHet mutation), HOMO (total count of individuals with Homozygous mutation), 
        PHENO=NA (total count of individuals with CompHet and/or Homozygous mutation in missing phenotype),
        PHENO=0 (total count of individuals with CompHet and/or Homozygous mutation in controls),
        PHENO=1 (total count of individuals with CompHet and/or Homozygous mutation in cases), and
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
        GENE (gene symbol), TOTALW (total weight of individuals with CompHet and/or Homozygous mutation),
        COMPHETW (total weight of individuals with CompHet mutation), HOMOW (total weight of individuals with Homozygous mutation), 
        PHENOW=NA (total weight of individuals with CompHet and/or Homozygous mutation in missing phenotype),
        PHENOW=0 (total weight of individuals with CompHet and/or Homozygous mutation in controls),
        PHENOW=1 (total weight of individuals with CompHet and/or Homozygous mutation in cases), and
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
