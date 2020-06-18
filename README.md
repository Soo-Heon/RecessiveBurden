# RecessiveBurden
Python codes to calculate recessive burden (compound heterozygous + homozygous) in annotated exome sequenced data

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
