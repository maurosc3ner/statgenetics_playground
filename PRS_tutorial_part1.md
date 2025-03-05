This tutorial is intended to perform QC using plink and R. It is an
adaptation of Choi et al., Tutorial: a guide to performing polygenic
risk score analyses. Nat Protoc (2020).
<https://doi.org/10.1038/s41596-020-0353-1>

## 0.1 Datasets

[Base
data](https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view?usp=sharing)

[Target
data](https://drive.google.com/file/d/1uhJR_3sn7RA8U5iYQbcmTp6vFdQiF4F2/view?usp=sharing):
Simulated data based on the 1000 Genomes Project European samples

# 1 Base data GWAS QC

## 1.1 MAF and genotyping

    ##    CHR     BP        SNP A1 A2      N         SE        P        OR      INFO
    ## 1    1 756604  rs3131962  A  G 388028 0.00301666 0.483171 0.9978869 0.8905579
    ## 2    1 768448 rs12562034  A  G 388028 0.00329472 0.834808 1.0006873 0.8958935
    ## 3    1 779322  rs4040617  G  A 388028 0.00303344 0.428970 0.9976036 0.8975083
    ## 4    1 801536 rs79373928  G  T 388028 0.00841324 0.808999 1.0020357 0.9089629
    ## 5    1 808631 rs11240779  G  A 388028 0.00242821 0.590265 1.0013083 0.8932125
    ## 6    1 809876 rs57181708  G  A 388028 0.00336785 0.714750 1.0012317 0.9235576
    ## 7    1 835499  rs4422948  G  A 388028 0.00237580 0.710884 0.9991198 0.9064377
    ## 8    1 838555  rs4970383  A  C 388028 0.00235773 0.150993 0.9966199 0.9077165
    ## 9    1 840753  rs4970382  C  T 388028 0.00207377 0.199967 0.9973457 0.9146026
    ## 10   1 846864   rs950122  C  G 388028 0.00255030 0.321085 0.9974727 0.9059456
    ##          MAF
    ## 1  0.3693896
    ## 2  0.3368458
    ## 3  0.3773680
    ## 4  0.4832122
    ## 5  0.4504096
    ## 6  0.4997439
    ## 7  0.4810160
    ## 8  0.3271640
    ## 9  0.4989362
    ## 10 0.3587834

## 1.2 Duplicate SNPs

    ##        CHR        BP         SNP A1 A2      N         SE        P        OR
    ## 116842   3 182254751   rs7622072  C  G 388028 0.00526681 0.126843 0.9919915
    ## 177369   5 139208711 rs113309830  C  G 388028 0.00377597 0.542249 0.9977015
    ##             INFO       MAF
    ## 116842 0.9032204 0.4727903
    ## 177369 0.9070873 0.4364399

How many duplicated SNPs are there? There are 2 duplicated snps

## 1.3 Ambiguous SNPs

    ## [1] 29876

How many non-ambiguous SNPs were there?

# 2 Target data QC (PLINK)

    plink \
        --bfile ../EUR/EUR \
        --maf 0.01 \
        --hwe 1e-6 \
        --geno 0.01 \
        --mind 0.01 \
        --write-snplist \
        --make-just-fam \
        --out EUR.QC

How many SNPs and samples were filtered? 14 samples removed due to
missing genotype data (–mind). 5353 variants removed due to missing
genotype data (–geno). 944 variants removed due to Hardy-Weinberg exact
test. 5061 variants removed due to minor allele threshold(s)

## 2.1 Inbreeding (PLINK)

    ./plink \
        --bfile ../EUR/EUR \
        --keep ../EUR/EUR.QC.fam \
        --extract ../EUR/EUR.QC.snplist \
        --indep-pairwise 200 50 0.25 \
        --out ../EUR/EUR.QC

Pruning complete. 272077 of 540534 variants removed.

    # Generate het file with F estimates per person
    ./plink \
        --bfile ../EUR/EUR \
        --extract ../EUR/EUR.QC.prune.in \
        --keep ../EUR/EUR.QC.fam \
        --het \
        --out ../EUR/EUR.QC

We will remove individuals with F coefficients that are more than 3
standard deviation (SD) units from the mean:

## 2.2 heterozygosity rate (R)

    ## [1] 0.003034763

    ## [1] 0.0138189

How many samples were excluded due to high heterozygosity rate? 2
samples were excluded (from 489 to 487)

## 2.3 Ambiguous SNPs

These were removed during the base data QC.

## 2.4 Mismatching SNPs (strand-flipping)

Identify SNPs that require recoding in the target (to ensure the coding
allele in the target data is the effective allele in the base summary
statistic)

    ## [1] "rs3131962"  "rs12562034" "rs4040617"  "rs79373928" "rs11240779"
    ## [6] "rs57181708"

## 2.5 Sex chromosomes (PLINK)

    plink \
        --bfile EUR \
        --extract EUR.QC.prune.in \
        --keep EUR.valid.sample \
        --check-sex \
        --out EUR.QC

How many samples were excluded due mismatched Sex information? 4 samples
were excluded

## 2.6 Relatedness

    plink \
        --bfile EUR \
        --extract EUR.QC.prune.in \
        --keep EUR.QC.valid \
        --rel-cutoff 0.125 \
        --out EUR.QC

How many related samples were excluded? 0 samples were excluded

## 2.7 Generate final QC’ed target data file

    plink \
        --bfile EUR \
        --make-bed \
        --keep EUR.QC.rel.id \
        --out EUR.QC \
        --extract EUR.QC.snplist \
        # --exclude EUR.mismatch \
        --a1-allele EUR.a1
