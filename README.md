# NGS Data Integration

Data integration is a key objective in biomedical research, as it allows the identification of hidden relationships and correlations between heterogeneous biomolecular data. This is a data integration prototype using R, currently supporting two NGS technologies (450K methylation array and RNA sequencing).

## Input Data Sources

Currently the process supports two technologies; 

The DNA methylation profiling was performed using the Infinium Human Methylation 450k array (Illumina) interrogating 485,577 CpG sites and data were analyzed using RnBeads (R package). RNA sequencing was performed using NextSeq 500 (Illumina) and data were analyzed using TopHat in Unix environment. The DNA methylation and expression levels were measured using b-values and Fragments Per Kilobase Million (FPKM), respectively. 

### 450K methylation data

The sample data (`betas_450k_Meth.csv`) is located within the `sample_data` folder, in a `csv` format, containing the following columns:

`ID`, `Chromosome`, `Start`, `End` and `Strand`

The first column (`ID`) corresponds to the identifier assigned to each methylation site (CpG site), the next four columns (`Chromosome`, `Start`, `End` and `Strand`) define the exact chromosomal position of the site

A snippet of the data is the following:

```
cg13869341,chr1,15865,15866,+
cg14008030,chr1,18827,18828,+
cg12045430,chr1,29407,29408,+
```

### RNA-Seq data

The sample data (`gene_expr_RNA_Seq.csv`) is located within the `sample_data` folder, in a `csv` format, containing the following columns:

`LOC`, `Chromosome`, `Start`, `End`, `gene.features.locus`, `Genes` and `Strand`

The first column (`LOC`) corresponds to the identifier assigned to each gene after the tuxedo protocol [1], columns `Chromosome`, `Start`, `End`, `gene.features.locus` and `Strand` define the exact chromosomal position of the site and column `Genes` contains the gene names that the particular loci has been annotated with.

A snippet of the data is the following:

```
XLOC_000001,chr1,11873,29370,chr1:11873-29370,DDX11L1,+
XLOC_000002,chr1,11873,29370,chr1:11873-29370,WASH7P,+
XLOC_000003,chr1,30365,30503,chr1:30365-30503,MIR1302-10,+
```

## Steps involved

The `R` script developed is based on the `GenomicRanges` package (`library(GenomicRanges)`)

### Stage 1: Find the overlap within the transcript

- **step A**. find the overlapping ranges between the CpG and the LOC
- **step B**. add a column on betas file (`betas_m`) and on expression file(`expr_m`) with the overlapping region in each case and create a total matrix (`info.within.1.2`)
- **step C**. save the total matrix within the transcript

### Stage 2: find the overlap within the TSS and the `+` strand

- **step A**. find the TSS of the 5'- 3' transcript
- **step B**. find the overlapping ranges between the CpG and the LOC
- **step C**. add a column on betas file (`betas_m`) and on expression file(`expr_m`) with the overlapping region in each case and create a total matrix (`info.tss.pos.1`)
- **step D**. save the total matrix within the TSS and the `+` strand

### Stage 3: find the overlap within the TSS the `-` strand

- **step A**. find the TSS of the 3'- 5' transcript
- **step B**. find the overlapping ranges between the CpG and the LOC
- **step C**. add a column on betas file (`betas_m`) and on expression file(`expr_m`) with the overlapping region in each case and create a total matrix (`info.tss.neg.1.2`)
- **step D**. save the total matrix within the TSS and the `-` strand

## References

[1] Cole Trapnell,	Adam Roberts,	Loyal Goff,	Geo Pertea,	Daehwan Kim,	David R Kelley, Harold Pimentel,	Steven L. Salzberg,	John L. Rinn	& Lior Pachter, "_Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks_", Nature Protocols 7, 562–578 (2012) doi:10.1038/nprot.2012.016 ¶6.

_Note: Testing data can be retrieved from [genome/gms repository](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data )_
