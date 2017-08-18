# NGS Data Integration process

Data integration is a key objective in biomedical research, as it allows the identification of hidden relationships and correlations between heterogeneous biomolecular data. This is a data integration prototype using R, currently supporting two NGS technologies (450K methylation array and RNA sequencing).

## Input Data Sources

Currently the process supports two NGS technologies; 450K methylation data (sample data produced by an Infinium Human Methylation 450k array) and RNA sequencing data (sample data produced by Illumina NextSeq 500).

### 450K methylation data

The sample data (`betas_450k_Meth.csv`) is located within the `sample_data` folder, in a `csv` format, containing the following columns:

`ID`, `Chromosome`, `Start`, `End`, `Strand` and `betas_m`

The first column (`ID`) corresponds to the identifier assigned to each methylation site, the next four columns (`Chromosome`, `Start`, `End` and `Strand`) define the exact chromosomal position of the site and the last column (`betas_m`) ...???

A snippet of the data is the following:

```
cg13869341,chr1,15865,15866,+,chr1_15865_15866
cg14008030,chr1,18827,18828,+,chr1_18827_18828
cg12045430,chr1,29407,29408,+,chr1_29407_29408
```

### RNA-Seq data

The sample data (`gene_expr_RNA_Seq.csv`) is located within the `sample_data` folder, in a `csv` format, containing the following columns:

`LOC`, `Chromosome`, `Start`, `End`, `gene.features.locus`, `Genes` and `Strand`

The first column (`LOC`) corresponds to the identifier assigned to each gene after the tuxedo protocol (???ref???), columns `Chromosome`, `Start`, `End`, `gene.features.locus` and `Strand` define the exact chromosomal position of the site and column `Genes` contains the gene names that the particular loci has been annotated with.

A snippet of the data is the following:

```
XLOC_000001,chr1,11873,29370,chr1:11873-29370,DDX11L1,+
XLOC_000002,chr1,11873,29370,chr1:11873-29370,WASH7P,+
XLOC_000003,chr1,30365,30503,chr1:30365-30503,MIR1302-10,+
```

## steps

##


_Note: Testing data can be retrieved from [genome/gms repository](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data )_
