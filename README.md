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

## find the overlap within the transcript
A. find the ranges
```library(GenomicRanges)
gr1<-GRanges(seqnames=expr_data$Chromosome,IRanges(start=expr_data$Start,end=expr_data$End))

gr2<-GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start,end=betas$End))


ranges <- merge(as.data.frame(gr2),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
print("1st run of ranges")
ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]
print("merge the ranges")```

B. find the overlaping region on betas file (`betas_m`) and on expression file(`expr_m`)

```info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
expr_data_<-within(expr_data, expr_m <- paste(Chromosome, Start, End,  sep='_'))


info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.within.1<- merge(info.all1, expr_data_, by.x = "expr_m", by.y = "expr_m")

colnames(info.within.1)
info.within.1$x <- paste(info.within.1$betas_m,info.within.1$expr_m, info.within.1$ID, info.within.1$LOC)
info.within.1.2<-info.within.1[!duplicated(info.within.1$x),]
colnames(info.within.1.2)

info.within.1.2$location <-  "body" ```
C. save the ranges within the transcript

 ```out.file <- "Body_overlap_regions.csv"
write.table(info.within.1.2,     file=out.file, col.names=T, row.names=F, quote=F, sep=",") ```

## find the overlap within the TSS + st

A. find the TSS
```positive<-subset(expr_data, expr_data[,7] == "+")

x<- as.data.frame(positive[,3]-2000)
file<- cbind(positive[,1],positive[,2],x, positive[,3], positive[,c(5:7)])
names(file)<- names(positive)
print("5-3 TSS ranges")```

B. find the ranges
```library(GenomicRanges)
gr1<-GRanges(seqnames=positive$Chromosome,IRanges(start=positive$Start,end=positive$End))

gr2<-GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start,end=betas$End))


ranges <- merge(as.data.frame(gr2),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
print("5-3: 2nd run of ranges")
ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]
```

C. find the overlaping region on betas file (`betas_m`) and on expression file(`expr_m`)

```info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
positive_<-within(positive, expr_m <- paste(Chromosome, Start, End,  sep='_'))


info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.tss.pos.1<- merge(info.all1, positive_, by.x = "expr_m", by.y = "expr_m")

print("5-3 complete overlap with TSS")


info.tss.pos.1$x <- paste(info.tss.pos.1$betas_m,info.tss.pos.1$expr_m, info.tss.pos.1$ID, info.tss.pos.1$LOC)
info.tss.pos.1.2<-info.tss.pos.1[!duplicated(info.tss.pos.1$x),]

info.tss.pos.1$location <-  "TSS" ```

D. save the ranges within the TSS (+ strd)

 ``` out.file <- "TSS5-3_overlap_regions.csv"
write.table(info.tss.pos.1.2,     file=out.file, col.names=T, row.names=F, quote=F, sep=",") ```

## find the overlap within the TSS - st

A. find the TSS
```negative<-subset(expr_data, expr_data[,7] == "-")

x<- as.data.frame(negative[,3]-2000)
file<- cbind(negative[,1],negative[,2],x, negative[,3], negative[,c(5:7)])
names(file)<- names(negative)
print("5-3 TSS ranges")```

B. find the ranges
```library(GenomicRanges)
gr1<-GRanges(seqnames=negative$Chromosome,IRanges(start=negative$Start,end=negative$End))

gr2<-GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start,end=betas$End))


ranges <- merge(as.data.frame(gr2),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
print("5-3: 2nd run of ranges")
ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]

```

C. find the overlaping region on betas file (`betas_m`) and on expression file(`expr_m`)

```info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
negative_<-within(negative, expr_m <- paste(Chromosome, Start, End,  sep='_'))


info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.tss.neg.1<- merge(info.all1, negative_, by.x = "expr_m", by.y = "expr_m")

print("5-3 complete overlap with TSS")
####add a column with location
info.tss.neg.1$location <-  "TSS"

info.tss.neg.1$x <- paste(info.tss.neg.1$betas_m,info.tss.neg.1$expr_m, info.tss.neg.1$ID, info.tss.neg.1$LOC)
info.tss.neg.1.2<-info.tss.neg.1[!duplicated(info.tss.neg.1$x),] ```

D. save the ranges within the transcript

 ``` out.file <- "TSS3-5_overlap_regions.csv"
write.table(info.tss.neg.1.2,     file=out.file, col.names=T, row.names=F, quote=F, sep=",")
 ```



_Note: Testing data can be retrieved from [genome/gms repository](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data )_
