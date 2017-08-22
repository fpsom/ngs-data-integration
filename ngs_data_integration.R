library(GenomicRanges)

# Import Data

betas <- read.csv("sample_data/betas_450k_Meth.csv", sep = ";")
expr_data <- read.csv("sample_data/gene_expr_RNA_Seq.csv", sep = ";")

# Output Data

info.within.1.2.out.file <- "Body_overlap_regions.csv"
info.tss.pos.1.2.out.file <- "TSS5-3_overlap_regions.csv"
info.tss.neg.1.2.out.file <- "TSS3-5_overlap_regions.csv"

# Data Integration: RNA-Seq (`expr_data`) and 450k Methylation data (`betas`)

## Find the overlap within the transcript

### Step A:
### find the overlapping ranges between the CpG and the LOC

print("Finding the overlap within the transcript")

gr_exp <- GRanges(seqnames=expr_data$Chromosome, IRanges(start=expr_data$Start, end=expr_data$End))

gr_betas <- GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start, end=betas$End))

print("   merging ranges by `seqnames`")

ranges <- merge(as.data.frame(gr_betas), as.data.frame(gr_exp), by="seqnames", suffixes=c("A", "B"))

ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]

### Step B:
### add a column on betas file (`betas_m`) and on expression file(`expr_m`) with the overlapping region in each case and create a total matrix (`info.within.1.2`)

info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
expr_data_<-within(expr_data, expr_m <- paste(Chromosome, Start, End,  sep='_'))

info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.within.1<- merge(info.all1, expr_data_, by.x = "expr_m", by.y = "expr_m")

colnames(info.within.1)
info.within.1$x <- paste(info.within.1$betas_m,info.within.1$expr_m, info.within.1$ID, info.within.1$LOC)
info.within.1.2<-info.within.1[!duplicated(info.within.1$x),]
colnames(info.within.1.2)

#### add a column with location
info.within.1.2$location <-  "body"

### Step C:
### save the total matrix within the transcript

print(paste0("   Saving output in file: ", info.within.1.2.out.file, sep=''))

write.table(info.within.1.2, file=info.within.1.2.out.file, col.names=T, row.names=F, quote=F, sep=",")

## Find the overlap within the TSS and the '+' strand

### Step A:
### find the TSS of the 5'-3' transcript

positive<-subset(expr_data, expr_data[,7] == "+")

print("Calculating 5'- 3' TSS ranges")

x <- as.data.frame(positive[,3]-2000)
file<- cbind(positive[,1],positive[,2],x, positive[,3], positive[,c(5:7)])
names(file)<- names(positive)

### Step B:
### find the overlapping ranges between the CpG and the LOC

gr_exp_pos<-GRanges(seqnames=positive$Chromosome,IRanges(start=positive$Start,end=positive$End))

print("   merging ranges 5'-3'")

ranges <- merge(as.data.frame(gr_betas),as.data.frame(gr_exp_pos),by="seqnames",suffixes=c("A","B"))

ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]

### Step C:
### add a column on betas file (`betas_m`) and on expression file(`expr_m`) with the overlapping region in each case and create a total matrix (`info.tss.pos.1`)

print("   calculating 5'-3' complete overlap with TSS")

info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
positive_<-within(positive, expr_m <- paste(Chromosome, Start, End,  sep='_'))

info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.tss.pos.1<- merge(info.all1, positive_, by.x = "expr_m", by.y = "expr_m")

info.tss.pos.1$x <- paste(info.tss.pos.1$betas_m,info.tss.pos.1$expr_m, info.tss.pos.1$ID, info.tss.pos.1$LOC)
info.tss.pos.1.2<-info.tss.pos.1[!duplicated(info.tss.pos.1$x),]


#### add a column with location
info.tss.pos.1$location <-  "TSS"

### Step D:
### save the total matrix within the TSS (+ strand)

print(paste0("   Saving output in file: ", info.tss.pos.1.2.out.file, sep=''))

write.table(info.tss.pos.1.2, file=info.tss.pos.1.2.out.file, col.names=T, row.names=F, quote=F, sep=",")


## Find the overlap within the TSS and the '+' strand

### Step A:
### find the TSS of the 3'-5' transcript

negative <- subset(expr_data, expr_data[,7] == "-")

x <- as.data.frame(negative[,3]-2000)
file <- cbind(negative[,1],negative[,2],x, negative[,3], negative[,c(5:7)])
names(file) <- names(negative)
print("Calculating 3'- 5' TSS ranges")

### Step B:
### find the overlapping ranges between the CpG and the LOC

gr_exp_neg <- GRanges(seqnames=negative$Chromosome,IRanges(start=negative$Start,end=negative$End))

print("   merging ranges 3'-5'")

ranges <- merge(as.data.frame(gr_betas), as.data.frame(gr_exp_neg), by="seqnames", suffixes=c("A","B"))

ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]

### Step C:
### add a column on betas file (`betas_m`) and on expression file(`expr_m`) with the overlapping region in each case and create a total matrix (`info.tss.neg.1.2`)

print("   calculating 3'-5' complete overlap with TSS")

info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
negative_<-within(negative, expr_m <- paste(Chromosome, Start, End,  sep='_'))

info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.tss.neg.1<- merge(info.all1, negative_, by.x = "expr_m", by.y = "expr_m")

####add a column with location
info.tss.neg.1$location <-  "TSS"

info.tss.neg.1$x <- paste(info.tss.neg.1$betas_m,info.tss.neg.1$expr_m, info.tss.neg.1$ID, info.tss.neg.1$LOC)
info.tss.neg.1.2<-info.tss.neg.1[!duplicated(info.tss.neg.1$x),]


### Step D:
### save the total matrix within the TSS (- strand)

print(paste0("   Saving output in file: ", info.tss.neg.1.2.out.file, sep=''))

write.table(info.tss.neg.1.2, file=info.tss.neg.1.2.out.file, col.names=T, row.names=F, quote=F, sep=",")
