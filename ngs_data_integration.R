#####################
# Import Data

library(GenomicRanges)

betas <- read.csv("sample_data/betas_450k_Meth.csv")
expr_data <- read.csv("sample_data/gene_expr_RNA_Seq.csv")


########################################################################
########################################################################
###########find the overlap within the transcript#######################
########################################################################

#####betas=betas1

print("betas")
gr1<-GRanges(seqnames=expr_data$Chromosome, IRanges(start=expr_data$Start, end=expr_data$End))

gr2<-GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start,end=betas$End))


ranges <- merge(as.data.frame(gr2),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
print("1st run of ranges")
ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]
print("merge the ranges")

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

info.within.1.2$location <-  "body"

##############################################

out.file <- "Body_overlap_regions.csv"
write.table(info.within.1.2,     file=out.file, col.names=T, row.names=F, quote=F, sep=",")

###########################################
#####################

########################################################################
########################################################################
###########find the overlap within the TSS + str#######################
########################################################################


positive<-subset(expr_data, expr_data[,7] == "+")

#positive
x<- as.data.frame(positive[,3]-2000)
file<- cbind(positive[,1],positive[,2],x, positive[,3], positive[,c(5:7)])
names(file)<- names(positive)
print("5-3 TSS ranges")

gr1<-GRanges(seqnames=positive$Chromosome,IRanges(start=positive$Start,end=positive$End))

gr2<-GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start,end=betas$End))


ranges <- merge(as.data.frame(gr2),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
print("5-3: 2nd run of ranges")
ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]


info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
positive_<-within(positive, expr_m <- paste(Chromosome, Start, End,  sep='_'))


info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.tss.pos.1<- merge(info.all1, positive_, by.x = "expr_m", by.y = "expr_m")

print("5-3 complete overlap with TSS")


info.tss.pos.1$x <- paste(info.tss.pos.1$betas_m,info.tss.pos.1$expr_m, info.tss.pos.1$ID, info.tss.pos.1$LOC)
info.tss.pos.1.2<-info.tss.pos.1[!duplicated(info.tss.pos.1$x),]


####add a column with location
info.tss.pos.1$location <-  "TSS"

########################################################################

out.file <- "TSS5-3_overlap_regions.csv"
write.table(info.tss.pos.1.2,     file=out.file, col.names=T, row.names=F, quote=F, sep=",")


########################################################################
########################################################################
###########find the overlap within the TSS - strd#######################
########################################################################
head(expr_data)
negative<-subset(expr_data, expr_data[,7] == "-")

#####betas=betas1
print("betas")

head(file)
#positive
x<- as.data.frame(negative[,3]-2000)
file<- cbind(negative[,1],negative[,2],x, negative[,3], negative[,c(5:7)])
names(file)<- names(negative)
print("5-3 TSS ranges")

gr1<-GRanges(seqnames=negative$Chromosome,IRanges(start=negative$Start,end=negative$End))

gr2<-GRanges(seqnames=betas$Chromosome, IRanges(start=betas$Start,end=betas$End))


ranges <- merge(as.data.frame(gr2),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
print("5-3: 2nd run of ranges")
ranges <- ranges[with(ranges, startB <= startA & endB >= endA),]


info<-within(ranges, betas_m <- paste(seqnames, startA, endA,  sep='_'))
info<- within(info, expr_m <- paste(seqnames, startB, endB,  sep='_'))
betas<-within(betas, betas_m <- paste(Chromosome, Start, End,  sep='_'))
negative_<-within(negative, expr_m <- paste(Chromosome, Start, End,  sep='_'))


info.all1<- merge(info, betas, by.x = "betas_m", by.y = "betas_m")
info.tss.neg.1<- merge(info.all1, negative_, by.x = "expr_m", by.y = "expr_m")

print("5-3 complete overlap with TSS")
####add a column with location
info.tss.neg.1$location <-  "TSS"

info.tss.neg.1$x <- paste(info.tss.neg.1$betas_m,info.tss.neg.1$expr_m, info.tss.neg.1$ID, info.tss.neg.1$LOC)
info.tss.neg.1.2<-info.tss.neg.1[!duplicated(info.tss.neg.1$x),]

########################################################################
out.file <- "TSS3-5_overlap_regions.csv"
write.table(info.tss.neg.1.2,     file=out.file, col.names=T, row.names=F, quote=F, sep=",")