# NGS Data Integration

## Project folder

`home/bio_tmp/NGSDataIntegration/`

## Folder Structure

```
|----/home/bio_tmp/NGSDataIntegration/

			                   |----- inputFiles/
										|----- README
										|----- ExomeSeq_Normal_R1.fastq
							            |----- ExomeSeq_Normal_R2.fastq
										|----- ExomeSeq_Tumor_R1.fastq
										|----- ExomeSeq_Tumor_R2.fastq
										|----- RNASeq_Normal_R1.fastq
										|----- RNASeq_Normal_R2.fastq
										|----- RNASeq_Tumor_R1.fastq
										|----- RNASeq_Tumor_R2.fastq
							   
			                   |----- outputFiles/
										|----- RNASeq_Normal_R1_val_1.fq
										|----- RNASeq_Normal_R2_val_2.fq
										|----- RNASeq_Normal_R1.fastq_trimming_report.txt  
										|----- RNASeq_Normal_R2.fastq_trimming_report.txt
					                    |----- nohupFastQC_RNASeq_Normal.out
										|----- RNASeq_Tumor_R1_val_1.fq
										|----- RNASeq_Tumor_R2_val_2.fq
										|----- RNASeq_Tumor_R1.fastq_trimming_report.txt  
										|----- RNASeq_Tumor_R2.fastq_trimming_report.txt
										|----- nohupFastQC_RNASeq_Tumor.out
										

							   
			                   |----- simulatedData/
										|----- rawData/
												|----- alignment_practical.tar
												|----- gerald_C1TD1ACXX_7_CGATGT.bam
												|----- gerald_C2DBEACXX_3.bam
												|----- gerald_C1TD1ACXX_7_ATCACG.bam
												|----- gerald_C1TD1ACXX_8_ACAGTG.bam
										|----- nohup.out
										|----- nohupRNASeqNormal.out
										|----- nohupRNASeqTumor.out
										|----- nohupExomeSeqNormal.out
										|----- nohupExomeSeqTumor.out										
```


### RNA-SEQ Pipeline (CERTH version)
Commands are executed in `outputFiles/`  directory.


### Step A0: Get FastQ files from BAM files
```
[???] nohup time java -Xmx2g -jar $PICARD SamToFastq INPUT=rawData/gerald_C1TD1ACXX_7_ATCACG.bam FASTQ=ExomeSeq_Tumor_R1.fastq SECOND_END_FASTQ=ExomeSeq_Tumor_R2.fastq > nohupExomeSeqTumor.out 2>&1&
```


#### Step A1: Join Lanes per Forward / Reverse
-----------------------------------------
_not necessary in this work_


#### Step A2: Adapter Trimming + FastQC
----------------------------------
```
[1:16.82elapsed]	nohup time ./../../HumanNGS/Apps/trim_galore/trim_galore --path_to_cutadapt ../../HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt --paired ../inputFiles/RNASeq_Normal_R1.fastq ../inputFiles/RNASeq_Normal_R2.fastq > nohupFastQC_RNASeq_Normal.out 2>&1&
```

```
[1:18.90elapsed]	nohup time ./../../HumanNGS/Apps/trim_galore/trim_galore --path_to_cutadapt ../../HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt --paired ../inputFiles/RNASeq_Tumor_R1.fastq ../inputFiles/RNASeq_Tumor_R2.fastq > nohupFastQC_RNASeq_Tumor.out 2>&1&
```


#### Step B: Run TopHat (assembly)
-----------------------------

```
[00:52:09 elapsed]	
nohup time tophat2 -p 8 -G /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o RNASeq_Normal /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome RNASeq_Normal_R1_val_1.fq RNASeq_Normal_R2_val_2.fq > nohupTopHatRNASeq_Normal.out 2>&1&
```
```
[00:59:15 elapsed]      
nohup time tophat2 -p 8 -G /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o RNASeq_Tumor /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome RNASeq_Tumor_R1_val_1.fq RNASeq_Tumor_R2_val_2.fq > nohupTopHatRNASeq_Tumor.out 2>&1&
```

#### Step C: Run Cufflinks
---------------------
```
[2:50.32elapsed]	
nohup time cufflinks -p 8 -o RNASeq_Normal_Cufflinks RNASeq_Normal/accepted_hits.bam > nohupCufflinks_RNASeq_Normal.out 2>&1&
```
```
[1:30.30elapsed]
nohup time cufflinks -p 8 -o RNASeq_Tumor_Cufflinks RNASeq_Tumor/accepted_hits.bam > nohupCufflinks_RNASeq_Tumor.out 2>&1&
```
#### Step D1: Prepare for cuffmerge
-------------------------------
```
[assemblies.txt]
./RNASeq_Normal_Cufflinks/transcripts.gtf
....
....
....
```


#### Step D2: Run Cuffmerge
---------------------
```
[4:56.32elapsed]	nohup time cuffmerge -g /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -s /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa -p 8 assemblies.txt > nohupCuffmerge.out 2>&1&
```

#### Step E: Run CuffDiff
--------------------
```
[1:20:07elapsed]	nohup time cuffdiff -o diff_out_perGroup -b /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa -p 8 -L N,T -u merged_asm/merged.gtf ./RNASeq_Normal/accepted_hits.bam ./RNASeq_Tumor/accepted_hits.bam > nohupCuffdiff.out 2>&1&
```

#### Step F: Create counts file (comes after step B)
---------------------------
```
[1:23:43elapsed]	nohup time htseq-count --format=bam RNASeq_Normal/accepted_hits.bam /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > Sample01_htseq_counts.txt > nohupHTSeqCount01.out 2>&1&
```
