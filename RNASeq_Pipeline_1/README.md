# NGS Data Integration

## Project folder

`home/bio_tmp/NGSDataIntegration/`

## Folder Structure

`|----/home/bio_tmp/NGSDataIntegration/
			                   |----- inputFiles/
							   
							   
			                   |----- outputFiles/

							   
			                   |----- simulatedData/
										|----- rawData/
												|----- 
												|----- 
												|----- 
										|----- nohup.out
										|----- nohupRNASeqNormal.out
										|----- nohupRNASeqTumor.out
										|----- nohupExomeSeqNormal.out
										|----- nohupExomeSeqTumor.out										
`


### RNA-SEQ Pipeline (CERTH version)
Commands are executed in `outputFiles/`  directory.


Step A1: Join Lanes per Forward / Reverse
-----------------------------------------
<not needed>


Step A2: Adapter Trimming + FastQC
----------------------------------
[1:16.82elapsed]	nohup time ./../../HumanNGS/Apps/trim_galore/trim_galore --path_to_cutadapt ../../HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt --paired ../inputFiles/RNASeq_Normal_R1.fastq ../inputFiles/RNASeq_Normal_R2.fastq > nohupFastQC_RNASeq_Normal.out 2>&1&

[1:18.90elapsed]	nohup time ./../../HumanNGS/Apps/trim_galore/trim_galore --path_to_cutadapt ../../HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt --paired ../inputFiles/RNASeq_Tumor_R1.fastq ../inputFiles/RNASeq_Tumor_R2.fastq > nohupFastQC_RNASeq_Tumor.out 2>&1&


Step B: Run TopHat (assembly)
-----------------------------

[]	nohup time tophat2 -p 8 -G /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o RNASeq_Normal /home/bio_tmp/HSapRefData/UCSCRefData/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome RNASeq_Normal_R1_val_1.fq RNASeq_Normal_R2_val_2.fq > nohupTopHatRNASeq_Normal.out 2>&1&


