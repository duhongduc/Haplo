############################################################################################    
#       S1: R codes for read FASTA file and find genes in sub sequences

#	      Duc Du, May 2023

#R version 4.2.3
#Copyright (C) 2020 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)
###########################################################################33

#set working directory and input file name below
#setwd("C:/Users/ducdh/haplo")

##########  INSTALL PACKAGES   ############### 

InstallPackages = FALSE #Change FALSE to TRUE to install packages for the first time

if (InstallPackages) {
  install.packages("devtools")
  install.packages("stringi")
  install.packages("rentrez")
}

library(stringi)
library(rentrez)
library(devtools)

devtools::install_github("BioGenies/tidysq")
library(tidysq)
library(seqinr)

##########  LOADING DATA FROM NCBI SEQUENCE DATABASE   ############### 

#Retrieve genome sequence data using rentrez
rm(list = ls(all.names = TRUE)) #clear up before work
fname <- rentrez::entrez_fetch(db = "nucleotide",
                               id = "LR757995", #You must enter the RefSeq or Genbank accession number of the dataset you want to read.
                               rettype = "fasta")

##########  FINDING GENES   ############### 

#Sequence creation
sq <- sq(fname)

#Subsetting sequences
sq[sq %has% "TATAAT"]

#Finding genes with specific motifs
gene1 <- find_motifs(sq, fname, "TATAAT")
gene1

#Export sq objects and output
write_fasta(sq, fname, "LR757995.fasta") #Write to FASTA file according to the original RefSeq or Genbank accession number
write.csv(gene1, "Gene1.csv", row.names=FALSE)