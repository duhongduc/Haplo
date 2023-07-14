# install.packages("rentrez")
library(rentrez)

##########  LOADING DATA FROM NCBI SEQUENCE DATABASE   ############### 

#Retrieve genome sequence data using rentrez
rm(list = ls(all.names = TRUE)) #clear up before work
myseq <- rentrez::entrez_fetch(db = "nucleotide", id = "LR757995", rettype = "fasta")
myseq

##########  FINDING GENES   ############### 

#Find position of every occurrence of "TATAAT"
unlist(gregexpr("TATAAT", myseq))