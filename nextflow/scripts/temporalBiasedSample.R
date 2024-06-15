# script that provides a list of IDs for an LCUBE subsample of sequences
# Input: a CSV file with columns "id", "lat", "lng", and "date"; a fasta file with sequences and ids that match the CSV file
# Output: a CSV file with columns "id", "lat", "lng", and "date" for the LCUBE subsample; a fasta file with sequences and ids that match the subsample CSV file

# Load necessary libraries
library(phangorn)
library(dplyr)
library(stringr)

# Load command line arguments into variables
args <- commandArgs(trailingOnly = TRUE)
in_csv <- args[2]
id_col <- args[3]
date_col <- args[4]
in_fasta <- args[5]
out_csv <- args[6]
out_fasta <- args[7]
n <- as.integer(args[8])
seed <- as.integer(args[9])
probs <- args[10]

print(args)

# # Adding some hard-coded testing comments because as an academic, I can get away with it
# in_csv<-"../test_files/sim_0613-1.csv"
# id_col <- "id"
# date_col <- "Collection_Date"
# in_fasta <- "../test_files/sim_0613-1.fa"
# out_fasta <- "../test_files/sim_0613-1_test1.fa"
# out_csv <- "../test_files/sim_0613-1_test1.csv"
# n <- 100
# seed <- 12

# Read in the CSV and FASTA files
metadata <- read.csv(in_csv, header = TRUE, stringsAsFactors = FALSE)
fasta <- read.dna(in_fasta, as.character=TRUE, format="fasta", as.matrix=TRUE)
rownames(fasta) <- trimws(rownames(fasta))
fa <- as.phyDat(fasta)
distance <- dist.hamming(fa)

# Set the seed for reproducibility
set.seed(seed)

# Get the number of sequences in the FASTA file
N <- length(fa)
megameta<-metadata
megameta$probability<- 1
megameta[megameta$location == probs, "probability"] <- 1-megameta[megameta$location == probs, date_col]

# get the subsample
p = rep(n/N, N)
print(nrow(metadata))
print(n)
locationSampleProbabilities=eval(parse(text=probs))
s = sample(row_number(metadata), n, prob=megameta[row_number(metadata), "probability"])


# write out the metadata from the subsample and the corresponding FASTA sequences
write.csv(metadata[s,], out_csv, row.names = FALSE)
write.dna(fasta[metadata[s, id_col],], file = out_fasta, format = "fasta")
