# script that provides a list of IDs for an LCUBE subsample of sequences
# Input: a CSV file with columns "id", "lat", "lng", and "date"; a fasta file with sequences and ids that match the CSV file
# Output: a CSV file with columns "id", "lat", "lng", and "date" for the LCUBE subsample; a fasta file with sequences and ids that match the subsample CSV file

# Load necessary libraries
library(BalancedSampling)
library(phangorn)
library(dplyr)
library(geonames)
library(stringr)

# Load command line arguments into variables
args <- commandArgs(trailingOnly = TRUE)
in_csv <- args[2]
in_fasta <- args[3]
out_csv <- args[4]
n <- as.integer(args[5])
seed <- as.integer(args[6])

# Read in the CSV and FASTA files
metadata <- read.csv(in_csv, header = TRUE, stringsAsFactors = FALSE)
fasta <- read.dna(in_fasta, as.character=TRUE, format="fasta", as.matrix=TRUE)
fa <- as.phyDat(fasta)
distance <- dist.hamming(fa)

# Set the seed for reproducibility
set.seed(seed)

# Get the number of sequences in the FASTA file
N <- length(fasta)

# get temporal distance matrix and collect lat long vals
temporal<- c()
for (i in 1:nrow(metadata)) {
  temporal<- c(temporal, as.numeric(difftime(as.Date(metadata$Collection_Date[i], format = "%Y-%m-%d"), as.Date("2020-01-01"), unit="days")) / 365)
}
metadata$timediff <- temporal
metadata$lat<-as.numeric(metadata$lat)
metadata$lng<-as.numeric(metadata$lng)

# bind distance matrices together
matr = as.matrix(distance)
matr<-cbind(matr, metadata$lat, metadata$lng, metadata$timediff)

# get the LCUBE subsample
p = rep(n/N, N)
s = lcube(p, matr, cbind(p))

# write out the metadata from the subsample and the corresponding FASTA sequences
write.csv(metadata[s,], out_csv, row.names = FALSE)
write.dna(fasta[s,], file = str_replace(in_fasta, ".fasta", "_lcube.fasta"), format = "fasta")
