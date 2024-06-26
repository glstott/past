# script that provides a list of IDs for an LCUBE subsample of sequences
# Input: a CSV file with columns "id", "lat", "lng", and "date"; a fasta file with sequences and ids that match the CSV file
# Output: a CSV file with columns "id", "lat", "lng", and "date" for the LCUBE subsample; a fasta file with sequences and ids that match the subsample CSV file

# Load necessary libraries
library(BalancedSampling)
library(phangorn)
library(dplyr)
# library(geonames)
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

# Adding some hard-coded testing comments because as an academic, I can get away with it
# in_csv<-"../test_files/HA_NorthAmerica_202401-20240507.csv"
# id_col <- "id"
# date_col <- "Collection_Date"
# in_fasta <- "../test_files/HA_NorthAmerica_202401-20240507.aligned.fasta"
# out_csv <- "../test_files/HA_NorthAmerica_202401-20240507_lcube.csv"
# n <- 100
# seed <- 12

# Read in the CSV and FASTA files
metadata <- read.csv(in_csv, header = TRUE, stringsAsFactors = FALSE)
fasta <- read.dna(in_fasta, as.character=TRUE, format="fasta", as.matrix=TRUE)
rownames(fasta) <- trimws(rownames(fasta))
fa <- as.phyDat(fasta)
distance <- dist.hamming(fa)
metadata <- metadata[match(names(fa), metadata[,id_col]),]

# Set the seed for reproducibility
set.seed(seed)

# Get the number of sequences in the FASTA file
N <- length(fa)

# get temporal distance matrix and collect lat long vals
temporal<- c()
for (i in 1:nrow(metadata)) {
  temporal<- c(temporal, as.numeric(difftime(as.Date(metadata[i, date_col], 
                                                     tryFormats = c("%Y-%m-%d", "%Y-%m", "%Y")), 
                                             as.Date("2020-01-01"), unit="days")) / 365)
}
metadata$timediff <- temporal
# metadata$lat<-as.numeric(metadata$lat)
# metadata$lng<-as.numeric(metadata$lng)

# bind distance matrices together
matr = as.matrix(distance)
print(matr[0:3,])
# selected_columns <- sample(ncol(matr), 4) # This chunk is commented out, but 
# matr[, selected_columns]                  # left in case we need to speed things up.
matr<-cbind(matr, select(metadata, -c(date_col, id_col, "location")))
matr<- scale(matr)

# scream if there is an NA
if(any(is.na(matr))) {
  stop("Error: The matrix contains NA values.")
}

# get the LCUBE subsample
p = rep(n/N, N)
s = lcube(p, matr, cbind(p))
# print("outputting results...")

# write out the metadata from the subsample and the corresponding FASTA sequences
write.csv(metadata[s,], out_csv, row.names = FALSE)
write.dna(fasta[metadata[s, id_col],], file = out_fasta, format = "fasta")
