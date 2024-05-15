# script that provides a list of IDs for an LCUBE subsample of sequences
# Input: a CSV file with a unique identifier ID column, date column (YYYY-MM-DD format), and any number of numeric columns; a fasta file with sequences and ids that match the CSV file
# Output: a CSV file for the LCUBE subsample; a fasta file with sequences and ids that match the subsample CSV file

# Load necessary libraries
library(BalancedSampling)
library(phangorn)
library(dplyr)
library(stringr)

# Load command line arguments into variables
args <- commandArgs(trailingOnly = TRUE)
in_csv <- args[1]
id_col <- args[2]
date_col <- args[3]
in_fasta <- args[4]
out_csv <- args[5]
out_fasta <- args[6]
n <- as.integer(args[7])
seed <- as.integer(args[8])

# Read in the CSV and FASTA files
metadata <- read.csv(in_csv, header = TRUE, stringsAsFactors = FALSE)
fasta <- read.dna(in_fasta, as.character=TRUE, format="fasta", as.matrix=TRUE)
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

# bind distance matrices together
matr = as.matrix(distance)
# selected_columns <- sample(ncol(matr), 4) # This chunk is commented out, but 
# matr[, selected_columns]                  # left in case we need to speed things up. I recommend a more useful selection of relevant distance measures.
matr<-cbind(select(metadata, -c(date_col, id_col)), matr)
matr<- scale(matr)

# scream if there is an NA
if(any(is.na(matr))) {
  stop("Error: The matrix contains NA values.")
}

# get the LCUBE subsample
p = rep(n/N, N)
s = lcube(p, matr, cbind(p))

# write out the metadata from the subsample and the corresponding FASTA sequences
write.csv(metadata[s,], out_csv, row.names = FALSE)
write.dna(fasta[s,], file = str_replace(in_fasta, ".fasta", ".lcube.fasta"), format = "fasta")