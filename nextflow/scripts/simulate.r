library(ape)
library(diversitree)
library(phangorn)
library(treeio)

# load command line arguments into variables outtree, 
#      outhist, taxa number, and seed
args <- commandArgs(trailingOnly = TRUE)
outtree <- args[2]
n <- args[3]
if (length(args) > 3) {
  seed <- as.integer(args[4])
}

#instantiate variables
lambda             <- 9
base_trans_freq    <- 0.5
mu                 <- 3
sample_proportions <- c(0.1, 0.01, 0.1, 0.05, 0.01)
set.seed(seed)

# a few generating functions to make our parameter list more human-readable, 
#   using rule of thumb from Liu et al.: low-0.1,0.3, high-0.7,0.9
low    <- function(){runif(1, 0.1, 0.3)}
med    <- function(){runif(1, 0.3, 0.6)}
high   <- function(){runif(1, 0.6, 0.9)}

# Generate rates for pairwise combos of continents: NA-1,SA-2, E-3, AS-4,AF-5
#   For now, we'll use the same rates for both directions of transition
#   12 13 14 15
rateNA_SA <- base_trans_freq * high()
rateNA_EU <- base_trans_freq * low()
rateNA_AS <- base_trans_freq * med()
rateNA_AF <- base_trans_freq * low()
#   23 24 25
rateSA_EU <- base_trans_freq * high()
rateSA_AS <- base_trans_freq * high()
rateSA_AF <- base_trans_freq * med()
#   34 35
rateEU_AS <- base_trans_freq * med()
rateEU_AF <- base_trans_freq * high()
#   45
rateAS_AF <- base_trans_freq * high()

# print each rate to the console
cat("rateNA_SA: ", rateNA_SA, "\n", "rateNA_EU: ", rateNA_EU, "\n",
    "rateNA_AS: ", rateNA_AS, "\n", "rateNA_AF: ", rateNA_AF, "\n",
    "rateSA_EU: ", rateSA_EU, "\n", "rateSA_AS: ", rateSA_AS, "\n",
    "rateSA_AF: ", rateSA_AF, "\n", "rateEU_AS: ", rateEU_AS, "\n",
    "rateEU_AF: ", rateEU_AF, "\n", "rateAS_AF: ", rateAS_AF, "\n")

# Place rates along with other necessary variables into rate matrix
p   <- c(
  lambda, lambda, lambda, lambda, lambda,
  mu, mu, mu, mu, mu,
  # 12 13 14 15
  rateNA_SA, rateNA_EU, rateNA_AS, rateNA_AF,
  # 21 23 24 25
  rateNA_SA, rateSA_EU, rateSA_AS, rateSA_AF,
  # 31 32 34 35
  rateNA_EU, rateSA_EU, rateEU_AS, rateEU_AF,
  # 41 42 43 45
  rateNA_AS, rateSA_AS, rateEU_AS, rateAS_AF,
  # 51 52 53 54
  rateNA_AF, rateSA_AF, rateEU_AF, rateAS_AF
)

print("Rate has been matrixed...")
print(n)
print(p)
# generate tree and estimate discrete history
phy <- tree.musse(p, max.taxa = n, x0=1, max.t=Inf, include.extinct=TRUE)
print("The MUSSE has run.")
# h <- history.from.sim.discrete(phy, 1:5)
# print("History simulated...")

# plot tree
phy2 <- phy

zapped_edge_lengths <- zapsmall(node.depth.edgelength(phy)[phy$orig$idx2[match(phy$tip.label, phy$orig$name2)]], digits=4)

# Include dates as part of tip names
phy2$tip.label <- paste(phy$tip.label, phy$tip.state,
                        zapped_edge_lengths, sep="_")

# write tree to file
write.tree(phy2, outtree)

# # write history to file
# write.table(h$history, outhist, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)