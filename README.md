# past

Phylogenetic Analysis Sampling Tool (PAST) is a tool which allows you to randomly select a balanced and well-spread subsample from a population of sequences prior to subsequent phylogenetic analysis. PAST makes use of the LCUBE algorithm developed by Grafström and Tillé, 2013.[1] 

## Why though?

With the explosion in available sequence data for phylogenetic analysis, we have moved from a problem of insufficient data to resolve evolutionary relationships to a problem of too much data to readily make use of more computationally intensive Bayesian methods. Furthermore, the data which we collect tend to be heavily influenced by whether it is sampled from a high-resource region or is from a host of concern/priority. As we down-sample our datasets, it's important to keep in mind the disparities between these hosts, locations, etc. Currently, there are sampling tools which select for a subset of distance measures, but there is not yet a tool which is able to create a balanced and well-spread sample across an arbitrary number of relevant variables. 

Note to @glstott: Include a demonstrative figure for balanced and well spread samples here...

### Existing landscape



## A rudimentary explanation of LCUBE

LCUBE is an extension of the cube method for sample selection, developed by Deville and Tillé in 2004, which generates a balanced random sample.[2] A balanced sample is one where the mean value for each balancing (also known as auxiliary) variable is the same in the population as it is in the selected sample. LCUBE adds the additional constraint of selecting a well-spread sample, i.e. LCUBE selects a balanced and well-spread sample, one that contains the extreme values for each balancing variable while preserving the mean values of the population for each balancing variable. 

### The cube representation of sample selection

Samples in a population can be represented geometrically by an n-dimensional space where n is the total size of the population. For a quick demonstrative example, we will use the case where our population is size 3 (represented by a pigeon, chicken, and human). If we desire a sample of size of size 2, then we could select any of the samples represented by a red dot along the vertices of the cube below.

<img src="figures/Sampling Explanation 1_n.png"  width="600"/>

Using this representation, we can then construct an intersecting plane which represents the system of balancing equations (seen below, along with a vector of inclusion probabilities). This system of equations produces a subspace, in this case, a plane, that represents what combination of individuals will maintain the population means. To select a balanced random sample, we would take a random walk beginning at the vector of inclusion probabilities and ending at one of the edges of the cube. However, this would lead to us having a sample with some inclusion probabilities being floating points instead of a binary decision. 

We then "round" these values using a linear program that optimizes the balancing cost of a given sample in a *landing phase*. 


<img src="figures/Sampling Explanation 3_n.png"  width="600"/>

### The localized part of LCUBE

<img src="figures/Sampling Explanation 4_n.png"  width="600"/>

LCUBE works by generating clusters of $p+1$ samples where $p$ is the number of balancing variables (sketched above). For each cluster of $p+1$ samples, it then performs a cube subsampling step. At each cube step, we only include the individuals who have a value of 1 at the end of the random walk. We wait to perform the *landing phase* until there are less than $p+1$ individuals left. At each cube step, at least one individual will be added to the sample since we are wandering into an edge along the cube. 


## Ok, but how do I use it?

### I/O

#### Inputs
* `in_csv`: Metadata CSV with unique identifier matching FASTA file's headers. This should have the following headers: "id", "lat", "lng", and "date"
* `in_fasta`: FASTA file
* `n`: subsample size
* `seed`: seed

#### Output

* `out_csv`: a text file without a header listing the IDs to be included in your balanced and well-spread subsample.

### Example usage

To get a sample of size 100 using seed 12...

```bash
Rscript lcube.r example.csv example.fasta out.txt 100 12
```


## References

1. Grafström, A., & Tillé, Y. (2013). Doubly balanced spatial sampling with spreading and restitution of auxiliary totals. Environmetrics, 24(2), 120–131. https://doi.org/10.1002/env.2194
2. Deville, J.-C., & Tillé, Y. (2004). Efficient Balanced Sampling: The Cube Method. Biometrika, 91(4), 893–912. http://www.jstor.org/stable/20441151

