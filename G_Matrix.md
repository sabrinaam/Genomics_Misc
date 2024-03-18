
### Load Data
```
id_geno2 <- scan("IDgenoC14_2.txt") # awk '{print $1}' genoC14.txt > IDgenoC14_2.txt
snps <- scan("genoC14_2.txt", what="list") # awk '{print $2}' genoC14.txt > genoC14_2.txt
```

### Creating W
- n: The number of SNP strings in the dataset.
- m: The number of SNPs per individual, determined by splitting the first SNP string into its individual characters (assumed to be digits here).
- W: Initializes a matrix to store integer representations of the SNPs, with n rows and m columns.

```
n <- length(snps)
m <- length(as.integer(strsplit(x = snps[1], split = '', fixed = T)[[1]] ))
W <- matrix (nrow=n, ncol=m, 10)
## The following loop converts SNP data from character strings to integer matrices, handling missing data represented by '5':
for (i in 1:n){
        W[i,] <- as.integer(strsplit(x=snps[i], split='', fixed=T)[[1]])
        W[i,] <- ifelse(W[i,] == 5, NA, W[i,])
        print(i)
    }
```
- The loop iterates over each SNP string, splits it into individual characters, converts them to integers, and stores them in the W matrix.
Replaces '5's with NA, assuming '5' indicates missing data.


### Imputation and MAF

- Imputes missing values (NAs) in the matrix W using random binomial draws, with the probability of success being the allele frequency p.
p is calculated as the sum of all alleles divided by twice the number of rows (individuals), representing the allele frequency.

```
set.seed(100)
p <- colSums(W, na.rm = TRUE) / (2*nrow(W)) # or colMean(X) / 2
for (j in 1:ncol(W)){
  W[,j] <- ifelse(is.na(W[,j]), rbinom(1, 2, prob = p[j]), W[,j])
}
```
- The following chunk updates allele frequencies and calculates the minor allele frequency (MAF).
Removes SNPs with MAF below 0.05, a common quality control step in genetic analyses.

```
p <- colSums(W) / (2*nrow(W))
maf <- ifelse(p > 0.5, 1-p, p)
index <- which(maf < 0.05)
length(index)
W2 <- W[, -index]
dim(W2) # 928 505367
```

### Constructing the Genomic Relationship Matrix (G1)

- Recalculates p using the filtered SNP matrix W2.
- Centers the SNP data by subtracting the mean.
- Computes the genomic relationship matrix G using VanRaden's method, which is the scaled transpose of the matrix multiplied by itself, adjusted by the allele frequencies.


```
p <- colSums(W2, na.rm = TRUE) / (2*nrow(W2)) # or colMean(W2) / 2
W_c <- scale(W2, center = TRUE, scale = FALSE)
G_C14 <- tcrossprod(W_c) / sum(2 * p * (1 - p))
```

### Save and Load

```
save(G_C14, file="G_C14.Rda")
load("G_C14.Rda")
dim(G_C14)
```
