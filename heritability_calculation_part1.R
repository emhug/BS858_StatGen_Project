'''
Statistical genetics final project
Part 1

1) Compute appropriate statistics to assess whether grip strength is heritable
in your population. In your summary of your findings, comment on
whether your results have any impact on whether or not you expect to be able to replicate the
SNP-hand grip association findings of the UKBstudy. 

In the files for this project, we were given the grip strength data from 1000 sibpairs.
These sibpairs can be used to estimate heritability of grip strength.
'''
# Set working directory
#setwd('~/Documents/Statistical Genetics/Final_Project')
# Read in the data
gripdata <- read.csv("project_sib_pheno_and_RV_data.csv", as.is=TRUE)

# Heritability based on sibling pair
# Subset only the sibling data
sibs <- subset(gripdata, GRIP1 & !is.na(GRIP2))
sibsgrip <- c(sibs$GRIP1, sibs$GRIP2) # Vector containing only grip measurements from siblings

# h^2 = 2*(Intraclass correlation coefficient)
# Calculate mean, standard deviation, and number of observations
meangrip <- mean(sibsgrip)
sdgrip <- sd(sibsgrip)
N <- nrow(sibs)
# Adjust grips based on mean
adjsib1 <- sibs$GRIP1 - meangrip # First sibling
adjsib2 <- sibs$GRIP2 - meangrip # Second sibling
# Compute intraclass correlation coefficient
icc <- sum(adjsib1*adjsib2)/sdgrip^2/(N-1)
h2 <- 2*icc
# Based on these calculations, heritability is 0.25748.

# Sanity check with ANOVA analysis to estimate heritability
# Subset the grips of each siblings and put into single dataset
sib1 <- sibs[, c("famid", "GRIP1")] # sib1 dataset
sib2 <- sibs[, c("famid", "GRIP2")] # sib2 dataset
names(sib1) <- c("famid", "GRIP") # Rename columns for binding
names(sib2) <- c("famid", "GRIP")
sib_anova <- rbind(sib1, sib2) # Combine the data into a single dataset

# Run the ANOVA on the dataset and print results
print(summary(aov(GRIP ~ as.factor(famid), data=sib_anova)))
# Based on the printed results,
MSB <- 48.89
MSW <- 37.70
k <- 2
# The formula for ICC based on ANOVA output is:
# (MSB - MSW) / (MSB + (k-1)*MSW)
anova_icc <- (MSB - MSW) / (MSB + (k-1)*MSW)
anova_h2 <- anova_icc*2 # Calculate heritability 
print(anova_h2)
# Based on the ANOVA method, the estimate of heritability is 0.25846.

'''
These two estimates of heritability are similar. 
  ICC-based h2: 0.25748
  ANOVA-based h2: 0.25846
Therefore the heritability of grip strength in this population is about 25.8%.'''
