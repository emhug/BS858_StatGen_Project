# Part 9 of the final project
# Test the loss of function rare variants in the sequenced individuals for association with handgrip strength.

# Set working directory
setwd("~/Documents/Statistical Genetics/Final_Project")

# Read in the data
data <- read.csv("project_sib_pheno_and_RV_data.csv", header=TRUE, as.is=TRUE)

# Subset for the first 5000 individuals
data <- data[1:5000, ]

# Calculate minor allele frequency
maf<-apply(data[,14:22],2,function(i)sum(i)/(2*length(i)))

# Determine weight of rare variant
weights<-1/sqrt(5000*maf*(1-maf))

# Determine the total weight of each individual
xxx<-sapply(1:9,function(i)data[,i+13]*weights[i])
data$MB <- apply(xxx,1,sum)
# Sanity check
# Mean and standard deviation
mean(data$MB)
sd(data$MB)

# testing association
# Madsen-Browning
mb_all <- glm(GRIP1 ~ age1+sex1+bmi1+hgt1+MB, data=data)
summary(mb_all)
# MB has a p-value of 0.0502 when population is included and a p-value of 0.0478 when population is not included

# CMC
data$CMC<-apply(data[,14:22],1,sum)
cmc <-glm(GRIP1 ~ age1+sex1+bmi1+hgt1+CMC, data=data)
summary(cmc)
# CMC has a p-value of 0.100 when population is included and a p-value of 0.0966 when population is not included