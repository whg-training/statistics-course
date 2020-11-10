
dataA = read.table( "datasetA.tsv", header = T )
dataA = read.table( "datasetB.tsv", header = T ) # ok, they are identical really

# We want to know about the effect of treatment on C!

# DATASET A:
# A is a therapy being tested
# B is measured cholesterol levels
# C is measured blood pressure

# DATASET B:
# A is a therapy being tested
# B is a measure of ethnic background / and or location of the individual (for example, computed using principal components from genome-wide genotypes)
# C is measured blood pressure


View(dataA)
View(dataB) # ok they are identical

# Fit linear models with A or A and B
lA = lm( C ~ A, data = dataA );
lAB = lm( C ~ A + B, data = dataA ); 

# Look at result
summary(lA)$coeff
summary(lAB)$coeff

#... so should we condition on B or not?
