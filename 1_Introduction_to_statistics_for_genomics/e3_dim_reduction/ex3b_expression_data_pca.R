## ex3b(i): Understanding gene expression data
## Looking at dimension reduction, using principal component analysis
# Two libraries needed (older versions of R maybe only need the first?)
install.packages("FactoMineR")
install.packages("factoextra")

library(FactoMineR)
library(factoextra)

# Get the data in
expression <- read.csv('microarray_all_expression.csv', header=TRUE, row.names = 1, check.names = FALSE)
treatment <- factor(expression$treatment, levels = c('Naive', 'IFN', 'LPS 2h', 'LPS 24h'))

# PCA
res.pca <- PCA(expression[, -1], scale.unit = FALSE, graph=FALSE)
eig.val <- get_eigenvalue(res.pca) # If we wanted to look at the eigenvalues, used to calculate % variance explained

# Plot a 'biplot'; the PCA with arrows showing which variables are driving the PCA
# First and second Component
fviz_pca_biplot(res.pca, axes = c(1,2), geom.ind = "point", pointshape = 21, pointsize = 2.5, 
select.var = list(contrib = 50), fill.ind = treatment, col.ind = "black", 
ggtheme = theme_classic(), repel = TRUE, title='PCA of naive and treated monocytes (1st & 2nd dimensions)') + 
ggpubr::fill_palette(c('#FF0000', '#008000', '#0000FF', '#FFD700'))

# Second and third Component
fviz_pca_biplot(res.pca, axes = c(2,3), geom.ind = "point", pointshape = 21, pointsize = 2.5, 
select.var = list(contrib = 50), fill.ind = treatment, col.ind = "black", 
ggtheme = theme_classic(), repel = TRUE, title='PCA of naive and treated monocytes (1st & 2nd dimensions)') + 
ggpubr::fill_palette(c('#FF0000', '#008000', '#0000FF', '#FFD700'))

# See how much variation is explained by each component (eigenvalue)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# See how much each sample contributes to these principal components
fviz_contrib(res.pca, choice = "ind", axes = 1:2)


## ex3b(ii) What is interesting about the top 50 components?
## Looking and testing for enrichment
# Optional installation, can use web plug in http://galahad.well.ox.ac.uk:3020/
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"

# Find the how much each genes contributes to component 1 and 2
vPC1 <- res.pca$var$coord[,1]
vPC2 <- res.pca$var$coord[,2]
vlabs <- rownames(res.pca$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- c("PC1", "PC2")

# Calculate the size of the arrow for component 1 and 2
vPCs$dist <- (vPCs$PC1**2 + vPCs$PC2**2)**0.5
# Order them by size
vPCs <- vPCs[order(-vPCs$dist),]
# Grab the top 50
interesting_vPCs <- head(vPCs,50)

# Get the IFN interesting genes
ifn <- rownames(interesting_vPCs[interesting_vPCs$PC1 > 0,])
for(i in ifn) {paste(i, '\n') %>% cat()}
eTerm <- xEnricherGenes(data=ifn, ontology="GOBP", RData.location=RData.location, test = "fisher") # Can also use "hypergeo", "fisher" or "binomial"
res <- xEnrichViewer(eTerm, top_num=30, sortBy="adjp",details=TRUE)
print(res)

# And the LPS interesing genes
lps <- interesting_vPCs[interesting_vPCs$PC1 < 0,]

# The genes linked with 24 hours
lps24 <- rownames(lps[lps$PC2 < 0,])
for(i in lps24) {paste(i, '\n') %>% cat()}
eTerm <- xEnricherGenes(data=lps24, ontology="GOBP", RData.location=RData.location, test = "fisher")
res <- xEnrichViewer(eTerm, top_num=30, sortBy="adjp",details=TRUE)
print(res)

# And two hours
lps2 <- rownames(lps[lps$PC2 > 0,])
for(i in lps2) {paste(i, '\n') %>% cat()}
eTerm <- xEnricherGenes(data=lps2, ontology="GOBP", RData.location=RData.location, test = "fisher")
res <- xEnrichViewer(eTerm, top_num=30, sortBy="adjp",details=TRUE)
print(res)


## ex3b(iii) Power calculations
## Example one: two independent groups (independent t-test)
# A proposed study wishes to investigate the effects of a hypertensive drug on patients with a genetic variant (experimental group) 
# compared to a patients without (control group). Previous studies show that the minimum clinically important difference is 15 mmHg 
# and a pilot study the pooled standard deviation (SD) is 20 mmHg and measurements were normally distributed. If we match the
# experimental group 1 to 1 with the control group, what sample size will we need for 80% power and alpha level p=0.05.

install.packages("pwr") # Useful website https://www.statmethods.net/stats/power.html
library(pwr)

eff = 15 / 20
sig = 0.05
pow = 0.8

print(pwr.t.test(type = "two.sample", d = eff, sig.level = sig, power = pow))

# After a further pilot study, we find the pooled standard deviation was 10mmHg
eff = 15 / 10
print(pwr.t.test(type = "two.sample", d = eff, sig.level = sig, power = pow))

## Example two: a post-hoc power calculation for two paired groups (paired t-test)
# Is  there  evidence  that  clofibrate  changes  the  mean cholesterol  level?  Cholesterol  is  to  be  measured 
# before and after receiving clofibrate. From previous studies,  a  mean  difference  of  40mg/dl  is  deemed 
# clinically significant, with pooled standard deviation of 50 with measurements normally distributed. What  sample 
# size  of  is  required for  alpha level p = 0.05 and  power 80%? 

eff = 40 / 50
sig = 0.05
pow = .8
  
print(pwr.t.test(type = "paired", d = eff, sig.level = sig, power = pow))

# How many samples would we need for alpha level p = 0.01
sig = 0.01
print(pwr.t.test(type = "paired", d = eff, sig.level = sig, power = pow))

# As a result of the poor response/dropout rate, only 12 patients were recruited. A mean difference of 50 and pooled 
# SD of 60  mg/dl were found. What is the power for this study?

eff = 50 / 60
sig = 0.05
num = 12
print(pwr.t.test(type = "paired", d = eff, sig.level = sig, n = num))

