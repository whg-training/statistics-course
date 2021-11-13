##R script to determine credible set of markers
##Input:
####1. regions determined by chromosomes and positions (Header: Locus SNP chr startPos endPos)
####2. Summary statistics (betas and SEs) from an association analysis (if the statistics are labelled differently from BETA and SE, please update their names in the script)
##Output:
####1. summary table with the number of selected set of SNPs and distance covered for each region
####2. association results for the selected set of SNPs in each region

## To run: Rscript --vanilla crediblesets.R alpha=0.99 data.file="results.file.name" regions.file="region.file.name" path="/well/abc/" 

 

#!/usr/bin/env Rscript
for (e in commandArgs(trailingOnly=TRUE))
{
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2]))
  {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    assign(ta[[1]][1],TRUE)
  }
}


if(!exists("alpha"))
{
  alpha <- 0.99
}

prior.var <- function()
{
	##theta's prior: N(0,W); return W
	return(0.04)
}

calc.BFs <- function(data)
{
	data$z = data$BETA / data$SE
	w = prior.var()

	## BF = data(P(y | H1) / P(y | H0)
	data$ABF = sqrt(data$SE^2 / (data$SE^2 + w)) * exp(data$z^2 / 2 * w / (data$SE^2 + w))
	return(data)
}

region.BF <- function(data)
{
	k = nrow(data)
	BFregion = 1/k * sum(data$ABF)
	return(BFregion)
}

posterior <- function(data)
{
	k = nrow(data)

	BFregion = region.BF(data)

	data$posterior = data$ABF / (k*BFregion)
	return(data)
}

get.posteriors <- function(data, region_l, region_u, region_chr)
{
	n = which(data$POS >= region_l & data$POS <= region_u & data$CHR == region_chr)
	data = data[n,]

	n=which(data$SE==0)
	if(length(n)!=0)
			data=data[-n,]
	data = calc.BFs(data)

	return(posterior(data))
}

credible.sets <- function(data, alpha, region_l, region_u, region_chr)
{
	data = get.posteriors(data, region_l, region_u, region_chr)
	data = data[sort.list(data$posterior, decreasing=TRUE),]

	sum_post=0
	credible=c()
	i = 1
	while(i <= nrow(data) & sum_post <= alpha)
	{
			sum_post = sum_post + data$posterior[i]
			credible = c(credible, as.character(data$SNP[i]))
			i = i+1
	}

	results = data[which(data$SNP %in% credible),]
   
	return(results)
}

determine.credible.sets <- function(alpha, regions.file, input.file)
{
         
	##provide file with the following columns: Locus, chr, startPos, endPos, SNP
	regions = read.table(regions.file, as.is=T, header=T)
	regions$n_snps = 0
	regions$distance = 0
	regions$start = 0
	regions$end = 0
	
	##Main data input file
	data = read.table(input.file, as.is=T, header=T)	
	## remove missing stats
	data<-subset(data, BETA!="NA" & SE!="NA" & P!="NA")
	## remove any duplicates in the data.
	data<- data[!duplicated(data[c("SNP")]),] 

	results<-vector("list", nrow(regions))
        
	for (i in 1:nrow(regions))
	{
		locus = as.character(paste(regions$Locus[i], regions$SNP[i], sep="_"))
		region_chr = as.numeric(regions$chr[i])
		region_l = as.numeric(regions$startPos[i])
		region_u = as.numeric(regions$endPos[i])

		results[[i]] = credible.sets(data, alpha, region_l, region_u, region_chr)
		results[[i]]$Locus=locus	   

		n_snps = nrow(results[[i]])
		distance = max(results[[i]]$POS) - min(results[[i]]$POS) 

		regions$n_snps[i] = n_snps
		regions$distance[i] = distance
		regions$start[i] = min(results[[i]]$POS)
		regions$end[i] = max(results[[i]]$POS)
		regions$total_posterior[i] = sum(results[[i]]$posterior)

	}

	resultCombined <-  do.call("rbind", results)
	regionsOrdered <- regions[order(regions$chr, regions$startPos),]
	mylist<-list(regionsOrdered, resultCombined) 
	return(mylist)
}

credibleSets=determine.credible.sets(alpha, regions.file, data.file)
write.table(credibleSets[[1]], file=paste(path, data.file, ".credibleset.summary", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(credibleSets[[2]], file=paste(path, data.file, ".credibleset.results", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

