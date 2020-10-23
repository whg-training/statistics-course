#practical 1 (with answers)

#Understanding the posterior distribution
#Source: Statistical rethinking, R. MacElreath (book "SR")

#Exercise 1: let's look at a common likelihood function: the Binomial distribution
#Context: you have studied the binomial distribution during your course. Try to reproduce the figures in the power point
# a) with n=15, p=0.2
# b) with n=15, p=0.8
# c) with n=15, p=0.5
# d) with n=40, p=0.2

#code for a
pa <- dbinom(x=seq(0,15,1),size=15,prob=0.2)#size=nb. trials, prob=probability of success,x=nb. of successes
names(pa) <- 0:15
barplot(pa,col="#98c8fd")

#how n and p influcence the shape of the distribution?

#Exercise 2
#Context: how much water covers the earth? imagine a version of the globe small enough to hold in your hands.
#Experience: you toss the globe up in the air, when you catch it you record whether or not the surface under
#your index finger of your right hand is water or land. You repeat the procedure n times.
#The true proportion of water is p
#A single toss of the globe has probability p of producing water and 1-p of producing land
#Each toss is independent

#Computing posterior distribution by grid approximation
#define number of points
n_pts<-100
#define grid
p_grid <-seq (0,1,length.out=n_pts)
#define prior
prior<-rep(1,n_pts)
#compute likelihood at each value in grid
likelihood <- dbinom(6,size=9,prob=p_grid)
#product prior * likelihood
unstd.posterior <- likelihood * prior
#standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

#plot the approximated posterior

plot(p_grid, posterior, type ="b",xlab="probability of water",ylab="posterior probability")
mtext(paste0(n_pts," points"))

#Next steps:
  
#1.change the number of points and plot the new graphs. What are the effect on the approximation of the posterior distribution?

#2.compare with the figure 2.5 p. 30 book SR (also in the power point presentation). Briefly commment your results.

#3. try to put different priors (one at a time). What are the results? Briefly comment your results
#a) alternative prior nb1
prior<-ifelse(p_grid < 0.5,0,1)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)
plot(p_grid, posterior, xlab="probability of water",ylab="posterior probability")
mtext(paste0(n_pts," points"))

#b) alternative prior nb2
prior<-exp(-5*abs(p_grid - 0.5))
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)
plot(p_grid, posterior, xlab="probability of water",ylab="posterior probability")
mtext(paste0(n_pts," points"))

#4. plot the prior, likelihood, and posterior (non standardised) in one plot with relevant colors so that you can compare and discuss the results with your colleague
#try to make a series of plots with varying priors (you can use the priors set in the previous exercise)
#example of plots

prior<-rep(1,n_pts)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

par(mfrow=c(2,2))
plot(p_grid, prior,col="red",xlab="probability of water",ylab="",ylim=c(0,1),type="p",pch=".",cex=2)
lines(p_grid,likelihood,col="green",type="h")
lines(p_grid, unstd.posterior)

prior2<-ifelse(p_grid < 0.5,0,1)
unstd.posterior <- likelihood * prior2
posterior <- unstd.posterior / sum(unstd.posterior)
plot(p_grid, prior2,col="red",xlab="probability of water",ylab="",ylim=c(0,1),type="p",pch=".",cex=2)
lines(p_grid,likelihood,col="green",type="h")
lines(p_grid, unstd.posterior)

prior3<-exp(-5*abs(p_grid - 0.5))
unstd.posterior <- likelihood * prior3
posterior <- unstd.posterior / sum(unstd.posterior)
plot(p_grid, prior3,col="red",xlab="probability of water",ylab="",ylim=c(0,1),type="p",pch=".",cex=2)
lines(p_grid,likelihood,col="green",type="h")
lines(p_grid, unstd.posterior)

