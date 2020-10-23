#practical 1 (without answers)

#Understanding the posterior distribution
#Source: Statistical rethinking, R. MacElreath (book "SR")

#Exercise 1 (SR, Ch. 2.4.1): let's look at a common likelihood function: the Binomial distribution
#Context: you have studied the binomial distribution during your course.
#Try to reproduce the histograms in the power point (slides 3-4)
# a) with n=15, p=0.2
# b) with n=15, p=0.8
# c) with n=15, p=0.5
# d) with n=40, p=0.2

#How n and p influcence the shape of the distribution?

#Exercise 2
#Context: how much water covers the earth? imagine a version of the globe small enough to hold in your hands.
#Experience: you toss the globe up in the air, when you catch it you record whether or not the surface under
#your index finger of your right hand is water or land. You repeat the procedure n times.
#The true proportion of water is p
#A single toss of the globe has probability p of producing water and 1-p of producing land
#Each toss is independent

#Computing posterior distribution by grid approximation
#Steps:

#1.define number of points
#2. define grid (the values should go from 0 to 1 given that we focus on probability to get "water")
#3. define prior (define a flat prior with constant value of 1)
#4. compute likelihood at each value in grid (tip: use command dbinom)
#5. product prior * likelihood
#6. standardize the posterior, so it sums to 1

#plot the approximated posterior

#Next steps:

#1.change the number of points and plot the new graphs. What are the effect on the approximation of the posterior distribution?

#2.compare with the figure 2.5 p. 30 book SR (also in the power point presentation). Briefly commment your results.

#3. try to put different priors (one at a time). What are the results? Briefly comment your results
#a) alternative prior nb1
prior<-ifelse(p_grid < 0.5,0,1)

#b) alternative prior nb2
prior<-exp(-5*abs(p_grid - 0.5))

#4. plot the prior, likelihood, and posterior (non standardised) in one plot with relevant colors so that you can compare and discuss the results with your colleague
#try to make a series of plots with varying priors (you can use the priors set in the previous exercise)
#example of plots

