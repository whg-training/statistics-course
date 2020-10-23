library(colorBrewer)

# ---------------------
# 1) modelling de novo mutation rate a single data point.
# --------------------
pois_data = c(42, 49, 40, 36, 50, 37, 39, 43, 54, 55)

prior_shape = 1
prior_rate = 1

x = seq(1,100,by=0.1)
y_prior = dgamma(x,shape=prior_shape,rate=prior_rate)

op <- par(mfrow = c(2,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,2,1) + 0.1)
plot(x,y_prior,ylab='density',type='l',ylim=c(0,0.3))

my_colours = colorRampPalette(c("blue", "red"))( 10)

for (i in 1:length(pois_data)){
  post_shape = prior_shape + sum(pois_data[1:i])
  post_rate = prior_rate + i
  
  y_post = dgamma(x,shape=post_shape,rate=post_rate)
  lines(x,y_post,col=my_colours[i])
  Sys.sleep(1) 
  }

our_prob = 1 - pgamma(45,shape=prior_shape+sum(pois_data),rate = prior_rate+length(pois_data))
our_prob

# change the prior
# --------------
prior_shape = 10
prior_rate = 0.2

y_prior = dgamma(x,shape=prior_shape,rate=prior_rate)

plot(x,y_prior,ylab='density',type='l',ylim=c(0,0.3))

my_colours = colorRampPalette(c("blue", "red"))( 10)
for (i in 1:length(pois_data)){
  post_shape = prior_shape + sum(pois_data[1:i])
  post_rate = prior_rate + i
  
  y_post = dgamma(x,shape=post_shape,rate=post_rate)
  lines(x,y_post,col=my_colours[i])
  Sys.sleep(1) 
}

our_prob = 1 - pgamma(45,shape=prior_shape+sum(pois_data),rate = prior_rate+length(pois_data))
our_prob


# ---------------------
# 2) height data
# --------------------
# Download the height data and then read it into R.
h_data = read.csv('~/Desktop/AzimPresentation/GMS_sampling/height_data.csv')
N = length(h_data$x)
mean(h_data$x)
sd(h_data$x)

sum(h_data$x <165) / N
sum(h_data$x > 180 & h_data$x < 190) / N



# ----------------------------
# 3) Monty Hall Problem Simulation.
# -----------------------------

# Do not swap:
# --------

# 
N = 1e4
playerdoor <- sample(3,N,replace=TRUE)
cardoor <- sample(3,N,replace=TRUE)
sum(cardoor == playerdoor)/N

# plot the results
op <- par(mfrow = c(2,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,2,1) + 0.1)
plot(cumsum(cardoor[1:100]==playerdoor[1:100])/(1:100),pch=20,ylim=c(0,1),ylab='estimated prob',xlab = 'simulation number')
abline(h=1/3,col='red')

plot(cumsum(cardoor==playerdoor)/(1:N),pch=20,ylim=c(0,1),ylab='estimated prob',xlab = 'simulation number')
abline(h=1/3,col='red')

# Always swap:
# -----------

# assume player always chooses door 1 
N = 1e4
doors <- 1:3
cardoor <- sample(doors,N,replace=TRUE)
playerdoor <- sample(doors,N,replace=TRUE)

# Monty opens a door (can’t be the player’s door or the car door) and player swaps.
result = rep(NA,N)
for (i in 1:N){
  # monty opens a door
  notpossible = unique(c(cardoor[i],playerdoor[i]))
  if (length(notpossible)==2) {
    montydoor = doors[-notpossible]
  } else {
    montydoor <- sample(doors[-notpossible],1)
  }
  
  # player switches.
  chosen = doors[-c(montydoor,playerdoor[i])]
  if (length(chosen)!=1){
    print(i)
    stop()
  }
  
  #get the result of the game
  result[i] <- (chosen == cardoor[i])
}

# plot the results

plot(cumsum(result[1:100])/(1:100),pch=20,ylim=c(0,1),ylab='estimated prob',xlab = 'simulation number')
abline(h=2/3,col='red')

plot(cumsum(result)/(1:N),pch=20,ylim=c(0,1),ylab='estimated prob',xlab = 'simulation number')
abline(h=2/3,col='red')


# function to play interactively
# ------------------

monty <- function() {
  
  doors <- 1:3# randomly pick where the car is
  cardoor <- sample(doors,1)
  
  # prompt player
  print("Monty Hall says ‘Pick a door, any door!’")
  
  # receive the player’s choice of door (should be 1,2, or 3)
  chosen <- scan(what = integer(), nlines = 1, quiet = TRUE)
  
  # pick Monty’s door (can’t be the player’s door or the car door)
  if (chosen != cardoor) montydoor <- doors[-c(chosen, cardoor)]
  else montydoor <- sample(doors[-chosen],1)
  
  # find out whether the player wants to switch doors
  print(paste("Monty opens door ", montydoor, "!", sep=""))
  print("Would you like to switch (y/n)?")
  reply <- scan(what = character(), nlines = 1, quiet = TRUE)
  
  # interpret what player wrote as "yes" if it starts with "y"
  if (substr(reply,1,1) == "y") chosen <- doors[-c(chosen,montydoor)]
  
  # announce the result of the game!
  if (chosen == cardoor) print("You won!")
  else print("You lost!")
  
}

# ---------------
# 4) Metropolis algorithm
# ----------------
pois_data = c(42, 49, 40, 36, 50, 37, 39, 43, 54, 55)

prior_shape = 1
prior_rate = 1

T = 1e4
lambda_old = 10
prop_sd = 1
lambda = rep(NA,T)
u = runif(T)
for (i in 1:T){
  lambda_star = rnorm(1,mean=lambda_old,sd=1)
  ratio = (lambda_star^(prior_shape+sum(pois_data)) * exp(-(prior_rate+length(pois_data))*lambda_star)) / 
                  (lambda_old^(prior_shape+sum(pois_data)) * exp(-(prior_rate+length(pois_data))*lambda_old))
  if (u[i] < ratio) {
    lambda[i] = lambda_star
  }else{
    lambda[i] = lambda_old
  }
  lambda_old = lambda[i]
}

# have to work in log likelihood scale.
# -------------------
T = 1e4
prior_shape = 1
prior_rate = 1


lambda_old = 10
lambda = rep(NA,T)
prop_sd = 4
u = runif(T)
for (i in 1:T){
  lambda_star = rnorm(1,mean=lambda_old,sd=prop_sd)
  ratio = ((prior_shape+sum(pois_data))*log(lambda_star) + (-(prior_rate+length(pois_data))*lambda_star)) - 
    ((prior_shape+sum(pois_data))*log(lambda_old) + (-(prior_rate+length(pois_data))*lambda_old)) 
  if (u[i] < exp(ratio)) {
    lambda[i] = lambda_star
  }else{
    lambda[i] = lambda_old
  }
  lambda_old = lambda[i]
}

plot(lambda,type='n',ylim=c(20,60))
lines(lambda)
autocorr.plot(lambda,auto.layout = F,lag.max = 100)


# choice of proposal: small moves.
# -------------------
T = 1e4
lambda_old = 10
lambda = rep(NA,T)
prop_sd = 0.1
u = runif(T)
for (i in 1:T){
  lambda_star = rnorm(1,mean=lambda_old,sd=prop_sd)
  ratio = ((prior_shape+sum(pois_data))*log(lambda_star) + (-(prior_rate+length(pois_data))*lambda_star)) - 
    ((prior_shape+sum(pois_data))*log(lambda_old) + (-(prior_rate+length(pois_data))*lambda_old)) 
  if (u[i] < exp(ratio)) {
    lambda[i] = lambda_star
  }else{
    lambda[i] = lambda_old
  }
  lambda_old = lambda[i]
}


plot(lambda,type='n',ylim=c(20,60))
lines(lambda)
autocorr.plot(lambda,auto.layout = F,lag.max = 100)


# choice of proposal: large moves.
# -------------------
T = 1e4
lambda_old = 10
lambda = rep(NA,T)
prop_sd = 100
u = runif(T)
for (i in 1:T){
  lambda_star = rnorm(1,mean=lambda_old,sd=prop_sd)
  if (lambda_star <= 0){
    lambda[i] = lambda_old
  }else{
    ratio = ((prior_shape+sum(pois_data))*log(lambda_star) + (-(prior_rate+length(pois_data))*lambda_star)) - 
      ((prior_shape+sum(pois_data))*log(lambda_old) + (-(prior_rate+length(pois_data))*lambda_old)) 
    if (u[i] < exp(ratio)) {
      lambda[i] = lambda_star
    }else{
      lambda[i] = lambda_old
    }
    lambda_old = lambda[i]
  }
}

plot(lambda,type='n',ylim=c(20,60))
lines(lambda)
autocorr.plot(lambda,auto.layout = F,lag.max = 100)


# has the chain reached its stationary distribution.
T = 1e3
lambda_old = 100
lambda_run_1 = rep(NA,T)
prop_sd = 4
u = runif(T)
for (i in 1:T){
  lambda_star = rnorm(1,mean=lambda_old,sd=prop_sd)
  ratio = ((prior_shape+sum(pois_data))*log(lambda_star) + (-(prior_rate+length(pois_data))*lambda_star)) - 
    ((prior_shape+sum(pois_data))*log(lambda_old) + (-(prior_rate+length(pois_data))*lambda_old)) 
  if (u[i] < exp(ratio)) {
    lambda_run_1[i] = lambda_star
  }else{
    lambda_run_1[i] = lambda_old
  }
  lambda_old = lambda_run_1[i]
}
plot(lambda_run_1,type='n')
lines(lambda_run_1)


# run 2
T = 1e3
lambda_old = 10
lambda_run_2 = rep(NA,T)
prop_sd = 4
u = runif(T)
for (i in 1:T){
  lambda_star = rnorm(1,mean=lambda_old,sd=prop_sd)
  ratio = ((prior_shape+sum(pois_data))*log(lambda_star) + (-(prior_rate+length(pois_data))*lambda_star)) - 
    ((prior_shape+sum(pois_data))*log(lambda_old) + (-(prior_rate+length(pois_data))*lambda_old)) 
  if (u[i] < exp(ratio)) {
    lambda_run_2[i] = lambda_star
  }else{
    lambda_run_2[i] = lambda_old
  }
  lambda_old = lambda_run_2[i]
}
lines(lambda_run_2,col=2)


# What is the accptence rate
library(coda)
mcmc.trace <- mcmc(lambda)
summary(mcmc.trace)
1-rejectionRate(mcmc.trace)
plot(mcmc.trace)
autocorr.plot(mcmc.trace)

final_chain = lambda_run_3[seq(1000,T,by=30)]

mean(final_chain)
sd(final_chain)

sum(final_chain>45)/length(final_chain)


# ----------------------
# 5) allele frequency MCMC
# ------------------
prior = function(p){
  if((p<0) || (p>1)){  # || here means "or"
    return(0)}
  else{
    return(1)}
}

likelihood = function(p, nAA, nAa, naa){
  return(p^(2*nAA) * (2*p*(1-p))^nAa * (1-p)^(2*naa))
}

p_sampler = function(nAA, nAa, naa, n_iter, p_start_val, pproposal_sd){
  p = rep(0,n_iter)
  p[1] = p_start_val
  for(i in 2:n_iter){
    currentp = p[i-1]
    newp = rnorm(1,currentp,pproposal_sd)
    A = prior(newp)*likelihood(newp,nAA,nAa,naa)/(prior(currentp) * likelihood(currentp,nAA,nAa,naa))
    if(runif(1)<A){
      p[i] = newp       # accept move with probabily min(1,A)
    } else {
      p[i] = currentp        # otherwise "reject" move, and stay where we are
    }
  }
  return(p)
}


z=p_sampler(50,21,29,10000,0.1,0.01)
plot(z,type='n')
lines(z)

x=seq(0,1,length=1000)
hist(z,prob=T,100,main='')
lines(x,dbeta(x,122, 80)) 


# -----------------------------------
# 7) allele frequency and the inbreading coef.
# ----------------------------------------
prior_p = function(p){
  if((p<0) || (p>1)){  # || here means "or"
    return(0)}
  else{
    return(1)}
}

prior_g = function(g){
  if((g<0) || (g>1)){  # || here means "or"
    return(0)}
  else{
    return(1)}
}

likelihood = function(p,g, nAA, nAa, naa){
  return((g*p+(1-g)*p^2)^nAA * ((1-g)*2*p*(1-p))^nAa * (g*(1-p)+(1-g)*(1-p)^2)^naa)
}

g_p_sampler = function(nAA, nAa, naa, n_iter, g_start_val, p_start_val, g_proposal_sd, p_proposal_sd){
  g = rep(0,n_iter)
  p = rep(0,n_iter)
  g[1] = g_start_val
  p[1] = p_start_val
  for(i in 2:n_iter){
    current_g = g[i-1]
    current_p = p[i-1]
    new_g = current_g + rnorm(1,0,g_proposal_sd)
    new_p = current_p + rnorm(1,0,p_proposal_sd)
    A = prior_p(new_p)*prior_g(new_g)*likelihood(new_p,new_g,nAA,nAa,naa)/(prior_p(current_p)*prior_g(current_g)*likelihood(current_p,current_g,nAA,nAa,naa))
    if(runif(1)<A){
      p[i] = new_p
      g[i] = new_g # accept move with probabily min(1,A)
    } else {
      p[i] = current_p        # otherwise "reject" move, and stay where we are
      g[i] = current_g
    }
  }
  return(list(g=g,p=p)) # return a "list" with two elements named f and p
}

# out = g_p_sampler(100,21,49,1e4,0.1,0.1,0.1,0.1)

out = g_p_sampler(50,21,29,1e4,0.1,0.1,0.1,0.1)

plot(out$g,type='n')
lines(out$g)
hist(out$g,30)


plot(out$p,type='n')
lines(out$p)
hist(out$p,30)

dev.off()
plot(out$g,out$p,type='n',xlab='g inbreeding coef',ylab = 'p')
lines(out$g,out$p)
