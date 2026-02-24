#####################################################
# Functions for luck analyses
    #Derived from Chrissy's code
#####################################################    

## Input needed
# Fmat, Umat and Amat
# mixdist

####################################################
#What I have deleted:
#offspring_weight
#contributions for mean/var/skewness of lifespan/LRO
#analysis for individuals who have reproduced at least once and related plotting functions
#sensitivity analysis

fundamental_matrix = function(P){
  solve( diag(ncol(P)) - P)
}

# Calculate the expected value of lifetime reproductive success
mean_LRO<- function(Umat, Fmat, mixdist=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  # offspring_weight removed — not needed for age-structured model
  
  ## Calculate Ex(R | birth size z)
  N<- exactLTRE::fundamental_matrix(Umat)
  expRCond_z<- rep(1,Nclasses)%*%Fmat%*%N
  
  if(!is.null(mixdist)){
    expR<- expRCond_z%*%mixdist
    return(expR)
  } else{
    return(expRCond_z)
  }
}

# Calculate the variance in lifetime reproductive success for a given population
# matrix, given a number of options: based on the Markov Chain method and moments
# repro_var can be: poisson, bernoulli, or fixed
# mixdist, if provided, should be a vector for combining across offspring types
var_LRO_mcr<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  
  
  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it
  
  # Z matrix operator:
  Zmat<- cbind(diag(1, Nclasses, Nclasses, names=FALSE), 0)
  # column matrix of 1's:
  OneVec<- rep(1, Nclasses+1)
  # The fundamental matrix: 
  Nmat<- exactLTRE::fundamental_matrix(Umat)
  # Take the sum of offspring:
  eff<- colSums(Fmat) #eff is the stage-specific offspring production
  
  # Calculate the raw moments of the reward matrix:
  R1<- OneVec %*% t(c(eff, 0)) # first moment is the stage-specific reproductive output
  if (repro_var %in% c("poisson", "Poisson")){
    R2<- R1 + (R1*R1)
    R3<- R1 + 3*(R1*R1) + (R1*R1*R1)
  } else if (repro_var %in% c("fixed", "Fixed")){
    R2<- R1*R1
    R3<- R1*R1*R1
  } else if (repro_var %in% c("bernoulli", "Bernoulli")){
    R2<- R1
    R3<- R1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }
  # Rtilde1 is the first raw moment of the reward matrix for only transient states
  Rtilde1<- Zmat%*%R1%*%t(Zmat)
  # Rtilde2 is the second raw moment of the reward matrix for only transient states
  Rtilde2<- Zmat%*%R2%*%t(Zmat)
  
  # Calculate the raw moments of LRO conditional on starting state:
  mu_prime1<- t(Nmat)%*%Zmat %*% t(Pmat*R1)%*%OneVec
  mu_prime2<- t(Nmat)%*% (Zmat %*% t(Pmat*R2)%*%OneVec + 2*t(Umat*Rtilde1) %*% mu_prime1)
  # mu_prime3<- t(Nmat)%*% (Zmat %*% t(Pmat*R3)%*%OneVec + 3*t(Umat*Rtilde2) %*% mu_prime1 + 3*t(Umat*Rtilde1) %*% mu_prime2)
  
  ## The first raw moment of LRO conditional on starting state is the expected value:
  expRCond_z<- mu_prime1
  
  ## Var(R | birth size z)
  varRCond_z<- mu_prime2 - (mu_prime1*mu_prime1)
  
  if(is.null(mixdist)){
    return(varRCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varR_within<- t(mixdist)%*%varRCond_z
    # variance in LRO due to differences among starting states:
    varR_between<- t(mixdist)%*%(expRCond_z^2) - (t(mixdist)%*%expRCond_z)^2
    # total variance in LRO, given the mixing distribution:
    varR<- varR_within + varR_between
    return(varR)
  }
}

var_LRO_ipmbook<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  
  #offspring_weights removed from the original function
  
  # Take the sum of offspring:
  betabar<- colSums(Fmat) #betabar is the average offspring production
  
  # Variance in per-capita number of offspring per time step:
  if (repro_var %in% c("poisson", "Poisson")){
    sigsqb<- betabar # poisson
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    sigsqb<- betabar*(1 - betabar) # bernoulli
  } else if (repro_var %in% c("fixed", "Fixed")){
    sigsqb<- 0
  }
  
  ## Calculate Ex(R | birth size z)
  expRCond_z<- mean_LRO(Umat, Fmat)
  
  ## Var(R | birth size z) (see ch. 3 of the IPM book, eq. 3.2.15 (p. 82 of pdf))
  N<- exactLTRE::fundamental_matrix(Umat)
  rbarPib = expRCond_z %*% Umat
  r2 = (sigsqb + (betabar)^2 + 2*betabar*rbarPib) %*% N
  varRCond_z = (r2 - expRCond_z^2)
  
  if(is.null(mixdist)){
    return(varRCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varR_within<- varRCond_z%*%mixdist 
    # variance in LRO due to differences among starting states:
    varR_between<- t(mixdist)%*%t(expRCond_z^2) - (t(mixdist)%*%t(expRCond_z))^2
    # total variance in LRO, given the mixing distribution:
    varR<- varR_within + varR_between
    return(varR)
  }
}

# Calculate skew in lifetime reproductive success:
skew_LRO<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  
  #offspring_weights removed from the original function
  
  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it
  
  # Z matrix operator:
  Zmat<- cbind(diag(1, Nclasses, Nclasses, names=FALSE), 0)
  # column matrix of 1's:
  OneVec<- rep(1, Nclasses+1)
  # The fundamental matrix: 
  Nmat<- exactLTRE::fundamental_matrix(Umat)
  # Take the sum of offspring:
  eff<- colSums(Fmat) #eff is the stage-specific offspring production
  
  # Calculate the raw moments of the reward matrix:
  R1<- OneVec %*% t(c(eff, 0)) # first moment is the stage-specific reproductive output
  if (repro_var %in% c("poisson", "Poisson")){
    R2<- R1 + (R1*R1)
    R3<- R1 + 3*(R1*R1) + (R1*R1*R1)
  } else if (repro_var %in% c("fixed", "Fixed")){
    R2<- R1*R1
    R3<- R1*R1*R1
  } else if (repro_var %in% c("bernoulli", "Bernoulli")){
    R2<- R1
    R3<- R1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }
  
  # Rtilde1 is the first raw moment of the reward matrix for only transient states
  Rtilde1<- Zmat%*%R1%*%t(Zmat)
  # Rtilde2 is the second raw moment of the reward matrix for only transient states
  Rtilde2<- Zmat%*%R2%*%t(Zmat)
  
  # Calculate the raw moments of LRO conditional on starting state:
  mu_prime1<- t(Nmat)%*%Zmat %*% t(Pmat*R1)%*%OneVec
  mu_prime2<- t(Nmat)%*% (Zmat %*% t(Pmat*R2)%*%OneVec + 2*t(Umat*Rtilde1) %*% mu_prime1)
  mu_prime3<- t(Nmat)%*% (Zmat %*% t(Pmat*R3)%*%OneVec + 3*t(Umat*Rtilde2) %*% mu_prime1 + 3*t(Umat*Rtilde1) %*% mu_prime2)
  
  ## The first raw moment of LRO conditional on starting state is the expected value:
  expRCond_z<- mu_prime1
  
  ## Second central moment: mu2 = Var(R | birth size z)
  varRCond_z<- mu_prime2 - (mu_prime1*mu_prime1)
  
  # third central moment of LRO, conditional on Z:
  mu3Cond_z<- mu_prime3 - 3*(mu_prime2*mu_prime1) + 2*(mu_prime1*mu_prime1*mu_prime1)
  # Skewness, conditional on Z:
  skewCond_z<- varRCond_z^(-3/2)*mu3Cond_z
  
  if(is.null(mixdist)){
    return(skewCond_z)
  } else{ # law of total cumulance for third central moment
    # expected value of third central moment:
    mu3_within<- t(mu3Cond_z) %*% mixdist
    # third central moment of the expected value:
    mu3_between<- t((mu_prime1-(t(mu_prime1)%*%mixdist)[1])^3) %*% mixdist
    # covariance of expected value and variance:
    covExpVar<- t(mu_prime1*varRCond_z)%*%mixdist - (t(mu_prime1)%*%mixdist) * (t(varRCond_z)%*%mixdist)
    # total third central moment:
    mu3_total<- mu3_within + mu3_between + 3*covExpVar
    
    # total variance:
    var_total<- t(varRCond_z)%*%mixdist + t(mu_prime1^2)%*%mixdist - (t(mu_prime1)%*%mixdist)^2
    
    # total skewness:
    skew_total<- (var_total)^(-3/2)*mu3_total
    return(skew_total)
  }
}

# Function for mean lifespan. I already wrote essentially the same function in 
# the exactLTRE package, but I want to add the mixing distribution option.
mean_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  ## Calculate Ex(R | current state)
  N<- exactLTRE::fundamental_matrix(Umat)
  expLCond_z<- rep(1,Nclasses)%*%N
  
  if(!is.null(mixdist)){
    expL<- expLCond_z%*%mixdist
    return(expL)
  } else{
    return(expLCond_z)
  }
}

# Calculate the variance in lifespan:
# note: this calculates the variance in the number of time steps!
var_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  # calculate the fundamental matrix
  N<- exactLTRE::fundamental_matrix(Umat)
  
  ## Calculate Ex(R | current state)
  expLCond_z<- mean_lifespan(Umat, mixdist = NULL)
  
  ## Var(L | current state) using eqn. 5.12 from Hal's book:
  eT<- matrix(data=1, ncol=Nclasses, nrow=1) # column vector of 1's
  varLCond_z<- eT %*% (2*N%*%N - N) - (expLCond_z)^2
  
  if(is.null(mixdist)){
    return(varLCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varL_within<- varLCond_z %*% mixdist 
    # variance in LRO due to differences among starting states:
    varL_between<- t(mixdist)%*%t(expLCond_z^2) - (t(mixdist)%*%t(expLCond_z))^2
    # total variance in lifespan, given the mixing distribution:
    varL<- varL_within + varL_between
    return(varL)
  }
}

# Function for skew in lifespan:
skew_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  # calculate the fundamental matrix
  N<- exactLTRE::fundamental_matrix(Umat)
  # row vector of 1's
  eT<- matrix(data=1, ncol=Nclasses, nrow=1)
  # identity matrix
  eye<- diag(x=1, nrow=Nclasses, ncol=Nclasses, names=FALSE)
  
  # calculate the moments of lifespan
  eta1<- mean_lifespan(Umat)
  eta2<- eta1%*%(2*N-eye)
  eta3<- eta1%*%(6*N%*%N-6*N+eye)
  
  # calculate the central moments:
  eta_hat2<- eta2-eta1*eta1
  eta_hat3<- eta3 - 3*eta1*eta2 + 2*eta1*eta1*eta1
  
  # calculate skew:
  skewCond_z<- eta_hat2^(-3/2)*eta_hat3
  
  if(is.null(mixdist)){
    return(skewCond_z)
  } else{ # law of total cumulance for third central moment
    # expected value of third central moment:
    eta_hat3_within<- eta_hat3 %*% mixdist
    # third central moment of the expected value:
    eta_hat3_between<- ((eta1-(eta1%*%mixdist)[1])^3) %*% mixdist
    # covariance of expected value and variance:
    covExpVar<- (eta1*eta_hat2)%*%mixdist - (eta1%*%mixdist) * (eta_hat2%*%mixdist)
    # total third central moment:
    eta_hat3_total<- eta_hat3_within + eta_hat3_between + 3*covExpVar
    
    # total variance:
    var_total<- eta_hat2%*%mixdist + (eta1)^2%*%mixdist - (eta1%*%mixdist)^2
    
    # total skewness:
    skew_total<- (var_total)^(-3/2)*eta_hat3_total
    return(skew_total)
  }
}



# Function for calculating the stable distribution:
stable_dist<- function(Amat){
  # Calculate the eigen values and vectors:
  eigz<- eigen(Amat)
  # Find the index of lambda:
  I<- which(Re(eigz$values)==max(Re(eigz$values)))
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,I]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
}

# Function for calculating the mixing distribution:
mixing_distro<- function(Amat, Fmat){
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,1]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
  # Multiply the stable distribution by Fmat to get a cohort of offspring:
  offspring<- Fmat%*%wmean
  # Rescale to sum to 1 (proportions):
  offspring<- offspring/sum(offspring)
  
  return(offspring)
}

# Function for coefficient of variation:
coeff_var<- function(mean, var){
  sqrt(var)/mean
}


lambda<- function(Amat){
  lambda<- Re(eigen(Amat, only.values = T)$values[1])
  return(lambda)
}


# Function for adding transparency to colors:
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# A function for calculating S, as an iteroparity score, from Demetrius' entropy:
# This function will calculate the age-based lx and mx for the provided mixing 
# distribution (i.e. cohort structure), and then calculates iteroparity following
# the calculation provided in Rage package: https://github.com/jonesor/Rage/blob/main/R/entropy_d.R
iteroparity_d<- function(Umat, Fmat, mixdist, nage=300){
  Amat<- Umat+Fmat
  e<- matrix(1, nrow=1, ncol=ncol(Umat))
  
  # age-specific survival and fertility:
  surv_age<- matrix(ncol=ncol(Umat), nrow=nage)
  surv_age[1,]<- 1
  repro_age<- matrix(ncol=ncol(Umat), nrow=nage)
  Umati<- diag(ncol(Umat)); # not Umat, because surv_age[2,] results from Umat, not Umat%*%Umat. 
  
  # repro_age[1] is calculated outside of the loop. The loop below starts at i=2. 
  repro_age[1,] = e%*%Fmat; 
  
  for (i in 2:nage){
    Umati<- Umat%*%Umati 
    surv_age[i,]<- e%*%Umati
    repro_age[i,]<- e%*%Fmat%*%Umati/surv_age[i,]
  }
  lx<- surv_age%*%mixdist
  repro_age[is.na(repro_age)]<- 0
  mx<- repro_age%*%mixdist
  
  lxmx<- lx*mx; px = lxmx/sum(lxmx); ## SPE: the scaling happens HERE.
  # if lxmx == 0, log(lxmx) == -Inf; for entropy calc below, these -Inf can be
  #  converted to 0, because lim(x * log(x)) as x->0 is 0
  log_px <- log(px)
  log_px[px==0]<- 0
  
  this_iteroparity<- -sum(px * log_px)
  return(list(iteroparity=this_iteroparity,lx=lx,mx=mx)) 
}

# Function to calculate typical "clutch size" (i.e., reproductive output per 
# time step) at the stable stage distribution:
clutch_size<- function(Amat, Fmat){
  
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,1]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
  # Multiply the stable distribution by Fmat to get a cohort of offspring:
  offspring<- Fmat%*%wmean
  # Calculate the proportion of the population that is reproductively active:
  repro_active<- which(colSums(Fmat)>0)
  # The clutch size is the sum of the offspring cohort, divided by the 
  # reproductively active individuals:
  clutchsize<- sum(offspring)/sum(wmean[repro_active])
  
  return(clutchsize)
}

## Function to simulate an individual life history for a given Markov Chain
# Assumes Poisson distributed reproduction
SimulateMarkovIndiv<- function(Umat, Fmat, maxage=100, mixdist=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  
  Nclasses<- dim(Umat)[1]
  
  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it
  
  B = apply(Pmat,2,cumsum); ## Cumulative sums of each column, used in coin-toss for state transitions
  # For simulating reproduction, average repro output by state:
  effVec<- c(colSums(Fmat), 0) ## column sums of Fmat, plus no reproduction for dead individuals
  
  ## Why cumulative sums? To decide your next state, generate a Unif(0,1) random number
  ## and compare it to the i^th column of B, where i is your current state. If the random
  ## number exceeds none of the entries in B, you go to state 1. If the random number
  ## exceeds exactly one of the entries in B, you go to state 2. And so on -- this gives 
  ## correct transition probabilities conditional on your current state. 
  
  nt = maxage             ## Length of the simulation 
  states=numeric(nt+1);   ## Vector to hold states over time 
  repro<- numeric(nt+1)   ## Vector to hold reproductive output over time
  rd=runif(nt);           ## Do all the Unif(0,1) coin tosses at once 
  if(is.null(mixdist)){
    states[1] = 1;          ## Start in state 1 
  } else{
    birth_prob<- cumsum(mixdist)
    states[1]<- sum(birth_prob<runif(1))+1  ## Coin toss for starting state
  }
  
  for(i in 1:nt) {
    repro[i]<- rpois(1, effVec[states[i]])         ## Reproduce based on current state, Poisson distro 
    b=B[,states[i]];             ## Relevant cumulative probabilities, given current state  
    states[i+1]=sum(rd[i]>b)+1   ## ``Coin toss'' based on current state 
  }
  output<- list(states=states, repro=repro)
  return(output)
}

################################################################################
## Alternative measures for generation time:
################################################################################
gen_time_R0<- function(Umat, Fmat){
  Amat <- Umat + Fmat
  Nmat <- fundamental_matrix(Umat)
  R0 <- max(Re(eigen(Fmat %*% Nmat, only.values = T)$values))
  
  this_lambda<- max(Re(eigen(Amat, only.values=T)$values))
  
  gen_time<- log(R0)/log(this_lambda)
  return(gen_time)
}

gen_time_mu1_v<- function(Umat, Fmat){
  Nclasses<- ncol(Umat)
  e<- matrix(1, nrow=1, ncol=Nclasses)
  
  # Calculate the fundamental matrix for survival:
  N_0<- fundamental_matrix(Umat)
  # Calculate the stable stage distribution:
  eigz <- eigen(Umat+Fmat)
  ilambda <- which(Re(eigz$values) == max(Re(eigz$values)))
  w <- Re(eigz$vectors[, ilambda]); w = w/sum(w); 
  
  eigtz <- eigen(t(Umat+Fmat))
  ilambda <- which(Re(eigtz$values) == max(Re(eigtz$values)))
  v <- Re(eigtz$vectors[, ilambda]); v = v/sum(v); 
  v = matrix(v,nrow=1); 
  
  
  # calculation from Ellner 2018 (Am Nat)
  numerator<- (v%*%Fmat)%*%(N_0%*%N_0)%*%(Fmat%*%w)
  denominator<- (v%*%Fmat)%*%N_0%*%(Fmat%*%w)
  gen_time<- numerator/denominator
  
  return(gen_time)
}

gen_time_Ta<- function(Umat, Fmat){
  Amat<- Umat + Fmat
  eigz <- eigen(Amat)
  ilambda <- which(Re(eigz$values) == max(Re(eigz$values)))
  lambda <- max(Re(eigz$values))
  w <- Re(eigz$vectors[, ilambda])
  v <- Re(eigen(t(Amat))$vectors[, ilambda])
  gentime <- lambda*sum(v*w)/sum(v*(Fmat %*% w)) 
  return(gentime)
}
