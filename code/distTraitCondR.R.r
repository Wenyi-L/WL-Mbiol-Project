#' Distribution of lifespan conditional on LRO
#'
#' Calculates the distribution of lifespan conditional on LRO in the
#' presence of environmental variation.
#'
#' @param Plist A list of survival/growth transition matrices.
#'   Plist\[\[q\]\]\[i,j\] is the probability of transitioning from
#'   state j to state i in environment q.
#' @param Flist A list of fecundity matrices.  Flist\[\[q\]\]\[i,j\]
#'   is the expected number of state i offspring from a state j parent
#'   in environment q
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param maxAge The maximum attainable age
#' @param percentileCutoff A value between 0 and 1.  Calculations are
#'   performed for values of LRO out to this percentile.  Optional.
#'   The default value is 0.99.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#'
#' @details Assumes a pre-breeding census (reproduction happens before
#'   survival and growth).  The details of this calculation can be
#'   found in Robin E. Snyder and Stephen P. Ellner.
#'   2024. "To prosper, live long: Understanding the sources of
#'   reproductive skew and extreme reproductive success in structured
#'   populations."  The American Naturalist 204(2),
#'   https://doi.org/10.1086/730557 and its online supplement.
#'
#' @return Returns a list containing the following:
#' * probLifespanCondR: A matrix whose \[i,j\]th entry is the
#'   probability that an individual with LRO = i-1 will have a
#'   lifespan of j years (i.e. age j-1).
#' * distKidsAtDeath: A vector whose jth entry is the probability of
#'   having an LRO of j-1.
#' * sdLifespanCondR: A vector whose jth entry is the standard
#'   deviation of lifespan conditional on having an LRO of j-1.
#' * CVLifespanCondR: A vector whose jth entry is the coefficient of
#'   variation of lifespan conditional on having an LRO of j-1.
#' * maxKidsIndex: calculations were performed out to values of LRO
#'   equal to maxKidsIndex-1.
#' * normalSurvProb: A vector whose jth entry is the unconditional
#'   probability of surviving j years (i.e. until age j-1).
#' @export
#'
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 12
#' maxLRO = 30
#' maxAge=20
#' out = distLifespanCondR2 (Plist, Flist, Q, c0, maxClutchSize,
#'   maxLRO, maxAge)
distLifespanCondR2 = function (Plist, Flist, Q,
                               c0, maxClutchSize, maxLRO, maxAge,
                               percentileCutoff = 0.99,
                               Fdist="Poisson") {
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001

  mT = maxLRO + 1
  mA = maxAge + 1

  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## We assume Poisson-distributed number of offspring.
  ## The columns of B should sum to 1 and they do.
  B = mk_B (F, maxClutchSize, Fdist)

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.
  out = make_AxT (B, M, mT)
  A = out$A
  mzA = bigmz*mT
  K = out$K

  #####################################################################
  ## making the bullet matrices
  #####################################################################

  ## Make "bullet" matrices for A.
  Fbullet = Qbullet = Gbullet = matrix (0, mzA, mzA)

  ## Sbullet updates survival
  Sbullet = diag (colSums(A))

  ## Fbullet updates number of kids.
  for (j in 1:mT) {
    for (i in j:(mT-1)) {
      for (z in 1:bigmz) {
        if (Fdist == "Poisson") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dpois (i-j, lambda=sum(F[,z]))
        } else if (Fdist == "Bernoulli") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dbinom (i-j, prob=sum(F[,z]), size=1)
        } else {
          stop ("Supported options for Fdist are 'Poisson' and 'Bernoulli'.")
        }
      }
    }
  }
  ## Make mT absorbing.  If you have #kids index < mT, the probability
  ## of going to mT while in stage z = 1 - all other kid transition
  ## probabilities with stage z.

  for (j in 1:(mT-1)) {
    for (z in 1:bigmz) {
      Fbullet[(mT-1)*bigmz + z, (j-1)*bigmz + z] =
        ##      1 - sum(Fbullet[1:((mT-1)*bigmz), (j-1)*bigmz + z])
        1 - sum(Fbullet[(0:(mT-2))*bigmz + z, (j-1)*bigmz + z])
    }
  }
  ## If you have #kids index = mT, you stay there.
  for (z in 1:bigmz)
    Fbullet[(mT-1)*bigmz + z, (mT-1)*bigmz + z] = 1

  ## Qbullet updates env.
  smallQbullet = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      smallQbullet[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = diag(mz)*Q[i,j]
    }
  }
  ## Qbullet is the same for each number of kids
  for (k in 1:mT)
    Qbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = smallQbullet

  ## Gbullet updates stage
  smallGbullet = matrix (0, bigmz, bigmz)
  for (q in 1:numEnv) {
    S = diag (colSums (Plist[[q]]))
    Sinv = diag (1 / diag(S))
    G = Plist[[q]] %*% Sinv
    smallGbullet[(q-1)*mz + 1:mz, (q-1)*mz + 1:mz] = G
  }
  for (k in 1:mT)
    Gbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = smallGbullet

  ## sanity check
  epsilon = 0.00001
  if (sum(abs(range(Qbullet %*% Gbullet %*% Sbullet %*% Fbullet - A)))
      > epsilon)
    ## warning("Qbullet %*% Gbullet %*% Sbullet %*% Fbullet differs
    ## substantially from A")
    cat ("Qbullet %*% Gbullet %*% Sbullet %*% Fbullet differs substantially from A")

  ####################################################################
  ## Now make the additionally extended space matrix that allows us to
  ## stop partway through a time step.  Our new prefix is es to
  ## indicate this additional extension.  ("es" = "extended space")
  ####################################################################

  esmzA = 4*mzA

  ##esbigF = esA = matrix (0, esmzA, esmzA)
  esA = matrix (0, esmzA, esmzA)
  esA[mzA + 1:mzA, 1:mzA] = Fbullet
  esA[2*mzA + 1:mzA, mzA + 1:mzA] = Sbullet
  esA[3*mzA + 1:mzA, 2*mzA + 1:mzA] = Gbullet
  esA[1:mzA, 3*mzA + 1:mzA] = Qbullet

  ##################################################################
  ## What is the distribution of #kids at death when reproduction
  ## comes before survival and growth?
  ##################################################################

  ## Create a survival probability vector
  surv = apply (esA, 2, sum);  die = 1-surv;

  ## And the fundamental matrix of esA
  fundesA <- solve(diag(ncol(esA))-esA)

  ## Make omega
  omega = matrix(0,nrow(fundesA),ncol(fundesA));

  for(j in 1:ncol(fundesA))
    omega[,j] = die*fundesA[,j]

  distAtBirth = matrix(0, bigmz, 4*mT)
  distAtBirth[,1] = m0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  ## distAtDeath = matrix(distAtDeath, nrow=bigmz) # state matrix
  ## distKidsAtDeath = apply(distAtDeath, 2, sum)[1:mT]
  distAtDeath = array (distAtDeath, dim=c(bigmz, mT, 4)) # state
  # matrix
  distKidsAtDeath = apply (distAtDeath, 2, sum)

  ## What threshold number of kids represents the 99th %ile of the LRO
  ## distribution?
  cumDistKidsAtDeath = cumsum (distKidsAtDeath)
  R90 = which (cumDistKidsAtDeath > 0.9)[1]
  R99 = which (cumDistKidsAtDeath > 0.99)[1]
  RCutoff = which (cumDistKidsAtDeath > percentileCutoff)[1]
  if (RCutoff > maxClutchSize)
    stop ("RCutoff = ", RCutoff, "but maxClutchSize = ", maxClutchSize, "\n")

  probLifespanCondR = matrix (0, RCutoff+2, mA)
  probThresholdOrMore = matrix (0, RCutoff+2, esmzA)
  probSurvAgeA = numeric (mA+1)
  probSuccessCondZ = matrix (0, RCutoff+1, esmzA)
  probSuccess = numeric (RCutoff+1)
  probLifespanCondR = matrix (0, RCutoff+1, mA)
  normalSurvProb = numeric (RCutoff+2)

  AaM0 = c(m0, rep(0, (mT-1)*bigmz))
  mzAZero = rep(0, mzA)
  extendedInit = c(AaM0, rep(mzAZero, 3))

  ## The case of R = 0 (kidsIndexThreshold = 1) #########
  probThresholdOrMore[1,] = 1
  esAa = diag (esmzA)
  esA4 = esA %*% esA %*% esA %*% esA
  for (a in 1:(mA+1)) {
    probSurvAgeA[a] = colSums (esAa) %*% extendedInit
    if (a >= 2)
      probLifespanCondR[1,a-1] =
        probSurvAgeA[a-1] - probSurvAgeA[a]
    ## Keep track of survival probabilities for unconditional kernel.
    ## we need to hit the last matrix 4 times to advance one time step
    esAa = esA4 %*% esAa
  }
  normalSurvProb = probSurvAgeA

  ## The case of R > 0
  ##for (kidsIndexThreshold in 2:(RCutoff+2)) {

  ## k is the actual number of kids for the threshold.  We start with
  ## the threshold for success at LRO = 1 kid.
  for (k in 1:(RCutoff+1)) {
    message ("'Success' = ", k, " offspring")

    ## Get esA conditional on success
    cutoff = k*bigmz
    transientStates = c(1:cutoff, mzA + 1:cutoff, 2*mzA + 1:cutoff,
                        3*mzA + 1:cutoff)
    out = makeCondKernel (esA, transientStates)
    ## condEsA = out$MCond
    probThresholdOrMore[k+1,] = out$q2Extended

    probSuccessCondZ[k,] =
      probThresholdOrMore[k,] - probThresholdOrMore[k+1,]
    probSuccess[k] = probSuccessCondZ[k,] %*% extendedInit

    message ("Pr(success) = ", probSuccess[k])

    condEsA = matrix (0, esmzA, esmzA)
    ## For any starting state with #kids > threshold, probSuccessCondZ
    ## is zero and I think condEsA should be zero too.  Only do the
    ## following for z s.t. probSuccessCondZ > 0.
    for (j in transientStates)
      condEsA[,j] = esA[,j] * probSuccessCondZ[k,] /
      probSuccessCondZ[k,j]

    ## Useful later
    condEsA4 = condEsA %*% condEsA %*% condEsA %*% condEsA

    ## Initial distribution conditional on having LRO = exactly
    ## kidsIndexThreshold-1
    initCondSuccess =
      probSuccessCondZ[k,] * extendedInit
    initCondSuccess = initCondSuccess / sum(initCondSuccess)

    ## What is the lifespan dist. given that you had at least
    ## numKidsThreshold kids?
    condEsAa = diag(esmzA)
    for (a in 1:(mA+1)) {
      probSurvAgeA[a] = colSums (condEsAa) %*% initCondSuccess
      ## Keep track of survival probabilities for unconditional kernel.
      ## we need to hit the last matrix 4 times to advance one time step
      condEsAa = condEsA4 %*% condEsAa
    }

    probLifespanCondR[k,] = -diff (probSurvAgeA)
  }  ## end loop over k

  ## Sanity check: passes if maxAge is large enough
  foo1 = apply (probLifespanCondR, 1, sum)
  foo2 = rep (1, nrow(probLifespanCondR))
  if (max(abs(range(foo1 - foo2))) > epsilon)
    ## warning ("Rows of probLifespanCondR should sum to 1 and they
    ## don't.")
    cat ("Rows of probLifespanCondR should sum to 1 and they don't.")

  meanLifespanCondR = sdLifespanCondR = CVLifespanCondR =
    numeric (RCutoff+1)

  for (k in 1:(RCutoff+1)) {
    meanLifespanCondR[k] = sum((1:mA)*probLifespanCondR[k,])
    varLifespanCondR = sum((1:mA)^2*probLifespanCondR[k,]) -
      meanLifespanCondR[k]^2
    sdLifespanCondR[k] = sqrt(varLifespanCondR)
    CVLifespanCondR[k] = sdLifespanCondR[k] / meanLifespanCondR[k]
  }

  return (list(probLifespanCondR=probLifespanCondR,
               distKidsAtDeath=distKidsAtDeath,
               sdLifespanCondR=sdLifespanCondR,
               CVLifespanCondR=CVLifespanCondR,
               maxKidsIndex=RCutoff+1,
               normalSurvProb=normalSurvProb)) }

#' Distribution of lifespan conditional on LRO
#'
#' Calculates the distribution of lifespan conditional on LRO in the presence of
#' environmental variation, assuming a post-breeding census.
#'
#' @param Plist A list of survival/growth transition matrices.
#'   Plist\[\[q\]\]\[i,j\] is the probability of transitioning from state j to
#'   state i in environment q.
#' @param Flist A list of fecundity matrices.  Flist\[\[q\]\]\[i,j\] is the
#'   expected number of state i offspring from a state j parent in environment q
#' @param Q The environment transition matrix.  Q\[i,j\] is the probability of
#'   transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution: c0\[j\] is
#'   the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param maxAge The maximum attainable age
#' @param percentileCutoff A value between 0 and 1.  Calculations are performed
#'   for values of LRO out to this percentile.  Optional. The default value is
#'   0.99.
#' @param Fdist The clutch size distribution.  The recognized options are
#'   "Poisson" and "Bernoulli".  Optional.  The default value is "Poisson".
#'
#' @details Assumes a post-breeding census (reproduction happens after survival
#'   and growth).  The details of this calculation can be found in Robin E.
#'   Snyder and Stephen P. Ellner. 2024. "To prosper, live long: Understanding
#'   the sources of reproductive skew and extreme reproductive success in
#'   structured populations."  The American Naturalist 204(2),
#'   https://doi.org/10.1086/730557 and its online supplement.
#'
#' @return Returns a list containing the following:
#' * probLifespanCondR: A matrix whose \[i,j\]th entry is the
#'   probability that an individual with LRO = i-1 will have a lifespan of j
#'   years (i.e. age j-1).
#' * distKidsAtDeath: A vector whose jth entry is the probability of
#'   having an LRO of j-1.
#' * sdLifespanCondR: A vector whose jth entry is the standard
#'   deviation of lifespan conditional on having an LRO of j-1.
#' * CVLifespanCondR: A vector whose jth entry is the coefficient of
#'   variation of lifespan conditional on having an LRO of j-1.
#' * maxKidsIndex: calculations were performed out to values of LRO
#'   equal to maxKidsIndex-1.
#' * normalSurvProb: A vector whose jth entry is the unconditional
#'   probability of surviving j years (i.e. until age j-1).
#' @export
#'
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 12
#' maxAge=20
#' maxLRO = 30
#' out = distLifespanCondR2PostBreeding (Plist, Flist, Q, c0,
#'   maxClutchSize, maxLRO, maxAge)
#'
distLifespanCondR2PostBreeding = function (Plist, Flist, Q,
                                           c0, maxClutchSize, maxLRO, maxAge,
                                           percentileCutoff = 0.99,
                                           Fdist="Poisson") {
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001
  epsilon = 0.00001

  mT = maxLRO + 1
  mA = maxAge + 1

  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz
  mzA = bigmz*mT

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  distKidsAtDeath = calcDistLROPostBreeding (Plist, Flist, Q, c0,
                                             maxClutchSize,
                                             maxLRO, Fdist)

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)
  
  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)
  
  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }
  
  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## We assume Poisson-distributed number of offspring.
  ## The columns of B should sum to 1 and they do.
  B = mk_BPostBreeding (M, F, maxClutchSize, Fdist)
  
  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.
  out = make_AxT (B, M, mT)
  A = out$A
  mzA = bigmz*mT
  
  ## What threshold number of kids represents the 99th %ile of the LRO
  ## distribution?
  cumDistKidsAtDeath = cumsum (distKidsAtDeath)
  R90 = which (cumDistKidsAtDeath > 0.9)[1]
  R99 = which (cumDistKidsAtDeath > 0.99)[1]
  RCutoff = which (cumDistKidsAtDeath > percentileCutoff)[1]
  if (RCutoff > maxClutchSize)
    stop ("RCutoff = ", RCutoff, "but maxClutchSize = ", maxClutchSize, "\n")

  probLifespanCondR = matrix (0, RCutoff+2, mA)
  probThresholdOrMore = matrix (0, RCutoff+2, mzA)
  probSurvAgeA = numeric (mA+1)
  probSuccessCondZ = matrix (0, RCutoff+1, mzA)
  probSuccess = numeric (RCutoff+1)
  probLifespanCondR = matrix (0, RCutoff+1, mA)
  normalSurvProb = numeric (RCutoff+2)

  AaM0 = c(m0, rep(0, (mT-1)*bigmz))
  mzAZero = rep(0, mzA)
  ## extendedInit = AaM0

  ## The case of R = 0 (kidsIndexThreshold = 1) #########
  probThresholdOrMore[1,] = 1
  Aa = diag (mzA)
    for (a in 1:(mA+1)) {
    probSurvAgeA[a] = colSums (Aa) %*% AaM0
    if (a >= 2)
      probLifespanCondR[1,a-1] =
        probSurvAgeA[a-1] - probSurvAgeA[a]
    ## Keep track of survival probabilities for unconditional kernel.
    ## we need to hit the last matrix 4 times to advance one time step
    Aa = A %*% Aa
  }
  normalSurvProb = probSurvAgeA

  ## The case of R > 0
  ##for (kidsIndexThreshold in 2:(RCutoff+2)) {

  ## k is the actual number of kids for the threshold.  We start with
  ## the threshold for success at LRO = 1 kid.
  for (k in 1:(RCutoff+1)) {
    message ("'Success' = ", k, " offspring")

    ## Get A conditional on success
    cutoff = k*bigmz
    transientStates = 1:cutoff
    out = makeCondKernel (A, transientStates)
    probThresholdOrMore[k+1,] = out$q2Extended

    probSuccessCondZ[k,] =
      probThresholdOrMore[k,] - probThresholdOrMore[k+1,]
    probSuccess[k] = probSuccessCondZ[k,] %*% AaM0

    message ("Pr(success) = ", probSuccess[k])

    condA = matrix (0, mzA, mzA)
    ## For any starting state with #kids > threshold, probSuccessCondZ
    ## is zero and I think condA should be zero too.  Only do the
    ## following for z s.t. probSuccessCondZ > 0.
    for (j in transientStates) {
      if (probSuccessCondZ[k,j] == 0) {
        condA[,j] = 0
      } else {
        condA[,j] = A[,j] * probSuccessCondZ[k,] /
          probSuccessCondZ[k,j]
      }
    }

    ## Initial distribution conditional on having LRO = exactly
    ## kidsIndexThreshold-1
    initCondSuccess =
      probSuccessCondZ[k,] * AaM0
    initCondSuccess = initCondSuccess / sum(initCondSuccess)

    ## What is the lifespan dist. given that you had at least
    ## numKidsThreshold kids?
    condAa = diag(mzA)
    for (a in 1:(mA+1)) {
      probSurvAgeA[a] = colSums (condAa) %*% initCondSuccess
      ## Keep track of survival probabilities for unconditional kernel.
      ## we need to hit the last matrix 4 times to advance one time step
      condAa = condA %*% condAa
    }

    probLifespanCondR[k,] = -diff (probSurvAgeA)
  }  ## end loop over k

  ## Sanity check: passes if maxAge is large enough
  foo1 = apply (probLifespanCondR, 1, sum)
  foo2 = rep (1, nrow(probLifespanCondR))
  
  if (max(abs(range(foo1 - foo2))) > epsilon)
    ##warning ("Rows of probLifespanCondR should sum to 1 and they
    ##don't.")
    cat ("Rows of probLifespanCondR should sum to 1 and they don't.")

  meanLifespanCondR = sdLifespanCondR = CVLifespanCondR =
    numeric (RCutoff+1)

  for (k in 1:(RCutoff+1)) {
    meanLifespanCondR[k] = sum((1:mA)*probLifespanCondR[k,])
    varLifespanCondR = sum((1:mA)^2*probLifespanCondR[k,]) -
      meanLifespanCondR[k]^2
    sdLifespanCondR[k] = sqrt(varLifespanCondR)
    CVLifespanCondR[k] = sdLifespanCondR[k] / meanLifespanCondR[k]
  }

  return (list(probLifespanCondR=probLifespanCondR,
               distKidsAtDeath=distKidsAtDeath,
               sdLifespanCondR=sdLifespanCondR,
               CVLifespanCondR=CVLifespanCondR,
               maxKidsIndex=RCutoff+1,
               normalSurvProb=normalSurvProb))
}

#' Distribution of lifetime reproductive output
#'
#' Calculates the distribution of lifetime reproductive output (LRO) the
#' presence of environmental variation.  Assumes a pre-breeding census
#' (reproductiion happens before survival and growth).
#'
#' Called by probTraitCondLRO.
#'
#' It is also possible to calculate the distribution of LRO
#' by cross-classifying states by stage and number of offspring
#' produced so far and calculating the state distribution at death.
#' If s(z') is the survival probability for cross-classified state z',
#' M(z', z) is the probability of transitioning from cross-classified
#' state z to cross-classified state z', and N is the fundamental
#' matrix for M: N = (I - M)^\{-1\}, then the state distribution at
#' death, conditional on starting in state z, is (1 - s(z')) * N(z',
#' z), where * denotes a Hadamard product (element-by-element
#' multiplication, not matrix multiplication).  (See, e.g., eq. 3.2.8
#' in Data-driven Modeling of Structured Populations: A Practical
#' Guide to the Integral Projection Model, by Stephen P. Ellner, Dylan
#' Z. Childs and Mark Rees, ed. 1, 2015.)  However, this implicitly assumes
#' that reproduction  happens after survival and growth, i.e. a
#' post-breeding census.
#'
#' @param Plist A list of survival/growth transition matrices.
#'   Plist\[\[q\]\]\[i,j\] is the probability of transitioning from state j to
#'   state i in environment q.
#' @param Flist A list of fecundity matrices.  Flist\[\[q\]\]\[i,j\] is the
#'   expected number of state i offspring from a state j parent in environment q
#' @param Q The environment transition matrix.  Q\[i,j\] is the probability of
#'   transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution: c0\[j\] is
#'   the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param Fdist The clutch size distribution.  The recognized options are
#'   "Poisson" and "Bernoulli".  Optional.  The default value is "Poisson".
#'
#' @return The distribution of lifetime reproductive output
#' @importFrom Matrix bdiag
#' @export
#'
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 10
#' maxLRO=20
#' out = calcDistLRO (Plist, Flist, Q, c0, maxClutchSize, maxLRO)
calcDistLRO = function (Plist, Flist, Q,
                        c0, maxClutchSize, maxLRO,
                        Fdist="Poisson") {
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001

  mT = maxLRO + 1

  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## The columns of B should sum to 1.
  B = mk_B (F, maxClutchSize, Fdist)

  ## Let user decide if the level of eviction is acceptable
  message("calcDistLRO: checking clutch size distribution sums.  min(colSums(B)) = ",
          min(colSums(B)), "and would ideally be 1.\n")

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.
  out = make_AxT (B, M, mT)
  message ("calcDistLRO: Making A...")
  A = out$A
  mzA = bigmz*mT
  K = out$K

  #####################################################################
  ## making the bullet matrices
  #####################################################################

  ## Make "bullet" matrices for A.
  Fbullet = matrix (0, mzA, mzA)

  ## Fbullet updates number of kids.
  for (j in 1:mT) {
    for (i in j:(mT-1)) {
      for (z in 1:bigmz) {
        if (Fdist == "Poisson") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dpois (i-j, lambda=sum(F[,z]))
        } else if (Fdist == "Bernoulli") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dbinom (i-j, prob=sum(F[,z]), size=1)
        } else {
          stop ("Supported options for Fdist are 'Poisson' and 'Bernoulli.'\n")
        }
      }
    }
  }
  ## Make mT absorbing.  If you have #kids index < mT, the probability
  ## of going to mT while in stage z = 1 - all other kid transition
  ## probabilities with stage z.

  for (j in 1:(mT-1)) {
    for (z in 1:bigmz) {
      Fbullet[(mT-1)*bigmz + z, (j-1)*bigmz + z] =
        1 - sum(Fbullet[(0:(mT-2))*bigmz + z, (j-1)*bigmz + z])
    }
  }
  ## If you have #kids index = mT, you stay there.
  for (z in 1:bigmz)
    Fbullet[(mT-1)*bigmz + z, (mT-1)*bigmz + z] = 1

  ## Mbullet updates everything else
  ## Why doesn't the calculation work when I use bdiag?  No idea.
  if (FALSE) {
    foo = list ()
    for (k in 1:mT)
      foo[[k]] = M
    Mbullet = bdiag (foo)
  } else {
    Mbullet = matrix (0, mzA, mzA)
    for (k in 1:mT)
      Mbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = M
  }


  ## sanity check
  epsilon = 0.00001
  if (sum(abs(range(Mbullet %*% Fbullet - A))) > epsilon)
    stop ("Mbullet %*% Fbullet is substantially different than A.")

  ####################################################################
  ## Now make the additionally extended space matrix that allows us to
  ## stop partway through a time step.  Our new prefix is es to
  ## indicate this additional extension.  ("es" = "extended space")
  ####################################################################

  ## We don't need a space that separately updates survival, growth,
  ## and environment.  We just need to do reproduction first, then
  ## the rest.
  esmzA = 2*mzA
  esA = matrix (0, esmzA, esmzA)
  esA[mzA + 1:mzA, 1:mzA] = Fbullet
  esA[1:mzA, mzA + 1:mzA] = Mbullet

  ##################################################################
  ## What is the distribution of #kids at death when reproduction
  ## comes before survival and growth?
  ##################################################################

  ## Create a survival probability vector
  surv = apply (esA, 2, sum);  die = 1-surv;

  ## And the fundamental matrix of esA
  message ("calcDistLRO: calculating fundamental matrix...\n")
  fundesA <- solve(diag(ncol(esA))-esA)

  ## Make omega
  omega = matrix(0,nrow(fundesA),ncol(fundesA));

  for(j in 1:ncol(fundesA))
    omega[,j] = die*fundesA[,j]

  ##  distAtBirth = matrix(0, bigmz, 4*mT)
  distAtBirth = matrix(0, bigmz, 2*mT)
  distAtBirth[,1] = m0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  distAtDeath = array (distAtDeath, dim=c(bigmz, mT, 2)) # state matrix
  distKidsAtDeath = apply (distAtDeath, 2, sum)

  return (distKidsAtDeath)
}

#' Distribution of lifetime reproductive output
#'
#' Calculates the distribution of lifetime reproductive output (LRO) the
#' presence of environmental variation.  Assumes a post-breeding census
#' (reproductiion happens after survival and growth).
#'
#' Called by probTraitCondLROPostBreeding.
#'
#' This function calculates the distribution of LRO by
#' cross-classifying states by stage and number of offspring produced
#' so far and calculating the state distribution at death.  If s(z')
#' is the survival probability for cross-classified state z', M(z', z)
#' is the probability of transitioning from cross-classified state z
#' to cross-classified state z', and N is the fundamental matrix for
#' M: N = (I - M)^\{-1\}, then the state distribution at death,
#' conditional on starting in state z, is (1 - s(z')) * N(z', z),
#' where * denotes a Hadamard product (element-by-element
#' multiplication, not matrix multiplication).  (See, e.g., eq. 3.2.8
#' in Data-driven Modeling of Structured Populations: A Practical
#' Guide to the Integral Projection Model, by Stephen P. Ellner, Dylan
#' Z. Childs and Mark Rees, ed. 1, 2015.)  This method implicitly
#' assumes that reproduction happens after survival and growth, i.e. a
#' post-breeding census.
#'
#' @param Plist A list of survival/growth transition matrices.
#'   Plist\[\[q\]\]\[i,j\] is the probability of transitioning from state j to
#'   state i in environment q.
#' @param Flist A list of fecundity matrices.  Flist\[\[q\]\]\[i,j\] is the
#'   expected number of state i offspring from a state j parent in environment q
#' @param Q The environment transition matrix.  Q\[i,j\] is the probability of
#'   transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution: c0\[j\] is
#'   the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param Fdist The clutch size distribution.  The recognized options are
#'   "Poisson" and "Bernoulli".  Optional.  The default value is "Poisson".
#'
#' @return The distribution of lifetime reproductive output
#' @export
#'
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 10
#' maxLRO = 30
#' out = calcDistLROPostBreeding (Plist, Flist, Q, c0, maxClutchSize,
#'   maxLRO)
calcDistLROPostBreeding = function (Plist, Flist, Q,
                                    c0, maxClutchSize, maxLRO,
                                    Fdist="Poisson") {
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz
  mT = maxLRO + 1

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## The columns of B should sum to 1.
  B = mk_BPostBreeding (M, F, maxClutchSize, Fdist)

  ## Let user decide if the level of eviction is acceptable
  message("calcDistLRO: checking clutch size distribution sums.  min(colSums(B)) = ",
          min(colSums(B)), "and would ideally be 1.\n")

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.
  out = make_AxT (B, M, mT)
  message ("calcDistLRO: Making A...")
  A = out$A
  mzA = bigmz*mT
  ## Create a survival probability vector
  surv = apply (A, 2, sum);  die = 1-surv;

  ## And the fundamental matrix of esA
  message ("calcDistLRO: Calculating the fundamental matrix of A...")
  fundA <- solve(diag(ncol(A)) - A)

  ## Make omega
  omega = matrix(0, nrow(A), ncol(A));

  for (j in 1:ncol(fundA))
    omega[,j] = die*fundA[,j]

  distAtBirth = matrix(0, bigmz, mT)
  distAtBirth[,1] = m0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  distAtDeath = matrix (distAtDeath, nrow=bigmz)
  distKidsAtDeath = apply (distAtDeath, 2, sum)
  distKidsAtDeath = distKidsAtDeath / sum(distKidsAtDeath)

  return (distKidsAtDeath)
}

#' Distribution of a trait conditional on LRO
#'
#' Calculates Prob(X | R), where X is a trait value and R is lifetime
#' reproductive output (LRO) in the presence of environmental
#' variation. 
#' @param PlistAllTraits A list of lists of survival/growth transition
#'   matrices. Plist\[\[x\]\]\[\[q\]\]\[i,j\] is the probability of
#'   transitioning from state j to state i in environment q when an
#'   individual has trait x.
#' @param FlistAllTraits A list of lists of fecundity  matrices.
#'   Flist\[\[x\]\]\[\[q\]\]\[i,j\] is the expected number of state i
#'   offspring from a state j, trait x parent in environment q.
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment
#'   i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param traitDist A vector whose jth entry is the unconditional
#' probability that the trait takes on the jth value
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#'
#' @details The details of this calculation can be found in Robin
#' E. Snyder and Stephen P. Ellner.  2018.  "Pluck or Luck: Does Trait
#' Variation or Chance Drive Variation in Lifetime Reproductive
#' Success?"  The American Naturalist 191(4).  DOI:
#' 10.1086/696125
#'
#' @return A list containing the following:
#' * probXCondR: a matrix whose i,jth component is the probability
#' that the trait has the jth value given that LRO = i-1
#' * probR: A vector whose jth entry is the unconditional probability
#' that LRO = j-1
#' * probRCondX: A matrix whose i,jth entry is the probability that
#' LRO = i-1 given that the trait takes the jth value.
#' @export
#'
#' @examples
#' PlistAllTraits = FlistAllTraits = list ()
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = 0.9*P1
#' PlistAllTraits[[1]] = list (P1, P2)
#' P1 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
#' P2 = 0.9*P1
#' PlistAllTraits[[2]] = list (P1, P2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = 0.8*F1
#' FlistAllTraits[[1]] = list (F1, F2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0.9*(0:2)
#' F2 = 1.1*F1
#' FlistAllTraits[[2]] = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' traitDist = rep(0.5, 2)
#' out = probTraitCondLRO (PlistAllTraits, FlistAllTraits, Q,
#'       c0, maxClutchSize=10, maxLRO=15, traitDist)
probTraitCondLRO = function (PlistAllTraits, FlistAllTraits, Q,
                             c0, maxClutchSize, maxLRO,
                             traitDist,
                             Fdist="Poisson") {
  numTraits = length (PlistAllTraits)
  mT = maxLRO + 1

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  mz = dim(PlistAllTraits[[1]][[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (x in 1:numTraits) {
      for (q in 1:numEnv)
        if (sum(colSums(FlistAllTraits[[x]][[q]]) > 1))
          stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
               q, " with trait ", x," sum to > 1 but clutch size is Bernoulli-distributed.")
    }
  }

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## P(LRO | trait value) for each trait value
  probRCondX = matrix (0, mT, numTraits)
  ## P(trait value | LRO)
  probXCondR = matrix (0, mT, numTraits)

  for (x in 1:numTraits) {
    cat ("Trait", x, "...\n")
    Plist = PlistAllTraits[[x]]
    Flist = FlistAllTraits[[x]]
    probRCondX[,x] = calcDistLRO (Plist, Flist, Q,
                                  c0, maxClutchSize, maxLRO,
                                  Fdist)
  }
  ## Sanity check
  epsilon = 0.00001
  foo1 = colSums(probRCondX)
  foo2 = rep (1, ncol(probRCondX))
  if (max(abs(range (foo1 - foo2))) > epsilon)
    ## warning ("Columns of probRCondX should sum to 1 and they
    ## don't.")
    cat ("Columns of probRCondX should sum to 1 and they don't.")

  ## Now use Bayes thm to get P(X | R). ################

  ## Get joint probability P(R, X)
  probRAndX = matrix (0, mT, numTraits)
  for (i in 1:numTraits)
    ## joint probability P(R, X) = P(R | X) P(X).
    probRAndX[,i] = probRCondX[,i] * traitDist[i]

  ## Marginalize over X to get P(R)
  probR = apply (probRAndX, 1, sum)
  if (abs(sum(probR) - 1) > epsilon)
    ## warning ("Pr(R) does not sum to 1.")
    cat ("Pr(R) does not sum to 1.")

  ## Condition on R
  for (x in 1:numTraits)
    probXCondR[,x] = probRAndX[,x]/probR

  return (list(probXCondR=probXCondR,
               probR=probR, probRCondX=probRCondX))
}

#' Distribution of a trait conditional on LRO
#'
#' Calculates Prob(X | R), where X is a trait value and R is lifetime
#' reproductive output (LRO), for a post-breeding census model.
#' @param PlistAllTraits A list of lists of survival/growth transition
#'   matrices. Plist\[\[x\]\]\[\[q\]\]\[i,j\] is the probability of
#'   transitioning from state j to state i in environment q when an
#'   individual has trait x.
#' @param FlistAllTraits A list of lists of fecundity  matrices.
#'   Flist\[\[x\]\]\[\[q\]\]\[i,j\] is the expected number of state i
#'   offspring from a state j, trait x parent in environment q.
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment
#'   i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param traitDist A vector whose jth entry is the unconditional
#' probability that the trait takes on the jth value
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#'
#' @details The details of this calculation can be found in Robin
#' E. Snyder and Stephen P. Ellner.  2018.  "Pluck or Luck: Does Trait
#' Variation or Chance Drive Variation in Lifetime Reproductive
#' Success?"  The American Naturalist 191(4).  DOI:
#' 10.1086/696125
#'
#' @return A list containing the following:
#' * probXCondR: a matrix whose i,jth component is the probability
#' that the trait has the jth value given that LRO = i-1
#' * probR: A vector whose jth entry is the unconditional probability
#' that LRO = j-1
#' * probRCondX: A matrix whose i,jth entry is the probability that
#' LRO = i-1 given that the trait takes the jth value.
#' @export
#'
#' @examples
#' PlistAllTraits = FlistAllTraits = list ()
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = 0.9*P1
#' PlistAllTraits[[1]] = list (P1, P2)
#' P1 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
#' P2 = 0.9*P1
#' PlistAllTraits[[2]] = list (P1, P2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = 0.8*F1
#' FlistAllTraits[[1]] = list (F1, F2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0.9*(0:2)
#' F2 = 1.1*F1
#' FlistAllTraits[[2]] = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' traitDist = rep(0.5, 2)
#' out = probTraitCondLROPostBreeding (PlistAllTraits, FlistAllTraits, Q,
#'       c0, maxClutchSize=10, maxLRO=15, traitDist)
probTraitCondLROPostBreeding = function (PlistAllTraits, FlistAllTraits, Q,
                             c0, maxClutchSize, maxLRO,
                             traitDist,
                             Fdist="Poisson") {
  numTraits = length (PlistAllTraits)
  mT = maxLRO + 1

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  mz = dim(PlistAllTraits[[1]][[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (x in 1:numTraits) {
      for (q in 1:numEnv)
        if (sum(colSums(FlistAllTraits[[x]][[q]]) > 1))
          stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
               q, " with trait ", x," sum to > 1 but clutch size is Bernoulli-distributed.")
    }
  }

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## P(LRO | trait value) for each trait value
  probRCondX = matrix (0, mT, numTraits)
  ## P(trait value | LRO)
  probXCondR = matrix (0, mT, numTraits)

  for (x in 1:numTraits) {
    cat ("Trait", x, "...\n")
    Plist = PlistAllTraits[[x]]
    Flist = FlistAllTraits[[x]]
    probRCondX[,x] = calcDistLROPostBreeding (Plist, Flist, Q,
                                              c0, maxClutchSize, maxLRO,
                                              Fdist)
  }
  ## Sanity check
  epsilon = 0.00001
  foo1 = colSums(probRCondX)
  foo2 = rep (1, ncol(probRCondX))
  if (max(abs(range (foo1 - foo2))) > epsilon)
    ##warning ("Columns of probRCondX should sum to 1 and they
    ##don't.")
    cat ("Columns of probRCondX should sum to 1 and they don't.")

  ## Now use Bayes thm to get P(X | R). ################

  ## Get joint probability P(R, X)
  probRAndX = matrix (0, mT, numTraits)
  for (i in 1:numTraits)
    ## joint probability P(R, X) = P(R | X) P(X).
    probRAndX[,i] = probRCondX[,i] * traitDist[i]

  ## Marginalize over X to get P(R)
  probR = apply (probRAndX, 1, sum)
  if (abs(sum(probR) - 1) > epsilon)
    ##warning ("Pr(R) does not sum to 1.")
    cat("Pr(R) does not sum to 1.")

  ## Condition on R
  for (x in 1:numTraits)
    probXCondR[,x] = probRAndX[,x]/probR

  return (list(probXCondR=probXCondR,
               probR=probR, probRCondX=probRCondX))
}

## Wrapper functions for when there's no env. var. #################

#' Distribution of lifespan conditional on LRO
#'
#' Calculates the distribution of lifespan conditional on LRO in the
#' absence of environmental variation.
#'
#' @param P The survival/growth transition matrix.
#'   P\[i,j\] is the probability of transitioning from
#'   state j to state i.
#' @param F The fecundity matrix.  F\[i,j\]
#'   is the expected number of state i offspring from a state j
#'   parent. 
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param maxAge The maximum attainable age
#' @param percentileCutoff A value between 0 and 1.  Calculations are
#'   performed for values of LRO out to this percentile.  Optional.
#'   The default value is 0.99.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#'
#' @details Assumes a pre-breeding census (reproduction happens before
#'   survival and growth).  The details of this calculation can be
#'   found in Robin E. Snyder and Stephen P. Ellner.
#'   2024. "To prosper, live long: Understanding the sources of
#'   reproductive skew and extreme reproductive success in structured
#'   populations."  The American Naturalist 204(2),
#'   https://doi.org/10.1086/730557 and its online supplement.  This
#'   is a wrapper function for distLifespanCondR2.
#'
#' @return Returns a list containing the following:
#' * probLifespanCondR: A matrix whose \[i,j\]th entry is the
#'   probability that an individual with LRO = i-1 will have a
#'   lifespan of j years (i.e. age j-1).
#' * distKidsAtDeath: A vector whose jth entry is the probability of
#'   having an LRO of j-1.
#' * sdLifespanCondR: A vector whose jth entry is the standard
#'   deviation of lifespan conditional on having an LRO of j-1.
#' * CVLifespanCondR: A vector whose jth entry is the coefficient of
#'   variation of lifespan conditional on having an LRO of j-1.
#' * maxKidsIndex: calculations were performed out to values of LRO
#'   equal to maxKidsIndex-1.
#' * normalSurvProb: A vector whose jth entry is the unconditional
#'   probability of surviving j years (i.e. until age j-1).
#' @export
#'
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' c0 = c(1,0,0)
#' maxClutchSize = 12
#' maxLRO = 30
#' maxAge=20
#' out = distLifespanCondR2NoEnv (P, F, c0, maxClutchSize,
#'   maxLRO, maxAge)
distLifespanCondR2NoEnv = function (P, F, c0, maxClutchSize,
                                   maxLRO, maxAge,
                                   percentileCutoff=0.99,
                                   Fdist="Poisson") {
  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  Plist = list (P)
  Flist = list (F)
  Q = matrix (1, 1, 1)
  
  out = distLifespanCondR2 (Plist, Flist, Q,
                            c0, maxClutchSize, maxLRO, maxAge,
                            percentileCutoff, Fdist)

  return (out)
}

#' Distribution of lifespan conditional on LRO
#'
#' Calculates the distribution of lifespan conditional on LRO in the presence of
#' environmental variation, assuming a post-breeding census.
#'
#' @param P The survival/growth transition matrix.
#'   P\[i,j\] is the probability of transitioning from
#'   state j to state i.
#' @param F The fecundity matrix.  F\[i,j\]
#'   is the expected number of state i offspring from a state j
#'   parent. 
#' @param c0 A vector specifying the offspring state distribution: c0\[j\] is
#'   the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param maxAge The maximum attainable age
#' @param percentileCutoff A value between 0 and 1.  Calculations are performed
#'   for values of LRO out to this percentile.  Optional. The default value is
#'   0.99.
#' @param Fdist The clutch size distribution.  The recognized options are
#'   "Poisson" and "Bernoulli".  Optional.  The default value is "Poisson".
#'
#' @details Assumes a post-breeding census (reproduction happens after survival
#'   and growth).  The details of this calculation can be found in Robin E.
#'   Snyder and Stephen P. Ellner. 2024. "To prosper, live long: Understanding
#'   the sources of reproductive skew and extreme reproductive success in
#'   structured populations."  The American Naturalist 204(2),
#'   https://doi.org/10.1086/730557 and its online supplement.  This
#'   is a wrapper function for distLifespanCondR2PostBreeding.
#'
#' @return Returns a list containing the following:
#' * probLifespanCondR: A matrix whose \[i,j\]th entry is the
#'   probability that an individual with LRO = i-1 will have a lifespan of j
#'   years (i.e. age j-1).
#' * distKidsAtDeath: A vector whose jth entry is the probability of
#'   having an LRO of j-1.
#' * sdLifespanCondR: A vector whose jth entry is the standard
#'   deviation of lifespan conditional on having an LRO of j-1.
#' * CVLifespanCondR: A vector whose jth entry is the coefficient of
#'   variation of lifespan conditional on having an LRO of j-1.
#' * maxKidsIndex: calculations were performed out to values of LRO
#'   equal to maxKidsIndex-1.
#' * normalSurvProb: A vector whose jth entry is the unconditional
#'   probability of surviving j years (i.e. until age j-1).
#' @export
#'
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' c0 = c(1,0,0)
#' maxClutchSize = 12
#' maxAge=20
#' maxLRO = 30
#' out = distLifespanCondR2PostBreedingNoEnv (P, F, c0, maxClutchSize,
#'   maxLRO, maxAge) 
distLifespanCondR2PostBreedingNoEnv = function (P, F,
                                           c0, maxClutchSize, maxLRO, maxAge,
                                           percentileCutoff = 0.99,
                                           Fdist="Poisson") {

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  Plist = list (P)
  Flist = list (F)
  Q = matrix (1, 1, 1)

  out = distLifespanCondR2PostBreeding (Plist, Flist, Q,
                                        c0, maxClutchSize, maxLRO, maxAge,
                                        percentileCutoff,
                                        Fdist)
  return (out)
}

#' Distribution of lifetime reproductive output
#'
#' Calculates the distribution of lifetime reproductive output (LRO) the
#' absence of environmental variation.  Assumes a pre-breeding census
#' (reproductiion happens before survival and growth).
#'
#' It is also possible to calculate the distribution of LRO
#' by cross-classifying states by stage and number of offspring
#' produced so far and calculating the state distribution at death.
#' If s(z') is the survival probability for cross-classified state z',
#' M(z', z) is the probability of transitioning from cross-classified
#' state z to cross-classified state z', and N is the fundamental
#' matrix for M: N = (I - M)^\{-1\}, then the state distribution at
#' death, conditional on starting in state z, is (1 - s(z')) * N(z',
#' z), where * denotes a Hadamard product (element-by-element
#' multiplication, not matrix multiplication).  (See, e.g., eq. 3.2.8
#' in Data-driven Modeling of Structured Populations: A Practical
#' Guide to the Integral Projection Model, by Stephen P. Ellner, Dylan
#' Z. Childs and Mark Rees, ed. 1, 2015.)  However, this implicitly assumes
#' that reproduction  happens after survival and growth, i.e. a
#' post-breeding census.
#'
#' @param P The survival/growth transition matrix.
#'   P\[i,j\] is the probability of transitioning from
#'   state j to state i.
#' @param F The fecundity matrix.  F\[i,j\]
#'   is the expected number of state i offspring from a state j
#'   parent. 
#' @param c0 A vector specifying the offspring state distribution: c0\[j\] is
#'   the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param Fdist The clutch size distribution.  The recognized options are
#'   "Poisson" and "Bernoulli".  Optional.  The default value is "Poisson".
#'
#' @return The distribution of lifetime reproductive output
#' @importFrom Matrix bdiag
#' @export
#'
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' c0 = c(1,0,0)
#' maxClutchSize = 10
#' maxLRO=20
#' out = calcDistLRONoEnv (P, F, c0, maxClutchSize, maxLRO)
calcDistLRONoEnv = function (P, F, c0, maxClutchSize, maxLRO,
                        Fdist="Poisson") {

  Plist = list (P)
  Flist = list(F)
  Q = matrix (1, 1, 1)

  out = calcDistLRO (Plist, Flist, Q, c0, maxClutchSize, maxLRO)

  return (out)
}

#' Distribution of lifetime reproductive output
#'
#' Calculates the distribution of lifetime reproductive output (LRO) the
#' absence of environmental variation.  Assumes a post-breeding census
#' (reproductiion happens after survival and growth).
#'
#' This function calculates the distribution of LRO by
#' cross-classifying states by stage and number of offspring produced
#' so far and calculating the state distribution at death.  If s(z')
#' is the survival probability for cross-classified state z', M(z', z)
#' is the probability of transitioning from cross-classified state z
#' to cross-classified state z', and N is the fundamental matrix for
#' M: N = (I - M)^\{-1\}, then the state distribution at death,
#' conditional on starting in state z, is (1 - s(z')) * N(z', z),
#' where * denotes a Hadamard product (element-by-element
#' multiplication, not matrix multiplication).  (See, e.g., eq. 3.2.8
#' in Data-driven Modeling of Structured Populations: A Practical
#' Guide to the Integral Projection Model, by Stephen P. Ellner, Dylan
#' Z. Childs and Mark Rees, ed. 1, 2015.)  This method implicitly
#' assumes that reproduction happens after survival and growth, i.e. a
#' post-breeding census.
#'
#' @param P The survival/growth transition matrix.
#'   P\[i,j\] is the probability of transitioning from
#'   state j to state i.
#' @param F The fecundity matrix.  F\[i,j\]
#'   is the expected number of state i offspring from a state j
#'   parent. 
#' @param c0 A vector specifying the offspring state distribution: c0\[j\] is
#'   the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param Fdist The clutch size distribution.  The recognized options are
#'   "Poisson" and "Bernoulli".  Optional.  The default value is "Poisson".
#'
#' @return The distribution of lifetime reproductive output
#' @export
#'
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 10
#' maxLRO = 30
#' out = calcDistLROPostBreedingNoEnv (P, F, c0, maxClutchSize,
#'   maxLRO)
calcDistLROPostBreedingNoEnv = function (P, F, c0, maxClutchSize, maxLRO,
                                         Fdist="Poisson") {
  Plist = list (P)
  Flist = list(F)
  Q = matrix (1, 1, 1)
  out = calcDistLROPostBreeding (Plist, Flist, Q, c0, maxClutchSize,
                                 maxLRO)
  return (out)
}

#' Distribution of a trait conditional on LRO
#'
#' Calculates Prob(X | R), where X is a trait value and R is lifetime
#' reproductive output (LRO) in the absence of environmental
#' variation. 
#' @param Plist A list of lists of survival/growth transition
#'   matrices. Plist\[\[x\]\]\[i,j\] is the probability of
#'   transitioning from state j to state i when an individual has trait x.
#' @param Flist A list of lists of fecundity  matrices.
#'   Flist\[\[x\]\]\[i,j\] is the expected number of state i
#'   offspring from a state j, trait x parent.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param traitDist A vector whose jth entry is the unconditional
#' probability that the trait takes on the jth value
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#'
#' @details The details of this calculation can be found in Robin
#' E. Snyder and Stephen P. Ellner.  2018.  "Pluck or Luck: Does Trait
#' Variation or Chance Drive Variation in Lifetime Reproductive
#' Success?"  The American Naturalist 191(4).  DOI:
#' 10.1086/696125
#'
#' @return A list containing the following:
#' * probXCondR: a matrix whose i,jth component is the probability
#' that the trait has the jth value given that LRO = i-1
#' * probR: A vector whose jth entry is the unconditional probability
#' that LRO = j-1
#' * probRCondX: A matrix whose i,jth entry is the probability that
#' LRO = i-1 given that the trait takes the jth value.
#' @export
#'
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
#' Plist = list(P1, P2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.9*(0:2)
#' Flist = list (F1, F2)
#' c0 = c(1,0,0)
#' traitDist = rep(0.5, 2)
#' out = probTraitCondLRONoEnv (Plist, Flist,
#'       c0, maxClutchSize=10, maxLRO=15, traitDist)
probTraitCondLRONoEnv = function (PlistAllTraits, FlistAllTraits,
                                  c0, maxClutchSize, maxLRO,
                                  traitDist,
                                  Fdist="Poisson") {
  P1list = list(Plist[[1]])
  P2list = list(Plist[[2]])
  PlistAllTraits = list(P1list, P2list)
  
  F1list = list(Flist[[1]])
  F2list = list(Flist[[2]])
  FlistAllTraits = list (F1list, F2list)

  Q = matrix (1, 1, 1)

  out = probTraitCondLRO (PlistAllTraits, FlistAllTraits, Q,
                          c0, maxClutchSize=10, maxLRO=15, traitDist)
  return (out)
}
