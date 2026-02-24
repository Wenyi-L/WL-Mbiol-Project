#' Convert a matrix to a vector
#'
#' Converts a matrix to a vector by stacking columns
#'
#' @param nmat The matrix to be converted
#'
#' @return The vectorized form of the matrix
#' @export
#'
#' @examples
#' foo = matrix (1:4, 2, 2)
#' vec (foo)
vec <- function(nmat){
  matrix(nmat,ncol=1)
}

#' Convert a vector to a matrix
#'
#' Converts a vector to a matrix with user-defined dimensions.
#'
#' @param nvec The vector to be turned into a matrix
#' @param nrow The number of rows in the matrix.  Either nrow or ncol
#'   must be specified, but not both.
#' @param ncol The number of columns in the matrix.  Either nrow or ncol
#'   must be specified, but not both.
#'
#' @return The matrix form of the vector
#' @export
#'
#' @details If the vector has length N and the desired number of rows
#'   or columns is M, then N must be evenly divisible by M.  The
#'   matrix is filled by columns.
#'
#' @examples
#' foo = 1:6
#' unvec (foo, nrow=2)
#' unvec (foo, ncol=2)
unvec <- function(nvec,nrow=NULL,ncol=NULL){
	if(is.null(nrow)) return(matrix(nvec,ncol=ncol));
	if(is.null(ncol)) return(matrix(nvec,nrow=nrow));
}

#' Convert 4-d array to matrix
#'
#' Converts a 4-d transition array, such as the K array produced by
#' make_AxT, to a cross-classified transition matrix ("megamatrix"),
#' such as the A matrix produced by make_AxT
#'
#' @param A4 The 4-d transition array
#'
#' @return The cross-classified transition matrix
#' @export
#'
#' @seealso unfold, make_AxT
#' @examples
#' M = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0.1*(0:2)
#' B = mk_B (F, Fdist="Bernoulli")
#' mT = 6
#' out = make_AxT (B, M, mT)
#' K = out$K
#' A = flatten (out$K) ## should equal out$A
flatten <- function(A4) {
	dim(A4)<-rep(sqrt(length(A4)),2)
	return(A4)
}

#' Convert a megamatrix to an array
#'
#' Converts a cross-classified transition matrix ("megamatrix") to a
#' 4-d transition array, such as the K array produced by make_AxT
#'
#' @param A2 The megamatrix to be transformed
#' @param dim A vector containing the dimensions of the array to be returned
#'
#' @return The transition array
#' @export
#'
#' @seealso flatten, make_AxT
#' @examples
#' M = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0.1*(0:2)
#' B = mk_B (F, Fdist="Bernoulli")
#' mT = 6
#' out = make_AxT (B, M, mT)
#' A = out$A
#' K = unfold (A, dim=c(3,6,3,6)  ## should equal out$K
unfold <- function(A2,dim) {
    dim(A2) <- dim;
    return(A2)
}

#' Make B matrix
#'
#' Calculates B, the clutch size distribution matrix
#'
#' @param maxKids The maximum clutch size.  Optional, with a default value of
#'   20.
#' @param F the fecundity matrix.  F\[i,j\] is the expected number of size i
#'   offspring from a size j parent.
#' @param Fdist the clutch size distribution.  Currently supported values are
#'   "Poisson" and "Bernoulli."  Optional, with a default value of "Poisson".
#'
#' @return The clutch size distribution matrix B, where B\[i,j\] is the
#'   probability that a size j parent produces i-1 offspring in a single
#'   reproductive bout
#' @export
#'
#' @examples
#' F = matrix (0, 3, 3); F[1,] = 0.1*(0:2)
#' out = mk_B (F, Fdist="Bernoulli")
mk_B = function (F, maxKids=20, Fdist="Poisson") {
  bigmz = ncol(F)
  B = matrix (0, maxKids+1, bigmz)
  if (Fdist == "Poisson") {
    for (z in 1:bigmz)
      B[,z] = dpois (0:maxKids, lambda=sum(F[,z]))
  } else if (Fdist == "Bernoulli") {
    for (z in 1:bigmz)
      B[,z] = dbinom (0:maxKids, size=1, prob=sum(F[,z]))
  } else {
    stop ("mk_B: I don't recognize that option for Fdist.\n")
  }

  return (B)
}

## For post-breeding census models.  Given the way we make A, the size
## x #kids transition matrix, the columns of B need to be the clutch
## size distribution conditional on surviving to reproduce.
#' Make B matrix for a post-breeding census
#'
#' Calculates B, the clutch size distribution matrix, for a
#' post-breeding census
#'
#' @param maxKids The maximum clutch size.  Optional, with a default value of
#'   20.
#' @param F the fecundity matrix.  F\[i,j\] is the expected number of size i
#'   offspring from a size j parent.
#' @param M The state transition matrix.  States can be a single
#'   quantity, such as size or life history stage, or can be
#'   cross-classified, such as size x environment.  M\[i,j\] is the
#'   probability of transitioning from state j to state i.
#' @param Fdist the clutch size distribution.  Currently supported values are
#'   "Poisson" and "Bernoulli."  Optional, with a default value of "Poisson".
#'
#' @return The clutch size distribution matrix B, where B\[i,j\] is the
#'   probability that a size j parent produces i-1 offspring in a single
#'   reproductive bout
#' @export
#'
#' @examples
#' F = matrix (c(0,0,0, 0.5,0,0, 0,0,0), 3, 3)
#' M = matrix (c(0,0.5,0, 0,0,0.5, 0,0,0), 3, 3)
#' out = mk_BPostBreeding (M, F, maxClutchSize, Fdist="Bernoulli")
mk_BPostBreeding = function (M, F, maxKids=20, Fdist="Poisson") {
  bigmz = ncol(F)
  surv = colSums(M); die = 1 - surv
  ## fecCondSurv = colSums(F) / surv
  ## zeroMass = c(1, rep(0, maxKids))
  B = matrix (0, maxKids+1, bigmz)
  if (Fdist == "Poisson") {
    for (z in 1:bigmz)
      if (surv[z] == 0) {
        B[,z] = rep(0, maxKids+1); B[1,z] = 1
      } else {
        ## B[,z] = die[z]*zeroMass +
        ## surv[z]*dpois (0:maxKids, lambda=)
        B[,z] = dpois (0:maxKids, lambda=sum(F[,z]/surv[z]))
      }
  }
  else if (Fdist == "Bernoulli") {
    for (z in 1:bigmz)
      if (surv[z] == 0) {
        B[, z] = rep (0, maxKids+1); B[1,z] = 1
      } else {
       ## B[,z] = die[z]*zeroMass +
        ## surv[z]*dbinom (0:maxKids, size=1, prob=sum(F[,z])/surv[z])
        B[,z] = dbinom (0:maxKids, size=1, prob=sum(F[,z]/surv[z]))
      }
  } else {
    stop ("mk_BPostBreeding: I don't recognize that option for Fdist.\n")
  }

  return (B)
}

######################################################################
# Function to take B and M matrices, and compute the transition
# probabilities from size-class i and j total kids, to all size classes
# and l total kids. This returns a vector of zeros if (l-j) is < 0
# or above the assumed maxnumber of kids per year,  ncol(B) - 1.
######################################################################
#' Helper function used for calculating size x #kids transition
#' matrices
#'
#' Function to take clutch size distribution matrix (B) and survival
#' and growth or survival and growth and environment transition matrix
#' (M), and compute the transition probabilities from size-class i and
#' j total kids, to all size classes and l total kids.
#' @param l the value of total \#kids in the next time step
#' @param i the current state
#' @param j the current total \#kids
#' @param B The clutch size distribution matrix, where B\[m,n\] is the
#'   probability that a size n parent produces m-1 offspring in a
#'   single reproductive bout
#' @param M The state transition matrix.  States can be a single
#'   quantity, such as size or life history stage, or can be
#'   cross-classified, such as size x environment.  M\[i,j\] is the
#'   probability of transitioning from state j to state i.
#' @return A vector whose kth entry is the probability that a size i
#'   individual with total number of offspring j will produce (l - j)
#'   offspring in the current reproductive bout and will transition to
#'   size k.  If l-j is < 0, the return value is a vector of zeros.
#' @details Called by make_AxT
#' @seealso make_AxT
#' @examples
#' M = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' B = mk_B (F)
#' p_xT (5, 2, 3, B, M)
p_xT <- function(l, i, j, B, M) {
  bigmz <- ncol(M); maxKids <- nrow(B)-1;
  newKids <- (l-j);
  if((newKids < 0) | (newKids > maxKids)) {
    return(rep(0, bigmz))
  }else{
    return(M[,i]*B[newKids+1,i])
  }
}

##############################################################################
## Function to make the 2D iteration matrix A for a size-kids model
## based on the M and B matrices summarizing a size-structured
## IPM. Apart from B and M the only input is mT, dimension for T
## (so range of T is 0 to mT-1). The 4-D iteration array K is also returned.
##
## Iteration matrix is modified so individuals who get to the maximum
## values of T in the matrix stay there, but continue to grow/shrink/die
#############################################################################

#' Make the size x \#kids transition matrix
#'
#' Makes the 2-dim iteration matrix A for a size x \#kids model, where
#' \#kids is the total number of offspring to date.
#' @param B The clutch size distribution matrix, where B\[m,n\] is the
#'   probability that a size n parent produces m-1 offspring in a
#'   single reproductive bout
#' @param M The state transition matrix.  States can be a single
#'   quantity, such as size or life history stage, or can be
#'   cross-classified, such as size x environment.  M\[i,j\] is the
#'   probability of transitioning from state j to state i.
#' @param mT The dimension of the \#kids part of the 4-d array.
#' I.e. the maximum \#kids is mT-1.
#' @details Called by distLifespanCondR2 and calcDistLRO.
#' @return A list containing
#' * A: The 2-dim transition matrix for states cross-classified by size
#' (or stage or...) and total number of offspring so far
#' * K: The 4-dim transition array, with dimensions (next size, next
#' #kids, current size, current #kids).  K is modified so individuals
#' who get to the maximum values of #kids in the matrix stay there,
#' but continue to grow/shrink/die.
#' @examples
#' M = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' B = mk_B (F)
#' mT = 20
#' out = make_AxT (B, M, mT)
make_AxT <- function(B, M, mT) {
  bigmz=ncol(M); Kvals=array(0,c(bigmz,mT,bigmz,mT));
  for(z in 1:bigmz){ # initial size
    for(k in 1:mT){ # initial T
      for(kp in 1:(mT-1)){ # final T
        Kvals[,kp,z,k]=p_xT(kp,z,k,B,M)
      }
      for (zp in 1:bigmz)
        ## make last kids class absorbing
        Kvals[zp,mT,z,k] = M[zp,z] - sum(Kvals[zp,,z,k])
    }
  }
  ## make kids-class mT absorbing: stay there with prob=1
  Kvals[1:bigmz,1:mT,1:bigmz,mT] <- 0;
  Kvals[1:bigmz,mT,1:bigmz,mT] <- M;
  A <- Kvals; dim(A) <- c(bigmz*mT,bigmz*mT);

  return(list(A=A,K=Kvals))
}

############################################################################
## Function to make megamatrix from a list of "fast" matrices
## (e.g. list of survival/growth transition matrices or fecundity
## matrices for different environments) and a "slow" matrix, e.g. the
## environment transition matrix.  The cross-classified states cycle
## through all the states in the relevant fast matrix for the first
## state in slowMatrix, then again for the second state in slowMatrix,
## etc.
##
## fastMatrixList: a list of matrices.  There should be a fast matrix
## for every state in slow matrix.  E.g. a fecundity matrix for each
## environment.
##
## slowMatrix: a matrix
############################################################################

#' Make megamtrix
#'
#' Calculates the transition matrix for cross-classified states from
#' the transition matrices for each sub-classification:
#' i.e. calculates the megamatrix.
#' @param fastMatrixList A list of transition matrices for the state
#'   that cycles more rapidly in the cross-classified state.  There
#'   should be one matrix for each state in the slow matrix.  e.g., if
#'   the cross-classified state vector is (env. 1, stage 1; env. 1,
#'   stage 2; ... env. 1, stage n; env. 2, stage 1, ...), then stage
#'   is the fast state, and there should be a survival/growth matrix
#'   for each environment or a fecundity matrix for each environment.
#' @param slowMatrix A transition matrix for the slow states in the
#'   cross-classified state.  e.g., if the cross-classified state
#'   vector is (env. 1, stage 1; env. 1, stage 2; ... env. 1, stage n;
#'   env. 2, stage 1, ...), then environment is the slow state and
#'   slowMatrix is the environment transition matrix.
#' @return M, the transition matrix for the cross-classified state
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' M = makeM (Plist, Q)
#' bigF = makeM(Flist, Q)
makeM = function (fastMatrixList, slowMatrix) {
  dimSlow = dim(slowMatrix)[1]
  dimFast = dim(fastMatrixList[[1]])[1]

  if (length(fastMatrixList) != dimSlow)
    stop("makeM: There must be a fast matrix for each state in the slow matrix.\n")
    if (dim(slowMatrix)[2] != dimSlow)
    stop ("makeM: Matrices must be square.\n")
  if (dim(fastMatrixList[[1]])[2] != dimFast)
    stop ("makeM: Matrices must be square.\n")

  M = matrix (NA, dimFast*dimSlow, dimFast*dimSlow)

  for (i in 1:dimSlow) {
    for (j in 1:dimSlow) {
      M[(i-1)*dimFast + 1:dimFast, (j-1)*dimFast + 1:dimFast] =
        fastMatrixList[[j]]*slowMatrix[i,j]
    }
  }

  return (M)
}

## Debugging: Try this for post-breeding? #############
p_xTPostBreeding <- function(l, i, j, B, M) {
  bigmz <- ncol(M); maxKids <- nrow(B) - 1;
  newKids <- (l-j);
  if((newKids < 0) | (newKids > maxKids)) {
    return(rep(0, bigmz))
  }else{
    return (M[,i] * B[newKids+1,i])
  }
}

make_AxTPostBreeding <- function(B, M, mT) {
  bigmz=ncol(M); Kvals=array(0,c(bigmz,mT,bigmz,mT));
  for(z in 1:bigmz){ # initial size
    for(k in 1:mT){ # initial T
      for(kp in 1:(mT-1)){ # final T
        Kvals[,kp,z,k] = p_xTPostBreeding (kp, z, k, B, M)
      }
      for (zp in 1:bigmz)
        ## make last kids class absorbing
        Kvals[zp,mT,z,k] = M[zp,z] - sum(Kvals[zp,,z,k])
    }
  }
  ## make kids-class mT absorbing: stay there with prob=1
  Kvals[1:bigmz,1:mT,1:bigmz,mT] <- 0;
  Kvals[1:bigmz,mT,1:bigmz,mT] <- M;
  A <- Kvals; dim(A) <- c(bigmz*mT,bigmz*mT);

  return(list(A=A,K=Kvals))
}
