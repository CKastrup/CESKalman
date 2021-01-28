# ===================================================================
# ------ FUNCTION FOR BUILDING A STATE SPACE REPRESENTATION ---------
# ===================================================================

# ===================================================================
#' @author Christian Sandholm Kastrup <CST@dreamgruppen.dk> , Anders Farver Kronborg <ANK@dreamgruppen.dk> and Peter Philip Stephensen <PSP@dreamgruppen.dk>
# by Christian Sandholm Kastrup
# Latest update: 17-12-2020

# ===================================================================
#' @title Function for building a state space representation to be used in the CESKalman function
#'
#' @description This function builds a state space representation to
#' be used in the CESKalman_Estimation, CESKalman, and CESKalman_Bootstrap functions.
#'
#' @details This function builds a state space representation to
#' be used in the CESKalman_Estimation, CESKalman, and CESKalman_Bootstrap functions. It is in particular suitable for production functions
#' with two different production factors where one of the factors is persistent, but can be used in any CES function with two factors.
#'
#' The estimated function is the Error-Correction model from Kronborg et al (2019):
#' \eqn{\Delta s_{t}=\alpha(s_{t-1}-(1-\sigma)p_{t-1}-\mu_{t-1})+\sum_{i=0}^{nlags}\kappa_{i}\Delta p_{t-i}+\sum_{i=1}^{nlags}\gamma_{i}\Delta s_{t-i}+\epsilon_t}.
#' \eqn{s_t} is the relative budget shares in logs (expenditure on factor 1 relative to factor 2), \eqn{p_t} is the relative prices in logs and \eqn{\mu_t=(\sigma-1)log(\Gamma_t)} where \eqn{\Gamma_t} is the relative
#' augmenting technical change. The process of \eqn{\mu_t}, the state variable, is an I(2) process: \eqn{\Delta \mu_t=\Delta \mu_{t-1}+\eta_t}.
#'
#' param can be a string of length 1 or 2. If lambda is freely estimated (set to NA) both param[1] (variance of error term in observation equation, \eqn{\epsilon_t})
#' and param[2] (variance of the error term in the relative augmenting technologies, \eqn{\eta_t}) are used. Else, only param[1] is used and \eqn{\lambda\Sigma^\eta=\Sigma^\epsilon}.
#' Note that the variances are specified such that \eqn{\Sigma^\epsilon=exp(param[1])} and \eqn{\Sigma^\eta=exp(param[2])} to ensure that the variances are positive.
#' We refer to Kronborg et al (2019) for further description of the methodology.
#'
#' X should be a matrix with all the explanatory variables to be used (all in logs). The dimension is TxK, where T
#' is time periods and K are the explanatory variables (not counting technical change as one). K=3 when nlags=0, =5 when nlags=1 and =7 when nlags=2. The column ordering should be:
#' 1: lags of relative budget shares in log-levels. 2: lags of the relative prices in log-levels. 3: First-difference of log relative prices.
#' 4: First-differences of log relative prices lagged one period. 5: First-differences of log relative budget shares lagged one period.
#' 6: First-differences of log relative prices lagged two periods. 7: First-differences of log relative budget shares lagged two periods.
#' When using the CESKalman function, data will be ordered in this way automatically.
#'
#' @param param Parameter values of the observation and state variance (see details)
#' @param X A matrix of all explanatory variables to be used (See details)
#' @param nlags Number of lags used of relative prices and relative budget shares. One of 0,1,2.
#' @param lambda The inverse of the signal-to-noise ratio. Can be any positive value larger than zero. If set to NA it is estimated
#' @param Leontief TRUE or FALSE (default). Should the elasticity of substitution be fixed at zero? If so, set =TRUE
#' @param sigma_init Initial value for the long run elasticity of substitution, \eqn{\sigma}
#' @param alpha_init Initial value for the adjustment parameter in the ECM, \eqn{\alpha}
#'
#' @usage build_SS(param,X,nlags,lambda,Leontief=FALSE,sigma_init,alpha_init)
#'
#' @return A list of class 'dlm' with all vectors and matrices to be used in the functions dlmMLE and dlmSmooth. We refer to Petris et al (2010)
#' for further description.
#'
#' @references Kronborg et al (2019) and Petris et al (2010)
#'
#' @export


build_SS<-function(param,X,nlags,lambda,Leontief=FALSE,sigma_init,alpha_init){


  if(nlags>2){
    stop("nlags has to be between zero and two")
  }

  if(!is.na(lambda)){
    if(lambda<0){
    stop("lambda must be NA or a positive value")
    }
  }

  if(!Leontief%in%c(FALSE,TRUE)){
    stop("Leontief must be TRUE or FALSE")
  }

  library("dlm")

  xData = X   ## The explanatory variables

  ## The first step is to build the state-space representation of the unobserved component, mu

    ## Observation equation
    FF   <- cbind(1, 1, 1)
    JFF   <- cbind(0,0,0) ## =0 for constant. Integer k>0 is element x[t,k]

    V     <- exp(param[1])

    ## State equations
    GG   <- rbind(c(1,1,1), c(0,1,1), c(0,0,0))

    R    <- rbind(0,0,1) ## S. 93 DLM vignette

    if(is.na(lambda)){
      W_struc    <- exp(param[2])*(R %*% t(R))
    }else{
      W_struc    <- (exp(param[1])/lambda)*(R %*% t(R))
    }

    ## Initial values
    m0 <- c(X[1,1]-(1-sigma_init)*X[1,2],0,0)

    C0    <- diag(c(5,5,5))


  ## Next, we add the other explanatory variables based on nlags and Leontief

  if(nlags==0&Leontief==FALSE){  ## No lags of first differences and unrestricted elasticity
  ## Observation equation
  FF   <- cbind(FF,1,1,1)
  JFF   <- cbind(JFF,1, 2,3) ## =0 for constant. Integer k>0 is element x[t,k]

  ## State equations
  GG   <- bdiag(GG,1,1,1)

  W    <- bdiag(W_struc,0,0,0)

  ## Initial values

  m0 <- c(m0,alpha_init,(sigma_init-1)*alpha_init,0)
  C0    <- bdiag(C0,5,5,5)
}

  if(nlags==1&Leontief==FALSE){  ## One lag of first differences and unrestricted elasticity
  ## Observation equation
  FF   <- cbind(FF,1,1,1,1,1)
  JFF   <- cbind(JFF,1, 2,3,4,5) ## =0 for constant. Integer k>0 is element x[t,k]

  ## State equations

  GG   <- bdiag(GG,1,1,1,1,1)

  W    <- bdiag(W_struc,0,0,0,0,0)

 ## Initial values

  m0 <- c(m0,alpha_init,(sigma_init-1)*alpha_init,0,0,0)
  C0    <- bdiag(C0,5,5,5,5,5)

  }

  if(nlags==2&Leontief==FALSE){  ## Two lags of first-differences and unrestricted elasticity
  ## Observation equation
  FF   <- cbind(FF,1,1,1,1,1,1,1)
  JFF   <- cbind(JFF,1, 2,3,4,5,6,7) ## =0 for constant. Integer k>0 is element x[t,k]

  ## State equations
  GG   <- bdiag(GG,1,1,1,1,1,1,1)

  W    <- bdiag(W_struc,0,0,0,0,0,0,0)

 ## Initial values

  m0 <- c(m0,alpha_init,(sigma_init-1)*alpha_init,0,0,0,0,0)
  C0    <- bdiag(C0,5,5,5,5,5,5,5)

}


# ## Next three is only used if the estimated LR elasticity is negative

  if(nlags==0&Leontief==TRUE){  ## Zero lags of first-differences and elasticity=0
  ## Observation equation
  FF   <- cbind(FF,1,1)
  JFF   <- cbind(JFF,1,2) ## =0 for constant. Integer k>0 is element x[t,k]

  ## State equation

  GG <- bdiag(GG,1,1)

  W    <- bdiag(W_struc,0,0)

  ## Initial values

  m0 <- c(m0,alpha_init,0)
  C0    <- bdiag(C0,5,5)

}

  if(nlags==1&Leontief==TRUE){ ## One lag of first-differences and elasticity=0
  ## Observation equation
  FF   <- cbind(FF,1,1,1,1)
  JFF   <- cbind(JFF,1, 2,3,4) ## =0 for constant. Integer k>0 is element x[t,k]

 # State equation
  GG   <- bdiag(GG,1,1,1,1)


  W    <- bdiag(W_struc,0,0,0,0)

 ## Initial values

  m0 <- c(m0,alpha_init,0,0,0)
  C0    <- bdiag(C0,5,5,5,5)

}

  if(nlags==2&Leontief==TRUE){ ## Number of lags of first-differences and elasticity=0
  ## Observation equation
  FF   <- cbind(FF,1,1,1,1,1,1)
  JFF   <- cbind(JFF,1, 2,3,4,5,6) ## =0 for constant. Integer k>0 is element x[t,k]

  ## State equation

  GG   <- bdiag(GG,1,1,1,1,1,1)


  W    <- bdiag(W_struc,0,0,0,0,0,0)

 ## Initial values

  m0 <- c(m0,alpha_init,0,0,0,0,0)
  C0    <- bdiag(C0,5,5,5,5,5,5)
}

build.dlm <- dlm(list(m0=m0,C0=C0,FF=FF,V=V,GG=GG,W=W,JFF=JFF,X=xData))

return(build.dlm)

}

