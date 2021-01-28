# ===================================================================
# ------ FUNCTION FOR ESTIMATING A CES PRODUCTION FUNCTION ----------
# ===================================================================

# ===================================================================
#' @author Christian Sandholm Kastrup <CST@dreamgruppen.dk> , Anders Farver Kronborg <ANK@dreamgruppen.dk> and Peter Philip Stephensen <psp@dreamgruppen.dk>
# by Christian Sandholm Kastrup
# Latest update: 03-12-2019

# ===================================================================
#' @title Function for estimating a CES production function with the Kalman filter.
#'
#' @description This function first estimates the variances with the dlmMLE function from dlm package. Second, the
#' dlmSmooth function is applied to estimate the state parameters. Third, several statistical tests
#' are performed using the bgtest and bptest from lmtest package and the jb.norm.test from package normtest. Also, it calculates the normalized innovations squared.
#'
#' @details This function first estimates the variances with the dlmMLE function from dlm package. Second, the
#' dlmSmooth function is applied to estimate the state parameters. Third, several statistical tests
#' are performed using the bgtest and bptest from lmtest package and the jaruqe.bera.test from package normtest. Also, it calculates the normalized innovations squared.
#'
#' \strong{NOTE:} This function can be used on its own, but we encourage to use the CESKalman instead as it is more robust.
#'
#' The estimated function is the Error-Correction model from Kronborg et al (2019):
#' \eqn{\Delta s_{t}=\alpha(s_{t-1}-(1-\sigma)p_{t-1}-\mu_{t-1})+\sum_{i=0}^{nlags}\kappa_{i}\Delta p_{t-i}+\sum_{i=1}^{nlags}\gamma_{i}\Delta s_{t-i}+\epsilon_t}.
#' \eqn{s_t} is the relative budget shares in logs (expenditure on factor 1 relative to factor 2), \eqn{p_t} is the relative prices in logs and \eqn{\mu_t=(\sigma-1)log(\Gamma_t)} where \eqn{\Gamma_t} is the relative
#' augmenting technical change. The process of \eqn{\mu_t}, the state variable, is an I(2) process: \eqn{\Delta \mu_t=\Delta \mu_{t-1}+\eta_t}.
#'
#' X should be a matrix with all the explanatory variables to be used (all in logs). The dimension is TxK, where T
#' is time periods and K are the explanatory variables (not counting technical change as one). K=3 when nlags=0, =5 when nlags=1 and =7 when nlags=2. The column ordering should be:
#' 1: lags of relative budget shares in log-levels. 2: lags of the relative prices in log-levels. 3: First-difference of log relative prices.
#' 4: First-differences of log relative prices lagged one period. 5: First-differences of log relative budget shares lagged one period.
#' 6: First-differences of log relative prices lagged two periods. 7: First-differences of log relative budget shares lagged two periods.
#' When using the CESKalman, CESKalman_Estimation or CESKalman_Bootstrap functions, data will be ordered in this way automatically.
#'
#' grid.param_init can be a vector of any length. It is the values to loop over for the variance of the error term in observation equation, \eqn{\epsilon_t}, in the numerical optimization. When lambda is set to NA, different values
#' for the signal-to-noise ratio in the range from 10-1000 found to be optimal in Kronborg et al (2019) are tried. The one that maximizes the likelihood is chosen.
#' @param Y First-difference of log relative budget shares
#' @param X A matrix of all explanatory variables to be used (See details)
#' @param grid.param_init A string of length 3. First and second argument are the, respectively, lower and upper bound of the grid of the initial variance of the error term used in the dlmMLE function. Third is the step size
#' @param nlags Number of lags used of relative prices and relative budget shares. One of 0,1,2
#' @param lambda The inverse of the signal-to-noise ratio. Can be any positive value. If set to FALSE (the default) it is estimated
#' @param Leontief TRUE or FALSE (default). Should the elasticity of substitution be fixed at zero? If so, set =TRUE
#' @param sigma_init Initial value for the long run elasticity of substitution, \eqn{\sigma}
#' @param alpha_init Initial value for the adjustment parameter, \eqn{\alpha}
#' @param cVal_NIS Critical value used to construct the confidence bands of the NIS test. Default is 0.10
#'
#' @usage CESKalman_Estimation(Y,data,grid.param_init=c(-9,-1,1),nlags,lambda,
#' @usage          Leontief=FALSE,alpha_init,sigma_init,cVal_NIS=0.10)
#'
#' @return
#' Returns a list with the following objects:
#'
#' \strong{Smooth:} The list returned from the dlmSmooth function.
#' \strong{Gamma:} The relative augmenting technologies in logarithms (input 1 relative to input 2).
#' \strong{residuals:} The residuals from the observation equation.
#' \strong{sigma:} The long run elasticity of substitution, \eqn{\sigma}.
#' \strong{alpha:} The adjustment parameter, \eqn{\alpha}.
#' \strong{LV:} The likelihood value.
#' \strong{est.lambda:} The estimated value of lambda. (only estimated if lambda was specified as FALSE, else just the chosen value).
#' \strong{BG_test:} p-value from the Breusch Godfrey test for autocorrelation from function bgtest.
#' \strong{BP_test:} p-value from the Breusch Pagan test for heteroscedasticity from function bptest.
#' \strong{JB_test:} p-value from the Jarque Bera test for normality from function jbtest.
#' \strong{NIS_test:} Value, lower and upper confidence bands of the Normalized Innovations Squared test.
#' \strong{data:} The data applied in estimation.
#'
#' @seealso
#' See dlm package for a decription of dlm objects, in particular the dlmSmooth function. See Kronborg et al (2019) for a description of the methodology.
#'
#' @export


# ===================================================================
# Starting function
# ===================================================================

CESKalman_Estimation <- function(Y,X,grid.param_init=c(-9,-1,3),nlags,lambda,Leontief=FALSE,alpha_init=-0.3,sigma_init=0.5,
                                 cVal_NIS=0.10){

  library("dlm")
  library("lmtest")
  library("normtest")
  library("dplyr")

  ## Outline:
  # 1. Loading the needed packages
  # 2. Creating the matrix of explanatory variables
  # 3. Estimate the model with MLE for a grid of parameter values
  # 4. Use the Kalman smoother
  # 5. Perform a set of misspecification tests
  # 6. Return the results.

if(!between(alpha_init,-1,0)){
  warning("alpha_init should be in the interval -1 to 0 in order to secure error-correction")
}

# ########################### Ordering data #######################
#
#   S = log((data[,1]*data[,3])/(data[,2]*data[,4]))
#   P = log(data[,1]/data[,2])
#
# if(Leontief==FALSE){
#     X <- cbind(S[(1+nlags):(length(S)-1)],P[(1+nlags):(length(S)-1)],diff(P[(1+nlags):length(P)]))
#     }else{
#       X <- cbind(S[(1+nlags):(length(S)-1)]-P[(1+nlags):(length(S)-1)],diff(P[(1+nlags):length(P)]))
#     }
#
#   if(nlags==1){
#     X <- cbind(X,diff(P[1:(length(P)-1)]),diff(S[1:(length(S)-1)]))
#   }
#
#   if(nlags==2){
#     X <- cbind(X,diff(P[nlags:(length(P)-1)]),diff(S[nlags:(length(S)-1)]),
#                diff(P[1:(length(P)-nlags)]),diff(S[1:(length(S)-nlags)]))
#   }

  ##################### Estimates with MLE #########################


  ## Creating a vector of parameters to loop over

  set.param_init = c()

  param_init = grid.param_init[1] ## First is equal to the first value of the grid

  while(param_init<=grid.param_init[2]){
    set.param_init = c(set.param_init,param_init)
    param_init = param_init + grid.param_init[3]  ## add the step size

  }


  param = c()
  for(n in 1:length(set.param_init)){ # Variance on structural share parameter

    if(is.na(lambda)){
      param[1]=set.param_init[n]

      for(i in 1:5){
        tryCatch({


          param[2]=c(c(set.param_init[n]-log(100)),c(set.param_init[n]-log(10)),c(set.param_init[n]-log(1000)),c(set.param_init[n]-log(500)),c(set.param_init[n]-log(50)))[i]

          MLE       <- dlmMLE(Y, param, build_SS,X=X,nlags=nlags,lambda=lambda,Leontief=Leontief,alpha_init=alpha_init,sigma_init=sigma_init, method = "L-BFGS-B", debug = FALSE,hessian=TRUE)
          MLE.build <- build_SS(MLE$par,X=X,nlags=nlags,lambda=lambda,Leontief=Leontief,alpha_init=alpha_init,sigma_init=sigma_init)
          MLE.param     <- MLE$par

          if(n&i==1){
            MLE.est=MLE
          }else{
            if(MLE$value<MLE.est$value & MLE$convergence == 0){
              MLE.est=MLE

            }
          }


        }, error=function(e){cat("ERROR :",conditionMessage(e), " - New initial parameter values are tried", "\n")})

      }
    }else{
      tryCatch({

        param=set.param_init[n]
        MLE       <- dlmMLE(Y, param, build_SS,X=X,nlags=nlags,lambda=lambda,Leontief=Leontief,alpha_init=alpha_init,sigma_init=sigma_init, method = "L-BFGS-B", debug = FALSE,hessian=TRUE)
        MLE.build     <- build_SS(MLE$par,X=X,nlags=nlags,lambda=lambda,Leontief=Leontief,alpha_init=alpha_init,sigma_init=sigma_init)
        MLE.param     <- MLE$par

        if(n==1){
          MLE.est=MLE
        }else{
          if(MLE$value<MLE.est$value & MLE$convergence==0){
            MLE.est=MLE
          }
        }

        if(MLE.est$convergence == 0){
          #cat("MLE, convergence")
          # break
        }
      }, error=function(e){cat("ERROR :",conditionMessage(e)," - New initial parameter values are tried", "\n")})
    }

    # if(MLE.est$convergence == 0){
    # #  cat("MLE, convergence")
    #   break
    # }

  }

  if(MLE.est$convergence==!0){
    stop("No Convergence!")
    cat("MLE,No convergence")
  }

  LV             <- -MLE.est$value

  ########## Kalman smoothing ####################################

  Smooth         <- dlm::dlmSmooth(Y, MLE.build, debug=FALSE, simplify=FALSE)

  Smooth_Struc <- Smooth$s[-1,1]+Smooth$s[-1,2]+Smooth$s[-1,3]

  if(Leontief==FALSE){
    Smooth_Short <- Y - Smooth_Struc -Smooth$s[2,4]*X[,1]- rowSums(sapply(0:nlags, function(l) Smooth$s[2,5+2*l]*X[,2*l+2]+Smooth$s[2,6+2*l]*X[,2*l+3]))
  }else{
    Smooth_Short <- Y - Smooth_Struc- rowSums(sapply(0:nlags, function(l) Smooth$s[2,4+2*l]*X[,2*l+1]+Smooth$s[2,5+2*l]*X[,2*l+2]))
  }

  if(is.na(lambda)){
  return_lambda = exp(MLE.param[1])/exp(MLE.param[2])
  }else{
  return_lambda=lambda
  }

  if(Leontief==FALSE){
    Adjust          <- Smooth$s[1,4]
    LR_Elasticity   <- (Smooth$s[1,5])/Adjust+1
  }else{
    Adjust          <- Smooth$s[1,4]
    LR_Elasticity   <- 0
  }



  res        = Smooth_Short

  ############### Diagnostics ########################

  JB_test    = jb.norm.test(res)$p.value

  trend = 1:length(Smooth_Struc)

  Y_OLS = Y-Smooth_Struc

  ## OLS estimate
  OLS = stats::lm(Y_OLS~X)

  BG_test    = lmtest::bgtest(OLS,order=1,type=c("Chisq"))$p.value
  BP_test    = lmtest::bptest(OLS)$p.value

  ## NIS test
  NIS = matrix(NA,nrow=length(Smooth_Struc),ncol=1)

  for(t in 1:length(Smooth_Struc)){
    NIS[t] <- res[t]^2/(exp(MLE.param[1]))
  }
  NIS_test <- c(mean(NIS[1:length(Smooth_Struc)])*length(Smooth_Struc),qchisq(cVal_NIS/2,length(Smooth_Struc)),qchisq((1-cVal_NIS/2),length(Smooth_Struc)))


  Smooth_Struc    <- Smooth_Struc/-Adjust

  Gamma = (1/(LR_Elasticity-1))*Smooth_Struc

  output = list(Smooth,Gamma,res,LR_Elasticity,Adjust,LV,return_lambda,BG_test,JB_test,BP_test,NIS_test,data)
  names(output)=c("Smooth","Gamma","residuals","sigma","alpha","LV","est.lambda","BG_test","JB_test","BP_test","NIS_test","data")

  class(output) = c("list","CESKalman")

  return(output)

}
