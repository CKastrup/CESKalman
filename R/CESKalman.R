# ===================================================================
# ------ FUNCTION FOR ESTIMATING A CES PRODUCTION FUNCTION ----------
# ===================================================================

# ===================================================================
#' @author Christian Sandholm Kastrup <CST@dreamgruppen.dk> , Anders Farver Kronborg <ANK@dreamgruppen.dk> and Peter Philip Stephensen <psp@dreamgruppen.dk>
# by Christian Sandholm Kastrup
# Latest update: 17-12-2020

# ===================================================================
#' @title Function for estimating a CES production function with the Kalman filter.
#'
#' @description The function applies the CESKalman_Estimation function to estimate the elasticity of substitution in a
#' CES production function using the Kalman filter. CESKalman is more robust than the CESKalman_Estimation as it loops over different values of
#' the signal-to-noise ratio, the number of lags, and initial parameter values. The combination that maximizes the likelihood
#' given the model is well specified based on the bgtest from lmtest package and the Normalized Innovations Squared is chosen.
#'
#' @details This function calls the function CESKalman_Estimation multiple times as it loops over a grid of values
#' of the signal-to-noise ratio, number of lags chosen to prevent autocorrelation and different combinations of initial parameter values of sigma and alpha.
#' The optimal signal-to-noise ratio and initial parameter values are chosen to maximize the likelihood, conditional on
#' the model being well specified based on a Breusch Godfrey test for autocorrelation and a Normialized Innovations Squared test for filter misspecification.
#'
#' The estimated function is the Error-Correction model from Kronborg et al (2019):
#' \eqn{\Delta s_{t}=\alpha(s_{t-1}-(1-\sigma)p_{t-1}-\mu_{t-1})+\sum_{i=0}^{nlags}\kappa_{i}\Delta p_{t-i}+\sum_{i=1}^{nlags}\gamma_{i}\Delta s_{t-i}+\epsilon_t}.
#' \eqn{s_t} is the relative budget shares in logs (expenditure on factor 1 relative to factor 2), \eqn{p_t} is the relative prices in logs and \eqn{\mu_t=(\sigma-1)log(\Gamma_t)} where \eqn{\Gamma_t} is the relative
#' augmenting technical change. The process of \eqn{\mu_t}, the state variable, is an I(2) process: \eqn{\Delta \mu_t=\Delta \mu_{t-1}+\eta_t}.
#'
#' data should be a matrix or ts matrix containing the data series and the dimension Tx4. First and second column is the price of the first and second factor, respectively. Third and forth is the quantity of factor 1 and factor 2, respectively.
#'
#' The function loops over a grid of different signal-to-noise ratios in the range specified in grid.lambda. In addition, it is possible to also estimate this parameter freely by setting lambda_est_freely=TRUE.
#' In this case, the preferred value from the grid search is compared to the free estimation and the one that maximizes likelihood given the estimation being well specified is chosen.
#'
#' If sigma is estimated to be negative it is restricted to zero, but we still allow for short run fluctuations of prices to influence the budget shares.
#'
#' Potential convergence errors of dlmMLE can most often be resolved by increasing the range of grid.param_init, e.g. to grid.param_init=c(-9,-1,1).
#' Note that the variances are specified such that \eqn{\Sigma^\epsilon=exp(param[1])} and \eqn{\Sigma^\eta=exp(param[2])} if lambda is freely estimated.
#' Else, only exp(param[1]) is used and \eqn{\lambda\Sigma^\eta=\Sigma^\epsilon}. When lambda is freely estimated the grid.param_init is still only over values of
#' exp(param[1]), i.e. the observation variance. But for any observation variance we try five different values of lambda in the interval 10-1000 as initial value of exp(param[2]).
#'
#' @param data A matrix or ts matrix with the data series (see details)
#' @param grid.param_init A string of length 3. First and second argument are the, respectively, lower and upper bound of the grid of the initial variance of the error term used in the dlmMLE function. Third is the step size
#' @param max_nlags Maximum number of lags used. One of 0,1,2
#' @param grid.lambda A string of length 3. First and second argument are the, respectively, lower and upper bound of the grid of the inverse of the signal-to-noise ratio. Third is the step size. If set to NA, no grid is used
#' @param grid.sigma_init A string of length 3. First and second argument are the lower and upper bound of initial values for the long run elasticity of substitution, \eqn{\sigma}, to loop over. Third is the step size
#' @param grid.alpha_init A string of length 3. First and second argument are the lower and upper bound of initial values for the adjustment parameter, \eqn{\alpha}, to loop over. Third is the step size
#' @param cVal_NIS Critical value to be used to construct the confidence bands of the NIS test
#' @param cVal_Auto Critical value to be used in the Breusch Godfrey test for autocorrelation
#' @param print_results Do you want the function to print results while estimating?
#' @param lambda_est_freely Should lambda also be estimated freely?
#'
#'
#' @usage CESKalman(data,grid.param_init=c(-9,-1,5),max_nlags=2,grid.lambda,
#'                 grid.alpha_init=c(-0.9,-0.1,0.3),grid.sigma_init=c(0,1.5,0.5),
#'                 cVal_NIS=0.10,cVal_Auto=0.10,print_results=TRUE,lambda_est_freely=TRUE)
#
#'
#'
#' @return
#' Returns a list of class CESKalman from the preferred call of function CESKalman_Estimation (highest likelihood and well specified). See CESKalman_Estimation for description. In addition, the data series applied is returned.
#'
#' @examples
#'
#' ## First, data is loaded with the Load_Data function (or any other data set)
#' data = Load_Data(Country="USA",tstart=1970,tend=2017)
#'
#' data = data.frame(data)
#'
#' data = cbind(data$q,data$w,data$K,data$L)
#'
#' ## We can then estimate with four different approaches that all depends on the value of lambda.
#'
#' # 1. Grid combined with free estimation (our preffered method)
#' Kalman = CESKalman(data=data,grid.lambda=c(10,500,20), lambda_est_freely=TRUE)
#'
#' # 2. Grid only
#' Kalman = CESKalman(data=data,grid.lambda=c(10,500,20), lambda_est_freely=FALSE)
#'
#' # 3. Free estimation only
#' Kalman = CESKalman(data=data,grid.lambda=NA, lambda_est_freely=TRUE)
#'
#' # 4. A fixed value of lambda, e.g. lambda=100
#' Kalman = CESKalman(data=data,grid.lambda=c(100,100,100), lambda_est_freely=FALSE) ## Only argument [1] used
#'
#' ## Output is displayed with e.g. the plot function:
#' plot(Kalman)
#'
#' sigma=Kalman$sigma # This is the elasticity
#' alpha=Kalman$alpha # This is the adjustment parameter
#' Gamma=Kalman$Gamma # This is the relative log augmenting technologies
#'
#'
#' @references Kronborg et al (2019)
#'
#' @export






# ===================================================================
# Starting function
# ===================================================================

CESKalman <- function(data,grid.param_init=c(-9,-1,5),max_nlags=2,grid.lambda,grid.alpha_init=c(-0.9,-0.1,0.3),grid.sigma_init=c(0,1.5,0.5),
                                 cVal_NIS=0.10,cVal_Auto=0.10,print_results=TRUE,lambda_est_freely=TRUE){


  ## Outline:
  # The function is used to determine the optimal lambda, number of lags and initial value of alpha and sigma.
  # The combination of alpha_init and sigma_init that maximizes the likelihood is chosen
  # Lags is included if autocorrelation is found
  # The lambda that is well specified and has the highest likelihood value of the filter is chosen

if(is.na(grid.lambda[1])&lambda_est_freely==FALSE){
  stop("grid.lambda must be non-NA and/or lambda_est freely must be TRUE")
}

  ## Loading the needed packages
  library("dlm")
  library("lmtest")
  library("normtest")

  ################# Creates a set of vectors with the values to loop over ################

  ## Lambda
  if(lambda_est_freely){
  set.lambda=NA   ## Always start with a free estimation
  }else{
  set.lambda=c()
}

  if(is.na(grid.lambda[1])==FALSE){   ## Creates a vector of lambdas to be looped over, if a grid is specified
    lambda=grid.lambda[1]
    while(lambda<=grid.lambda[2]){
  set.lambda = c(set.lambda,lambda)
  lambda = lambda + grid.lambda[3]
    }
  }


  ## Alpha

  alpha_init    =  grid.alpha_init[1] ## First equal to the lower bound
  set.alpha_init = c()
  while(alpha_init<=grid.alpha_init[2]){   ## Creates a vector of alphas to be looped over
   set.alpha_init = c(set.alpha_init,alpha_init)
   alpha_init = alpha_init+grid.alpha_init[3]
  }

 ## sigma
  sigma_init    =  grid.sigma_init[1] ## First equal to the lower bound
  set.sigma_init = c()
  while(sigma_init<=grid.sigma_init[2]){   ## Creates a vector of alphas to be looped over
    set.sigma_init = c(set.sigma_init,sigma_init)
    sigma_init = sigma_init+grid.sigma_init[3]
  }


  ########################### Ordering data #######################

  S = log((data[,1]*data[,3])/(data[,2]*data[,4]))
  P = log(data[,1]/data[,2])
  
## Start looping and choosing the optimal one



for(lambda in set.lambda){  ## Outer loop is over different signal-to-noise ratios, lambda
for(Leontief in c(FALSE,TRUE)){   ## If the elasticity is found to be negative, restrict it to zero, else break.

  for(nlags in 0:max_nlags){  ##  Include lags until no autocorrelation is present.


    if(Leontief==FALSE){
      X <- cbind(S[(1+nlags):(length(S)-1)],P[(1+nlags):(length(S)-1)],diff(P[(1+nlags):length(P)]))
    }else{
      X <- cbind(S[(1+nlags):(length(S)-1)]-P[(1+nlags):(length(S)-1)],diff(P[(1+nlags):length(P)]))
    }

    if(nlags==1){
      X <- cbind(X,diff(P[1:(length(P)-1)]),diff(S[1:(length(S)-1)]))
    }

    if(nlags==2){
      X <- cbind(X,diff(P[nlags:(length(P)-1)]),diff(S[nlags:(length(S)-1)]),
                 diff(P[1:(length(P)-nlags)]),diff(S[1:(length(S)-nlags)]))
    }

      Y = diff(S[(1+nlags):length(S)])   # Correct the length of Y if more lags are included

    for(alpha_init in set.alpha_init){   ## Loop over all values in the grid of alpha
      for(sigma_init in set.sigma_init){ ## Loop over all values in the grid of sigma. Choose the one that maximizes the likelihood

 ## tmp1 is the newest call. The preferred one is called tmp
        tmp1 = CESKalman_Estimation(Y=Y,X=X,grid.param_init=grid.param_init,nlags=nlags,lambda=lambda,Leontief=Leontief,alpha_init=alpha_init,sigma_init=sigma_init,
                                      cVal_NIS=cVal_NIS)

        tmp1$sigma_init = sigma_init
        tmp1$alpha_init = alpha_init
        tmp1$nlags = nlags

        if(sigma_init==set.sigma_init[1]&alpha_init==set.alpha_init[1]){ ## If it is the first one, set it equal to the preferred one

            tmp <- tmp1

        } else{

          ## ...If not, test if the preferred elasticity is negative and the new estimate is higher. Then update

          if(tmp$sigma<0){
            if(tmp$sigma<tmp1$sigma){
              tmp = tmp1
            }
          }else{

            ##.. Else, i.e. the elasticity is positive, choose the one with the highest likelihood value.

            if(tmp$LV < tmp1$LV){

              tmp <- tmp1

            }
          }
        }


      } ## Stop looping over sigma

    } ## Stop looping over alpha



    if(tmp$BG_test>cVal_Auto){
      break  ## Breaks the loop over lags, i.e. no more lags are included if no autocorrelation is found
    }
  }   ## Stop the loop over nlags


  if(tmp$sigma>0){
    break   ## Breaks the loop over the elasticity, i.e. if estimated to be positive we continue, else, it is restricted to zero.
  }
  } ##Stops the loop of Leontief

  ## Print some results
  if(print_results==TRUE){
  cat("     lambda=",round(tmp$est.lambda,2),"  Likelihood=",round(tmp$LV,2),"  sigma=",round(tmp$sigma,2),"  NIS=",(tmp$NIS_test[1]>tmp$NIS_test[2]&tmp$NIS_test[1]<tmp$NIS_test[3]),"  BG_test=",(tmp$BG_test>cVal_Auto),"\n")
  }


  ## Choosing the best lambda
if(lambda_est_freely){
  if(is.na(lambda)){  ## If first one, it is always the preferred one
    tmp_pref = tmp
  }else{
    if(all(tmp_pref$NIS_test[1]<tmp_pref$NIS_test[3] & tmp_pref$NIS_test[1]>tmp_pref$NIS_test[2] & tmp_pref$BG_test>cVal_Auto)==TRUE){  ## Preffered is within the NIS confidence bands and the BG test for autocorrelation
      if(all(tmp$NIS_test[1]<tmp$NIS_test[3] & tmp$NIS_test[1]>tmp$NIS_test[2] & tmp$BG_test>cVal_Auto &tmp$LV>tmp_pref$LV)==TRUE){  ## Last try is within the confidence bands and likelihood is higher
        tmp_pref = tmp  ## Then replace!

      }
    }else{   ## If it is misspecified
      if(all(tmp$NIS_test[1]<tmp$NIS_test[3] & tmp$NIS_test[1]>tmp$NIS_test[2] & tmp$BG_test>cVal_Auto)==TRUE){
        tmp_pref = tmp
      }
    }

  }
}else{
  if(lambda==set.lambda[1]){  ## If first one, it is always the preferred one
    tmp_pref = tmp
  }else{
    if(all(tmp_pref$NIS_test[1]<tmp_pref$NIS_test[3] & tmp_pref$NIS_test[1]>tmp_pref$NIS_test[2] & tmp_pref$BG_test>cVal_Auto)==TRUE){  ## Preffered is within the NIS confidence bands and the BG test for autocorrelation
      if(all(tmp$NIS_test[1]<tmp$NIS_test[3] & tmp$NIS_test[1]>tmp$NIS_test[2] & tmp$BG_test>cVal_Auto &tmp$LV>tmp_pref$LV)==TRUE){  ## Last try is within the confidence bands and likelihood is higher
        tmp_pref = tmp  ## Then replace!

      }
    }else{   ## If preferred is misspecified
   #   if(all(tmp$NIS_test[1]<tmp$NIS_test[3] & tmp$NIS_test[1]>tmp$NIS_test[2] & tmp$BG_test>cVal_Auto)==TRUE){
      if(tmp$LV>tmp_pref$LV){    # Replace if the likelihood is higher
        tmp_pref = tmp
      }
    }

  }
}


} ## Stops looping over lambda

  tmp_pref$data = data
  
  return(tmp_pref) ## Return the preferred one.

}




