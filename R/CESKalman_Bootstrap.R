
# ===================================================================
#' @author Christian Sandholm Kastrup <CST@dreamgruppen.dk>, Anders Farver Kronborg <ANK@dreamgruppen.dk> and Peter Philip Stephensen <PSP@dreamgruppen.dk>
# by Christian Sandholm Kastrup
# Latest update: 03-12-2019

# ===================================================================
#' @title Function for bootstrapping confidence intervals of the elasticity and adjustment parameter
#'
#' @description The function uses a recursive-design bootstrapping procedure to produce confidence intervals of
#' the elasticity of substitution and adjustment parameter. The CESKalman_Estimation function is applied and
#' estimated for every new draw.
#'
#' @details This function takes an object of class CESKalman as input and calls the function CESKalman_Estimation for every new draw. The parameter estimates
#' from CESKalman are used as initial values of sigma and alpha. Also, the number of lags and lambda are determined
#' by the CESKalman function.
#'
#' grid.param_init can be a vector of any length. Values to loop over for the variance of error term in observation equation, \eqn{\epsilon_t}, in the numerical optimization. When lambda is FALSE, differen values
#' for the signal-to-noise ratio in the range from 10-1000 found to be optimal in Kronborg et al (2019) are tried. The one that maximizes the likelihood is chosen.
#'
#' Only draws with no autocorrelation in the estimated model and the NIS test within the limits specified in cVal_NIS is accepted.
#'
#' @param Estimation An object of class CESKalman
#' @param grid.param_init A vector of initial Parameter values of the observation variance to loop over (see details)
#' @param ndraw The number of draws in the bootstrapping procedure, 1000 is default.
#' @param print_results Should results be printed while bootstrapping?
#' @param cVal_auto Critical value for the Breusch-Godfrey test for autocorrelation
#' @param cVal_NIS Critical value for the NIS test 
#' 
#' @usage CESKalman_Bootstrap(Estimation,grid.param_init=c(-9,-1,3),ndraw=1000,print_results=TRUE,cVal_Auto=0.1,cVal_NIS=0.1)
#
#'
#'
#' @return
#' Returns a matrix of dimension ndrawX3. First column is the draws of sigma, second is the draws of alpha and third is the likelihood value.
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
#' ## Next, the CESKalman function is called
#' Kalman = CESKalman(data=data,grid.lambda=c(10,500,20), lambda_est_freely=TRUE)
#'
#' ## Bootstrapping confidence intervals
#' Bootstrap = CESKalman_Bootstrap(Estimation=Kalman)
#'
#' @references Kronborg et al (2019) and Petris et al (2010)
#'
#' @export






# ===================================================================
# Starting function
# ===================================================================

CESKalman_Bootstrap <- function(Estimation,grid.param_init=c(-9,-1,5),ndraw=1000,cVal_NIS=0.10,cVal_Auto=0.10,print_results=TRUE){


  ## Start bootstrapping

  boot_alpha = matrix(NA,ncol=1,nrow=ndraw)
  boot_sigma  = matrix(NA,ncol=1,nrow=ndraw)
  boot_likelihood = matrix(NA,ncol=1,nrow=ndraw)


  nlags = Estimation$nlags

  res      = Estimation$residuals

  alpha = Estimation$alpha

  sigma = Estimation$sigma

  Gamma = c(rep(NA,nlags),Estimation$Gamma,NA)

  data = Estimation$data

  if(sigma==0){kappa=Estimation$Smooth$s[1,5]}else{kappa=Estimation$Smooth$s[1,6]}
  if(sigma==0&nlags>0){kappa1=Estimation$Smooth$s[1,6]}
  if(sigma>0&nlags>0){kappa1=Estimation$Smooth$s[1,7]}
  if(sigma==0&nlags>1){kappa2=Estimation$Smooth$s[1,8]}
  if(sigma>0&nlags>1){kappa2=Estimation$Smooth$s[1,9]}

  if(sigma==0&nlags>0){gamma1=Estimation$Smooth$s[1,7]}
  if(sigma>0&nlags>0){gamma1=Estimation$Smooth$s[1,8]}
  if(sigma==0&nlags>1){gamma2=Estimation$Smooth$s[1,9]}
  if(sigma>0&nlags>1){gamma2=Estimation$Smooth$s[1,10]}

  if(sigma==0){Leontief=TRUE}else{Leontief=FALSE}

  ########################### Ordering data #######################

  S = log((data[,1]*data[,3])/(data[,2]*data[,4]))

  P = log(data[,1]/data[,2])



  DP = c(0,diff(P))

  nacc <- 1
  k    <- 1
  while(k <= ndraw){


    ## Draw residuals
    draw_res = c(rep(0,1+nlags),sample(res,replace=TRUE))

    S_boot = matrix(NA,ncol=1,nrow=length(S))

    S_boot[1:(1+nlags)] = S[1:(1+nlags)]


    ## Recursive design
    for(t in (2+nlags):length(S)){
      if(nlags==0){phi = kappa*(P[t]-P[t-1])}

      if(nlags==1){phi = kappa*(P[t]-P[t-1])+kappa1*(P[t-1]-P[t-2])+gamma1*(S_boot[t-1]-S_boot[t-2])}

      if(nlags==2){phi = kappa*(P[t]-P[t-1])+kappa1*(P[t-1]-P[t-2])+kappa2*(P[t-2]-P[t-3])+gamma1*(S_boot[t-1]-S_boot[t-2])+gamma2*(S_boot[t-2]-S_boot[t-3])}

     S_boot[t] = S_boot[t-1]+alpha*(S_boot[t-1]-(1-sigma)*P[t-1]-(sigma-1)*Gamma[t-1])+phi+draw_res[t]

    }

    ## Creating the resulting data to be used in estimation

   Y_boot = diff(S_boot[(1+nlags):length(S_boot)])   # Correct the length of Y if more lags are included


    if(Leontief==FALSE){
        X_boot <- cbind(S_boot[(1+nlags):(length(S)-1)],P[(1+nlags):(length(S)-1)],diff(P[(1+nlags):length(P)]))
        }else{
          X_boot <- cbind(S_boot[(1+nlags):(length(S)-1)]-P[(1+nlags):(length(S)-1)],diff(P[(1+nlags):length(P)]))
        }

      if(nlags==1){
        X_boot <- cbind(X_boot,diff(P[1:(length(P)-1)]),diff(S_boot[1:(length(S)-1)]))
      }

      if(nlags==2){
        X_boot <- cbind(X_boot,diff(P[nlags:(length(P)-1)]),diff(S_boot[nlags:(length(S)-1)]),
                   diff(P[1:(length(P)-nlags)]),diff(S_boot[1:(length(S)-nlags)]))
      }


    tmp = CESKalman_Estimation(Y=Y_boot,X=X_boot,grid.param_init=grid.param_init,nlags=nlags,lambda=Estimation$est.lambda,Leontief=FALSE,alpha_init=Estimation$alpha,sigma_init=Estimation$sigma,
                                cVal_NIS=cVal_NIS)
if(tmp$BG_test>cVal_Auto&tmp$NIS_test[1]>tmp$NIS_test[2]&tmp$NIS_test[1]<tmp$NIS_test[3]){
      boot_alpha[k] = tmp$alpha
      boot_sigma[k]  = tmp$sigma
      boot_likelihood[k] = tmp$LV
      k <- k + 1


}
    
    nacc = nacc+1
    nacc_ratio= k/nacc

if(print_results){
    if(k %% 10 == 0){
      cat("\014")
      cat(paste("This is draw number ",k,sep=""), "\n")
      cat(paste("Acceptance ratio ",round(nacc_ratio,digits=2)), "\n")
    }
}
  }
if(print_results){
  cat("\014")
  cat("Bootstrapping done!", "\n")
  cat("\n")
  cat("Elasticity: ", "\n", "\n")
  cat("Estimate           : ", Estimation$sigma, "\n")
  cat("Bootstrap median   : ", median(boot_sigma), "\n")
  cat("Bootstrap S.d.     : ", sd(boot_sigma), "\n")
  cat("5 pct. quantile    : ", quantile(boot_sigma, probs = 0.05), "\n")
  cat("95 pct. quantile   : ", quantile(boot_sigma, probs = 0.95), "\n", "\n")

  cat("Adjustment: ", "\n", "\n")
  cat("Estimate           : ", Estimation$alpha, "\n")
  cat("Bootstrap median   : ", median(boot_alpha), "\n")
  cat("Bootstrap S.d.     : ", sd(boot_alpha), "\n")
  cat("5 pct. quantile    : ", quantile(boot_alpha, probs = 0.05), "\n")
  cat("95 pct. quantile   : ", quantile(boot_alpha, probs = 0.95), "\n", "\n")


  par(mfrow=c(1,2))
  hist(boot_sigma, main="Elasticity",breaks=20)
  abline(v=Estimation$sigma, lwd=2, lty=1)
  abline(v=quantile(boot_sigma, probs = 0.05), lwd=2, lty=2)
  abline(v=quantile(boot_sigma, probs = 0.95), lwd=2, lty=2)

  hist(boot_alpha, main="Adjustment coefficient",breaks=15)
  abline(v=Estimation$alpha, lwd=2, lty=1)
  abline(v=quantile(boot_alpha, probs = 0.05), lwd=2, lty=2)
  abline(v=quantile(boot_alpha, probs = 0.95), lwd=2, lty=2)
}

  boot = cbind(boot_sigma,boot_alpha,boot_likelihood)
  colnames(boot) = c("sigma","alpha","likelihood")

  return(boot)

}




