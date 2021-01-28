# ===================================================================
# ------ PLOT FUNCTION FOR CESKalman OBJECT -------------------------
# ===================================================================

# ===================================================================
#' @author Christian Sandholm Kastrup <CST@dreamgruppen.dk>, Anders Farver Kronborg <ANK@dreamgruppen.dk> and Peter Philip Stephensen <PSP@dreamgruppen.dk>
# by Christian Sandholm Kastrup
# Latest update: 17-12-2020

# ===================================================================
#' @title Function for plotting the model fit and trend/price decomposition of a CESKalman object
#'
#' @description The first plot contains the model fit of the relative quantities. The second is a decomposition of the fitted values in a price and a trend component.
#' The third plot is the demeaned-data series applied in estimation.
#'
#' @details
#' \strong{NOTE:}
#' As we are primarily interested to explain the model fit of the relative quantities, the estimated equation is first translated into quantities, denoted \eqn{x_t}:
#' \eqn{\Delta x_{t}=\alpha(x_{t-1}+\sigma p_{t-1}-\mu_{t-1})+(\kappa_{i}-1)\Delta p_{t}+\epsilon_t}. Note that nlags=0 is currently only available option.
#'
#' This constitutes a differential equation. The fitted values can be decomposed into two terms, \eqn{\hat{x}_t=\hat{x}_t^{Trend}+\hat{x}_t^{Price}}.
#' The trend component is \eqn{\hat{x}_t^{Trend}\equiv \hat{x}_0(1+\hat{\alpha})^t-\hat{\alpha}\sum_{s=1}^{t}\hat{\mu}_{s-1} *(1+\hat{\alpha})^{t-s}}.
#' The price component is \eqn{\hat{x}_t^{Price}\equiv \sum_{s=1}^t(\hat{\alpha}\hat{\sigma}p_{s-1}+(\hat{\kappa}-1)\Delta p_s)(1+\hat{\alpha})^{t-s}}
#
#' For now, the function can only handle yearly data, if you want to specify the start and end date. If using data
#' at another frequency, you will have to use the default t0 and tEnd. In addition, the function is only suitable when nlags=0.
#'
#' @usage plot.CESKalman(Kalman,t0=1,tEnd=nrow(Kalman$data))
#'
#' @param Kalman an object of class CESKalman returned from the CESKalman function
#' @param t0 Start date. Only available for yearly data, else it is an index
#' @param tEnd End date. Only available for yearly data, else it is an index
#'
#'
#' @examples
#' ## First, data is loaded with the Load_Data function (or any other data set)
#' data = Load_Data(Country="USA",tstart=1970,tend=2017)
#'
#' data = data.frame(data)
#'
#' data = cbind(data$q,data$w,data$K,data$L)
#'
#' ## The second step is a call to the CESKalman function
#' Kalman = CESKalman(data,grid.lambda=c(10,1000,20),lambda_est_freely = F,max_nlags=0)
#'
#' ## Lastly we can plot the fit of the model:
#' plot(Kalman)
#'
#' @references Kronborg et al (2019)
#'
#' @export


# ===================================================================
# Starting function
# ===================================================================

plot.CESKalman <- function(Kalman,t0=1,tEnd=nrow(Kalman$data),main=""){
  par(mar=c(2,2,2,2))
#  par(mar=c(4,4,4,4))
  par(mfrow=c(3,1),oma = c(3, 0, 3, 0))

  nlags = Kalman$nlags  ## Correct for lags
  data = Kalman$data[(1+nlags):nrow(Kalman$data),]

  x = log((data[,3])/(data[,4]))
  x_0 = log((data[,3])/(data[,4]))[1]
  alpha_hat = Kalman$alpha
  p = log(data[,1]/data[,2])
  if(Kalman$sigma==0){phi_hat = c(NA,(Kalman$Smooth$s[1,5]-1)*diff(p))}else{phi_hat = c(NA,(Kalman$Smooth$s[1,6]-1)*diff(p))} ## We loose the first observation due to the difference
  beta_hat = -Kalman$sigma
  mu_hat = c((Kalman$sigma-1)*Kalman$Gamma,NA)  ## Note that since Gamma is lagged in the ECM, no observation is available at time T
  ## Fitted values

  x_hat = matrix(NA,nrow=length(x))
  x_hat_trend = matrix(NA,nrow=length(x))
  x_hat_price = matrix(NA,nrow=length(x))

    ## Initialization
  x_hat[1] = x_0

  ## Note that t=1 is time 0 and t=2 is time 1
  for(t in 2:length(x)){
    x_hat[t] = x_hat[t-1]+alpha_hat*(x_hat[t-1]-beta_hat*p[t-1]-mu_hat[t-1])+phi_hat[t]
    x_hat_trend[t] = x_hat[1]*(1+alpha_hat)^(t-1)+sum(-alpha_hat*mu_hat[(t-1):1]*(1+alpha_hat)^(0:(t-2)) )
    x_hat_price[t] = sum((-alpha_hat*beta_hat*p[(t-1):1]+phi_hat[t:2])*(1+alpha_hat)^(0:(t-2)) )
  }
  cor=round(cor(x_hat_trend[-1],x_hat_price[-1]),2)
  x_hat = x_hat_trend+x_hat_price

  x_hat   = ts(x_hat[-1],start=t0+1+nlags,end=tEnd,frequency = 1)
  x_hat_trend = ts(x_hat_trend[-1],start=t0+1+nlags,end=tEnd,frequency = 1)
  x_hat_price = ts(x_hat_price[-1],start=t0+1+nlags,end=tEnd,frequency = 1)
  x      = ts(x[-1],start=t0+1+nlags,end=tEnd,frequency = 1)
  p      = ts(p[-1],start=t0+1+nlags,end=tEnd,frequency=1)

#
#   Constant = mean(s-s_hat_price-s_hat_trend,na.rm = T)
#   plot(s-s_hat_price,lty=1,type="l",ylim=c(min(s-s_hat_price,s_hat_trend+Constant,na.rm = T),max(s-s_hat_price,s_hat_trend+Constant,na.rm = T)),ylab="",xlab="",main="Trend and cycle",lwd=2)
#   lines(s_hat_trend+Constant,lty=2,lwd=2)
#
#
#   par(new = TRUE)
#   plot(s-s_hat_price-s_hat_trend-Constant,type="h",axes=F,ylab="",xlab="",lty=1)
#   axis(4)
#   abline(h=0,col="grey40")
#
#   legend("topleft",c("Trend+cycle","Trend","Cycle rhs."),lty=c(1,2,1),cex=1,lwd=c(2,2,1))
#
#
#   Constant = mean(s-s_hat_price-s_hat_trend,na.rm = T)
#   plot(s-s_hat_trend,lty=1,type="l",ylim=c(min(s-s_hat_trend,s_hat_price+Constant,na.rm = T),max(s-s_hat_trend,s_hat_trend+Constant,na.rm = T)),ylab="",xlab="",main="Prices and cycle",lwd=2)
#   lines(s_hat_price+Constant,lty=2,lwd=2)
#
#
#   par(new = TRUE)
#   plot(s-s_hat_price-s_hat_trend-Constant,type="h",lty=1,axes=F,ylab="",xlab="")
#   axis(4)
#   abline(h=0,col="grey40")
#
#   legend("topleft",c("Prices+cycle","Prices","Cycle rhs."),lty=c(1,2,1),cex=1,lwd=c(2,2,1))
  dat = cbind(x,x_hat)

    plot(x,lty=1,type="l",ylim=c(min(dat),max(dat,na.rm = T)),ylab="",xlab="",main="Model fit of relative quantities",lwd=2)
  lines(x_hat,lty=1,lwd=2,col="red")


  par(new = TRUE)
  plot(x-x_hat,type="h",axes=F,ylab="",xlab="",lty=1)
  axis(4)
  abline(h=0,col="grey40")

  legend("topleft",c("Relative quantities","Fit","Residuals rhs."),lty=c(1,1,1),,col=c("black","red","black"),cex=1,lwd=c(2,2,1))

dat = cbind(x_hat,x_hat_trend,x_hat_price)

    plot(x_hat,lty=1,type="l",ylim=c(min(dat,na.rm = T),max(dat,na.rm = T)),ylab="",xlab="",main="Decomposition of fitted relative quantities",lwd=2)
    lines(x_hat_trend,lty=1,lwd=2,col="red")
    lines(x_hat_price,lty=1,lwd=2,col="blue")

    legend("topleft",c("Fit","Trend","Price-effects"),lty=c(1,1,1),col=c("black","red","blue"),cex=1,lwd=c(2,2,2))

    dat = cbind(x-mean(x),p-mean(p))

    plot(x-mean(x),lty=1,type="l",ylab="",xlab="",main="Demeaned data",lwd=2,ylim=c(min(dat,na.rm = T),max(dat,na.rm = T)))

    par(new = TRUE)
    plot(p-mean(p),type="l",ylim=rev(range(dat)),axes=F,ylab="",xlab="",lty=1,lwd=2,col="red")
    axis(4)
   # abline(h=0,col="grey40")

    # lines(-p+mean(p),type="l",lwd=2,col="red")

    legend("topleft",c("Relative quantities","Relative prices rhs"),lty=1,lwd=2,col=c("black","red"))


    lambda = round(Kalman$est.lambda)
    sigma= round(Kalman$sigma,2)
    alpha = round(Kalman$alpha,2)

    mtext(main,side=3,outer=TRUE,cex=1.5,line=0.5)
  #  mtext(paste("nlags=",nlags,sep=""),side=4,outer=TRUE,cex=1.5,line=0.5)
  #  mtext(paste("correlation=",cor,sep=""),side=2,outer=TRUE,cex=1.5,line=0.5)
    mtext(bquote(lambda == .(lambda) ~~~ sigma == .(sigma) ~~~ alpha == .(alpha)), outer=TRUE, side=1, cex=1.5, line=0.5)
    #      mtext(expression(paste(lambda,"=",round(Kalman$est.lambda)," ",sigma,"=",round(Kalman$sigma,2)," ", alpha,"=",round(Kalman$alpha,2),sep=""), outer=TRUE, side=1, cex=1.5, line=0.5)

}
