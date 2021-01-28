---
title: "'Estimating the Constant Elasticity of Substitution Subject to Time-varying,"
  Nonlinear, Technical Change: The CESKalman R-package'
author: "Christian S. Kastrup, Anders F. Kronborg and Peter P. Stephensen"
date: "\today{}"
output:
  html_document:
    df_print: paged
linestretch: 1.5
bibliography: Library_Seminar.bib
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 1. Introduction
The elasticity of substitution between capital and labor is a crucial parameter in economic models. It affects many outcomes such as the distribution of income, effect of labor market reforms and industry allocations of production factors. An issue when estimating this parameter is how to specify the process of technical change: If misspecified, it is known to potentially bias the estimate of the elasticity. In this note, we show how this elasticity can be estimated subject to a flexible data-driven process of technical change using the **CESKalman** R-package.

The methodology applied in the package relies on @Kronborg2019 who show that a Kalman filter can be applied to estimate the elasticity of substitution. In the paper, we show that the method is superior to the common linear trend assumption and nests several different types of technical change, such as the commonly used Box-Cox transformation. The process of technical change is specified in a flexible data-driven way as an I(2) process. This allows for so-called "medium run" fluctuations, that is, several periods with shifts in the direction of technical change. The value of the inverse of the signal-to-noise ratio (the ratio of measurement error variance to variance of technical change error, $\lambda$) determines the degree of smoothness of the process of technical change: Low values of $\lambda$ subscribes all short run fluctuations of relative expenditure shares not described by the relative prices to technical change. Oppositely, $\lambda\rightarrow\infty$ corresponds to the commonly used linear trend assumption (e.g. @Antras2004 and @Leon-Ledesma2015). Section 2 describes in short the methodology. The interested reader is referred to @Kronborg2019 for further information and motivation. 

In Section 3 we describe how the package can be used to estimate the elasticity of substitution in any CES function with two input factors. The main function is **CESKalman** and returns a list of class CESKalman containing the estimated value of the elasticity along with a process of relative technical change. The level of smoothness is determined by a maximum likelihood estimation combined with a grid search procedure. The preferred value of the ratio is the one that maximizes the likelihood conditional on being well specified based on a Breusch-Godfrey test for autocorrelation and a Normalized Innovations Squarred test (NIS). We motivate the role for medium run fluctuations based on US data from the PWT database (@Feenstra2015) available in the package through the function **Load_Data**. We argue that these medium run fluctuations are persistent and that failure of allowing for such shifts in the direction of technical change may bias the estimate of the elasticity of substitution. Next, we show how the **CESKalman** function can be applied to estimate the elasticity while allowing for medium run fluctuations. We estimate the US elasticity to 0.57 in our preferred specification. Technical change is labor-augmenting in the long run, but capital-augmenting in the 90'es possibly reflecting the IT-revolution, consistent with @Klump2008 for the Euro area. In addition, we show how the output can be displayed with the **plot.CESKalman** function and how to bootstrap standard errors with **CESKalman_Bootstrap**. Lastly, we analyse how the value of $\lambda$ affects the elasticity and illustrate how the smoothness of technical change is determined by $\lambda$.

# 2. Methodology
The CES production function consisting of capital and labor as input is given by:
$$Y_t=[(\Gamma_t^KK_t)^{\frac{\sigma-1}{\sigma}}+(\Gamma_t^LL_t)^{\frac{\sigma-1}{\sigma}}]^{\frac{\sigma}{\sigma-1}}$$
$Y_t$ denotes total output, $K_t$ is capital and $L_t$ is labor. $\Gamma_t^K$ and $\Gamma_t^L$ are, respectively, capital and labor augmenting technical change. We denote the ratio of these as the relative augmenting technical change, $\Gamma_t\equiv\Gamma_t^K/\Gamma_t^L$. If this ratio is increasing it implies that technical change is relatively augmenting capital, and labor if decreasing. $\sigma$ is the elasticity of substitution, the parameter of primary interest. Assuming profit-maximizing firms the relative expenditure shares for the two factors is found from the cost minimization problem:

$$log(\frac{q_tK_t}{w_tL_t})=(\sigma-1)log(\Gamma_t)+(1-\sigma)log(\frac{q_t}{w_t})\quad(1)$$

$q_t$ and $w_t$ denotes, respectively, the user cost of capital and the wage. From the equation three different cases arise depending on the value of $\sigma$: If $\sigma=1$, the production function is Cobb-Douglas and technical change nor the relative prices change the relative expenditure share. If $\sigma<1$, capital and labor are gross complements. In this case, labor-augmenting technical change will imply a *decreasing* share of labor. Oppositely, an increase in the wage will imply an *increase* in the share of labor. If $\sigma>1$, they are gross substutites and the opposite holds true. This highlights the well-known result that the direction of technical change depends crucially on the elasticity of substitution. Consequently, they have to be estimated jointly. 

## 2.1 Adjustment costs
Equation (1) above is a static equation as it does not contain any dynamics and is therefore the long run equilibrium. However, in macroeconomic modelling, adjustment costs are often imposed to ensure that there are lags in the response of quantities to changes in the relative price (e.g. @Christiano2005 and @Smets2007). For this reason, the long run equation above may be subject to a small sample bias due to the relatively short sample in most time series. Therefore, we apply an error-correction model to allow for short run dynamics. We estimate the equation:
$$\Delta s_t=\alpha(s_{t-1}-(1-\sigma)p_{t-1}-\mu_{t-1})+\sum_{i=0}^{k}\kappa_i\Delta p_{t-i}+\sum_{i=1}^{k}\gamma_i\Delta s_{t-i}+\varepsilon_t,\quad \varepsilon_t = N(0,\Sigma^\varepsilon)\quad(2)$$
$s_t\equiv log(\frac{q_tK_t}{w_tL_t})$ is the relative expenditure share, $p_t\equiv log(\frac{q_t}{w_t})$ is the relative prices and $\mu_t\equiv (\sigma-1)log(\Gamma_t)$, the state parameter in the Kalman filter reflecting technical change.[^note1] The speed of adjustment to the long run is $\alpha$ and $\kappa$ and $\gamma$ are short run elasticities with respect to prices and expenditure shares, respectively.

## 2.2 The process of technical change
The process of technical change is chosen such that it satisfies three conditions: First, a trend in $\mu_t$ should be allowed since the relative prices often contains a trend.[^note2] Second, we deviate from the linear trend assumption to account for medium run fluctuations of factor expenditures. Third, $\mu_t$ should be a persistent process, i.e. it cannot reflect the year-to-year noise in the data. A process that fits these conditions is an I(2) process:
$$\Delta \mu_t = \Delta \mu_{t-1} +\eta_t,\quad \eta_t=N(0,\Sigma^\eta) $$
The degree of smoothness is determined by the so-called inverse signal-to-noise ratio, $\lambda\equiv \Sigma^\varepsilon/\Sigma^\eta$. When $\lambda\rightarrow0$ it implies that all noise in the data not captured by the relative prices is subscribed to technical change. Oppositely, when $\lambda\rightarrow\infty$, $\mu_t$ is a linear trend. This illustrates the flexibility of $\mu_t$ and that it nests several other assumptions, e.g. the linear trend or Box-Cox transformation. 

## 2.3 The CESKalman function
The function **CESKalman** estimates equation (2) with $\mu_t$ an I(2) process. The variances $\Sigma^\varepsilon$ and $\Sigma^\eta$ are estimated with maximum likelihood. The function allows both for a free estimation of both parameters and for a fixed value of $\lambda$ such that the estimated parameter is $\Sigma^\varepsilon$ and $\lambda\Sigma^\eta=\Sigma^\varepsilon$. Since the Kalman filter is a Markov Chain, starting values of $\alpha$ and $\sigma$ are needed.[^note3] Therefore, we perform a grid search over different values of these parameters and choose the combination that maximizes the likelihood. If $\sigma$ is estimated to be negative, it is restricted to zero. If $\varepsilon_t$ is either autocorrelated based on a Breush-Godfrey test for autocorrelation or the variances do not satisfy a Normalized Innovations Squared test (NIS), additional lags of the short run parameters are included in equation (2). To summarize, the procedure is:

1. Estimate the variances in equation (2). A grid of different values are used as initial parameter values, specified in grid.param_init. This object is a string of length 3 containing respectively, the lower bound, upper bound and the increment (step size). The grid search continues until convergence.
2. Perform a grid search procedure of initial values of $\sigma$ and $\alpha$, specified in grid.sigma_init and grid.alpha_init. The combination that maximizes the likelihood is chosen. Note that step 1 is performed for every combination.
3. If the estimated elasticity is negative, it is restricted to zero and the procedure returns to step 1 with this restriction imposed.
4. If autocorrelation or a violation of the NIS test is found, lags are included and the procedure returns to step 1.
5. Try different values of $\lambda$. Our preferred method is first a free estimation and next a grid search procedure specified in grid.lambda_init. For every new value of $\lambda$, step 1-4 is caried out.


# 3. The R-package
In this Section, we show how to apply the CESKalman R-package. In Section 3.1 we show how the package is installed. In Section 3.2 we illustrate the role for medium run fluctuations in technical change based on a data set available through the R-package. Section 3.3 describes how to estimate the elasticity of substitution subject to a flexible approach of technical change with the CESKalman package. It also contains a description of the **plot.CESKalman** function which can be used to plot the fit of model along with the decomposition of the fitted relative quantities in prices and technical change. Lastly, the function **CESKalman_Bootstrap** is used to bootstrap standard errors of the parameters. Section 3.4 contains a short analysis of how the degree of smoothness (value of $\lambda$), affects technical change and the elasticity.

## 3.1 Installation of the package
The first step in our analysis is the installation of the package. The package is installed from GitHub: (mangler at finde ud af hvordan man gør det) 

```{r,results='hide',warning=FALSE,message=FALSE}
## Installing the package
library(devtools)
install_github("CKastrup/CESKalman")
library(CESKalman)

```

## 3.2 The role for medium run fluctuations
Having installed the package, we are ready to investigate if medium run fluctuations in factor effeciency is present in the US data from the **Load_Data** function. Data is available from 1970 to 2017 for 16 OECD countries.[^note4] The data series are primarily from the PWT database, but in a few cases, the labor income share is unobserved in the first years of the sample and set equal to the first available observation in the PWT data. When the OECD database contains longer data series on the labor income share, this dataset is used. The function **Load_Data** is only made available to the user in order to illustrate how the package works. We do not take responsibility of potential errors in the function.

The function is applied as:
```{r,warning=FALSE,message=FALSE}
US = Load_Data(Country="USA",tstart=1970,tend=2017) 
US = ts(US,start=1970,end=2017,freq=1)
```

The dataframe US contains many different variables, such as GDP ("q_gdp"), employment ("L"), capital ("K"), a static Hall-Jorgenson user cost of capital without taxes ("q") and the wage ("w"). Some "stylized facts" about the role of the relative prices and technical change are shown in the following. We do not attempt to provide any proof of the relative importance of technical change and the relative prices, but only to motivate that persistent deviations from the long run are present in the data, which motivates the role of non-linear technical change. The following code plots the relative prices, relative quantities, and factor producivities:

```{r,warning=FALSE,message=FALSE}
par(mfrow=c(2,2),mar=c(2,4,2,4),oma = c(0, 0, 3, 0))

plot(US[,'q']/US[,'w'],main="Relative prices, q/w",ylab="Level",xlab="",lwd=2)
  par(new = TRUE)
  plot(100*diff(log(US[,'q']/US[,'w'])),type="l",axes=F,ylab="",xlab="",lty=2,lwd=2)
  axis(4)
  abline(h=0,col="grey40")
  abline(h=mean(100*diff(log(US[,'q']/US[,'w']))),lty=2)
  mtext("Growth rate, %", side=4, line=2.5,cex=0.7)
  legend("topright",legend=c("Level","Growth rate, rhs"),lty=c(1,2),lwd=2,bty="n")
  
  plot(US[,'K']/US[,'L'],main="K/L-ratio",ylab="Level",xlab="",lwd=2)
  par(new = TRUE)
  plot(100*diff(log(US[,'K']/US[,'L'])),type="l",axes=F,ylab="",xlab="",lty=2,lwd=2)
  axis(4)
  abline(h=0,col="grey40")
    abline(h=mean(100*diff(log(US[,'K']/US[,'L']))),lty=2)
  mtext("Growth rate, %", side=4, line=2.5,cex=0.7)
  
  plot(US[,'q_gdp']/US[,'K'],main="Capital productivity, Y/K",ylab="Level",xlab="",lwd=2)
  par(new = TRUE)
  plot(100*diff(log(US[,'q_gdp']/US[,'K'])),type="l",axes=F,ylab="",xlab="",lty=2,lwd=2)
  axis(4)
  abline(h=0,col="grey40")
    abline(h=mean(100*diff(log(US[,'q_gdp']/US[,'K']))),lty=2)
  mtext("Growth rate, %", side=4, line=2.5,cex=0.7)
  
  plot(US[,'q_gdp']/US[,'L'],main="Labor productivity, Y/L",ylab="Level",xlab="",lwd=2)
  par(new = TRUE)
  plot(100*diff(log(US[,'q_gdp']/US[,'L'])),type="l",axes=F,ylab="",xlab="",lty=2,lwd=2)
  axis(4)
  abline(h=0,col="grey40")
    abline(h=mean(100*diff(log(US[,'q_gdp']/US[,'L']))),lty=2)
  mtext("Growth rate, %", side=4, line=2.5,cex=0.7)
  
  mtext("Figure 1: US data on capital and labor",side=3,outer=TRUE,cex=1.5,line=0.5)
  
  
```

In Figure 1 we see that the relative prices, that is, the user cost relative to the wage, has been decreasing through the sample with `r -round(mean(100*diff(log(US[,'q']/US[,'w']))),2)`% on average. This holds in particular true in the 90'es and start 00'es whereafter the relative growth rate has been `r round(-mean(100*diff(log(US[32:47,'q']/US[32:47,'w']))),2)`% from 2003-2017. 

The capital to labor ratio has been increasing through the sample with an average growth rate of `r round(mean(100*diff(log(US[,'K']/US[,'L']))),2)`%. This may reflect that labor is a scarce factor of production and cannot be expanded beyond the population size. Oppositely, capital can be accumulated through further investments. 

In the bottom two graphs of Figure 1 capital and labor productivity, respectively, i.e. the Y/K and Y/L, ratios are displayed. These may reflect changes in the composition of factors due to technical change or relative prices. In the long run, technical change appears to be Harrod-neutral, that is, a positive growth rate of labor-augmenting technical change, and only minor increases in the growth rate of capital-augmenting technical change. This is consistent with labor-augmenting technical change increasing the wage in the long run as seen above. We interpret this in accordance with @Acemoglu2002 and @Acemoglu2003, who argue that since labor is the scarce factor of production, technical change will be directed at this factor. However, there are shifts is the direction of technical change, for an example in the 90'es and after the financial crisis in 2008 where capital-augmenting technical change may have been high.[^note5] Such shifts motivate the role for medium run fluctuations and why technical change should be allowed to evolve non-linearly.


## 3.3 The CESKalman package: Estimation of the elasticity of substitution and technical change
In this Section, we show how to use the CESKalman package to estimate the elasticity and a process of technical change. Prior to estimation, data should be ordered in the following way:
```{r,warning=FALSE,message=FALSE}
data = cbind(US[,"q"],US[,"w"],US[,"K"],US[,"L"])
```

### 3.3.1 Estimation: The CESKalman function

As explained above, the **CESKalman** function is the main function of the package and used to estimate the elasticity along with a process of technical change. It has the advantage that it loops over lots of different combinations of initial parameter values, add lags to ensure no autocorrelation, and chooses $\lambda$ to ensure the best fit of the model, given it is well specified. Our preferred method for choosing $\lambda$ is a free ML estimation combined with a grid searching procedure to address potential issues of non-concavity of the objective function.[^note6] The parameter grid.lambda should contain the, respectively, lower and upper bound of the grid and the step size.[^note7] In the example below with grid.lambda=c(100,1000,100), the lowest value of $\lambda$ applied is 100 and the upper is 1000. The parameter lambda_est_freely should be TRUE or FALSE which indicates if a free estimation of $\lambda$ should be included also.

In the following we try with a grid of $\lambda$ in the interval between 100 and 1000. The step size is chosen to 100 and the grid is combined with a free estimation. By default the function displays output for every $\lambda$ applied. The output is the value of $\lambda$, the likelihood, the value of the elasticity and if NIS=TRUE and BG_test=TRUE it implies that the model is well-specified on a confidence level set by c_val_Auto and c_val_NIS, 10% by default.[^note8] The function is applied as:
```{r,warning=FALSE,message=F}
Estimation = CESKalman(data=data,grid.lambda = c(100,1000,100),lambda_est_freely = T)


```
We see that the free estimation implies a value of $\lambda$ close to one, which will ascribe a large part of fluctuations in the relative expenditure shares not described by prices to technical change. However, this estimation with a low value of $\lambda$ results in a misspecified model which is a result consistent with the estimations in @Kronborg2019. This highlights the importance of combining the free estimation with the grid search. When $\lambda=100$, the likelihood is conciderably higher than the free estimation and the model is well specified. The elasticity is 0.46. As the value of $\lambda$ is increased, the likelihood increases along with the value of $\sigma$. It is a general result that as we converge to a linear trend ($\lambda\rightarrow\infty$), more of the variation in expenditure shares has to be explained by the relative prices which increases the elasticity. The optimal value of $\lambda$ is 600 where the elasticity is 0.57. When $\lambda=700$ or above, the returned likelihood is decreasing, but the elasticity is higher. Since all estimations in the grid search procedure are well specified, the preferred estimation returned is the one with the highest likelihood value, $\lambda=600$. 

In Figure 2 below we plot the process of technical change, $log(\Gamma_t)$, returned by the function. In the long run, technical change is relatively augmenting labor, as predicted from the moviation in Section 3.2. However, there are important periods of medium run fluctuations, for an example in the 90'es potentially due to the IT-revolution (@Klump2008 on EU data) and a decline in the growth rate after the financial crisis.[^note9]
```{r,warning=FALSE,message=FALSE}
plot(ts(Estimation$Gamma,start=1970+Estimation$nlags,end=2016,freq=1),ylab="",xlab="",lwd=2
     ,main="Figure 2: The process of technical change")

```

### 3.3.2 Displaying output: The plot function
The primary variable we want the model to fit is the relative quantities.[^note10] The fit of the model can be investigated by using the **plot.CESKalman** function. The plot contains three different graphs. The first shows the model fit of the relative quantities. The second shows a decomposition of the fitted values in a trend-component reflecting technical change and a price-component reflecting short and long run variation in prices. The two components are defined as differential equations such that the lag structure of the relative quantities in equation (2) is taken into account. The last plot shows the relative prices and relative quantities applied in estimation. The relative price is shown on the right axis with a negative sign in front. If the two variables are possitively correlated, it reflects that the estimated elasticity should be positive.

Unfortunately, the function is currently only prepared to handle estimations where no lags of the growth rate in expenditure shares are included (i.e. k=0 in equation (2)). Since one lag is included in the estimation above, the function cannot be used. We can however still illustrate the functionality of the plot by using a value of $\lambda$ which implies that no lags are added. Such a value is $\lambda=6.25$ which is applied below.[^note11] The resulting plot is shown in Figure 3. From the first graph it can be seen that the fit of the model is good and that the residuals do not appear to be autocorrelated. In the second graph, the decomposition of the fit of relative quantities in a trend and price component is shown. It can be seen that the price component explains the largest part of the fitted values, whereas the trend component appear at first to be constant throughout the sample. However, by a closer expection it can be seen that the trend component is decreasing in the 90'es and after the financial crisis. Since $\sigma<1$, this reflects capital-augmenting technical change. The last graph can be used to investigate the correlation between demand and prices. The correlation between these are high which motivates the importance of the price component and why  positive (and relatively high) long run elasticity is found.

```{r,warning=FALSE,message=FALSE}
Estimation_2 = CESKalman(data=data,grid.lambda = c(6.25,6.25,100),lambda_est_freely = F)
plot(Estimation_2,t0=1970,tEnd=2017,main="Figure 3: The plot.CESKalman function")
```
 
### 3.3.3 Bootstrapping
Since $\alpha$ and $\sigma$ is specified as constant in the Kalman filter, we do not obtain standard errors of these directly. Therefore, a residual-based recursive bootstrapping procedure is used to obtain standard errors. The function applied is **CESKalman_Bootstrap**. The function returns the bootstrapped estimates of $\alpha$ and $\sigma$. When print_results=TRUE (the default), a set of summary statistics about the bootstrapped values and a density plot is returned. As default, ndraw=1000. The function is applied as:
```{r,warning=FALSE,message=F}
Bootstrap = CESKalman_Bootstrap(Estimation)
```
The bootstrapped standard error of $\sigma$ is `r round(sd(Bootstrap[,"sigma"]),2)` implying that the estimate of $\sigma$ is significantly different from zero and one. Thus, we are both able to reject that production is Leontief (perfect complements) or Cobb-Douglas. 
 
## 3.4 How the degree of smoothness affects the process of technical change
In this Section, we show how the degree of smoothness affects the value of the elasticity and the process of technical change. In the code below we estimate the elasticity with $\lambda=100$, $\lambda=1000$ and $\lambda=10000$. The elasticity increases from 0.46 to 0.78 when increasing $\lambda$ from 100 to 1000. When $\lambda=10000$, the elasticity decreases to 0.68 which is due to the inclusion of an additional lag (k=2 in equation (2)) to ensure no autocorrelation in the residuals. In general, it seems that as $\lambda$ is increased, so is the risk of autocorrelation. Therefore, more lags are added. As an example, when $\lambda=6.25$ as above, zero lags are included, one lag is included when $\lambda\in(100,600,1000)$ and two lags when $\lambda=10000$. Thus, holding the number of lags constant, the elasticity is increasing in $\lambda$.
```{r,warning=FALSE,message=F}
Estimation_100 = CESKalman(data=data,grid.lambda = c(100,100,100),lambda_est_freely = F)
Estimation_1000 = CESKalman(data=data,grid.lambda = c(1000,1000,100),lambda_est_freely = F)
Estimation_10000 = CESKalman(data=data,grid.lambda = c(10000,10000,100),lambda_est_freely = F)
```

Below we compare the process of technical change for $\lambda\in(6.25,100,600,1000,10000)$. We see that the process of technical change is most volatile for low values of $\lambda$ and may feature several periods of changes in the direction of technical change. This is most easily seen when $\lambda=6.25$ with particular large changes in the direction of technical change in the 90'es and after the financial crisis. When $\lambda=10000$, the process of technical change is almost a linear trend, which shows that as $\lambda\rightarrow\infty$, the process converges to a linear trend.

```{r,warning=FALSE,message=FALSE}

ylim = c(Estimation_2$Gamma,Estimation_100$Gamma,Estimation$Gamma,
         Estimation_1000$Gamma,Estimation_10000$Gamma)

plot(ts(Estimation_2$Gamma,start=1970+Estimation_2$nlags,end=2016,freq=1),ylab="",xlab="",lwd=2
     ,ylim=range(ylim),main="Figure 4: Different degrees of smoothness")
lines(ts(Estimation_100$Gamma,start=1970+Estimation_100$nlags,end=2016,freq=1),lty=2,lwd=2)
lines(ts(Estimation$Gamma,start=1970+Estimation$nlags,end=2016,freq=1),lty=3,lwd=2)
lines(ts(Estimation_1000$Gamma,start=1970+Estimation_1000$nlags,end=2016,freq=1),lty=4,lwd=2)
lines(ts(Estimation_10000$Gamma,start=1970+Estimation_10000$nlags,end=2016,freq=1),lty=5,lwd=2)

legend("topleft",legend=c(expression(paste(lambda,"=6.25")),expression(paste(lambda,"=100")),
                          expression(paste(lambda,"=600")),expression(paste(lambda,"=1000")),
                          expression(paste(lambda,"=10000"))),lty=c(1,2,3,4,5),lwd=2)
```

# 4. Concluding remarks
In this note, we have shown how the R-package CESKalman can be used to estimate the elasticity of substitution subject to a flexible, data-driven, process of technical change. The package uses the fact that the process have a natural state-space representation and uses a Kalman filter to estimate the elasticity and a process of technical change. The method uses several grid searching procedures of the model parameters which ensures that the model with the highest likelihood is chosen, conditional on being well specified. 

Using a dataset made available to the user of the package, we show that the US elasticity of substitution between capital and labor is 0.57 in our preferred specification and that technical change is labor augmenting in the long run, but with several periods of changes in the direction of technical change.



# References
[^note1]: Note that $\mu_t$ and $log(\Gamma_t)$ are of opposite signs when $\sigma<1$. That is, when $\sigma<1$, an increase in $\mu_t$ reflects an increase in labor-augmenting technical change.
[^note2]: If no trend is included in the direction of technical change, but $\sigma\neq1$, it implies that the relative expenditure shares are trending due to a trend in the relative prices, which is theoretically inconsistent.
[^note3]: The initial values of short run dynamics are set to zero. The initial value of $\mu_{t-1}$ is chosen by assuming that the economy is in the long run equilibrium, that is, the initial value is $\mu_{t-1}=s_{t-1}-(1-\sigma)p_{t-1}$ with $\sigma$ being the initial value of $\sigma$. The variances of all parameters are set relatively high to secure a fast convergence as argued in @Durbin2012.
[^note4]: See ?Load_Data for description of the methodology and available countries.
[^note5]: obviously, we can only speculate if the increase in the Y/K ratio is due to technical change or the relative prices
[^note6]: As we shall see below, it is also possible to fix $\lambda$ at a certain value and thereby not apply the free estimation and the grid search procedure.
[^note7]: If grid.lambda is set to NA, no grid search is performed.
[^note8]: The graphical output can be turned off by setting print_results=FALSE.
[^note9]: If $\lambda$ is set to 6.25, technical change would be augmenting capital after the financial crisis as shown below. 
[^note10]: We only use expenditure shares in estimation instead of quantities in order to account for potential measurement errors.
[^note11]: In fact, @Ravn2002 argue that the optimal value in the HP-filter is 6.25 on yearly US data.
