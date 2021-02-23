# The CESKalman package
Welcome to the Github page of the R-package CESKalman. The package is used to estimate CES production functions in time series data. 

The main advantage of this package relative to other approches in the literature (linear trend assumption, Box-Cox transformation ect.) is that the process of technical change, known to potentially bias the elasticity estimate, is estimated in a flexible manner as an I(2) process. The inverse signal-to-noise ratio (the ratio of measurement error variance to variance of the process of technical change), determines the smoothness of the process: High values imply that technical change converges to a linear trend. Oppositely, low values imply that all year-to-year noise in the data not explained by the relative prices are subscribed to technical change. The optimal value is found as the value that maximizes the likelihood based on first a free estimation and next a grid search procedure to address potential non-concavity of the objective function.

For information on the methodology and how to apply the package, the reader is referred to the [Reference_manual](https://github.com/CKastrup/CESKalman/blob/master/Reference_manual.pdf) and the [Vignette](https://github.com/CKastrup/CESKalman/blob/master/VIGNETTE.pdf) made available to the user from this Github page along with the [working paper version](https://dreamgruppen.dk/media/9274/w2019_06.pdf). We are happy to take comments, questions, or recommendations on mail CST@dreammodel.dk. 

We hope that the package can be an improvement to future researchers estimating the elasticity of substitution.

Christian S. Kastrup (mail: CST@dreammodel.dk), Anders F. Kronborg (mail: ANK@dreammodel.dk) and Peter P. Stephensen (mail: PSP@dreammodel.dk)

# Usage
How to apply the package is described in detail in the Reference_manual and the Vignette. The package is installed as:

```{r,results='hide',warning=FALSE,message=FALSE}
## Installing the package
library(devtools)
install_github("CKastrup/CESKalman")
library(CESKalman)

```

A tutorial dataset is made available in the R-package via the Load_Data function to illustrate how the package can be applied. This dataset is based on the PWT database (Feenstra et al., 2015):

```{r,warning=FALSE,message=FALSE}
US = Load_Data(Country="USA",tstart=1970,tend=2017) 
```

The returned object 'US' returns a time series object with variables such as the user cost, 'q', the wage, 'w', the capital stock, 'K', and the labor supply 'L'.

The main function in the package is *CESKalman* and applied in the following manner:

```{r,warning=FALSE,message=F}
data = cbind(US[,'q'],US[,'w'],US[,'K'],US[,'L'])
Estimation = CESKalman(data=data,grid.lambda = c(100,1000,100),lambda_est_freely = T)

Estimation$sigma ## This is the elasticity
Estimation$alpha ## This is the adjustment parameter in the ECM
Estimation$Gamma ## This is the relative log technical change
```

The function loops over a grid of different values of the inverse signal-to-noise ratio, $\lambda$. In the above estimation, it is first estimated freely since lambda_est_freely=TRUE, next a grid search procedure is perform on the interval from 100 to 1000 with 100 as step size. The value that maximizes the likelihood conditional on being well specified is returned as the preferred estimation.

A plot function is also available in the package. The plot contains three different graphs: The first evaluates the fit of the model, the second shows the process of relative augmenting technical change and the last displays the data series applied in estimation, but with their means subtracted. The function is applied as: 

```{r,warning=FALSE,message=FALSE}
plot(Estimation,main="The plot.CESKalman function",t0=1970,tEnd=2017)
```

Since the elasticity is specified as a state variable in the Kalman filter with a zero variance (else it would be time-varying), standard errors are obtained from a recursive design residual-based bootstrapping procedure:

```{r,warning=FALSE,message=FALSE}
Bootstrap = Bootstrap(Estimation)

sd_sigma = sd(Bootstrap[,"sigma"])
sd_alpha = sd(Bootstrap[,"alpha"])
```
