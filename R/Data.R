# ===================================================================
# ------ DATA FROM THE PWT AND OECD USED IN FACTOR ESTIMATIONS ------
# ===================================================================

# ===================================================================
#' @author Christian Sandholm Kastrup <CST@dreamgruppen.dk>, Anders Farver Kronborg <ANK@dreamgruppen.dk> and Peter Philip Stephensen <PSP@dreamgruppen.dk>
# by Christian Sandholm Kastrup
# Latest update: 05-12-2019

# ===================================================================
#' @title Data series for 16 OECD countries obtained from the PWT and OECD database
#'
#' @description This function loads data from the PWT and OECD databases and is suitable for estimating the elasticity of substitution between capital and labor at the country level.
#'
#' @details This function loads various types of data and returns a time series matrix for the chosen country.
#' \strong{Note:} We do not take responsibility for potential errors and the data is primarily included for illustration.
#'
#' The data series are obtained from the PWT version 9.1 and the OECD data base.
#'
#' Available countries are Australia (AUS), Austria (AUT), Belgium (BEL), Canada (CAN), Denmark (DNK), Finland (FIN)
#' France (FRA), Germany (DEU), Italy (ITA), Japan (JPN), the Netherlands (NLD), New Zealand (NZL), Norway (NOR),
#' Sweden (SWE), Great Britain (GBR) and the United States (USA).
#'
#' The user cost of capital is calculated as: \eqn{q_{t}=(v_{it}/q_{it})*(irr_{t}+delta_{t})} with i denoting investments
#'
#' @param Country The ISO code of the country (see details)
#' @param tstart Initial time period. Earliest is 1970 (default)
#' @param tend Last time period. Latest possible is 2017 (default)
#'
#' @usage Load_Data(Country,tstart=1970,tend=2017)
#'
#' @return
#' Returns a time series matrix with the objects:
#' \strong{q_gdp:} Real GDP in national currency.
#' \strong{v_gdp:} Nominal GDP in national currency.
#' \strong{emp:} Employment in national currency.
#' \strong{avh:} Average number of yearly working hours.
#' \strong{labsh:} Labor share in nominal GDP.
#' \strong{irr:}  Real Internal Rate of Return (see PWT version 9.1).
#' \strong{delta:} Residually calculated depreciation rate.
#' \strong{v_i:} Nominal investments in national currency.
#' \strong{q_i:} Real investments in national currency.
#' \strong{K:} Net capital stock.
#' \strong{L:} Number of yearly labor hours.
#' \strong{w:} The wage.
#' \strong{pi:} Inflation rate.
#' \strong{q:} A simplified Hall-Jorgenson user cost (see details).
#' \strong{KL:} Capital/labor ratio.
#' \strong{markup:} The calculated profits rate.
#' \strong{GDP_US:} Real GDP in USD. Can be used in weighted regressions.
#'
#' @examples
#' Data = Load_Data(Country="USA",tstart=1970,tend=2017)
#'
#' Data = data.frame(Data)
#'
#' ## Create the data object needed in the CESKalman function:
#' data = cbind(Data$q,Data$w,Data$K,Data$L)
#'
#'
#' @seealso
#'
#' https://stats.oecd.org/Index.aspx?DataSetCode=STAN and  https://www.rug.nl/ggdc/productivity/pwt/
#'
#' @references Feenstra et al (2015)
#'
#' @export


# ===================================================================
# Starting function
# ===================================================================

Load_Data = function(Country,tstart=1970,tend=2017){

  load("Data.RData")

  tmp = get(paste(Country))

  tmp = window(tmp,start=tstart,end=tend)

    return(data.frame(tmp))

}
