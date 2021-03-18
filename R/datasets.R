#' Simulated CFA data with LI
#'
#' CFA data simulated based on 18 items, 3 factors and local independence;
#'  factorial correlation \eqn{\Phi=.3}.
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{A dataset with simulated responses of 1000 individuals to 18 items}
#' \item{\code{qlam}}{Loading pattern and values used to simulated the data}
#' }
#'
"sim18cfa0"


#' Simulated CFA data with LD
#'
#' CFA data simulated based on 18 items, 3 factors and local dependence; factorial correlation \eqn{\Phi=.3}.
#'
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{A dataset with simulated responses of 1000 individuals to 18 items}
#' \item{\code{qlam}}{Loading pattern and values used to simulated the data}
#' \item{\code{LD}}{Local dependence between items (LD effect = .3)}
#' }
#'
"sim18cfa1"

#' Simulated CCFA data with LI and missingness
#'
#' Categorical CFA data simulated based on 18 items, 3 factors, and 4 categories
#' with local independence and  10 percent missingness at random; factorial correlation \eqn{\Phi=.3}.
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{A dataset with simulated responses of 1000 individuals to 18 items}
#' \item{\code{qlam}}{Loading pattern and values used to simulated the data}
#' }
#'
"sim18ccfa40"

#' Simulated CCFA data with LD and missingness
#'
#' Categorical CFA data simulated based on 18 items, 3 factors, and 4 categories
#' with local dependence and  10 percent missingness at random; factorial correlation \eqn{\Phi=.3}.
#'
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{A dataset with simulated responses of 1000 individuals to 18 items}
#' \item{\code{qlam}}{Loading pattern and values used to simulated the data}
#' \item{\code{LD}}{Local dependence between items (LD effect = .3)}
#' }
#'
"sim18ccfa41"

#' Simulated MCFA data with LD and Missingness
#'
#' CFA data mixed with continuous and categorical responses simulated based on 3 factors,
#' 6 4-category items, 12 continuous items, local dependence, and 10 percent missigness at random;
#' factorial correlation \eqn{\Phi=.3}.
#'
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{A dataset with simulated responses of 1000 individuals to 18 items}
#' \item{\code{qlam}}{Loading pattern and values used to simulated the data}
#' \item{\code{LD}}{Local dependence between items (LD effect = .3)}
#' }
#'
"sim18mcfa41"


#' Simulated CCFA data (dichotomous) with LD and a minor factor/trait
#'
#' Categorical CFA data simulated based on 24 items, 4 factors, 2 categories
#'  and local dependence; factorial correlation \eqn{\Phi=.3}.
#'  The last factor/trait is minor (measured by cross-loadings only).
#'
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{A dataset with simulated responses of 1000 individuals to 24 items}
#' \item{\code{qlam}}{Loading pattern and values used to simulated the data}
#' \item{\code{LD}}{Local dependence between items (LD effect = .3)}
#' }
#'
"sim24ccfa21"

#' National Longitudinal Survey of Youth 1997
#'
#' A data set consisted of 3,458 individual responses to 27 mixed-type items, with a 1.12 percentage of missing data
#'
#'
#' @format A list with components:
#' \describe{
#' \item{\code{dat}}{The response data}
#' \item{\code{Q}}{Intial design matrix with three factors and two to three specified loadings per factor}
#' \item{\code{cati}}{Indices of categorical (polytomous) items}
#' }
#'
"nlsy27"

