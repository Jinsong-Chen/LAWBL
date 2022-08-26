#' LAWBL: Latent (Variable) Analysis with Bayesian Learning
#'
#' @description This package is to provide a variety of models to analyze latent variables based on Bayesian learning.
#'
#'
#' @details
#'
#' \emph{LAWBL} represents a partially confirmatory / exploratory approach to model latent variables based on Bayesian learning.
#' Built on the power of statistical learning, it can address psychometric challenges such as parameter specification, local dependence,
#' and factor extraction. Built on the scalability and flexibility of Bayesian inference and resampling techniques, it can accommodate
#' modeling frameworks such as factor analysis, item response theory, cognitive diagnosis modeling and causal or explanatory modeling.
#' The package can also handle different response formats or a mix of them, with or without missingness. The variety of models provide
#' a partial approach covering a wide range of the exploratory-confirmatory continuum under the context of latent variable modeling.
#'
#' Towards the confirmatory end, this package includes the Partially Confirmatory Factor Analysis (PCFA) model for continuous data
#' (Chen, Guo, Zhang, & Pan, 2020), the generalized PCFA (GPCFA) model covering continuous, categorical, and mixed-type data,
#' and the partially confirmatory item response model (PCIRM) for continuous and dichotomous data with intercept terms (Chen, 2020).
#' For PCFA, GPCFA, and PCIRM, there are two major model variants with different constraints for identification. One assumes local
#' independence (LI) with a more exploratory tendency, which can be also called the E-step. The other allows local dependence (LD)
#' with a more confirmatory tendency, which can be also called the C-step.
#'
#' Towards the exploratory end, the Bayesian regularized EFA (BREFA) with factor extraction and parameter estimation in one step (Chen 2021)
#' is offered. It's further improved as the Fully and partially EFA with better performance and partial knowledge.
#'
#' Parameters are obtained by sampling from the posterior distributions with the Markov chain Monte Carlo (MCMC) techniques. Different
#' Bayesian learning methods are used to regularize the loading pattern, local dependence, and/or factor identification.
#'
#'
#' @name LAWBL-package
#' @aliases LAWBL-package
#' @docType package
#'
#' @author {Jinsong Chen, \email{jinsong.chen@live.com}}
#'
#' @note This package is under development. You are very welcome to send me any comments or
#'  suggestions for improvements, and to share with me any problems you may encounter with
#'   the use of this package.
#'
#' @import stats
#' @import MASS
#' @import coda
#'
#' @references
#'
#' Chen, J. (2020). A partially confirmatory approach to the multidimensional item response theory with the Bayesian Lasso.
#' \emph{Psychometrika}. 85(3), 738-774. DOI:10.1007/s11336-020-09724-3.
#'
#' Chen, J., Guo, Z., Zhang, L., & Pan, J. (2021). A partially confirmatory approach to scale development
#'  with the Bayesian Lasso. \emph{Psychological Methods}. 26(2), 210-235. DOI: 10.1037/met0000293.
#'
#' Chen, J. (2021). A generalized partially confirmatory factor analysis framework with mixed Bayesian Lasso methods.
#'  \emph{Multivariate Behavioral Research}. DOI: 10.1080/00273171.2021.1925520.
#'
#' Chen, J. (2021). A Bayesian regularized approach to exploratory factor analysis in one step.
#' \emph{Structural Equation Modeling: A Multidisciplinary Journal}. DOI: 10.1080/10705511.2020.1854763.
#'
#' Chen, J. (2022). Partially confirmatory approach to factor analysis with Bayesian learning: A LAWBL tutorial.
#'  \emph{Structural Equation Modeling: A Multidisciplinary Journal}. DOI: 10.1080/00273171.2021.1925520.
#'
#' Chen, J. (In Press). Fully and partially exploratory factor analysis with bi-level Bayesian regularization.
#' \emph{Behavior Research Methods}.
#'
#' @keywords package
NULL

