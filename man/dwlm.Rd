\name{dwlm}
\alias{dwlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solves the doubly weighted simple linear model
}
\description{
Fits the simple linear model using weights on both the predictor and the response
}
\usage{
dwlm(x, y, weights.x = rep(1, length(x)), 
	weights.y = rep(1, length(y)), subset = rep(TRUE, length(x)), 
	sigma2.x = rep(0, length(x[subset])), 
	from = min((y[subset] - mean(y[subset]))/(x[subset] - mean(x[subset]))), 
	to = max((y[subset] - mean(y[subset]))/(x[subset] - mean(x[subset]))), 
	n = 1000, max.iter = 100, tol = .Machine$double.eps^0.25, 
	method = c("MLE", "SSE", "R"), trace = FALSE, coef.H0 = c(0,1), alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	the predictor values
}
  \item{y}{
	the response values
}
  \item{weights.x}{
	the weight attached to the predictor values
}
  \item{weights.y}{
	the weight attached to the response values
}
  \item{subset}{
	a logical vector or a numeric vector with the positions to be considered
}
  \item{sigma2.x}{
	numeric, the variance due to heterogeneity in the predictor value
}
  \item{from}{
	numeric, the lowest value of the slope to search for a solution
}
  \item{to}{
	numeric, the highest value of the slope to search for a solution
}
  \item{n}{
	integer, the number of slices the search interval (from, to) is divided in
}
  \item{max.iter}{
	integer, the maximum number of allowed iterations
}
  \item{tol}{
	numeric, the maximum allowed error tolerance
}
  \item{method}{
	string, the selected method (MSE, SSE, R) as described in the references.
}
  \item{trace}{
	logical, flag to keep track of the solution
}
  \item{coef.H0}{
	numeric vector, the coeffients to test against to for significant difference
}
  \item{alpha}{
	numeric, the significance level for estimating the Degrees of Equivalence (DoE)
}
}
\value{
	A list with the following elements:
	\item{x}{original pedictor values}
	\item{y}{original response values}
	\item{weights.x}{original predictor weigths}
	\item{weights.y}{original response weights}
	\item{subset}{original subset parameter}
	\item{coef.H0}{original parameter value for hypothesis testing against to}
	\item{coefficients}{estimated parameters for the linear model solution}
	\item{cov.mle}{Maximum Likelihood Estimafor for the covariance matrix}
	\item{cov.lse}{Least Squares Estiimator for the  covariance matrix}
	\item{x.hat}{fitted true predictor value, this is a latent (unobserved) variable}
	\item{y.hat}{fitted true response value, this is a latent (unobserved) variable}
	\item{df.residuals}{degrees of freedom}
	\item{MSE}{mean square error of the solution}
	\item{DoE}{pointwise degrees of equivalente between the observed and the latent variables}
	\item{u.DoE.mle}{uncerainty of the pointwise degrees of equivalence using maximum likelihood solution}
	\item{u.DoE.lse}{uncertainty of the pointwise degrees of equivalence using least squares solution}
	\item{dm.dXj}{partial gradient of the slope with respect to the jth predictor variable}
	\item{dm.dYj}{partial gradient of the slope with respect to the jth response variable}
	\item{dc.dXj}{partial gradient of the intercept with respect to the jth predictor variable}
	\item{dc.dYj}{partial gradient of the intercept with respect to the jth response variable}
	\item{curr.iter}{number of iterations}
	\item{curr.tol}{reached tolerance}
}
\references{
Reed, B.C. (1989) "Linear least-squares fits with errors in both coordinates", American Journal of Physics, 57, 642. \url{https://doi.org/10.1119/1.15963} 

Reed, B.C. (1992) "Linear least-squares fits with errors in both coordinates. II: Comments on parameter variances", American Journal of Physics, 60, 59. \url{https://doi.org/10.1119/1.17044}

Ripley and Thompson (1987) "Regression techniques for the detection of analytical bias", Analysts, 4. \url{https://doi.org/10.1039/AN9871200377}
}
\author{
Hugo Gasca-Aragon

Maintainer: Hugo Gasca-Aragon <hugo_gasca_aragon@hotmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{lm}}
}
\examples{
# Example ISO 28037 Section 7
X.i<- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
Y.i<- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
sd.X.i<- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
sd.Y.i<- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)

# values obtained on sep-26, 2016
dwlm.res <- dwlm(X.i, Y.i, 1/sd.X.i^2, 1/sd.Y.i^2, 
	from = 0, to=3, coef.H0=c(0, 2), tol = 1e-10)
dwlm.res$coefficients
dwlm.res$cov.mle[1, 2]
dwlm.res$curr.tol
dwlm.res$curr.iter
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dwlm }
% __ONLY ONE__ keyword per line
