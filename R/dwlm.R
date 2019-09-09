#
# Doubly Weighted Functional Linear Regression
#
# [X_i,Y_i]' ~ N([x.i, alpha+beta*x.i]', Diagonal[sigma_xi^2, sigma_yi^2])
#
# (X, Y): are the version of (x, y=alpha+beta*x) observed with measurement error.
#
# where:
#	 mu_i are fixed unknown.
#	 sigma_xi^2,sigma_yi^2 are known.
#
# Gasca (2011)
#
# References: Deming(1943), York(1966), Williamson(1968), Reed(1992)
#
# extended to estimate the degrees of equivalence.
# extended to compute the deleted residuals (Leave One Out) regression.
# extended to compute the multiple deleted residuals (Leave Many Out) regression.
# 2016-08-13: modified to return the coef.H0 parameter in the result.
# 2017-05-25: extended to compute the inverse prediction with its uncertainty
# 2019-08-02: edited to upload it to CRAN
#

dwlm <- function(x, y, weights.x = rep(1, length(x)), 
	weights.y = rep(1, length(y)), subset = rep(TRUE, length(x)), 
	sigma2.x = rep(0, length(x[subset])), 
	from = min((y[subset] - mean(y[subset]))/(x[subset] - mean(x[subset]))), 
	to = max((y[subset] - mean(y[subset]))/(x[subset] - mean(x[subset]))), 
	n = 1000, max.iter = 100, tol = .Machine$double.eps^0.25, 
	method = c("MLE", "SSE", "R"), trace = FALSE, coef.H0 = c(0, 1), 
	alpha = 0.05) {
#
# sigma2.x: variance of the material (heterogeneity).
# method: 	"MLE"=Maximum Likelihood Estimator, "SSE"=Least Squares Estimator, "R"=SSE'.
#		uses Newton Raphson method as approximation method for MLE and SSE methods.
#		uses binary search method for R method.
#
	log <- data.frame(beta = rep(NA, max.iter), R = rep(NA, max.iter), 
		SSE = rep(NA, max.iter), LLH = rep(NA, max.iter))
	b.interval <- from + c(0:(n - 1))/(n - 1)*(to - from)
	method <- method[1]

	old.x <- x
	old.y <- y
	old.w.x.i <- weights.x
	old.w.y.i <- weights.y
	old.lambda <- old.w.y.i/old.w.x.i

	x <- x[subset]
	y <- y[subset]
	weights.x <- weights.x[subset]
	weights.y <- weights.y[subset]

	p.b.R <- rep(0, n)
	SSE <- p.b.R
	llh <- p.b.R
	w.x.i <- weights.x
	w.y.i <- weights.y
	lambda <- w.y.i/w.x.i

	# this locate a rough estimate
	for (kk in 1:n) {
		b <- b.interval[kk]
		w.i <- w.x.i*w.y.i/(w.x.i + b^2*w.y.i)
		x.bar <- sum(w.i*x)/sum(w.i)
		y.bar <- sum(w.i*y)/sum(w.i)
		u.i <- x - x.bar
		v.i <- y - y.bar

		# LS solution: Reed 1992, verified by Gasca 2011
		p.b.R[kk] <- b^2*sum(w.i^2*u.i*v.i/w.x.i) -
			sum(w.i^2*u.i*v.i/w.y.i) +
			b*sum(w.i^2*(u.i^2/w.y.i - v.i^2/w.x.i-2*(b*u.i - v.i)^2*sigma2.x))
		SSE[kk] <- sum(w.i*(v.i - b*u.i)^2)

		a <- y.bar - b*x.bar
		x.hat <- (x + lambda*b*(y - a))/(1 + lambda*b^2)
		y.hat <- a + b*x.hat

		# MLE solution: Gasca 2011, llh=Log LikeliHood
		llh[kk]<- -length(x)*(log(2*pi)) + sum(log(w.x.i))/2 + 
			sum(log(w.y.i))/2 - sum(w.x.i*(x - x.hat)^2)/2 - 
			sum(w.y.i*(y - y.hat)^2)/2
	}

	# now we improve the certainty of the estimate
	if (method == "MLE") {
		idx <- c(1:length(b.interval))[llh == max(llh)]
	} else 
	if (method == "SSE") {
		idx <- c(1:length(b.interval))[SSE == min(SSE)]
	} else 
	if (method == "R") {
		idx <- c(1:length(b.interval))[abs(p.b.R) == min(abs(p.b.R))]
		if (p.b.R[idx]*p.b.R[idx + 1] < 0) {
			idx <- idx
		} else
		if (p.b.R[idx - 1]*p.b.R[idx] < 0) {
			idx <- idx - 1
		} else {
			stop("Error: method R cannot find root.")
		}
	} else {
		stop("Error: method must be [MLE, SSE, R]")
	}

	if (length(idx) > 1) { stop("multiple extrema points detected") }

	if (idx == 1) {
		b.l <- b.interval[idx]
		b.r <- b.interval[idx + 1]
		R.l <- p.b.R[idx]
		R.r <- p.b.R[idx + 1]
		SSE.l <- SSE[idx]
		SSE.r <- SSE[idx + 1]
		llh.l <- llh[idx]
		llh.r <- llh[idx + 1]
	} else {
		if (idx == length(b.interval)) {
			b.l <- b.interval[idx - 1]
			b.r <- b.interval[idx]
			R.l <- p.b.R[idx - 1]
			R.r <- p.b.R[idx]
			SSE.l <- SSE[idx - 1]
			SSE.r <- SSE[idx]
			llh.l <- llh[idx - 1]
			llh.r <- llh[idx]
		} else {
			b.l <- b.interval[idx - 1]
			b.r <- b.interval[idx + 1]
			R.l <- p.b.R[idx - 1]
			R.r <- p.b.R[idx + 1]
			SSE.l <- SSE[idx - 1]
			SSE.r <- SSE[idx + 1]
			llh.l <- llh[idx - 1]
			llh.r <- llh[idx + 1]
		}
	}

	curr.iter <- 1
	curr.tol <- Inf
	while (curr.iter < max.iter & curr.tol > tol) {
		new.b <- (b.l + b.r)/2
		w.i <- w.x.i*w.y.i/(w.x.i + new.b^2*w.y.i)
		x.bar <- sum(w.i*x)/sum(w.i)
		y.bar <- sum(w.i*y)/sum(w.i)
		u.i <- x - x.bar
		v.i <- y - y.bar

		# Reed solution 1992, verified by Gasca 2011
		# extended to include the variance of x.i for the ultrastructural model
		new.p.b.R <- new.b^2*sum(w.i^2*u.i*v.i/w.x.i) - 
			sum(w.i^2*u.i*v.i/w.y.i) + 
			new.b*sum(w.i^2*(u.i^2/w.y.i - v.i^2/w.x.i - 2*sigma2.x*(new.b*u.i - v.i)^2))

		new.SSE <- sum(w.i*(v.i - new.b*u.i)^2)

		a <- y.bar - new.b*x.bar
		lambda <- w.y.i/w.x.i # var.x/var.y
		x.hat <- (x+lambda*new.b*(y - a))/(1 + lambda*new.b^2)
		y.hat <- a + new.b*x.hat

		# MLE solution: Gasca 2011
		new.llh <- -length(x)*(log(2*pi)) + sum(log(w.x.i))/2 + 
			sum(log(w.y.i))/2 - sum(w.x.i*(x - x.hat)^2)/2 - 
			sum(w.y.i*(y - y.hat)^2)/2

		log[curr.iter, ] <- c(new.b, new.p.b.R, new.SSE, new.llh)

		# the search is wrong
		if (method == "MLE") {
			q.b <- -((llh.l - new.llh)*(b.l^2 - b.r^2) - (llh.l - llh.r)*
				(b.l^2 - new.b^2))/((llh.l - llh.r)*(b.l - new.b) - 
				(llh.l - new.llh)*(b.l - b.r))/2
			if (new.b < q.b) { 
				llh.l <- new.llh
				b.l <- new.b
			} else {
				llh.r <- new.llh
				b.r <- new.b
			}
		} else 
		if (method == "SSE") {
			q.b <- -((SSE.l - new.SSE)*(b.l^2 - b.r^2) - (SSE.l - SSE.r)*
				(b.l^2 - new.b^2))/((SSE.l - SSE.r)*(b.l - new.b) - 
				(SSE.l - new.SSE)*(b.l - b.r))/2
			if (new.b < q.b) { 
				SSE.l <- new.SSE
				b.l <- new.b
			} else {
				SSE.r <- new.SSE
				b.r <- new.b
			}
		} else {
			if (R.r*new.p.b.R < 0) { 
				R.l <- new.p.b.R
				b.l <- new.b
			} else {
				R.r <- new.p.b.R
				b.r <- new.b
			}
		}
		curr.tol <- abs(b.r - b.l)/((b.l + b.r)/2)
		curr.iter <- curr.iter + 1
	}

	if (curr.iter >= max.iter & curr.tol > tol) {
		warning("Estimates did not converge.")
	}

	# once found we can estimate the slope, we can estimate the rest of the parameters.

	# update the values at the minumim SSE/llh within the control interval
	slope <- new.b
	w.i <- w.x.i*w.y.i/(w.x.i + slope^2*w.y.i)
	x.bar <- sum(w.i*x)/sum(w.i)
	y.bar <- sum(w.i*y)/sum(w.i)
	u.i <- x - x.bar
	v.i <- y - y.bar
	intercept <- y.bar - slope*x.bar
	w <- sum(w.i)

	# standard error of the estimates, first order approximation by the delta method
	# Reed solution 1992, verified by Gasca 2011.
	HH <- -2*slope/w*sum(w.i^2*v.i/w.x.i)
	JJ <- -2*slope/w*sum(w.i^2*u.i/w.x.i)
	AA <- 4*slope*sum(w.i^3*u.i*v.i/w.x.i^2) - w*HH*JJ/slope
	BB <- -sum(w.i^2*(4*slope*w.i/w.x.i*(u.i^2/w.y.i - v.i^2/w.x.i) - 
		2*v.i*HH/w.x.i + 2*u.i*JJ/w.y.i))
	CC <- -sum(w.i^2/w.y.i*(4*slope*w.i*u.i*v.i/w.x.i + v.i*JJ + u.i*HH))
	DDj <- w.i^2*v.i/w.x.i - w.i/w*sum(w.i^2*v.i/w.x.i)
	EEj <- 2*(w.i^2*u.i/w.y.i - w.i/w*sum(w.i^2*u.i/w.y.i))
	FFj <- w.i^2*v.i/w.y.i - w.i/w*sum(w.i^2*v.i/w.y.i)
	GGj <- w.i^2*u.i/w.x.i - w.i/w*sum(w.i^2*u.i/w.x.i)
	A <- sum(w.i^2*u.i*v.i/w.x.i)
	B <- sum(w.i^2*(u.i^2/w.y.i - v.i^2/w.x.i))

	dm.dXj <- rep(0, length(old.x))
	dm.dYj <- rep(0, length(old.x))
	dc.dXj <- rep(0, length(old.x))
	dc.dYj <- rep(0, length(old.x))
	dm.dXj[subset] <- -(slope^2*DDj + slope*EEj - FFj)/(2*slope*A + B - 
		AA*slope^2 + BB*slope - CC)
	dm.dYj[subset] <- -(slope^2*GGj - 2*slope*DDj - EEj/2)/(2*slope*A + B - 
		AA*slope^2 + BB*slope - CC)
	dc.dXj[subset] <- (HH - slope*JJ - x.bar)*dm.dXj[subset] - slope*w.i/w
	dc.dYj[subset] <- (HH - slope*JJ - x.bar)*dm.dYj[subset] + w.i/w

	se.slope.mR <- sqrt(sum(dm.dYj[subset]^2/w.y.i + dm.dXj[subset]^2/w.x.i))
	se.intercept.mR <- sqrt(sum(dc.dYj[subset]^2/w.y.i + 
		dc.dXj[subset]^2/w.x.i))
	cov.slope.intercept.mR <- sum(dc.dXj[subset]*dm.dXj[subset]/w.x.i + 
		dc.dYj[subset]*dm.dYj[subset]/w.y.i)

#	var.slope.Reed<- sum(w.i*(slope*U.i-V.i)^2)/(length(w.X.i)-2)*sum(dm.dyj^2/w.Y.i+dm.dxj^2/w.X.i)
#	var.intercept.Reed<- sum(w.i*(slope*U.i-V.i)^2)/(length(w.X.i)-2)*sum(dc.dyj^2/w.Y.i+dc.dxj^2/w.X.i)
#	se.slope<- sqrt(var.slope.Reed)
#	se.intercept<- sqrt(var.intercept.Reed)

	# uncertainty estimation using Williamson.1968 results, 
	# NOT verified algebraically but numerically equivalent to the Reed solution in the test.
	Z.i <- w.i*(u.i/w.y.i + slope*v.i/w.x.i)
	Z.bar <- sum(w.i*Z.i)/sum(w.i)
	Q.inv <- sum(w.i*(u.i*v.i/slope + 4*(Z.i - Z.bar)*(Z.i - u.i)))
	var.slope.Williamson <- 1/Q.inv^2*sum(w.i^2*(u.i^2/w.y.i + v.i^2/w.x.i))
	var.intercept.Williamson <- 1/(sum(w.i)) + 
		2*(x.bar + 2*Z.bar)*Z.bar/Q.inv + 
		(x.bar + 2*Z.bar)^2*var.slope.Williamson
	se.slope.W <- sqrt(var.slope.Williamson)
	se.intercept.W <- sqrt(var.intercept.Williamson)

#	if (abs(se.slope.W-se.slope.mR)/(se.slope.W+se.slope.mR)>max.tol) warning(paste("se{Slope}: Williamson solution(",se.slope.W,") and Reed solution(",se.slope.mR,") difference exceeded max.tolerance", max.tol))
#	if (abs(se.intercept.W-se.intercept.mR)/(se.intercept.W+se.intercept.mR)>max.tol) warning(paste("se{Intercept}: Williamson solution(",se.slope.W,") and Reed solution(",se.slope.mR,") difference exceeded max.tolerance", max.tol))

	dof <- length(x) - 2

	# the respective estimate (x.hat, y.hat) of the expected values (x.i=E[X.i], y.i=E[Y.i]) are:
	x.hat <- (old.x + old.lambda*slope*(old.y - intercept))/
		(1 + old.lambda*slope^2)
	y.hat <- intercept + slope*x.hat

	x.bar <- sum(w.i*x.hat[subset])/sum(w.i)
	y.bar <- sum(w.i*y.hat[subset])/sum(w.i)

	u.i <- x.hat[subset] - x.bar
	v.i <- y.hat[subset] - y.bar

	var.beta.mle <- 1/sum(w.i*u.i^2)
	var.alpha.mle <- 1/sum(w.i) + x.bar^2*var.beta.mle
	cov.alpha.beta.mle <- -x.bar*var.beta.mle
	
	# build coefficients table
	coef<- matrix(c(intercept, sqrt(var.alpha.mle), 
		(intercept - coef.H0[1])/sqrt(var.alpha.mle), 
		2*(1 - pt(abs(intercept - coef.H0[1])/sqrt(var.alpha.mle), dof)), 
		slope, sqrt(var.beta.mle), (slope - coef.H0[2])/sqrt(var.beta.mle), 
		2*(1 - pt(abs(slope - coef.H0[2])/sqrt(var.beta.mle), dof))), 2, 4,
		byrow = TRUE)
	rownames(coef) <- c("intercept", "slope")
	colnames(coef) <- c("Estimate", "Std. error", "t.value(H0)", 
		"Pr(t.value>|t| | H0)")

	# degrees of equivalence estimation for all the observed points 
	# using only the subset of observations to estimate the intercept and slope
	DoE <- (old.x - x.hat)*sqrt(1 + 1/(old.lambda^2*slope^4))

	# first we obtain the uncertainty for the deleted observations
	# DoE.i(X.i, Y.i, alpha(i), beta(i))
	dDoE.i.dx <- old.lambda*slope^2/(1 + old.lambda*slope^2)*
		sqrt(1 + 1/old.lambda^2/slope^4)
	dDoE.i.dy <- -1/slope*dDoE.i.dx
	dDoE.i.da <- -dDoE.i.dy
	dDoE.i.db <- (old.lambda*(intercept + 2*slope*old.x - old.y))/
		(1 + old.lambda*slope^2 - 2*old.lambda^2*slope^2*
			(intercept + slope*old.x - old.y)/(1 + old.lambda*slope^2)^2)*
			sqrt(1 + 1/old.lambda^2/slope^4) - 
			2*(intercept + slope*old.x - old.y)/((1 + old.lambda*slope^2)*
			slope^2*sqrt(1 + old.lambda^2*slope^4))

	var.DoE.i.mle <- dDoE.i.dx^2/old.w.x.i + dDoE.i.dy^2/old.w.y.i + 
		dDoE.i.da^2*var.alpha.mle + dDoE.i.db^2*var.beta.mle + 
		2*dDoE.i.da*dDoE.i.db*cov.alpha.beta.mle
	u.DoE.i.mle <- sqrt(qchisq(1-alpha, 2))/2*sqrt(var.DoE.i.mle)

	var.DoE.i.lse <- dDoE.i.dx^2/old.w.x.i + dDoE.i.dy^2/old.w.y.i + 
		dDoE.i.da^2*se.intercept.mR^2 + dDoE.i.db^2*se.slope.mR^2 + 
		2*dDoE.i.da*dDoE.i.db*cov.slope.intercept.mR
	var.DoE.i.lse[var.DoE.i.lse < 0] <- 0
	u.DoE.i.lse <- sqrt(qchisq(1-alpha, 2))/2*sqrt(var.DoE.i.lse)


	# second we obtain the uncertainties for the non-deleted observations
	# the uncertainty for the included observations needs to take additional dependences into account.
	# DoE({(X,Y)}, alpha({(X,Y)}), beta({(X,Y)}))
	dDoE.dXj <- diag(sqrt(1 + lambda^2*slope^4)/(1 + lambda*slope^2)) +
		((lambda*(intercept + 2*slope*x - y)/(1 + lambda*slope^2) - 
			2*(lambda^2*slope^2*(intercept + slope*x - y))/
			(1 + lambda*slope^2)^2)*sqrt(1 + 1/lambda^2/slope^4) - 
			2*(intercept + slope*x - y)/(1 + lambda*slope^2)/
			slope^2/sqrt(1 + lambda^2*slope^4))%*%t(dm.dXj[subset]) +
		(lambda*slope/(1 + lambda*slope^2)*sqrt(1 + 1/lambda^2/slope^4))%*%
			t(dc.dXj[subset])
	dDoE.dYj <- -diag(sqrt(1 + lambda^2*slope^4)/(1 + lambda*slope^2)) +
		((lambda*(intercept + 2*slope*x - y)/(1 + lambda*slope^2) - 
			2*(lambda^2*slope^2*(intercept + slope*x - y))/
			(1 + lambda*slope^2)^2)*sqrt(1 + 1/lambda^2/slope^4) - 
			2*(intercept + slope*x - y)/(1 + lambda*slope^2)/
			slope^2/sqrt(1 + lambda^2*slope^4))%*%t(dm.dYj[subset]) +
		(lambda*slope/(1 + lambda*slope^2)*sqrt(1 + 1/lambda^2/slope^4))%*%
			t(dc.dYj[subset])
	var.DoE <- dDoE.dXj^2%*%(1/w.x.i)+dDoE.dYj^2%*%(1/w.y.i)
	u.DoE.mle <- rep(0, length(old.x))
	u.DoE.mle[subset] <- sqrt(var.DoE)
	u.DoE.mle[!subset] <- u.DoE.i.mle[!subset]

	u.DoE.lse <- rep(0, length(old.x))
	u.DoE.lse[subset] <- sqrt(var.DoE)
	u.DoE.lse[!subset] <- u.DoE.i.lse[!subset]

	res <- list(X = old.x, Y = old.y, weights.x = old.w.x.i, 
			weights.y = old.w.y.i, subset = subset, coef.H0 = coef.H0,
			coefficients = coef, 
			cov.mle = matrix(c(var.alpha.mle, cov.alpha.beta.mle, 
				cov.alpha.beta.mle, var.beta.mle), 2, 2),
			cov.lse = matrix(c(se.intercept.mR^2, cov.slope.intercept.mR, 
				cov.slope.intercept.mR, se.slope.mR^2), 2, 2),
			x.hat = x.hat, y.hat = y.hat, 
			residuals = old.y - intercept - slope*old.x,
			df.residuals = dof, 
			MSE = sum((y - intercept-slope*x)^2)/dof,
			DoE = DoE, 
			u.DoE.mle = u.DoE.mle,
			u.DoE.lse = u.DoE.lse,
			dm.dXj = dm.dXj,
			dm.dYj = dm.dYj,
			dc.dXj = dc.dXj,
			dc.dYj = dc.dYj,
			curr.iter = curr.iter,
			curr.tol = curr.tol )

	class(res) <- "dwlm"
	if (trace) print(log)
	return( res )	
}
