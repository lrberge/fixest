# load("data/trade.RData") ; load("data/base_did.RData")
# roxygen2::roxygenise(roclets = "rd")

femlm_only_clusters <- function(env){
	# Estimation with only the cluster coefficients

	#
	# 1st step => computing the dummies
	#

	res = get("res", env)
	nobs = get("nobs", env)
	family = get("family", env)
	offset.value = get("offset.value", env)
	model0 = get("model0", env)

	if(family == "negbin"){
		coef[[".theta"]] = model0$theta
	} else {
		coef = list()
	}

	# indicator of whether we compute the exp(mu)
	useExp = family %in% c("poisson", "logit", "negbin")
	useExp_clusterCoef = family %in% c("poisson")

	# mu, using the offset
	mu_noDum = offset.value
	if(length(mu_noDum) == 1) mu_noDum = rep(mu_noDum, nobs)

	# we create the exp of mu => used for later functions
	exp_mu_noDum = NULL
	if(useExp_clusterCoef){
		exp_mu_noDum = rpar_exp(mu_noDum, env)
	}

	dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef)

	exp_mu = NULL
	if(useExp_clusterCoef){
		# despite being called mu, it is in fact exp(mu)!!!
		exp_mu = exp_mu_noDum*dummies
		mu = rpar_log(exp_mu, env)
	} else {
		mu = mu_noDum + dummies
		if(useExp){
			exp_mu = rpar_exp(mu, env)
		}
	}

	#
	# 2nd step => saving information
	#

	fixef_sizes = get("fixef_sizes", env)
	famFuns = get("famFuns", env)
	lhs = get("lhs", env)

	# The log likelihoods
	loglik = famFuns$ll(lhs, mu, exp_mu, env, coef)
	ll_null = model0$loglik

	# degres de liberte
	df_k = sum(fixef_sizes - 1) + 1
	pseudo_r2 = 1 - (loglik - df_k)/(ll_null - 1)

	# Calcul residus
	expected.predictor = famFuns$expected.predictor(mu, exp_mu, env)
	residuals = lhs - expected.predictor

	# calcul squared corr
	if(sd(expected.predictor) == 0){
		res$sq.cor = NA
	} else {
		res$sq.cor = stats::cor(lhs, expected.predictor)**2
	}

	# calcul r2 naif
	res$naive.r2 = 1 - sum(residuals**2) / sum((lhs - mean(lhs))**2)

	ssr_null = cpp_ssr_null(lhs)

	res$loglik = loglik
	res$n = length(lhs)
	res$nparams = df_k
	res$ll_null = ll_null
	res$ssr_null = ssr_null
	res$pseudo_r2 = pseudo_r2
	res$expected.predictor = expected.predictor
	res$residuals = residuals
	res$family = family
	res$method = "femlm"

	#
	# Information on the dummies

	if(useExp_clusterCoef){
		dummies = rpar_log(dummies, env)
	}

	res$sumFE = dummies

	if(family == "negbin"){
		theta = coef[".theta"]
		res$theta = theta
	}

	res$convStatus = TRUE

	class(res) = "fixest"

	return(res)
}

femlm_hessian <- function(coef, env){
	# Computes the hessian

	verbose = get("verbose", env)
	if(verbose >= 2) ptm = proc.time()
	params <- get("params", env)
	names(coef) <- params
	nonlinear.params <- get("nonlinear.params", env)
	k <- length(nonlinear.params)
	isNL <- get("isNL", env)
	hessian.args = get("hessian.args", env)
	famFuns = get("famFuns", env)
	family = get("family", env)
	y = get("lhs", env)
	isFixef = get("isFixef", env)
	nthreads = get("nthreads", env)
	mu_both = get_savedMu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	jacob.mat = get_Jacobian(coef, env)

	ll_d2 = famFuns$ll_d2(y, mu, exp_mu, coef)
	if(isFixef){
		dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
		jacob.mat = jacob.mat + dxi_dbeta
	} else dxi_dbeta = 0

	# hessVar = crossprod(jacob.mat, jacob.mat * ll_d2)
	hessVar = cpppar_crossprod(jacob.mat, ll_d2, nthreads)

	if(isNL){
		# we get the 2nd derivatives
		z = numDeriv::genD(evalNLpart, coef[nonlinear.params], env=env, method.args = hessian.args)$D[, -(1:k), drop=FALSE]
		ll_dl = famFuns$ll_dl(y, mu, exp_mu, coef=coef, env=env)
		id_r = rep(1:k, 1:k)
		id_c = c(sapply(1:k, function(x) 1:x), recursive=TRUE)
		H = matrix(0, nrow=k, ncol=k)
		H[cbind(id_r, id_c)] = H[cbind(id_r, id_c)] = colSums(z*ll_dl)
	} else H = 0

	# on ajoute la partie manquante
	if(isNL) hessVar[1:k, 1:k] = hessVar[1:k, 1:k] + H

	if(family == "negbin"){
		theta = coef[".theta"]
		ll_dx_dother = famFuns$ll_dx_dother(y, mu, exp_mu, coef, env)

		if(isFixef){
			dxi_dother = deriv_xi_other(ll_dx_dother, ll_d2, env, coef)
		} else {
			dxi_dother = 0
		}

		# calcul des derivees secondes vav de theta
		h.theta.L = famFuns$hess.thetaL(theta, jacob.mat, y, dxi_dbeta, dxi_dother, ll_d2, ll_dx_dother)
		hessVar = cbind(hessVar, h.theta.L)
		h.theta = famFuns$hess_theta_part(theta, y, mu, exp_mu, dxi_dother, ll_dx_dother, ll_d2, env)
		hessVar = rbind(hessVar, c(h.theta.L, h.theta))

	}

	if(anyNA(hessVar)){
		stop("NaN in the Hessian, can be due to a possible overfitting problem.\nIf so, to have an idea of what's going on, you can reduce the value of the argument 'rel.tol' of the nlminb algorithm using the argument 'opt.control = list(rel.tol=?)' with ? the new value.")
	}

	# warn_0_Hessian = get("warn_0_Hessian", env)
	# if(!warn_0_Hessian && any(diag(hessVar) == 0)){
	# 	# We apply the warning only once
	# 	var_problem = params[diag(hessVar) == 0]
	# 	warning("Some elements of the diagonal of the hessian are equal to 0: likely presence of collinearity. FYI the problematic variables are: ", paste0(var_problem, collapse = ", "), ".", immediate. = TRUE)
	# 	assign("warn_0_Hessian", TRUE, env)
	# }

	# if(verbose >= 2) cat("Hessian: ", (proc.time()-ptm)[3], "s\n", sep="")
	- hessVar
}

femlm_gradient <- function(coef, env){
	# cat("gradient:\n") ; print(as.vector(coef))

	params = get("params", env)
	names(coef) = params
	nonlinear.params = get("nonlinear.params", env)
	linear.params = get("linear.params", env)
	famFuns = get("famFuns", env)
	family = get("family", env)
	y = get("lhs", env)
	mu_both = get_savedMu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# calcul de la jacobienne
	res <- list() #stocks the results

	# cat("\tgetting jacobian")
	# ptm = proc.time()
	jacob.mat = get_Jacobian(coef, env)
	# cat("in", (proc.time()-ptm)[3], "s.\n")

	# cat("\tComputing gradient ")
	# ptm = proc.time()
	# res = famFuns$grad(jacob.mat, y, mu, env, coef)
	res = getGradient(jacob.mat, y, mu, exp_mu, env, coef)
	# cat("in", (proc.time()-ptm)[3], "s.\n")
	names(res) = c(nonlinear.params, linear.params)

	if(family=="negbin"){
		theta = coef[".theta"]
		res[".theta"] = famFuns$grad.theta(theta, y, mu, exp_mu, env)
	}

	return(-unlist(res[params]))
}

femlm_scores <- function(coef, env){
	# Computes the scores (Jacobian)
	params = get("params", env)
	names(coef) <- params
	famFuns = get("famFuns", env)
	family = get("family", env)
	y = get("lhs", env)
	mu_both = get_savedMu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	jacob.mat = get_Jacobian(coef, env)
	scores = getScores(jacob.mat, y, mu, exp_mu, env, coef)

	if(family=="negbin"){
		theta = coef[".theta"]
		score.theta = famFuns$scores.theta(theta, y, mu, exp_mu, env)
		scores = cbind(scores, score.theta)
	}

	return(scores)
}

femlm_ll <- function(coef, env){
	# Log likelihood

	# misc funs
	iter = get("iter", env) + 1
	assign("iter", iter, env)
	pastLL = get("pastLL", env)
	verbose = get("verbose", env)
	ptm = proc.time()
	if(verbose >= 1){
		coef_names = sapply(names(coef), charShorten, width = 10)
		coef_line = paste0(coef_names, ": ", signif(coef), collapse = " -- ")
		cat("\nIter", iter, "- Coefficients:", coef_line, "\n")
	}

	# we save the coefs => for debugging if optimization fails
	assign("coef_evaluated", coef, env)

	# computing the LL
	famFuns = get("famFuns", env)
	family = get("family", env)
	y <- get("lhs", env)

	if(anyNA(coef)) stop("Divergence... (some coefs are NA)\nTry option verbose=2 to figure out the problem.")

	mu_both = get_mu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# for the NEGBIN, we add the coef
	ll = famFuns$ll(y, mu, exp_mu, env, coef)

	if(!is.finite(ll)){
		stop("Divergence... (evaluation of the log-likelihood is not finite)\nTry option verbose=2 to figure out the problem.")
	}

	# We save the ll and ssr with FEs of the first iteration
	isFixef = get("isFixef", env)
	if(iter == 1 && isFixef && all(head(coef, length(coef) - (family == "negbin")) == 0)){
		assign("ll_fe_only", ll, env)
		ep = famFuns$expected.predictor(mu, exp_mu)
		# assign("ssr_fe_only", drop(crossprod(y - ep)), env)
		assign("ssr_fe_only", cpp_ssq(y - ep), env)
	}

	evolutionLL = ll - pastLL
	assign("evolutionLL", evolutionLL, env)
	assign("pastLL", ll, env)

	if(iter == 1) evolutionLL = "--"
	if(verbose >= 1) cat("LL = ", ll, " (", (proc.time()-ptm)[3], "s)\tevol = ", evolutionLL, "\n", sep = "")

	# To cope with difficult convergence
	# if(iter > 20 && ll > pastLL && evolutionLL / (0.1 + abs(ll)) < 1e-8){
	# 	print( evolutionLL / (0.1 + ll))
	# 	# Convergence
	# 	assign("coef_evaluated", coef, env)
	# 	assign("convergence_ll", TRUE, env)
	#
	# 	stop("Convergence.")
	# }

	return(-ll) # je retourne -ll car la fonction d'optimisation minimise
}

evalNLpart = function(coef, env){
	# cat("Enter evalNLpart : ", as.vector(coef), "\n")
	# fonction qui evalue la partie NL
	isNL = get("isNL", env)
	if(!isNL) return(0)

	envNL = get("envNL", env)
	nonlinear.params <- get("nonlinear.params", env)
	nl.call <- get("nl.call", env)
	nbSave = get("nbSave", env)
	nbMaxSave = get("nbMaxSave", env)
	savedCoef = get("savedCoef", env)
	savedValue = get("savedValue", env)

	if(!is.null(names(coef))){
		coef = coef[nonlinear.params]
	} else if (length(coef) != length(nonlinear.params)){
		stop("Problem with the length of the NL coefficients.")
	}

	if(nbMaxSave == 0){
		for(var in nonlinear.params) assign(var, coef[var], envNL)
		y_nl <- eval(nl.call, envir = envNL)

		# we check problems
		if(anyNA(y_nl)){
			stop("Evaluation of non-linear part returns NAs. The coefficients were: ", paste0(nonlinear.params, " = ", signif(coef[nonlinear.params], 3)), ".")
		}

		return(y_nl)
	}

	for(i in nbSave:1){
		#les valeurs les + recentes sont en derniere position
		if(all(coef == savedCoef[[i]])){
			return(savedValue[[i]])
		}
	}

	# Si la valeur n'existe pas, on la sauvegarde
	# on met les valeurs les plus recentes en derniere position
	for(var in nonlinear.params) assign(var, coef[var], envNL)
	y_nl = eval(nl.call, envir = envNL)

	# we check problems
	if(anyNA(y_nl)){
		stop("Evaluation of non-linear part returns NAs. The coefficients were: ", paste0(nonlinear.params, " = ", signif(coef[nonlinear.params], 3)), ".")
	}

	if(nbSave < nbMaxSave){
		savedCoef[[nbSave + 1]] = coef
		savedValue[[nbSave + 1]] = y_nl
		assign("nbSave", nbSave + 1, env)
	} else if(nbMaxSave > 1){
		tmp = list()
		tmp[[nbSave]] = coef
		tmp[1:(nbSave-1)] = savedCoef[2:nbSave]
		savedCoef = tmp

		tmp = list()
		tmp[[nbSave]] = y_nl
		tmp[1:(nbSave-1)] = savedValue[2:nbSave]
		savedValue = tmp
	} else{
		savedCoef = list(coef)
		savedValue = list(y_nl)
	}

	# cat("computed NL part:", as.vector(coef), "\n")

	assign("savedCoef", savedCoef, env)
	assign("savedValue", savedValue, env)
	return(y_nl)
}

get_mu = function(coef, env, final = FALSE){
	# This function computes the RHS of the equation
	# mu_L => to save one matrix multiplication
	isNL = get("isNL", env)
	isLinear = get("isLinear", env)
	isFixef = get("isFixef", env)
	nobs = get("nobs", env)
	params = get("params", env)
	family = get("family", env)
	offset.value = get("offset.value", env)
	names(coef) = params

	# UseExp: indicator if the family needs to use exp(mu) in the likelihoods:
	#     this is useful because we have to compute it only once (save computing time)
	# useExp_clusterCoef: indicator if we use the exponential of mu to obtain the cluster coefficients
	#     if it is TRUE, it will mean that the dummy will be equal
	#     to exp(sumFE) despite being named sumFE
	useExp = family %in% c("poisson", "logit", "negbin")
	useExp_clusterCoef = family %in% c("poisson")

	# For managing mu:
	coefMu = get("coefMu", env)
	valueMu = get("valueMu", env)
	valueExpMu = get("valueExpMu", env)
	wasUsed = get("wasUsed", env)
	if(wasUsed){
		coefMu = valueMu = valueExpMu = list()
		assign("wasUsed", FALSE, env)
	}

	if(length(coefMu)>0){
		for(i in 1:length(coefMu)){
			if(all(coef==coefMu[[i]])){
				return(list(mu = valueMu[[i]], exp_mu = valueExpMu[[i]]))
			}
		}
	}

	if(isNL){
		muNL = evalNLpart(coef, env)
	} else muNL = 0

	if(isLinear){
		linear.params = get("linear.params", env)
		linear.mat = get("linear.mat", env)
		nthreads = get("nthreads", env)
		# mu_L = c(linear.mat %*% coef[linear.params])
		mu_L = cpppar_xbeta(linear.mat, coef[linear.params], nthreads)
	} else mu_L = 0

	mu_noDum = muNL + mu_L + offset.value

	# Detection of overfitting issues with the logit model:
	if(family == "logit" && FALSE){
		warn_overfit_logit = get("warn_overfit_logit", env)

		if(!warn_overfit_logit && max(abs(mu_noDum)) >= 100){
			# overfitting => now finding the precise cause
			if(!isNL || (isLinear && max(abs(mu_L)) >= 100)){
				# we create the matrix with the coefficients to find out the guy
				mat_L_coef = linear.mat * matrix(coef[linear.params], nrow(linear.mat), length(linear.params), byrow = TRUE)
				max_var = apply(abs(mat_L_coef), 2, max)
				best_suspect = linear.params[which.max(max_var)]
				stop("in femlm(): Likely presence of a separation problem (coef->Inf). One suspect variable is: ", best_suspect, ".", immediate. = TRUE, call. = FALSE)
			} else {
				warning("in femlm(): Likely presence of an overfitting problem due to the non-linear part.", immediate. = TRUE, call. = FALSE)
			}

			assign("warn_overfit_logit", TRUE, env)
		}
	}


	# we create the exp of mu => used for later functions
	exp_mu_noDum = NULL
	if(useExp_clusterCoef){
		exp_mu_noDum = rpar_exp(mu_noDum, env)
	}

	if(isFixef){
		# we get back the last dummy
		sumFE = getDummies(mu_noDum, exp_mu_noDum, env, coef, final)
	} else {
		if(useExp_clusterCoef){
			sumFE = 1
		} else {
			sumFE = 0
		}
	}

	# We add the value of the dummy to mu and we compute the exp if necessary
	exp_mu = NULL
	if(useExp_clusterCoef){
		# despite being called sumFE, it is in fact exp(sumFE)!!!
		exp_mu = exp_mu_noDum*sumFE
		mu = rpar_log(exp_mu, env)
	} else {
		mu = mu_noDum + sumFE
		if(useExp){
			exp_mu = rpar_exp(mu, env)
		}
	}

	if(isFixef){
		# BEWARE, if useExp_clusterCoef, it is equal to exp(sumFE)
		attr(mu, "sumFE") = sumFE
	}

	if(length(mu)==0) mu = rep(mu, nobs)

	# we save the value of mu:
	coefMu = append(coefMu, list(coef))
	valueMu = append(valueMu, list(mu))
	valueExpMu = append(valueExpMu, list(exp_mu))
	assign("coefMu", coefMu, env)
	assign("valueMu", valueMu, env)
	assign("valueExpMu", valueExpMu, env)

	return(list(mu = mu, exp_mu = exp_mu))
}

get_savedMu = function(coef, env){
	# This function gets the mu without computation
	# It follows a LL evaluation
	coefMu = get("coefMu", env)
	valueMu = get("valueMu", env)
	valueExpMu = get("valueExpMu", env)
	assign("wasUsed", TRUE, env)

	if(length(coefMu)>0) for(i in 1:length(coefMu)) if(all(coef==coefMu[[i]])){
		# cat("coef nb:", i, "\n")
		return(list(mu = valueMu[[i]], exp_mu = valueExpMu[[i]]))
	}

	stop("Problem in \"get_savedMu\":\n gradient did not follow LL evaluation.")
}

get_Jacobian = function(coef, env){
	# retrieves the Jacobian of the "rhs"
	params <- get("params", env)
	names(coef) <- params
	isNL <- get("isNL", env)
	isLinear <- get("isLinear", env)
	isGradient = get("isGradient", env)

	if(isNL){
		nonlinear.params = get("nonlinear.params", env)
		jacob.mat = get_NL_Jacobian(coef[nonlinear.params], env)
	} else jacob.mat = c()

	if(isLinear){
		linear.mat = get("linear.mat", env)
		if(is.null(dim(jacob.mat))){
			jacob.mat = linear.mat
		} else {
			jacob.mat = cbind(jacob.mat, linear.mat)
		}
	}

	return(jacob.mat)
}

get_NL_Jacobian = function(coef, env){
	# retrieves the Jacobian of the non linear part
	#cat("In NL JAC:\n")
	#print(coef)
	nbSave = get("JC_nbSave", env)
	nbMaxSave = get("JC_nbMaxSave", env)
	savedCoef = get("JC_savedCoef", env)
	savedValue = get("JC_savedValue", env)

	nonlinear.params <- get("nonlinear.params", env)
	coef = coef[nonlinear.params]

	if(nbSave>0) for(i in nbSave:1){
		#les valeurs les + recentes sont en derniere position
		if(all(coef == savedCoef[[i]])){
			# cat("Saved value:", as.vector(coef), "\n")
			return(savedValue[[i]])
		}
	}

	#Si la valeur n'existe pas, on la sauvegarde
	#on met les valeurs les plus recentes en derniere position
	isGradient <- get("isGradient", env)
	if(isGradient){
		call_gradient <- get("call_gradient", env)
		#we send the coef in the environment
		for(var in nonlinear.params) assign(var, coef[var], env)
		jacob.mat <- eval(call_gradient, envir=env)
		jacob.mat <- as.matrix(as.data.frame(jacob.mat[nonlinear.params]))
	} else {
		jacobian.method <- get("jacobian.method", env)
		jacob.mat <- numDeriv::jacobian(evalNLpart, coef, env=env, method=jacobian.method)
	}

	#Controls:
	if(anyNA(jacob.mat)){
		qui <- which(apply(jacob.mat, 2, function(x) anyNA(x)))
		variables <- nonlinear.params[qui]
		stop("ERROR: The Jacobian of the nonlinear part has NA!\nThis concerns the following variables:\n", paste(variables, sep=" ; "))
	}

	#Sauvegarde
	if(nbSave<nbMaxSave){
		savedCoef[[nbSave+1]] = coef
		savedValue[[nbSave+1]] = jacob.mat
		assign("JC_nbSave", nbSave+1, env)
	} else if(nbMaxSave>1){
		tmp = list()
		tmp[[nbSave]] = coef
		tmp[1:(nbSave-1)] = savedCoef[2:nbSave]
		savedCoef = tmp

		tmp = list()
		tmp[[nbSave]] = jacob.mat
		tmp[1:(nbSave-1)] = savedValue[2:nbSave]
		savedValue = tmp
	} else{
		savedCoef = list(coef)
		savedValue = list(jacob.mat)
	}

	# print(colSums(jacob.mat))

	# cat("computed NL Jacobian:", as.vector(coef), "\n")
	# print(savedCoef)

	assign("JC_savedCoef", savedCoef, env)
	assign("JC_savedValue", savedValue, env)
	return(jacob.mat)
}

get_model_null <- function(env, theta.init){
	# I have the closed form of the ll0
	famFuns = get("famFuns", env)
	family = get("family", env)
	N = get("nobs", env)
	y = get("lhs", env)
	verbose = get("verbose", env)
	ptm = proc.time()

	# one of the elements to be returned
	theta = NULL

	if(family == "poisson"){
		# There is a closed form

		if("lfactorial" %in% names(env)){
			lfact = get("lfactorial", env)
		} else {
			lfact = sum(rpar_lgamma(y + 1, env))
			assign("lfactorial", lfact, env)
		}

	    w = env$weights.value
	    isWeight = length(w) > 1
	    if(isWeight){
	        sy = sum(y * w)
	        N = sum(w)
	    } else {
	        sy = sum(y)
	    }

	    constant = log(sy / N)
	    loglik =  sy*log(sy) - sy*log(N) - sy - lfact


	} else if(family == "gaussian"){
		# there is a closed form
		constant = mean(y)
		ss = sum( (y - constant)**2 )
		sigma = sqrt( ss / N )
		loglik = -1/2/sigma^2*ss - N*log(sigma) - N*log(2*pi)/2

	} else if(family == "logit"){
		# there is a closed form

	    w = env$weights.value
		isWeight = length(w) > 1
	    if(isWeight){
	        sy = sum(y * w)
	        N = sum(w)
	    } else {
	        sy = sum(y)
	    }

		constant = log(sy) - log(N - sy)
		loglik = sy*log(sy) - sy*log(N-sy) - N*log(N) + N*log(N-sy)

	} else if(family=="negbin"){

		if("lgamma" %in% names(env)){
			lgamm = get("lgamma", env)
		} else {
			lgamm = sum(rpar_lgamma(y + 1, env))
			assign("lgamma", lgamm, env)
		}

		sy = sum(y)
		constant = log(sy / N)

		mean_y = mean(y)
		invariant = sum(y*constant) - lgamm

		if(is.null(theta.init)){
			theta.guess = max(mean_y**2 / max((var(y) - mean_y), 1e-4), 0.05)
		} else {
			theta.guess = theta.init
		}

		# I set up a limit of 0.05, because when it is too close to 0, convergence isnt great

		opt <- nlminb(start=theta.guess, objective=famFuns$ll0_theta, y=y, gradient=famFuns$grad0_theta, lower=1e-3, mean_y=mean_y, invariant=invariant, hessian = famFuns$hess0_theta, env=env)

		loglik = -opt$objective
		theta = opt$par
	}

	if(verbose >= 2) cat("Null model in ", (proc.time()-ptm)[3], "s. ", sep ="")

	return(list(loglik=loglik, constant=constant, theta = theta))
}

getGradient = function(jacob.mat, y, mu, exp_mu, env, coef, ...){
	famFuns = get("famFuns", env)
	nthreads = get("nthreads", env)
	ll_dl = famFuns$ll_dl(y, mu, exp_mu, coef=coef, env=env)

	cpppar_xwy(jacob.mat, ll_dl, 1, nthreads)
}

getScores = function(jacob.mat, y, mu, exp_mu, env, coef, ...){
	famFuns = get("famFuns", env)
	isFixef = get("isFixef", env)

	ll_dl = famFuns$ll_dl(y, mu, exp_mu, coef=coef, env=env)
	scores = jacob.mat* ll_dl

	if(isFixef){
		ll_d2 = famFuns$ll_d2(y, mu, exp_mu, coef=coef, env=env)
		dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
		scores = scores + dxi_dbeta * ll_dl
	}

	return(as.matrix(scores))
}

getDummies = function(mu, exp_mu, env, coef, final = FALSE){
	# function built to get all the dummy variables
	# We retrieve past dummies (that are likely to be good
	# starting values)

	sumFE = get("saved_sumFE", env)
	family = get("family", env)
	fixef.tol = get("fixef.tol", env)
	verbose = get("verbose", env)
	if(verbose >= 2) ptm = proc.time()

	#
	# Dynamic precision
	#

	iterCluster = get("iterCluster", env)
	evolutionLL = get("evolutionLL", env)
	nobs = get("nobs", env)
	iter = get("iter", env)
	iterLastPrecisionIncrease = get("iterLastPrecisionIncrease", env)

	nbIterOne = get("nbIterOne", env)
	if(iterCluster <= 2){
		nbIterOne = nbIterOne + 1
	} else { # we reinitialise
		nbIterOne = 0
	}
	assign("nbIterOne", nbIterOne, env)

	# nber of times LL almost didn't increase
	nbLowIncrease = get("nbLowIncrease", env)
	if(evolutionLL/nobs < 1e-8){
		nbLowIncrease = nbLowIncrease + 1
	} else { # we reinitialise
		nbLowIncrease = 0
	}
	assign("nbLowIncrease", nbLowIncrease, env)

	if(!final && fixef.tol > .Machine$double.eps*10000 && iterCluster <= 2 && nbIterOne >= 2 && ((nbLowIncrease >= 2 && (iter - iterLastPrecisionIncrease) >= 3) || (iter - iterLastPrecisionIncrease) >= 15)){
		fixef.tol = fixef.tol/10
		if(verbose >= 2) cat("Precision increased to", fixef.tol, "\n")
		assign("fixef.tol", fixef.tol, env)
		assign("iterLastPrecisionIncrease", iter, env)

		# If the precision increases, we must also increase the precision of the dummies!
		if(family %in% c("negbin", "logit")){
			assign("NR.tol", fixef.tol / 100, env)
		}

		# we also set acceleration to on
		assign("useAcc", TRUE, env)
	} else if(final){
		# we don't need ultra precision for these last dummies
		fixef.tol = fixef.tol * 10**(iterLastPrecisionIncrease != 0)
		if(family %in% c("negbin", "logit")){
			assign("NR.tol", fixef.tol / 100, env)
		}
	}

	iterMax = get("fixef.iter", env)
	fixef_sizes = get("fixef_sizes", env)
	Q = length(fixef_sizes)

	# whether we use the eponentiation of mu
	useExp_clusterCoef = family %in% c("poisson")
	if(useExp_clusterCoef){
		mu_in = exp_mu * sumFE
	} else {
		mu_in = mu + sumFE
	}

	#
	# Computing the optimal mu
	#

	useAcc = get("useAcc", env)
	carryOn = FALSE

	# Finding the complexity of the problem

	firstRunCluster = get("firstRunCluster", env)
	if(firstRunCluster && Q >= 3){
		# First iteration: we check if the problem is VERY difficult (for Q = 3+)
		useAcc = TRUE
		assign("useAcc", TRUE, env)
		res = convergence(coef, mu_in, env, iterMax = 15)
		if(res$iter == 15){
			assign("difficultConvergence", TRUE, env)
			carryOn = TRUE
		}
	} else if(useAcc){
		res = convergence(coef, mu_in, env, iterMax)
		if(res$iter <= 2){
			# if almost no iteration => no acceleration next time
			assign("useAcc", FALSE, env)
		}
	} else {
		res = convergence(coef, mu_in, env, iterMax = 15)
		if(res$iter == 15){
			carryOn = TRUE
		}
	}

	if(carryOn){
		# the problem is difficult => acceleration on
		useAcc = TRUE
		assign("useAcc", TRUE, env)

		res = convergence(coef, res$mu_new, env, iterMax)
	}

	mu_new = res$mu_new
	iter = res$iter

	#
	# Retrieving the value of the dummies
	#

	if(useExp_clusterCoef){
		sumFE = mu_new / exp_mu
	} else {
		sumFE = mu_new - mu
	}

	# Warning messages if necessary:
	if(iter == iterMax) {
		fixef.iter.limit_reached = get("fixef.iter.limit_reached", env)
		fixef.iter.limit_reached = fixef.iter.limit_reached + 1
		assign("fixef.iter.limit_reached", fixef.iter.limit_reached, env)
	}

	assign("iterCluster", iter, env)

	# we save the dummy:
	assign("saved_sumFE", sumFE, env)

	if(verbose >= 2){
		acc_info = ifelse(useAcc, "+Acc. ", "-Acc. ")
		cat("Cluster Coef.: ", (proc.time()-ptm)[3], "s (", acc_info, "iter:", iter, ")\t", sep = "")
	}

	# we update the flag
	assign("firstRunCluster", FALSE, env)

	sumFE
}

deriv_xi = function(jacob.mat, ll_d2, env, coef){
	# Derivative of the cluster coefficients

	# data:
	iterMax = get("deriv.iter", env)
	fixef_sizes = get("fixef_sizes", env)
	Q = length(fixef_sizes)

	verbose = get("verbose", env)
	if(verbose >= 2) ptm = proc.time()

	#
	# initialisation of dxi_dbeta
	#

	if(Q >= 2){
		# We set the initial values for the first run
		if(!"sum_deriv" %in% names(env)){
			# init of the sum of the derivatives => 0
			dxi_dbeta = matrix(0, nrow(jacob.mat), ncol(jacob.mat))
		} else {
			dxi_dbeta = get("sum_deriv", env)
		}
	} else {
		# no need if only 1, direct solution
		dxi_dbeta = NULL
	}

	#
	# Computing the optimal dxi_dbeta
	#

	accDeriv = get("accDeriv", env)
	carryOn = FALSE

	# Finding the complexity of the problem

	firstRunDeriv = get("firstRunDeriv", env)
	if(firstRunDeriv){
		# set accDeriv: we use information on cluster deriv
		iterCluster = get("iterCluster", env)
		diffConv = get("difficultConvergence", env)
		if(iterCluster < 20 & !diffConv){
			accDeriv = FALSE
			assign("accDeriv", FALSE, env)
		}
	}

	if(firstRunDeriv && accDeriv && Q >= 3){
		# First iteration: we check if the problem is VERY difficult (for Q = 3+)
		assign("accDeriv", TRUE, env)
		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, iterMax = 15)
		if(res$iter == 15){
			assign("derivDifficultConvergence", TRUE, env)
			carryOn = TRUE
		}
	} else if(accDeriv){
		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, iterMax)
		if(res$iter <= 10){
			# if almost no iteration => no acceleration next time
			assign("accDeriv", FALSE, env)
		}
	} else {
		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, iterMax = 50)
		if(res$iter == 50){
			carryOn = TRUE
		}
	}

	if(carryOn){
		# the problem is difficult => acceleration on
		accDeriv = TRUE
		assign("accDeriv", TRUE, env)
		res = dconvergence(res$dxi_dbeta, jacob.mat, ll_d2, env, iterMax)
	}

	dxi_dbeta = res$dxi_dbeta
	iter = res$iter

	if(iter == iterMax){
		deriv.iter.limit_reached = get("deriv.iter.limit_reached", env)
		deriv.iter.limit_reached = deriv.iter.limit_reached + 1
		assign("deriv.iter.limit_reached", deriv.iter.limit_reached, env)
	}

	assign("firstRunDeriv", FALSE, env)
	assign("sum_deriv", dxi_dbeta, env)

	if(verbose >= 2){
		acc_info = ifelse(accDeriv, "+Acc. ", "-Acc. ")
		cat("  Derivatives: ", (proc.time()-ptm)[3], "s (", acc_info, "iter:", iter, ")\n", sep = "")
	}

	return(dxi_dbeta)
}

deriv_xi_other = function(ll_dx_dother, ll_d2, env, coef){
	# derivative of the dummies wrt an other parameter
	fixef_id.matrix_cpp = get("fixef_id.matrix_cpp", env)
	fixef_sizes = get("fixef_sizes", env)
	fixef_id = get("fixef_id_list", env)
	deriv.tol = get("deriv.tol", env)
	Q = length(fixef_id)
	iterMax = get("deriv.iter", env)

	if(Q == 1){
		dum = fixef_id[[1]]
		k = max(dum)
		S_Jmu = cpp_tapply_vsum(k, ll_dx_dother, dum)
		S_mu = cpp_tapply_vsum(k, ll_d2, dum)
		dxi_dother = - S_Jmu[dum] / S_mu[dum]
	} else {
		# The cpp way:

		N = length(ll_d2)

		# We set the initial values for the first run
		if(!"sum_deriv_other" %in% names(env)){
			init = rep(0, N)
		} else {
			init = get("sum_deriv_other", env)
		}

		dxi_dother <- cpp_partialDerivative_other(iterMax, Q, N, epsDeriv = deriv.tol, ll_d2, ll_dx_dother, init, fixef_id.matrix_cpp, fixef_sizes)

		# we save the values
		assign("sum_deriv_other", dxi_dother, env)

	}

	as.matrix(dxi_dother)
}

####
#### Convergence ####
####

convergence = function(coef, mu_in, env, iterMax){
	# computes the new mu wrt the cluster coefficients

	fixef_sizes = get("fixef_sizes", env)
	Q = length(fixef_sizes)
	useAcc = get("useAcc", env)
	diffConv = get("difficultConvergence", env)
	mem.clean = get("mem.clean", env)

	if(mem.clean){
	    gc()
	}

	if(useAcc && diffConv && Q > 2){
		# in case of complex cases: it's more efficient
		# to initialize the first two clusters

		res = conv_acc(coef, mu_in, env, iterMax, only2 = TRUE)
		mu_in = res$mu_new
	}

	if(Q == 1){
		mu_new = conv_single(coef, mu_in, env)
		iter = 1
	} else if(Q >= 2){
		# Dynamic setting of acceleration

		if(!useAcc){
			res = conv_seq(coef, mu_in, env, iterMax = iterMax)
		} else if(useAcc){
			res = conv_acc(coef, mu_in, env, iterMax = iterMax)
		}

		mu_new = res$mu_new
		iter = res$iter
	}

	# we return a list with: new mu and iterations
	list(mu_new = mu_new, iter = iter)
}

conv_single = function(coef, mu_in, env){
	# convergence for a single cluster
	# it returns: the new mu (NOT sumFE)

	# Loading all the required variables
	lhs = get("lhs", env)
	fixef_sizes = get("fixef_sizes", env)
	fixef_id_vector = get("fixef_id_vector", env)
	fixef_table_vector = get("fixef_table_vector", env)
	sum_y_vector = get("sum_y_vector", env)
	fixef_cumtable_vector = get("fixef_cumtable_vector", env)
	fixef_order_vector = get("fixef_order_vector", env)
	nthreads = get("nthreads", env)

	NR.tol = get("NR.tol", env)
	family = get("familyConv", env)

	family_nb = switch(family, poisson=1, negbin=2, logit=3, gaussian=4, lpoisson=5)
	theta = ifelse(family == "negbin", coef[".theta"], 1)

	mu_new = update_mu_single_cluster(family = family_nb, nb_cluster = fixef_sizes, theta = theta, diffMax_NR = NR.tol, mu_in = mu_in, lhs = lhs, sum_y = sum_y_vector, dum = fixef_id_vector, obsCluster = fixef_order_vector, table = fixef_table_vector, cumtable = fixef_cumtable_vector, nthreads = nthreads)

	return(mu_new)
}

conv_seq = function(coef, mu_in, env, iterMax){
	# convergence of cluster coef without acceleration
	# Now all in cpp

	# Loading all the required variables
	lhs = get("lhs", env)
	fixef_sizes = get("fixef_sizes", env)
	dum_vector = get("fixef_id_vector", env)
	fixef_table_vector = get("fixef_table_vector", env)
	sum_y_vector = get("sum_y_vector", env)
	fixef_cumtable_vector = get("fixef_cumtable_vector", env)
	fixef_order_vector = get("fixef_order_vector", env)
	nthreads = get("nthreads", env)

	fixef.tol = get("fixef.tol", env)
	NR.tol = get("NR.tol", env)
	family = get("familyConv", env)

	family_nb = switch(family, poisson=1, negbin=2, logit=3, gaussian=4, lpoisson=5)
	theta = ifelse(family == "negbin", coef[".theta"], 1)

	Q = length(fixef_sizes)

	if(family == "lpoisson"){
		# we transform the mu_in into a non exponential form
		mu_in = log(mu_in)
	}

	if(Q == 2 & family == "poisson"){
		# Required Variables
		setup_poisson_fixedcost(env)
		info = get("fixedCostPoisson", env)

		res = cpp_conv_seq_poi_2(n_i = info$n_i, n_j = info$n_j, n_cells = info$n_cells, index_i = info$index_i, index_j = info$index_j, order = info$order, dum_vector = dum_vector, sum_y_vector = sum_y_vector, iterMax = iterMax, diffMax = fixef.tol, exp_mu_in = mu_in)

	} else if(Q == 2 & family == "gaussian"){
		# Required variables
		setup_gaussian_fixedcost(env)
		info = get("fixedCostGaussian", env)
		invTableCluster_vector = get("fixef_invTable", env)

		res = cpp_conv_seq_gau_2(n_i = info$n_i, n_j = info$n_j, n_cells = info$n_cells, r_mat_row = info$mat_row, r_mat_col = info$mat_col, r_mat_value_Ab = info$mat_value_Ab, r_mat_value_Ba = info$mat_value_Ba, dum_vector = dum_vector, lhs = lhs, invTableCluster_vector = invTableCluster_vector, iterMax = iterMax, diffMax = fixef.tol, mu_in = mu_in)

	} else {
		res = cpp_conv_seq_gnl(family = family_nb, iterMax = iterMax, diffMax = fixef.tol, diffMax_NR = NR.tol, theta = theta, lhs = lhs, nb_cluster_all = fixef_sizes, mu_init = mu_in, dum_vector = dum_vector, tableCluster_vector = fixef_table_vector, sum_y_vector = sum_y_vector, cumtable_vector = fixef_cumtable_vector, obsCluster_vector = fixef_order_vector, nthreads = nthreads)
	}

	if(family == "lpoisson"){
		# we transform the mu_in into an exponential form
		res$mu_new = exp(res$mu_new)
	}

	return(res)
}

conv_acc = function(coef, mu_in, env, iterMax, only2 = FALSE){
	# convergence of cluster coef without acceleration
	# Now all in cpp

	# Loading all the required variables
	lhs = get("lhs", env)
	fixef_sizes = get("fixef_sizes", env)
	dum_vector = get("fixef_id_vector", env)
	fixef_table_vector = get("fixef_table_vector", env)
	sum_y_vector = get("sum_y_vector", env)
	fixef_cumtable_vector = get("fixef_cumtable_vector", env)
	fixef_order_vector = get("fixef_order_vector", env)
	nthreads = get("nthreads", env)

	fixef.tol = get("fixef.tol", env)
	NR.tol = get("NR.tol", env)
	family = get("familyConv", env)

	family_nb = switch(family, poisson=1, negbin=2, logit=3, gaussian=4, lpoisson=5)
	theta = ifelse(family == "negbin", coef[".theta"], 1)

	if(only2){
		# means we compute the CC of the first two FE
		# we recreate the values we send
		fixef_sizes = fixef_sizes[1:2]
		nb_keep = sum(fixef_sizes)
		fixef_table_vector = fixef_table_vector[1:nb_keep]
		sum_y_vector = sum_y_vector[1:nb_keep]
		fixef_cumtable_vector = fixef_cumtable_vector[1:nb_keep]

		dum_vector = dum_vector[1:(2*length(lhs))]
		fixef_order_vector = fixef_order_vector[1:(2*length(lhs))]
	}

	Q = length(fixef_sizes)

	if(family == "lpoisson"){
		# we transform the mu_in into a non exponential form
		mu_in = log(mu_in)
	}

	if(Q == 2 & family == "poisson"){
		# Required Variables
		setup_poisson_fixedcost(env)
		info = get("fixedCostPoisson", env)

		res = cpp_conv_acc_poi_2(n_i = info$n_i, n_j = info$n_j, n_cells = info$n_cells, index_i = info$index_i, index_j = info$index_j, order = info$order, dum_vector = dum_vector, sum_y_vector = sum_y_vector, iterMax = iterMax, diffMax = fixef.tol, exp_mu_in = mu_in)

	} else if(Q == 2 & family == "gaussian"){
		# Required variables
		setup_gaussian_fixedcost(env)
		info = get("fixedCostGaussian", env)
		invTableCluster_vector = get("fixef_invTable", env)

		res = cpp_conv_acc_gau_2(n_i = info$n_i, n_j = info$n_j, n_cells = info$n_cells, r_mat_row = info$mat_row, r_mat_col = info$mat_col, r_mat_value_Ab = info$mat_value_Ab, r_mat_value_Ba = info$mat_value_Ba, dum_vector = dum_vector, lhs = lhs, invTableCluster_vector = invTableCluster_vector, iterMax = iterMax, diffMax = fixef.tol, mu_in = mu_in)

	} else {
		res = cpp_conv_acc_gnl(family = family_nb, iterMax = iterMax, diffMax = fixef.tol, diffMax_NR = NR.tol, theta = theta, lhs = lhs, nb_cluster_all = fixef_sizes, mu_init = mu_in, dum_vector = dum_vector, tableCluster_vector = fixef_table_vector, sum_y_vector = sum_y_vector, cumtable_vector = fixef_cumtable_vector, obsCluster_vector = fixef_order_vector, nthreads = nthreads)
	}

	if(family == "poisson" && res$any_negative_poisson){
		# we need to switch to log poisson
		assign("familyConv", "lpoisson", env)

		verbose = get("verbose", env)
		if(verbose >= 3) cat("Switch to log-poisson (to cope with high valued FEs).\n")

		res = conv_acc(coef, mu_in, env, iterMax, only2)

		# we switch back to original poisson
		assign("familyConv", "poisson", env)
	}

	if(family == "lpoisson"){
		# we transform the mu_in into an exponential form
		res$mu_new = exp(res$mu_new)
	}

	return(res)
}

####
#### Convergence Deriv cpp ####
####

dconvergence = function(dxi_dbeta, jacob.mat, ll_d2, env, iterMax){

	fixef_sizes = get("fixef_sizes", env)
	Q = length(fixef_sizes)
	accDeriv = get("accDeriv", env)
	derivDiffConv = get("derivDifficultConvergence", env)
	mem.clean = get("mem.clean", env)

	if(mem.clean){
	    gc()
	}

	if(accDeriv && derivDiffConv && Q > 2){
		# in case of complex cases: it's more efficient
		# to initialize the first two clusters

		res = dconv_acc(dxi_dbeta, jacob.mat, ll_d2, env, iterMax, only2 = TRUE)
		dxi_dbeta = res$dxi_dbeta
	}


	if(Q == 1){
		# calculer single en cpp
		dxi_dbeta = dconv_single(jacob.mat, ll_d2, env)
		iter = 1

	} else {

		# The convergence algorithms
		if(accDeriv){
			res = dconv_acc(dxi_dbeta, jacob.mat, ll_d2, env, iterMax)
			dxi_dbeta = res$dxi_dbeta
			iter = res$iter
		} else {
			res = dconv_seq(dxi_dbeta, jacob.mat, ll_d2, env, iterMax)
			dxi_dbeta = res$dxi_dbeta
			iter = res$iter
		}
	}

	return(list(dxi_dbeta = dxi_dbeta, iter = iter))
}

dconv_single = function(jacob.mat, ll_d2, env){

	# data:
	jacob_vector = as.vector(jacob.mat)
	n_vars = ncol(jacob.mat)
	nb_cluster_all = get("fixef_sizes", env)
	dum_vector = get("fixef_id_vector", env)
	nb_coef = nb_cluster_all[[1]]

	dxi_dbeta = update_deriv_single(n_vars, nb_coef, ll_d2, jacob_vector, dum_vector)

	return(dxi_dbeta)
}

dconv_seq = function(dxi_dbeta, jacob.mat, ll_d2, env, iterMax){

	# Parameters
	jacob_vector = as.vector(jacob.mat)
	n_vars = ncol(jacob.mat)
	nb_cluster_all = get("fixef_sizes", env)
	dum_vector = get("fixef_id_vector", env)
	deriv_init_vector = as.vector(dxi_dbeta)
	deriv.tol = get("deriv.tol", env)

	Q = length(nb_cluster_all)

	if(Q == 2){
		setup_poisson_fixedcost(env)
		info = get("fixedCostPoisson", env)

		res <- cpp_derivconv_seq_2(iterMax = iterMax, diffMax = deriv.tol, n_vars = n_vars, nb_cluster_all = nb_cluster_all, n_cells = info$n_cells, index_i = info$index_i, index_j = info$index_j, order = info$order, ll_d2 = ll_d2, jacob_vector = jacob_vector, deriv_init_vector = deriv_init_vector, dum_vector = dum_vector)
	} else {
		res <- cpp_derivconv_seq_gnl(iterMax = iterMax, diffMax = deriv.tol, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
	}

	return(list(dxi_dbeta = res$dxi_dbeta, iter = res$iter))
}

dconv_acc = function(dxi_dbeta, jacob.mat, ll_d2, env, iterMax, only2 = FALSE){

	# Parameters
	jacob_vector = as.vector(jacob.mat)
	n_vars = ncol(jacob.mat)
	nb_cluster_all = get("fixef_sizes", env)
	dum_vector = get("fixef_id_vector", env)
	deriv_init_vector = as.vector(dxi_dbeta)
	deriv.tol = get("deriv.tol", env)

	if(only2){
		# we update everything needed
		nb_cluster_all = nb_cluster_all[1:2]
		dum_vector = dum_vector[1:(2*nrow(jacob.mat))]
	}

	Q = length(nb_cluster_all)

	if(Q == 2){
		setup_poisson_fixedcost(env)
		info = get("fixedCostPoisson", env)

		res <- cpp_derivconv_acc_2(iterMax = iterMax, diffMax = deriv.tol, n_vars = n_vars, nb_cluster_all = nb_cluster_all, n_cells = info$n_cells, index_i = info$index_i, index_j = info$index_j, order = info$order, ll_d2 = ll_d2, jacob_vector = jacob_vector, deriv_init_vector = deriv_init_vector, dum_vector = dum_vector)
	} else {
		res <- cpp_derivconv_acc_gnl(iterMax = iterMax, diffMax = deriv.tol, n_vars = n_vars, nb_cluster_all = nb_cluster_all, ll_d2 = ll_d2, jacob_vector = jacob_vector, deriv_init_vector = deriv_init_vector, dum_vector = dum_vector)
	}

	return(list(dxi_dbeta = res$dxi_dbeta, iter = res$iter))
}


####
#### Misc FE ####
####

setup_poisson_fixedcost = function(env){

	# We set up only one
	if("fixedCostPoisson" %in% names(env)){
		return(NULL)
	}

	ptm = proc.time()

	fixef_id = get("fixef_id_list",env)

	dum_A = as.integer(fixef_id[[1]])
	dum_B = as.integer(fixef_id[[2]])

	myOrder = order(dum_A, dum_B)
	index_i = dum_A[myOrder] - 1L
	index_j = dum_B[myOrder] - 1L

	n_cells = get_n_cells(index_i, index_j)

	res = list(n_i = max(dum_A), n_j = max(dum_B), n_cells = n_cells, index_i = index_i, index_j = index_j, order = myOrder - 1L)

	assign("fixedCostPoisson", res, env)

	verbose = get("verbose", env)
	if(verbose >= 2) cat("Poisson fixed-cost setup: ", (proc.time()-ptm)[3], "s\n", sep = "")
}

setup_gaussian_fixedcost = function(env){

	# We set up only one
	if("fixedCostGaussian" %in% names(env)){
		return(NULL)
	}

	ptm = proc.time()

	invTableCluster_vector = get("fixef_invTable", env)
	dum_vector = get("fixef_id_vector", env) # already minus 1
	fixef_id = get("fixef_id_list", env)

	dum_i = as.integer(fixef_id[[1]])
	dum_j = as.integer(fixef_id[[2]])

	n_i = max(dum_i)
	n_j = max(dum_j)

	myOrder = order(dum_i, dum_j)
	index_i = dum_i[myOrder] - 1L
	index_j = dum_j[myOrder] - 1L

	n_cells = get_n_cells(index_i, index_j)
	res = cpp_fixed_cost_gaussian(n_i, n_cells, index_i, index_j, myOrder - 1L, invTableCluster_vector, dum_vector)
	res$n_i = n_i
	res$n_j = n_j
	res$n_cells = n_cells

	assign("fixedCostGaussian", res, env)

	verbose = get("verbose", env)
	if(verbose >= 2) cat("Gaussian fixed-cost setup: ", (proc.time()-ptm)[3], "s\n", sep = "")
}

####
#### Parallel Functions ####
####

# In this section, we create all the functions that will be parallelized

rpar_exp = function(x, env){
	# fast exponentiation
	nthreads = get("nthreads", env)

	if(nthreads == 1){
		# simple exponentiation
		return(exp(x))
	} else {
		# parallelized one
		return(cpppar_exp(x, nthreads))
	}

}

rpar_log = function(x, env){
	# fast log
	nthreads = get("nthreads", env)

	if(nthreads == 1){
		# simple log
		return(log(x))
	} else {
		# parallelized one
		return(cpppar_log(x, nthreads))
	}

}

rpar_lgamma = function(x, env){
	# fast lgamma
	nthreads = get("nthreads", env)

	if(nthreads == 1){
		# lgamma via cpp is faster
		return(cpp_lgamma(x))
	} else {
		# parallelized one
		return(cpppar_lgamma(x, nthreads))
	}

}

rpar_digamma = function(x, env){
	nthreads = get("nthreads", env)

	if(nthreads == 1){
		# digamma via R is as fast => no need cpp
		return(digamma(x))
	} else {
		# parallelized one
		return(cpppar_digamma(x, nthreads))
	}

}

rpar_trigamma = function(x, env){
	nthreads = get("nthreads", env)

	if(nthreads == 1){
		# trigamma via R is as fast => no need cpp
		return(trigamma(x))
	} else {
		# parallelized one
		return(cpppar_trigamma(x, nthreads))
	}

}

rpar_log_a_exp = function(a, mu, exp_mu, env){
	# compute log_a_exp in a fast way
	nthreads = get("nthreads", env)

	if(nthreads == 1){
		# cpp is faster
		return(cpp_log_a_exp(a, mu, exp_mu))
	} else {
		# parallelized one
		return(cpppar_log_a_exp(nthreads, a, mu, exp_mu))
	}
}

