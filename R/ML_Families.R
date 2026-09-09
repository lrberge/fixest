

#-------------------------------#
#### ===== WEIGHTS UTIL ===== ####
#_______________________________#

# All the ML families below support weighted estimations. The weights are read
# from the estimation environment at call time -- and not captured once when the
# family is created -- because split-sample estimations reshape the environment
# and subset the weights accordingly. Every caller must therefore pass the `env`
# of the estimation at hand.
#
# The convention is the one used everywhere else in fixest: `weights.value` is
# the scalar 1 when the estimation is unweighted, and a vector of length nobs
# otherwise.

ml_weights = function(env){

	if(missing(env) || !is.environment(env)){
		stop("Internal error in the ML families: the estimation environment must be ",
		     "passed in order to access the weights.")
	}

	get0("weights.value", envir = env, ifnotfound = 1)
}

# sum(w * x), avoiding the copy of x when there is no weight
wsum = function(w, x){
	if(length(w) == 1) sum(x) else sum(w * x)
}

# w * x, avoiding the copy of x when there is no weight
wmul = function(w, x){
	if(length(w) == 1) x else w * x
}

# the "number of observations": the sum of the weights when weighted
wn = function(w, y){
	if(length(w) == 1) length(y) else sum(w)
}

#---------------------------#
#### ===== POISSON ===== ####
#___________________________#

ml_poisson = function(){
	# Cette fonction renvoie une famille de fonctions
	ll = function(y, mu, exp_mu, env, ...){

		w = ml_weights(env)

		# computing the lfactorial is costly
		# note that the cached value is the weighted one: the cache lives in the
		# estimation env, hence it is never shared across estimations
		if("lfactorial" %in% names(env)){
			lfact = get("lfactorial", env)
		} else {
			lfact = wsum(w, rpar_lgamma(y + 1, env))
			assign("lfactorial", lfact, env)
		}

		wsum(w, y*mu - exp_mu) - lfact
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, env, ...){
		wmul(ml_weights(env), c(y - exp_mu))
	}

	# Derivee seconde
	ll_d2 = function(y, mu, exp_mu, env, ...){
		wmul(ml_weights(env), -c(exp_mu))
	}

	expected.predictor = function(mu, exp_mu, ...){
		exp_mu
	}

	linearFromExpected = function(exp_pred){
		log(exp_pred)
	}

	closedFormDummies = function(dum, y, mu, env, sum_y, orderCluster, tableCluster, ...){
		# We send only the dummies (not the vector of dummies)
		# sum_y is already weighted

		mu_dum = cpp_tapply_vsum(length(tableCluster), wmul(ml_weights(env), c(exp(mu))), dum)
		log(sum_y) - log(mu_dum)
	}

	return(list(ll = ll, expected.predictor = expected.predictor, ll_dl = ll_dl,
	            ll_d2 = ll_d2, closedFormDummies = closedFormDummies,
	            linearFromExpected = linearFromExpected))
}

#------------------------#
#### ===== LOGIT ==== ####
#------------------------#


ml_logit = function(){

	ll = function(y, mu, exp_mu, env, ...){
		wsum(ml_weights(env), y*mu - rpar_log_a_exp(1, mu, exp_mu, env))
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, env, ...){
		wmul(ml_weights(env), c(y - 1/(1+1/exp_mu)))
	}

	# Derivee seconde
	ll_d2 = function(y, mu, exp_mu, env, ...){
		wmul(ml_weights(env), - c(1/ ( (1+exp_mu) * (1+1/exp_mu) )))
	}

	expected.predictor = function(mu, exp_mu, ...){
		1/(1+1/exp_mu)
	}

	linearFromExpected = function(exp_pred){
		log(exp_pred) - log(1 - exp_pred)
	}

	guessDummy = function(sum_y, n_group, mu, ...){
		# guess for the dummy:
		log(sum_y) - log(n_group-sum_y) - mu
	}

	guessExpDummy = function(sum_y, n_group, exp_mu, ...){
		#guess for the dummy:
		sum_y / (n_group-sum_y) / exp_mu
	}

	return(list(ll = ll, expected.predictor = expected.predictor, ll_dl = ll_dl,
	            ll_d2 = ll_d2, guessDummy = guessDummy,
	            guessExpDummy = guessExpDummy, linearFromExpected = linearFromExpected))
}


#*************************#
#### ===== NEGBIN ==== ####
#*************************#

ml_negbin = function(){
	# Cette fonction renvoie une famille de fonctions
	ll = function(y, mu, exp_mu, env, coef, ...){

		w = ml_weights(env)

		# computing the lgamma is costly
		if("lgamma" %in% names(env)){
			lgamm = get("lgamma", env)
		} else {
			lgamm = wsum(w, rpar_lgamma(y + 1, env))
			assign("lgamma", lgamm, env)
		}

		theta = coef[".theta"]
		N = wn(w, y)

		wsum(w, rpar_lgamma(theta+y, env) + y*mu - (theta+y)*rpar_log_a_exp(theta, mu, exp_mu, env)) - lgamm + N*(- lgamma(theta) + theta*log(theta))
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, coef, env, ...){
		theta = coef[".theta"]

		wmul(ml_weights(env), c(y - (theta+y) / (theta/exp_mu + 1)))
	}

	# Derivee croisee
	ll_dx_dother = function(y, mu, exp_mu, coef, env, ...){
		#Means the second derivative of the LL wrt to the linear part and theta
		theta = coef[".theta"]

		wmul(ml_weights(env), c( -1/(theta/exp_mu + 1)^2 + y/( (theta/exp_mu + 1) * (theta + exp_mu) ) ))
	}

	# Derivee seconde:
	ll_d2 = function(y, mu, exp_mu, coef, env, ...){
		theta = coef[".theta"]

		wmul(ml_weights(env), - theta * (theta+y) / ( (theta/exp_mu + 1) * (theta + exp_mu) ))
	}

	grad.theta = function(theta, y, mu, exp_mu, env, ...){
		w = ml_weights(env)
		N = wn(w, y)

		wsum(w, rpar_digamma(theta+y, env) - rpar_log_a_exp(theta, mu, exp_mu, env) - (theta+y)/(theta+exp_mu) ) + N*(- psigamma(theta) + log(theta) + 1 )
	}

	scores.theta = function(theta, y, mu, exp_mu, env){
		wmul(ml_weights(env), rpar_digamma(theta+y, env) - rpar_log_a_exp(theta, mu, exp_mu, env) - (theta+y)/(theta+exp_mu) + (- psigamma(theta) + log(theta) + 1))
	}

	hess.theta = function(theta, y, mu, exp_mu, env){
		w = ml_weights(env)
		N = wn(w, y)

		wsum(w, rpar_trigamma(theta+y, env) - 1/(theta+exp_mu) + y/(theta+exp_mu)^2 - 1/( (theta/exp_mu + 1) * (theta + exp_mu) ) ) + N*(- psigamma(theta, 1) + 1/theta)
	}

	hess.thetaL = function(theta, jacob.mat, y, dxi_dbeta, dxi_dother, ll_d2, ll_dx_dother){
		# ll_d2 and ll_dx_dother are already weighted
		res = crossprod((jacob.mat+dxi_dbeta), dxi_dother*ll_d2+ll_dx_dother)
		return(res)
	}

	hess_theta_part = function(theta, y, mu, exp_mu, dxi_dother, ll_dx_dother, ll_d2, env){
		# La derivee vav de theta en prenant en compte les dummies
		# ll_d2 and ll_dx_dother are already weighted
		d2ll_d2theta = hess.theta(theta, y, mu, exp_mu, env)
		res = sum(dxi_dother^2*ll_d2 + 2*dxi_dother*ll_dx_dother) + d2ll_d2theta
		return(res)
	}

	expected.predictor = function(mu, exp_mu, ...){
		exp_mu
	}

	linearFromExpected = function(exp_pred){
		log(exp_pred)
	}

	guessDummy = function(sum_y, n_group, mu, ...){
		# guess for the dummy:
		log(sum_y) - log(n_group) - mu
	}

	guessExpDummy = function(sum_y, n_group, exp_mu, ...){
		# guess for the dummy:
		sum_y / n_group / exp_mu
	}

	ll0_theta = function(theta, y, mean_y, invariant, env){
		# La fonction minimise, on renvoie "-"
		# mean_y and invariant are already weighted (see get_model_null)
		w = ml_weights(env)
		N = wn(w, y)

		ll = wsum(w, rpar_lgamma(theta+y, env)) + invariant + N*(- mean_y*log(theta+mean_y) - lgamma(theta) + theta*log(theta) - theta*log(theta+mean_y))

		-ll
	}

	grad0_theta = function(theta, y, mean_y, env, ...){
		# La fonction minimise, on renvoie "-"
		w = ml_weights(env)
		N = wn(w, y)

		grad.theta = wsum(w, rpar_digamma(theta+y, env)) + N*(- mean_y/(theta+mean_y) - psigamma(theta) + log(theta) + 1 - log(theta+mean_y) - theta/(theta+mean_y))

		- grad.theta
	}

	hess0_theta = function(theta, y, mean_y, env, ...){
		# La fonction minimise, on renvoie "-"
		w = ml_weights(env)
		N = wn(w, y)

		hess.theta = wsum(w, rpar_trigamma(theta+y, env)) + N*(- psigamma(theta, 1) + 1/theta - 1/(theta+mean_y))

		- as.matrix(hess.theta)
	}

	return(list(ll = ll, expected.predictor = expected.predictor,
	            linearFromExpected = linearFromExpected, ll0_theta = ll0_theta,
	            grad0_theta = grad0_theta, hess0_theta = hess0_theta, hess.theta = hess.theta,
	            hess.thetaL = hess.thetaL, hess_theta_part = hess_theta_part,
	            grad.theta = grad.theta, scores.theta = scores.theta, ll_dl = ll_dl,
	            ll_d2 = ll_d2, guessDummy = guessDummy, ll_dx_dother = ll_dx_dother,
	            guessExpDummy = guessExpDummy))
}


#---------------------------#
#### ===== GAUSSIAN ==== ####
#---------------------------#


ml_gaussian = function(){

	# the weighted variance of the residuals: sum(w*(y-mu)^2) / sum(w)
	sigma2 = function(y, mu, w){
		wsum(w, (y-mu)^2) / wn(w, y)
	}

	ll = function(y, mu, exp_mu, env, ...){
		w = ml_weights(env)
		s2 = sigma2(y, mu, w)
		n = wn(w, y)
		-1/2/s2*wsum(w, (y-mu)^2) - n*log(sqrt(s2)) - n*log(2*pi)/2
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, env, ...){
		w = ml_weights(env)
		s2 = sigma2(y, mu, w)
		wmul(w, (y-mu)/s2)
	}

	# Derivee seconde
	ll_d2 = function(y, mu, exp_mu, env, ...){
		w = ml_weights(env)
		s2 = sigma2(y, mu, w)
		wmul(w, rep(-1/s2, length(y)))
	}

	expected.predictor = function(mu, exp_mu, ...){
		mu
	}

	linearFromExpected = function(exp_pred){
		exp_pred
	}

	closedFormDummies = function(dum, y, mu, env, sum_y, orderCluster, tableCluster, ...){
		# We send only the dummies (not the vector of dummies)
		# sum_y is already weighted, hence the divisor is the sum of the weights

		w = ml_weights(env)
		sum_w = if(length(w) == 1) tableCluster else cpp_tapply_vsum(length(tableCluster), w, dum)
		(sum_y - cpp_tapply_vsum(length(tableCluster), wmul(w, mu), dum)) / sum_w
	}

	return(list(ll = ll, expected.predictor = expected.predictor,
	            linearFromExpected = linearFromExpected, ll_dl = ll_dl,
	            ll_d2 = ll_d2, closedFormDummies = closedFormDummies))
}
