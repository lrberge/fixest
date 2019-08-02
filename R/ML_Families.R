

#***************************#
#### ===== POISSON ===== ####
#***************************#

ml_poisson = function(){
	# Cette fonction renvoie une famille de fonctions
	ll = function(y, mu, exp_mu, env, ...){

		# computing the lfactorial is costly
		if("lfactorial" %in% names(env)){
			lfact = get("lfactorial", env)
		} else {
			lfact = sum(rpar_lgamma(y + 1, env))
			assign("lfactorial", lfact, env)
		}

		sum(y*mu - exp_mu) - lfact
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, ...){
		c(y - exp_mu)
	}

	# Derivee seconde
	ll_d2 = function(y, mu, exp_mu, ...){
		-c(exp_mu)
	}

	expected.predictor = function(mu, exp_mu, ...){
		exp_mu
	}

	linearFromExpected = function(exp_pred){
		log(exp_pred)
	}

	closedFormDummies = function(dum, y, mu, env, sum_y, orderCluster, tableCluster, ...){
		# We send only the dummies (not the vector of dummies)

		mu_dum = cpp_tapply_vsum(length(tableCluster), c(exp(mu)), dum)
		log(sum_y) - log(mu_dum)
	}

	return(list(ll=ll, expected.predictor=expected.predictor, ll_dl=ll_dl, ll_d2=ll_d2, closedFormDummies=closedFormDummies, linearFromExpected=linearFromExpected))
}

#************************#
#### ===== LOGIT ==== ####
#************************#


ml_logit = function(){

	ll = function(y, mu, exp_mu, env, ...){
		sum(y*mu - rpar_log_a_exp(1, mu, exp_mu, env))
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, ...){
		c(y - 1/(1+1/exp_mu))
	}

	# Derivee seconde
	ll_d2 = function(y, mu, exp_mu, ...){
		- c(1/ ( (1+exp_mu) * (1+1/exp_mu) ))
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

	return(list(ll=ll, expected.predictor=expected.predictor, ll_dl=ll_dl, ll_d2=ll_d2, guessDummy=guessDummy, guessExpDummy=guessExpDummy, linearFromExpected=linearFromExpected))
}


#*************************#
#### ===== NEGBIN ==== ####
#*************************#

ml_negbin = function(){
	# Cette fonction renvoie une famille de fonctions
	ll = function(y, mu, exp_mu, env, coef, ...){
		# computing the lgamma is costly
		if("lgamma" %in% names(env)){
			lgamm = get("lgamma", env)
		} else {
			lgamm = sum(rpar_lgamma(y + 1, env))
			assign("lgamma", lgamm, env)
		}

		theta = coef[".theta"]
		N = length(y)

		sum(rpar_lgamma(theta+y, env) + y*mu - (theta+y)*rpar_log_a_exp(theta, mu, exp_mu, env)) - lgamm + N*(- lgamma(theta) + theta*log(theta))
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, coef, env, ...){
		theta = coef[".theta"]

		c(y - (theta+y) / (theta/exp_mu + 1))
	}

	# Derivee croisee
	ll_dx_dother = function(y, mu, exp_mu, coef, env, ...){
		#Means the second derivative of the LL wrt to the linear part and theta
		theta = coef[".theta"]

		c( -1/(theta/exp_mu + 1)^2 + y/( (theta/exp_mu + 1) * (theta + exp_mu) ) )
	}

	# Derivee seconde:
	ll_d2 = function(y, mu, exp_mu, coef, env, ...){
		theta = coef[".theta"]

		- theta * (theta+y) / ( (theta/exp_mu + 1) * (theta + exp_mu) )
	}

	grad.theta = function(theta, y, mu, exp_mu, env, ...){
		N = length(y)

		sum( rpar_digamma(theta+y, env) - rpar_log_a_exp(theta, mu, exp_mu, env) - (theta+y)/(theta+exp_mu) ) + N*(- psigamma(theta) + log(theta) + 1 )
	}

	scores.theta = function(theta, y, mu, exp_mu, env){
		rpar_digamma(theta+y, env) - rpar_log_a_exp(theta, mu, exp_mu, env) - (theta+y)/(theta+exp_mu) + (- psigamma(theta) + log(theta) + 1)
	}

	hess.theta = function(theta, y, mu, exp_mu, env){
		N = length(y)

		sum( rpar_trigamma(theta+y, env) - 1/(theta+exp_mu) + y/(theta+exp_mu)^2 - 1/( (theta/exp_mu + 1) * (theta + exp_mu) ) ) + N*(- psigamma(theta, 1) + 1/theta)
	}

	hess.thetaL = function(theta, jacob.mat, y, dxi_dbeta, dxi_dother, ll_d2, ll_dx_dother){
		res = crossprod((jacob.mat+dxi_dbeta), dxi_dother*ll_d2+ll_dx_dother)
		return(res)
	}

	hess_theta_part = function(theta, y, mu, exp_mu, dxi_dother, ll_dx_dother, ll_d2, env){
		# La derivee vav de theta en prenant en compte les dummies
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
		N = length(y)

		ll = sum(rpar_lgamma(theta+y, env)) + invariant + N*(- mean_y*log(theta+mean_y) - lgamma(theta) + theta*log(theta) - theta*log(theta+mean_y))

		-ll
	}

	grad0_theta = function(theta, y, mean_y, env, ...){
		# La fonction minimise, on renvoie "-"
		N = length(y)

		grad.theta = sum(rpar_digamma(theta+y, env)) + N*(- mean_y/(theta+mean_y) - psigamma(theta) + log(theta) + 1 - log(theta+mean_y) - theta/(theta+mean_y))

		- grad.theta
	}

	hess0_theta = function(theta, y, mean_y, env, ...){
		# La fonction minimise, on renvoie "-"
		N = length(y)

		hess.theta = sum(rpar_trigamma(theta+y, env)) + N*(- psigamma(theta, 1) + 1/theta - 1/(theta+mean_y))

		- as.matrix(hess.theta)
	}

	return(list(ll=ll, expected.predictor=expected.predictor, linearFromExpected=linearFromExpected, ll0_theta=ll0_theta, grad0_theta=grad0_theta, hess0_theta=hess0_theta, hess.theta=hess.theta, hess.thetaL=hess.thetaL, hess_theta_part=hess_theta_part, grad.theta=grad.theta, scores.theta=scores.theta, ll_dl=ll_dl, ll_d2=ll_d2, guessDummy=guessDummy, ll_dx_dother=ll_dx_dother, guessExpDummy=guessExpDummy))
}


#***************************#
#### ===== GAUSSIAN ==== ####
#***************************#


ml_gaussian = function(){

	ll = function(y, mu, exp_mu, env, ...){
		sigma = sqrt(mean((y-mu)^2))
		n = length(y)
		-1/2/sigma^2*sum((y-mu)^2) - n*log(sigma) - n*log(2*pi)/2
	}

	# Derivee
	ll_dl = function(y, mu, exp_mu, env, ...){
		sigma = sqrt(mean((y-mu)^2))
		(y-mu)/sigma^2
	}

	# Derivee seconde
	ll_d2 = function(y, mu, exp_mu, env, ...){
		sigma = sqrt(mean((y-mu)^2))
		rep(-1/sigma^2, length(y))
	}

	expected.predictor = function(mu, exp_mu, ...){
		mu
	}

	linearFromExpected = function(exp_pred){
		exp_pred
	}

	closedFormDummies = function(dum, y, mu, env, sum_y, orderCluster, tableCluster, ...){
		# We send only the dummies (not the vector of dummies)

		(sum_y - cpp_tapply_vsum(length(tableCluster), mu, dum)) / tableCluster
	}

	return(list(ll=ll, expected.predictor=expected.predictor, linearFromExpected=linearFromExpected, ll_dl=ll_dl, ll_d2=ll_d2, closedFormDummies=closedFormDummies))
}

