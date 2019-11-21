#------------------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Oct 19 09:20:25 2019
# ~: Functions facilitating control within funs
# NOTA to myself: Do not develop these functions here!
# but in lbmisctools!
#------------------------------------------------------#



control_variable = function(x, myType, prefix, name, charVec, mustBeThere = FALSE){
	# Format of my types:
	#   - single => must be of lenght one
	#   - Vector => must be a vector
	#   - Matrix => must be a matrix
	#   - GE/GT/LE/LT: greater/lower than a given value
	#   - predefinedType => eg: numeric, integer, etc
	#   - match.arg => very specific => should match the charVec
	#   - noNA: NAs not allowed
	#   - null: null type allowed
	# If there is a parenthesis => the class must be of specified types:
	# ex: "(list, data.frame)" must be a list of a data.frame
	
	ignore.case = TRUE
	
	if(missing(prefix)){
		msg = deparse(sys.calls()[[sys.nframe()-1]])[1] # call can have svl lines
		if(length(msg) > 1) browser()
		nmax = 40
		if(nchar(msg) > nmax) msg = paste0(substr(msg, 1, nmax-1), "...")
		prefix = paste0(msg, ": ")
	}
	
	x_name = deparse(substitute(x))
	if(missing(name)) name = x_name
	
	firstMsg = paste0(prefix, "The argument '", name, "' ")
	
	if(missing(x)){
		if(mustBeThere){
			stop(firstMsg, "is missing => it must be provided.", call. = FALSE)
		} else {
			return(invisible(NULL))
		}
	}
	
	
	# simple function to extract a pattern
	# ex: if my type is VectorIntegerGE1 => myExtract("GE[[:digit:]]+","VectorIntegerGE1") => 1
	myExtract = function(expr, text, trim=2){
		start = gregexpr(expr,text)[[1]] + trim
		length = attr(start,"match.length") - trim
		res = substr(text,start,start+length-1)
		as.numeric(res)
	}
	
	#
	# General types handling
	#
	
	loType = tolower(myType)
	
	# null type is caught first
	if(grepl("null", loType)){
		if(is.null(x)) return(invisible(NULL))
	}
	
	if(grepl("single", loType)){
		if(length(x) != 1) stop(firstMsg,"must be of length one.", call. = FALSE)
	}
	
	if(grepl("vector", loType)){
		if(!isVector(x)) stop(firstMsg,"must be a vector.", call. = FALSE)
		if(is.list(x)) stop(firstMsg,"must be a vector (and not a list).", call. = FALSE)
	}
	
	res = checkTheTypes(loType, x)
	if(!res$OK) stop(firstMsg,res$message, call. = FALSE)
	
	# # INTEGER is a restrictive type that deserves some explanations (not included in getTheTypes)
	# if(grepl("integer",loType)){
	#     if(grepl("single",loType)){
	#         if(!is.numeric(x)) stop(firstMsg,"must be an integer (right now it is not even numeric).", call. = FALSE)
	#         if(!(is.integer(x) || x%%1==0)) stop(firstMsg,"must be an integer.", call. = FALSE)
	#     } else {
	#         if(!is.numeric(x)) stop(firstMsg,"must be composed of integers (right now it is not even numeric).", call. = FALSE)
	#         if(!(is.integer(x) || all(x%%1==0))) stop(firstMsg,"must be composed of integers.", call. = FALSE)
	#     }
	# }
	
	if(grepl("nona", loType)){
		if(any(is.na(x))){
			stop(firstMsg,"contains NAs, this is not allowed.", call. = FALSE)
		}
	}
	
	# GE: greater or equal // GT: greater than // LE: lower or equal // LT: lower than
	if(is.numeric(x)){
		x = x[!is.na(x)]
		
		if(grepl("ge[[:digit:]]+",loType)){
			n = myExtract("ge[[:digit:]]+", loType)
			if( !all(x>=n) ) stop(firstMsg,"must be greater than, or equal to, ", n, ".", call. = FALSE)
		}
		if(grepl("gt[[:digit:]]+",loType)){
			n = myExtract("gt[[:digit:]]+", loType)
			if( !all(x>n) ) stop(firstMsg,"must be strictly greater than ", n, ".", call. = FALSE)
		}
		if(grepl("le[[:digit:]]+",loType)){
			n = myExtract("le[[:digit:]]+", loType)
			if( !all(x<=n) ) stop(firstMsg,"must be lower than, or equal to, ", n, ".", call. = FALSE)
		}
		if(grepl("lt[[:digit:]]+",loType)){
			n = myExtract("lt[[:digit:]]+", loType)
			if( !all(x<n) ) stop(firstMsg,"must be strictly lower than ", n, ".", call. = FALSE)
		}
	}
	
	#
	# Specific Types Handling
	#
	
	if(grepl("match.arg", loType)){
		if(ignore.case){
			x = toupper(x)
			newCharVec = toupper(charVec)
		} else {
			newCharVec = charVec
		}
		
		if( is.na(pmatch(x, newCharVec)) ){
			n = length(charVec)
			if(n == 1){
				msg = paste0("'",charVec,"'")
			} else {
				msg = paste0("'", paste0(charVec[1:(n-1)], collapse="', '"), "' or '",charVec[n],"'")
			}
			stop(firstMsg,"must be one of:\n",msg,".", call. = FALSE)
		} else {
			qui = pmatch(x, newCharVec)
			return(charVec[qui])
		}
	}
}

matchTypeAndSetDefault = function(myList, myDefault, myTypes, prefix){
	# Cette fonction:
	#   i) verifie que tous les elements de la liste sont valides
	#   ii) mes les valeurs par defauts si elles certaines valeurs sont manquantes
	#   iii) Envoie des messages d'erreur si les typages ne sont pas bons
	# En fait cette fonction "coerce" myList en ce qu'il faudrait etre (donne par myDefault)
	
	# 1) check that the names of the list are valid
	if(is.null(myList)) myList = list()
	list_names = names(myList)
	
	if(length(list_names)!=length(myList) || any(list_names=="")){
		stop(prefix,"The elements of the list should be named.", call. = FALSE)
	}
	
	obj_names = names(myDefault)
	
	isHere = pmatch(list_names,obj_names)
	
	if(anyNA(isHere)){
		if(sum(is.na(isHere))==1) stop(prefix, "The following argument is not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
		else stop(prefix, "The following arguments are not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
	}
	
	# 2) We set the default values and run Warnings
	res = list()
	for(i in 1:length(obj_names)){
		obj = obj_names[i]
		qui = which(isHere==i) # qui vaut le numero de l'objet dans myList
		type = myTypes[i] # we extract the type => to control for "match.arg" type
		if(length(qui)==0){
			# we set to the default if it's missing
			if(type == "match.arg") {
				res[[obj]] = myDefault[[i]][1]
			} else {
				res[[obj]] = myDefault[[i]]
			}
		} else {
			# we check the integrity of the value
			val = myList[[qui]]
			if(type == "match.arg"){
				# If the value is to be a match.arg => we use our controls and not
				# directly the one of the function match.arg()
				charVec = myDefault[[i]]
				control_variable(val, "singleCharacterMatch.arg", prefix, obj, charVec)
				val = match.arg(val, charVec)
			} else {
				control_variable(val, type, prefix, obj)
			}
			
			res[[obj]] = val
		}
	}
	
	return(res)
}



checkTheTypes = function(str, x){
	# This function takes in a character string describing the types of the
	# element x => it can be of several types
	
	# types that are controlled for:
	allTypes = c("numeric", "integer", "character", "logical", "list", "data.frame", "matrix", "factor", "formula")
	
	OK = FALSE
	message = c()
	
	for(type in allTypes){
		
		if(grepl(type, str)){
			
			# we add the type of the control
			message = c(message, type)
			
			if(type == "numeric"){
				if(!OK & is.numeric(x)){
					OK = TRUE
				}
			} else if(type == "integer"){
				if(is.numeric(x) && (is.integer(x) || all(x%%1==0))){
					OK = TRUE
				}
			} else if(type == "character"){
				if(is.character(x)){
					OK = TRUE
				}
			} else if(type == "logical"){
				if(is.logical(x)){
					OK = TRUE
				}
			} else if(type == "list"){
				if(is.list(x)){
					OK = TRUE
				}
			} else if(type == "data.frame"){
				if(is.data.frame(x)){
					OK = TRUE
				}
			} else if(type == "matrix"){
				if(is.matrix(x)){
					OK = TRUE
				}
			} else if(type == "factor"){
				if(is.factor(x)){
					OK = TRUE
				}
			}  else if(type == "formula"){
				if(length(class(x)) == 1 && class(x) == "formula"){
					OK = TRUE
				}
			}
		}
		
		if(OK) break
	}
	
	if(length(message) == 0){
		OK = TRUE #ie there is no type to be searched
	} else if(length(message) >= 3){
		n = length(message)
		message = paste0("must be of type: ",  paste0(message[1:(n-1)], collapse = ", "), " or ", message[n], ".")
	} else {
		message = paste0("must be of type: ",  paste0(message, collapse = " or "), ".")
	}
	
	
	return(list(OK=OK, message=message))
}


check_arg = function(x, type, message, call_depth = 0){
	# function that makes it easy to check arguments:
	#    provides precise and meaningful error messages
	# only for scalars or vectors
	# must be of single types
	# possible types: logical/numeric/integer/character
	# osf => one sided formula
	# tsf: two sided formula
	
	type = tolower(type)
	
	# function for greater than/ lower than, etc
	myExtract = function(expr, trim=2){
		# type is global
		start = gregexpr(expr, type)[[1]] + trim
		length = attr(start, "match.length") - trim
		res = substr(type, start, start + length - 1)
		as.numeric(res)
	}
	
	# Creating the message with the most possible precision
	
	if(missing(message)){
		while(TRUE){
			# while: Trick to exist the condition at any time
			message = paste0("Argument '", deparse(substitute(x)),"' must be ")
			
			if(grepl("osf", type)){
				message = paste0(message, "a one sided formula.")
				break
			} else if(grepl("tsf", type)){
				message = paste0(message, "a two sided formula.")
				break
			} else if(grepl("logical", type)){
				my_type = "logical"
			} else if(grepl("character", type)){
				my_type = "character"
			} else if(grepl("integer", type)){
				my_type = "integer"
			} else if(grepl("numeric", type)){
				my_type = "numeric"
			} else {
				# no requirement
				my_type = ""
			}
			
			if(grepl("single", type)){
				message = paste0(message, "a single ", my_type)
			} else {
				# This is a vector
				message = paste0(message, "a ", my_type, " vector")
				message = gsub("a int", "an int", message)
				message = gsub(" +", " ", message)
			}
			
			if(grepl("integer|numeric", type)){
				
				first_msg = ifelse(grepl("vector", type), " of values", "")
				
				if(grepl(expr <- "ge[[:digit:]]+", type)){
					n = myExtract(expr)
					message = paste0(message, first_msg, " greater than, or equal to, ", n)
					first_msg = " and"
				}
				
				if(grepl(expr <- "gt[[:digit:]]+", type)){
					add_and = TRUE
					n = myExtract(expr)
					message = paste0(message, first_msg, " strictly greater than ", n)
					first_msg = " and"
				}
				
				if(grepl(expr <- "le[[:digit:]]+", type)){
					n = myExtract(expr)
					message = paste0(message, first_msg, " lower than, or equal to, ", n)
				}
				
				if(grepl(expr <- "lt[[:digit:]]+", type)){
					n = myExtract(expr)
					message = paste0(message, first_msg, " strictly lower than ", n)
				}
			}
			
			message = paste0(message, ". REASON")
			break
		}
	}
	
	stop_now = function(...){
		# message is a global
		
		reason = paste0(...)
		
		# The original call
		my_call = deparse(sys.calls()[[sys.nframe()-(2 + call_depth)]])[1] # call can have svl lines
		nmax = 40
		if(nchar(my_call) > nmax) my_call = paste0(substr(my_call, 1, nmax-1), "...")
		
		# The formatted message
		msg_split = strsplit(message, " ?REASON ?")[[1]]
		
		msg_new = c(msg_split[1], reason, msg_split[-1])
		msg_new = paste(msg_new, collapse = " ")
		
		stop("in ", my_call, ":\n", msg_new, call. = FALSE)
	}
	
	#
	# SPECIAL ARGS
	#
	
	if(missing(x)){
		if(grepl("mbt", type)){
			stop_now("But it is missing.")
		} else {
			return(NULL)
		}
	}
	
	if(grepl("null", type)){
		if(is.null(x)){
			return(NULL)
		}
	}
	
	#
	# FORMULAS
	#
	
	if(grepl("(o|t)sf", type)){
		if(!"formula" %in% class(x)){
			stop_now("But it is currently not a formula (instead it is of class '", class(x)[1], "').")
		}
		
		if(grepl("osf", type) && length(x) == 3){
			stop_now("But it is currently two-sided.")
		}
		
		if(grepl("tsf", type) && length(x) == 2){
			stop_now("But it is currently only one-sided.")
		}
		
		return(NULL)
	}
	
	
	isSingle = FALSE
	if(grepl("single|scalar", type)){
		isSingle = TRUE
		if(length(x) == 0){
			stop_now("But it is of length 0.")
		} else if(length(x) != 1){
			stop_now("But it is of length ", length(x), ".")
		}
	}
	
	if(grepl("character", type) && !is.character(x)){
		stop_now("But it is not of type character.")
	}
	
	if(grepl("logical", type) && !is.logical(x)){
		stop_now("But it is not logical.")
	}
	
	if(grepl("numeric|integer", type) && !is.numeric(x) && !is.logical(x)){
		# logicals are OK as numerics
		stop_now("But it is not numeric.")
	}
	
	if(!grepl("naok", type) && anyNA(x)){
		if(isSingle){
			stop_now("But it is equal to NA.")
		} else {
			stop_now("But it contains NAs.")
		}
	}
	
	if(grepl("integer", type) && !all(x %% 1 == 0)){
		stop_now("But it is not integer (although numeric).")
	}
	
	# Greater than, lower than
	
	if(grepl(expr <- "ge[[:digit:]]+", type)){
		n = myExtract(expr)
		if( any(x < n) ) stop_now("But it is lower than ", n, ".")
	}
	
	if(grepl(expr <- "gt[[:digit:]]+", type)){
		n = myExtract(expr)
		if( any(x == n) ) stop_now("But it is equal to ", n, " (while it should be *striclty* greater).")
		if( any(x < n) ) stop_now("But it is lower than ", n, ".")
	}
	
	if(grepl(expr <- "le[[:digit:]]+", type)){
		n = myExtract(expr)
		if( any(x > n) ) stop_now("But it is greater than ", n, ".")
	}
	
	if(grepl(expr <- "lt[[:digit:]]+", type)){
		n = myExtract(expr)
		if( any(x == n) ) stop_now("But it is equal to ", n, " (while it should be *striclty* lower).")
		if( any(x > n) ) stop_now("But it is greater than ", n, ".")
	}
	
}

# Avoids the problem of multiple lines deparse
deparse_long = function(x){
	dep_x = deparse(x)
	if(length(dep_x) == 1){
		return(dep_x)
	} else {
		return(paste(gsub("^ +", "", dep_x), collapse = ""))
	}
}



####
#### Control Utilities ####
####

enumerate_items = function (x, type, verb = FALSE, addS = FALSE, past = FALSE, or = FALSE, start_verb = FALSE, quote = FALSE){
	# function that enumerates items and add verbs
	# in argument type, you can have a mix of the different arguments, all separated with a "."
	# to add an additional, regular, verb, you can use the underscore: e.g. "s._look.or"
	
	if(!missing(type)){
		args = strsplit(type, "\\.")[[1]]
		verb = intersect(c("is", "has", "contain"), args)
		if(length(verb) == 0){
			qui = grepl("_", args)
			if(any(qui)){
				verb = gsub("^_|s$", "", args[qui])
			} else {
				verb = "no"
			}
		}
		addS = "s" %in% args
		past = "past" %in% args
		or = "or" %in% args
		start_verb = "start" %in% args
		quote = "quote" %in% args
	} else {
		verb = match.arg(as.character(verb), c("is", "has", "no", "contain", "FALSE"))
		if(verb == "FALSE") verb = "no"
	}
	
	n = length(x)
	
	if(past){
		if(verb %in% c("no", "is", "has")){
			verb_format = switch(verb, is = ifelse(n == 1, " was", " were"), no = "", has=" had")
		} else {
			verb_format = paste0(" ", verb, "ed")
		}
	} else {
		if(verb %in% c("no", "is", "has")){
			verb_format = switch(verb, is = ifelse(n == 1, " is", " are"), no = "", has = ifelse(n == 1, " has", " have"))
		} else {
			verb_format = ifelse(n == 1, paste0(" ", verb, "s"), paste0(" ", verb))
		}
		
	}
	
	if (addS) {
		startWord = ifelse(n == 1, " ", "s ")
	} else {
		startWord = ""
	}
	
	if(quote){
		x = paste0("'", x, "'")
	}
	
	if (n == 1) {
		if(!start_verb){
			res = paste0(startWord, x, verb_format)
		} else {
			res = paste0(startWord, gsub(" ", "", verb_format), " ", x)
		}
		
	} else {
		and_or = ifelse(or, " or ", " and ")
		if(!start_verb){
			res = paste0(startWord, paste0(x[-n], collapse = ", "), and_or, x[n], verb_format)
		} else {
			res = paste0(startWord, gsub(" ", "", verb_format), " ", paste0(x[-n], collapse = ", "), and_or, x[n])
		}
		
	}
	
	res
}

ifsingle = function(x, yes, no){
	if(length(x) == 1){
		return(yes)
	} else {
		return(no)
	}
}



isVector = function(x){
	# it seems that when you subselect in data.table
	# sometimes it does not yield a vector
	# so i cannot use is.vector to check the consistency
	
	if(is.vector(x)){
		return(TRUE)
	} else {
		if(is.null(dim(x)) && !is.list(x)){
			return(TRUE)
		}
	}
	return(FALSE)
}













































































