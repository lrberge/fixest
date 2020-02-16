
***********************
** STATA COMPARISONS **
***********************

/*
OLS:
reghdfe
*/

/*
POISSON:
ppmlhdfe
*/

/*
NEGBIN:
nbreg (in lack of FE methods)
*/

/*
LOGIT:
logit (in lack of FE methods)
*/


// ssc install mat2txt
// ssc install reghdfe, replace
// ssc install ppmlhdfe, replace

// reghdfe, compile

  /**********************/
 /* General information*/
/**********************/

* set the right path
cd "_STATA"

set more off, perm
set matsize 11000
set maxvar 100000


  /*************/
 /**** OLS ****/
/*************/


*******************
***** REGHDFE *****
*******************

*
* 1 cluster
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/5{
	if `size' == 5 {
		insheet using "../DATA/base_10M.csv", clear
		rename x1 X1
	}

	forvalues rep = 1/10{

		if `size' < 5 {
			use "base_s`size'_r`rep'.dta", clear
		}

		timer clear 1

		timer on 1
		reghdfe ln_y X1, a(dum_1)
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}


** SAVE **
mat2txt , matrix(A) saving("reghdfe_G1.txt") replace

*
* 2 clusters
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/5{
	if `size' == 5 {
		insheet using "../DATA/base_10M.csv", clear
		rename x1 X1
	}
	forvalues rep = 1/10{

		if `size' < 5 {
			use "base_s`size'_r`rep'.dta", clear
		}

		timer clear 1

		timer on 1
		reghdfe ln_y X1, a(dum_1 dum_2)
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}


** SAVE **
mat2txt , matrix(A) saving("reghdfe_G2.txt") replace

*
* 3 clusters
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/5{
	if `size' == 5 {
		insheet using "../DATA/base_10M.csv", clear
		rename x1 X1
	}
	forvalues rep = 1/10{

		if `size' < 5 {
			use "base_s`size'_r`rep'.dta", clear
		}

		timer clear 1

		timer on 1
		reghdfe ln_y X1, a(dum_1 dum_2 dum_3)
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("reghdfe_G3.txt") replace


***            *** 
* Difficult data *
***            ***

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 3/6{

	use "base_diff_e`size'.dta", clear
	
	forvalues g = 1/3{
	
		forvalues rep = 1/10{

			timer clear 1

			timer on 1
			if `g' == 1{
				reghdfe y x1 x2, a(id_indiv)
			} 
			else if `g' == 2{
				reghdfe y x1 x2, a(id_indiv id_firm)
			} 
			else {
				reghdfe y x1 x2, a(id_indiv id_firm id_year)
			}
			timer off 1

			// to create the value r(t1)
			timer list 1
			matrix A = A \ (`size', r(t1))
		}
	}
}

** SAVE **
mat2txt , matrix(A) saving("reghdfe_diff.txt") replace



  /***********/
 /* POISSON */
/***********/



*********************
***** XTPOISSON *****
*********************

*
* 1 cluster
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/4{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1
		xtset dum_1

		timer on 1
		xtpoisson y X1, fe
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("xtpoisson_G1.txt") replace

*
* 2 clusters
*

* OK: s1, s2 , s3 // s4: just one estimation exceeds 5600s

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/3{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1
		xtset dum_1

		timer on 1
		xtpoisson y X1 i.dum_2, fe
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("xtpoisson_G2.txt") replace

*
* 3 clusters
*

*OK: s1, s2, s3 // s4 => circa 1000s

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/3{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1
		xtset dum_1

		timer on 1
		xtpoisson y X1 i.dum_2 i.dum_3, fe
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("xtpoisson_G3.txt") replace


********************
***** PPMLHDFE *****
********************

*
* 1 cluster
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/4{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1

		timer on 1
		ppmlhdfe y X1, a(dum_1)
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("ppmlhdfe_G1.txt") replace

*
* 2 clusters
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/4{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1

		timer on 1
		ppmlhdfe y X1, a(dum_1 dum_2)
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("ppmlhdfe_G2.txt") replace

*
* 3 clusters
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/4{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1

		timer on 1
		ppmlhdfe y X1, a(dum_1 dum_2 dum_3)
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("ppmlhdfe_G3.txt") replace



  /*********************/
 /* NEGATIVE BINOMIAL */
/*********************/

*****************
***** NBREG *****
*****************

* We cannot take advantage of the xtnbreg because it is conditional, thus a different model

*
* 1 cluster
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/2{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1

		bysort dum_1: egen sum_y = sum(y)
		drop if sum_y == 0

		timer on 1
		nbreg y X1 i.dum_1, diff
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("nbreg_G1.txt") replace

*
* 2 clusters
*

* OK: s1, s2 // s3 is doable but very long => to do

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/2{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1

		* dropping the problems
		bysort dum_1: egen sum_y_1 = sum(y)
		drop if sum_y_1 == 0

		bysort dum_2: egen sum_y_2 = sum(y)
		drop if sum_y_2 == 0

		timer on 1
		nbreg y X1 i.dum_1 i.dum_2, diff
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("nbreg_G2.txt") replace

*
* 3 clusters
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/2{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1

		* dropping the problems
		bysort dum_1: egen sum_y_1 = sum(y)
		drop if sum_y_1 == 0

		bysort dum_2: egen sum_y_2 = sum(y)
		drop if sum_y_2 == 0

		bysort dum_3: egen sum_y_3 = sum(y)
		drop if sum_y_3 == 0

		timer on 1
		nbreg y X1 i.dum_1 i.dum_2 i.dum_3, diff
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("nbreg_G3.txt") replace



  /*******************/
 /*      LOGIT      */
/*******************/

*****************
***** LOGIT *****
*****************

** xtlogit is a conditionnal estimator
* we clean only 1/0 values

*
* 1 cluster
*

* initialisation of the matrix
matrix A = (0, 0)

* size 3: 8057s

forvalues size = 1/2{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1
		gen sign_y = 1 if y>0
		replace sign_y = 0 if y ==0

		* dropping the problems
		bysort dum_1: egen sum_y_1 = sum(sign_y)
		drop if sum_y_1 == 0
		bysort dum_1: egen sum_y1_1 = sum(1-sign_y)
		drop if sum_y1_1 == 0


		timer on 1
		logit sign_y X1 i.dum_1, diff
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("logit_G1.txt") replace

*
* 2 clusters
*

* OK: s1, s2 // s3 is doable but very long => to do

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/2{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1
		gen sign_y = 1 if y>0
		replace sign_y = 0 if y ==0

		* dropping the problems
		bysort dum_1: egen sum_y_1 = sum(sign_y)
		drop if sum_y_1 == 0
		bysort dum_1: egen sum_y1_1 = sum(1-sign_y)
		drop if sum_y1_1 == 0
		bysort dum_2: egen sum_y_2 = sum(sign_y)
		drop if sum_y_2 == 0
		bysort dum_2: egen sum_y1_2 = sum(1-sign_y)
		drop if sum_y1_2 == 0

		timer on 1
		logit sign_y X1 i.dum_1 i.dum_2, diff
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("logit_G2.txt") replace

*
* 3 clusters
*

* initialisation of the matrix
matrix A = (0, 0)

forvalues size = 1/2{
	forvalues rep = 1/10{
		use "base_s`size'_r`rep'.dta", clear

		timer clear 1
		gen sign_y = 1 if y>0
		replace sign_y = 0 if y ==0

		* dropping the problems
		bysort dum_1: egen sum_y_1 = sum(sign_y)
		drop if sum_y_1 == 0
		bysort dum_1: egen sum_y1_1 = sum(1-sign_y)
		drop if sum_y1_1 == 0
		bysort dum_2: egen sum_y_2 = sum(sign_y)
		drop if sum_y_2 == 0
		bysort dum_2: egen sum_y1_2 = sum(1-sign_y)
		drop if sum_y1_2 == 0
		bysort dum_3: egen sum_y_3 = sum(sign_y)
		drop if sum_y_3 == 0
		bysort dum_3: egen sum_y1_3 = sum(1-sign_y)
		drop if sum_y1_3 == 0

		timer on 1
		logit sign_y X1 i.dum_1 i.dum_2 i.dum_3, diff
		timer off 1

		// to create the value r(t1)
		timer list 1
		matrix A = A \ (`size', r(t1))
	}
}

** SAVE **
mat2txt , matrix(A) saving("logit_G3.txt") replace


