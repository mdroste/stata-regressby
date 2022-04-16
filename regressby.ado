*===============================================================================
* PROGRAM: regressby.ado
* PURPOSE: Performs fast grouped univariate OLS regressions.
*          The following commands are equivalent:
*			 regressby y x, by(byvars)
*			 statsby, by(byvars) clear: reg y x	
*		   Except regressby will run 10-100x faster.
*          Also computes standard errors in a variety of flavors: usual
*          asymptotic standard errors, robust standard errors, and clustered
*          standard errors.
* AUTHORS: Michael Stepner, Michael Droste, Wilbur Townsend
*===============================================================================


*-------------------------------------------------------------------------------
* Stata wrapper
*-------------------------------------------------------------------------------

program define regressby

	version 12.0
	syntax varlist(min=2 numeric) [aweight], by(varlist) [vce(string) covs save(string)] 
	
	* Preserve dataset in case we crash
	preserve

	* Restrict sample with if/in conditions
	marksample touse, strok novarlist
	qui drop if `touse'==0

	* Parse VCE option, if specified
	if `"`vce'"' != "" {
		my_vce_parse , vce(`vce') 
		local vcetype    "robust"
		local clusterby  "`r(clustervar)'"
		if "`vcetype'"=="robust" local robust "robust"
		if "`clusterby'"!="" local robust = ""
	}

	* Check to make sure save data file path is valid
	if ("`replace'"=="") & (`"`savegraph'"'!="") {
		if regexm(`"`savegraph'"',"\.[a-zA-Z0-9]+$") confirm new file `"`save'"'
		else confirm new file `"`save'.dta"'
	}

	* Error checking: can't specify both robust and clusterby
	if "`robust'"!="" & "`clusterby'"!="" {
		di as error "Error: can't specify both clustered and robust standard errors at once! Choose one."
		exit
	}
	
	* Display type of standard error chosen
	if "`robust'"=="" & "`clusterby'"=="" {
		di "Running regressby with normal OLS standard errors."
	}
	if "`robust'"!="" {
		di "Running regressby with robust standard errors."
	}
	if "`clusterby'"!="" {
		di "Running regressby with cluster-robust standard errors (clustered by `clusterby')."
	}
	
	* Construct analytical weight variable
	if ("`weight'"!="") {
		local wt [`weight'`exp']
		tempvar tmpwt
		gen `tmpwt' `exp'
		local weightby `tmpwt'
		di "Using analytical weights, weight `exp'."
	}

	* Display weighting scheme, if applicable
	if "`weightby'"!="" {
		foreach v in `varlist' {
			qui replace `v' = `v' * sqrt(`weightby')
		}
		qui replace `weightby' = sqrt(`weightby')
	}
	
	* Convert string by-vars to temporary numeric variables
	foreach var of varlist `by' {
		cap confirm numeric variable `var', exact
		if _rc==0 {  // numeric var
			local bynumeric `bynumeric' `var'
		}
		else {  // string var
			tempvar `var'N
			encode `var', gen(``var'N')
			local bynumeric `bynumeric' ``var'N'
			local bystr `bystr' `var'  // list of string by-vars
		}
	}
	
	* Sort using by-groups
	sort `by' `clusterby'
	
	* Generate a single by-variable counting by groups
	tempvar grp
	egen `grp'=group(`bynumeric')
	qui drop if mi(`grp')
	
	* Drop observations missing independent or dependent variables
	* Also count number of variables here including constant, awkward and should be replaced
	local num_x = 0
	foreach v in `varlist'{
		qui drop if mi(`v')
		local num_x = `num_x' + 1
	}
	local num_x = `num_x' - 1
	if "`nocons'"=="" local num_x = `num_x' + 1

	* Drop observations missing weight, if weights are specified
	if "`weightby'"!="" {
		drop if `weightby'==.
	}

	* XX revisit this later to handle missing data
	

	* Perform regressions on each by-group, store in dataset
	mata: _regressby("`varlist'", "`grp'", "`bynumeric'","`clusterby'","`robust'","`weightby'")
	
	* Convert string by-vars back to strings, from numeric
	foreach var in `bystr' {
		decode ``var'N', gen(`var')
	}
	order `by'


	* XX find out if it is faster to compute R2 in Mata or Stata
	if "`nocov'"!="" {
		cap drop _cov_*
	}

	* XX optionally save out to dta and just restore with a message
	if "`save'"=="" {
		restore, not
	}
	if "`save'"!="" {
		save `save', replace
		restore
	}

end

*-------------------------------------------------------------------------------
* Mata program: _regressby3
* Inputs:
*  	- A y-var and x-var for an OLS regression
*  	- A group var, for which each value represents a distinct by-group.	
*		This var must be in ascending order.
*	- A list of numeric by-variables, whose groups correspond to th group var.
* Outputs:
*  	- dataset of coefficients from OLS regression for each by-group
*-------------------------------------------------------------------------------

version 13.1
set matastrict on

mata:
void _regressby(string scalar regvars, string scalar grpvar, string scalar byvars, string scalar clusterby, string scalar robust, string scalar weightby) {

// Convert variable names to column indices
real rowvector regcols, bycols, clustercol, weightcol
real scalar ycol, xcol, grpcol
regcols 	= st_varindex(tokens(regvars))
bycols 		= st_varindex(tokens(byvars))
clustercol 	= st_varindex(tokens(clusterby))
weightcol 	= st_varindex(tokens(weightby))
grpcol 		= st_varindex(grpvar)

// Fetch number of groups
real scalar numgrp, startobs, curgrp
numgrp 		= _st_data(st_nobs(),grpcol)
startobs 	= 1  
curgrp	 	= _st_data(1,grpcol)

// Preallocate matrices for output
real matrix groups, coefs, ses, covs, nobs
groups 		= J(numgrp, cols(bycols), .)
coefs 		= J(numgrp, cols(regcols), .)
Vs          = J(numgrp, cols(regcols)^2, .)
nobs 		= J(numgrp, 1, .)

// Preallocate regression objects
real matrix XX, Xy, XX_inv, V, Z, M, y, x, w
real scalar N, k, cov, p, nc
real vector beta, e, s2, cvar, xi, ei

// -----------------------------------------------------------------------------
// Iterate over groups
// -----------------------------------------------------------------------------

// Iterate over groups 1 to Ng-1
for (obs=1; obs<=st_nobs()-1; obs++) {
	if (_st_data(obs,grpcol)!=curgrp) {
		st_view(M, (startobs,obs-1), regcols, 0)
		st_subview(y, M, ., 1)
		st_subview(X, M, ., (2\.))
		N    = rows(X)
	// Augment x with either column of 1's or weights
	// TODO -- noconstant option needs to be specified here and also accounted for in df
		if (weightby!="") {
			st_view(w, (startobs,obs-1), weightcol, 0)
			X = X,w
		}
		if (weightby=="") {
			X = X,J(N,1,1)
		}
		// Define matrix products
		XX 		= quadcross(X,X)
		Xy 		= quadcross(X,y)
		XX_inv 	= invsym(XX)
		// ------------ COMPUTE COEFFICIENTS --------------------
		beta 	= (XX_inv*Xy)'
        e 		= y - X*beta'
		p    	= cols(X)
		k    	= p - diag0cnt(XX_inv)
		// ------------ COMPUTE STANDARD ERRORS -----------------
		if (robust == "" & clusterby=="") {
			V 	= quadcross(e,e)/(N-k)*cholinv(XX)
		}
		if (robust != "") {
			V   = (N/(N-k))*XX_inv*quadcross(X, e:^2, X)*XX_inv
		}
		if (clusterby != "") {
			st_view(cvar,(startobs,obs-1),clustercol,0)
			info = panelsetup(cvar, 1)
			nc  = rows(info)
			Z   = J(k, k, 0)
			if (nc>2) {
				for (i=1; i<=nc; i++) {
					xi = panelsubmatrix(X,i,info)
					ei = panelsubmatrix(e,i,info)
					Z  = Z + xi'*(ei*ei')*xi
				}
				V   = ((N-1)/(N-k))*(nc/(nc-1))*XX_inv*Z*XX_inv
			}
		}
		// ------------ STORE OUTPUT ----------------------------
		coefs[curgrp,.] 	= beta
	    Vs[curgrp,.]        = rowshape(V, 1)
		nobs[curgrp,1]  	= N
		groups[curgrp,.] 	= st_data(startobs,bycols)
		// ------------ WRAP UP BY ITERATING COUNTERS -----------
		curgrp	 = _st_data(obs,grpcol)
		startobs = obs
	}
}

// Iterate over last group manually
obs=st_nobs()
if (_st_data(obs,grpcol)==curgrp) {  // last observation is not a group to itself
	// increment obs, since code is written as processing the observation that is 1 past the last in the group
	++obs
	// compute OLS coefs: beta = inv(X'X) * X'y. --> see Example 4 of -help mf_cross-
	st_view(M, (startobs,obs-1), regcols, 0)
	st_subview(y, M, ., 1)
	st_subview(X, M, ., (2\.))
	N    = rows(X)
	// Augment X with either column of 1's (unweighted) or weights (weighted)
	// TODO -- noconstant option needs to be specified here and also accounted for in df
	if (weightby!="") {
		st_view(w, (startobs,obs-1), weightcol, 0)
		X = X,w
	}
	if (weightby=="") {
		X = X,J(N,1,1)
	}
	// Define matrix products
	XX 		= quadcross(X,X)
	Xy 		= quadcross(X,y)
	XX_inv 	= invsym(XX)
	beta 	= (XX_inv*Xy)'
    e 		= y - X*beta'
	p    	= cols(X)
	k    	= p - diag0cnt(XX_inv)
	// USUAL OLS STANDARD ERRORS
	if (robust == "" & clusterby == "") {
		V 	= quadcross(e,e)/(N-k)*cholinv(XX)
	}
	// ROBUST STANDARD ERRORS
	if (robust != "") {
		V    = (N/(N-k))*XX_inv*quadcross(X, e:^2, X)*XX_inv
	}
	// CLUSTERED STANDARD ERRORS
	if (clusterby != "") {
		st_view(cvar,(startobs,obs-1),clustercol,0)
		info = panelsetup(cvar, 1)
		nc  = rows(info)
		Z   = J(k, k, 0)
		if (nc>2) {
			for (i=1; i<=nc; i++) {
				xi = panelsubmatrix(X,i,info)
				ei = panelsubmatrix(e,i,info)
				Z  = Z + xi'*(ei*ei')*xi
			}
			V   = ((N-1)/(N-k))*(nc/(nc-1))*XX_inv*Z*XX_inv
		}
	}
	// STORE REGRESSION OUTPUT
	coefs[curgrp,.] 	= beta
    Vs[curgrp,.]        = rowshape(V, 1)
	nobs[curgrp,1]  	= N
	groups[curgrp,.] 	= st_data(startobs,bycols)
}
	
else {
	display("{error} last observation is in a singleton group")
	exit(2001)
}

// -----------------------------------------------------------------------------
// Gather output and pass back into Stata
// -----------------------------------------------------------------------------

// Store group identifiers in dataset
stata("qui keep in 1/"+strofreal(numgrp, "%18.0g"))
stata("keep "+byvars)
st_store(.,tokens(byvars),groups)

// Store coefficients in dataset:

// ... Number of observations,
(void) st_addvar("long", "N")
st_store(., ("N"), nobs)

// ... And then looping over covariates,
covariates = (cols(regcols)>1) ? tokens(regvars)[|2 \ .|], "cons" : ("cons")
for (k=1; k<=length(covariates); k++) {
    covName = covariates[k]
    // ... Coefficients and standard errors,
    (void) st_addvar("float", "_b_"+covName)
    (void) st_addvar("float", "_se_"+covName)
    st_store(., "_b_"+covName,  coefs[., k])
    st_store(., "_se_"+covName, sqrt(Vs[., k + cols(regcols)*(k - 1)]))
    // ... And the sampling covariances.
    for (j=1; j<k; j++) {
        otherCovName = covariates[j]
        (void) st_addvar("float", "_cov_"+covName+"_"+otherCovName)
        st_store(., "_cov_"+covName+"_"+otherCovName, Vs[., k + cols(regcols) * (j - 1)])
    }
}

}

end

*-------------------------------------------------------------------------------
* Auxiliary Stata programs for parsing vce (standard error options)
* Source: https://blog.stata.com/2015/12/08/programming-an-estimation-command-in-stata-using-a-subroutine-to-parse-a-complex-option/
*-------------------------------------------------------------------------------

program define my_vce_parse, rclass
    syntax  [, vce(string) ]
    local case : word count `vce'
    if `case' > 2 {
        my_vce_error , typed(`vce')
    }
    local 0 `", `vce'"' 
    syntax  [, Robust CLuster * ]
    if `case' == 2 {
        if "`robust'" == "robust" | "`cluster'" == "" {
            my_vce_error , typed(`vce')
        }
        capture confirm numeric variable `options'
        if _rc {
            my_vce_error , typed(`vce')
        }
        local clustervar "`options'" 
    }
    else {    // case = 1
        if "`robust'" == "" {
            my_vce_error , typed(`vce')
        }
    }
    return clear    
    return local clustervar "`clustervar'" 
end
 
program define my_vce_error
    syntax , typed(string)
    display `"{red}{bf:vce(`typed')} invalid"'
    error 498
end
