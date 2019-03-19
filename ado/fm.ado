/* MIT License

Copyright 2005 Mitchell A. Petersen

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*  Fama-MacBeth estimation of coefficients and standard errors
  run N time series regressions

* Send questions and corrections to: mpetersen@northwestern.edu
 Program takes list of variables as inputs (first is dep variable,
 rest are independent)

* The program requires the firm and time * variable to be set w/ tsset
  command (e.g. tsset firmid timeid).* Additions and edits

* Compliant with outreg -- June, 2007 (Jingling Guan) [eclass]

* Allows specification of byfm variable (e.g. if year then run
  regression for each of T years) June, 2007 (Jingling Guan)

*  Macros: 
	We can use ` ' or $_, they mean the same thing
	Thus `i' and $_i are the same and `b1' and $_b1 are the same

*/

#delimit ;
           
program define fm, eclass sortpreserve byable(recall);
	* version 8;

	syntax varlist [if] [in], [ timeseries(integer 0) ] [byfm(varname)];	/* if specify regpanel=1, then run N time series regressions */

	quietly capture tsset;			/* Checks that panel variable has been set by tsset */
		if _rc ~=0 {;
			di as err "panel variable not set, use -tsset panelvar timevar ...-";
			exit 111;
			};

	local panel_ = r(panelvar);		/* capture the name of the panel variable from tsset */
	local time_ = r(timevar);		/* capture the name of the time series variable from tsset */

   	capture confirm variable `byfm';
		if !_rc {;
		   	gen regvar_ = `byfm';	 /* by variable has been specified */
	                };
		   else {;
			gen regvar_ = `r(timevar)';
	   		};

	quietly save junk_fm, replace;		/* This saves the raw data set for later use */

	* ---------------  generate matrices coef and vc of the right size to -----------;
	* ----------------- save the coefficients and var-covar matrices ----------------;
	qui reg `varlist' `in' `if';
	matrix coef = e(b);	
	matrix vc = e(V);	
	local nobs_=e(N);

	* --------------- do N/T time series/cross sectional regressions ----------------;
	* ---------- and save the coefficients for each regression in variables b_ ------;
	quietly statsby "reg `varlist' `in' `if'" _b iobs_=e(N), by(regvar_);  * replaces data in memory;


	* ----------------- Calculate the coefficients for indep vars -------------------;	
	tokenize `varlist';  			/* this puts varlist in to the macros `1' `2' etc */
	macro shift;				/* drops first arguement and shifts the rest up one */
	local i = 1;
	* loop as long as macro `1' is not empty;
	while "`1'" ~= "" {;
		quietly sum b_`1';
		local b`i' = r(mean);
		matrix coef[1,`i']= r(mean);			/* save the slope coefficients in matrix coef */
		matrix vc[`i',`i']= (r(sd)/sqrt(r(N)))^2;	/* save the variances in matrix vc */
		macro shift;
		local i = `i' + 1;
		};

	* ------------------------ Now calculate the constant ----------------------------;
	quietly sum b_cons;
	local bcons = r(mean);
	local nreg_ = r(N);			/* nreg_ saves the number of regressions run */

	* save the slope coefficient in matrix coef;
	matrix coef[1,`i']= r(mean);				/* save the slope coefficients in matrix coef */
	matrix vc[`i',`i']= (r(sd)/sqrt(r(N)))^2;		/* save the variances in matrix vc */

	drop _all;
	use junk_fm;
	erase junk_fm.dta;		
	

	* ---------------------------------------------------------------------------------;
	* ------------------ Calculate R2 based on the FM coefficient estimates -----------;
	* ---------------------------------------------------------------------------------;
	gen yhat_ = $_bcons;
	tokenize `varlist';
	macro shift;
	local i = 1;
	while "`1'" ~= ""{;
		quietly replace yhat_ = yhat_ + `1' * `b`i'';
		macro shift;
		local i = `i' + 1;
		};	

	tokenize `varlist'; 			* puts dependent variable in first arguement;
	quietly reg `1' yhat_ `if' `in';
	local e_r2 = e(r2);
	* --------------------- end of R2 calculation -------------------------------------;

	* ---------------------------------------------------------------------------------;
	* --------------------- upload the statistics to e() ------------------------------;
	* ---------------------------------------------------------------------------------;
	qui reg `varlist' `in' `if', robust;
	tempname b V;
	matrix `b' = coef;
	matrix `V' = vc;
	ereturn post `b' `V'; 		* post coef and var-covariance matrix V;

	* scalars;
	ereturn scalar N = $_nobs_;
	ereturn scalar r2 = $_e_r2;
	ereturn scalar N_g = $_nreg_;
	ereturn scalar df_m = wordcount("`varlist'")-1;
	ereturn scalar df_r = $_nreg_;

	ereturn local title "Fama-MacBeth Estimation";
	ereturn local method "Fama-MacBeth";
*	ereturn local vcetype "Fama-MacBeth";
	ereturn local cmd "fm";
	ereturn local depvar "`1'";
	* ------------------------------ end of uploading ---------------------------------;

	* ---------------------------------------------------------------------------------;
	* ------------------------ Print out regression results ---------------------------;
	* ---------------------------------------------------------------------------------; 
	dis " ";
	dis in green e(title);	
	
	
	
	dis " ";
	dis in green "panel variable:  " in yellow "`panel_'"
	    in green _column(48) "Number of obs" _column(70) "= " %7.0f in yellow e(N);

   	capture confirm variable `byfm';
		if !_rc {;
			dis in green " time variable:  " in yellow "`time_'"
			    in green _column(48) "Number of " "`byfm'" "(s)" _column(70) "= " %7.0f in yellow $_nreg_;
	                };
		   else {;
			dis in green " time variable:  " in yellow "`time_'"
			    in green _column(48) "Number of " "`time_'" "(s)" _column(70) "= " %7.0f in yellow $_nreg_;
	   		};
	dis in green _column(48) "R-squared             = " %7.4f in yellow e(r2);
	
	ereturn display;

  	capture confirm variable `byfm';
		if !_rc {;
			dis _column(10) in yellow %7.0f $_nreg_ in green " " "`byfm'" " regressions";
	                };
		   else {;
			dis _column(10) in yellow %7.0f $_nreg_ in green " " "`time_'" " regressions";
	   		};
	
	* ------------------------------ end of printing ----------------------------------;

	drop yhat_ regvar_;
	
	end;

