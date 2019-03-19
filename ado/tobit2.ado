/* MIT License

Copyright 2008 by Jingling Guan and Mitchell A. Petersen

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


* tobit2.ado ---- written by Jingling Guan/Mitchell Petersen -- February 2008

* Program calculates clustered standard errors in both a firm and time
*  dimension for tobit models as described by Thompson in "A Simple
*  Formula for Standard Errors that Cluster by Both Firm and Time"
*  and and Cameron, Gelbach, and Miller, 2006, "Robust Inference
*  with Multi-way Clustering" Additions and edits 

* Compliant with outreg 

* Checks for multiple observations per fcluster-tcluster


#delimit ;

program define tobit2, eclass sortpreserve byable(recall);
        syntax [varlist] [in] [if], fcluster(varname)  tcluster(varname) ll(real); 
        tokenize `varlist';
        marksample touse;
        local depv `"`1'"';

* ----------------------------------------------------------------;
* ------ generate 2 dependent variable to run intreg  ------------;
* ----------------------------------------------------------------;

capture confirm existence `ll';
if !_rc {;
   scalar cutoff = `ll';
   qui gen y2a = `depv'; 
   qui replace y2a = . if `depv'<=cutoff;  
   qui gen y3a = `depv'; 
   qui replace y3a = cutoff if `depv'<=cutoff; 
   };
else {;
     di as err "Please specify lower  limit using ll( )"; };   
/*
else{;
	capture confirm existence `ul';
	if !_rc {;
	 	scalar cutoff = `ul';
		qui gen y2a = `depv';
		qui replace y2a = cutoff if `depv'>=cutoff;   
		qui gen y3a = `depv';
		qui replace y3a = . if `depv'>=cutoff;
	                };
	else {;
		di as err "Please specify either ll or ul";
		exit 111;
	   		};
	};
*/	


* ---------------------------------------------------------------- ;
* ------ generate the new variables list to run intreg ----------- ;
* ---------------------------------------------------------------- ;
tokenize `varlist';
macro shift;

qui reg `varlist';
scalar numindvar = e(df_m);

local varlist0 y2a y3a;

local i=1;

while `i' <= numindvar {;
   local varlistnew `varlist0' ``i'';
   local varlist0 `varlistnew';
   local i=`i'+1;
   };


* ---------------------------------------------------------------- ;
* ------- Tobit Clustering by First Variable (e.g. Firm) --------- ;
* ---------------------------------------------------------------- ;
	quietly intreg `varlistnew' if `touse', robust cluster(`fcluster');
	matrix vcf = e(V);
	local nfcluster=e(N_clust);

* ---------------------------------------------------------------- ;
* -------- Tobit Clustering by Second Variable (e.g. Time) ------- ;
* ---------------------------------------------------------------- ;
	quietly intreg `varlistnew' if `touse', robust cluster(`tcluster');
	matrix vct = e(V);
	local ntcluster=e(N_clust);
* ---------------------------------------------------------------- ;
* ------------------  Tobit with "No Clustering"  ---------------- ;
* ---------------------------------------------------------------- ;
	capture confirm string variable `fcluster';
	if !_rc {;
		gen bc1 = `fcluster'; /* string variable */
		};
		else {;
		gen bc1 = string(`fcluster'); /* numeric */
		};
	capture confirm string variable `tcluster';
	if !_rc {;
		gen bc2 = `tcluster'; /* string variable */
		};
		else {;
		gen bc2 = string(`tcluster'); /* numeric */
		};
	gen bc3 = bc1 + "_" + bc2;
	* --------------------------------------------------------- ;
	*   Check for multiple observations per fcluster-tcluster   ;
	* --------------------------------------------------------- ;
	bysort bc3: gen unique_obs = _n==1;	* =1 if only one obs per fcluster-tcluster;
	qui sum unique_obs;		

	if r(mean)==1 {;
	   	quietly intreg `varlistnew' if `touse', robust;
	   	local mcluster=0;
		};
	   else {;
	   	quietly intreg `varlistnew' if `touse', robust cluster(bc3);
	   	local mcluster =1 ;
		};

	drop bc1 bc2 bc3 unique_obs;

	local nparm = e(df_m)+1;
	matrix coef = e(b);

	matrix vwhite = e(V);

	matrix vc = vcf + vct - vwhite;
*	matrix vc1 = vc[1..numindvar+1,1..numindvar+1];
	
* ---------------------------------------------------------------- ;
* ------------------- Print out Tobit Results -------------------- ;
* ---------------------------------------------------------------- ;
	tokenize `varlist';  		/* this puts varlist in to the macros `1' `2' etc */
	macro shift;			/* drops first arguement (dep var) and shifts the rest up one */
	
	dis " ";
	dis in green "Tobit with 2D clustered SEs"
		_column (56) "Number of obs = " %7.0f in yellow e(N);
*	dis in green _column(56) "F(" %3.0f e(df_m) "," %6.0f e(df_r) ") =" %8.2f in yellow e(F);
*	dis in green _column(56) "Prob > F      ="  %8.4f in yellow 1-F(e(df_m),e(df_r),e(F));
	dis in green "Number of clusters (`fcluster') = " _column(31) %5.0f in yellow $_nfcluster;
*          in green _column(56) "R-squared     =" %8.4f in yellow e(r2);
   	dis in green "Number of clusters (`tcluster') = " _column(31) %5.0f in yellow $_ntcluster;
*	    in green _column(56) "Root MSE      =" %8.4f in yellow e(rmse);

* ---------------------------------------------------------------- ;
* -------------------- upload tobit Results into e()-------------- ;
* ---------------------------------------------------------------- ;

* save statistics from the last tobit (clustered by fcluster+tcluster);
* scalars;
	scalar e_N=e(N);
	scalar e_df_m = e(df_m);
	scalar e_df_r = e(df_r);

	scalar e_ll = e(ll);
	scalar e_ll_0 = e(ll_0);
	
	scalar e_N_unc = e(N_unc);
	scalar e_N_lc = e(N_lc);
	scalar e_N_rc = e(N_rc);

* prepare matrices to upload into e();
	ereturn clear;
	tempname b V;
	matrix `b' = coef;
	matrix `V' = vc;

* post the resuls in e();
	ereturn post `b' `V';
	ereturn scalar N = e_N;
	ereturn scalar df_m = e_df_m;
	ereturn scalar df_r = e_df_r;
 
	ereturn scalar ll = e_ll;
	ereturn scalar ll_0 = e_ll_0;
	ereturn local title "Tobit with clustered SEs";
	ereturn local method "2-dimension clustered SEs";
	ereturn local depvar "`depv'";
	ereturn local cmd "cluster2";

	ereturn scalar N_lc = e_N_lc;
	ereturn scalar N_rc = e_N_rc;
	ereturn scalar llopt = cutoff;

* end of uploading;

drop y2a y3a;
* ==================================================================;

* display coefficients and se;
	ereturn display;

	dis " ";
	if $_mcluster==1 {;
		dis "     SE clustered by " "`fcluster'" " and " "`tcluster'" " (multiple obs per " "`fcluster'" "-" "`tcluster'" ")";
		};
	   else {;
		dis "     SE clustered by " "`fcluster'" " and " "`tcluster'";
		};		
	dis " "; 

	scalar drop e_N e_df_m e_df_r e_ll e_ll_0 e_N_lc e_N_rc e_N_unc cutoff;
	matrix drop coef vc vcf vct;

end;

