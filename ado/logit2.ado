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

* logit2.ado ---- written by Jingling Guan/Mitchell Petersen -- February 2008

* Program calculates clustered standard errors in both a firm and time
  dimension for lobit models as described by Thompson in "A Simple
  Formula for Standard Errors that Cluster by Both Firm and Time" and
  and Cameron, Gelbach, and Miller, 2006, "Robust Inference with
  Multi-way Clustering" Additions and edits Compliant with outreg
  Checks for multiple observations per fcluster-tcluster

*/

#delimit ;
program define logit2, eclass sortpreserve byable(recall);
	syntax [varlist] [in] [if], fcluster(varname) tcluster(varname);
	tokenize `varlist';
	marksample touse;
	local depv `"`1'"';
* ---------------------------------------------------------------- ;
* ------- Logit Clustering by First Variable (e.g. Firm) --------- ;
* ---------------------------------------------------------------- ;
	quietly logit `varlist' if `touse', robust cluster(`fcluster');
	matrix vcf = e(V);
	local nfcluster=e(N_clust);
* ---------------------------------------------------------------- ;
* -------- Logit Clustering by Second Variable (e.g. Time) ------- ;
* ---------------------------------------------------------------- ;
	quietly logit `varlist' if `touse', robust cluster(`tcluster');
	matrix vct = e(V);
	local ntcluster=e(N_clust);
* ---------------------------------------------------------------- ;
* ------------------  Logit with "No Clustering"  ---------------- ;
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
	   	quietly logit `varlist' if `touse', robust;
	   	local mcluster=0;
		};
	   else {;
	   	quietly logit `varlist' if `touse', robust cluster(bc3);
	   	local mcluster =1 ;
		};
	drop bc1 bc2 bc3 unique_obs;

	local nparm = e(df_m)+1;
	matrix coef = e(b);
	matrix vc = vcf+vct-e(V);
	
* ---------------------------------------------------------------- ;
* ------------------- Print out Logit Results -------------------- ;
* ---------------------------------------------------------------- ;
	tokenize `varlist';  		/* this puts varlist in to the macros `1' `2' etc */
	macro shift;			/* drops first arguement (dep var) and shifts the rest up one */
	
	dis " ";
	dis in green "Logit with 2D clustered SEs";
	dis in green _column (56) "Number of obs = " %7.0f in yellow e(N);
	dis in green "Number of clusters (`fcluster') = " _column(31) %5.0f in yellow $_nfcluster
            in green _column(56) "Wald chi2(" %3.0f e(df_m) ") =" %7.2f in yellow e(chi2);
	dis in green "Number of clusters (`tcluster') = " _column(31) %5.0f in yellow $_ntcluster
            in green _column(56) "Prob > chi2   =" %8.4f in yellow chi2tail(e(df_m),e(chi2));
   	dis in green "Log pseudolikelihood = " %10.5f in yellow e(ll)
	    in green _column(56) "Pseudo R2     =" %8.4f in yellow e(r2_p);

* ---------------------------------------------------------------- ;
* -------------------- upload Logit Results into e()-------------- ;
* ---------------------------------------------------------------- ;

* save statistics from the last logit (clustered by fcluster+tcluster);
* scalars;
	scalar e_N=e(N);
	scalar e_df_m = e(df_m);
	scalar e_df_r = e(df_r);
	scalar e_F = e(F);
	scalar e_r2 = e(r2);
	scalar e_rmse = e(rmse);
	scalar e_mss = e(mss);
	scalar e_rss = e(rss);
	scalar e_r2_a = e(r2_a);
	scalar e_ll = e(ll);
	scalar e_ll_0 = e(ll_0);

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
 	ereturn scalar F= e_F;
	ereturn scalar r2= e_r2;
	ereturn scalar rmse = e_rmse;
	ereturn scalar mss = e_mss;
	ereturn scalar rss = e_rss;
	ereturn scalar r2_a = e_r2_a;
	ereturn scalar ll = e_ll;
	ereturn scalar ll_0 = e_ll_0;
	ereturn local title "Logit with clustered SEs";
	ereturn local method "2-dimension clustered SEs";
	ereturn local depvar "`depv'";
	ereturn local cmd "cluster2";
* end of uploading;
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

	scalar drop e_N e_df_m e_df_r e_F e_r2 e_rmse e_mss e_rss e_r2_a e_ll e_ll_0;
	matrix drop coef vc vcf vct;

end;
