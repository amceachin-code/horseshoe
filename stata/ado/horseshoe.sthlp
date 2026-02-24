{smcl}
{* *! version 1.1.0  21feb2026}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[BAYES] bayes" "help bayes"}{...}
{viewerjumpto "Syntax" "horseshoe##syntax"}{...}
{viewerjumpto "Description" "horseshoe##description"}{...}
{viewerjumpto "Convergence" "horseshoe##convergence"}{...}
{viewerjumpto "Reproducibility" "horseshoe##reproducibility"}{...}
{viewerjumpto "Options" "horseshoe##options"}{...}
{viewerjumpto "Postestimation" "horseshoe##postestimation"}{...}
{viewerjumpto "Stored results" "horseshoe##results"}{...}
{viewerjumpto "Examples" "horseshoe##examples"}{...}
{viewerjumpto "References" "horseshoe##references"}{...}
{viewerjumpto "Author" "horseshoe##author"}{...}

{title:Title}

{p2colset 5 20 22 2}{...}
{p2col:{cmd:horseshoe} {hline 2}}Bayesian linear regression with the horseshoe prior{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:horseshoe}
{depvar}
{indepvars}
{ifin}
[{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opth pen:alized(matname)}}1 x {it:p} matrix of 0/1 penalization flags;
default is all 1s (all penalized){p_end}

{syntab:MCMC}
{synopt:{opt nmcmc(#)}}number of post-burnin MCMC draws to store; default is
{cmd:nmcmc(1000)}{p_end}
{synopt:{opt burn:in(#)}}number of burnin iterations to discard; default is
{cmd:burnin(500)}{p_end}
{synopt:{opt thin(#)}}thinning interval; default is {cmd:thin(1)}{p_end}
{synopt:{opt seed(#|string)}}set RNG seed for reproducibility{p_end}

{syntab:Prior}
{synopt:{opt lambdas:cale(#)}}scale of the half-Cauchy prior on local shrinkage
lambda_j; default is {cmd:lambdascale(1)}{p_end}
{synopt:{opt taus:cale(#)}}scale of the half-Cauchy prior on global shrinkage
tau; default is {cmd:tauscale(1)}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set credible interval level; default is
{cmd:level(95)}{p_end}
{synopt:{opt v:erbose}}print progress every 500 iterations{p_end}

{syntab:Output}
{synopt:{opt sav:ing(filename)}}save full posterior draws to {it:filename}.dta{p_end}
{synopt:{opt replace}}overwrite existing {opt saving()} file{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:horseshoe} fits a Bayesian linear regression model using the horseshoe
prior (Carvalho, Polson, and Scott 2010), estimated via the auxiliary-variable
Gibbs sampler of Makalic and Schmidt (2015). The horseshoe prior is a
continuous shrinkage prior that adapts to sparsity: it aggressively shrinks
noise coefficients toward zero while allowing true signals to remain large.

{pstd}
The key feature of this implementation is {bf:per-parameter penalization
control} via the {opt penalized()} option. This lets you specify which
coefficients receive the horseshoe prior (shrinkage toward zero) and which
receive a flat prior (no shrinkage). This is particularly useful when you want
to regularize high-dimensional interaction terms while leaving treatment main
effects unpenalized.

{pstd}
The model is:

{p 8 8 2}
y = X * beta + epsilon,   epsilon ~ N(0, sigma^2)

{pstd}
For penalized coefficients:

{p 8 8 2}
beta_j | lambda_j, tau, sigma^2 ~ N(0, lambda_j^2 * tau^2 * sigma^2){break}
lambda_j ~ C+(0, lambda_scale){break}
tau ~ C+(0, tau_scale)

{pstd}
For unpenalized coefficients:

{p 8 8 2}
beta_j ~ flat (improper) prior

{pstd}
The half-Cauchy priors are represented via scale-mixture of inverse-gamma
distributions using auxiliary variables nu_j and xi, following Makalic and
Schmidt (2015). This yields fully conjugate conditional posteriors, enabling
efficient Gibbs sampling.

{pstd}
{bf:Important:} The outcome {it:depvar} should be centered (mean-subtracted)
if no intercept column is included among {it:indepvars}. If you include an
intercept column, leave it unpenalized by setting its entry in
{opt penalized()} to 0.

{pstd}
{bf:Propriety requirement:} When using unpenalized parameters, the number of
observations must exceed the number of unpenalized parameters (n > p_free),
and the corresponding columns of X must have full column rank.


{marker convergence}{...}
{title:Convergence assessment}

{pstd}
{cmd:horseshoe} automatically computes the effective sample size (ESS) for all
sampled parameters using the initial positive sequence estimator (Geyer 1992).
This is the same estimator used by Stan's {cmd:monitor()} and R's
{cmd:coda::effectiveSize()}.

{pstd}
ESS measures how many independent draws the MCMC chain is equivalent to after
accounting for autocorrelation. Guidelines:

{p 8 8 2}
{bf:ESS > 400} per parameter is generally adequate for reliable posterior
summaries (means, SDs, and 95% credible intervals).{break}
{bf:ESS 100-400} may be sufficient for posterior means but credible interval
tails may be imprecise.{break}
{bf:ESS < 100} suggests poor mixing. Consider increasing {opt nmcmc()} or
{opt burnin()}, or check for near-collinearity in the design matrix.

{pstd}
The ESS for tau^2 (global shrinkage) is often the bottleneck because global
parameters typically mix slowest. If {cmd:e(ess_tau2)} is much lower than
{cmd:e(ess_min)} for the betas, the chain may need more iterations.

{pstd}
{bf:Choosing MCMC settings:}

{p 8 8 2}
The defaults {cmd:nmcmc(1000)} {cmd:burnin(500)} are adequate for moderate p
(up to ~100 predictors).{break}
For p > 200, consider {cmd:nmcmc(2000)} {cmd:burnin(1000)}.{break}
Thinning ({cmd:thin()}) is usually unnecessary unless storing draws for
external analysis and disk space is a concern. It does not improve ESS per
unit of compute time.

{pstd}
A warning is displayed automatically if the minimum ESS across all parameters
falls below 100.


{marker reproducibility}{...}
{title:Reproducibility}

{pstd}
Use the {opt seed()} option to set the random number generator seed before
sampling. This ensures that results are exactly reproducible across runs:

{phang2}{cmd:. horseshoe y x1 x2 x3, seed(12345)}{p_end}

{pstd}
After estimation, {cmd:e(seed)} stores the seed you specified and
{cmd:e(rng_state)} stores the full RNG state captured just before the sampler
began. You can restore the exact same state by passing the RNG state string
back to {opt seed()}:

{phang2}{cmd:. local saved_state = e(rng_state)}{p_end}
{phang2}{cmd:. horseshoe y x1 x2 x3, seed(`saved_state')}{p_end}

{pstd}
This produces bitwise-identical results to the original run.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth penalized(matname)} specifies a 1 x {it:p} Stata matrix of 0/1 flags
controlling which coefficients receive the horseshoe prior. A value of 1 means
the coefficient is penalized (horseshoe shrinkage); a value of 0 means it
receives a flat prior (no shrinkage). The matrix must have exactly as many
columns as there are independent variables, and must contain only 0s and 1s.
If not specified, all coefficients are penalized (standard horseshoe).

{dlgtab:MCMC}

{phang}
{opt nmcmc(#)} sets the number of post-burnin MCMC draws to store. Larger
values give more precise posterior summaries but take longer. Default is 1000.

{phang}
{opt burnin(#)} sets the number of initial MCMC iterations to discard as
burnin. These draws are used to let the Markov chain converge to the
stationary distribution. Default is 500.

{phang}
{opt thin(#)} sets the thinning interval. Only every {it:thin}-th post-burnin
iteration is stored. Use thinning > 1 to reduce autocorrelation in the stored
chain at the cost of longer runtime. Default is 1 (no thinning).

{phang}
{opt seed(#|string)} sets the random number generator seed before sampling.
Accepts either a non-negative integer (passed to {cmd:set seed}) or a full RNG
state string (passed to {cmd:set rngstate}), as returned in {cmd:e(rng_state)}.
If not specified, sampling uses whatever RNG state is current.

{dlgtab:Prior}

{phang}
{opt lambdascale(#)} sets the scale of the half-Cauchy prior on each local
shrinkage parameter lambda_j. Larger values allow individual coefficients to
escape shrinkage more easily (wider local tails). The standard horseshoe uses
{cmd:lambdascale(1)}. Default is 1.

{phang}
{opt tauscale(#)} sets the scale of the half-Cauchy prior on the global
shrinkage parameter tau. Smaller values enforce more aggressive overall
shrinkage. The standard horseshoe uses {cmd:tauscale(1)}. Default is 1.

{dlgtab:Reporting}

{phang}
{opt level(#)} specifies the credible interval level, as a percentage.
The default is {cmd:level(95)}, producing 95% equal-tailed credible intervals
from the 2.5th and 97.5th posterior percentiles.

{phang}
{opt verbose} requests that the sampler print progress messages every 500
iterations, including the current values of tau and sigma. Also enables
warnings for Cholesky ridge adjustments and extreme lambda^2 values.

{dlgtab:Output}

{phang}
{opt saving(filename)} saves the full posterior draws to {it:filename}.dta.
The saved dataset contains {it:nmcmc} observations with variables:
{cmd:draw} (draw number), {cmd:sigma2} (posterior draw of sigma^2),
{cmd:tau2} (posterior draw of tau^2), and {cmd:b_}{it:varname} for each
coefficient (using the original variable name, prefixed with {cmd:b_}).
This allows custom posterior analysis beyond the summaries stored in {cmd:e()}.

{phang}
{opt replace} permits {opt saving()} to overwrite an existing file.


{marker postestimation}{...}
{title:Postestimation commands}

{pstd}
The following standard postestimation commands are available after
{cmd:horseshoe}:

{synoptset 20}{...}
{synopt:Command}Description{p_end}
{synoptline}
{synopt:{helpb predict}}fitted values, residuals, standard errors of prediction{p_end}
{synopt:{helpb test}}Wald tests of linear hypotheses{p_end}
{synopt:{helpb lincom}}linear combinations of coefficients{p_end}
{synopt:{helpb nlcom}}nonlinear combinations of coefficients{p_end}
{synopt:{helpb margins}}marginal means and marginal effects{p_end}
{synoptline}

{pstd}
These commands work because {cmd:horseshoe} posts both {cmd:e(b)} (posterior
means) and {cmd:e(V)} (posterior variance-covariance matrix from the MCMC
draws) to Stata's {cmd:e()} class.

{pstd}
{bf:Replay:} Running {cmd:horseshoe} without arguments after estimation
re-displays the results table (standard Stata behavior for e-class commands).

{dlgtab:predict}

{p 8 17 2}
{cmd:predict}
{newvar}
[{cmd:,} {it:statistic}]

{synoptset 20 tabbed}{...}
{synopt:{opt xb}}linear prediction X*b; the default{p_end}
{synopt:{opt r:esiduals}}residuals y - X*b{p_end}
{synopt:{opt stdp}}standard error of the linear prediction{p_end}
{synoptline}

{pstd}
{opt xb} computes the linear prediction using the posterior mean coefficients
in {cmd:e(b)}.

{pstd}
{opt residuals} computes observed minus fitted: {it:depvar} - xb.

{pstd}
{opt stdp} computes the standard error of the linear prediction using the
posterior variance-covariance matrix {cmd:e(V)}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:horseshoe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(p)}}number of predictors{p_end}
{synopt:{cmd:e(p_pen)}}number of penalized predictors{p_end}
{synopt:{cmd:e(nmcmc)}}number of stored MCMC draws{p_end}
{synopt:{cmd:e(sigma2)}}posterior mean of sigma^2{p_end}
{synopt:{cmd:e(tau2)}}posterior mean of tau^2{p_end}
{synopt:{cmd:e(lambda_scale)}}lambda_scale used{p_end}
{synopt:{cmd:e(tau_scale)}}tau_scale used{p_end}
{synopt:{cmd:e(burnin)}}number of burnin iterations{p_end}
{synopt:{cmd:e(thin)}}thinning interval{p_end}
{synopt:{cmd:e(n_chol_fallback)}}number of iterations requiring Cholesky ridge adjustment{p_end}
{synopt:{cmd:e(ess_min)}}minimum effective sample size across all parameters{p_end}
{synopt:{cmd:e(ess_sigma2)}}effective sample size for sigma^2{p_end}
{synopt:{cmd:e(ess_tau2)}}effective sample size for tau^2{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:horseshoe}{p_end}
{synopt:{cmd:e(cmdline)}}full command as typed{p_end}
{synopt:{cmd:e(predict)}}{cmd:horseshoe_p}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(seed)}}seed specified in {opt seed()} (empty if not specified){p_end}
{synopt:{cmd:e(rng_state)}}full RNG state captured before sampling{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1 x p posterior means of beta{p_end}
{synopt:{cmd:e(V)}}p x p posterior variance-covariance matrix of beta{p_end}
{synopt:{cmd:e(b_sd)}}1 x p posterior standard deviations of beta{p_end}
{synopt:{cmd:e(b_lower)}}1 x p lower credible interval bounds{p_end}
{synopt:{cmd:e(b_upper)}}1 x p upper credible interval bounds{p_end}
{synopt:{cmd:e(penalized)}}1 x p penalization flags (0/1){p_end}
{synopt:{cmd:e(ess_beta)}}1 x p effective sample sizes for each beta_j{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
{bf:Example 1: Standard horseshoe (all coefficients penalized)}

{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. horseshoe price mpg weight length turn, nmcmc(2000) burnin(500) verbose}{p_end}

{pstd}
{bf:Example 2: Selective penalization}

{pstd}
Suppose you have treatment dummies that should NOT be shrunk, and interaction
terms that should be penalized:

{phang2}{cmd:. * Create penalization flag: 0 for treatment (cols 1-4), 1 for interactions (cols 5-20)}{p_end}
{phang2}{cmd:. matrix pen = J(1, 20, 1)}{p_end}
{phang2}{cmd:. matrix pen[1, 1] = 0}{p_end}
{phang2}{cmd:. matrix pen[1, 2] = 0}{p_end}
{phang2}{cmd:. matrix pen[1, 3] = 0}{p_end}
{phang2}{cmd:. matrix pen[1, 4] = 0}{p_end}
{phang2}{cmd:. horseshoe y treat1 treat2 treat3 treat4 x1-x16, penalized(pen)}{p_end}

{pstd}
{bf:Example 3: Reproducible results with seed()}

{phang2}{cmd:. horseshoe price mpg weight length, seed(12345)}{p_end}
{phang2}{cmd:. matrix list e(b)}{p_end}
{phang2}{cmd:. * Run again with same seed — results will be identical}{p_end}
{phang2}{cmd:. horseshoe price mpg weight length, seed(12345)}{p_end}

{pstd}
{bf:Example 4: Saving posterior draws for custom analysis}

{phang2}{cmd:. horseshoe price mpg weight length, saving(posterior_draws.dta) replace seed(42)}{p_end}
{phang2}{cmd:. * Load draws and compute custom summaries}{p_end}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. use posterior_draws.dta, clear}{p_end}
{phang2}{cmd:. summarize sigma2 tau2}{p_end}
{phang2}{cmd:. * Posterior probability that b_mpg > 0}{p_end}
{phang2}{cmd:. count if b_mpg > 0}{p_end}
{phang2}{cmd:. display "Pr(beta_mpg > 0) = " r(N) / _N}{p_end}
{phang2}{cmd:. restore}{p_end}

{pstd}
{bf:Example 5: Aggressive shrinkage with modified scales}

{phang2}{cmd:. horseshoe y x1-x100, lambdascale(1) tauscale(0.05) nmcmc(2000)}{p_end}

{pstd}
{bf:Example 6: Checking convergence diagnostics}

{phang2}{cmd:. horseshoe price mpg weight length}{p_end}
{phang2}{cmd:. display "Min ESS = " e(ess_min)}{p_end}
{phang2}{cmd:. display "ESS sigma^2 = " e(ess_sigma2)}{p_end}
{phang2}{cmd:. display "ESS tau^2 = " e(ess_tau2)}{p_end}
{phang2}{cmd:. matrix list e(ess_beta)}{p_end}

{pstd}
{bf:Example 7: Predict and postestimation}

{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. horseshoe price mpg weight length, nmcmc(1000) seed(42)}{p_end}
{phang2}{cmd:. * Fitted values}{p_end}
{phang2}{cmd:. predict yhat, xb}{p_end}
{phang2}{cmd:. * Residuals}{p_end}
{phang2}{cmd:. predict resid, residuals}{p_end}
{phang2}{cmd:. * Standard errors of prediction}{p_end}
{phang2}{cmd:. predict se_yhat, stdp}{p_end}
{phang2}{cmd:. * Test a linear hypothesis}{p_end}
{phang2}{cmd:. test mpg = weight}{p_end}
{phang2}{cmd:. * Linear combination}{p_end}
{phang2}{cmd:. lincom mpg + weight}{p_end}
{phang2}{cmd:. * Replay results}{p_end}
{phang2}{cmd:. horseshoe}{p_end}


{marker mata}{...}
{title:CATE extraction (horseshoe_cates.ado)}

{pstd}
A companion command {cmd:horseshoe_cates} extracts individual-level conditional
average treatment effects (CATEs) from posterior draws saved by
{cmd:horseshoe, saving()}, for multi-arm treatment designs.

{pstd}
{bf:Stata command interface:}

{phang2}{cmd:. horseshoe_cates, draws(posterior.dta) xtest(X_test) xsd(X_sd) nx(54) nd(4)}{p_end}

{pstd}
Results are returned as {cmd:r()} matrices: {cmd:r(cate_hat)}, {cmd:r(cate_lo)},
{cmd:r(cate_hi)}.

{pstd}
{bf:Mata function interface} (for pipeline scripts):

{phang2}{cmd:. run "horseshoe_cates.ado"}{p_end}
{phang2}{cmd:. mata: _hs_extract_cates(beta_all, X_test, X_sd, n_x, n_d, lo_prob, hi_prob)}{p_end}

{pstd}
See the header of {cmd:horseshoe_cates.ado} for full argument documentation.


{marker dependencies}{...}
{title:Dependencies}

{pstd}
The horseshoe package consists of the following files, all of which must be
on the adopath:

{p 8 8 2}
{cmd:horseshoe.ado} — main estimation command{break}
{cmd:horseshoe_p.ado} — predict postestimation command{break}
{cmd:horseshoe_cates.ado} — CATE extraction command and Mata functions{break}
{cmd:_hs_utils.ado} — shared Mata utilities (required by horseshoe.ado and horseshoe_cates.ado){break}
{cmd:horseshoe.sthlp} — this help file


{marker references}{...}
{title:References}

{phang}
Carvalho, C. M., N. G. Polson, and J. G. Scott. 2010.
The horseshoe estimator for sparse signals.
{it:Biometrika} 97(2): 465-480.

{phang}
Geyer, C. J. 1992.
Practical Markov chain Monte Carlo.
{it:Statistical Science} 7(4): 473-483.

{phang}
Makalic, E., and D. F. Schmidt. 2015.
A simple sampler for the horseshoe estimator.
{it:arXiv:1508.03884v4}.

{phang}
Bhattacharya, A., D. Chakraborty, and B. K. Mallick. 2016.
Fast sampling with Gaussian scale-mixture priors in high-dimensional regression.
{it:Biometrika} 103(4): 985-991. (Related work; this implementation uses the
standard Makalic and Schmidt sampler, not the Bhattacharya et al. fast algorithm.)


{marker author}{...}
{title:Author}

{pstd}
Andrew McEachin{p_end}
