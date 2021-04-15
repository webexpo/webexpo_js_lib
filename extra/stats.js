//==============================================================================
// Statistical fcts
//
// Version 0.3 (Apr 2021)
//                      
// Author: Patrick Bélisle


// requires genericFcts.js
// requires log-arithmetic.js
// requires dataPreparation/*.js 
//  (in particular S.js, but others as well, at least in the current vesion)


// Change log
// ======================
//
// Version 0.3 (Apr 2021)
// ----------------------
//   Added the following functions:
//     - pbeta, qbeta, rbeta
//
//   The functions in logf.js & stat-BetaDistrn.js
//   were included in the present file.
//
// Version 0.2 (Mar 2021)
// ----------------------
//   - Slighlty modified calculation of phi
//      [that fct is not using zygotine any more either]
//   - Moved fct varphi from normal-OneSubjectEstimates.js to this file
//   - Phi_interval (not used yet) was renamed PhiInterval
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code


BetaParmsFromMoments = function(m, v)
{
  alpha = m*m*(1-m)/v-m;
  beta = alpha*(1/m-1);
  
  return {alpha: alpha, beta: beta};
} // end of BetaParmsFromMoments


BetaParmsFromQuantiles = function(lcl, ucl, level)
{
  if (typeof level == 'undefined') level = 0.95;

  const epsilon = {derivative: 1e-6, NR: 1e-6};
  const Z = qnorm((1 + level) / 2);
  
  var  m = (lcl + ucl) / 2,
      sd = (ucl - lcl) / (2*Z);
      
  var tmp = BetaParmsFromMoments(m, sd**2);
  
  var alpha = tmp.alpha,
      beta  = tmp.beta,
      area, area1, area2, m,
      target = [(1-level)/2, (1+level)/2],
      cont = true,
      converged,
      mode, mode1, mode2, 
      iter = 0;
      

    
  
  while (cont)
  {
    var Beta  = new logf_beta(alpha,                      beta);
    var Beta1 = new logf_beta(alpha,                      beta + epsilon.derivative);
    var Beta2 = new logf_beta(alpha + epsilon.derivative, beta);  
    
    F1 = Beta.cdf(lcl);             
    F2 = Beta.cdf(ucl);
    
    F1pa = (Beta2.cdf(lcl) - F1) / epsilon.derivative;
    F1pb = (Beta1.cdf(lcl) - F1) / epsilon.derivative;
    
    F2pa = (Beta2.cdf(ucl) - F2) / epsilon.derivative;
    F2pb = (Beta1.cdf(ucl) - F2) / epsilon.derivative;
    
    J1 = [F1pa, F1pb];
    J2 = [F2pa, F2pb];
    J = [J1, J2];
        
    change = NewtonRaphsonChange(J, [F1, F2], target);
    alpha -= change[0];
    beta  -= change[1];
    
    converged = change.map(Math.abs).max() < epsilon.NR;
    cont = !converged && ++iter < 300;
    
    if (alpha < 0) alpha = (alpha + change[0]) / 2;
    if (beta  < 0)  beta = ( beta + change[1]) / 2;
  }
  
  
  if (!converged) alpha = NaN, beta = NaN;
  
  return {converged: converged, alpha: alpha, beta: beta};
} // end of BetaParmsFromQuantiles


Array.prototype.CI = function(level=0.95)
{
  var CI = [];
  var len = this.length;
  var s = this.sort(function(a, b){return a-b}); // sorted values
  
  // CI lower limit
  var j = (1-level) / 2 * len;
  var i = Math.floor(j);
  var z = s[i-1];
  if (j > i) z += (j-i) * (s[i] - z);
  CI.push(z);
  
  // CI upper limit
  var j = (1+level) / 2 * len;
  var i = Math.floor(j);
  var z = s[i-1];
  if (j > i) z += (j-i) * (s[i] - z);
  CI.push(z);
  
  return CI;
} // end of Array.CI


lgamma = function(z)
{
  // z: number or array
  if (typeof z == 'number') return zNum.logGamma(z);
  else return z.map(z => zNum.logGamma);
} // end of lgamma


logPhi = function(z)
{
  return zygotine.S.normal.cdf(z, 0, 1, true, true);
} // end of logPhi


pbeta = function(x, alpha, beta)
{
  var Beta = new logf_beta(alpha, beta);
  return Beta.cdf(x);
} // end of pbeta


Array.prototype.phi = function(log_output=false)
{
  const k = Math.log(2*Math.PI) / 2;
  var u = this.map(z => -z*z/2 - k);
  
  if (log_output) return u;
  else            return u.map(Math.exp);
} // end of Array.phi


phi_diff = function(interval, mean=0, sd=1, log_output=false)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mean & sd: two numbers
  
  // Returns phi((this-mean)/sd) - phi((a-mean)/sd)  if log_output = false
  // Returns the above as a log-notation object      if log_output = true
  
  if (typeof interval.gt == 'number')
  {
    let tmp = phi_diff({gt: [interval.gt], lt: [interval.lt]}, mean, sd, log_output);
    
    if (log_output) return {gt: tmp.gt[0], lt: tmp.lt[0]};
    else            return tmp[0];
  }
  else
  {
    const k = Math.log(2*Math.PI) / 2;
  
    a = interval.gt.map(x => -0.5*((x-mean)/sd)**2 - k);
    b = interval.lt.map(x => -0.5*((x-mean)/sd)**2 - k);

    a = {x: a, s: [-1].rep(a.length)};
    b = {x: b, s:  [1].rep(b.length)};
  
    var d = b.plus_log(a);
  
    if (log_output) return d;
    else            return d.as_real();
  }
} // end of phi_diff



Array.prototype.Phi = function(mean=0, sd=1, lower_tail=true, log_p=false)
{
  return this.map(x => zygotine.S.normal.cdf(x, mean, sd, lower_tail, log_p));
} // end of Array.Phi


PhiInterval = function(interval, mean=0, sd=1, log_output=false)
{
  // interval: a literal with dimension gt & lt (two arrays of same length)
  // mean, sd: scalar
  
  // Returns pnorm(b, mean, sd) - pnorm(this, mean, sd)  [in R notation],
  //   or its log() value if log_output = true
  
  
  var a, b;
  var interval_dim_are_numbers = typeof interval.gt == 'number';
  
  if (interval_dim_are_numbers) 
  {
    a = [interval.gt];
    b = [interval.lt];
    // now a & b are arrays
  }
  else
  {
    a = interval.gt.slice();
    b = interval.lt.slice();    
  }
    
  
  var logPhi_a,
      logPhi_b;
  
  if (Array.isArray(mean))
  {
    logPhi_a = a.map((x, i) => zygotine.S.normal.cdf(x, mean[i], sd, true, true));  
    logPhi_b = b.map((x, i) => zygotine.S.normal.cdf(x, mean[i], sd, true, true));
  }
  else
  {
    logPhi_a = a.map(x => zygotine.S.normal.cdf(x, mean, sd, true, true));
    logPhi_b = b.map(x => zygotine.S.normal.cdf(x, mean, sd, true, true));
  }
    
  
  logPhi_a = {x: logPhi_a, s: [-1].rep(b.length)}; 
  logPhi_b = {x: logPhi_b, s:  [1].rep(b.length)};
  
  var PhiInterval = logPhi_b.plus_log(logPhi_a).x;
  
  if (interval_dim_are_numbers)
  {
    // return numbers, as the input interval dimensions consisted in numbers
    if (!log_output) return Math.exp(PhiInterval[0]);
    else             return PhiInterval[0];
  }
  else if (!log_output) return PhiInterval.map(Math.exp);
  else                  return PhiInterval;
} // end of PhiInterval


pnorm = function(x, mean=0, sd=1, lower_tail=true, log_p=false)
{
  // x: number or array
  
  if (typeof x == 'number')
  {
    return zygotine.S.normal.cdf(x, mean, sd, lower_tail, log_p);
  }
  else
  {
    return x.map(z => zygotine.S.normal.cdf(z, mean, sd, lower_tail, log_p));
  }
} // end of pnorm


qbeta = function(p, alpha, beta)
{
  var Beta = new logf_beta(alpha, beta);
  return Beta.quantile(p);
} // end of pbeta


function qnorm(p, mean=0, sd=1, lower_tail=true, log_p=false)
{
  // p: either a number or an array
  // mean, sd: if p is a number, (mean, sd) must be numbers
  //            if p is an array, (mean, sd) can be numbers or arrays of same length as p
  // lower_tail, log_p: boolean ()
  
  if (typeof p == 'number')
  {
    return zygotine.S.normal.icdf(p, mean, sd, lower_tail, log_p);
  }
  else
  {
    // p is an array
    let z = [];
    let size = p.length;
    
    // Turn mean & sd into arrays if they're not already
    if (typeof mean == 'number')    mean = as_array(mean, mean, size);
    if (typeof   sd == 'number') sd = as_array(sd, sd, size);
    
    for (let i=0; i<size; i++)
    {
      let x = zygotine.S.normal.icdf(p[i], mean[i], sd[i], lower_tail, log_p);
      z.push(x);
    }
    
    return z;
  }
} // end of qnorm


rbeta = function(alpha, beta, n, other, iter)
{
  var Beta = new logf_beta(alpha, beta);
//    Beta.other = other;
//    Beta.iter = iter;
  
  if (typeof n == 'undefined')
  {
    if (typeof alpha == 'number') n = 1;
    else n = alpha.length;
  }
  
  var U = runif(n);
  
  tmp = Beta.quantile(U);
  
  if (tmp.error)
  { 
    console.error("Error in rbeta");
    return;
  }
  else return tmp.x;
} // end of rbeta


// ICI cette fct ne fait pas ce que j'attends
//ResetSeed = function()
//{
//  var seed = new Date().getTime();
//  
//  zygotine.S.prng.init(12); // Reset Seed [= function of date and time]
//} // end of ResetSeed


rnorm = function(size, mean=0, sd=1)
{
  // size: an integer
  // (mean, sd): either numbers or arrays of length = size
  
  var u = runif(size);
  
  return qnorm(u, mean, sd);
} // end of rnorm


Array.prototype.rnorm_gt = function(mean, sd)
{
  // mean & sd: scalar
  
  var z = [];
  
  for (let i=0; i<this.length; i++)
  {
    let x = zO.rNormLeftCensored(mean, sd, this[i]);
    z.push(x);
  }
  
  return z;
} // end of rnorm_gt


rnorm_interval = function(mean, sd, gt, lt)
{
  var z = [];
  
  for (let i=0; i<gt.length; i++)
  {
    let x = zO.rNormIntervalCensored(mean, sd, gt[i], lt[i]);
    z.push(x);
  }
  
  return z;
} // end of rnorm_interval


Array.prototype.rnorm_lt = function(mean, sd)
{
  // mean & sd: scalar
  
  var z = [];
  
  for (let i=0; i<this.length; i++)
  {
    let x = zO.rNormRightCensored(mean, sd, this[i]);
    z.push(x);
  }
  
  return z;
} // end of rnorm_lt


runif = function(size, min, max)
{
  // size: integer
  // min: Number or array of length = size
  // max: Number or array of length = size
  
  if (size == 1)
  {
    if (typeof min === 'undefined') min = 0.0;
    if (typeof max === 'undefined') max = 1.0;
    
    return zS.uniform.sample(1, min, max);
  }
  else
  {
    // Fct is called to generate more than one random value => min & max should be (or converted to) arrays
    
    let u = [];
  
    let a = as_array(min, 0, size); // defaults to rep(0, size) [in R notation]
    let b = as_array(max, 1, size); // defaults to rep(1, size) [in R notation]
    
    for (let i=0; i<size; i++) u.push(zS.uniform.sample(1, a[i], b[i]));
    return u;
  }
} // end of runif


Array.prototype.sample_discrete = function(p, log_p=false)
{
  // this & p: two arrays of same length
  // Sample a value from 'this' with corresponding weights found in p
  
  var sample_discrete;
  
  if (log_p)
  {
    let p_max = p.max();
    p = p.map(p => p - p_max).map(Math.exp);    
  }

  var p_sum = p.sum();
  p = p.map(p => p/p_sum);  // standardize (make them sum to 1)
  
  var U = runif(1);
  
  for (let i=0; i<this.length; i++)
  {
    if (U <= p[i])
    {
      sample_discrete = this[i];
      break;
    }
    
    U -= p[i];
  }
  
  if (typeof sample_discrete == 'undefined') sample_discrete = this[this.length-1]; // to prevent numeric imprecision
  
  return sample_discrete;
} // end of Array.sample_discrete


Array.prototype.varphi = function()
{
  // this: array
  
  var log_phi = this.phi(true);
  var log_Phi = this.Phi(0, 1, true, true);
  
  var vphi       = log_phi.minus(log_Phi).map(Math.exp);
  var vphi_prime = this.times(vphi).plus(vphi.sq()).times(-1);
  
  var varphi = {f: vphi, fp: vphi_prime};
  
  return varphi;
} // end of Array.varphi


////////////////////////////////////////////////////////////////////////////////
// The functions below are useful, for example, for use as prior distributions 
// for parameters to include in posterior or full conditional 
// posterior distributions
//

  // The function g (see logf_beta.prototype below) can be given (it is optional)
  // in case exp(logf) would lead to NaN when calculing cdf: indeed, if we take the case
  // the Beta distribution, logf(x) = NaN for x = 0 when alpha = 1 
  // as (alpha-1)*log(x) = 0 * log(0) = NaN, and so does exp(logf(0)). However, in the *real*
  // distribution function x**(alpha-1) = 0**0 is well defined (and = 1); in order to
  // go around that computational problem, one provides g [which could have been called "alternate f"]
  // which will be used by the algorithm when the cdf using f [the exponential of logf] leads to NaN.


logf_beta = function(alpha, beta)
{
  this.alpha = alpha;
  this.beta  =  beta;
  
  this.domain = [0, 1];
  this.label = 'Beta distrn';
  
    mode = this.mode();
    this.g_compute = this.use_g();
    this.compute_M(mode);
    this.total_area = this.area();
} // end of logf_beta


logf_beta.prototype = 
{
  g:           function(x){return x**(this.alpha-1) * (1-x)**(this.beta-1) / this.M},
  logf:        function(x){return (this.alpha - 1) * Math.log(x) + (this.beta - 1) * Math.log(1-x)},
  logf_prime:  function(x){return (this.alpha - 1)/x - (this.beta - 1)/(1-x)},
  logf_second: function(x){return -(this.alpha - 1)/x**2 - (this.beta - 1)/(1-x)**2},
  mode:        function() {return (this.alpha - 1) / (this.alpha + this.beta - 2)},
  use_g:       function() {return this.alpha == 1 || this.beta == 1}
} // end of logf_beta


logf_lnorm = function(mean, sd) 
{
  this.mean = mean;
  this.prec = 1/sd**2;
  this.label = 'logNormal distrn';
} // end of logf_lnorm


logf_lnorm.prototype = 
{
  logf:        function(x){return -Math.log(x) - this.prec/2 * (Math.log(x) - this.mean)**2},
  logf_prime:  function(x){return - (this.prec * (Math.log(x) - this.mean) + 1) / x},
  logf_second: function(x){return - (this.prec * (this.mean + 1 - Math.log(x)) - 1) / x**2}
} // end of logf_lnorm


logf_norm = function(mean, sd) 
{
  this.mean = mean;
  this.prec = 1/sd**2;
  this.label = 'Normal distrn';
} // end of logf_norm


logf_norm.prototype =
{
  logf:        function(x){return -this.prec/2 * (x-this.mean)**2},
  logf_prime:  function(x){return -this.prec * (x - this.mean)},
  logf_second: function(x){return -this.prec}
} // end of logf_norm.prototype


logf_unif = function(){}; // end of logf_unif


logf_unif.prototype =
{
  logf: function(x){return 0},
  logf_prime: function(x){return 0},
  logf_second: function(x){return 0}
} // end of logf_unif.prototype



// _______________________________________________________________________________
//
// ___ Used in Between-Worker ____________________________________________________
//
// Save the functions below as they are still used by the Between-Worker algorithm
// (as well as riskband)
// but new versions are available above for new uses


logPhiInterval = function(a, b, mean, sd)
{
  // a, b, mean, sd: all Numbers
  // Return log(pnorm(b, mean, sd) - pnorm(a, mean, sd)), in R notation
  
  var logp1 = {x: [], s: [1]},
      logp2 = {x: [], s: [1]};
      
  logp1.x = [a].Phi(mean, sd, true, true);
  logp2.x = [b].Phi(mean, sd, true, true);
  
  log_pdiff = logp2.minus_log(logp1);
  
  return log_pdiff.x[0];
} // end of logPhiInterval
