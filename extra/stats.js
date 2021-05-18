//==============================================================================
// Statistical fcts
//                      
// Author: Patrick Bélisle


// Version 0.4 (May 2021)
// [distributed]

  // requires genericFcts.js
  // requires log-arithmetic.js
  // requires dataPreparation/*.js 
  //  (in particular S.js, but others as well, at least in the current version)


// Change log
// ===============================
//
//
// Version 0.4 (May 2021)
// ----------------------
//   Added the following function:
//     - CrI
//     - desc_stats
//     - rbinom
//  
//   The following function(s) was (were) modified [now calling as_log]
//     - phi_diff [it also now accepts mu argument as an array]
//
//   The definition of the following functions were moved from Array.prototype.* to classical fct defns
//     - phi
//     - Phi
//
//   The Jacobian in the Newton-Raphson algorithm is now defined through MyMatrix
//
//   rbeta is now based on variable type DensityFct
//
//   The code in rnorm_gt, rnorm_lt & rnorm_interval was revisited.
//
//   Functions Phi & pnorm now accept arrays in mean & sd arguments.
//
//   The property .domain was added to all logf* objects.
//   The property .quantile (used in OneSubjectEstimates) was added to a few logf* objects.
//
//
// Version 0.3 (Apr 2021)
// ----------------------
//   Added the following functions:
//     - pbeta, qbeta, rbeta
//
//   The functions in logf.js & stat-BetaDistrn.js
//   were included in the present file.
//
//
// Version 0.2 (Mar 2021)
// ----------------------
//   - Slighlty modified calculation of phi
//      [that fct is not using zygotine any more either]
//   - Moved fct varphi from normal-OneSubjectEstimates.js to this file
//   - Phi_interval (not used yet) was renamed PhiInterval
//
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
      
  var J = create_matrix(2, 2);

  
  while (cont)
  {    
    var Beta = BetaDensity(alpha, beta);
    var Beta1 = BetaDensity(alpha, beta + epsilon.derivative);
    var Beta2 = BetaDensity(alpha + epsilon.derivative, beta);      
    
    
    F1 = Beta.cdf(lcl);             
    F2 = Beta.cdf(ucl);
    
    F1pa = (Beta2.cdf(lcl) - F1) / epsilon.derivative;
    F1pb = (Beta1.cdf(lcl) - F1) / epsilon.derivative;
    
    F2pa = (Beta2.cdf(ucl) - F2) / epsilon.derivative;
    F2pb = (Beta1.cdf(ucl) - F2) / epsilon.derivative;
    
    J.m = [[F1pa, F1pb], [F2pa, F2pb]]; // Jacobian
        
    change = NewtonRaphsonChange(J, [F1, F2], target);
    alpha -= change[0];
    beta  -= change[1];
    
    converged = max_abs(change) < epsilon.NR;
    cont = !converged && ++iter < 300;
    
    if (alpha < 0) alpha = (alpha + change[0]) / 2;
    if (beta  < 0)  beta = ( beta + change[1]) / 2;
  }
  
  
  if (!converged) alpha = NaN, beta = NaN;
  
  return {converged: converged, alpha: alpha, beta: beta};
} // end of BetaParmsFromQuantiles


CrI = function(arr, level=0.95)
{
  var CrI = [];
  var len = arr.length;
  var x = arr.slice(); // copy by value, not by reference
  var s = x.sort(function(a, b){return a-b}); // sorted values
  
  // CrI lower limit
  var j = (1-level) / 2 * len;
  var i = Math.floor(j);
  var z = s[i-1];
  if (j > i) z += (j-i) * (s[i] - z);
  CrI.push(z);
  
  // CrI upper limit
  var j = (1+level) / 2 * len;
  var i = Math.floor(j);
  var z = s[i-1];
  if (j > i) z += (j-i) * (s[i] - z);
  CrI.push(z);
  
  return CrI;
} // end of CrI


desc_stats = function(obj, CrI_level=0.95, ndigits=4)
{
  console.log("*** Desc stats (CrI level = %s %) ***", 100*CrI_level);
  
  for (const property in obj) 
  {
    let tmp = {};
    
    tmp.mean   = round(mean(obj[property]), ndigits);
    tmp.median = round(median(obj[property]), ndigits);
    tmp.CrI    = round(CrI(obj[property], CrI_level), ndigits);
  
    console.log(`${property}:`, tmp);
  }
} // end of desc_stats


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


phi = function(arr, log_output=false)
{
  const k = Math.log(2*Math.PI) / 2;
  var u = arr.map(z => -z*z/2 - k);
  
  if (log_output) return u;
  else            return u.map(Math.exp);
} // end of phi


phi_diff = function(interval, mean=0, sd=1, log_output=false)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mean: either a) an array of same length as interval.gt
  //           or b) as number
  // sd: a number
  
  // Returns phi((this-mean)/sd) - phi((a-mean)/sd)  if log_output = false
  // Returns the above as a log-notation object      if log_output = true
  
  if (typeof interval.gt == 'number')
  {
    let tmp = phi_diff({gt: [interval.gt], lt: [interval.lt]}, mean, sd, log_output);
    
    if (log_output)
    {
      tmp.x = tmp.x[0];
      tmp.s = tmp.s[0];
      return tmp;
    }
    else            return tmp[0];
  }
  else
  {
    const k = Math.log(2*Math.PI) / 2;
    var a, b;
    
    if (Array.isArray(mean))
    {
      a = as_log(interval.gt.map((x,i) => -0.5*((x-mean[i])/sd)**2 - k), true, -1);
      b = as_log(interval.lt.map((x,i) => -0.5*((x-mean[i])/sd)**2 - k), true);
    }
    else
    {
      a = as_log(interval.gt.map(x => -0.5*((x-mean)/sd)**2 - k), true, -1);
      b = as_log(interval.lt.map(x => -0.5*((x-mean)/sd)**2 - k), true);
    }

    var d = b.plus(a);
  
    if (log_output) return d;
    else            return d.as_real(false);
  }
} // end of phi_diff


Phi = function(arr, mean=0, sd=1, lower_tail=true, log_p=false)
{
  // arr: an array
  // mean: either a) an array of same length as arr
  //           or b) a number
  // sd: either a) an array of same length as arr
  //         or b) a number  
  
  if (typeof arr == 'number')
  {
    return zygotine.S.normal.cdf(x, mean, sd, lower_tail, log_p);
  }
  else if (Array.isArray(mean))
  {
    if (Array.isArray(sd))     return arr.map((z,i) => zygotine.S.normal.cdf(z, mean[i], sd[i], lower_tail, log_p));
    else                       return arr.map((z,i) => zygotine.S.normal.cdf(z, mean[i], sd, lower_tail, log_p));
  }
  else if (Array.isArray(sd))  return arr.map((z,i) => zygotine.S.normal.cdf(z, mean, sd[i], lower_tail, log_p));
  else                         return arr.map(z => zygotine.S.normal.cdf(z, mean, sd, lower_tail, log_p));
} // end of Phi


PhiInterval = function(interval, mean=0, sd=1, log_output=false)
{
  // interval: a literal with dimension gt & lt (two arrays of same length)
  // mean: an array (of same length as interval.gt & interval.lt)
  //    or a scalar
  // sd: scalar
  
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
    
  
  logPhi_a = as_log(logPhi_a, true, -1);
  logPhi_b = as_log(logPhi_b, true);
  
  var PhiInterval = logPhi_b.plus(logPhi_a).x;
  
  if (interval_dim_are_numbers)
  {
    // return numbers, as the input interval dimensions consisted in numbers
    if (!log_output) return Math.exp(PhiInterval[0]);
    else             return PhiInterval[0];
  }
  else if (!log_output) return PhiInterval.map(Math.exp);
  else                  return PhiInterval;
} // end of PhiInterval


pnorm = function(arr, mean=0, sd=1, lower_tail=true, log_p=false)
{
  // arr: an array
  // mean: either a) an array of same length as arr
  //           or b) a number
  // sd: either a) an array of same length as arr
  //         or b) a number  
  
  if (typeof arr == 'number')
  {
    return zygotine.S.normal.cdf(x, mean, sd, lower_tail, log_p);
  }
  else if (Array.isArray(mean))
  {
    if (Array.isArray(sd))     return arr.map((z,i) => zygotine.S.normal.cdf(z, mean[i], sd[i], lower_tail, log_p));
    else                       return arr.map((z,i) => zygotine.S.normal.cdf(z, mean[i], sd, lower_tail, log_p));
  }
  else if (Array.isArray(sd))  return arr.map((z,i) => zygotine.S.normal.cdf(z, mean, sd[i], lower_tail, log_p));
  else                         return arr.map(z => zygotine.S.normal.cdf(z, mean, sd, lower_tail, log_p));
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
    if (typeof mean == 'number') mean = as_array(mean, mean, size);
    if (typeof   sd == 'number')   sd = as_array(sd, sd, size);
    
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
  var Beta = Object.create(DensityFct);
    var f = new logf_beta(alpha, beta);
    Beta.alpha  = f.alpha;
    Beta.beta   = f.beta;
    Beta.mode   = f.mode;
    
    Beta.g_compute = f.g_compute; 
    Beta.g         = f.g;
    
    Beta.logf        = f.logf;
    Beta.logf_prime  = f.logf_prime;
    Beta.logf_second = f.logf_second;
    
    Beta.domain = f.domain;
    Beta.label  = f.label;
    
  
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


rbinom = function(N, p, len=1)
{
  // N: a number (the maximum value taken by the binomial random variable for which we are generating random values)
  // p: a number or an array, with probability-ies of 'success'
  // len: a number (IMPORTANT: this parameter value is ignored when p is an array)


  if (Array.isArray(p))
  {
    let out = [];
    
    for (let i=0; i<p.length; i++)
    {
      let out_i = 0;
      let U = runif(N);
      for (let j=0; j<N; j++) out_i += U[j] <= p[i] ? 1 : 0;
      out.push(out_i);      
    }
    
    return out;    
  }
  else
  {
    // p is a number
    
    if (len == 1)
    {
      let out = 0;
      let U = runif(N);
      for (let i=0; i<N; i++) out += U[i] <= p ? 1 : 0;
      
      return out;
    }
    else
    {
      // len > 1
      let out = [];
      
      for (let i=0; i<len; i++)
      {
        let out_i = 0;
        let U = runif(N);
        for (let j=0; j<N; j++) out_i += U[j] <= p ? 1 : 0;
        out.push(out_i);      
      }
      
      return out;
    }
  }
} // end of rbinom


// ICI cette fct ne fait pas ce que j'attendais
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


rnorm_gt = function(mean, sd, gt)
{
  // mean & sd: scalars
  // gt: array or number
  
  var gt_is_array = Array.isArray(gt);
  if (!gt_is_array) gt = [gt];
  
  var z = [];
  
  for (let i=0; i<gt.length; i++)
  {
    let a = zygotine.S.normal.cdf(gt[i], mean, sd, false, true); // pnorm[log] - upper tail
    let u = runif(1);
    let logp = a + Math.log(u);
    let x = qnorm(logp, mean, sd, false, true);
    z.push(x);
  }
  
  if (gt_is_array) return z;
  else             return z[0];
} // end of rnorm_gt



rnorm_interval = function(mean, sd, gt, lt)
{
  // mean & sd: scalars
  // gt & lt: either a) two arrays of same length 
  //              or b) two numbers
  
  var lt_is_array = Array.isArray(lt);
  if (!lt_is_array) lt = [lt], gt = [gt];
  
  var z = [];
  
  for (let i=0; i<lt.length; i++)
  {
    let lower_tail = gt[i] < mean;
    
    if (lower_tail)
    {
      a = zygotine.S.normal.cdf(gt[i], mean, sd, true, true); // pnorm[log]
      b = zygotine.S.normal.cdf(lt[i], mean, sd, true, true);
    }
    else
    {
      a = zygotine.S.normal.cdf(lt[i], mean, sd, false, true);
      b = zygotine.S.normal.cdf(gt[i], mean, sd, false, true);
    }
  
  
    let u = runif(1);
    let logp = b + Math.log(1 - u*(1 - Math.exp(a-b)));
    let x = qnorm(logp, mean, sd, lower_tail, true);
    z.push(x);
  }
  
  if (lt_is_array) return z;
  else             return z[0];
} // end of rnorm_interval


rnorm_lt = function(mean, sd, lt)
{
  // mean & sd: scalars
  // lt: array or number
  
  var lt_is_array = Array.isArray(lt);
  if (!lt_is_array) lt = [lt];
  
  var z = [];
  
  for (let i=0; i<lt.length; i++)
  {
    let a = zygotine.S.normal.cdf(lt[i], mean, sd, true, true); // pnorm[log]
    let u = runif(1);
    let logp = a + Math.log(u);
    let x = qnorm(logp, mean, sd, true, true);
    z.push(x);
  }
  
  if (lt_is_array) return z;
  else             return z[0];
} // end of rnorm_lt


runif = function(size, min=0, max=1)
{
  // size: integer
  // min: Number or array of length = size
  // max: Number or array of length = size
  
  if (size == 1)
  {    
    return zS.uniform.sample(1, min, max);
  }
  else if (typeof min == 'number' && typeof max == 'number')
  {
    return zS.uniform.sample(size, min, max);
  }
  else
  {
    // Fct is called to generate more than one random value
    // and at least one of min & max is an array
    
    var u = [],
        m, M; // min & max
        
    if (!Array.isArray(min)) m = min;
    if (!Array.isArray(max)) M = max;
    
    for (let i=0; i<size; i++)
    {
      if (Array.isArray(min)) m = min[i];
      if (Array.isArray(max)) M = max[i];
       
      u.push(zS.uniform.sample(1, m, M));
    }
    
    return u;
  }
} // end of runif


sample_discrete = function(values_list, p, log_p=false)
{
  // values_list & p: two arrays of same length
  // Sample a value from 'values_list' with corresponding weights found in p
  
  var sample_discrete;
  
  if (log_p)
  {
    let p_max = max(p);
    p = p.map(p => p - p_max).map(Math.exp);    
  }

  var p_sum = sum(p);
  p = p.map(p => p/p_sum);  // standardize (make them sum to 1)
  
  var U = runif(1);
  
  for (let i=0; i<values_list.length; i++)
  {
    if (U <= p[i])
    {
      sample_discrete = values_list[i];
      break;
    }
    
    U -= p[i];
  }
  
  if (typeof sample_discrete == 'undefined') sample_discrete = values_list[values_list.length-1]; // to prevent numeric imprecision
  
  return sample_discrete;
} // end of sample_discrete


varphi = function(arr)
{  
  var log_phi = phi(arr, true);
  var log_Phi = Phi(arr, 0, 1, true, true);
  
  var vphi       = substract(log_phi, log_Phi).map(Math.exp);
  var vphi_prime = vphi.map((v,i) => -(arr[i]*v + v**2));
  
  return {f: vphi, fp: vphi_prime};
} // end of varphi


// Deprecated
//logPhiInterval = function(a, b, mean, sd)
//{
//  // a, b, mean, sd: all Numbers
//  // Return log(pnorm(b, mean, sd) - pnorm(a, mean, sd)), in R notation
//    
//  var logp1 = as_log(Phi([a], mean, sd, true, true),    true, -1);
//  var logp2 = as_log(Phi([b], mean, sd, true, true),    true);
//  
//  log_pdiff = logp2.plus(logp1);
//  
//  return log_pdiff.x[0];
//} // end of logPhiInterval


////////////////////////////////////////////////////////////////////////////////
// The functions below are useful, for example, for use as prior distributions 
// for parameters to include in posterior or full conditional 
// posterior distributions
//

  // The function g (see logf_beta.prototype below) can be given (it is optional)
  // in case exp(logf) would lead to NaN when calculing cdf: indeed, if we take the case
  // of the Beta distribution, logf(x) = NaN for x = 0 when alpha = 1 
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
  
    //this.mode = this.mode();
    this.g_compute = this.use_g();
    //this.compute_M(mode);
    //this.total_area = this.area();
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
  this.mean   = mean;
  this.sd     = sd;
  this.sd2    = sd**2;
  this.prec   = 1/sd**2;
  this.label  = 'logNormal distrn';
  this.domain = [0, Infinity];
} // end of logf_lnorm


logf_lnorm.prototype = 
{
  logf:        function(x){return -Math.log(x) - this.prec/2 * (Math.log(x) - this.mean)**2},
  logf_prime:  function(x){return - (this.prec * (Math.log(x) - this.mean) + 1) / x},
  logf_second: function(x){return - (this.prec * (this.mean + 1 - Math.log(x)) - 1) / x**2},
  mode:        function(){return Math.exp(this.mean-this.sd2)},
  quantile:    function(p){return Math.exp(qnorm(p, this.mean, this.sd))}
} // end of logf_lnorm


logf_norm = function(mean, sd) 
{
  this.mean   = mean;
  this.sd     = sd;
  this.prec   = 1/sd**2;
  this.label  = 'Normal distrn';
  this.domain = [-Infinity, Infinity];
} // end of logf_norm


logf_norm.prototype =
{
  logf:        function(x){return -this.prec/2 * (x-this.mean)**2},
  logf_prime:  function(x){return -this.prec * (x - this.mean)},
  logf_second: function(x){return -this.prec},
  mode:        function(){return this.mean},
  quantile:    function(p){return qnorm(p, this.mean, this.sd)}
} // end of logf_norm.prototype


logf_unif = function(lower, upper)
{
  if (typeof lower == 'undefined') lower = -Infinity;
  if (typeof upper == 'undefined') upper =  Infinity;
  
  this.domain = [lower, upper];
}; // end of logf_unif


logf_unif.prototype =
{
  logf: function(x){return 0},
  logf_prime: function(x){return 0},
  logf_second: function(x){return 0}
} // end of logf_unif.prototype


////////////////////////////////////////////////////////////////////////////////
// Defining Density Functions


function BetaDensity(alpha, beta)
{
  var Beta = Object.create(DensityFct);
    var f = new logf_beta(alpha, beta);
    Beta.alpha  = f.alpha;
    Beta.beta   = f.beta;
    Beta.mode   = f.mode;
    
    Beta.g_compute = f.g_compute; 
    Beta.g         = f.g;
    
    Beta.logf        = f.logf;
    Beta.logf_prime  = f.logf_prime;
    Beta.logf_second = f.logf_second;
    
    Beta.domain = f.domain;
    Beta.label  = f.label;
    
  return Beta;
} // 