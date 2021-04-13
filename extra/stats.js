//==============================================================================
// Statistical fcts
//
// Version 0.1 (Mar 2021)
//
// Author: Patrick Bélisle


// requires genericFcts.js
// requires log-arithmetic.js
// requires dataPreparation/*.js (in particular S.js, but others as well, at least in the current vesion)


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


logPhiInterval = function(a, b, mu, sigma)
{
  // a, b, mu, sigma: all Numbers
  // Return log(pnorm(b, mu, sigma) - pnorm(a, mu, sigma)), in R notation
  
  var logp1 = {x: [], s: [1]},
      logp2 = {x: [], s: [1]};
      
  logp1.x = [a].map(z => (z-mu)/sigma).Phi(true, true);
  logp2.x = [b].map(z => (z-mu)/sigma).Phi(true, true);
  
  log_pdiff = log_diff(logp1, logp2);
  
  return log_pdiff.x[0];
} // end of logPhiInterval


Array.prototype.phi = function(log_p=false)
{
  return this.map(z => zygotine.S.normal.pdf(z, 0, 1, log_p));
} // end of Array.phi


Array.prototype.Phi = function(lower_tail=true, log_p=false)
{
  return this.map(z => zygotine.S.normal.cdf(z, 0, 1, lower_tail, log_p));
} // end of Array.Phi


pnorm = function(x, mu=0, sd=1, lower_tail=true, log_p=false)
{
  // x: number or array
  
  if (typeof x == 'number')
  {
    return zygotine.S.normal.cdf(x, mu, sd, lower_tail, log_p)
  }
  else
  {
    return x.map(z => zygotine.S.normal.cdf(z, mu, sd, lower_tail, log_p));
  }
} // end of pnorm


function qnorm(p, mu=0, sigma=1, lower_tail=true, log_p=false)
{
  // p: either a number or an array
  // mu, sigma: if p is a number, (mu, sigma) must be numbers
  //            if p is an array, (mu, sigma) can be numbers or arrays of same length as p
  // lower_tail, log_p: boolean ()
  
  if (typeof p == 'number')
  {
    return zygotine.S.normal.icdf(p, mu, sigma, lower_tail, log_p);
  }
  else
  {
    // p is an array
    let z = [];
    let size = p.length;
    
    // Turn mu & sigma into arrays if they're not already
    if (typeof    mu == 'number')    mu = as_array(mu, mu, size);
    if (typeof sigma == 'number') sigma = as_array(sigma, sigma, size);
    
    for (let i=0; i<size; i++)
    {
      let x = zygotine.S.normal.icdf(p[i], mu[i], sigma[i], lower_tail, log_p);
      z.push(x);
    }
    
    return z;
  }
} // end of qnorm


rnorm = function(size, mu=0, sigma=1)
{
  // size: an integer
  // (mu, sigma): either numbers or arrays of length = size
  
  var u = runif(size);
  return qnorm(u, mu, sigma);
} // end of rnorm


Array.prototype.rnorm_gt = function(mu, sigma)
{
  // mu & sigma: scalar
  
  var z = [];
  
  for (let i=0; i<this.length; i++)
  {
    let x = zO.rNormLeftCensored(mu, sigma, this[i]);
    z.push(x);
  }
  
  return z;
} // end of rnorm_gt


rnorm_interval = function(mu, sigma, gt, lt)
{
  var z = [];
  
  for (let i=0; i<gt.length; i++)
  {
    let x = zO.rNormIntervalCensored(mu, sigma, gt[i], lt[i]);
    z.push(x);
  }
  
  return z;
} // end of rnorm_interval


Array.prototype.rnorm_lt = function(mu, sigma)
{
  // mu & sigma: scalar
  
  var z = [];
  
  for (let i=0; i<this.length; i++)
  {
    let x = zO.rNormRightCensored(mu, sigma, this[i]);
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
    
    let u = [], a = [], b = [];
  
    a = as_array(min, 0, size); // defaults to rep(0, size) [in R notation]
    b = as_array(max, 1, size); // defaults to rep(1, size) [in R notation]
    
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
  
  u = runif(1);
  
  for (let i=0; i<this.length; i++)
  {
    if (u <= p[i])
    {
      sample_discrete = this[i];
      break;
    }
    
    u -= p[i];
  }
  
  if (typeof sample_discrete == 'undefined') sample_discrete = this[this.length-1]; // to prevent numeric imprecision
  
  return sample_discrete;
} // end of Array.sample_discrete



