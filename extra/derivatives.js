// Functions to calculate derivatives often involved in log(posterior)
//
// Author: Patrick Bélisle
//                                                                                                       
// Version 0.3 (May 2021)
// [distributed]                                              


// Change Log
// ========================
//
// Version 0.3 (May 2021)
// ----------------------
//   //   A correction was made in the function d2logPhiInterval_dsigma2
//
//   Commented out the two (unused) functions dprod_dtheta & dratio_dtheta
//
//   Changed the following fct calls:
//     - lo_div_  =>  ratio
//     - lo_mult  =>  mult
//     - lo_mult_ =>  mult
//
//   Now using mult, ratio & changed_sign in d2logPhi_dmudsigma
//
//   The calls to the following fcts were changed as they are not defined through Array.prototype anymore.
//     - dot_product
//     - sum
//
//   The definition of the following functions were moved from Array.prototype.* to classical fct defns:
//     - d2logPhi_dmudsigma
//
//   dlogPhi_dmu & d2logPhi_dmu2 were joined in one function -> dklogPhi_dmuk
//    (which is accepting mu as an array or a scalar)
//
//   dlogPhi_dsigma & d2logPhi_dsigma2 were also joined in one function -> dklogPhi_dsigmak
//    (which is accepting mu as an array or a scalar)
//
//   dlogPhiInterval_dsigma, d2logPhiInterval_dsigma2 are now accepting mu as an array or a scalar
//
//   dsumLogB_dtheta          was renamed dLogB_dtheta           with an additional argument return_sum=true   
//   d2sumLogB_dtheta1dtheta2 was renamed d2LogB_dtheta1dtheta2  with an additional argument return_sum=true
//   d2sumLogB_dtheta2        was renamed d2LogB_dtheta2         with an additional argument return_sum=true
//
//
// Version 0.2 (Apr 2021)
// ----------------------
//   Added the following functions:
//     - dlogPhi_dsigma
//     - d2logPhi_dsigma2
//     - dlogPhiInterval_dsigma
//     - d2logPhiInterval_dsigma2
//     - d2logPhi_dmudsigma
//     - d2logPhiInterval_dmudsigma
//
//   Changed calls to log_notation for calls to as_log


dLogB_dtheta = function(B, Bp, return_sum=true)
{
  // B & Bp are either
  //    a) both arrays (on natural scale) or
  // or b) B is an array (on log-scale) and Bp is a log-notation object

  // Returns d(sum(log(B_i)))/dtheta 
  //   given that the               B_i   terms are given in B, 
  //   and their first derivatives  B_i'  are given in Bp

  //  d(log(B))/dtheta = B'/B 

  
  if (Array.isArray(Bp))
  { 
    // B & Bp are arrays [on natural scale]
    if (return_sum) return ratio_sum(Bp, B);
    else            return ratio(Bp, B);
  }
  else
  {
    // Bp is a log-notation object, B an array (on log-scale)
    if (return_sum) return Bp.ratio_sum(B);
    else            return Bp.div(B).as_real(false);
  }
} // end of dLogB_dtheta


d2LogB_dtheta2 = function(B, Bp, Bs, return_sum=true)
{
  // B, Bp, Bs are either
  //    a) all arrays (on natural scale) or
  // or b) B is an array (on log-scale) and Bp & Bs are log-notation objects

  // Returns d2(sum(log(B_i)))/dtheta2 
  //   given that the            B_i terms are given in B, 
  //   their first derivatives   B_i'  given in Bp
  //   and their 2nd derivatives B_i'' given in Bs
  
  //  d2(log(B))/dtheta2 = d/dtheta (B'/B) = = (B''B - B'**2)/B**2 = B''/B - (B'/B)**2
  
    
  if (Array.isArray(Bp))
  { 
    // B, Bp & Bs are arrays
    if (return_sum) return ratio_sum(Bs, B) - ratio_sqSum(Bp, B);
    else
    {
      let a = ratio(Bs, B);
      let b = ratio(Bp, B);

      return a.map((a,i) => a - b[i]**2)
    }
  }
  else
  {
    // Bs is a log-notation object, B & Bp are arrays
    if (return_sum)  return Bs.ratio_sum(B) - Bp.ratio_sum(B, [2,2]); // scalar
    else
    {
      let a = Bs.div(B).as_real(false);                // arrays
      let b = Bp.div(B).as_real(false).map(z => z**2);
      
      return a.map((a,i) => a - b[i]);
    }
  }
} // end of d2LogB_dtheta2


d2LogB_dtheta1dtheta2 = function(B, Bp1, Bp2, Bm, return_sum=true)
{
  // B is an array
  // all other arguments are either:
  //     a) arrays (on natural scale, as well as B)
  //  or b) log-notation objects (and then B is an array on log-scale)
  
 
  if (Array.isArray(Bp1))
  {
    // scenario a
    if (return_sum) return ratio_sum(Bm, B) - sum(Bp1.map((b,i) => b * Bp2[i] / B[i]**2));
    else
    {
      let a = ratio(Bm, B);
      let b = Bp1.map((b,i) => b * Bp2[i] / B[i]**2);
      
      return a.map((a, i) => a - b[i]); // an array
    }
  }
  else
  {
    // scenario b
    if (return_sum)
    {
      let Bm1 = Bm.ratio_sum(B);                    // scalar
      let Bm2 = Bp1.mult(Bp2).ratio_sum(B, [1,2]);  // scalar
       
      return Bm1 - Bm2;
    }
    else
    {
      let Bm1 = Bm.div(B);               // log-notation objects
      let Bm2 = Bp1.mult(Bp2).div(B, 2); 
      
      return Bm1.minus(Bm2).as_real(false); // an array
    }
  } 
} // end of d2LogB_dtheta1dtheta2


dklogPhi_dmuk = function(arr, mu, sigma, lower_tail=true, return_sum=true)
{
  // arr: an array, with values y_i
  // mu:  an array of same length as arr
  //   or a scalar
  // sigma: scalar
  //
  // If return_sum = true, return the sum of the terms
  //   otherwise return an array with all the terms

  var mu_sign = lower_tail ? -1 : 1;
  var z;
  
  if (Array.isArray(mu)) z = arr.map((x, i) => mu_sign * (mu[i] - x) / sigma);
  else                   z = arr.map(x      => mu_sign * (mu    - x) / sigma);
  
  
  var vphi = varphi(z);
  
  if (return_sum)
  { 
    let order1 = mu_sign * sum(vphi.f) / sigma;
    let order2 = sum(vphi.fp) / sigma**2;
    
    return {order1: order1, order2: order2};
  }
  else
  {            
    let order1 = vphi.f.map(v => mu_sign * v / sigma);
    let order2 = vphi.fp.map(v => v / sigma**2)
    
    return {order1: order1, order2: order2};
  }
} // end of dklogPhi_dmuk


d2logPhi_dmudsigma = function(arr, mu, sigma, lower_tail=true, return_sum=true)
{
  // arr: an array
  // mu:  an array of same length as arr
  //   or a scalar
  // sigma: scalar
  
  
  const log_sigma = Math.log(sigma);
  const mu_sign = lower_tail ? -1 : 1;
  
  var z;
  
  if (Array.isArray(mu)) z = arr.map((x, i)      => mu_sign * (mu[i] - x) / sigma)
  else                   z = arr.map(x           => mu_sign * (mu    - x) / sigma);  
  
  var B = Phi(z, 0, 1, true, true);
  
  var logphi = phi(z, true); // array
  
  // First-order derivatives
  
  var Bp1 = as_log(z).mult(logphi).div(log_sigma).changed_sign(); // dlogPhi/dsigma
  var Bp2 = as_log(logphi, true, mu_sign).div(log_sigma);         // dlogPhi/dmu
  
  // Mixed derivative
  
  var Bm = as_log(z.map(z => z**2 - 1)).mult(logphi).div(2*log_sigma); // log-notation object
  if (mu_sign < 0) Bm.s = Bm.s.map(s => -s);
                               
  return d2LogB_dtheta1dtheta2(B, Bp1, Bp2, Bm, return_sum);         
} // end of d2logPhi_dmudsigma


dklogPhi_dsigmak = function(arr, mu, sigma, lower_tail=true)
{
  // arr: an array
  // mu: an array (of same length as arr) [for individual means]
  //     or a scalar
  // sigma: scalar
    
  const mu_sign = lower_tail ? -1 : 1;
 
  if (Array.isArray(mu))  z = arr.map((x, i) => mu_sign * (mu[i] - x)/sigma);
  else                    z = arr.map(x      => mu_sign * (mu    - x)/sigma);
  
  var vphi = varphi(z).f;
  
  order1 = - dot_product(z, vphi) / sigma;
  
  order2 = z.map((z, i) => vphi[i] * (z**3 + vphi[i]*z**2 - 2*z));
  order2 = - sum(order2) / sigma**2;
  
  return {order1: order1, order2: order2};
} // end of dklogPhi_dsigmak


dlogPhiInterval_dmu = function(interval, mu, sigma, return_sum=true)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu, sigma:  scalars
  
  const log_sigma = Math.log(sigma);
    
  var num = phi_diff(interval, mu, sigma, true);      // log-notation object
  var denom = PhiInterval(interval, mu, sigma, true); // array
      
  let tmp = num.div(denom).div(log_sigma).as_real(return_sum); // an array or a scalar
  
  if (return_sum) return -1 * tmp; // a number
  else return tmp.map(x => -x);    // an array
} // end of dlogPhiInterval_dmu


d2logPhiInterval_dmu2 = function(interval, mu, sigma, return_sum=true)
{
  // interval: a literal with dimensions a & b (two arrays of same length, the respective lower & upper interval endpoints)
  // mu: either a) an array of same length as interval.gt
  //         or b) a scalar
  // sigma:  a scalar
    
  const log_sigma = Math.log(sigma);
  
  var B = PhiInterval(interval, mu, sigma, true); // array 
  
  var tmp = dPhiInterval_dmu(interval, mu, sigma, log_sigma);    
  
  var Bs1 = as_log(tmp.z1).mult(tmp.logphi1).div(2*log_sigma); // in log-notation
  var Bs2 = as_log(tmp.z2).mult(tmp.logphi2).div(2*log_sigma);
    
  var Bs = Bs1.minus(Bs2);
  
  return d2LogB_dtheta2(B, tmp.Bp, Bs, return_sum); // scalar
} // end of d2logPhiInterval_dmu2


d2logPhiInterval_dmudsigma = function(interval, mu, sigma, return_sum=true)
{ 
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  //    mu: an array (of same length as interval.gt) OR a scalar
  // sigma:  scalar
  
  const log_sigma = Math.log(sigma);

  var B = PhiInterval(interval, mu, sigma, true); // array

  var tmp = dPhiInterval_dmu(interval, mu, sigma, log_sigma);
    var Bp_mu = tmp.Bp; // log-notation object
    
  
  var z2m1_1 = as_log(tmp.z1.map(z => z**2 - 1)).mult(tmp.logphi1);
  var z2m1_2 = as_log(tmp.z2.map(z => z**2 - 1)).mult(tmp.logphi2);
  
  var Bm = z2m1_1.minus(z2m1_2).div(2*log_sigma);
  
  var Bp_sigma = dPhiInterval_dsigma(interval, mu, sigma, log_sigma).Bp; // log-notation object

  return d2LogB_dtheta1dtheta2(B, Bp_mu, Bp_sigma, Bm, return_sum);
} // end of d2logPhiInterval_dmudsigma


dlogPhiInterval_dsigma = function(interval, mu, sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu, sigma:  scalars
  
  var B  =         PhiInterval(interval, mu, sigma, true); // array 
  var Bp = dPhiInterval_dsigma(interval, mu, sigma, Math.log(sigma)).Bp; // log-notation object 
  
  return dLogB_dtheta(B, Bp); // scalar  
} // end of dlogPhiInterval_dsigma
 
  
d2logPhiInterval_dsigma2 = function(interval, mu, sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu:    an array of same length as interval.gt & interval.lt
  //     or a scalar
  // sigma:  a scalar

  const log_sigma  = Math.log(sigma);
  
  var B = PhiInterval(interval, mu, sigma, true); // array
  
  var tmp = dPhiInterval_dsigma(interval, mu, sigma, log_sigma);

    var lz1 = tmp.lz1;
    var lz2 = tmp.lz2;


  var tmp1 = lz1.sq().minus(2).mult(tmp.logphi1).mult(lz1);
  var tmp2 = lz2.sq().minus(2).mult(tmp.logphi2).mult(lz2);  
    
  var Bs = tmp1.minus(tmp2).div(2*log_sigma);

  
  return d2LogB_dtheta2(B, tmp.Bp, Bs);
} // end of d2logPhiInterval_dsigma2


dkPhi_dmuk = function(arr, mu, sigma, lower_tail=true, log=false)
{
  // Returns the first two order derivatives
  // mu & sigma: two scalars
  
  var mu_sign = lower_tail ? -1 : 1;
  var phi, logphi;
  
  var z = arr.map(y => -mu_sign*(y-mu)/sigma);
  
  if (log) logphi = phi(z, true); // arrays
  else        phi = phi(z);     

  
  if (log)
  {
    let log_sigma = Math.log(sigma);
    order1 = {x: logphi.map(f => f - log_sigma), s: rep(mu_sign, this.length)};
    
    let logz = as_log(z);
    order2 = as_log(logphi.map(t => t-2*log_sigma), true, -1).mult(logz);
  }
  else
  {
    order1 = phi.map(f => mu_sign*f/sigma);
    order2 = phi.map((f, i) => -f * z[i]/sigma**2);
  }
  
  return {order1: order1, order2: order2, phi: phi, logphi: logphi}
} // end of dkPhi_dmuk


d2Phi_dmudsigma = function(arr, mu, sigma, lower_tail=true, log=false)
{
  // Returns the mixed derivative
  // mu & sigma: two scalars

  const log_sigma = Math.log(sigma);  
  var mu_sign = lower_tail ? -1 : 1;

  var z = arr.map(x => mu_sign * (mu - x)/sigma);
  var log_z2m1 = as_log(z.map(z => z**2 - 1));
  
  var logphi = phi(z, true); // array
  var dPhi2_dmudsigma = as_log(logphi.map(t => t - 2*log_sigma), true, mu_sign).mult(log_z2m1);
  
  if (log) return dPhi2_dmudsigma;
  else     return dPhi2_dmudsigma.as_real();
} // end of dPhi2_dmudsigma


dkPhi_dsigmak = function(arr, mu, sigma, lower_tail=true, log=false)
{
  // Returns the first two order derivatives
  // mu & sigma: two scalars
  
  const log_sigma = Math.log(sigma);  
  var mu_sign = lower_tail ? -1 : 1;
  
  var z = arr.map(y => -mu_sign*(y-mu)/sigma);
  var log_z_z2m2 = as_log(z.map(z => z*(z**2 - 2)));
  var logz = as_log(z);

  var phi = phi(z);
  var logphi = phi(z, true);
  
  var order1 = as_log(logphi.map(f => f -   log_sigma), true, -1).mult(logz);
  var order2 = as_log(logphi.map(f => f - 2*log_sigma), true, -1).mult(log_z_z2m2);
  
  if (!log)
  {
    order1 = order1.as_real();
    order2 = order2.as_real();
  }    

  return {order1: order1, order2: order2, z: z, phi: phi, logphi: logphi};
} // end of dkPhi_dsigmak


dPhiInterval_dmu = function(interval, mu, sigma, log_sigma)
{            
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu:       an array of same length as 'interval.gt' OR a scalar
  // sigma & log_sigma:  two scalars
   
  var z1, z2;
  
  if (Array.isArray(mu))
  {
    z1 = interval.gt.map((x, i) => (x-mu[i])/sigma);
    z2 = interval.lt.map((x, i) => (x-mu[i])/sigma);
  }
  else
  {     
    z1 = interval.gt.map(x => (x-mu)/sigma);
    z2 = interval.lt.map(x => (x-mu)/sigma);
  }
  
  
  var logphi1 = phi(z1, true); // arrays
  var logphi2 = phi(z2, true);
  
  var Bp1 = as_log(logphi1.map(l => l-log_sigma), true);
  var Bp2 = as_log(logphi2.map(l => l-log_sigma), true, -1);
  
  var Bp = Bp1.plus(Bp2);
  
  return {Bp: Bp, z1: z1, z2: z2, logphi1: logphi1, logphi2: logphi2};
} // end of dPhiInterval_dmu


dPhiInterval_dsigma = function(interval, mu, sigma, log_sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu:       an array of same length as 'interval.gt' OR a scalar
  // sigma & log_sigma:  two scalars

  var z1, z2;
  
  if (Array.isArray(mu))
  {
    z1 = interval.gt.map((x, i) => (x-mu[i])/sigma);
    z2 = interval.lt.map((x, i) => (x-mu[i])/sigma);
  }
  else
  {     
    z1 = interval.gt.map(x => (x-mu)/sigma);
    z2 = interval.lt.map(x => (x-mu)/sigma);
  }

    
  var logphi1 = phi(z1, true); // arrays
  var logphi2 = phi(z2, true);
  
  var lz1 = as_log(z1); // log-notation objects
  var lz2 = as_log(z2);
  
  var Bp1 = lz1.mult(logphi1);
  var Bp2 = lz2.mult(logphi2);

  var Bp = Bp1.minus(Bp2).div(log_sigma); // log-notation object
  
  return {Bp: Bp, lz1: lz1, lz2: lz2, logphi1: logphi1, logphi2: logphi2};
} // end of dPhiInterval_dsigma


//Object.prototype.dprod_dtheta = function(v, u_prime, v_prime)
//{
//  // this [u], v, u_prime & v_prime are either
//  //     a) all arrays (NOT on log-scale)
//  //  or b) all log-notation objects
//  
//  // IMPORTANT / WARNING
//  //  => This fct was not used nor tested yet
//  
//  if (Array.isArray(this))
//  {
//    return u_prime.times(v).plus(this.times(v_prime));
//  }
//  else if (Array.isArray(this.x))
//  {
//    var t1 = {x: u_prime.x.plus(v.x),    s: u_prime.s.times(v.s)};
//    var t2 = {x: this.x.plus(v_prime.x), s: this.s.times(v_prime.s)};
//    
//    return t1.plus(t2);
//  } 
//  else
//  {
//    // this.x is a number
//    
//    var t1 = {x: u_prime.x + v.x,     s: u_prime.s * v.s};
//    var t2 = {x: this.x + v_prime.x,  s: this.s * v_prime.s};
//    
//    return t1.plus(t2);
//  } 
//} // end of Object.dprod_dtheta


//Object.prototype.dratio_dtheta = function(v, u_prime, v_prime)
//{
//  // this [u], v, u_prime & v_prime are either
//  //     a) all arrays (NOT on log-scale)
//  //  or b) all log-notation objects
//  
//  // IMPORTANT / WARNING
//  //  => This fct was not used nor tested yet
//  
//  if (Array.isArray(this))
//  {
//    return u_prime.times(v).minus(this.times(v_prime)).divided_by(v.sq());
//  }
//  else if (Array.isArray(this.x))
//  {
//    var up_v = {x: u_prime.x.plus(v.x),    s: u_prime.s.times(v.s)};
//    var u_vp = {x: this.x.plus(v_prime.x), s: this.s.times(v_prime.s)};
//    
//    return up_v.minus(u_vp).div(v.x, 2);
//  }  
//  else
//  {
//    // this.x is a number
//    var up_v = {x: u_prime.x + v.x,     s: u_prime.s * v.s};
//    var u_vp = {x: this.x + v_prime.x,  s: this.s * v_prime.s};
//    
//    var tmp = up_v.minus(u_vp);
//    
//    return Math.exp(tmp.x - 2*v.x) * tmp.s;
//  }
//} // end of Object.dratio_dtheta