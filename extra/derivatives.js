// Functions to calculate derivatives often involved in log(posterior)
//
// Author: Patrick Bélisle
//
// Version 0.2 (Apr 2021)                                              


// Change Log
// ========================
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



Array.prototype.dlogPhi_dmu = function(mu, sigma, lower_tail=true, return_sum=true)
{
  // this: an array, with values y_i
  // mu & sigma: scalars
  // If return_sum = true, return the sum of the terms
  //   otherwise return an array with all the terms

  var mu_sign = lower_tail ? -1 : 1;
  var z = this.map(x => mu_sign * (mu - x)/sigma);
  
  var varphi = z.varphi().f;
  
  if (return_sum) return mu_sign * varphi.sum() / sigma;
  else            return varphi.times(mu_sign/sigma);
} // end of Array.dlogPhi_dmu


Array.prototype.d2logPhi_dmu2 = function(mu, sigma, lower_tail=true, return_sum=true)
{
  // this: an array, with values y_i
  // mu & sigma: scalars
  // If return_sum = true, return the sum of the terms
  //   otherwise return an array with all the terms

  var mu_sign = lower_tail ? -1 : 1;
  var z = this.map(x => mu_sign * (mu - x)/sigma);
  
  var vphi_prime = z.varphi().fp;
  
  if (return_sum) return vphi_prime.sum() / sigma**2;
  else            return vphi_prime.divided_by(sigma**2);
  
  //if (return_sum) return mu_sign * vphi_prime.sum() / sigma**2;
  //else            return vphi_prime.divided_by(mu_sign * sigma**2);
} // end of Array.d2logPhi_dmu2


Array.prototype.d2logPhi_dmudsigma = function(mu, sigma, lower_tail=true)
{
  // this: an array
  // mu, sigma:  scalars
  
  const log_sigma = Math.log(sigma);
  const mu_sign = lower_tail ? -1 : 1;
  
  var z = this.map(x => mu_sign * (mu - x)/sigma);
  var B = z.Phi(0, 1, true, true);
  
  var logphi = z.phi(true); // array
  
  // First-order derivatives
  
  var lz = z.as_log();
  
  var Bp1 = {x: lz.x.plus(logphi).minus(log_sigma), s: lz.s.times(-1)};  // dlogPhi/dsigma
  var Bp2 = {x: logphi.minus(log_sigma), s: [mu_sign].rep(this.length)}; // dlogPhi/dmu
  
  // Mixed derivative
  
  var Bm = z.sq().minus(1).as_log().lo_mult_(logphi).lo_div_(2*log_sigma); // log-notation object 
  if (mu_sign < 0) Bm.s = Bm.s.times(-1);
                               
  return d2sumLogB_dtheta1dtheta2(B, Bp1, Bp2, Bm);         
} // end of Array.d2logPhi_dmudsigma


Array.prototype.dlogPhi_dsigma = function(mu, sigma, lower_tail=true)
{
  // this: an array
  // mu, sigma: scalars
  
  // NOTE: there is another version of this fct available at the bottom of this file (used by MAPInits-BW.js)
  
  const mu_sign = lower_tail ? -1 : 1;
  
  z = this.map(x => mu_sign * (mu - x)/sigma);
  var vphi = z.varphi().f;
  
  return - z.dot_product(vphi) / sigma;
} // end of Array.dlogPhi_dsigma


Array.prototype.d2logPhi_dsigma2 = function(mu, sigma, lower_tail=true)
{
  // this: an array
  // mu, sigma: scalars
  
  // NOTE: There is an alternative version of this fct (used by MAPInits-BW.js) at the bottom of this file
  
  const mu_sign = lower_tail ? -1 : 1;
  
  z = this.map(x => -mu_sign * (x-mu)/sigma);
  var vphi = z.varphi().f;
  
  var a = z.times(vphi);
  var b = a.plus(z.sq()).map(x => x - 2);
  
  return -1 * a.dot_product(b) / sigma**2;
} // end of Array.d2logPhi_dsigma2


dlogPhiInterval_dmu = function(interval, mu, sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu, sigma:  scalars
    
  var num = phi_diff(interval, mu, sigma, true); // log-notation object
  var denom = PhiInterval(interval, mu, sigma, true); // array
  
  var l = num.x.minus(denom);
    
  return -1 * l.sum_signedExp(num.s) / sigma;
} // end of dlogPhiInterval_dmu


d2logPhiInterval_dmu2 = function(interval, mu, sigma)
{
  // interval: a literal with dimensions a & b (two arrays of same length, the respective lower & upper interval endpoints)
  // mu, sigma:  scalars
    
  const log_sigma = Math.log(sigma);
  
  var B = PhiInterval(interval, mu, sigma, true); // array 
  
  var tmp = dPhiInterval_dmu(interval, mu, sigma, log_sigma);
  
    var Bp = tmp.Bp;
    var z1 = tmp.z1;
    var z2 = tmp.z2;
    var logphi1 = tmp.logphi1;
    var logphi2 = tmp.logphi2;
    
  
  var Bs1 = z1.as_log().lo_mult_(logphi1).lo_div_(2*log_sigma); // in log-notation
  var Bs2 = z2.as_log().lo_mult_(logphi2).lo_div_(2*log_sigma);
  
// ICI now included in the two lines above
//  Bs1.x = Bs1.x.plus(logphi1);
//  Bs2.x = Bs2.x.plus(logphi2);
  
//  Bs1.x = Bs1.x.map(x => x - 2*log_sigma);
//  Bs2.x = Bs2.x.map(x => x - 2*log_sigma);
  
  var Bs = Bs1.minus_log(Bs2);
  
  return d2sumLogB_dtheta2(B, Bp, Bs); // scalar
} // end of d2logPhiInterval_dmu2


d2logPhiInterval_dmudsigma = function(interval, mu, sigma)
{ 
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  //    mu: an array (of same length as interval.gt) OR a scalar
  // sigma:  scalar
  
  // NOTE: There is an alternative version of this fct at the bottom of this file (it is used by MAPInits-BW.js)
  
  const log_sigma = Math.log(sigma);

  var B = PhiInterval(interval, mu, sigma, true); // array

  var tmp = dPhiInterval_dmu(interval, mu, sigma, log_sigma);
  
    var Bp_mu = tmp.Bp; // log-notation object
    
    var z1 = tmp.z1;  // arrays
    var z2 = tmp.z2;
    var logphi1 = tmp.logphi1;
    var logphi2 = tmp.logphi2;  
  
  
  var z2m1_1 = z1.sq().minus(1).as_log().lo_mult_(logphi1);
  var z2m1_2 = z2.sq().minus(1).as_log().lo_mult_(logphi2);
  
  var Bm = z2m1_1.minus_log(z2m1_2).lo_div_(2*log_sigma);
  
  var Bp_sigma = dPhiInterval_dsigma(interval, mu, sigma, log_sigma).Bp; // log-notation object

  
  return d2sumLogB_dtheta1dtheta2(B, Bp_mu, Bp_sigma, Bm);
} // end of d2logPhiInterval_dmudsigma


dlogPhiInterval_dsigma = function(interval, mu, sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  // mu, sigma:  scalars
  
  // NOTE: There is an alterative version of this fct at the bottom of this file
  //       [it is used by MAPInits-BW.js]
  
  var B  =         PhiInterval(interval, mu, sigma, true); // array 
  var Bp = dPhiInterval_dsigma(interval, mu, sigma, Math.log(sigma)).Bp; // log-notation object 
  
  return dsumLogB_dtheta(B, Bp); // scalar  
} // end of dlogPhiInterval_dsigma
 
  
d2logPhiInterval_dsigma2 = function(interval, mu, sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  
  // mu, sigma:  scalars

  const log_sigma  = Math.log(sigma);
  
  var B = PhiInterval(interval, mu, sigma, true); // array
  
  var tmp = dPhiInterval_dsigma(interval, mu, sigma, log_sigma);
  
    var Bp = tmp.Bp;
    var lz1 = tmp.lz1;
    var lz2 = tmp.lz2;
    var logphi1 = tmp.logphi1;
    var logphi2 = tmp.logphi2;

    
  var tmp1 = lz1.x.times(2).as_log(true).minus_log(2).lo_mult_(logphi1).lo_mult(lz1);
  var tmp2 = lz2.x.times(2).as_log(true).minus_log(2).lo_mult_(logphi2).lo_mult(lz2);
  
// ICI maintenant inclues dans les 2 lignes ci-dessus  
//  tmp1.x = tmp1.x.plus(logphi1); 
//  tmp2.x = tmp2.x.plus(logphi2);
  
//  tmp1 = tmp1.lo_mult(lz1);
//  tmp2 = tmp2.lo_mult(lz2);
  
  var Bs = tmp1.minus_log(tmp2).lo_div_(2*log_sigma);
  // Bs.x = Bs.x.minus(2*log_sigma);  ICI maintenant inclus dans ligne ci-dessus
  
  return d2sumLogB_dtheta2(B, Bp, Bs);
  
} // end of d2logPhiInterval_dsigma2


Array.prototype.dkPhi_dmuk = function(mu, sigma, lower_tail=true, log=false)
{
  // Returns the first two order derivatives
  // mu & sigma: two scalars
  
  var mu_sign = lower_tail ? -1 : 1;
  var phi, logphi;
  
  var z = this.map(y => -mu_sign*(y-mu)/sigma);
  
  if (log) logphi = z.phi(true); // arrays
  else        phi = z.phi();     

  
  if (log)
  {
    let log_sigma = Math.log(sigma);
    order1 = {x: logphi.map(f => f - log_sigma), s: [mu_sign].rep(this.length)};
    
    let logz = z.as_log();
    order2 = logphi.map(t => t-2*log_sigma).as_log(true, false).lo_mult(logz);
  }
  else
  {
    order1 = phi.map(f => mu_sign*f/sigma);
    order2 = phi.times(z).map(t => -t/sigma**2);
  }
  
  return {order1: order1, order2: order2, phi: phi, logphi: logphi}
} // end of Array.dkPhi_dmuk


Array.prototype.d2Phi_dmudsigma = function(mu, sigma, lower_tail=true, log=false)
{
  // Returns the mixed derivative
  // mu & sigma: two scalars

  const log_sigma = Math.log(sigma);  
  var mu_sign = lower_tail ? -1 : 1;

  var z = this.map(x => mu_sign * (mu - x)/sigma);
  var log_z2m1 = z.map(z => z**2 - 1).as_log();
  
  var logphi = z.phi(true); // array
  var dPhi2_dmudsigma = logphi.map(t => t - 2*log_sigma).as_log(true, !lower_tail).lo_mult(log_z2m1);

  
  if (log) return dPhi2_dmudsigma;
  else     return dPhi2_dmudsigma.as_real();
} // end of Array.dPhi2_dmudsigma


Array.prototype.dkPhi_dsigmak = function(mu, sigma, lower_tail=true, log=false)
{
  // Returns the first two order derivatives
  // mu & sigma: two scalars
  
  const log_sigma = Math.log(sigma);  
  var mu_sign = lower_tail ? -1 : 1;
  
  var z = this.map(y => -mu_sign*(y-mu)/sigma);
  var log_z_z2m2 = z.map(z => z*(z**2 - 2)).as_log();
  var logz = z.as_log();

  var phi = z.phi();
  var logphi = z.phi(true);
  
  
  var order1 = logphi.map(f => f - log_sigma).as_log(true, false).lo_mult(logz);
  
  var order2 = logphi.map(f => f - 2*log_sigma).as_log(true, false).lo_mult(log_z_z2m2);
  
  
  if (!log)
  {
    order1 = order1.as_real();
    order2 = order2.as_real();
  }    

  return {order1: order1, order2: order2, z: z, phi: phi, logphi: logphi};
} // end of Array.dkPhi_dsigmak


dPhiInterval_dmu = function(interval, mu, sigma, log_sigma)
{            
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  //         mu:         an array of same length as 'this' & b, OR a scalar
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
  
  
  var logphi1 = z1.phi(true); // arrays
  var logphi2 = z2.phi(true);
  
  var Bp1 = {x: logphi1.map(l => l-log_sigma), s:  [1].rep(interval.gt.length)};
  var Bp2 = {x: logphi2.map(l => l-log_sigma), s: [-1].rep(interval.gt.length)};
  
  var Bp = Bp1.plus_log(Bp2);
  
  return {Bp: Bp, z1: z1, z2: z2, logphi1: logphi1, logphi2: logphi2};
} // end of dPhiInterval_dmu


dPhiInterval_dsigma = function(interval, mu, sigma, log_sigma)
{
  // interval: a literal with dimensions gt & lt [two arrays of same length, the lower & upper interval endpoints, respectively]
  //         mu:         an array of same length as 'this' & b, OR a scalar
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

    
  var logphi1 = z1.phi(true); // arrays
  var logphi2 = z2.phi(true);
  
  var lz1 = z1.as_log(); // log-notation objects
  var lz2 = z2.as_log();
  
  var Bp1 = {x: lz1.x.plus(logphi1), s: lz1.s};
  var Bp2 = {x: lz2.x.plus(logphi2), s: lz2.s};
  
  var Bp = Bp1.minus_log(Bp2).lo_div_(log_sigma); // log-notation object
  // Bp.x = Bp.x.minus(log_sigma); // ICI maintenant inclus dans ligne ci-dessus
  
  return {Bp: Bp, lz1: lz1, lz2: lz2, logphi1: logphi1, logphi2: logphi2};
} // end of dPhiInterval_dsigma


Object.prototype.dprod_dtheta = function(v, u_prime, v_prime)
{
  // this [u], v, u_prime & v_prime are either
  //     a) all arrays (NOT on log-scale)
  //  or b) all log-notation objects
  
  // IMPORTANT / WARNING
  //  => This fct was not used nor tested yet
  
  if (Array.isArray(this))
  {
    return u_prime.times(v).plus(this.times(v_prime));
  }
  else if (Array.isArray(this.x))
  {
    var t1 = {x: u_prime.x.plus(v.x),    s: u_prime.s.times(v.s)};
    var t2 = {x: this.x.plus(v_prime.x), s: this.s.times(v_prime.s)};
    
    return t1.plus_log(t2);
  } 
  else
  {
    // this.x is a number
    
    var t1 = {x: u_prime.x + v.x,     s: u_prime.s * v.s};
    var t2 = {x: this.x + v_prime.x,  s: this.s * v_prime.s};
    
    return t1.plus_log(t2);
  } 
} // end of Object.dprod_dtheta


Object.prototype.dratio_dtheta = function(v, u_prime, v_prime)
{
  // this [u], v, u_prime & v_prime are either
  //     a) all arrays (NOT on log-scale)
  //  or b) all log-notation objects
  
  // IMPORTANT / WARNING
  //  => This fct was not used nor tested yet
  
  if (Array.isArray(this))
  {
    return u_prime.times(v).minus(this.times(v_prime)).divided_by(v.sq());
  }
  else if (Array.isArray(this.x))
  {
    var up_v = {x: u_prime.x.plus(v.x),    s: u_prime.s.times(v.s)};
    var u_vp = {x: this.x.plus(v_prime.x), s: this.s.times(v_prime.s)};
    
    return up_v.minus_log(u_vp).lo_div_(v.x, 2);
  }  
  else
  {
    // this.x is a number
    var up_v = {x: u_prime.x + v.x,     s: u_prime.s * v.s};
    var u_vp = {x: this.x + v_prime.x,  s: this.s * v_prime.s};
    
    var tmp = up_v.minus_log(u_vp);
    
    return Math.exp(tmp.x - 2*v.x) * tmp.s;
  }
} // end of Object.dratio_dtheta


dsumLogB_dtheta = function(B, Bp)
{
  // B & Bp are either
  //   -- both arrays (on natural scale) or
  //   -- B is an array (on log-scale) and Bp is a log-notation objects

  // Returns d(sum(log(B_i)))/dtheta 
  //   given that the               B_i   terms are given in B, 
  //   and their first derivatives  B_i'  are given in Bp

  //  d(log(B))/dtheta = B'/B 

  
  if (Array.isArray(Bp))
  { 
    // B, Bp & Bs are arrays
    return Bp.division_sum(B);
  }
  else
  {
    // Bp is a log-notation object, B an array
    return Bp.x.minus(B).sum_signedExp(Bp.s);
  }
} // end of dsumLogB_dtheta


d2sumLogB_dtheta2 = function(B, Bp, Bs)
{
  // B, Bp, Bs are either
  //   -- all arrays (on natural scale) or
  //   -- B is an array (on log-scale) and Bp & Bs are log-notation objects

  // Returns d2(sum(log(B_i)))/dtheta2 
  //   given that the            B_i terms are given in B, 
  //   their first derivatives   B_i'  given in Bp
  //   and their 2nd derivatives B_i'' given in Bs
  
  //  d2(log(B))/dtheta2 = d/dtheta (B'/B) = = (B''B - B'**2)/B**2 = B''/B - (B'/B)**2
  
    
  if (Array.isArray(Bp))
  { 
    // B, Bp & Bs are arrays
    return Bs.division_sum(B) - Bp.division_sqSum(B);
  }
  else
  {
    // Bs is a log-notation object, B & Bp are arrays
    return Bs.x.minus(B).sum_signedExp(Bs.s) - Bp.x.division_sqSum(B, true);
  }
} // end of d2sumLogB_dtheta2


d2sumLogB_dtheta1dtheta2 = function(B, Bp1, Bp2, Bm)
{
  // B is an array
  // all other arguments are either 
  //  - arrays (on natural scale, as well as B)
  //  - log-notation objects (and then B is an array on log scale)
  
 
  if (Array.isArray(Bp1))
  {
    return Bm.division_sum(B) - Bp1.times(Bp2).division_sum(B.sq());
  }
  else
  {
    let Bm1 = Bm.x.minus(B).sum_signedExp(Bm.s);              // scalar    
    let Bm2 = Bp1.lo_mult(Bp2).lo_div_(B, 2).sum_signedExp(); // scalar

    // ICI code un peu efface, maintenant inclus dans les deux lignes ci-dessus
//    Bm2.x = Bm2.x.minus(B.times(2));   
//    Bm2 = Bm2); // scalar
    
    return Bm1 - Bm2;
  } 
} // end of d2sumLogB_dtheta1dtheta2



// *!*!*!*!*!*!*!*!*!*!**!*!*!*!*!*!*!*!*!*!*!*!*!*!*!**!*!*!*!*!*!*!*!*!*!*!*!*
//
// Functions that are used in MAPInits-BW.js


function dlogPhi_dsigma(a, sigma)
{
  // a: array
  // sigma: scalar
  // was named g2
  
  z = a.map(z => z/sigma);
  var vphi = z.varphi();
  
  return z.dot_product(vphi.f) / sigma;
} // end of dlogPhi_dsigma


function d2logPhi_dsigma2(z, sigma)
{
  // z: array
  // sigma: scalar
  // was named dg2_dsigma
  
  var z = z.map(z => z/sigma);
  var vphi = z.varphi();
  
  var z2 = z.sq();
  var d2logPhi_dsigma2 = - (z2.dot_product(vphi.fp) + 2*z.dot_product(vphi.f)) / sigma**2;
  
  return d2logPhi_dsigma2;
} // end of d2logPhi_dsigma2