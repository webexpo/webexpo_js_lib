// Functions to calculate derivatives often involved in log(posterior)
//
// Author: Patrick Bélisle
//
// Version 0.1 (Mar 2021)


// ICI importer (de normal-OneSubjectEstimates.js) les fonctions suivantes
// en les rendant des Array.prototype et en les validant/testant au fur et a mesure
//   -> 4 choses a tester left_tail=T|F  * return_sum=T|F
// De meme, pour les fcts djlogPhi_dtheta, generaliser la fct avec l'option left_tail=false
//

//   dlogPhi_dsigma
//   d2logPhi_dsigma2
//   dlogPhiInterval_dsigma
//   d2logPhiInterval_dsigma2


Array.prototype.dlogPhi_dmu = function(mu, sigma, left_tail=true, return_sum=true)
{
  // this: an array, with values y_i
  // mu & sigma: scalars
  // If return_sum = true, return the sum of the terms
  //   otherwise return an array with all the terms

  var mu_sign = left_tail ? -1 : 1;
  var z = this.map(x => -mu_sign*(x-mu)/sigma);
  
  var varphi = z.varphi().f;
  
  if (return_sum) return mu_sign * varphi.sum() / sigma;
  else            return varphi.times(mu_sign/sigma);
} // end of Array.dlogPhi_dmu


Array.prototype.d2logPhi_dmu2 = function(mu, sigma, left_tail=true, return_sum=true)
{
  // this: an array, with values y_i
  // mu & sigma: scalars
  // If return_sum = true, return the sum of the terms
  //   otherwise return an array with all the terms

  var mu_sign = left_tail ? -1 : 1;
  var z = this.map(x => -mu_sign*(x-mu)/sigma);
  
  var vphi_prime = z.varphi().fp;
  
  if (return_sum) return mu_sign * vphi_prime.sum() / sigma**2;
  else            return vphi_prime.divided_by(mu_sign * sigma**2);
} // end of Array.d2logPhi_dmu2


Array.prototype.dlogPhiInterval_dmu = function(b, mu, sigma)
{
  // this: an array with lower interval endpoints
  //    b: an array (of same length as 'this') with upper interval endpoints
  // mu, sigma:  scalars
    
  var num = b.phi_diff(this, mu, sigma, true);      // log-notation object
  var denom = this.PhiInterval(b, mu, sigma, true); // array
  
  var l = num.x.minus(denom);
  
  return -1 * l.sum_signedExp(num.s) / sigma;
} // end of dlogPhiInt_dmu


Array.prototype.d2logPhiInterval_dmu2 = function(b, mu, sigma)
{
  // this: an array with lower interval endpoints
  //    b: an array (of same length as 'this') with upper interval endpoints
  // mu, sigma:  scalars
    
  const log_sigma = Math.log(sigma);
  
  var B = this.PhiInterval(b, mu, sigma, true); // array 
  
  var z1 = this.map(x => (x-mu)/sigma);
  var z2 =    b.map(x => (x-mu)/sigma);
  
  var logphi1 = z1.phi(true); // arrays
  var logphi2 = z2.phi(true);
  
  
  var Bp1 = {x: logphi1.map(l => l-log_sigma), s:  [1].rep(b.length)};
  var Bp2 = {x: logphi2.map(l => l-log_sigma), s: [-1].rep(b.length)};
  
  var Bp = Bp1.plus_log(Bp2);
  
  
  var Bs1 = log_notation(z1); // in log-notation
  var Bs2 = log_notation(z2);
  
  Bs1.x = Bs1.x.plus(logphi1);
  Bs2.x = Bs2.x.plus(logphi2);
  
  Bs1.x = Bs1.x.map(x => x - 2*log_sigma);
  Bs2.x = Bs2.x.map(x => x - 2*log_sigma);
  
  var Bs = Bs1.minus_log(Bs2);
  
  return d2sumLogB_dtheta2(B, Bp, Bs); // scalar
} // end of d2logPhiInterval_dmu2


dsumLogB_dtheta = function(B, Bp)
{
  // B & Bp are either
  //   -- abothll arrays (on natural scale) or
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
    return Bs.division_sum(B) + Bp.x.division_sqSum(B);
  }
  else
  {
    // Bs is a log-notation object, B & Bp are arrays
    
    return Bs.x.minus(B).sum_signedExp(Bs.s) - Bp.x.division_sqSum(B, true);
  }
} // end of d2sumLogB_dtheta2