// Code to find initial values in the context of Normal (or log-Normal observations)
// in presence of censored observations
//
// Author: Patrick Bélisle

// Version 0.1 (Mar 2021)


function d2logPhi_dmu2(z, sigma, mu_sign=-1)
{
  // z: array
  // sigma: scalar
  // was name dg1_dmu
  
  let vphi = z.map(z => z/sigma).varphi();
  
  return mu_sign * vphi.fp.sum() / sigma**2;
} // end of d2logPhi_dmu2


function d2logPhi_dmudsigma(z, sigma, mu_sign=-1)
{
  // z: array
  // sigma: scalar
  // was named dg2_dmu
  
  var  z = z.map(z => z/sigma);
  var vphi = z.varphi();
  
  return mu_sign * (vphi.f.sum() + z.dot_product(vphi.fp)) / sigma**2;
} // end of d2logPhi_dmudsigma


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


function d2logPhiInt_dmu2(interval, mu, sigma)
{
  // was named dgd_dmu
  
  var z1 = interval.gt.minus(mu),
      z2 = interval.lt.minus(mu);
  
  var z1num = z1,
      z2num = z2;
  
  z1 = z1.divided_by(sigma);
  z2 = z2.divided_by(sigma);
  
  var log_sigma = Math.log(sigma);
  
  var log_u = log_dphi(z1, z2);
  var log_v = log_dPhi(z1, z2);  // array
  
  var log_u1 = z1.phi(true),
      log_u2 = z2.phi(true);
  
  var log_up1 = {x: [], s: []},
      log_up2 = {x: [], s: []};
  
  log_up1.x = log_u1.plus(z1.log_abs()).minus(log_sigma);  
  log_up2.x = log_u2.plus(z2.log_abs()).minus(log_sigma);
  
  log_up1.s = z1.sign();
  log_up2.s = z2.sign();
    
  var log_vp1 = dPhi_dmu(z1num, sigma, true),
      log_vp2 = dPhi_dmu(z2num, sigma, true);
  
  var log_up = log_diff(log_up1, log_up2),
      log_vp = log_diff(log_vp1, log_vp2);
  
  var log_d1 = {x: [], s: log_up.s};
  var log_d2 = {x: [], s: log_u.s.times(log_vp.s)};
  
      log_d1.x = log_up.x.plus(log_v);
      log_d2.x = log_u.x.plus(log_vp.x);

  
  var log_d = log_diff(log_d2, log_d1);
      log_d.x.minus(log_v.times(2));
      
  var d = log_d.x.map(Math.exp).times(log_d.s);
  
  return d.sum() / sigma; // scalar
} // end of d2logPhiInt_dmu2


function d2logPhiInt_dmudsigma(interval, mu, sigma)
{
  // was named dgm_dmu
  
  var x1 = interval.gt.minus(mu),
      x2 = interval.lt.minus(mu);

  var z1 = x1.divided_by(sigma),
      z2 = x2.divided_by(sigma);
      
  var u1 = dphi_dsigma(x1, sigma, true),
      u2 = dphi_dsigma(x2, sigma, true);
      
  var u = log_diff(u1, u2);
  var v = log_dPhi(z1, z2);

  var up1 = d2phi_dmudsigma(x1, sigma, true),
      up2 = d2phi_dmudsigma(x2, sigma, true);

  var u_prime = log_diff(up1, up2);

  var vp1 = dPhi_dmu(x1, sigma, true),
      vp2 = dPhi_dmu(x2, sigma, true);

  var v_prime = log_diff(vp1, vp2);

  var d1 = log_mult(u, v_prime);

  var d2 = {x: u_prime.x.plus(v), s: u_prime.s};

  var d = log_diff(d1, d2);
      d = d.x.minus(v.times(2)).map(Math.exp).times(d.s);

  return d.sum();
} // end of d2logPhiInt_dmudsigma


function d2logPhiInt_dmudsigma(interval, mu, sigma)
{
  // was named dgm_dsigma
  
  var x1 = interval.gt.minus(mu),
      x2 = interval.lt.minus(mu);

  var z1 = x1.divided_by(sigma),
      z2 = x2.divided_by(sigma);

  var log_u1 = dphi_dsigma(x1, sigma, true),
      log_u2 = dphi_dsigma(x2, sigma, true);

  var log_u = log_diff(log_u1, log_u2);
  var log_v = log_dPhi(z1, z2);

  var log_up1 = d2phi_dsigma2(x1, sigma, true),
      log_up2 = d2phi_dsigma2(x2, sigma, true);

  var log_up = log_diff(log_up1, log_up2);

  var log_vp1 = dPhi_dsigma(x1, sigma, true),
      log_vp2 = dPhi_dsigma(x2, sigma, true);

  var log_vp = log_diff(log_vp1, log_vp2);

  var log_d1 = {x: log_up.x.plus(log_v), s: log_up.s},
      log_d2 = log_mult(log_u, log_vp);

  var log_d = log_diff(log_d2, log_d1);
      log_d.x = log_d.x.minus(log_v.times(2));

  var d = log_d.x.map(Math.exp).times(log_d.s);

  return d.sum();
} // end of d2logPhiInt_dmudsigma


function d2phi_dmudsigma(a, sigma, log=false)
{
  // a: array
  // sigma: scalar
  // was named dgm1_dmu

  var d2phi_dmudsigma;
  var z = a.divided_by(sigma);

  if (log)
  {
    var b = z.sq().minus(1); // array
    var x = z.phi(true).plus(b.log_abs()).minus(2*Math.log(sigma));

    d2phi_dmudsigma = {x: x, s: b.sign()};
  }
  else
  {
    var z2m1 = z.sq().minus(1);

    d2phi_dmudsigma = z.phi().times(z2m1).divided_by(sigma**2);
  }


  return d2phi_dmudsigma;
} // end of d2phi_dmudsigma


function d2phi_dsigma2(a, sigma, log=false)
{
  // was named dgm1_dsigma
  var d2phi_dsigma2;

  var z = a.divided_by(sigma);
  var z2m2 = z.sq().minus(2);


  if (log)
  {
    var log_sigma = Math.log(sigma);
    var x =  z.phi(true).plus(z.log_abs()).plus(z2m2.log_abs()).minus(2*log_sigma);

    d2phi_dsigma2 = {x: x, s: z.sign().times(z2m2.sign())};
  }
  else
  {
    d2phi_dsigma2 = z.phi().times(z).times(z2m2).divided_by(sigma**2);
  }


  return d2phi_dsigma2;
} // end of d2phi_dsigma2


function dPhi_dmu(a, sigma, log=false)
{
  // a: array
  // sigma: scalar
  
  var dPhi_dmu;
  var z = a.divided_by(sigma);
  
  if (log)
  {
    var x = z.phi(true).minus(Math.log(sigma));
        
    dPhi_dmu = {x: x, s: [-1].rep(z.length)};
  }
  else
  {
    dPhi_dmu = z.phi().map(u => -u/sigma);    
  }

  return dPhi_dmu;
} // end of dPhi_dmu


function dPhi_dsigma(a, sigma, log=false)
{
  // a: array
  // sigma: scalar

  var dPhi_dsigma;

  var z = a.divided_by(sigma);

  if (log)
  {
    var log_sigma = Math.log(sigma);
    var x = z.phi(true).plus(a.log_abs()).minus(2*log_sigma);

    dPhi_dsigma = {x: x, s: a.sign().times(-1)};
  }
  else
  {
    dPhi_dsigma = z.phi().times(z).map(u => -u/sigma);
  }


  return dPhi_dsigma;
} // end of dPhi_dsigma


function dlogPhi_dmu(z, sigma)
{
  // z: array
  // sigma: scalar
  // was named g1
  
  var vphi = z.map(z => z/sigma).varphi();
  
  return vphi.f.sum() / sigma;
} // end of dlogPhi_dmu


function dlogPhi_dsigma(a, sigma)
{
  // z: array
  // sigma: scalar
  // was named g2
  
  z = a.map(z => z/sigma);
  var vphi = z.varphi();
  
  return z.dot_product(vphi.f) / sigma;
} // end of dlogPhi_dsigma


function dlogPhiInt_dmu(interval, mu, sigma)
{
  // interval:   object with vectors {l, u}, of same length
  // mu, sigma:  scalars
  // was named gd
  
  var z1 = interval.gt.minus(mu).divided_by(sigma),
      z2 = interval.lt.minus(mu).divided_by(sigma);

  var num1_log = z1.phi(true),
      num2_log = z2.phi(true);
  
  var   num_log = log_diff_exp(num1_log, num2_log);
  var denom_log = log_dPhi(z1, z2);
  
  var log_ratios = num_log.x.minus(denom_log);
  var     ratios = log_ratios.map(Math.exp);

  ratios = ratios.times(num_log.s);  
  
  return ratios.sum() / sigma;
} // end of dlogPhiInt_dmu


function dlogPhiInt_dsigma(interval, mu, sigma)
{
  // was named gm
  
  var z1 = interval.gt.minus(mu)
      z2 = interval.lt.minus(mu);
  
  var u1 = dphi_dsigma(z1, sigma, true),
      u2 = dphi_dsigma(z2, sigma, true);
  
  var u = log_diff(u1, u2);
  
  z1 = z1.divided_by(sigma);
  z2 = z2.divided_by(sigma);
  
  var log_v = log_dPhi(z1, z2); // array
  var log_ratios = u.x.minus(log_v);
  var ratios = log_ratios.map(Math.exp);
  ratios = ratios.times(u.s);  
  
  return ratios.sum(); // scalar
} // end of dlogPhiInt_dsigma


function dphi_dsigma(a, sigma, log=false)
{
  // a: array
  // sigma : scalar
  // was named gm1
  
  var dphi_dsigma;
  var z = a.divided_by(sigma);
  
  if (log)
  {
    var phi = z.phi(true).plus(a.log_abs()).minus(2*Math.log(sigma));
    
    dphi_dsigma = {x: phi, s: a.sign()};
  }
  else
  {
    var phi = z.phi();
    dphi_dsigma = phi.times(z).divided_by(sigma);
  }
  
  // when log = TRUE, return a list with dimensions (x, s)
  // when log = FALSE, return an array

  return dphi_dsigma;
} // end of dphi_dsigma


function log_dphi(z1, z2)
{
  var log1 = z1.phi(true),
      log2 = z2.phi(true);
  
  return log_diff_exp(log1, log2);
} // end of log_dphi


function log_dPhi(z1, z2)
{  
  // z1, z2: arrays of same length
  
  var denom1_log = z1.Phi(true, true),
      denom2_log = z2.Phi(true, true);
  
  var denom_log = log_diff_exp(denom1_log, denom2_log);
  
  var w =  denom_log.x.whichInf();
  
  if (w.length > 0)
  {
    // Take the problem from the opposite side (for thosw that return a result = Infinity)
    var opposite = {z1: [], z2: []};
    w.forEach(function(w){opposite.z1.push(z1[w]); opposite.z2.push(z2[w]);});
    
    var opposite_denom1_log = opposite.z1.Phi(false, true),
        opposite_denom2_log = opposite.z2.Phi(false, true);
    
    var opposite_denom_log = log_diff_exp(opposite_denom1_log, opposite_denom2_log);
    
    w.forEach(function(w, j){denom_log.x[w] = opposite_denom_log.x[j]});
  }
  
  return denom_log.x; // array
} // end of log_dPhi


Array.prototype.varphi = function()
{
  // this: array
  
  var log_phi = this.phi(true);
  var log_Phi = this.Phi(true, true);
  
  var vphi       = log_phi.minus(log_Phi).map(Math.exp);
  var vphi_prime = this.times(vphi).plus(vphi.sq()).times(-1);
  
  var varphi = {f: vphi, fp: vphi_prime};
  
  return varphi;
} // end of Array.varphi


// -------------- Main fct -----------------------------------------------------


function OneSubjectEstimates(data, epsilon=1e-6, max_niter=250)
{
  // data: an object with one subjects's data, including censored data, in
  // the form {y, gt, lt, interval.gt, interval.lt}
  
  // If this fct does not converge, the user, depending on his needs,
  // may want to use first_guess.mu & first_guess.sigma as initial values
    
  var estimates = {converged: false, mu: [], sigma: []};
  var tmp = {mu: [], sigma: []};
  var y_sum;
 
  
  if (data.y.length > 0)
  {
    var ybar = data.y.mean(),
        ysd  = data.y.sd();
      
    tmp.mu.push(ybar);
    if (!isNaN(ysd)) tmp.sigma.push(ysd);
      
    y_sum = data.y.sum();
  }
  else
  {
    y_sum = 0;
  }
    
    
  if (data.interval.gt.length > 0)
  {
    var l = data.interval.gt.median(),
        u = data.interval.lt.median();
      
    var mu = (l + u) / 2;
    var sigma = (u - l) / 4;
      
    tmp.mu.push(mu);
    tmp.sigma.push(sigma);
  }
    
    
  if (tmp.mu.length > 0)
  {
    var mu = tmp.mu.median();
      
    if (data.gt.length > 0)
    {
      var l = data.gt.median();
      var sigma = Math.abs(mu - l) / 2;
      tmp.sigma.push(sigma);
    }
      
    if (data.lt.length > 0)
    {
      var u = data.lt.median();
      var sigma = Math.abs(u - mu) / 2;
      tmp.sigma.push(sigma);
    }       
  }
  else if (data.lt.length > 0 && data.gt.length > 0)
  {
    var l = data.gt.median(),
        u = data.lt.median();
      
    var mu = (l + u) / 2;
    var sigma = Math.abs(u - l) / 4;
    
    tmp.mu.push(mu);
    tmp.sigma.push(sigma);
  }
    
    
  var mu    = tmp.mu.median(),
      sigma = tmp.sigma.median();
      
  var first_guess = {mu: mu, sigma: sigma};


  // Run a univariate N-R to find the best sigma given (fixed) mu
    
    // But first make sure that we are on the right side of sigma-hat
      
    var aligned = {y: [], lt: [], gt: []};
      
    if (data.y.length  > 0) aligned.y  =  data.y.minus(mu);
    if (data.lt.length > 0) aligned.lt = data.lt.minus(mu);
    if (data.gt.length > 0) aligned.gt = data.gt.minus(mu).times(-1);
      
    var ay2sum = aligned.y.sqSum();
      
      
    cont = true;
      
    while (cont)
    {
      df2_dsigma = data.y.length / sigma**2;
      if (data.y.length  > 0) df2_dsigma -= 3*ay2sum / sigma**4;
      if (data.lt.length > 0) df2_dsigma -= d2logPhi_dsigma2(aligned.lt, sigma);
      if (data.gt.length > 0) df2_dsigma -= d2logPhi_dsigma2(aligned.gt, sigma);
      if (data.interval.gt.length > 0) df2_dsigma -= d2logPhiInt_dmudsigma(data.interval, mu, sigma);
    
      if (df2_dsigma > 0) sigma /= 2;
      else cont = false;
    }
    

  cont = true;
  iter = 0;
    
    
  while (cont)
  { 
    f2 = -data.y.length / sigma;
    if (data.y.length  > 0) f2 += ay2sum / sigma**3; 
    if (data.lt.length > 0) f2 -= dlogPhi_dsigma(aligned.lt, sigma)
    if (data.gt.length > 0) f2 -= dlogPhi_dsigma(aligned.gt, sigma)
    if (data.interval.gt.length > 0) f2 -= dlogPhiInt_dsigma(data.interval, mu, sigma);  
   
    df2_dsigma = data.y.length / sigma**2;
    if (data.y.length  > 0) df2_dsigma -= 3*ay2sum / sigma**4;
    if (data.lt.length > 0) df2_dsigma -= d2logPhi_dsigma2(aligned.lt, sigma);
    if (data.gt.length > 0) df2_dsigma -= d2logPhi_dsigma2(aligned.gt, sigma);
    if (data.interval.gt.length > 0) df2_dsigma -= d2logPhiInt_dmudsigma(data.interval, mu, sigma);
  
    if (df2_dsigma > 0)
    {
      sigma /= 2;
    }
    else
    {
      change = f2 / df2_dsigma;
      sigma -= change;
      if (sigma < 0) sigma = (sigma + change) / 2;

      iter++;
      converged = Math.abs(change) < epsilon;
      cont = !converged & iter < max_niter; 
    } 
  }  
        

  // Run multivariate Newton-Raphson algorithm
    
  var cont = true,
      converged = false,
      iter = 0,
      theta = [mu, sigma];
    
  var aligned, aligned_sum,
      f1, f2, 
      df1_dmu,
      df2_dmu, df2_dsigma; 

  var f, J;
  var max_change, previous_sigma;
            
    
  while (cont)
  {
    mu    = theta[0],
    sigma = theta[1];
     
    aligned = data.y.minus(mu); // array
    aligned_sum = y_sum - mu * data.y.length;
      
    f1 = aligned_sum / sigma**2;                           // scalar
    f2 = -data.y.length/sigma + aligned.sqSum()/sigma**3;  // scalar
      
      
    // Define Jacobian (first part)
      
    df1_dmu    =  - data.y.length / sigma**2;
    df2_dmu    = -2 * aligned_sum / sigma**3;
    df2_dsigma =    data.y.length / sigma**2;
      if (data.y.length > 0) df2_dsigma -= 3 * aligned.sqSum() / sigma**4;
      
       
    if (data.lt.length > 0)
    {
      let z = data.lt.minus(mu); // array
      f1 -= dlogPhi_dmu(z, sigma);
      f2 -= dlogPhi_dsigma(z, sigma);
        
      // Complete defn of Jacobian
        
      df1_dmu    -=      d2logPhi_dmu2(z, sigma);
      df2_dmu    -= d2logPhi_dmudsigma(z, sigma);
      df2_dsigma -=   d2logPhi_dsigma2(z, sigma);
    } 
      
      
    if (data.gt.length > 0)
    {
      let z = data.gt.minus(mu).times(-1);
      f1 += dlogPhi_dmu(z, sigma);
      f2 -= dlogPhi_dsigma(z, sigma);
        
      // Complete defn of Jacobian
        
      df1_dmu +=       d2logPhi_dmu2(z, sigma, 1);
      df2_dmu -=  d2logPhi_dmudsigma(z, sigma, 1);
      df2_dsigma -= d2logPhi_dsigma2(z, sigma);
    }
      
      
    if (data.interval.gt.length > 0)
    {
      f1 -=    dlogPhiInt_dmu(data.interval, mu, sigma);
      f2 -= dlogPhiInt_dsigma(data.interval, mu, sigma);
      
      // Complete defn of Jacobian
        
      df1_dmu    -=      d2logPhiInt_dmu2(data.interval, mu, sigma);
      df2_dmu    -= d2logPhiInt_dmudsigma(data.interval, mu, sigma);
      df2_dsigma -= d2logPhiInt_dmudsigma(data.interval, mu, sigma);
    }
          
      
    J = [[df1_dmu, df2_dmu], [df2_dmu, df2_dsigma]];
    f = [f1, f2];
      
    change = NewtonRaphsonChange(J, f); 
      
    if (change.any_NaN())
    {
      cont = false;
    }
    else
    {
      let max_change = change.max_abs();
      
      let previous_sigma = theta[1];
      theta = theta.minus(change);
      if (theta[1] < 0) theta[1] = previous_sigma / 2;

      iter++;
      converged = max_change < epsilon;
      cont = !converged & iter < max_niter;
    }
  }
    
  mu = theta[0];
  sigma = theta[1];
    
    
  var OneSubjectEstimates = {converged: converged, mu: mu, sigma: sigma, n_iter: iter, first_guess: first_guess};
    
  return OneSubjectEstimates;
} // end of OneSubjectEstimates