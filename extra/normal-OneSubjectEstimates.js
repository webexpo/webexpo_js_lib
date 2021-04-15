// Code to find initial values in the context of Normal (or log-Normal observations)
// in presence of censored observations
//
// Author: Patrick Bélisle

// Version 0.2 (Apr 2021)  -- not distributed


// Change log
// ======================
//
// Version 0.2 (Apr 2021)
// ----------------------
//   - moved function varphi to file stats.js
//   - calls to log_diff -> now calling minus_log
//   - calls to log_mult -> now calling lo_mult
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code

// Requires: derivatives.js  


function OneSubjectEstimates(data, mu_prior, sigma_prior, epsilon=1e-6, max_niter=250)
{
  // data: an object with one subjects's data, including censored data, in
  // the form {y, gt, lt, interval.gt, interval.lt}
  
  // ICI decrire mu_priorf
  
  // If this fct does not converge, the user, depending on his needs,
  // may want to use first_guess.mu & first_guess.sigma as initial values
  
  
  // Default uniform prior for mu & sigma
  
  if (typeof    mu_prior == 'undefined')    mu_prior = new logf_unif();
  if (typeof sigma_prior == 'undefined') sigma_prior = new logf_unif();
  
    
  var estimates = {converged: false, mu: [], sigma: []};
  var tmp = {mu: [], sigma: []};
  var sum_yi, sum_yi2;
 
  
  if (data.y.length > 0)
  {
    var ybar = data.y.mean(),
        ysd  = data.y.sd();
      
    tmp.mu.push(ybar);
    if (!isNaN(ysd)) tmp.sigma.push(ysd);
      
    sum_yi  = data.y.sum();
    sum_yi2 = data.y.sqSum();
  }
  else
  {
    sum_yi  = 0;
    sum_yi2 = 0;
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
      
    var y_var = sum_yi2 - 2*mu*sum_yi + data.y.length*mu**2;
      
      
    cont = true;
      
    while (cont)
    {
      d2l_dsigma2 = data.y.length / sigma**2;
      if (data.y.length  > 0) d2l_dsigma2 -= 3*y_var / sigma**4;
      if (data.lt.length > 0) d2l_dsigma2 += data.lt.d2logPhi_dsigma2(mu, sigma);
      if (data.gt.length > 0) d2l_dsigma2 += data.gt.d2logPhi_dsigma2(mu, sigma, false);
      if (data.interval.gt.length > 0) d2l_dsigma2 += d2logPhiInterval_dmudsigma(data.interval, mu, sigma);
    
      if (d2l_dsigma2 > 0) sigma /= 2;
      else cont = false;
    }
    

  cont = true;
  iter = 0;
    
    
  while (cont)
  { 
    dl_dsigma = -data.y.length / sigma;
    if (data.y.length  > 0) dl_dsigma += y_var / sigma**3; 
    if (data.lt.length > 0) dl_dsigma += data.lt.dlogPhi_dsigma(mu, sigma);
    if (data.gt.length > 0) dl_dsigma += data.gt.dlogPhi_dsigma(mu, sigma, false);
    if (data.interval.gt.length > 0) dl_dsigma += dlogPhiInterval_dsigma(data.interval, mu, sigma);  
   
    d2l_dsigma2 = data.y.length / sigma**2;
    if (data.y.length  > 0) d2l_dsigma2 -= 3*y_var / sigma**4;
    if (data.lt.length > 0) d2l_dsigma2 += data.lt.d2logPhi_dsigma2(mu, sigma);
    if (data.gt.length > 0) d2l_dsigma2 += data.gt.d2logPhi_dsigma2(mu, sigma, false);
    if (data.interval.gt.length > 0) d2l_dsigma2 += d2logPhiInterval_dmudsigma(data.interval, mu, sigma);
  
    if (d2l_dsigma2 > 0)
    {
      sigma /= 2;
    }
    else
    {
      change = dl_dsigma / d2l_dsigma2;
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
    
  var aligned_sum,
      dl_dmu, dl_dsigma, 
      d2l_dmu2,
      d2l_dmudsigma, d2l_dsigma2; 

  var f, J;
  var max_change, previous_sigma;
            
    
  while (cont)
  {
    mu    = theta[0],
    sigma = theta[1];
     
    aligned_sum = sum_yi - mu * data.y.length;
    y_var = sum_yi2 - 2*mu*sum_yi + data.y.length*mu**2;
      
    dl_dmu    = aligned_sum / sigma**2                +    mu_prior.logf_prime(mu);      // scalar
    dl_dsigma = -data.y.length/sigma + y_var/sigma**3 + sigma_prior.logf_prime(sigma);   // scalar
    
          
    // Define Jacobian (first part)
      
    d2l_dmu2      =  - data.y.length / sigma**2 +    mu_prior.logf_second(mu);
    d2l_dmudsigma = -2 * aligned_sum / sigma**3;
    d2l_dsigma2   =    data.y.length / sigma**2 + sigma_prior.logf_second(sigma);
      if (data.y.length > 0) d2l_dsigma2 -= 3 * y_var / sigma**4;
      
       
    if (data.lt.length > 0)
    {
      dl_dmu    +=    data.lt.dlogPhi_dmu(mu, sigma);
      dl_dsigma += data.lt.dlogPhi_dsigma(mu, sigma);
        
      // Complete defn of Jacobian
        
      d2l_dmu2       +=      data.lt.d2logPhi_dmu2(mu, sigma);
      d2l_dmudsigma  += data.lt.d2logPhi_dmudsigma(mu, sigma);
      d2l_dsigma2    +=   data.lt.d2logPhi_dsigma2(mu, sigma);
    } 
      
  
    if (data.interval.gt.length > 0)
    {
      dl_dmu    +=    dlogPhiInterval_dmu(data.interval, mu, sigma);
      dl_dsigma += dlogPhiInterval_dsigma(data.interval, mu, sigma);
      
      // Complete defn of Jacobian
        
      d2l_dmu2      +=      d2logPhiInterval_dmu2(data.interval, mu, sigma);
      d2l_dmudsigma += d2logPhiInterval_dmudsigma(data.interval, mu, sigma);
      d2l_dsigma2   +=   d2logPhiInterval_dsigma2(data.interval, mu, sigma);
    }
    
    
    if (data.gt.length > 0)
    {
      dl_dmu    +=    data.gt.dlogPhi_dmu(mu, sigma, false);
      dl_dsigma += data.gt.dlogPhi_dsigma(mu, sigma, false);
        
      // Complete the definition of the Jacobian
        
      d2l_dmu2      +=      data.gt.d2logPhi_dmu2(mu, sigma, false);
      d2l_dmudsigma += data.gt.d2logPhi_dmudsigma(mu, sigma, false);
      d2l_dsigma2   +=   data.gt.d2logPhi_dsigma2(mu, sigma, false);
    }
          
      
    J = [[d2l_dmu2, d2l_dmudsigma], [d2l_dmudsigma, d2l_dsigma2]];
    f = [dl_dmu, dl_dsigma];
      
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