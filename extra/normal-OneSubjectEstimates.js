// Code to find initial values in the context of Normal (or log-Normal observations)
// in presence of censored observations
//
// Author: Patrick Bélisle

// Version 0.3 (May 2021)
// [distributed]


// Change log
// ======================
//
//
// Version 0.3 (May 2021)
// ----------------------
//   Changed the following fct calls, as they are not defined through Array.prototype anymore
//     - any_NaN
//     - max_abs
//     - mean
//     - sd
//     - sqSum
//
//   The Jacobian in the Newton-Raphson algorithm is now defined through MyMatrix
//
//   Changed the calls to (deprecated) dlogPhi_dsigma & d2logPhi_dsigma2 for calls to dklogPhi_dsigmak
//
//   Added functions 
//     OneSubjectEstimates_log_post, OneSubjectEstimates_fixed_mu & OneSubjectEstimates_fixed_sigma
//
//   Moved the N-R section itself to a separate function (NR1Subject)
//
//
// Version 0.2 (Apr 2021)
// ----------------------
//   - moved function varphi to file stats.js
//   - calls to log_diff -> now calling minus_log
//   - calls to log_mult -> now calling lo_mult
//
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code

// Requires: derivatives.js  


function NR1Subject(data, mu, sigma, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter)
{
  if (typeof sum_yi  == 'undefined') sum_yi  =   sum(data.y);
  if (typeof sum_yi2 == 'undefined') sum_yi2 = sqSum(data.y);

  // Run multivariate Newton-Raphson algorithm
    
  var cont = true,
      converged = false,
      iter = 0;
    
  var aligned_sum,
      dl_dmu, dl_dsigma, 
      d2l_dmu2,
      d2l_dmudsigma, d2l_dsigma2; 

  var f;
  var J = create_matrix(2, 2);
       
    
  while (cont)
  {     
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
      let dkl_dmuk    =    dklogPhi_dmuk(data.lt, mu, sigma);
      let dkl_dsigmak = dklogPhi_dsigmak(data.lt, mu, sigma);

      dl_dmu    += dkl_dmuk.order1;
      dl_dsigma += dkl_dsigmak.order1;
        
      // Complete defn of Jacobian
        
      d2l_dmu2       += dkl_dmuk.order2;
      d2l_dmudsigma  += d2logPhi_dmudsigma(data.lt, mu, sigma);
      d2l_dsigma2    += dkl_dsigmak.order2;
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
      let dkl_dmuk    =    dklogPhi_dmuk(data.gt, mu, sigma, false);
      let dkl_dsigmak = dklogPhi_dsigmak(data.gt, mu, sigma, false);
    
      dl_dmu    += dkl_dmuk.order1;
      dl_dsigma += dkl_dsigmak.order1;
        
      // Complete the definition of the Jacobian
        
      d2l_dmu2      += dkl_dmuk.order2;
      d2l_dmudsigma += d2logPhi_dmudsigma(data.gt, mu, sigma, false);
      d2l_dsigma2   += dkl_dsigmak.order2;
    }
          
      
    J.m = [[d2l_dmu2, d2l_dmudsigma], [d2l_dmudsigma, d2l_dsigma2]];
    f = [dl_dmu, dl_dsigma];
      
    change = NewtonRaphsonChange(J, f);
    
      
    if (any_NaN(change))
    {
      cont = false;
    }
    else
    {
      let previous_theta = [mu, sigma];
      
      mu    -= change[0];
      sigma -= change[1];
      
      let max_change = max_abs(change);
      converged = max_change < epsilon;
              
                   
        if (sigma < 0)
        {
          let previous_sigma = sigma + change[1];
          
          let tmp = OneSubjectEstimates_fixed_mu(data, mu, previous_sigma/2, sigma_prior, [0, Infinity]);
          if (tmp.converged) sigma = tmp.sigma;
          else sigma = previous_sigma / 2;
        }
          
        // We ignored the restricted domain of both mu & sigma
        //if      (mu < domain.mu[0]) mu = (domain.mu[0] + mu + change[0]) / 2;
        //else if (mu > domain.mu[1]) mu = (domain.mu[1] + mu + change[0]) / 2;
        //
        //if      (sigma < domain.sigma[0]) sigma = (domain.sigma[0] + sigma + change[1]) / 2;
        //else if (sigma > domain.sigma[1]) sigma = (domain.sigma[1] + sigma + change[1]) / 2;
      

      iter++;
      max_change = NewtonRaphson_max_change([mu, sigma], previous_theta);
      let stuck = max_change < epsilon;
      
      cont = !converged & !stuck & iter < max_niter;
    }
  }
  
  return {converged: converged, iter: iter, mu: mu, sigma: sigma};
} // end of NR1Subject


function OneSubjectEstimates(data, mu_prior, sigma_prior, domain, theta_start, loop_over_mu=true, epsilon=1e-5, max_niter=100)
{
  // data: an object with one subjects's data, including censored data, in
  // the form {y, gt, lt, interval.gt, interval.lt}
  
  // domain: a literal object with mu & sigma dimensions/properties, arrows of length 2
  
  // theta_start: an optional literal, with dimensions/proporties mu & sigma,
  // giving  starting values in N-R algorithm; when missing, the mode will be used
  // -- if the prior distribution (mu_prior or sigma_prior) does have a mode() function ---
  // or the midpoint of the domain.
  
  // If the first attempt to estimate Subject Estimates does not work, the algorithm
  // will loop over either mo or sigma (based on the value of loop_over_mu),
  // the other parameter taking the value maximizing log posterior given the (fixed) value 
  // of the other parameter.
    
  // If the fct OneSubjectEstimates does not converge, the user, depending on his needs,
  // may want to use best_guess.mu & best_guess.sigma as initial values
  
    
  // Default uniform prior for mu & sigma
  
  if (typeof    mu_prior == 'undefined')    mu_prior = new logf_unif();
  if (typeof sigma_prior == 'undefined') sigma_prior = new logf_unif();
  
  if (typeof domain == 'undefined') domain = {mu: [-Infinity, Infinity], sigma: [0, Infinity]};
  
  const Q = 100;
  
  const NR2_parms = {epsilon: 1e-4, max_niter: 50};
    

  var converged = false,
      cont = true,
      iter = 0,
      mu, sigma,
      sum_yi, sum_yi2,
      best_guess = {logpost: -Infinity, mu: undefined, sigma: undefined};
      
      
  if (typeof theta_start == 'undefined')
  {
    if (typeof mu_prior.mode == 'function') mu = mu_prior.mode();
    else mu = (domain.mu[0] + domain.mu[1]) / 2;
    
    if (typeof sigma_prior.mode == 'function') sigma = sigma_prior.mode();
    else sigma = (domain.sigma[0] + domain.sigma[1]) / 2;
    
    theta_start = {mu: mu, sigma: sigma};
  }
  else mu = theta_start.mu, sigma = theta_start.sigma;     
  

  if (data.y.length == 0) sum_yi = 0,           sum_yi2 = 0;
  else                    sum_yi = sum(data.y), sum_yi2 = sqSum(data.y);
  
  
  while (cont)
  {
    // find the 2nd parm value that maximizes log posterior GIVEN the value of the 1st parm
    if (loop_over_mu)
    {
      tmp1 = OneSubjectEstimates_fixed_mu(data, mu, sigma, sigma_prior, domain.sigma, NR2_parms.epsilon, NR2_parms.max_niter);  
      if (tmp1.converged) sigma = tmp1.sigma;
    }
    else
    {
      tmp1 = OneSubjectEstimates_fixed_sigma(data, mu, sigma, mu_prior, domain.mu, NR2_parms.epsilon, NR2_parms.max_niter);
      if (tmp1.converged) mu = tmp1.mu;
    }
       
    
    if (tmp1.converged)
    {    
      let tmp = NR1Subject(data, mu, sigma, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter);
      converged = tmp.converged;
      
      if (converged) mu = tmp.mu, sigma = tmp.sigma, cont = false;
      else
      {
        let logpost = OneSubjectEstimates_log_post(data, mu, sigma, mu_prior, sigma_prior); 
        if (logpost > best_guess.logpost) best_guess = {logpost: logpost, mu: mu, sigma: sigma};
      }
    }
    
    
    if (!converged)
    {
      // will try another approach
      let my_start = {mu: mu, sigma: sigma};
      
      let tmp = OneSubjectEstimates_latentY(data, mu_prior, sigma_prior, domain, my_start, sum_yi, sum_yi2, epsilon, max_niter);
    
      if (tmp.converged) converged = true, mu = tmp.mu, sigma = tmp.sigma;
    }
    
    
    if (!converged && iter == 0)
    {      
      let tmp = OneSubjectEstimates_latentY_dataDriven(data, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter);
    
      if (tmp.converged) converged = true, mu = tmp.mu, sigma = tmp.sigma;
    }
        
     
    cont = !converged && iter < Q;
    if (!cont) break; 
      

    p = (0.5 + iter++) / Q;
        
    if (loop_over_mu)
    {
      if (typeof mu_prior.quantile == 'function') mu = mu_prior.quantile(p);
      else mu = domain.mu[0] + p * (domain.mu[1] - domain.mu[0]);
    }
    else
    {
      // looping over sigma
      if (typeof sigma_prior.quantile == 'function') sigma = sigma_prior.quantile(p);
      else sigma = domain.sigma[0] + p * (domain.sigma[1] - domain.sigma[0]);
    }
  }
  
    
  if (!converged)
  {
    // will try another approach
    let tmp = OneSubjectEstimates_latentY_alt(data, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter);
    
    if (tmp.converged) converged = true, mu = tmp.mu, sigma = tmp.sigma;
  }
   

  return {converged: converged, mu: mu, sigma: sigma, best_guess: best_guess};
} // end of OneSubjectEstimates


function OneSubjectEstimates_latentY(data, mu_prior, sigma_prior, domain, theta_start, sum_yi, sum_yi2, epsilon, max_niter)
{
  var latent_y = data.y.slice();
  
  var mu    = theta_start.mu;
  var sigma = theta_start.sigma;
  

  for (let i=0; i<data.interval.gt.length; i++)
  {
    let y = rnorm_interval(mu, sigma, data.interval.gt[i], data.interval.lt[i]);
    latent_y.push(y);
  }

  
  for (let i=0; i<data.gt.length; i++)
  {
    let y = rnorm_gt(mu, sigma, data.gt[i]);
    latent_y.push(y);
  }

  
  for (let i=0; i<data.lt.length; i++)
  {
    let y = rnorm_lt(mu, sigma, data.lt[i]);
    latent_y.push(y);
  }

  
  
  mu = mean(latent_y);
  if (latent_y.length > 1) sigma = sd(latent_y);

  tmp = NR1Subject(data, mu, sigma, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter);
   
  return tmp;
} // end of OneSubjectEstimates_latentY


function OneSubjectEstimates_latentY_alt(data, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter)
{
  // This is an alternative to the above OneSubjectEstimates_latentY function.
  // If there is a good number of 'y' points, the algorithm should converge.

  var latent_y = data.y.slice();
  
  var mu = mean(latent_y);
  var sigma = sd(latent_y);
  

  for (let i=0; i<data.interval.gt.length; i++)
  {
    let y = rnorm_interval(mu, sigma, data.interval.gt[i], data.interval.lt[i]);
    latent_y.push(y);
  }
  
  mu = mean(latent_y);
  sigma = sd(latent_y);

  
  for (let i=0; i<data.gt.length; i++)
  {
    let y = rnorm_gt(mu, sigma, data.gt[i]);
    latent_y.push(y);
  }

  
  for (let i=0; i<data.lt.length; i++)
  {
    let y = rnorm_lt(mu, sigma, data.lt[i]);
    latent_y.push(y);
  }


  mu = mean(latent_y);
  sigma = sd(latent_y);

  tmp = NR1Subject(data, mu, sigma, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter);
   
  return tmp;
} // end of OneSubjectEstimates_latentY_alt


function OneSubjectEstimates_latentY_dataDriven(data, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter)
{           
  // Data-driven inits (to Inits)
  
  const SIGMA_MIN = 0.001;
  var latent_y = data.y.slice();
  var mu, sigma;
  

  if (latent_y.length > 0)
  {
    mu = mean(latent_y);
    sigma = sd(latent_y);
    
    if (sigma < SIGMA_MIN) sigma = SIGMA_MIN;
  }


  if (data.interval.gt.length > 0)
  {
    if (typeof mu == 'undefined')
    {
      for (let i=0; i<data.interval.gt.length; i++)
      {
        let midpoint = (data.interval.gt[i] + data.interval.lt[i]) / 2;
        latent_y.push(midpoint);
      } 
    }
    else
    {
      for (let i=0; i<data.interval.gt.length; i++)
      {  
        let endpoints = [data.interval.gt[i], data.interval.lt[i]];
        let my_sigma = max(endpoints.map(e => e - mu).map(Math.abs)) / 2;
        if (my_sigma < SIGMA_MIN) my_sigma = SIGMA_MIN;
        
        let y = rnorm_interval(mu, my_sigma, data.interval.gt[i], data.interval.lt[i]);
        latent_y.push(y);
      }
    }

    // Update (mu, sigma) estimates
    mu = mean(latent_y);
    sigma = sd(latent_y);
    
    if (sigma < SIGMA_MIN) sigma = SIGMA_MIN;
  }
  
  
  if (data.gt.length > 0 || data.lt.length > 0)
  {
    if (data.gt.length > 0)
    {
      let gt_mu;
      let gt_sigma;
    
      if (typeof mu == 'undefined')
      {
        gt_mu = mean(data.gt);
        gt_sigma = sd(data.gt);
        if (gt_sigma < SIGMA_MIN) gt_sigma = SIGMA_MIN;
      }
      else gt_mu = mu, gt_sigma = sigma;
      
    
      for (let i=0; i<data.gt.length; i++)
      {
        if (gt_mu < data.gt[i]) 
        {
          let d = data.gt[i] - gt_mu;
          let my_sigma = d / 2;
          if (my_sigma < gt_sigma) my_sigma = gt_sigma;
          if (my_sigma < SIGMA_MIN) my_sigma = SIGMA_MIN;
        
          y = rnorm_gt(gt_mu, my_sigma, data.gt[i]);
        }
        else if (typeof mu != 'undefined') y = rnorm_gt(gt_mu, gt_sigma, data.gt[i]);
        else y = gt_mu;

        latent_y.push(y);
      }
    }
    
    
    if (data.lt.length > 0)
    {
      let lt_mu;
      let lt_sigma;
    
      if (typeof mu == 'undefined')
      {
        lt_mu = mean(data.lt);
        lt_sigma = sd(data.lt);
        if (lt_sigma < SIGMA_MIN) lt_sigma = SIGMA_MIN;
      }
      else lt_mu = mu, lt_sigma = sigma;
      
    
      for (let i=0; i<data.lt.length; i++)
      {
        if (lt_mu > data.lt[i])
        { 
          let d = lt_mu - data.lt[i];
          let my_sigma = d / 2;
          if (my_sigma < lt_sigma) my_sigma = lt_sigma;
          if (my_sigma < SIGMA_MIN) my_sigma = SIGMA_MIN;
        
          y = rnorm_lt(lt_mu, my_sigma, data.lt[i]);
        }
        else if (typeof mu != 'undefined') y = rnorm_lt(lt_mu, lt_sigma, data.lt[i]);
        else y = lt_mu;
        
        latent_y.push(y);
      }
    }    
    
    // Update (mu, sigma) estimates
    mu = mean(latent_y);
    sigma = sd(latent_y);
    
    if (sigma < SIGMA_MIN) sigma = SIGMA_MIN;
  }
  
  
  tmp = OneSubjectEstimates_fixed_mu(data, mu, sigma, sigma_prior);
  if (tmp.converged) sigma = tmp.sigma;

  tmp = NR1Subject(data, mu, sigma, mu_prior, sigma_prior, domain, sum_yi, sum_yi2, epsilon, max_niter);
   
  return tmp;
} // end of OneSubjectEstimates_latentY_dataDriven


function OneSubjectEstimates_fixed_mu(data, mu, sigma, sigma_prior, sigma_range, epsilon=1e-4, max_niter=100)
{
  // Run a univariate N-R to find the best sigma given (fixed) mu
  
  var sum_yi, sum_yi2;
      
  if (data.y.length > 0)
  {
    sum_yi  = sum(data.y); 
    sum_yi2 = sqSum(data.y);
  }
  else sum_yi = 0, sum_yi2 = 0;


  var y_var = sum_yi2 - 2*mu*sum_yi + data.y.length*mu**2,
      cont = true,
      iter = 0;
      
  if (typeof sigma_range == 'undefined') lower = 0,              upper = Infinity;
  else                                   lower = sigma_range[0], upper = sigma_range[1];
    
    
  while (cont)
  { 
    dl_dsigma   = -data.y.length / sigma;
    d2l_dsigma2 =  data.y.length / sigma**2;
    
    
    if (data.y.length  > 0) 
    {
      dl_dsigma   +=   y_var / sigma**3; 
      d2l_dsigma2 -= 3*y_var / sigma**4;
    }
    
    
    if (data.gt.length > 0)
    {
      let tmp = dklogPhi_dsigmak(data.gt, mu, sigma, false);
      
      dl_dsigma   += tmp.order1;
      d2l_dsigma2 += tmp.order2
    }
    
      
    if (data.lt.length > 0)
    {
      let tmp = dklogPhi_dsigmak(data.lt, mu, sigma);
      
      dl_dsigma   += tmp.order1;
      d2l_dsigma2 += tmp.order2;
    }
    
    
    if (data.interval.gt.length > 0)
    { 
      dl_dsigma   +=  dlogPhiInterval_dsigma(data.interval, mu, sigma);
      d2l_dsigma2 += d2logPhiInterval_dsigma2(data.interval, mu, sigma);
    }
    
    dl_dsigma   += sigma_prior.logf_prime(sigma);
    d2l_dsigma2 += sigma_prior.logf_second(sigma);
    
    
    if (dl_dsigma > 0) lower = sigma;
    else upper = sigma;
    
    
    let previous_sigma = sigma;
    
    change = - dl_dsigma / Math.abs(d2l_dsigma2);
    sigma -= change;
    
    
    if      (sigma < lower) sigma = (lower + sigma + change) / 2;
    else if (sigma > upper) sigma = (upper + sigma + change) / 2;

    iter++;
    converged = Math.abs(change) < epsilon;
    let stuck = Math.abs(sigma - previous_sigma) < epsilon;
    
    cont = !converged & !stuck & iter < max_niter;
  }  

  return {converged: converged, sigma: sigma, 
          sum_yi: sum_yi, sum_yi2: sum_yi2};
} // end of OneSubjectEstimates_fixed_mu


function OneSubjectEstimates_fixed_sigma(data, mu, sigma, mu_prior, mu_range, epsilon=1e-4, max_niter=100)
{
  // Run a univariate N-R to find the best mu given (fixed) sigma
  
  var sum_yi,
      n = data.y.length;
      
  if (n > 0) sum_yi = sum(data.y); 
  else       sum_yi = 0;


  var cont = true,
      iter = 0,
      sigma2 = sigma**2;
      
  if (typeof mu_range == 'undefined') lower = -Infinity,   upper = Infinity;
  else                                lower = mu_range[0], upper = mu_range[1];
    
    
  while (cont)
  { 
    dl_dmu = (sum_yi - n * mu) / sigma2;
    d2l_dmu2 = -n / sigma2;
    
  
    if (data.gt.length > 0)
    {
      let tmp = dklogPhi_dmuk(data.gt, mu, sigma, false)
      
      dl_dmu   += tmp.order1;
      d2l_dmu2 += tmp.order2
    }
    
      
    if (data.lt.length > 0)
    {
      let tmp = dklogPhi_dmuk(data.lt, mu, sigma);
      
      dl_dmu   += tmp.order1;
      d2l_dmu2 += tmp.order2;
    }
    
    
    if (data.interval.gt.length > 0)
    { 
      dl_dmu   +=   dlogPhiInterval_dmu(data.interval, mu, sigma);
      d2l_dmu2 += d2logPhiInterval_dmu2(data.interval, mu, sigma);
    }
    
    dl_dmu   += mu_prior.logf_prime(mu);
    d2l_dmu2 += mu_prior.logf_second(mu);
    
    
    if (dl_dmu > 0) lower = mu;
    else upper = mu;
    
    
    change = - dl_dmu / Math.abs(d2l_dmu2);
    mu -= change;
    
    if      (mu < lower) mu = (lower + mu + change) / 2;
    else if (mu > upper) mu = (upper + mu + change) / 2;

    iter++;
    converged = Math.abs(change) < epsilon;
    cont = !converged & iter < max_niter; 
  }  

  return {converged: converged, mu: mu};
} // end of OneSubjectEstimates_fixed_sigma


function OneSubjectEstimates_log_post(data, mu, sigma, mu_prior, sigma_prior)
{
  var n = data.y.length;
  var log_post = 0;
  
  if (n > 0)
  {
    log_post -= n * Math.log(sigma);
    log_post -= sum(data.y.map(y => (y-mu)**2)) / 2 / sigma**2;
  }

  
  if (data.lt.length > 0) log_post += sum(Phi(data.lt, mu, sigma,  true, true));
  if (data.gt.length > 0) log_post += sum(Phi(data.gt, mu, sigma, false, true));
  
  if (data.interval.gt.length > 0) log_post += sum(PhiInterval(data.interval, mu, sigma, true));
  
  
  log_post += mu_prior.logf(mu);
  log_post += sigma_prior.logf(sigma);

  return log_post;
} // end of OneSubjectEstimates_log_post