// Riskband/Banerjee algorithm
// Author -- Patrick B�lisle

// Version 0.2 (May 2021)
// [distributed]

// Requires S.js


// Change Log
// ======================
//
// Version 0.2 (7 May 2021)
// ------------------------
//  Changed the call to logPhiInterval [deprecated] to a call to PhiInterval
//
//  Changed the following fct calls, as they are not defined through Array.prototype anymore:
//    - any_duplicated_cutoff
//    - diff
//    - dot_product
//    - mean
//    - rnorm_gt, rnorm_lt, rnorm_interval
//    - sd
//    - sqSum
//    - times [renamed prod]
//
//   Dropped a call to .plus (renamed add(), but not used in this new version)
//
//   The properties a & b in literal 'candidates' were renamed gt & lt 
//   (making it directly usable as an argument in call to PhiInterval)
//
//   Slightly modified the function log_prior_density to cause an immediate error message when the values specified in A
//   do not all cross the domain mu x sigma. In the original version of this program, the model was run before the error
//   msg was displayed, even though it did not make sense (as one of the regions cannot have a non-null prior probability
//   if it does NOT cross the mu x sigma domain, hence being mathematically impossible/empty).
//
//   Moved the call to function PleaseBePatient() -- displaying a message informing that MCMC has started -- 
//   a little bit further down in the main fct.
//
//   Added a protection against negative values in data when assuming log-normal distribution of the data.
//
//   Added a call to desc_stats.

//
// Version 0.1 (Mar 2021)
// ----------------------
//   Original code

                                        
any_duplicated_cutoff = function(arr)
{
  var any_duplicate = false;
  
  for (let i=1; i<arr.length && !any_duplicate; i++)
  {
    if (arr[i] == arr[i-1]) any_duplicate = true;
  }
  
  return any_duplicate;
} // end of any_duplicated_cutoff


log_prior_density = function(A, Z, domain, prior_density, region_prior_prob)
{
  const mu_lower = domain.mu[0];
  const mu_upper = domain.mu[1];
  
  const sigma_lower = domain.sigma[0];
  const sigma_upper = domain.sigma[1];
  
  const R = A.length + 1;


  // Make sure that each region delimited by terms in 'A' 
  // does intersect with the domain (mu_lower, mu_upper) x (sigma_lower, sigma_upper)
  
  var intersects = {b: [], r: [], t: [], l: []}; // whether each cut-off in A cross the bottom(b), right(r), top(t) & left(l)
                                                 // borders of the rectangle domain of (mu, sigma)
 
  var tmp;                            
  tmp = A.map(a => a - Z*sigma_lower);
  intersects.b = tmp.map(mu => mu_lower < mu & mu < mu_upper);
  
  tmp = A.map(a => (a - mu_upper) / Z);
  intersects.r = tmp.map(sigma => sigma_lower < sigma & sigma < sigma_upper);
  
  tmp = A.map(a => a - Z*sigma_upper);
  intersects.t = tmp.map(mu => mu_lower < mu & mu < mu_upper);
  
  tmp = A.map(a => (a - mu_lower) / Z);
  intersects.l = tmp.map(sigma => sigma_lower < sigma & sigma < sigma_upper);
  
    // List cut-offs that do not cross the domain
    var A_problematic = [];
    
    for (let i=0; i<R; i++)
    {
      if (intersects.b[i] == 0 & intersects.r[i] == 0) A_problematic.push(A_orig[i]);
    }
    
    if (A_problematic.length > 0) return {error: true, A: A_problematic};
    
  
  var S = RegionAreas(intersects, A, domain, Z); // Region surface/areas
  var log_S = S.map(Math.log);
  

  // Compute prior density for each region
  
  if (prior_density.equal_region_probs) log_region_prior_prob = rep(0, R);
  else                                  log_region_prior_prob = region_prior_prob.map(Math.log);   
  
  
  var log_prior_dens = [0];
  
  for (let j=1; j<R; j++)
  {
    let log_prior_dens_j = log_region_prior_prob[j] + log_S[0] - log_region_prior_prob[0] - log_S[j];
    log_prior_dens.push(log_prior_dens_j);
  }


  return {prior: log_prior_dens, error: false};  // no need to standardize log_prior_dens
} // end of log_prior_density


muGen = function(y_bar, mu_sd, mu_range, sigma, Z, A, prior_density, log_prior_dens)
{
  var lower_limit,
      upper_limit;

  if (prior_density.uniform)
  {
    lower_limit = mu_range[0];
    upper_limit = mu_range[1];
  }
  else
  {
    mu_cutoffs = A.map(a => a - Z*sigma);
    candidates = sampling_j_and_intervals(mu_cutoffs, mu_range);

    // Sample a segment
  
    if (candidates.j.length == 1)
    { 
      lower_limit = mu_range[0];
      upper_limit = mu_range[1];
    }
    else
    { 
      let logp = PhiInterval(candidates, y_bar, mu_sd, true);
  
      // Add log(prior_density) to each segment
      candidates.j.forEach((j, i) => logp[i] += log_prior_dens[j]); 

      j = sample_discrete(seq(candidates.j.length), logp, true);
            
      lower_limit = candidates.gt[j];
      upper_limit = candidates.lt[j];
    }
  }


  // Sample a value within the sampling range (range of the above-selected segment if that is the case)
  
  mu = rnorm_interval(y_bar, mu_sd, lower_limit, upper_limit);
  
  return mu;
} // end of muGen


RegionAreas = function(intersects, A, domain, Z)
{
  var S = [];
  
  const mu_lower = domain.mu[0];
  const mu_upper = domain.mu[1];
  
  const sigma_lower = domain.sigma[0];
  const sigma_upper = domain.sigma[1];
  
  const mu_width    = mu_upper    - mu_lower;
  const sigma_width = sigma_upper - sigma_lower;
  const R = intersects.b.length + 1;
  
  
  for (let r=0; r<R; r++)
  {
    var area = 0;
    
    if (r == 0)
    {
      if (intersects.b[r] & intersects.l[r])
      {
        // region 72 in R code
        let H = (A[r]-mu_lower)/Z - sigma_lower;
        let B = Z * H;
        area = B * H / 2;
      }
      else if (intersects.b[r] & intersects.t[r])
      {
        // region 63 in R code
        let b =  domain.sigma.map(x => -Z*x + A[r] - mu_lower);
        let H = sigma_width;
        area = H * (b[0] + b[1]) / 2;
      }
      else if (intersects.r[r] & intersects.l[r])
      {
        // region 71 in R code
        let h = domain.mu.map(mu => (A[r] - mu)/Z - sigma_lower);
        let B = mu_width;
        area = B * (h[0] + h[1]) / 2;
      }
      else if (intersects.r[r] & intersects.t[r])
      {
        // region 62 in R code
        let B = mu_width;
        let H = sigma_width;
        let h = sigma_upper - (A[r] - mu_upper) / Z;
        let b = Z * h;
        area = B*H - b*h/2;
      }
    }
    else if (r == R-1)
    {
      if (intersects.b[r-1] & intersects.l[r-1])
      {
        // region 34 in R code
        let B = mu_width;
        let H = sigma_width;
        let h = (A[r-1] - mu_lower)/Z - sigma_lower;
        let b = h*Z;
        area = B*H - b*h/2; 
      }
      else if (intersects.b[r-1] & intersects.t[r-1])
      {
        // region 7 in R code
        let H = sigma_width;
        let B = mu_upper - A[r-1] + Z*sigma_upper;
        let b = B - Z*H;
        area = (B+b)/2 * H;
      }
      else if (intersects.r[r-1] & intersects.l[r-1])
      {
        // region 31 in R code
        let h = domain.mu.map(mu => sigma_upper - (A[r-1] - mu)/Z);
        let B = mu_width;
        area = B * (h[0] + h[1]) / 2;
      }
      else if (intersects.r[r-1] & intersects.t[r-1])
      {
        // region 4 in R code
        let H = sigma_upper - (A[r-1] - mu_upper)/Z;
        let B = Z*H;
        area = B*H/2;
      }
    }
    else if (intersects.b[r-1] & intersects.l[r-1] & intersects.b[r] & intersects.l[r])
    { 
      // region 45 in R code
      let h = [A[r-1], A[r]].map(a => (a - mu_lower) / Z - sigma_lower);
      let b = h.map(h => h * Z);
      area = diff(prod(b,h)) / 2;
    }
    else if (intersects.b[r-1] & intersects.l[r-1] & intersects.b[r] & intersects.t[r])
    {
      // region 36 in R code
      let b = A[r] - Z*sigma_upper - mu_lower;
      let h = sigma_upper - (A[r-1]-mu_lower) / Z;
      let H = sigma_width;
      let B = A[r] - A[r-1];
      area = B*H - h*(B-b)/2;
    }
    else if (intersects.b[r-1] & intersects.l[r-1] & intersects.r[r] & intersects.l[r])
    {
      // region 44 in R code
      let B = mu_width;
      let H = (A[r] - A[r-1]) / Z; 
      let b = mu_upper - A[r-1] + Z*sigma_lower;
      let h = b/Z;
      area = B*H - b*h/2;
    }
    else if (intersects.b[r-1] & intersects.l[r-1] & intersects.r[r] & intersects.t[r])
    {
      // region 35 in R code
      let H = sigma_width; 
      let B = mu_width;
      let sigma = [];
      sigma.push((A[r-1] - mu_lower)/Z);
      sigma.push((A[r]   - mu_upper)/Z);
      let h = [];
      h.push(sigma[0] - sigma_lower);
      h.push(sigma_upper - sigma[1]);
      let b = h.map(h => h*Z);
      area = B*H - dot_product(b, h) / 2;
    }
    else if (intersects.b[r-1] & intersects.t[r-1] & intersects.b[r] & intersects.t[r])
    {
      // region 9 in R code
      let H = sigma_width;
      let B = A[r] - A[r-1];
      area = B*H; 
    }
    else if (intersects.b[r-1] & intersects.t[r-1] & intersects.r[r] & intersects.t[r])
    {
      // region 8 in R code
      let H = sigma_width ;
      let B = A[r] - A[r-1];
      let b_prime = Z*H;
      let h = sigma_upper - (A[r]-mu_upper) / Z;
      let b = Z*h;
      area = H*(B+b) - (H*b_prime + b*h)/2;
    }
    else if (intersects.r[r-1] & intersects.l[r-1] & intersects.r[r] & intersects.l[r])
    {
      // region 41 in R code 
      let B = mu_width;
      let H = (A[r] - A[r-1]) / Z;
      area = B*H;
    }
    else if (intersects.r[r-1] & intersects.l[r-1] & intersects.r[r] & intersects.t[r])
    {
      // region 32 in R code
      let B = mu_width;
      let H = (A[r] - A[r-1]) / Z;
      let b = A[r] - Z*sigma_upper - mu_lower;
      let h = b/Z;
      area = B*H - b*h/2;
    }
    else if (intersects.r[r-1] & intersects.t[r-1] & intersects.r[r] & intersects.t[r])
    {
      // region 5 in R code
      let b = [A[r-1], A[r]].map(a => mu_upper - a + Z*sigma_upper)
      let h = b.map(b => b/Z);
      area = - diff(prod(b, h)) / 2;  
    }
    
        
    S.push(area);
  }
  
  return S;
} // end of RegionAreas

                          
RiskbandSigmaGen = function(n, beta, sigma_range, mu, Z, A, prior_density, log_prior_dens)
{
  var lower_limit = sigma_range[0],
      upper_limit = sigma_range[1],
      area;

  var o = new RiskbandSigmaGenObject(n, beta);
    o.M = o.logf(o.mode);

  
  
  if (!prior_density.uniform)
  {
    var sigma_cutoffs = A.map(a => (a - mu) / Z);
    var candidates = sampling_j_and_intervals(sigma_cutoffs, sigma_range);
    
    if (candidates.j.length > 1)
    {
      let log_prob = [];
      let r = [];
      let areas = [];
      
      for (let i=0; i<candidates.j.length; i++)
      {
        let area = o.fCum(candidates.gt[i], candidates.lt[i]).result;
        areas.push(area);
        log_prob.push(Math.log(area) + log_prior_dens[candidates.j[i]]);
        r.push(candidates.j[i]);
      }
  
      
      let selected_index = sample_discrete(seq(candidates.j.length), log_prob, true);
  
      area = areas[selected_index];
      lower_limit = candidates.gt[selected_index];
      upper_limit = candidates.lt[selected_index];
    }
  }
   
   
  if (typeof area == 'undefined') area = o.fCum(lower_limit, upper_limit).result; 
       
  
  var u = runif(1);
  var target = area * u;
  
    
  var cont = true,
      iter = 0,
      NR = {epsilon: 1e-8, max_niter: 100},
      sigma = (lower_limit + upper_limit) / 2; // N-R starting point
    
  
  while (cont) 
  {
    iter++;
    let previous_sigma = sigma;
    let F = o.fCum(lower_limit, sigma).result;
    change = (F - target) / o.f(sigma);
    sigma -= change;
    
    if (F > target) upper_limit = previous_sigma;
    else lower_limit = previous_sigma, target -= F;
    
    if      (sigma < lower_limit) sigma = (previous_sigma + lower_limit) / 2;
    else if (sigma > upper_limit) sigma = (previous_sigma + upper_limit) / 2;
    
    converged = Math.abs(change) < NR.epsilon;
    cont = !converged && (iter < NR.max_niter);
  }
  
  
  return sigma;
} // end of RiskbandSigmaGen


RiskbandSigmaGenObject = function(n, beta) 
{
  this.n = n;
  this.b = beta;
  
  this.mode = Math.sqrt(2*beta/n);
} // end of RiskbandSigmaGenObject


RiskbandSigmaGenObject.prototype = 
{
  f: function(sigma) 
  {
      var tmp = Math.exp(this.logf(sigma) - this.M);
      return tmp;
  },

  fCum: function(lowerBound, upperBound) 
  {
      var self = this;
      var F = function(xA) 
      {
          return xA.map(self.f, self);
      };

      return new zygotine.Num.NumericIntegration(F, lowerBound, upperBound).compute();
  },

  logf: function(sigma) 
  {
      var A = this;
      return -A.n * Math.log(sigma) - A.b / Math.pow(sigma, 2);
  }
} // end of RiskbandSigmaGenObject.prototype


sampling_j_and_intervals = function(cutoffs, range)
{
  // Return the indices j along with their corresponding interval
  //  crossing the allowed domain (specified through 'range')
  //
  // cutoffs: array with cut-off points
  // range: an array of length with [min,max] for the random variable

  // For example, if cutoffs = [a, b, c]
  // then the possible values for j are
  //   j = 0 with domain = (-Inf, a)
  //       1 with domain   (a, b)
  //       2 with domain   (b, c)
  //   and 3 with domain   (c, Inf)
  // all of which will be checked against the actual allowed 'range'
  
  var candidates = {j: [], gt: [], lt: []};
  var a = range[0],
      b;
      
  var j = 0;
  
  var cont  = true,
      cont2 = true;
  
      
  // find first j crossing the domain
  
  if (cutoffs[cutoffs.length-1] <= range[0])
  {
    candidates.j.push(cutoffs.length-1);
    candidates.gt.push(range[0]);
    candidates.lt.push(range[1]);
  
    cont2 = false;
  }
  else
  {
    while (cont)
    {
      b = cutoffs[j];
    
      if (b > range[0])
      {
        cont = false;
      
        candidates.j.push(j);
        candidates.gt.push(a); 
      
        if (b > range[1] )
        {
          b = range[1];
          cont2 = false;
        }
      
        candidates.lt.push(b);
        a = b;
      }
    
      j++;
    }
  }
  
  
  // keep going until the end of domain is covered
  
  while (cont2)
  {
    b = cutoffs[j];
    
    if (b > range[1] | typeof b == 'undefined')
    {
      b = range[1];
      cont2 = false;
    }
    
    candidates.j.push(j);
    candidates.gt.push(a);
    candidates.lt.push(b);
    a = b;
    
    if (j++ == cutoffs.length) cont2 = false;
  }
  
  return candidates;
} // end of sampling_j_and_interval


function SimulatedValues(data, mu, sigma)
{
  // Generating fake values for the censored data
  
  var unobserved = [];
  
  if (data.gt.length > 0)
  {
    let z = rnorm_gt(mu, sigma, data.gt);
    unobserved = unobserved.concat(z);
  }
  
  if (data.lt.length > 0)
  {
    let z = rnorm_lt(mu, sigma, data.lt);
    unobserved = unobserved.concat(z);  
  }
  
  if (data.interval.gt.length > 0)
  {
    let z = rnorm_interval(mu, sigma, data.interval.gt, data.interval.lt);
    unobserved = unobserved.concat(z);
  }
  
  var unobserved = {sum_yi: sum(unobserved), sum_yi2: sqSum(unobserved)};
  
  return unobserved;
} // end of SimulatedValues


// =============================================================================
// ==  Main function                                                          ==


run_Riskband = function(specificParams)
{
  var R = document.riskband_form.R.value;
  var logNormalDistrn = document.getElementById('logN').checked;
    
  console.clear();
  ClearErrorMsg();
  
  
    ////////////////////////////////////////////////////////////////////////////
    // Read Problem Description from Html Form
    
    // Read data
    
    var data = ReadData(document.riskband_form, specificParams);
  
    if (data.N == 0)
    {
      return { ok: false, errmsg: ErrorMsg(3) }
    }
    
    // Read region cut-offs values (A) & prior probabilities
    // from calling html page
    
    var A = Read_A_fromHtml(R),
        region_prior_prob = [];
    
    var mu_lower  = Number(document.riskband_form.mu_lower.value),
        mu_upper  = Number(document.riskband_form.mu_upper.value),
        gsd_lower = Number(document.riskband_form.gsd_lower.value),
        gsd_upper = Number(document.riskband_form.gsd_upper.value);
        
    var prior_density = {equal_region_probs: document.getElementById("rp_equalwts").checked, 
                         uniform:            document.getElementById("rp_unif").checked};

       
    if (!is_sorted(A))
    {
      return { ok: false, errmsg: ErrorMsg(4) }
    }
    else if (any_duplicated_cutoff(A))
    {
      return { ok: false, errmsg: ErrorMsg(5) }
    }
    

    if (!prior_density.equal_region_probs && !prior_density.uniform)
    {
      // user-defined region probs
      
      for (let i=0; i<R; i++)
      {
        let html_varname = "rpp" + i;
        let region_prior_prob_i = document.getElementById(html_varname).value;
        region_prior_prob.push(region_prior_prob_i);
      }
      
      region_prior_prob = region_prior_prob.map(Number); // convert to numeric values
      
      let tot_prob = sum(region_prior_prob);
      
      if (tot_prob != 1)
      {
        return { ok: false, errmsg: ErrorMsg(0) }
      }
      
      let any_negative_value = region_prior_prob.filter(p => p < 0).length > 0;
      
      if (any_negative_value)
      {
        return { ok: false, errmsg: ErrorMsg(1) }
      }
    }
    
    
    mcmc = MCMCParms(document.riskband_form);

    
    ////////////////////////////////////////////////////////////////////////////
    
      
    var gm_min = Math.exp(mu_lower),
        gm_max = Math.exp(mu_upper),
      
        sigma_lower = Math.log(gsd_lower),
        sigma_upper = Math.log(gsd_upper); 
    

  const Z = 1.644854; // 95th percentile of N(0,1)
  
  const domain = {mu:    [   mu_lower,    mu_upper],
                  sigma: [sigma_lower, sigma_upper]};

  
  A_orig = A.slice();
  
  if (logNormalDistrn)
  { 
    if (data.any_le0())
    { 
      return { ok: false, errmsg: ErrorMsg(6) }
    }
  
    data.TakeLog();
    A = A.map(Math.log); // A is on the same scale as data
  }

 
  var log_prior_dens;
  
  if (!prior_density.uniform)
  { 
    let tmp = log_prior_density(A, Z, domain, prior_density, region_prior_prob);
    
    if (tmp.error)
    { 
      return { ok: false, errmsg: ErrorMsg(2, tmp.A) }
    }
    
    log_prior_dens = tmp.prior;
  }
  // else: we do not need to calculate log_prior_dens (leave it undefined) 
    
  
  // MCMC sampling
  
  PleaseBePatient();
  
  
  var burnin = {mu: [], sigma: []};
  var sample = {mu: [], sigma: []};
  

  var sqrt_N = Math.sqrt(data.N);
  
  var sum_yi    = sum(data.y),
      sum_yi2   = sqSum(data.y),
      simulated = {sum_yi: 0, sum_yi2: 0};
  
  
  // Find inits for (mu, sigma)
  
  var mu, sigma;
  

                   
  // var tmp = OneSubjectEstimates(data, undefined, undefined, domain, alt_start);
  var tmp = OneSubjectEstimates(data, undefined, undefined, domain);
  console.log("Inits", tmp);

  
  
  if (tmp.converged)
  {
    mu = tmp.mu;
    sigma = tmp.sigma;
  }
  else if (!isNaN(tmp.best_guess.mu))
  {
    mu = tmp.best_guess.mu;
    sigma = tmp.sigma;
  }
  else if (data.N > 0)
  {
    // we are most likely in the presence of data that are all left-censored or all right-censored
    
    if (data.lt.length > 0)
    {
      last_chance = data.lt;
    }
    else if (data.gt.length > 0)
    {
      last_chance = data.gt;
    }
    
    mu = mean(last_chance);
    sigma = sd(last_chance);
    if (isNaN(sigma)) sigma = -Infinity; // possible if there was only one observation
  }
  

  // Make sure that initial values for (mu, sigma) are in the domain of possible values
  
  if      (mu < mu_lower) mu = mu_lower;
  else if (mu > mu_upper) mu = mu_upper;
  
  if      (sigma < sigma_lower) sigma = sigma_lower;
  else if (sigma > sigma_upper) sigma = sigma_upper;
  
  
  // Prepare for MCMC
  
  var cont = true;  
  
  for (let iter=0; iter < mcmc.niter + mcmc.burnin & cont; iter++)
  {
    // Generate unobserved values for censored data
    
    if (data.any_censored)
    {
      let tmp = SimulatedValues(data, mu, sigma);
      simulated.sum_yi  = tmp.sum_yi;
      simulated.sum_yi2 = tmp.sum_yi2;    
    }

  
    // --- Sample from f(mu | sigma) ---------------------------------------------
    
    y_bar = (sum_yi + simulated.sum_yi) / data.N;
    mu_sd = sigma / sqrt_N;
      
    mu = muGen(y_bar, mu_sd, domain.mu, sigma, Z, A, prior_density, log_prior_dens);
    

    // --- Sample from f (sigma | mu) --------------------------------------------
     
    beta = (sum_yi2 + simulated.sum_yi2 - 2*mu*(sum_yi + simulated.sum_yi) + data.N*(mu**2)) / 2;
      
    sigma = RiskbandSigmaGen(data.N, beta, domain.sigma, mu, Z, A, prior_density, log_prior_dens);
    
    if (isNaN(sigma))
    {
      console.log("Generated sigma = NA on iteration = %d", iter);
      cont = false;
    }
          
  
    // Save/monitor parameter values
    
    if (iter < mcmc.burnin)
    {
      if (mcmc.monitor_burnin)
      {
        burnin.mu.push(mu);
        burnin.sigma.push(sigma);
      }
    }
    else
    {
      sample.mu.push(mu);
      sample.sigma.push(sigma);      
    }
  }
  
  
  desc_stats(sample);
  

  WorkComplete(sample);
  console.log("Goodbye!");
  return { ok: true, burnin, sample }
} // end of run_Riskband