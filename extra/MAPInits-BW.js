// =============================================================================
// Code written to calculate Initial Values for the Between-Worker model
// using the MAP (Maximum Posterior Estimate) as initial values 
// [based on inits-BetweenWorker.R version 0.25]                                                             
//
// Author: Patrick Bélisle                                                
//                              
// Version 0.4 (May 2021)
// [distributed]


// Change Log
// =============================================================================
//
//
// Version 0.4 (May 2021)
// ----------------------
//   The following fct calls were modified (as they are not defined through Array.prototype anymore):
//     - allTheSame [renamed all_the_same]
//     - divided_by [renamed ratio]
//     - mean
//     - sqSum
//
//   Using built-in .some() rather than home-made .any_* in a few occurences.
//
//   Locally-defined functions dlogPhi_dsigma & d2logPhi_dsigma2 were dropped
//   in favor of dklogPhi_dsigmak defined in derivatives.js (with different parametrization than that of previous local versions)
//
//   The Jacobian in the Newton-Raphson algorithm is now defined through MyMatrix
//
//   Slightly modified MAP_Inits to give sensible results when only one of the gt/lt endpoints is available
//     (in addition to at least one 'y' data point)
//
//   Changes (& correction!) in the function dl_dsigmaw (renamed dl_dsigma2_fct)
//
//   Modified logN fct (added f property)
//
//   Added mu_prior & sw_prior as arguments to OneSubjectEstimates
//
//   Added the function ThisWorkerData (easying the rerun of individual parm estimates for subjects
//    for which OneSubjectEstimates did not converge on first pass -- hoping for convergence the second time,
//     using parms obtained from other workers/subjects)
//
//
// Version 0.3 (Mar 2021)
// ----------------------
//   The name of a few functions called (from derivatives.js)
//   were changed (for internal/library consistency)
//
//
// Version 0.2 (Feb 2021)
// ----------------------
//   Now use the library lib/normal-OnesubjectEstimates rather than 
//   the previously embedded fct OneWorkerEstimates

         
dl_dsigmaw_fct = function(sigma_within, mu_worker, mu_overall, data, id, parms)
{
  var f  = - id.y.length / sigma_within,
      fp =   id.y.length / sigma_within**2; 
  
  
  if (id.y.length > 0)
  {
    let tmp = sqSum(data.y.map((z, i) => z - mu_worker[id.y[i]] - mu_overall));
    
    f  +=     tmp / sigma_within**3;
    fp -= 3 * tmp / sigma_within**4;
  }
  
  if (id.lt.length > 0)
  {
    let mu_j = elements(mu_worker, id.lt).map(m => m + mu_overall);
    let d = dklogPhi_dsigmak(data.lt, mu_j, sigma_within);
    
    f  += d.order1;
    fp += d.order2;
  }  
  
  if (id.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.gt).map(m => m + mu_overall);
    let d = dklogPhi_dsigmak(data.gt, mu_j, sigma_within, false);
    
    f  += d.order1;
    fp += d.order2;
  }
  
  if (id.interval.length > 0)
  {
    let mu_j = elements(mu_worker, id.interval).map(m => m + mu_overall);
    
    f  +=   dlogPhiInterval_dsigma(data.interval, mu_j, sigma_within);
    fp += d2logPhiInterval_dsigma2(data.interval, mu_j, sigma_within);
  }


  if (!parms.uupOnSds)
  {
    let logn = logN(sigma_within, parms);
    
    f  += logn.fp;
    fp += logn.fs;
  }
  
  
  return {f: f, fp: fp};
} // end of dl_dsigmaw_fct


function dl_dtheta_fct(theta, i, data, id, parms, epsilon=1e-6)
{
  // This function was written for validation ends only
  // -- it approximates the first derivatives of log posterior
  
  let f1 = log_posterior(theta, data, id, parms);
  theta[i] += epsilon;
  let f2 = log_posterior(theta, data, id, parms);
  let fp = (f2 - f1) / epsilon;
  
  return fp;
} // end of dl_dtheta_fct


function d2l_dtheta1dtheta2(theta, i, j, data, id, parms, epsilon=1e-5)
{
  let f1a = log_posterior(theta, data, id, parms);
  theta[i] += epsilon;
  let f1b = log_posterior(theta, data, id, parms);
  let fp1 = (f1b - f1a) / epsilon;
  
  theta[i] -= epsilon;
  theta[j] += epsilon;
  
  let f2a = log_posterior(theta, data, id, parms);
  theta[i] += epsilon;
  let f2b = log_posterior(theta, data, id, parms);
  let fp2 = (f2b - f2a) / epsilon;
  
  let fs = (fp2 - fp1) / epsilon;
  return fs;
} // end of d2l_dtheta1dtheta2


function d2l_dtheta2(theta, i, data, id, parms, epsilon=1e-5)
{
  // This function was written for validation ends only
  // -- it approximates the first derivatives of log posterior
  
  let f1 = log_posterior(theta, data, id, parms);
  theta[i] += epsilon;
  let f2 = log_posterior(theta, data, id, parms);
  theta[i] += epsilon;
  let f3 = log_posterior(theta, data, id, parms);
  
  let fp1 = (f2 - f1) / epsilon;
  let fp2 = (f3 - f2) / epsilon;
  
  let fs = (fp2 - fp1) / epsilon;
  
  return fs;
} // end of d2l_dtheta2


function f12(mu_worker, mu_overall, sigma_within, sigma_between, count, y_sum, data, id)
{   
  var tmp = y_sum.map((s, i) => s - count.y[i] * (mu_worker[i] + mu_overall));  // array of length #workers
      
      
  var dl_dmui = tmp.map(m => m / sigma_within**2);              // arrays of length #workers
  var d2l_dmui2 = count.y.map(u => -u /sigma_within**2); 
  var d2l_dmuidsigmaw = tmp.map(u => -2*u / sigma_within**3);
  
  
  if (data.lt.length > 0)
  {
    let mu_j = elements(mu_worker, id.lt).map(m => m + mu_overall);
    
    let tmp   =      dklogPhi_dmuk(data.lt, mu_j, sigma_within, true, false)
    let mixed = d2logPhi_dmudsigma(data.lt, mu_j, sigma_within, true, false);
    
    for (let i=0; i<data.lt.length; i++)
    {
      let j = id.lt[i];
      
      dl_dmui[j]         += tmp.order1[i];
      d2l_dmui2[j]       += tmp.order2[i];
      d2l_dmuidsigmaw[j] += mixed[i];
    }
  }
  
  
  if (data.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.gt).map(m => m + mu_overall);
    
    let tmp   =      dklogPhi_dmuk(data.gt, mu_j, sigma_within, false, false);
    let mixed = d2logPhi_dmudsigma(data.gt, mu_j, sigma_within, false, false);
    
    for (let i=0; i<data.gt.length; i++)
    {
      let j = id.gt[i];
      
      dl_dmui[j]         += tmp.order1[i];
      d2l_dmui2[j]       += tmp.order2[i];
      d2l_dmuidsigmaw[j] += mixed[i];
    }
  }
  
  
  if (data.interval.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.interval).map(m => m + mu_overall);
    
    let d1 = dlogPhiInterval_dmu(data.interval, mu_j, sigma_within, false);
    let d2 = d2logPhiInterval_dmu2(data.interval, mu_j, sigma_within, false);
    let mixed = d2logPhiInterval_dmudsigma(data.interval, mu_j, sigma_within, false);
       
    
    for (let i=0; i<data.interval.gt.length; i++)
    {
      let j = id.interval[i];
      
      dl_dmui[j]         += d1[i];
      d2l_dmui2[j]       += d2[i];
      d2l_dmuidsigmaw[j] += mixed[i]; 
    }     
  }
  
  
  var d2l_dmuidmu    = d2l_dmui2.slice();
  var dl_dmu         = sum(dl_dmui);
  var d2l_dmu2       = sum(d2l_dmui2);
  var d2l_dmudsigmaw = sum(d2l_dmuidsigmaw);
        
  dl_dmui = substract(dl_dmui, mu_worker.map(m => m/sigma_between**2));
  d2l_dmui2 = d2l_dmui2.map(d => d - 1/sigma_between**2); 
  
  
  var f12 = {dl_dmui: dl_dmui,
             dl_dmu: dl_dmu,
              
             d2l_dmui2: d2l_dmui2, 
             d2l_dmuidmu: d2l_dmuidmu,
             d2l_dmuidsigmaw: d2l_dmuidsigmaw,
              
             d2l_dmu2: d2l_dmu2,
             d2l_dmudsigmaw: d2l_dmudsigmaw};
  
  return f12;
} // end of f12


function logN(sigma, m)
{
  var log_sigma = Math.log(sigma);
  
  var f  =  - log_sigma - m.logSigmaWithinPrec/2 * (log_sigma - m.logSigmaWithinMu)**2;
  var fp = -1/sigma * (1 + m.logSigmaWithinPrec  * (log_sigma - m.logSigmaWithinMu));
  var fs = (1 - m.logSigmaWithinPrec * (m.logSigmaWithinMu + 1 - log_sigma)) / sigma**2;
  
  var logN = {f: f, fp: fp, fs: fs};
  
  return logN;
} // end of logN


function log_posterior(theta, data, id, parms)
{
  // function used only to validate intermediate results
  
  // theta = [mu_worker, mu_overall, sw, sb]
  
  let my_theta = theta.slice();
  
  let sb = my_theta.pop();
  let sw = my_theta.pop();
  let mu_overall = my_theta.pop();
  let mu_worker  = my_theta.slice();
  let n_workers = mu_worker.length;  
  
  var logf = 0;
  
  
 if (data.y.length > 0)
  {
    logf -= data.y.length * Math.log(sw);
    
    let mu_j = elements(mu_worker, id.y).map(m => m + mu_overall);
    
    let beta = 0;
    for (let i=0; i<data.y.length; i++) beta += (data.y[i] - mu_j[i])**2;
    
    logf -= beta / 2 / sw**2;
  }
  
  
  if (data.lt.length > 0)
  {
    let mu_j = elements(mu_worker, id.lt).map(m => m + mu_overall);
    logf += sum(Phi(data.lt, mu_j, sw, true, true)) 
  }
  
  
  if (data.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.gt).map(m => m + mu_overall);
    logf += sum(Phi(data.gt, mu_j, sw, false, true));
  }
  
  
  if (data.interval.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.interval).map(m => m + mu_overall);
    logf += sum(PhiInterval(data.interval, mu_j, sw, true));
  }  
  
  // prior on mu_worker
  
  logf += - 1/2 * sum(mu_worker.map(m => m**2))/sb**2 - n_workers * Math.log(sb);
  
  // prior on sigma_within & sigma_between

  if (!parms.uupOnSds)
  {
    let log_sw = Math.log(sw);
    logf += - log_sw - parms.logSigmaWithinPrec/2 * (log_sw - parms.logSigmaWithinMu)**2;
    
    let log_sb = Math.log(sb);
    logf += - log_sb - parms.logSigmaBetweenPrec/2 * (log_sb - parms.logSigmaBetweenMu)**2;
  }  
  
  
  return logf;
} // end of log_posterior


function optimal_sigma_between(mu_worker, parms, max_niter=100, epsilon=1e-5)
{
  var n = mu_worker.length;
  var beta = sum(mu_worker.map(m => m**2));
  var sb = Math.sqrt(beta/n); // starting value
  
  var cont = true,
      iter = 0;

  
  while (cont)
  {
    log_sb = Math.log(sb);
  
    fp =     beta/sb**3 - 1/sb * (n + 1 + parms.logSigmaBetweenPrec * (log_sb - parms.logSigmaBetweenMu));
    fs = - 3*beta/sb**4 + (n + 1 - parms.logSigmaBetweenPrec * (parms.logSigmaBetweenMu + 1 - log_sb)) / sb**2;

    change = fp / fs;
    sb -= change;

    
    if (sb < 0) sb = (sb + change) / 2;
  
    let converged = Math.abs(change) < epsilon;
    cont = !converged && ++iter < max_niter;    
  }
  
  if (!converged) sb =  Math.sqrt(beta/n);
  
  return {converged: converged, sigma_between: sb};
} // end of optimal_sigma_between


function optimal_sigma_within(sigma_within, mu_worker, mu_overall, data, id, parms, max_niter=100, epsilon=1e-5)
{
  var cont = true,
      iter = 0;
  
  
  while (cont)
  {
    let previous_sw = sigma_within;
    let tmp = dl_dsigmaw_fct(sigma_within, mu_worker, mu_overall, data, id, parms);
    
    var fp = tmp.f;
    var fs = tmp.fp;
    
    var change = - fp / Math.abs(fs);
    sigma_within -= change;
    
    if (parms.uupOnSds)
    {
      if      (sigma_within > parms.sigmaWithinRange[1]) sigma_within = (parms.sigmaWithinRange[1] + sigma_within + change) / 2;
      else if (sigma_within < parms.sigmaWithinRange[0]) sigma_within = (parms.sigmaWithinRange[0] + sigma_within + change) / 2;
    }
    else if (sigma_within < 0) sigma_within = previous_sw / 2;
    
    
    converged = Math.abs(sigma_within - previous_sw) < epsilon;
    cont = !converged && ++iter < max_niter;
  }


  return {converged: converged, sigma_within: sigma_within};
} // end of optimal_sigma_within


function sigma_within_RespectingRange(sigma_within, sp, previous_sigma_within)
{
  var out_of_range = false;
  
  if (sigma_within < sp.sigmaWithinRange[0])
  {
    if (typeof previous_sigma_within == 'undefined')
    {
      sigma_within = sp.sigmaWithinRange[0];
    }
    else
    {
      sigma_within = (previous_sigma_within + sp.sigmaWithinRange[0]) / 2;
    }
    
    out_of_range = true;
  }
  else if (sigma_within > sp.sigmaWithinRange[1])
  {
    if (typeof previous_sigma_within == 'undefined')
    {
      sigma_within = sp.sigmaWithinRange[1];
    }
    else
    {
      sigma_within = (previous_sigma_within + sp.sigmaWithinRange[1]) / 2;
    }
    
    out_of_range = true;
  }
  
  var sigma_within_RespectingRange = {sigma_within: sigma_within, out_of_range: out_of_range};
  
  return sigma_within_RespectingRange;
}  // end of sigma_within_RespectingRange


function sigma_within_fullcond_log(sw, mu_worker, mu_overall, data, id, sp)
{
  // this fct was written only for investigation
  var logf = 0;
  
  if (data.y.length > 0)
  {
    logf -= data.y.length * Math.log(sw);
    
    let mu_j = elements(mu_worker, id.y).map(m => m + mu_overall);
    
    let beta = 0;
    for (let i=0; i<data.y.length; i++) beta += (data.y[i] - mu_j[i])**2;
    
    logf -= beta / 2 / sw**2;
  }
  
  
  if (data.lt.length > 0)
  {
    let mu_j = elements(mu_worker, id.lt).map(m => m + mu_overall);
    logf += sum(Phi(data.lt, mu_j, sw, true, true)) 
  }
  
  
  if (data.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.gt).map(m => m + mu_overall);
    logf += sum(Phi(data.gt, mu_j, sw, false, true));
  }
  
  
  if (data.interval.gt.length > 0)
  {
    let mu_j = elements(mu_worker, id.interval).map(m => m + mu_overall);
    logf += sum(PhiInterval(data.interval, mu_j, sw, true));
  }
  
  
  // Add prior
  logf += logN(sw, sp).f;
  
  return round(logf, 2);
} // end of sigma_within_fullcond_log


function ThisWorkerData(obj)
{
  var my_data = {y: [], lt: [], gt: [], interval: {gt: [], lt: []}};
    
  my_data.y  = obj.filter(m => m.type === "uncensored").map(m => m.a);
  my_data.lt = obj.filter(m => m.type === "lessThan").map(m => m.a);
  my_data.gt = obj.filter(m => m.type === "greaterThan").map(m => m.a);
    
  my_data.interval.obj = obj.filter(m => m.type === "interval").map(m => m.a);
  my_data.interval.obj = obj.filter(m => m.type === "interval").map(m => m.b);
    
  return my_data;
} // end of ThisWorkerData


////////////////////////////////////////////////////////////////////////////////


function MAP_Inits(wd, sp, epsilon=1e-6) 
{
  // Find the Maximum A Posteriori (MAP) Estimator 
  //   and use it as initial values for the Between-Workers algorithm
  
  // console.clear();

  
  var workers = wd.workerIds;       
  
  var mu_prior = new logf_unif(),
      sw_prior;
  
  var domain = {mu: [sp.muOverallLower, sp.muOverallUpper], sigma: [0, Infinity]}; 
  var theta_start = {mu: (sp.muOverallLower + sp.muOverallUpper) / 2};
      
         
  if (sp.uupOnSds)
  {    
    domain.sigma = [sp.sigmaWithinRange[0], sp.sigmaWithinRange[1]];
    theta_start.sigma = (sp.sigmaWithinRange[0] + sp.sigmaWithinRange[1]) / 2;
    
    sw_prior = new logf_unif();
  }
  else
  {
    sw_prior = new logf_lnorm(sp.logSigmaWithinMu, 1/Math.sqrt(sp.logSigmaWithinPrec)); 
    theta_start.sigma = sw_prior.mode();
  }


  var count = {y : [], lt :[], gt: [], interval: []};
  var y_sum = [];
  var individual_estimates = {mu : [], sw : [], worker: [], worker_index : []};
  var availableIndividualEstimates = [];
  
  var data = {y: [], lt: [], gt: [], interval: {gt: [], lt: []}};
  var   id = {y: [], lt: [], gt: [], interval: []};


  for (let i = 0; i < workers.length; i++) 
  {
    var worker = workers[i];
    var wk = wd._measureByWorker[worker];    
    var this_worker_data = ThisWorkerData(wk);
    
     
    count.y.push(this_worker_data.y.length);
    y_sum.push(sum(this_worker_data.y));
            
    count.gt.push(this_worker_data.gt.length); 
    count.lt.push(this_worker_data.lt.length);
                            
    count.interval.push(this_worker_data.interval.gt.length);
    
      // Append to combined data
      
      data.y  = data.y.concat(this_worker_data.y);
      data.lt = data.lt.concat(this_worker_data.lt);
      data.gt = data.gt.concat(this_worker_data.gt);
      data.interval.gt = data.interval.gt.concat(this_worker_data.interval.gt);
      data.interval.lt = data.interval.lt.concat(this_worker_data.interval.lt);
      
      id.y  = id.y.concat(rep(i, this_worker_data.y.length));
      id.lt = id.lt.concat(rep(i, this_worker_data.lt.length));
      id.gt = id.gt.concat(rep(i, this_worker_data.gt.length));
      id.interval = id.interval.concat(rep(i, this_worker_data.interval.gt.length));
      
                  
    var estimable_parms = this_worker_data.y.length > 1 || this_worker_data.interval.gt.length > 0 || (this_worker_data.gt.length > 0 && this_worker_data.lt.length > 0); // for this worker 
    estimable_parms = true; // we want to *at least* try to get estimates for each worker       
            
    if (estimable_parms)
    {
      var run_1WorkerEstimate = true;
      var thisWorkerEstimates = {converged: false, mu: [], sigma: []};
      var ybar;
      var concatenated_l, 
          concatenated_u;
              
              
      if (this_worker_data.y.length <= 1)
      {
        run_1WorkerEstimate = false;
      }
      else
      {
        run_1WorkerEstimate = !all_the_same(this_worker_data.y);
      }
      
      run_1WorkerEstimate = true; // force estimation for each worker
              
    
      if (!run_1WorkerEstimate && this_worker_data.y.length > 0)
      {
        ybar = mean(this_worker_data.y);
        
        concatenated_l = this_worker_data.gt.concat(this_worker_data.interval.gt);
        concatenated_u = this_worker_data.lt.concat(this_worker_data.interval.lt);

        run_1WorkerEstimate = concatenated_l.some(e => e > ybar);
        if (!run_1WorkerEstimate) run_1WorkerEstimate = concatenated_u.some(e => e < ybar);
      }    
              
              
      if (!run_1WorkerEstimate)
      {
        thisWorkerEstimates.converged = false;
        thisWorkerEstimates.mu = NaN;
        
        let l_median;
        let u_median;
        let n_criteria = 0;
        
        if (concatenated_l.length > 0)
        {
          n_criteria++;
          l_median = median(concatenated_l);
        }
        
        if (concatenated_u.length > 0)
        {
          n_criteria++;
          u_median = median(concatenated_u);
        }
        
        if (this_worker_data.y.length > 0)
        {
          n_criteria++;
          thisWorkerEstimates.mu = ybar;
          
          if (n_criteria == 1 && !all_the_same(this_worker_data.y))
          {
            thisWorkerEstimates.converged = true;
            thisWorkerEstimates.sigma = sd(this_worker_data.y); 
          }
        }
        
        
        if (n_criteria > 1)
        {
          thisWorkerEstimates.converged = true;
          
          if (isNaN(thisWorkerEstimates.mu)) thisWorkerEstimates.mu = (l_median + u_median) / 2;
          
          if (concatenated_l.length > 0)
          {
            if (concatenated_u.length > 0)  thisWorkerEstimates.sigma = Math.abs(u_median - l_median) / 6;
            else                            thisWorkerEstimates.sigma = Math.abs(ybar     - l_median) / 3; 
          }
          else
          {
            thisWorkerEstimates.sigma = Math.abs(u_median - ybar) / 3; 
          }
        }  
      }
      else
      {
        if (individual_estimates.mu.length > 0)
        {
          theta_start.mu    = median(individual_estimates.mu);
          theta_start.sigma = median(individual_estimates.sw);
        }
        

        let tmp = OneSubjectEstimates(this_worker_data, mu_prior, sw_prior, domain, theta_start, false);
        
        thisWorkerEstimates.converged = tmp.converged;
        thisWorkerEstimates.mu        = tmp.mu;
        thisWorkerEstimates.sigma     = tmp.sigma;
      }
      
      
      availableIndividualEstimates.push(thisWorkerEstimates.converged);
      
      if (thisWorkerEstimates.converged)
      {
        individual_estimates.mu.push(thisWorkerEstimates.mu);
        individual_estimates.sw.push(thisWorkerEstimates.sigma);
        individual_estimates.worker_index.push(i);
        individual_estimates.worker.push(worker);
      }
    }
    else
    {
      availableIndividualEstimates.push(false);
    }          
  } 
  
  
  // Run OneSubjectEstimates again for workers for which first pass did not converge
  // (in case the 'educated' starting values may help convergence in 2nd pass)
  
  any_converged = false;
  all_converged = true;
  
  for (let i=0; i<availableIndividualEstimates.length; i++)
  {
    if (availableIndividualEstimates[i]) any_converged = true;
    else all_converged = false; 
  }
  
  
  if (!all_converged && any_converged)
  {
    let tmp_mu = individual_estimates.mu.slice();
    let tmp_sw = individual_estimates.sw.slice();
  
    for (let i=0; i<availableIndividualEstimates.length; i++)
    {
      if (!availableIndividualEstimates[i])
      {
        var worker = workers[i];
        var wk = wd._measureByWorker[worker];
        var this_worker_data = ThisWorkerData(wk);
        
        let tmp = OneSubjectEstimates(this_worker_data, mu_prior, sw_prior, domain, theta_start, false);
        
        if (tmp.converged)
        {
          individual_estimates.mu.push(tmp.mu);
          individual_estimates.sw.push(tmp.sigma);
          individual_estimates.worker_index.push(i);
          individual_estimates.worker.push(worker);
          
          theta_start.mu    = median(individual_estimates.mu);
          theta_start.sigma = median(individual_estimates.sw);          
        }
      }
    }
  }
  

  // From the above-obtained individual estimates, derive initial values for the *real thing* (also by Newton-Raphson)
  
  
  var mu_overall    =   mean(individual_estimates.mu);
  var sigma_between =     sd(individual_estimates.mu);
  var sigma_within  = median(individual_estimates.sw);
  
  //console.log("Indiv Estimates", individual_estimates);
  //console.log("available Individual Estimates", availableIndividualEstimates);
  
  
  var mu_worker_centered = individual_estimates.mu.map(m => m - mu_overall);
  var mu_worker = rnorm(workers.length, 0, sigma_within);
  
  for (let j=0; j<individual_estimates.worker_index.length; j++)
  {
    let i = individual_estimates.worker_index[j];
    mu_worker[i] = mu_worker_centered[j];
  }
    
    
  // Run a univariate Newton-Raphson algorithm to find an initial value
  // for sigma_within, with fixed values for mu_overall & mu_worker
  // as calculated above

  cont = true;
  var iter = 0;
  var max_niter = 250;
  
  
  if (sp.uupOnSds)
  {
    let tmp = sigma_within_RespectingRange(sigma_within, sp);
    sigma_within = tmp.sigma_within;
  }
  

  var converged = false;
  var safer_sigma_within = sigma_within; // store a copy of current value in case the algorithm below does not converge
  

  while (cont)
  {
    let tmp = dl_dsigmaw_fct(sigma_within, mu_worker, mu_overall, data, id, sp);       
    var change = - tmp.f / Math.abs(tmp.fp);

    
    if (!isNaN(change))
    {
      sigma_within -= change;
      
      if (sigma_within < 0)
      {
        let previous_sw = sigma_within + change;
        sigma_within = previous_sw / 2;
        change = sigma_within - previous_sw;
      }
    
      converged = Math.abs(change) < epsilon; 
      
      iter++;
      cont = !converged & iter < max_niter
    
      if (sp.uupOnSds)
      {
        let tmp = sigma_within_RespectingRange(sigma_within, sp);
      
        if (tmp.out_of_range)
        {
          sigma_within = tmp.sigma_within;
          cont = iter < max_niter;
        }
      }
    }
    else
    {
      cont = false;
    }
  }
  

  if (!converged) sigma_within = safer_sigma_within; 
      
      
  // Take optimal sigma_within value (given other parm values fixed)
  // to ease later convergence
  
  let tmp = optimal_sigma_within(sigma_within, mu_worker, mu_overall, data, id, sp);
  if (tmp.convergence) sigma_within = tmp.sigma_within;

  
  // Make sure that sigma_within lies on
  // the right side of its optimal value (that is, shows a negative slope)

  var sw_lower_limit = 0;
  cont = true;
  
  if (sp.uupOnSds)
  {
    if (sigma_within  > sp.sigmaWithinRange[1])  sigma_within  = sp.sigmaWithinRange[1];
    if (sigma_between > sp.sigmaBetweenRange[1]) sigma_between = sp.sigmaBetweenRange[1];
    
    sw_lower_limit = sp.sigmaWithinRange[0];
  }
      

  while (cont)
  {
    let tmp = dl_dsigmaw_fct(sigma_within, mu_worker, mu_overall, data, id, sp);
    cont = tmp.fp > 0;

    if (cont) 
    {
      sigma_within = (sigma_within + sw_lower_limit) / 2;
      if (Math.abs(sigma_within - sw_lower_limit) < epsilon) cont = false;
    }
  }

  
  // Newton-Raphson algorithm to find optimal theta
  
  var theta = mu_worker.slice(); // copy by value
  theta.push(mu_overall, sigma_within);
  var w = seq(mu_worker.length);
  
  
  var J = create_matrix(1, 1); // not the true dimensions -- this is just to create the MyMatrix object
  var J_bottom_left = create_matrix(1, workers.length);
  
  cont = true;
  iter = 0;
  
  
  if (sp.uupOnSds) sigma_between = sd(mu_worker);
  else 
  {
    let tmp = optimal_sigma_between(mu_worker, sp);
    
    if (tmp.converged) sigma_between = tmp.sigma_between;
    else sigma_between = sd(mu_worker);
  }
  
    
  while (cont)
  {       
    let tmp = optimal_sigma_within(sigma_within, mu_worker, mu_overall, data, id, sp);
    if (tmp.converged) sigma_within = tmp.sigma_within; // else we keep last iteration result
    
    // We do not optimize for sigma_between at this stage, as it would inevitably converge towards 0
          
    tmp = f12(mu_worker, mu_overall, sigma_within, sigma_between, count, y_sum, data, id);
    
    var dl_dmui         = tmp.dl_dmui,
        dl_dmu          = tmp.dl_dmu,
        d2l_dmui2       = tmp.d2l_dmui2,
        d2l_dmu2        = tmp.d2l_dmu2,
        d2l_dmuidmu     = tmp.d2l_dmuidmu;
          
//        d2l_dmuidsigmaw = tmp.d2l_dmuidsigmaw,
//        d2l_dmudsigmaw  = tmp.d2l_dmudsigmaw;
        

   
   
    // Start building Jacobian (J)
     
    J.m = diag(d2l_dmui2, false); // matrix of dimension #workers x #workers (top left part of J)
     
      J_bottom_left.m[0] = d2l_dmuidmu;  // matrix of dimensions 1 x #workers
      
      var J_top_right = J_bottom_left.transpose();
     
    J = J.rbind(J_bottom_left).cbind(J_top_right); // missing bottom-right corner (2x2) to be filled below  
    
    J.m[mu_worker.length].push(d2l_dmu2); // Complete Jacobian

    
    let previous_theta = theta.slice();
    
    var dl_dtheta = dl_dmui.slice();
    dl_dtheta.push(dl_dmu);
    
    var change = NewtonRaphsonChange(J, dl_dtheta);  
    theta = substract(theta, change);
        
    mu_worker = elements(theta, w);
    mu_overall = theta[mu_worker.length];
    
    
    iter++;
    let max_change = NewtonRaphson_max_change(theta, previous_theta);
    
    converged = max_change < epsilon;
    cont = !converged & iter < max_niter;   
  }


  var inits = {mu_overall: mu_overall, mu_worker: mu_worker, 
               sigma_within: sigma_within, sigma_between: sigma_between,
               converged: converged, niter: iter};

  return inits;
} // end of MAP_Inits