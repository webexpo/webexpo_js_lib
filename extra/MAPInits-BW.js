// =============================================================================
// Code written to calculate Initial Values for the Between-Worker model
// using the MAP (Maximum Posterior Estimate) as initial values                                                              
                                                        
// Version 0.2 (Feb 2021) 
// [based on inits-BetweenWorker.R version 0.25]                            

// -----------------------------------------------------------------------------
// PB, Feb 2021

// Change Log
//
// Version 0.2 (Feb 2021)
// ----------------------
//   Now use the library lib/normal-OnesubjectEstimates rather than 
//   the previously embedded fct OneWorkerEstimates

         
dl_dsigmaw = function(sigma_within, mu_worker, mu_overall, data, id, parms)
{ 
  var f  = - id.y.length / sigma_within,
      fp =   id.y.length / sigma_within**2; 
  
  
  if (id.y.length > 0)
  {
    let mu_w = mu_worker.indexed(id.y); 
    let tmp = data.y.minus(mu_w).minus(mu_overall).sqSum();
    
    f  +=     tmp / sigma_within**3;
    fp -= 3 * tmp / sigma_within**4;
  }
  
  if (id.lt.length > 0)
  {
    let mu_w = mu_worker.indexed(id.lt);
    let tmp = data.lt.minus(mu_w).minus(mu_overall);
    
    f  -=   dlogPhi_dsigma(tmp, sigma_within);
    fp -= d2logPhi_dsigma2(tmp, sigma_within);
  }  
  
  if (id.gt.length > 0)
  {
    let mu_w = mu_worker.indexed(id.gt);
    let tmp = mu_w.plus(mu_overall).minus(data.gt);
    
    f  -=   dlogPhi_dsigma(tmp, sigma_within);
    fp -= d2logPhi_dsigma2(tmp, sigma_within);
  }
  
  if (id.interval.length > 0)
  {
    let mu_w = mu_worker.indexed(id.interval);
    let tmp = mu_w.plus(mu_overall);
    
    f  -= dlogPhiInt_dsigma(data.interval, tmp, sigma_within);
    fp -=        d2logPhiInt_dmudsigma(data.interval, tmp, sigma_within);
  }


  if (!parms.uupOnSds)
  {
    let logn = logN(sigma_within, parms);
    
    f  += logn.f;
    fp += logn.fp;
  }
  

  var dl_dsigmaw = {f: f, fp: fp};
  
  return dl_dsigmaw;
} // end of dl_dsigmaw


function f12(mu_worker, mu_overall, sigma_within, sigma_between, count, y_sum, data, id)
{
  var tmp = mu_worker.plus(mu_overall).times(count.y);
      tmp = y_sum.minus(tmp);                             // array of length #workers
      
  var f1 = tmp.divided_by(sigma_within**2);               // arrays of length #workers
  var df1_dmui = count.y.map(u => -u /sigma_within**2); 
  var df1_dsigmaw = tmp.map(u => -2*u / sigma_within**3);
  

  for (let i=0; i<mu_worker.length; i++)
  {
    let my_mu = mu_worker[i] + mu_overall;
    
    if (count.lt[i] > 0)
    {
      let lt = data.lt.filter((e,index) => id.lt[index]===i).minus(my_mu);
      
      f1[i]          -=        dlogPhi_dmu(lt, sigma_within);
      df1_dmui[i]    -=      d2logPhi_dmu2(lt, sigma_within);
      df1_dsigmaw[i] -= d2logPhi_dmudsigma(lt, sigma_within);
    }
    
    if (count.gt[i] > 0)
    {
      let gt = data.gt.filter((e,index) => id.gt[index]===i).minus(my_mu).times(-1);
      
      f1[i]          +=        dlogPhi_dmu(gt, sigma_within);
      df1_dmui[i]    +=      d2logPhi_dmu2(gt, sigma_within, 1);
      df1_dsigmaw[i] += d2logPhi_dmudsigma(gt, sigma_within);
    }
    
    if (count.interval[i] > 0)
    {
      let gt = data.interval.gt.filter((e,index) => id.interval[index]===i),
          lt = data.interval.lt.filter((e,index) => id.interval[index]===i);     
      let my_interval = {gt: gt, lt: lt};
      
      f1[i]          -=        dlogPhiInt_dmu(my_interval, my_mu, sigma_within);
      df1_dmui[i]    -=      d2logPhiInt_dmu2(my_interval, my_mu, sigma_within);
      df1_dsigmaw[i] -= d2logPhiInt_dmudsigma(my_interval, my_mu, sigma_within);
    }
  } 
  
  
  var df1_dmu = df1_dmui.slice();
  var f2 = f1.sum();
  var df2_dmu = df1_dmui.sum();
  var df2_dsigmaw = df1_dsigmaw.sum();
        
  f1 = f1.minus(mu_worker.divided_by(sigma_between**2));
  df1_dmui = df1_dmui.minus(1/sigma_between**2); 
  
  
  var f12 = {f1: f1, df1_dmui: df1_dmui, df1_dmu: df1_dmu, df1_dsigmaw: df1_dsigmaw,
             f2: f2, df2_dmu: df2_dmu, df2_dsigmaw: df2_dsigmaw};
  
  return f12;
} // end of f12


function logN(sigma, m)
{
  var log_sigma = Math.log(sigma);
  var f = -1/sigma * (1 + m.logSigmaWithinPrec * (log_sigma - m.logSigmaWithinMu));
  var fp = (1 - m.logSigmaWithinPrec * (m.logSigmaWithinMu + 1 - log_sigma)) / sigma**2;
  
  var logN = {f: f, fp: fp};
  
  return logN;
} // end of logN


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



////////////////////////////////////////////////////////////////////////////////

function MAP_Inits(wd, sp, epsilon=1e-6) 
{
  // Find the Maximum A Posteriori (MAP) Estimator 
  //   and use it as initial values for the Between-Workers algorithm

  
  var workers = wd.workerIds;       

  var sigmaBetweenRange = {l : 0, u : Infinity};
  var sigmaWithinRange  = {l : 0, u : Infinity};
         
  if (sp.uupOnSds)
  {
    sigmaBetweenRange = {l : sp.sigmaBetweenRange[0], u : sp.sigmaBetweenRange[1]};
    sigmaWithinRange  = {l : sp.sigmaWithinRange[0],  u : sp.sigmaWithinRange[1]};
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
    
    var this_worker_data = {y: [], lt: [], gt: [], interval: {gt: [], lt: []}};
    
    this_worker_data.y  = wk.filter(m => m.type === "uncensored").map(m => m.a);
    this_worker_data.lt = wk.filter(m => m.type === "lessThan").map(m => m.a);
    this_worker_data.gt = wk.filter(m => m.type === "greaterThan").map(m => m.a);
    
    this_worker_data.interval.gt = wk.filter(m => m.type === "interval").map(m => m.a);
    this_worker_data.interval.lt = wk.filter(m => m.type === "interval").map(m => m.b);
     
    count.y.push(this_worker_data.y.length);
    y_sum.push(this_worker_data.y.sum());
            
    count.gt.push(this_worker_data.gt.length); 
    count.lt.push(this_worker_data.lt.length);
                            
    count.interval.push(this_worker_data.interval.gt.length);
    
      // Append to combined data
      
      data.y  = data.y.concat(this_worker_data.y);
      data.lt = data.lt.concat(this_worker_data.lt);
      data.gt = data.gt.concat(this_worker_data.gt);
      data.interval.gt = data.interval.gt.concat(this_worker_data.interval.gt);
      data.interval.lt = data.interval.lt.concat(this_worker_data.interval.lt);
      
      id.y  = id.y.concat([i].rep(this_worker_data.y.length));
      id.lt = id.lt.concat([i].rep(this_worker_data.lt.length));
      id.gt = id.gt.concat([i].rep(this_worker_data.gt.length));
      id.interval = id.interval.concat([i].rep(this_worker_data.interval.gt.length));
      
                  
    var estimable_parms = this_worker_data.y.length > 1 || this_worker_data.interval.gt.length > 0 || (this_worker_data.gt.length > 0 && this_worker_data.lt.length > 0); // for this worker        
            
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
        run_1WorkerEstimate = !this_worker_data.y.allTheSame();
      }
              
    
      if (!run_1WorkerEstimate && this_worker_data.y.length > 0)
      {
        ybar = this_worker_data.y.mean();
        
        concatenated_l = this_worker_data.gt.concat(this_worker_data.interval.gt);
        concatenated_u = this_worker_data.lt.concat(this_worker_data.interval.lt);

        if (concatenated_l.length > 0)                         run_1WorkerEstimate = concatenated_l.any_gt(ybar);
        if (!run_1WorkerEstimate && concatenated_u.length > 0) run_1WorkerEstimate = concatenated_u.any_lt(ybar);
      }    
              
              
      if (!run_1WorkerEstimate)
      {
        let l_median = concatenated_l.median();
        let u_median = concatenated_u.median();
                
        thisWorkerEstimates.converged = true;
        thisWorkerEstimates.mu = ybar;
        thisWorkerEstimates.sigma = Math.abs(u_median - l_median) / 6;
      }
      else
      {
        let tmp = OneSubjectEstimates(this_worker_data);
        
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


  // From the above-obtained individual estimates, derive initial values for the *real thing* (also by Newton-Raphson)
  
  var mu_overall    = individual_estimates.mu.mean();
  var sigma_between = individual_estimates.mu.sd();
  var sigma_within  = individual_estimates.sw.median();
  
  var mu_worker_centered = individual_estimates.mu.minus(mu_overall);
  var mu_worker = [];
  var j = 0;
  
  for (i=0; i<availableIndividualEstimates.length; i++)
  {
    var worker_mean =  availableIndividualEstimates[i] ? mu_worker_centered[j++] : 0;
    mu_worker.push(worker_mean);
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

  while (cont)
  {
    let tmp = dl_dsigmaw(sigma_within, mu_worker, mu_overall, data, id, sp);       
    var change = tmp.f / tmp.fp;

    
    if (!isNaN(change))
    {
      sigma_within -= change;
    
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
    let tmp = dl_dsigmaw(sigma_within, mu_worker, mu_overall, data, id, sp);
    
    cont = tmp.fp > 0;
    if (cont) sigma_within = (sigma_within + sw_lower_limit) / 2;
  }

  sigma_within = (sigma_within + sw_lower_limit) / 2;

  
  // Newton-Raphson algorithm to find optimal theta
  
  var theta = mu_worker.slice(); // copy by value
  theta.push(mu_overall, sigma_within);
  var w = seq(mu_worker.length);
  
  
  cont = true;
  iter = 0;
    
    
  while (cont)
  {
    //  f1 = dl/dmu_i, i = 1, 2, ..., n.workers
    //  f2 = dl/dmu 
    
    let tmp = f12(mu_worker, mu_overall, sigma_within, sigma_between, count, y_sum, data, id);
    
    var f1          = tmp.f1,
        f2          = tmp.f2,
        df1_dmui    = tmp.df1_dmui,
        df1_dmu     = tmp.df1_dmu,
        df1_dsigmaw = tmp.df1_dsigmaw,
        f2          = tmp.f2,
        df2_dmu     = tmp.df2_dmu,
        df2_dsigmaw = tmp.df2_dsigmaw;
   
   
    // Start building Jacobian (J)
     
    var J = df1_dmui.diag(); // matrix of dimension #workers x #workers (top left part of J)
     
    var J_bottom_left = Array(2).fill(0).map(x => Array(0).fill(0));
    J_bottom_left[0] = df1_dmu;
    J_bottom_left[1] = df1_dsigmaw;
     
    var J_top_right = J_bottom_left.transpose();
    J = J.rbind(J_bottom_left);
    J = J.cbind(J_top_right);
    
    
    // f3 = dl/dsigmaw
    
    tmp = dl_dsigmaw(sigma_within, mu_worker, mu_overall, data, id, sp);
    
    var f3          = tmp.f,
        df3_dsigmaw = tmp.fp;
        
    // Complete Jacobian
    J[mu_worker.length].push(df2_dmu, df2_dsigmaw);
    J[mu_worker.length+1].push(df2_dsigmaw, df3_dsigmaw);
    
    
    var f = f1.slice();
    f.push(f2, f3);
    
    var change = NewtonRaphsonChange(J, f);
    theta = theta.minus(change);
    
    var previous_sigma_within = sigma_within;
    
    mu_worker = theta.indexed(w);
    mu_overall = theta[mu_worker.length];
    sigma_within = theta[mu_worker.length+1];
    
    var max_change = change.max_abs();
    iter++;
    converged = max_change < epsilon;
    cont = !converged & iter < max_niter;

    // Safety against sigma_within falling out of range
    
    if (sp.uupOnSds)
    {            
      let tmp = sigma_within_RespectingRange(sigma_within, sp, previous_sigma_within);
     
      if (tmp.out_of_range)
      {
        sigma_within = tmp.sigma_within;
        theta[mu_worker.length+1] = sigma_within;
        cont = iter < max_niter;
      } 
    }
  }


  sigma_between = mu_worker.sd();

  var inits = {mu_overall: mu_overall, mu_worker: mu_worker, 
               sigma_within: sigma_within, sigma_between: sigma_between,
               converged: converged};

  return inits;
} // end of MAP_Inits