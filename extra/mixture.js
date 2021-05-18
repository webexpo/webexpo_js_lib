// Mixture algorithm
//
// Author -- Patrick Bélisle
//                      
// Version 0.4 (May 2021) 
//   [distributed]                                              
                                                                              

// Change Log
// =========================
//
// Version 0.4 (May 2021)
// ----------------------
//   - Slight correction in fct latent_data_lt (changed data.lt -> lt)
//
//   Changed the following fct calls, as they are not defined through Array.prototype anymore:
//     - CI [renamed sample_quantiles]
//     - divided_by [renamed ratio]
//     - phi
//     - Phi
//
//   Calls to d2sumLogB_dtheta2        were changed to calls to d2LogB_dtheta2
//   Calls to d2sumLogB_dtheta1dtheta2 were changed to calls to d2LogB_dtheta1dtheta2
//
//   Added fcts
//     - log_post, omega_start, optimal_sigma
//
//
// Version 0.3 (Apr 2021)
// ----------------------
//   sigma is now sampled from its full conditional distribution through means
//   of a two-step algorithm (sampling from the log-normal part of the distrn fct,
//   the sampled value being accepted with probability proportional to the second
//   part of the distrn fct) => hence no need for ref_points nor rgen_icdf anymore [except in fct rbeta, 
//   used for the generation of random values for omega].
//   Not only that makes the algorithm more simple, but it also makes it yet (~5 times) faster! 
//
// Version 0.2 (Apr 2021)
// ----------------------
//   We are now generating random values for the unobserved/censored scores
//   at each iteration, making it much easier to sample values from the full conditional distributions
//   for mu, sigma & omega.
//   That usually comes at the cost of higher autocorrelation but it does not really show here,
//   and the resulting algorithm is *much* faster (~60 times) than it was before.
//     - Dropped the fcts muGen* and omegaGen*, made unnecessary due to the recognizable form of
//       the full conditional distrns of mu & omega when all the values "y" are known [which is what
//        we assume at each iteration].
//
// Version 0.1 (Mar 2021)
// ----------------------
//   Original code.


Inits = function(data, parms)
{
  const n = {y:  data.y.length, 
             lt: data.lt.length, 
             interval: data.interval.lt.length, 
             gt: data.gt.length};
             
  const m = n.y + n.interval + n.gt;
             
  const max_niter = 100;
  const epsilon = 1e-4;
               
  var    mu_prior = new logf_norm(parms.mu.mean, parms.mu.sd);
  var sigma_prior = new logf_lnorm(parms.ls.mean, parms.ls.sd);
  
  var simplified_data = {y: data.y, gt: data.gt, interval: data.interval, lt: []}; // dropping 'lt'
  
  var tmp = OneSubjectEstimates(simplified_data, mu_prior, sigma_prior);
  
  var observed_ltprop = n.lt / data.N;
  var omega;


  
  if (tmp.converged)  mu = tmp.mu,             sigma = tmp.sigma;
  else                mu = tmp.best_guess.mu,  sigma = tmp.best_guess.sigma;
//  else
//  {
//    // Take prior mode of either mu or sigma (and take the other parm value that maximizes posterior)
//    
//    mu = parms.mu.mean;
//    sigma = Math.exp(parms.ls.mean - parms.ls.sd**2);
//    
//    let tmp1 =    OneSubjectEstimates_fixed_mu(simplified_data, mu, sigma, sigma_prior);
//    let tmp2 = OneSubjectEstimates_fixed_sigma(simplified_data, mu, sigma, mu_prior);
//    
//      tmp1.mu = mu;
//      tmp2.sigma = sigma;
//    
//    // Consider also the 'first guess' obtained from OneSubjectEstimates
//    
//    mu    = tmp.first_guess.mu;
//    sigma = tmp.first_guess.sigma;
//    
//      let tmp3 =    OneSubjectEstimates_fixed_mu(simplified_data, mu, sigma, sigma_prior);
//      let tmp4 = OneSubjectEstimates_fixed_sigma(simplified_data, mu, sigma, mu_prior);
//      
//        tmp3.mu = mu;
//        tmp4.sigma = sigma;
//
//    
//    let best_log_post = - Infinity;
//
//    
//    if (tmp1.converged)
//    {
//      let o = omega_start(data, tmp1.mu, tmp1.sigma, observed_ltprop, parms);
//      let my_log_post = log_post(simplified_data, tmp1.mu, tmp1.sigma, o, mu_prior, sigma_prior, parms);
//
//      if (my_log_post > best_log_post)
//      {
//        best_log_post = my_log_post;
//        mu = tmp1.mu;
//        sigma = tmp1.sigma;
//        omega = o;
//      }
//    }
//
//    if (tmp2.converged)
//    {
//      let o = omega_start(data, tmp2.mu, tmp2.sigma, observed_ltprop, parms);
//      let my_log_post = log_post(simplified_data, tmp2.mu, tmp2.sigma, o, mu_prior, sigma_prior, parms);
//
//      if (my_log_post > best_log_post)
//      {
//        best_log_post = my_log_post;
//        mu = tmp2.mu;
//        sigma = tmp2.sigma;
//        omega = o;
//      }
//    }
//    
//    if (tmp3.converged)
//    {
//      let o = omega_start(data, tmp3.mu, tmp3.sigma, observed_ltprop, parms);
//      let my_log_post = log_post(simplified_data, tmp3.mu, tmp3.sigma, o, mu_prior, sigma_prior, parms);
//
//      if (my_log_post > best_log_post)
//      {
//        best_log_post = my_log_post;
//        mu = tmp3.mu;
//        sigma = tmp3.sigma;
//        omega = o;
//      }
//    }
//    
//    if (tmp4.converged)
//    {
//      let o = omega_start(data, tmp4.mu, tmp4.sigma, observed_ltprop, parms);
//      let my_log_post = log_post(simplified_data, tmp4.mu, tmp4.sigma, o, mu_prior, sigma_prior, parms);
//
//      if (my_log_post > best_log_post)
//      {
//        best_log_post = my_log_post;
//        mu = tmp4.mu;
//        sigma = tmp4.sigma;
//        omega = o;
//      }
//    }  
//  }

  
  if (typeof omega == 'undefined') omega = omega_start(data, mu, sigma, observed_ltprop, parms);
  
  
  // We are now ready to run a multivariate Newton-Raphson algorithm to find the Max Posterior Estimates,
  // which will serve as inits
  
    // in case it does not converge, we save a copy of the above initial values
    var inits = {mu: mu, sigma: sigma, omega: omega, converged: false, iter: -1};
    
  
  var cont = true,
      iter = 0,
      theta = [mu, sigma, omega],
      omega_a = parms.omega.alpha - 1,
      omega_b = parms.omega.beta  - 1;
      
  var J = create_matrix(1,1); // incomplete dimension -- call create_matrix simply in order to create the MyMatrix object

  
  
  while (cont)
  {  
    let sigma2 = sigma**2;
    let sigma3 = sigma**3;
    let log_sigma = Math.log(sigma);
    let y_var = data.sum_yi2 - 2*mu*data.sum_yi + n.y*mu**2;
    
    dl_dmu    = (data.sum_yi - n.y*mu) / sigma2 - parms.mu.prec * (mu - parms.mu.mean);
    
    dl_dsigma = y_var / sigma3 - (n.y + 1)/sigma  - (log_sigma - parms.ls.mean) / sigma * parms.ls.prec;
    dl_domega = - (m + omega_b) / (1-omega);
      if (omega_a > 0) dl_domega += omega_a / omega;
    
    d2l_dmu2    = - n.y / sigma2 - parms.mu.prec;
    d2l_dsigma2 = -3*y_var/sigma**4 + (n.y+1) / sigma2 - (parms.ls.mean + 1 - log_sigma) / sigma2 * parms.ls.prec;
    d2l_domega2 = - (m + omega_b) / (1-omega)**2;
      if (omega_a > 0) d2l_domega2 -= omega_a / omega**2;
    
    d2l_dmudsigma    = -2*(data.sum_yi - n.y*mu) / sigma3;
    d2l_dmudomega    = 0;
    d2l_dsigmadomega = 0;
   
    
    if (n.lt > 0)
    {
      let z = data.lt.map(y => (y-mu)/sigma);
      let z2m1 = z.map(x => x**2 - 1);
      let Phi_z = Phi(z);
      let Phi_minus1 = Phi_z.map(f => f-1);
      let phi_z = phi(z);
      let B = Phi_z.map(f => omega + (1-omega)*f);
      let B2 = B.map(b => b**2);
      let zphi = z.map((z,i) => z*phi_z[i]);
      
      Bp_mu = phi_z.map(f => (omega - 1) * f / sigma);
      Bs_mu = zphi.map(z => z * (omega - 1) / sigma2);
      
      Bp_sigma = zphi.map(z => z * (omega - 1) / sigma);
      Bs_sigma = zphi.map((f, i) => f * (z[i]**2 - 2) * (omega-1) / sigma2);
            
      dl_dmu    += ratio_sum(Bp_mu, B);
      dl_dsigma -= (1-omega) / sigma * ratio_sum(zphi, B);
      dl_domega -= ratio_sum(Phi_minus1, B);
      
      d2l_dmu2    += d2LogB_dtheta2(B, Bp_mu,    Bs_mu);
      d2l_dsigma2 += d2LogB_dtheta2(B, Bp_sigma, Bs_sigma);
      d2l_domega2 -= ratio_sqSum(Phi_minus1, B);
      
      Bm = z.map((z,i) => (z**2 - 1) * phi_z[i] * (omega-1) / sigma2);
      d2l_dmudsigma += d2LogB_dtheta1dtheta2(B, Bp_mu, Bp_sigma, Bm);
      
      d2l_dmudomega    += ratio_sum(phi_z, B2) / sigma;
      d2l_dsigmadomega += ratio_sum(zphi, B2)  / sigma;
    }
    
    
    if (n.interval > 0)
    {
      dl_dmu        +=        dlogPhiInterval_dmu(data.interval, mu, sigma);
      dl_dsigma     +=     dlogPhiInterval_dsigma(data.interval, mu, sigma);
 
      d2l_dmu2      +=      d2logPhiInterval_dmu2(data.interval, mu, sigma);
      d2l_dsigma2   +=   d2logPhiInterval_dsigma2(data.interval, mu, sigma);
      
      d2l_dmudsigma += d2logPhiInterval_dmudsigma(data.interval, mu, sigma);
    }
    
    
    if (n.gt > 0)
    {
      let z = data.gt.map(y => (mu-y)/sigma);
      let my_phi = phi(z);
      let vphi = varphi(z).f;
      let zvarphi = z.map((z,i) => z * vphi[i]);

      // mu
      dl_dmu   += sum(vphi) / sigma;
      d2l_dmu2 -= sum(zvarphi)/sigma2 + sqSum(vphi)/sigma2;
      
      // sigma
      dl_dsigma   -= sum(zvarphi) / sigma;
      d2l_dsigma2 -= sum(z.map((z, i) => (z**2 - 2) * zvarphi[i])) / sigma2 + sqSum(zvarphi)/sigma2;
      
      // mixed mu-sigma
      let B = Phi(z);
      let Bp1 = my_phi.map(x => x/sigma);
      let Bp2 = my_phi.map((x,i) => - x * z[i] / sigma);
      let Bm = z.map((z,i) => (z**2 - 1) * my_phi[i] / sigma2);
      d2l_dmudsigma += d2LogB_dtheta1dtheta2(B, Bp1, Bp2, Bm);
    }
    
    
    // Define Jacobian
    
    J.m = diag([d2l_dmu2, d2l_dsigma2, d2l_domega2], false);
    
      J.m[1][0] = d2l_dmudsigma;
      J.m[2][0] = d2l_dmudomega;
      J.m[2][1] = d2l_dsigmadomega;
    
    J = J.sym_filled(); 
    
    f = [dl_dmu, dl_dsigma, dl_domega]; 
    change = NewtonRaphsonChange(J, f);

    previous_theta = theta.slice();
    
    theta = substract(theta, change);
    
      mu = theta[0];
      sigma = theta[1];
      omega = theta[2];
    
    // Make sure that each parameter remains in its mathematical domain
    
    if (sigma < 0)
    {
      let previous_sigma = sigma + change[1];
    
      if (omega >= 0 && omega <= 1)
      {
        let tmp = optimal_sigma(data, n, parms, mu, omega, previous_sigma/2);
        if (tmp.converged) sigma = tmp.sigma;
      }
    
      if (sigma < 0) sigma = previous_sigma / 2;
      theta[1] = sigma;
    }
    
    if (omega < 0)
    {
      omega = (omega + change[2]) / 2; // previous value / 2
      theta[2] = omega;
    }
    else if (omega > 1)
    {
      omega = (1 + omega + change[2]) / 2; // midpoint between 1 & previous value
      theta[2] = omega;
    }
    
    
    iter++;
    let max_change = NewtonRaphson_max_change(theta, previous_theta);

    converged = max_change < epsilon;    
    cont = !converged && iter < max_niter;
  }
  
  
  if (converged) inits = {mu: mu, sigma: sigma, omega: omega, converged: true, iter: iter};
  else inits.iter = iter;

  
  return inits;
} // end of Inits


latent_data = function(any, data, mu, sigma)
{
  var n0 = 0,
      sum_yi  = 0,
      sum_yi2 = 0;

  if (any.lt)
  {       
    let tmp = latent_data_lt(data.lt, mu, sigma);
      sum_yi  += tmp.sum_yi;
      sum_yi2 += tmp.sum_yi2;
      n0 = tmp.n0;
  }
  
  if (any.interval)
  { 
    let tmp = latent_data_interval(data.interval, mu, sigma);
      sum_yi  += tmp.sum_yi;
      sum_yi2 += tmp.sum_yi2;
  }
  
  if (any.gt)
  {      
    let tmp = latent_data_gt(data.gt, mu, sigma);
      sum_yi  += tmp.sum_yi;
      sum_yi2 += tmp.sum_yi2;
  }
    
    
  return {n0: n0, sum_yi: sum_yi, sum_yi2: sum_yi2};
} // end of latent_data


latent_data_gt = function(gt, mu, sigma)
{
  var sum_yi = 0, sum_yi2 = 0;
  var p_gt = Phi(gt, mu, sigma, false);
  
  for (let i=0; i<gt.length; i++)
  {
    let U = runif(1);
    let yi = qnorm(U*p_gt[i], mu, sigma, false);
    sum_yi  += yi;
    sum_yi2 += yi**2;
  }

  
  return {sum_yi: sum_yi, sum_yi2: sum_yi2};
} // end of latent_data_gt


latent_data_interval = function(interval, mu, sigma)
{
  var Phi1 = Phi(interval.gt, mu, sigma);
  var Phi2 = Phi(interval.lt, mu, sigma);
  
  var u = runif(interval.gt.length, Phi1, Phi2);
  
  var y = qnorm(u, mu, sigma);
  
  return {sum_yi: sum(y), sum_yi2: sqSum(y)};
} // end of latent_data_interval


latent_data_lt = function(lt, mu, sigma)
{ 
  var n0 = 0, sum_yi = 0, sum_yi2 = 0;
  var p_lt = Phi(lt, mu, sigma);
  var pr1 = p_lt.map(f => (1-omega)*f / (omega + (1-omega)*f));
  
  for (let i=0; i<lt.length; i++)
  {
    let U = runif(1);
    
    if (U < pr1[i])
    {
      let U = runif(1);
      let yi = qnorm(U*p_lt[i], mu, sigma);
      sum_yi  += yi;
      sum_yi2 += yi**2;
    }
    else n0++;
  }
  
    
  return {n0: n0, sum_yi: sum_yi, sum_yi2: sum_yi2};
} // end of latent_data_lt


function log_post(data, mu, sigma, omega, mu_prior, sigma_prior, parms)
{
  var log_post = 0;
  var n = data.y.length;
  var n0 = data.lt.length;
  
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
  
  var omega_a = parms.omega.alpha + n0;
  var omega_b = parms.omega.beta + data.N - n0;

  return log_post;
} // end of log_post


Mixture = function(data, parms, mcmc)
{
  const any = {lt:       data.lt.length          > 0, 
               interval: data.interval.lt.length > 0, 
               gt:       data.gt.length          > 0,
               censored: data.lt.length > 0 || data.interval.lt.length > 0 || data.gt.length > 0};
   
  const N = data.N;
  
  var delta,
      lambda, lambda_star,
      mu_mean, mu_sd,
      n0 = 0,
      omega_alpha, omega_beta,
      sum_yi  = data.sum_yi, 
      sum_yi2 = data.sum_yi2,
      y_len;


  // Compute Initial Values
    
  var inits = Inits(data, parms);
    console.log("*** Inits", inits);
  
    mu    = inits.mu;
    sigma = inits.sigma;
    omega = inits.omega;
    
    
  var sample = {mu: [], sigma: [],  omega: []};
  var burnin = {mu: [], sigma: [],  omega: []};
  
  
  // ICI tempo
  //mcmc.burnin = 0;
  //mcmc.niter = 1;
    
  
  var debug = {node: 'omega',  iter: -1,
               // debug.iter < 0 => no debugging
               
               icdf             : false,
               mode             : false,
               remote_left      : false,
               remote_right     : false,
               safe_mode        : false,
               safe_remoteLeft  : false,
               safe_remoteRight : false,
               
               bazooka          : false,
               bulldozer        : false,
               flags            : false};
  
  
  if      (debug.safe_remoteLeft)  debug.code = 1;
  else if (debug.safe_remoteRight) debug.code = 2;
  else if (debug.safe_mode)        debug.code = 3;
  else if (debug.icdf)             debug.code = 4;
  else if (debug.mode)             debug.code = 5;
  else if (debug.remote_right)     debug.code = 6;
  else if (debug.remote_left)      debug.code = 7;
  
  if (debug.bulldozer) debug.flags = true;
    
 
  // Start MCMC sampling
      
  for (let iter=0; iter < mcmc.niter + mcmc.burnin; iter++)
  {  
    if (debug.flags || iter == debug.iter || (iter%10000 == 0 && iter > 0)) console.log("Iteration #", iter);


    // Generate random values for unobserved/censored data
    
    if (any.censored)
    { 
      let tmp = latent_data(any, data, mu, sigma);  
          
      sum_yi  = data.sum_yi  + tmp.sum_yi;
      sum_yi2 = data.sum_yi2 + tmp.sum_yi2;
      n0 = tmp.n0;
    }
    
    y_len = N - n0;
    

    // Sample from f(mu | other)
                
    lambda = 1 / sigma**2;
    lambda_star = y_len * lambda + parms.mu.prec;
    delta = lambda * sum_yi + parms.mu.mean * parms.mu.prec;
    
    mu_mean = delta / lambda_star;
    mu_sd   = 1 / Math.sqrt(lambda_star); 
    
    mu = rnorm(1, mu_mean, mu_sd);
        
  
    // Sample sigma from f(sigma | mu, omega)
       
        if      (debug.node == 'sigma' && debug.iter == iter && debug.bazooka) debug_code = -1;
        else if (debug.node == 'sigma' && debug.iter == iter)                  debug_code = debug.code;
        else                                                                   debug_code = debug.bulldozer ? -1 : 0;
      
        if (debug.flags || debug_code != 0) console.log("**** sampling sigma *****\nmu, omega==", mu, omega);
    
  
    sigma = sigmaGen(mu, parms.ls, y_len, sum_yi, sum_yi2);
    
    
    // Sample from f(omega | other)
    
    omega_alpha = parms.omega.alpha     + n0;
    omega_beta  = parms.omega.beta  + N - n0;
    
    omega = rbeta(omega_alpha, omega_beta, 1, [mu, sigma], iter);
    
    
    // Save/monitor parameter values
    
    if (iter < mcmc.burnin)
    {
      if (mcmc.monitor_burnin)
      {
        burnin.mu.push(mu);
        burnin.sigma.push(sigma);
        burnin.omega.push(omega);
      }
    }
    else
    {
      sample.mu.push(mu);                         
      sample.sigma.push(sigma);                         
      sample.omega.push(omega);      
    }                                              
  } 
  
  
  desc_stats(sample);
    
  return {sample: sample, burnin: burnin};
} // end of Mixture


function omega_start(data, mu, sigma, observed_ltprop, parms, epsilon=1e-5, max_niter=100)
{
  predicted_ltprop = mean(Phi(data.lt, mu, sigma));
  
  if (observed_ltprop > predicted_ltprop) omega = (observed_ltprop - predicted_ltprop) / (1 - predicted_ltprop);
  else omega = 0;
  
  omega0 = omega;
    
  var omega_a = parms.omega.alpha - 1,
      omega_b = parms.omega.beta  - 1,
      m = data.interval.gt.length + data.gt.length + data.y.length;
  
  
  // Run a univariate Newton-Raphson algorithm to find the optimal value for omega with fixed (mu, sigma)
 
  var cont = true, 
      iter = 0;
  
  const z = data.lt.map(y => (y-mu)/sigma);
  const Phi_z = Phi(z);
  const Phi_minus1 = Phi_z.map(f => f-1); 
  
 
  while (cont)
  {
    dl_domega = - (m + omega_b) / (1-omega);
      if (omega_a > 0) dl_domega += omega_a / omega;
    
    d2l_domega2 = - (m + omega_b) / (1-omega)**2;
      if (omega_a > 0) d2l_domega2 -= omega_a / omega**2;
      
      
    let B = Phi_z.map(f => omega + (1-omega)*f);
    
    dl_domega   -=   ratio_sum(Phi_minus1, B); 
    d2l_domega2 -= ratio_sqSum(Phi_minus1, B);    

  
    change = - dl_domega / Math.abs(d2l_domega2);
    omega -= change;
    
    if      (omega < 0) omega =     (omega + change) / 2;
    else if (omega > 1) omega = (1 + omega + change) / 2;
      
    iter++;
    converged = Math.abs(change) < epsilon;
    cont = !converged && iter < max_niter;
  }
  
  if (!converged) omega = omega0;
  
  return omega;
} // end of omega_start


function optimal_sigma(data, n, parms, mu, omega, sigma, max_niter=100, epsilon=1e-5)
{
  // Find sigma that maximizes log posterior when (mu, omega) are fixed
  
  var cont = true,
      converged = false,
      iter = 0;
  
  
  while (cont)
  {  
    let sigma2 = sigma**2;
    let sigma3 = sigma**3;
    let sigma4 = sigma**4;
    let log_sigma = Math.log(sigma);
    let y_var = data.sum_yi2 - 2*mu*data.sum_yi + n.y*mu**2;
    

    dl_dsigma   =    y_var / sigma3 - (n.y + 1) / sigma  -   (log_sigma - parms.ls.mean)   / sigma  * parms.ls.prec;
    d2l_dsigma2 = -3*y_var / sigma4 + (n.y + 1) / sigma2 - (parms.ls.mean + 1 - log_sigma) / sigma2 * parms.ls.prec;

    
    if (n.lt > 0)
    {
      let z = data.lt.map(y => (y-mu)/sigma);
      let Phi_z = Phi(z);
      let phi_z = phi(z);
      let B = Phi_z.map(f => omega + (1-omega)*f);
      let zphi = z.map((z,i) => z*phi_z[i]);
      
      Bp = zphi.map(z => z * (omega - 1) / sigma);
      Bs = zphi.map((f, i) => f * (z[i]**2 - 2) * (omega-1) / sigma2);
            
      dl_dsigma   -= (1-omega) / sigma * ratio_sum(zphi, B);
      d2l_dsigma2 += d2LogB_dtheta2(B, Bp, Bs);
    }
    
    
    if (n.interval > 0)
    {
      dl_dsigma   +=   dlogPhiInterval_dsigma(data.interval, mu, sigma);
      d2l_dsigma2 += d2logPhiInterval_dsigma2(data.interval, mu, sigma);
    }
    
    
    if (n.gt > 0)
    {
      let z = data.gt.map(y => (mu-y)/sigma);
      let vphi = varphi(z).f;
      let zvarphi = z.map((z,i) => z * vphi[i]);

      dl_dsigma   -= sum(zvarphi) / sigma;
      d2l_dsigma2 -= sum(z.map((z, i) => (z**2 - 2) * zvarphi[i])) / sigma2 + sqSum(zvarphi)/sigma2;
    }
    
    change = - dl_dsigma / Math.abs(d2l_dsigma2);
    sigma -= change;
    
    if (sigma < 0) sigma = (sigma + change) / 2;
    
    converged = Math.abs(change) < epsilon;
    cont = !converged && ++iter < max_niter;
  }
  
  return {converged: converged, sigma: sigma};
} // end of optimal_sigma


sigmaGen = function(mu, ls_parms, y_len, sum_yi, sum_yi2)
{
  if (y_len == 0)
  {
    let log_sigma = rnorm(1, ls_parms.mean, ls_parms.sd);
    return Math.exp(log_sigma);  
  }


  tau_mean = ls_parms.mean * (1 + y_len);
  tau_sd   = ls_parms.sd   * (1 + y_len);
  
  beta = sum_yi2 - 2 * sum_yi * mu + y_len * mu**2;
  tau_e = 1 / (1 + y_len);
  
  wt_mode = Math.sqrt(beta/y_len);
  log_M = -beta/2/wt_mode**2 - y_len * Math.log(wt_mode);
  
  cont = true;
  
  while (cont)
  {
    let log_tau = rnorm(1, tau_mean, tau_sd);
    sigma = Math.exp(tau_e * log_tau);
    let log_U = Math.log(runif(1));
    let accepted = log_U <= -beta/2/sigma**2 - y_len * log_tau / (1 + y_len) - log_M;
    cont = !accepted;
  }

  return sigma;
} // end of sigmaGen



// =============================================================================
// ==   Main function


Run = function()
{
  var mu, 
      sigma, 
      omega = 0.1; // just to start the search for mode on first call to omegaGen
      
  var sample = {mu: [], sigma: [], omega: []};

  console.clear();
  DisplayLocalTime();
  HideResultsTextbox();
  
  var t0 = performance.now();

  
  //////////////////////////////////////////////////////////////////////////////
  // Read Data and relevant info from Html Form
  
  // Read data
    
  var data = ReadData(document.form);

  if (data.N == 0)
  {
    ErrorMsg(0);
    return;                                       
  }
  
  
  // Regroup parameters 
   
  var parms = {mu: {mean: Number(document.form.mu0.value), 
                    sd:   Number(document.form.sigma0.value),
                    prec: 0}, // temp
               ls: {mean: Number(document.form.logSigmaMu.value), 
                    sd:   Number(document.form.logSigmaSD.value),
                    prec: 0}, // temp
               omega: {alpha: NaN,
                       beta: NaN,
                       uniform_prior: false}
               };
               
  parms.mu.prec  = 1 / parms.mu.sd**2;
  parms.ls.prec  = 1 / parms.ls.sd**2;

  
  var mcmc = MCMCParms(document.form);
  
  var logNormalDistrn = document.form.logNormalDistrn.value == 1;
  
  
  parms.omega.uniform_prior = document.form.uniformOmegaPrior.value == 1;


  if (parms.omega.uniform_prior)
  {
    parms.omega.alpha = 1;
    parms.omega.beta  = 1;
  }
  else
  {
    if (document.form.omegaPriorLCL.value.length == 0 || document.form.omegaPriorUCL.value.length == 0)
    { 
      ErrorMsg(1);
      return;
    }
    
    let lcl = Number(document.form.omegaPriorLCL.value);
    let ucl = Number(document.form.omegaPriorUCL.value);

    
    // Validate LCL & UCL entries
    
    let level = Number(document.form.level.value);
    
    if (lcl >= ucl)               
    {
      ErrorMsg(2);
      return;
    }
    else if (lcl <= 0 || ucl >= 1)
    { 
      ErrorMsg(3);
      return;
    }
    
    let tmp = BetaParmsFromQuantiles(lcl, ucl, level);
    
      if (tmp.converged)
      {
        parms.omega.alpha = tmp.alpha;
        parms.omega.beta  = tmp.beta;
      }
      else
      {
        ErrorMsg(6);
        return;
      }
  }
  

  //////////////////////////////////////////////////////////////////////////////
  
 
  if (logNormalDistrn) 
  {
    // First check that no data value is negative

    if (data.any_le0())
    { 
      ErrorMsg(4);
      return; 
    }
  
    data.TakeLog();
  }
  
  data.sum_yi  = sum(data.y);
  data.sum_yi2 = sqSum(data.y);
  
  // Prepare Model Output(s)
  
  var out = {mixture: {}};
  
  // Run Model(s) 

  out.mixture = Mixture(data, parms, mcmc);
 

  WorkComplete(out.mixture.sample, out.mixture.burnin, mcmc);
     
  DisplayRuntime(t0, undefined, mcmc.burnin + mcmc.niter);
  
  
  // ___________________________________________________________________________
  

   console.log("Goodbye!");
} // end of Run