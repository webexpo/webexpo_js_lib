// Functions to find reference points (mode & left- and right-side remote limits)
// of a distribution function
// AND to sample from any general density function 
//  (on a finite domain, defined through ref_points)                 
//
// Author: Patrick Bélisle                   

// Version 0.4 (May 2021)
//   [distributed]


// Change log
// ======================
//
// Version 0.4 (May 2021)
// ----------------------
//   All Object.prototype.* function definitions were embedded 
//   in variable type declaration DensityFcts
//
// Version 0.3 (Apr 2021)
// ----------------------
//   The fct quantile does not have a 'range' argument anymore
//
// Version 0.2 (Apr 2021)
// ----------------------
//   - renamed this file (was called dens-gen-icdf.js before)
//   - added fcts 
//       compute_M
//       error_log
//       ProtectAgainstNaN_logfp
//       quantile
//       QuickLookAtBorderAsARemoteEndpoint
//       rgen
//       rgen_icdf 
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code


const DensityFct =
{
  // This variable type defn was used in: 
  // - stats.js/rbeta


  approx_mode: function(x)
  {
    // x: starting point
    
    var logfp = this.logf_prime(x);
    var logfs = this.logf_second(x);
    var mode;
    
    
    if (this.dens_approx == 1)
    {
      // Normal distrn approx
      
      let error = logfs > 0;
        
      if (!error)
      {
        let sigma = Math.sqrt(-1/logfs);
        mode = x + logfp*sigma**2;
      }
      
      return {error: error, mode: mode};
    }
    else if (this.dens_approx == 2)
    {
      //  log-Normal distrn approx
            
      let lambda = -logfs * x**2 - x * logfp;
      let error = lambda < 0;
      
      if (!error)
      { 
        let mu = Math.log(x) + (x * logfp + 1) / lambda;
        mode = Math.exp(mu - 1/lambda);
      }
  
      return {error: error, mode: mode};
    }
    else if (this.dens_approx == 3)
    {
      // Beta distrn approx
      
      let B = -1/(1-x)**2 - 1/x/(1-x);
      let b = (logfs + logfp/x) / B;
      let a = x * (logfp + b/(1-x));    
      let error = a < 0 || b < 0;
  
      if (!error)
      {
        let alpha = a + 1;
        let beta  = b + 1;
        mode = (alpha - 1) / (alpha + beta - 2);
      }
      
      return {error: error, mode: mode};
    }
    else
    { 
      // square root inverted gamma
      
      let beta = -logfp/4 * x**3 - logfs/4 * x**4;
      let a = -3/2*logfp*x - logfs*(x**2)/2;
      let error = a < 1 || beta < 0;
        
      if (!error)
      {
        let alpha = (a-1) / 2;
        mode = Math.sqrt(beta/(alpha+1/2));
      }
      
      return {error: error, mode: mode};
    }
  }, // end of approx_mode


  area: function(upper_bound, lower_bound) 
  {
    // Before using this function, make sure that either
    //     i) this.M (the height at the mode) is defined
    // or ii) the function this.mode is defined
    
  
    if (typeof lower_bound == 'undefined')
    { 
      lower_bound = typeof this.remote_left  == 'undefined' ? this.domain[0] : this.remote_left;
    }
    
    if (typeof upper_bound == 'undefined')
    { 
      upper_bound = typeof this.remote_right == 'undefined' ? this.domain[1] : this.remote_right;
    }  
    
    
    if (typeof this.g_compute == 'undefined') this.g_compute = false;
    
  
    if (isNaN(this.M))
    {
      let mode = this.mode();
      
      if (this.g_compute)
      {
        this.M = 1; // necessary for the line below to give the right value
        this.M = this.g(mode);
      }
      else
      {
        this.M = this.logf(mode);
      }
    }
      
  
    var self = this;
    
    var F = function(xA) 
    {
      return xA.map(self.f, self);
    }
    
    // with user-specified density function
    var G = function(xA) 
    {
      return xA.map(self.g, self);
    }
  
  
    if (!this.g_compute) return new zNum.NumericIntegration(F, lower_bound, upper_bound).compute().result;
    else                 return new zNum.NumericIntegration(G, lower_bound, upper_bound).compute().result;
  }, // end of Object.area


  cdf: function(x)
  {
    if (typeof this.total_area == 'undefined') this.total_area = this.area();
    
    return this.area(x) / this.total_area;
  }, // end of cdf
  

  compute_M: function(mode = this.mode())
  {
    // this function works iff the function this.mode is defined
    
    if (typeof this.g_compute == 'undefined') this.g_compute = false;
    
    if (this.g_compute) this.M = 1, this.M = this.g(mode); // need to set it to 1 first!
    else this.M = this.logf(mode);
  }, // end of compute_M
  
  
  f: function(x) 
  {
    if (!this.g_compute) return Math.exp(this.logf(x) - this.M);
    else                 return this.g(x);
  }, // end of f


  higher_local_max: function(x_local_min, hs, epsilon)
  {
    var change,
        fp,
        local_mode = [];
    
    // Find a starting point (for N-R to work smoothly) on each side of x_local_min
    
    // i) on right side
    
    var x = x_local_min;
    var fs = hs;
    var cont = true;
  
    while (cont)
    {
      change = 1/fs;
      x += change;
      if (x > this.domain[1]) x = (x - change + this.domain[1]) / 2;
      
      fp = this.logf_prime(x);
      fs = this.logf_second(x);
      cont = fp > 0 && fs > 0;
    }  
    
    var start_right = x,
        fp_right    = fp;
        
    
    // ii) on left side
    
    x = x_local_min;
    fs = hs;
    cont = true;
    
    while (cont)
    {
      change = -1/fs;
      x += change;
      if (x < this.domain[0]) x = (x - change + this.domain[0]) / 2;
      
      fp = this.logf_prime(x);
      fs = this.logf_second(x);
      cont = fp < 0 && fs > 0;
    }
  
    
    // Run Newton-Raphson to find local max on each side
    
    // i) left side
    
    cont = true;
    
    while (cont)
    {
      change = fp / this.logf_second(x);
      x -= change;
      if (x < this.domain[0]) x = (x + change + this.domain[0]) / 2;
      
      fp = this.logf_prime(x);
      cont = Math.abs(fp) > epsilon;
  
    }
    
    local_mode.push(x);
    
    // ii) right side
    
    cont = true;
    x = start_right;
    fp = fp_right;
    
    while (cont)
    {
      change = fp / this.logf_second(x);
      x -= change;
      if (x > this.domain[1]) x = (x + change + this.domain[1]) / 2;
      
      fp = this.logf_prime(x);
      cont = Math.abs(fp) > epsilon;
    }
    
    local_mode.push(x);
    
    
    // Compute height at each local max and pick the higher one
    
    h = local_mode.map(this.logf);
    var j = h[0] > h[1] ? 0 : 1;  
    var out = {x: local_mode[j], f: h[j], remote_start: []};
    
    // Also find alternate left- and side- starting points for the search for remote endpoints
  
    fs = local_mode.map(this.logf_second);
  
    // left side
    x = local_mode[0] + 1/fs[0];
    if (x < this.domain[0]) x = (local_mode[0] + this.domain[0]) / 2;
    out.remote_start.push(x);
    
    // right side
    x = local_mode[1] - 1/fs[1];
    if (x > this.domain[1]) x = (local_mode[1] + this.domain[1]) / 2;
    out.remote_start.push(x);
    
    
    return out;  
  }, // end of higher_local_max
  
  
  mode_SafeStartingPoint: function(x, logfp, logfs, debug=false)
  {
    var cont = true;
    var iter = 0;
    
    while (cont)
    {
      if (iter++ == 250) 
      {
        error_log("mode_SafeStartingPoint", this, iter);
        return {x: NaN, error: true};
      }
      
      change = -logfp/logfs;  
      x -= change;
    
      if (x < this.domain[0]) x = (x + change + this.domain[0]) / 2;
      
      logfp = this.logf_prime(x);
      logfs = this.logf_second(x);
      
      cont = logfs > 0;
    }
    
    return {x: x, logfp: logfp, error: false};
  }, // end of mode_SafeStartingPoint
  

  ProtectAgainstNaN_logfp: function(x, previous_x)
  {
    // x: the new tentative/suggested x
    
    var change = previous_x - x;
    var logfp = this.logf_prime(x);
    var cont = isNaN(logfp);
    var iter = 0;
  
    
    while (cont)
    {
      if (iter++ == 100)
      {
        x = previous_x;
        change = 0;
        logfp = this.logf_prime(x);
        break;
      }
      
      change /= 2;
      x = previous_x - change;
      logfp = this.logf_prime(x);
      cont = isNaN(logfp);
    }
    
    return {x: x, logfp: logfp, change: change};
  }, // end of ProtectAgainstNaN_logfp 
  
  
  quantile: function(p, xstart, debug=false, epsilon=1e-8)
  {
    // p: an array or a number
    
    // The object funs [this] includes:
    //   - The functions:
    //       logf, logf_prime, logf_second
    //
    // xstart: the starting point in the search for the distrn quantile of order p -- it is a good idea to start
    //           at the mode of the density function (for a question of quick & smooth convergence) 
    //           when you don't have an informed guess. 
    //
    //  If xstart is left undefined, funs must also include
    //   mode: a function returning the mode of the density function this.logf
    //
    // IMPORTANT: the property this.M must be defined beforehand, serving as a standardizing constant
    //            in the calculation of total area under the distribution density function on its domain
    //            (defined as points within this.remote_left & this.remote right or within this.domain limits)
    //
    // See top of rgen for more info.
   
    
    if (typeof xstart == 'undefined')
    {
      let mode = this.mode();
      xstart = mode;
              
      this.M = this.logf(mode);
    } 
    
   
    var lower_bound = typeof this.remote_left  == 'undefined' ? this.domain[0] : this.remote_left;
    var upper_bound = typeof this.remote_right == 'undefined' ? this.domain[1] : this.remote_right; 
  
    if (typeof this.total_area == 'undefined') this.total_area = this.area(upper_bound, lower_bound);
    
    
      if (debug)
      { 
        console.log("range", lower_bound, upper_bound);
        console.log("domain", this.domain);
        console.log("area", this.total_area);
      }
    
    
    var p_array,
        q = [],
        x = xstart; // start at the mode for a faster & safer convergence
        
  
    
    var array_input = Array.isArray(p);
    
    if (array_input) p_array = p;
    else p_array = [p];
    
    
    var area_xstart = this.area(xstart, lower_bound);
    var u_xstart = area_xstart / this.total_area;
    var r = rank(p_array);
    
    var p_split = split(p_array, u_xstart);
      p_split.le.sort(function(a, b){return b-a}); // desc order
      p_split.gt.sort(function(a, b){return a-b}); // asc order
      
    
    for (let s=-1; s<=1; s+=2)
    {
      // when s=-1, we treat values below xstart, in descending order
      // when s=+1, we treat values above xstart, in  ascending order
      
      x = xstart;
      var myp;
      
      if (s < 0) myp = p_split.le;
      else
      {       
        previous_target = area_xstart;
        myp = p_split.gt;
      }
      
      
      for (let i=0; i<myp.length; i++)
      {
        if (s < 0) my_range = [lower_bound, x];
        else       my_range = [x, upper_bound];
        
        // Sample a value -- from Inverse CDF, by a mixed Newton-Raphson/bisectional algorithm
      
        let target = this.total_area * myp[i];
        
        if (s > 0)
        {
          let tmp = target;
          target -= previous_target;
          previous_target = tmp;
        }
        
    
        let cont = true,
            F, change,
            iter = 0;  
    
        while (cont) 
        {
          if (iter++ == 250)
          { 
            error_log(this.label + " quantile", this, iter, xstart, undefined, p, [lower_bound, upper_bound]);
            return {error: true};
          }
          
          F = this.area(x, my_range[0]);
          change = (F - target) / this.f(x);
          
          if (F < target)
          {
            my_range[0] = x;
            target -= F;  
          }
          else if (F > target) my_range[1] = x;
      
          let previous_x = x;
          x -= change;
          
          if      (x > my_range[1]) x = (x + change + my_range[1]) / 2, change = previous_x - x;
          else if (x < my_range[0]) x = (x + change + my_range[0]) / 2, change = previous_x - x;
          
          converged = Math.abs(change) < epsilon; 
          cont = !converged;
        }
    
        if (s > 0) q.push(x);
        else q.unshift(x);
      }  
    }
    
    // At this point, the values in 'q' are in ascending order:
    // we now reorganize them in the same order as their p-counterpart
  
    if (array_input)
    {
      let rearranged_q = [];
      for (let i=0; i<r.length; i++) rearranged_q.push(q[r[i]]);
      return {error: false, x: rearranged_q};
    }
    else return {error: false, x: q[0]};
  }, // end of quantile 
  
  
  QuickLookAtBorderAsARemoteEndpoint: function(target, remote_start, left_side=true)
  {
    var cont = true,  
        last_remote = left_side ? remote_start.left : remote_start.right,
        side = left_side ? 0 : 1;
        
    
    if (last_remote == this.domain[side])
    {
      let logf = this.logf(last_remote);
      if (!isNaN(logf)) cont = logf < target;
    }
  
    return {cont: cont, x: last_remote};
  }, // end of QuickLookAtBorderAsARemoteEndpoint 
  
  
  ref_points: function(xstart, remote_start, debug_code=0, f_ratio_remote=1e-8, epsilon=1e-6) 
  {  
    const ALT_EPSILON = 1e-4;
    
    var x = xstart,
        found_mode_on_domain_border = false,
        mode_border_side; // Irrelevant when mode is NOT on the border of the domain (most of the time, that is)
        
    var out = {error: false, mode: NaN, M:NaN, remote: {left: NaN, right: NaN}};
    
    // Default values
    
    if (typeof remote_start == 'undefined') remote_start = {left: this.domain[0], right: this.domain[1]};
    
    
    // Find x with log.f.prime > 0
  
    // If the previous mode was at one end of the domain, first look if it is still the case 
    // for the current MCMC iteration
      
    if (x == this.domain[0] || x == this.domain[1])
    {
      let logfp_sign_if_mode = x == this.domain[0] ? -1 : 1;
      
      logfp = this.logf_prime(x);
  
      found_mode_on_domain_border = Math.sign(logfp) == logfp_sign_if_mode;
      if (found_mode_on_domain_border) mode_border_side = x == this.domain[0] ? 0 : 1;
      else x -= logfp_sign_if_mode * ALT_EPSILON;
    }
  
  
    if (!found_mode_on_domain_border && typeof this.dens_approx != 'undefined')
    { 
      let tmp = this.approx_mode(x);
      if (!tmp.error) x = tmp.mode;
    }
   
    
    var change,
        logfp = this.logf_prime(x),
        logfs = this.logf_second(x);
        
    var debug_mode                   = debug_code == 5 || debug_code < 0,
        debug_remote_left            = debug_code == 7 || debug_code < 0,
        debug_remote_left_safeStart  = debug_code == 1 || debug_code < 0,
        debug_remote_right           = debug_code == 6 || debug_code < 0,
        debug_remote_right_safeStart = debug_code == 2 || debug_code < 0,
        debug_icdf                   = debug_code == 4 || debug_code < 0
        ; 
        
    var debug_remote = debug_remote_left || debug_remote_left_safeStart || debug_remote_right || debug_remote_right_safeStart;
      
    
    if (debug_code != 0)
    { 
      console.log("Previous iter ref points:", remote_start.left, xstart, remote_start.right);
      console.log("logfp, logfs (at previous mode)", logfp, logfs);
    }
    
    
    // Reset starting points for remote limits if last limits were equal to the domain limit (in which case we'd be stuck at that value)
    if (remote_start.left  == this.domain[0]) remote_start.left  = -Infinity;
    if (remote_start.right == this.domain[1]) remote_start.right =  Infinity;
  
    
      // Protect against non-convergence
    
      if (isFinite(remote_start.left))
      {
        let logfp = this.logf_prime(remote_start.left);
        if (isNaN(logfp)) remote_start.left = -Infinity;
      }
    
      if (isFinite(remote_start.right))
      {
        let logfp = this.logf_prime(remote_start.right);
        if (isNaN(logfp)) remote_start.right = Infinity;
      }
      
            
        
    if (logfp > 0)
    { 
         if (!isFinite(remote_start.left))   remote_start.left  = x;
    }
    else if (!isFinite(remote_start.right))  remote_start.right = x;
  
    
    if (debug_remote) console.log("1) Remote starting points", remote_start);
    
    if (debug_mode)
    {
      console.log("Previous mode = ", x); 
      console.log("logfs = ", logfs);
    }
    
    
    if (!found_mode_on_domain_border && logfs > 0)
    { 
      if (debug_mode) console.log("Je cherche un point initial plus safe pour le mode, je pars avec x, logfp, logfs", x, logfp, logfs);
      let tmp = this.mode_SafeStartingPoint(x, logfp, logfs, debug_code == 3);
      
      if (tmp.error) return {error: true};
      else
      { 
        x     = tmp.x;
        logfp = tmp.logfp;
        logfs = this.logf_second(x);
      }
      
      if (debug_mode) console.log("et je reviens avec x = ", x);
    }
  
  
    // Find mode
    
    cont = !found_mode_on_domain_border;
    iter = 0;
    bound = this.domain.slice();
    if (debug_mode || debug_icdf) console.log("Starting search mode at x =", x);
    
    if (cont)
    {
      logfp = this.logf_prime(x);
      logfs = this.logf_second(x);
    }
      
    
    while (cont) 
    {    
      if (iter++ == 250) return {error: true};
      
      if (logfp > 0) bound[0] = x;
      else           bound[1] = x;
      
      change = - logfp / Math.abs(logfs);    
      previous_x = x;
      x -= change;
      if (debug_mode) console.log("x, logfp, logfs, change", previous_x, logfp, logfs, change);
      
      if      (x > bound[1]) x = (x + change + bound[1]) / 2, change = previous_x - x, logfp = this.logf_prime(x);
      else if (x < bound[0]) x = (x + change + bound[0]) / 2, change = previous_x - x, logfp = this.logf_prime(x);
    
      
      cont = Math.abs(change) > epsilon;
      
      if (!cont)
      {
        if      (Math.abs(x - this.domain[0]) < epsilon) x = this.domain[0];
        else if (Math.abs(x - this.domain[1]) < epsilon) x = this.domain[1];
        else if (Math.abs(logfp) > ALT_EPSILON && (bound[1] - bound[0]) > ALT_EPSILON)
        {
          change = - Math.sign(logfp) * epsilon;
          x -= change; 
          cont = true;
        }
      }
      
      if (cont)
      { 
        let tmp = this.ProtectAgainstNaN_logfp(x, previous_x);
        x = tmp.x;
        logfp = tmp.logfp;
        
        logfs = this.logf_second(x);
      }
    }
      
    if (isNaN(x))
    {
      error_log("mode = NaN", this, undefined, xstart, remote_start);
      out.error = true;
      return out;
    }
    
    
    var logfs = this.logf_second(x);
    var target;
    
    out.remote =  {left:  this.domain[0], 
                   right: this.domain[1]}; // to be updated below
  
    
    if (logfs > 0)
    {
      // We have found a mimimum point (between two local modes):
      // we redo the search on both sides of this dip
      
      let tmp = this.higher_local_max(x, logfs, epsilon);
      
      out.mode = tmp.x;
      out.M    = tmp.f;
      
      x = tmp.x;
      logfs = this.logf_second(x);
      
      if (!isFinite(remote_start.left)) remote_start.left = tmp.remote_start[0];
      else remote_start.left = remote_start.left < tmp.remote_start[0] ? remote_start.left : tmp.remote_start[0];
  
      if (!isFinite(remote_start.right)) remote_start.right = tmp.remote_start[1];
      else remote_start.right = remote_start.right > tmp.remote_start[1] ? remote_start.right : tmp.remote_start[1];
    }
    else
    {  
      out.mode = x;
      out.M    = this.logf(out.mode);
      
      let alt_left_start = x + 1/logfs;
      if (alt_left_start < this.domain[0]) alt_left_start = (x + this.domain[0]) / 2;
      
      if      (!isFinite(remote_start.left))        remote_start.left = alt_left_start;
      else if (remote_start.left > alt_left_start)  remote_start.left = alt_left_start;
      
      let alt_right_start = x - 1/logfs;
      if (alt_right_start > this.domain[1]) alt_right_start = (x + this.domain[1]) / 2;
      
      if  (!isFinite(remote_start.right))             remote_start.right = alt_right_start;
      else if (remote_start.right < alt_right_start)  remote_start.right = alt_right_start;
    }
    
    if (debug_code != 0 || debug_mode) console.log("Found mode at x = %f, logf=%f  logfs=%f", x, this.logf(out.mode), this.logf_second(out.mode)); 
    if (debug_remote_left || debug_remote_left_safeStart) console.log("2) Remote starting points", remote_start);
    
    var target = out.M + Math.log(f_ratio_remote);
      
      
    // Find right-side remote point
    
    
    if (out.mode == this.domain[1]) x = this.domain[1], cont = false;
    else cont = true;
    
    if (cont)
    {
      let tmp = this.QuickLookAtBorderAsARemoteEndpoint(target, remote_start, false);
      if (!tmp.cont) x = tmp.x, cont = false;
    }
    
    
    if (cont)
    { 
      bound = [x, this.domain[1]];
    
      if (debug_remote_right_safeStart || debug_remote_right) console.log("Start search for remote right at x ", remote_start.right);
      tmp = this.remoteRef_SafeStartingPoint(remote_start.right, target, false, debug_remote_right_safeStart);
    
      if (tmp.error)
      {
        error_log("remote_SafeStartingPoint/right", this, undefined, xstart, remote_start);
        out.error = true;
        return out;
      }
      else x = tmp.x;
  
      if (debug_remote_right) console.log("or rather at x ", x);
  
      change = out.mode - x; // to step back to mode if log_fp(x) = NaN
      cont = true;
      iter = 0;
    }
    
      
    while (cont) 
    {
      if (iter++ == 250)
      {
        error_log("search for right-remote", this, undefined, xstart, remote_start);
        return {error: true};
      }
    
      let tmp = this.ProtectAgainstNaN_logfp(x, x+change);
        x = tmp.x;
        logfp = tmp.logfp;
      
      
      if (logfp > 0) logfp *= -1; // we are on the wrong side of a local max => force the search on the opposite side
      let g = this.logf(x) - target;
      
      if (g < 0)
      {
        if (x < bound[1]) bound[1] = x;
      }
      else if (x > bound[0]) bound[0] = x;
      
      change = g / logfp;
      x -= change;
      
      
      if (x > bound[1]) 
      {
        // Use a quadratic approx to the curve to find the next point of the N-R search
        x += change; // back to previous x
  
        logfp = this.logf_prime(x);
        logfs = this.logf_second(x);
        
        change = - quadratic_soln(logfs/2, logfp, g).right;
        x -= change;
  
        // in case we are still out of bounds
        if (x > bound[1])
        {
          let previous_x = x + change;
          x = (previous_x + bound[1]) / 2.0;
          change = x - previous_x;
        }      
      }
      else if (x < bound[0])
      {
        let previous_x = x + change;
        x = (previous_x + bound[0]) / 2;
        change = previous_x - x;
      }
  
      
      cont = Math.abs(change) > epsilon;
      
      if (!cont && g > 0.01)
      {
        if (Math.abs(x - this.domain[1]) < epsilon) x = this.domain[1];
        else
        {
          let previous_x = x;
          x = (bound[0] + bound[1]) / 2;
          change = previous_x - x;
          cont = true;
        }
      }
    }
  
    out.remote.right = x;
    if (debug_code != 0) console.log("Found remote right at x ", x);
  
  
    // Find left-side remote point
    
    if (out.mode == this.domain[0]) x = this.domain[0], cont = false;
    else cont = true;
    
    if (cont)
    {
      let tmp = this.QuickLookAtBorderAsARemoteEndpoint(target, remote_start);
      if (!tmp.cont) x = tmp.x, cont = false;
    }
    
    
    if (cont)
    {
      if (debug_remote_left || debug_remote_left_safeStart) console.log("Start search for remote left at x ", remote_start.left);
      tmp = this.remoteRef_SafeStartingPoint(remote_start.left, target, true, debug_remote_left_safeStart);
    
      if (tmp.error)
      {
        error_log("remoteLeft_SafeStartingPoint/left", this, undefined, xstart, remote_start);
        out.error = true;
        return out;
      }
      else x = tmp.x;
    
      if (debug_remote_left || debug_remote_left_safeStart) console.log("or rather at x ", x);
  
      bound = [this.domain[0], out.mode];
      cont = true;
      iter = 0;  
    
      if (debug_remote_left)
      { 
        console.log("Start search for left remote", x);
        console.log("\tlogf=%f\n\ttarget=%f", round(this.logf(x), 3), round(target, 3));
      }
    }
    
      
    if (x != this.domain[0])
    {
      change = out.mode - x; // to step back to mode if log_fp(x) = NaN
    
      while (cont) 
      {
        if (iter++ == 250)
        {
          error_log("search for left-remote", this, iter, xstart, remote_start);
          return {error: true};
        }
  
        
        let tmp = this.ProtectAgainstNaN_logfp(x, x+change);
          x = tmp.x;
          logfp = tmp.logfp;
          
        
        if (debug_remote_left) console.log("x, logfp ==", x, logfp);
        if (logfp < 0) logfp *= -1; // we are on the wrong side of a local max => force the search on the opposite side
        let logf = this.logf(x); 
        let g = logf - target;
        change = g / logfp;
          if (debug_remote_left) console.log("\tlogf=", logf);
      
        if (g < 0)
        {
          if (x > bound[0]) bound[0] = x;
        }
        else if (x < bound[1]) bound[1] = x;
      
        x -= change;
        if (debug_remote_left) console.log("Searching for left remote, x = ", x);
      
      
        if (x < bound[0]) 
        {
          // Use a quadratic approx to the curve to find the next point of the N-R search
          x += change; // back to previous x
  
          logfp = this.logf_prime(x);
          logfs = this.logf_second(x);
          change = - quadratic_soln(logfs/2, logfp, g).left;
          if (isNaN(change)) change = g/logfp;        
          x -= change;
        
          // in case we are still out of bounds
          if (x < bound[0])
          {
            let previous_x = x + change;
            x = (previous_x + bound[0]) / 2.0;
            change = x - previous_x;
          }
        }
        else if (x > bound[1])
        {
          let previous_x = x + change;
          x = (previous_x + bound[1]) / 2.0;
          change = x - previous_x;
        }
  
        converged = Math.abs(change) < epsilon;
        cont = !converged;
  
        if (!cont && g > 0)
        {
          if (Math.abs(x - this.domain[0]) < epsilon) x = this.domain[0];
          else
          {
            let previous_x = x + change;
            
            if (isFinite(bound[0])) x = (bound[0] + bound[1]) / 2;
            else x = previous_x - 1/logfp;
            
            change = previous_x - x;
            cont = true;
          }
        }
      }
    }
  
    out.remote.left = x;
    if (debug_code != 0) console.log("Found remote left at x ", x);
    
    out.error = false;
    return out;
  }, // end of ref_points
  
  
  remoteRef_SafeStartingPoint: function(x_start, target, left_side, debug=false, epsilon=1e-6)
  {
    // Looking for a point within domain with logf_prime not equal to NaN
  
    var cont = true,
        x = x_start,
        iter = 0,
        change = this.mode - x; // to step back to mode if log_fp(x) = NaN
   
        
    var expected_slope_sign = left_side ? 1 : -1;
        
    while (cont)
    {
      if (iter++ == 250) return {error: true};
      
      let tmp = this.ProtectAgainstNaN_logfp(x, x+change);
        x = tmp.x;
        logfp = tmp.logfp;
        change = tmp.change;
          
      logf  = this.logf(x);
  
      
      if (debug) console.log("trying x = %f (logf=%f, logfp=%f)", x, round(logf, 3), round(logfp, 3));
      
      if (logf < target || Math.sign(logfp) == expected_slope_sign) return {x: x, error: false};
      
      // At this point, logfp is not of expected sign
      logfp *= -1;
      
      change = (logf - target) / logfp;
      x -= change;
      if (debug) console.log("x = %f", x);
      
      if (left_side)
      {
        if (x < this.domain[0])
        {
          x += change;
          let previous_x = x;
          
          let logfs = this.logf_second(x);
          let target_slope = (logf - target) / (x - this.domain[0]);
          change = (logfp-target_slope)/logfs;
          x -= change;
          
          if (x < this.domain[0]) x = (x + change + this.domain[0]) / 2;
          if ((x - this.domain[0]) < epsilon) return {x: this.domain[0], error: false};
          change = previous_x - x;
        } 
      }
      else if (x > this.domain[1])
      {
        x += change;
        let previous_x = x;
        
        let logfs = this.logf_second(x);
        let target_slope = (logf - target) / (this.domain[1] - x);
        change = (logfp-target_slope)/logfs;
        change = - Math.abs(change);
        x -= change;
        
        if (debug) console.log("search for right-remote x=%f logfp=%f target_slope=%f, previous_x", x, logfp, target_slope, x+change);
        if (x > this.domain[1]) x = (previous_x + this.domain[1]) / 2;
        if ((this.domain[1] - x) < epsilon) return {x: this.domain[1], error: false};
        change = previous_x - x;
      }
    }
  }, // end of remoteRef_SafeStartingPoint 
  
  
  rgen: function(U, xstart, epsilon=1e-8)
  {
    // The object funs [this] includes:
    //   - The functions:
    //       logf, logf_prime, logf_second
    //
    // xstart: the starting point in the search for the distrn quantile of order U -- it is a good idea to start
    //           at the mode of the density function (for a question of smooth convergence) 
    //           when you don't have an informed guess. 
    //
    //  If range & xstart are left undefined, funs must also include
    //   - mode: a function returning the mode of the density function this.logf
    //   - domain: an array giving the domain endpoints of the random variable (e.g. [0, Infinity], or [0, 1], etc.)
    //  If range IS defined, then the value of logf at its max (mode) must be specified through the value this.M
    //    (as a way of standardizing the area under the curve).
    
    // Example of use:
    //  To sample a value from a beta distrn, one could do the following:
    //    var Beta = new logf_beta(64.32, 15.09);  // where logf_beta is defined in file logf.js
    //
    //    for (let i=0; i<10000; i++) sample.omega[i] = Beta.rgen(); 
    //  One can easily change the Beta distrn parms and go on:
    //    Beta.alpha = 5.4;
    //    Beta.beta  = 8.44;
    //    [etc.]
    
    // Difference between rgen & rgen_icdf
    // ------------------------------------
    // If the reference points (mode & remote endpoints) were not calculated beforehand,
    //   use rgen_icdf, UNLESS they do no need to be calculated:
    //   they may not need to be calculated indeed if funs.domain is finite and if the function mode() is provided
    //   in funs, which is often not possible or not done in complicated distrn fcts)
    
    
    if (typeof U == 'undefined') U = runif(1);
  
    return this.quantile(U, xstart, false, epsilon);                       
  }, // end of Object.rgen
  
  
  rgen_icdf: function(U, xstart, remote_start, debug_code, epsilon=1e-8)
  {
    // The object 'this' is the usual set of functions [funs] describing a distribution function:
    //   - f, logf, logf_prime, logf_second
    //     (and possibly g, see details in logf section of file stats.js)
    // xstart: (a scalar) 
    //         - A starting value in the search for the mode of the distribution;
    //           it will be corrected by the code below if funs.dens_approx is defined and > 0 
    // remote_start: an object with properties (left, right)
    //           --- e.g. last iteration remote endpoints in an MCMC context ---
    //           which will serve as tentative initial values in the search for 
    //           the distrn remote endpoints.
  
  
    var debug_icdf = debug_code == 4 || debug_code < 0;  
    
    // Default values
      
    let ref = this.ref_points(xstart, remote_start, debug_code);
                                                    
    this.M = ref.M;
      this.remote_left  = ref.remote.left;
      this.remote_right = ref.remote.right;
    
    if (debug_code != 0)
    {
      console.log("found mode at x = %f (height=%f)", ref.mode, ref.M);
      console.log("and remote limits at", ref.remote);
    }
    
    var tmp = this.quantile(U, ref.mode, false, epsilon); 
    
    return {x: tmp.x, 
            ref: {mode: ref.mode, remote: ref.remote}
           };
  } // end of rgen_icdf
} // end of variable type declaration DensityFcts


error_log = function(msg, funs, iter, xstart, remote_start, U, range)
{
  console.error("[%s], iter # %d -- Error in %s", funs.label, iter, msg);
  console.error("\tother parms = ", funs.other);
  if (typeof U            != 'undefined') console.error("\tU = ", U);
  if (typeof remote_start != 'undefined') console.error("\tlast ref points", remote_start);
  if (typeof xstart       != 'undefined') console.error("\tstarting point [last/last mode]", xstart);
  if (typeof range        != 'undefined') console.error("\trange", range);
} // end of error_log