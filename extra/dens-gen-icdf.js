// Functions to find reference points (mode & left- and right-side remote limits)
// of a distribution function
//
// Author: Patrick Bélisle

// Version 0.1 (Mar 2021)


higher_local_max = function(funs, dip_x, range, hs, epsilon)
{
  var change,
      fp,
      local_mode = [];
  
  // Find a starting point (for N-R to work smoothly) on each side of dip_x
  
  // i) on right side
  
  var x = dip_x;
  var fs = hs;
  var cont = true;

  while (cont)
  {
    change = 1/fs;
    x += change;
    if (x > range[1]) x = (x - change + range[1]) / 2;
    
    fp = funs.logf_prime(x);
    fs = funs.logf_second(x);
    cont = fp > 0 && fs > 0;
  }  
  
  var start_right = x,
      hp_right    = fp;
  
  // ii) on left side
  
  x = dip_x;
  fs = hs;
  cont = true;
  
  while (cont)
  {
    change = -1/fs;
    x += change;
    if (x < range[0]) x = (x - change + range[0]) / 2;
    
    fp = funs.logf_prime(x);
    fs = funs.logf_second(x);
    cont = fp < 0 && fs > 0;
  }
    
  
  // Run pure Newton-Raphson to find local max on each side
  
  // i) left side
  
  cont = true;
  hp = fp;
  
  while (cont)
  {
    change = hp / funs.logf_second(x);
    x -= change;
    hp = funs.logf_prime(x);
    cont = Math.abs(hp) > epsilon;
  }
  
  local_mode.push(x);
  
  // ii) right side
  
  cont = true;
  x = start_right;
  hp = hp_right;
  
  while (cont)
  {
    change = hp / funs.logf_second(x);
    x -= change;
    hp = funs.logf_second(x);
    cont = abs(hp) > epsilon;
  }
  
  local_mode.push(x);
  
  
  // Compute height at each local max and pick the higher one
  
  h = local_mode.map(funs.logf);
  var j = h[0] > h[1] ? 0 : 1;  
  var out = {x: local_mode[j], f: h[j]};
  
  return out;  
} // end of higher_local_max


ref_points = function(funs, x_start, x_range) 
{
  // Fonction largement inspiree de zygotine.O.ReferencePoints
  //  (qui traitait les cas avec domaine positif seulement)

  // Find x with log.f.prime > 0
  
  var x = x_start,
      change,
      logfp = funs.logf_prime(x_start),
      left_start;
      
  if (logfp > 0) left_start = x;
  

  // Find mode
  
  cont = true;
  var bound = x_range.slice(); // points between which we know the mode is to be found
  
  while (cont) 
  {    
    if (logfp > 0)
    {
      if (x > bound[0]) bound[0] = x;
    }
    else if (x < bound[1]) bound[1] = x;
    
    change = logfp / funs.logf_second(x);
    x -= change;
    
    if      (x > bound[1]) x = (x + change + bound[1]) / 2;
    else if (x < bound[0]) x = (x + change + bound[0]) / 2;
    
    cont = Math.abs(change) > this.epsilon;
    if (cont) logfp = funs.logf_prime(x);
  }
  
  var logfs = funs.logf_second(x);
  
  if (logfs > 0)
  {
    // We have found a mimimum point (between two local modes):
    // we redo the search on both sides of this dip
    
    let tmp = higher_local_max(funs, x, x_range, logfs, this.epsilon);
    
    this.mode = tmp.x;
    this.M    = tmp.f;
    
    x = tmp.x;
    logfs = funs.logf_second(x);
  }
  else
  {  
    this.mode = x;
    this.M    = funs.logf(this.mode);
  }
  
  if (typeof left_start == 'undefined') left_start = x + 1/logfs;

  // Find right-side remote point
  
  cont = true;
  var target = this.M + Math.log(this.fRatioRemote);
  bound = [x, x_range[1]];
  
  x -= 1 / logfs;
    
  while (cont) 
  {
    let logfp = funs.logf_prime(x);
    if (logfp > 0) logfp *= -1; // we are on the wrong side of a local max => force the search on the opposite side
    let g = funs.logf(x) - target;
    
    if (g < 0)
    {
      if (x < bound[1]) bound[1] = x;
    }
    else if (x > bound[0]) bound[0] = x;
    
    change = g / logfp;
    x -= change;
    
    
    if (x > bound[1]) 
    {
      let previous_x = x + change;
      x = (previous_x + bound[1]) / 2;
      change = x - previous_x;
    }
    else if (x < bound[0])
    {
      let previous_x = x + change;
      x = (previous_x + bound[0]) / 2;
      change = x - previous_x;
    }
    
    cont = Math.abs(change) > this.epsilon;
  }

  this.remote_right = x;


  // Find left-side remote point
  
  x = left_start;

  bound = [x_range[0], this.mode];
  cont = true;
  
    
  while (cont) 
  {
    let logfp = funs.logf_prime(x);
    if (logfp < 0) logfp *= -1; // we are on the wrong side of a local max => force the search on the opposite side
    let g = funs.logf(x) - target;
    change = g / logfp;
    
    if (g < 0)
    {
      if (x > bound[0]) bound[0] = x;
    }
    else if (x < bound[1]) bound[1] = x;
    
    x -= change;
    
    if (x < bound[0]) 
    {
      let previous_x = x + change;
      x = (previous_x + bound[0]) / 2.0;
      change = x - previous_x;
    }
    else if (x > bound[1])
    {
      let previous_x = x + change;
      x = (previous_x + bound[1]) / 2.0;
      change = x - previous_x;
    }

    cont = Math.abs(change) > this.epsilon;
  }

  this.remote_left = x;
  
  return this;
} // end of ref_points


ref_points.prototype = 
{
  fRatioRemote: 1e-8,
  epsilon: 1e-6
} // end of ref_points.prototype