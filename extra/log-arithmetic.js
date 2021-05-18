// Functions related to calculation with log-objects
// {x: [], s: []}  where x is the log(abs(z)) and s = sign(z)
//
// Author: Patrick Bélisle
//
// Version 0.3 (May 2021)
// [distributed]


// Change log
// ======================
//
// Version 0.3 (May 2021)
// ----------------------
//   Added the logNotation object definition
//     The functions listed below were embedded in logNotation object definition
//     - Object.as_real
//     - Object.lo_div_           [renamed div]
//     - Object.lo_mult           [renamed mult]
//     - Object.lo_mult_          [embedded in the above mult]
//     - Object.lo_opp_sign       [renamed changed_sign]
//     - Object.minus_log         [renamed minus]
//     - Object.plus_log          [renamed plus]
//     - Object.sum_signedExp     [renamed real_sum]
//
//   The function as_log is now calling Object.create(logNotation)
//   The function log_notation was dropped (use as_log instead)
//
//   The type of the following functions was changed from Array.prototype to classical:
//     - as_log
//     - as_real
//
//   The call to fct dot_product was changed as it is not defined through Array.prototype anymore.
//
//
// Version 0.2 (Apr 2021)
// ----------------------
//   Added the following functions:
//     - as_real
//     - minus_log
//     - lo_div_
//     - lo_mult_
//     - lo_opp_sign
//     - plus_log
//     - sum_signedExp 
//
//   Renamed fcts:
//     log_mult -> lo_mult [lo stands for 'log object']
//
//   Dropped function () in favor of:
//     log_diff -> minus_log
//     log_sum -> plus_log
//
//   Added the argument is_already_on_log_scale to function log_notation
//
//   Added the possibility to use an array or a number as 2nd argument to function plus_log
//     (along with options is_already_on_log_scale & positive_sign)
//
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code


// Requires: genericFcts.js


as_log = function(x, is_already_on_log_scale=false, s=1)
{
  // To convert an array with real numbers into a log-notation object
  // (use the function as_real to get back to real numbers)
  //
  // x: an array
  // s: either i) a number (+/- 1)
  //          ii) an array of numbers (+/-1) of same length as 'x'
  //      or iii) a boolean value (true => sign = +, false => sign = -)
  //
  //   NOTE: The argument s is irrelevant when is_already_on_log_scale is false
  
  var o = Object.create(logNotation);
  
  
  if (is_already_on_log_scale)
  {
    if (typeof s == 'boolean') s = s ? 1 : -1;
    if (typeof s == 'number' && Array.isArray(x)) s = rep(s, x.length);
  
    o.x = x;
    o.s = s;
  }
  else if (typeof x == 'number')
  {
    o.x = Math.log(Math.abs(x));
    o.s = Math.sign(x);
  }
  else
  {
    o.x = x.map(Math.abs).map(Math.log);
    o.s = x.map(Math.sign);
  }
  
  return o;
} // end of as_log


as_real = function(x, s, return_sum=true)
{
  // x & s: either a) two arrays of same length
  //            or b) two numbers

  if (typeof x == 'number') return Math.exp(x) * s;
  else
  {
    let e = x.map(Math.exp);
    
    if (return_sum)
    {
      out = 0;
      for (let i=0; i<e.length; i++) out += e[i] * s[i];
    }
    else
    {
      out = [];
      for (let i=0; i<e.length; i++) out.push(e[i] * s[i]);
    }
    
    return out;
  }
} // end of as_real


function log_diff_exp(l1, l2)
{
  // Return log(abs(exp(l2) - exp(l1))) along with its sign , that is: {x, s}
  //
  // l1, l2: two arrays
  
  var ldiff = substract(l2, l1);
  var lde = {x: [], s: ldiff.map(Math.sign)};
  
  
  for (let i=0; i < l1.length; i++)
  {
    if (lde.s[i] > 0)
    {
      let x = l2[i] + Math.log(1 - Math.exp(l1[i]-l2[i]));
      lde.x.push(x);
    }
    else if (lde.s[i] < 0)
    {
      let x = l1[i] + Math.log(1 - Math.exp(l2[i]-l1[i]));
      lde.x.push(x);
    }
    else
    {
      lde.x.push(-Infinity);
    }
  }
  
  return lde;
} // end of log_diff_exp


const logNotation =
{
  x: [],
  s: [],


  as_real: function(return_sum=true)
  {
    // To convert 'this' log-notation object back to real numbers
    return as_real(this.x, this.s, return_sum);
  }, // end of as_real
  
  
  changed_sign: function()
  {
    // this: a log-notation object
    // Returns the object with its sign changed, that is, multiplies the object by -1
    
    if (typeof this.x == 'number') this.s *= -1;
    else                           this.s = this.s.map(s => -s);
    
    return this;
  }, // end of changed_sign
  
  
  div: function(a, pow=1, a_is_already_on_log_scale=true)
  {
    // this: a log-notation object
    //    a: a number or an array (of same length as this.x)
    //       to substract from this.x
    //  pow: scalar (most of the times, an integer, 1 or 2 [for division by square terms])
    //  Substracts 'a' from this.x; in other words, divide the log-object by 'a' 
    //    [or exp(a), if a_is_already_on_log_scale is true]
    
    
    if (typeof this.x == 'number')
    {
      if (!a_is_already_on_log_scale) a = Math.log(a);
      
      this.x -= pow*a;
      return this;
    }
    else if (Array.isArray(a))
    {
      if (!a_is_already_on_log_scale) a = a.map(Math.log);
      
      this.x = this.x.map((x, i) => x - pow*a[i]);
      return this;
    }
    else if (typeof a == 'number')
    {
      if (!a_is_already_on_log_scale) a = Math.log(a);
      
      this.x = this.x.map(x => x - pow*a);
      return this;
    }
  }, // end of div
  

  minus: function(l, l_notLogObject_is_already_on_log_scale=false)
  {
    // this: a log-notation object with dimensions (x, s) of same length
    //    l: either a) a log-notation with same dimensions, of same length,
    //              b) an array, of same length as this.x
    //              c) a number
    //
    // Note that the option l_notLogObject_is_already_on_log_scale is 
    //   irrelevant/ignored when l is a log-notation object
    
    if (typeof l == 'number' || Array.isArray(l)) 
    {
      return this.plus(l, l_notLogObject_is_already_on_log_scale, false);
    }
    else return this.plus(l.changed_sign());
  }, // end of minus
  
  
  mult: function(a, a_is_already_on_log_scale=true)
  {
    // this: a logNotation object
    // a: either a) a logNotation object with properties (x, s) of same length as in 'this'
    //              [in which case a_is_already_on_log_scale is irrelevant/ignored]
    //     or    b) a (positive*) number
    //     or    c) an array of (positive*) numbers
    //
    //   *: in its current version, mult multiplies by positive numbers only, thus avoiding the need
    //      to pay special attention to this.s (not modified by the multiplication by positive numbers) 
    var o = Object.create(logNotation);
    o.s = this.s;
    
  
    if (typeof a == 'number')
    {
      if (!a_is_already_on_log_scale) a = Math.log(a);
    
      if (typeof this.x == 'number')  o.x = this.x + a;
      else                            o.x = this.x.map(x => x + a); // this.x is an array
           
      return o;
    }
    else if (Array.isArray(a))
    {
      if (!a_is_already_on_log_scale) a = a.map(Math.log);
      
      o.x = this.x.map((x, i) => x + a[i]);
      return o;
    }
    else 
    {
      // a is a logNotation object
      o.x = this.x.map((x,i) => x + a.x[i]);
      o.s = this.s.map((s,i) => s * a.s[i]);
      return o;
    }
  }, // end of mult
  
  
  plus: function(l, l_notLogObject_is_already_on_log_scale=false, positive_sign=true)
  {
    // this: a log-notation object with dimensions (x, s) of same length
    //    l: either a) a log-notation with same dimensions, of same length,
    //              b) an array, of same length as this.x
    //              c) a number
    //
    // Note that the options l_notLogObject_is_already_on_log_scale & positive_sign are 
    //   irrelevant/ignored when l is a log-notation object
    
    // If l is an array and l_notLogObject_is_already_on_log_scale=true, then positive_sign should read 'not reverse sign'
    //  => in other words, in that case, if positive_sign = true, then l is ADDED to this,
    //                                      otherwise it is substracted from it
  
  
    if (typeof this.x == 'number')
    {
      if (typeof l == 'number')
      {
        if (l_notLogObject_is_already_on_log_scale) l_x = l;
        else                                        l_x = Math.log(l);
      
        l_s = positive_sign ? 1 : -1;
      }
      else
      {
        l_x = l.x;
        l_s = l.x;
      }
    
      let e_max = this.x > l_x ? this.x : l_x;
      let u = Math.exp(this.x - e_max)*this.s + Math.exp(l_x - e_max)*l_s;
      
      let x = Math.log(Math.abs(u)) + e_max;
      
      return as_log(x, true, Math.sign(u));
    }
    else if (typeof l == 'number')
    {
      if (l_notLogObject_is_already_on_log_scale)
      {
        let s = positive_sign ? 1 : -1;
        s = rep(s, this.x.length);
        l = rep(l, this.x.length);
        
        new_l = {x: l, s: s};
      }
      else
      {       
        let s = Math.sign(l);
        let x = rep(Math.log(Math.abs(l)), this.x.length);
        if (!positive_sign) s *= -1;
        s = rep(s, this.x.length);
                             
        new_l = {x: x, s: s};
      }
      
      return this.plus(new_l); 
    }
    else if (Array.isArray(l))
    {
      if (l_notLogObject_is_already_on_log_scale)
      {
        let s = positive_sign ? 1 : -1;
        s = rep(s, this.x.length);
        
        new_l = {x: l, s: s};
      }
      else
      {       
        let s = l.map(Math.sign);
        let x = l.map(Math.abs).map(Math.log);
        if (!positive_sign) s = s.map(s => -s); 
                            
        new_l = {x: x, s: s};
      }
      
      return this.plus(new_l);
    }
    
    
    var combo = {l1x: this.x, l1s: this.s, 
                 l2x:    l.x, l2s:    l.s};
  
    var same_sign = combo.l1s.map((s,i) => s == combo.l2s[i]);
    
    var x = rep(0, this.x.length), 
        s = rep(0, this.x.length);  // Templates for x & s object output dimensions  
    
    var w = which(same_sign);
    
    if (w.length > 0)
    {
      let tmp = {l1x: [], l2x: [], l2s: []};
      w.forEach(function(w){tmp.l1x.push(combo.l1x[w]); tmp.l2x.push(combo.l2x[w]); tmp.l2s.push(combo.l2s[w]);})
  
      let e_max = pmax(tmp.l1x, tmp.l2x); // array
      
      tmp.l1x = substract(tmp.l1x, e_max); // max element on each row is now 0
      tmp.l2x = substract(tmp.l2x, e_max);
      
      let log_psumExp = log_psum_exp(tmp.l1x, tmp.l2x);
      e_max = add(e_max, log_psumExp);
      
      let which_inf = which_Inf(e_max);
      
      // Prevent NaN's
      if (which_inf.length > 0) which_inf.forEach(function(w){e_max[w] = - Infinity;});
  
      w.forEach(function(w, j){x[w] = e_max[j]; s[w] = tmp.l2s[j];});
    }
    
    
    w = which(same_sign, false);
    
    if (w.length > 0)
    {
      let tmp = {l1x: [], l2x: [], l2s: []};  
      w.forEach(function(w){tmp.l1x.push(combo.l1x[w]); tmp.l2x.push(combo.l2x[w]); tmp.l2s.push(combo.l2s[w]);})
    
      let lde = log_diff_exp(tmp.l1x, tmp.l2x);
      
      w.forEach(function(w, j){x[w] = lde.x[j]; s[w] = tmp.l2s[j] * lde.s[j];});    
    }
    
    return as_log(x, true, s);
  }, // end of plus
    
  
  ratio_sum: function(d, pow=[1,1])
  {
    // d: an array (on log-scale) of same length as this.x
    // => Returns the sum of this.s[i] * exp(this.x[i]*pow[0] - d.x[i]*pow[1]))
    // Note: pow[0] MUST be 1 or 2.
    
    if (pow[0] == 1)
    {
      return as_real(this.x.map((x,i) => x - d[i]*pow[1]), this.s);
    }
    else if (pow[0] == 2)
    {
      return sum(this.x.map((x,i) => x*2 - d[i]*pow[1]).map(Math.exp));
    }
  }, // end of ratio_sum
  
  
  real_sum: function()
  {
    // this: a log-notation object
    // Returns sum(exp(this.x[i])*this.s[i])
  
    return sum_signedExp(this.x, this.s);
  }, // end of real_sum
  
  sq: function()
  {
    var o = Object.create(logNotation);
    o.x = this.x.map(x => 2*x);
    o.s = this.s.map(x => 1);
    
    return o;
  }
} // end of logNotation variable type declaration


sum_signedExp = function(x, s)
{
  // x & s: two arrays of same length 
  // Returns sum(exp(x[i])*s[i])
  // See also fct real_sum in logNotation objects
  
  var e = x.slice();
  var e_max = max(e);
  
  for (let i=0; i<e.length; i++) e[i] -= e_max;
  var low_sum = dot_product(e.map(Math.exp), s);
  
  return low_sum*Math.exp(e_max);
} // end of sum_signedExp