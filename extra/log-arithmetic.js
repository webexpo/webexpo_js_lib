// Functions related to calculation with log-objects
// {x: [], s: []}  where x is the log(abs(z)) and s = sign(z)
//
// Author: Patrick Bélisle
//
// Version 0.2 (Apr 2021)


// Change log
// ======================
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


Array.prototype.as_log = function(is_already_on_log_scale=false, positive_sign=true)
{
  // To convert an array with real numbers into a log-notation object
  // (use the function as_real to get back to real numbers)
  
  // The argument positive_sign is irrelevant when is_already_on_log_scale is false

  if (is_already_on_log_scale)
  {
    let sign = positive_sign ? 1 : -1;
  
    return {x: this, s: [sign].rep(this.length)};
  }
  else
  {
    return {x: this.map(Math.abs).map(Math.log),
            s: this.map(Math.sign)};
  }
} // end of Array.as_log


Number.prototype.as_log = function(is_already_on_log_scale=false)
{
  if (is_already_on_log_scale)
  {
    return {x: this, s: 1};
  }
  else
  {
    return {x: Math.log(Math.abs(this)), 
            s: Math.sign(this)};
  }
} // end of Number.as_log


Array.prototype.as_real = function(s)
{
  // see function below
  return this.map(Math.exp).times(s); // array  
} // end of Array.as_real


Object.prototype.as_real = function()
{
  // To convert a log-notation object (obtained with as_log) back to real numbers
  
  if (typeof this.x == 'number') return this.s * Math.exp(this.x);
  else                           return this.x.as_real(this.s);
} // end of Object.as_real


Object.prototype.lo_div_ = function(a, pow=1, a_is_already_on_log_scale=true)
{
  // this: a log-notation object
  //    a: a number or an array (of same length as this.x)
  //       to add to this x
  //  pow: scalar (most of the times, an integer, 1 or 2[for division by square terms])
  //  Adds a to this.x; in other words, divide the log-object by 'a' 
  //    [or exp(a), if a_is_already_on_log_scale is true
  
    // The final underscore in the function name is to mimic that of lo_mult_,
    // where the underscore is to underline its difference with the function lo_mult
  
  
  if (typeof this.x == 'number')
  {
    if (!a_is_already_on_log_scale) a = Math.log(a);
    return {x: this.x - pow*a, s: this.s};
  }
  else if (Array.isArray(a))
  {
    if (!a_is_already_on_log_scale) a = a.map(Math.log);
    
    return {x: this.x.map((x, i) => x - pow*a[i]), s: this.s};
  }
  else if (typeof a == 'number')
  {
    if (!a_is_already_on_log_scale) a = Math.log(a);
    
    return {x: this.x.map(x => x - pow*a), s: this.s};
  }
} // end of Object.lo_div_


function lo_mult(l1, l2)
{
  var l = {x: l1.x.plus(l2.x), s: l1.s.times(l2.s)};
  return l;
} // end of lo_mult;


Object.prototype.lo_mult = function(l)
{
  // this & l: two objects with dimension (x, s), of same length in the two objects
  
  return {x: this.x.plus(l.x), s: this.s.times(l.s)};
} // end of Object.lo_mult


Object.prototype.lo_opp_sign = function()
{
  // this: a log-notation object
  // Returns the object with its sign changed, that is, multiplies the object by -1
  
  if (typeof this.x == 'number') return {x: this.x, s: -this.s};
  else                           return {x: this.x, s: this.s.map(s => -s)};
} // end of Object.lo_opp_sign


Object.prototype.lo_mult_ = function(a, a_is_already_on_log_scale=true)
{
  // this: a log-notation object
  //    a: a number or an array (of same length as this.x)
  //       to add to this.x
  //  Adds 'a' to this.x; in other words, multiply the log-object by 'a' 
  //    [or exp(a), if a_is_already_on_log_scale is true]
  
    // It was not embedded in lo_mult as it does not deal with multiplication by negative numbers,
    //   hence avoiding the need to work on this.s -- (lo_mult can deal by multiplication by negative numbers, but you have to transform your multiplication factors/numbers into a log-object beforehand)
  
  
  if (typeof this.x == 'number')
  {
    if (!a_is_already_on_log_scale) a = Math.log(a);
    return {x: this.x + a, s: this.s};
  }
  else if (Array.isArray(a))
  {
    if (!a_is_already_on_log_scale) a = a.map(Math.log);
    
    return {x: this.x.map((x, i) => x + a[i]), s: this.s};
  }
  else if (typeof a == 'number')
  {
    if (!a_is_already_on_log_scale) a = Math.log(a);
    
    return {x: this.x.map(x => x + a), s: this.s};
  }
} // end of Object.lo_mult_


function log_diff_exp(l1, l2)
{
  // Return log(abs(exp(l2) - exp(l1))) along with its sign , that is: {x, s}
  //
  // l1, l2: two arrays
  
  var ldiff = l2.minus(l1);
  var lde = {x: [], s: ldiff.sign()};
  
  
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


function log_notation(z, is_already_on_log_scale=false)
{
  return z.as_log(is_already_on_log_scale);
} // end of log_notation


Object.prototype.minus_log = function(l, l_notLogObject_is_already_on_log_scale=false)
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
    return this.plus_log(l, l_notLogObject_is_already_on_log_scale, false);
  }
  else return this.plus_log(l.lo_opp_sign());
} // end of Object.minus_log


// ICI deprecated
//Object.prototype.minus_log = function(l, l_notLogObject_is_already_on_log_scale=false)
//{
  // this: a log-notation object with dimensions (x, s) of same length
  //    l: either a) a log-notation with same dimensions, of same length,
  //              b) an array, of same length as this.x
  //              c) a number
  //
  // Note that the option l_notLogObject_is_already_on_log_scale is 
  //   irrelevant/ignored when l is a log-notation object
  
//  if (typeof l == 'number' || Array.isArray(l)) 
//  {
//    return this.plus_log(l, l_notLogObject_is_already_on_log_scale, false);
//  }
  

//  a = Object.assign({}, l); // copy object 'l' by value, not by reference
//  a.s = a.s.map(z => -z);
  
//  return this.plus_log(a);
//} // end of Object.minus_log


Object.prototype.plus_log = function(l, l_notLogObject_is_already_on_log_scale=false, positive_sign=true)
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
    
    return {x: Math.log(Math.abs(u)) + e_max, s: Math.sign(u)};
  }
  else if (typeof l == 'number')
  {
    if (l_notLogObject_is_already_on_log_scale)
    {
      let s = positive_sign ? 1 : -1;
      s = [s].rep(this.x.length);
      l = [l].rep(this.x.length);
      
      new_l = {x: l, s: s};
    }
    else
    {       
      let s = Math.sign(l);
      let x = [Math.log(Math.abs(l))].rep(this.x.length);
      if (!positive_sign) s *= -1;
      s = [s].rep(this.x.length);
                           
      new_l = {x: x, s: s};
    }
    
    return this.plus_log(new_l); 
  }
  else if (Array.isArray(l))
  {
    if (l_notLogObject_is_already_on_log_scale)
    {
      let s = positive_sign ? 1 : -1;
      s = [s].rep(this.x.length);
      
      new_l = {x: l, s: s};
    }
    else
    {       
      let s = l.map(Math.sign);
      let x = l.map(Math.abs).map(Math.log);
      if (!positive_sign) s = s.map(s => -s); 
                          
      new_l = {x: x, s: s};
    }
    
    return this.plus_log(new_l);
  }
  
  
  var combo = {l1x: this.x, l1s: this.s, 
               l2x:    l.x, l2s:    l.s};

  var same_sign = combo.l1s.equals(combo.l2s);  
  
  var x = [0].rep(this.x.length), 
      s = [0].rep(this.x.length);  // Templates for x & s object output dimensions  
  
  var w = same_sign.which();
  
  if (w.length > 0)
  {
    let tmp = {l1x: [], l2x: [], l2s: []};
    w.forEach(function(w){tmp.l1x.push(combo.l1x[w]); tmp.l2x.push(combo.l2x[w]); tmp.l2s.push(combo.l2s[w]);})

    let e_max = tmp.l1x.pmax(tmp.l2x); // array
    
    tmp.l1x = tmp.l1x.minus(e_max); // max element on each row is now 0
    tmp.l2x = tmp.l2x.minus(e_max);
    
    let log_sum_exp = tmp.l1x.log_sum_exp(tmp.l2x);
    e_max = e_max.plus(log_sum_exp);
    
    let whichInf = e_max.whichInf();
    
    // Prevent NaN's
    if (whichInf.length > 0) whichInf.forEach(function(w){e_max[w] = - Infinity;});

    w.forEach(function(w, j){x[w] = e_max[j]; s[w] = tmp.l2s[j];});
  }
  
  
  w = same_sign.which(false);
  
  if (w.length > 0)
  {
    let tmp = {l1x: [], l2x: [], l2s: []};  
    w.forEach(function(w){tmp.l1x.push(combo.l1x[w]); tmp.l2x.push(combo.l2x[w]); tmp.l2s.push(combo.l2s[w]);})
  
    let lde = log_diff_exp(tmp.l1x, tmp.l2x);
    
    w.forEach(function(w, j){x[w] = lde.x[j]; s[w] = tmp.l2s[j] * lde.s[j];});    
  }
  
  return {x: x, s: s};
} // end of Object.plus_log


Array.prototype.sum_signedExp = function(s)
{
  // this & s: two arrays of same length 
  // Returns sum(exp(this[i])*s[i])
  // See also Object.prototype.sum_signedExp
  
  var e = this.slice();
  var e_max = e.max();
  
  for (let i=0; i<this.length; i++) e[i] -= e_max;
  var low_sum = e.map(Math.exp).dot_product(s);
  
  return low_sum*Math.exp(e_max);
} // end of Array.sum_signedExp


Object.prototype.sum_signedExp = function()
{
  // this: a log-notation object (see the function as_log)
  // Returns sum(exp(this.x[i])*this.s[i])
  // See also Array.prototype.sum_signedExp

  return this.x.sum_signedExp(this.s);
} // end of Object.sum_signedExp