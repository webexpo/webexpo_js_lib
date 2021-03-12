// Functions of general interest
// Author: Patrick Bélisle

// Version 0.1 (Mar 2021)


Array.prototype.all = function()
{
  if (this.length > 0)
  {
    const bool = (element) => element;
    return this.every(bool);
  }
  else
  {
    return undefined;
  }
} // end of Array.all


Array.prototype.allTheSame = function()
{
  if (this.length > 0)
  {
    const bool = (element) => element;
    let head = this[0];
    return this.map(z => z == head).every(bool);
  }
  else
  {
    return undefined;
  }
} // end of Array.allTheSame


Array.prototype.any = function()
{
  if (this.length > 0)
  {
    const bool = (element) => element;
    return this.some(bool);
  }
  else
  {
    return undefined;
  }
} // end of Array.any


Array.prototype.any_NaN = function()
{
  if (this.length > 0)
  {
    const bool = (element) => element;
    return this.map(isNaN).some(bool);
  }
  else
  {
    return undefined;
  }
} // end of Array.any_NaN


Array.prototype.any_gt = function(k)
{
  // k: scalar comparison value
 
  if (this.length > 0)
  {
    return this.map(x => x > k).any();
  }
  else
  {
    return undefined;
  }
} // end of Array.any_gt


Array.prototype.any_lt = function(k)
{
  // k: scalar comparison value
 
  if (this.length > 0)
  {
    return this.map(x => x < k).any();
  }
  else
  {
    return undefined;
  }
} // end of Array.any_lt


as_array = function(o, default_value, size)
{
  // o is either an array (in which case the function returns that array untouched)
  //             or a number to be repeated 'size' times
  //   (if o is undefined, then the repeated number is 'default_value')
  
  let a = [];
  
  if (typeof o === 'undefined')
  {
    for (let i = 0; i<size; i++) a.push(default_value);
  }
  else if (typeof o === 'number')
  {
    for (let i = 0; i<size; i++) a.push(o);  
  }
  else
  {
    a = o;
  }
  
  return a;
} // end of as_array


Array.prototype.cumsum = function()
{
  var cumsum = [],
      tot = 0;
  
  for (let i=0; i< this.length; i++) 
  {
    tot += this[i];
    cumsum.push(tot);
  }
  
  return cumsum;
} // end of Array.cumsum


Array.prototype.diff = function()
{
  // Returns the equivalent of diff(a) in R
  // -> as an array if length(a) >= 3, as a number if length(a) = 2
  
  var diff = [];
  for (let i=1; i<this.length; i++) diff.push(this[i]-this[i-1]);
  
  if (diff.length == 1) diff = diff[0];
  
  return diff;
} // end of Array.diff


Array.prototype.divided_by = function(a)
{
  // a: either s scalar or an array of same length as 'this'

  if (Array.isArray(a))
  {
    // Return the elementwise division of the two arrays (of same length) this & a  (this / a)
     
    return this.map(function (x, i){return x / a[i];});
  }
  else
  {
    // Divide each element of 'this' by a
    return this.map(x => x / a);
  }
} // end of Array.divided_by


Array.prototype.dot_product = function(a)
{
  // Returns dot-product (or scalar product) of this & a
  
  return this.reduce(function(r, x, i){return r + x*a[i]}, 0);
} // end of Array.dot_product


Array.prototype.equals = function(a)
{
  // Return the elementwise comparison of the two arrays (of same length) this & a
  // a: array
  
  return this.map(function (x, i){return x == a[i];});
} // end of Array.equals


Array.prototype.indexed = function(j)
{
  // this: array
  // j: array of integers (>= 0)
  
  // Return this[j] (as R would do)
  //   [see also Array.selected]
  
  var indexed = [];
  var z = this.slice(); // hard copy
  
  j.forEach(function(j){indexed.push(z[j]);});
  
  return indexed;
} // end of Array.indexed


Array.prototype.log_abs = function()
{
  // Return log(|this|) -- an array
  
  return this.map(x => Math.log(Math.abs(x)));
} // end of Array.log_abs


Array.prototype.log_sum_exp = function(a)
{
  // Return log(exp(this) + exp(a)), where this & a are two arrays of same length
  
  var z1 = this.map(Math.exp);
  var z2 =    a.map(Math.exp);
  
  z1 = z1.plus(z2);
  
  return z1.map(Math.log);
} // end of Array.log_sum_exp


Array.prototype.max = function() 
{
  return Math.max.apply(null, this);
} // end of Array.max


Array.prototype.max_abs = function() 
{
  return Math.max.apply(null, this.map(Math.abs));
} // end of Array.max_abs


Array.prototype.mean = function()
{
  // Return the mean of the elements of 'this'
  
  if (this.length > 0)
  {
    return this.sum() / this.length;
  }
  else
  {
    return undefined;
  }
} // end of Array.mean


Array.prototype.median = function(na_rm=true)
{ 
  var median;
  var a = [];
  
  if (na_rm)
  {
    for (let i=0; i<this.length; i++)
    {
      if (!isNaN(this[i]) & typeof(this[i]) != 'undefined') a.push(this[i]);
    }
  }
  else
  {
    a = this.slice();
  }
  
  
  a.sort(function(a, b){return a-b});
  var len = a.length;
    
           
  if (len%2 == 1)
  {
    let j = (len - 1) / 2;
    median = a[j];
  }
  else
  {
    let j = len / 2;
    median = (a[j-1] + a[j]) / 2;
  }
          
  return median;
} // end of Array.median


Array.prototype.min = function() 
{
  return Math.min.apply(null, this);
} // end of Array.min


Array.prototype.minus = function(a)
{
  // a: either s scalar or an array of same length as 'this'

  if (Array.isArray(a))
  {
    // Return the elementwise difference between the two arrays (of same length) this & a (this - a)
    
    return this.map(function (x, i){return x - a[i];});
  }
  else
  {
    // Substract the scalar 'a' to each element of 'this'
    return this.map(x => x - a);
  }
} // end of Array.minus


Array.prototype.plus = function(a)
{
  // a: either s scalar or an array of same length as 'this'

  if (Array.isArray(a))
  {
    // Return the elementwise sum of the two arrays (of same length) this & a
    
    return this.map(function (x, i){return x + a[i];});
  }
  else
  {
    // Add the scalar 'a' to each element of 'this'
    return this.map(x => x + a);
  }
} // end of Array.plus


Array.prototype.pmax = function(a)
{
  // a: array
  
  var pmax = [];
  
  for (let i = 0; i < this.length; i++)
  {
    let x = this[i] >= a[i] ? this[i] : a[i];
    pmax.push(x);
  }
  
  return pmax;
} // end of Array.pmax


Array.prototype.rep = function(n_rep)
{
  // this: array
  // n_rep: a positive integer or an array of positive integers of same length as 'this'
  
  var rep = [];

  
  if (n_rep.length == 1 || typeof n_rep == 'number')
  {
    let m = Array.isArray(n_rep) ? n_rep[0] : n_rep;
  
    for (let i=0; i<this.length; i++)
    {
      for (let j=0; j<m; j++) rep.push(this[i]);
    }
  }
  else if (n_rep.length == this.length)
  {
    for (let i=0; i<this.length; i++)
    {
      for (let j=0; j<n_rep[i]; j++) rep.push(this[i]);
    }
  }
  
  return rep;
} // end of Array.rep


Array.prototype.sd = function()
{  
  var sd;
  var n = this.length;
  
  
  if (n > 1)
  {
    let xbar = this.sum() / n;
    let sqSum = this.sqSum();
    
    sd = Math.sqrt((sqSum - n * xbar**2) / (n - 1));
  }
  else
  {
    sd = NaN;
  }
  
  return sd;
} // end of Array.sd


Array.prototype.selected = function(cond)
{
  // this: array
  // cond: array of boolean
  
  // Return this[cond] (as R would do)
  //   [see also Array.indexed]
  
  return this.indexed(cond.which());
} // end of Array.selected


seq = function(stop, from=0)
{
  // return seq(from=from, to=stop-1), in R terms
  
  var seq = [];
  for (let i=from; i<stop; i++) seq.push(i);
  
  return seq;
} // end of seq


Array.prototype.sign = function()
{  
  return this.map(Math.sign);
} // end of Array.sign


Array.prototype.sorted = function()
{
  var sorted = true;
  
  for (let i=1; i<this.length && sorted; i++)
  {
    if (this[i] < this[i-1]) sorted = false;
  }
  
  return sorted;
} // end of Array.sorted


Array.prototype.sq = function () 
{
  return this.map(function (x){return Math.pow(x, 2);});
} // end of Array.sq


Array.prototype.sqSum = function()
{
  return this.reduce(function (a, b) { return a + (b * b); }, 0);
} // end of Array.sqSum


Array.prototype.sum = function()
{
  // Return the sum of the elements of 'this'
  return this.reduce(function (a, b) { return a + b; }, 0); // scalar
} // end of Array.sum


Array.prototype.times = function(a)
{
  // a: either a scalar or an array of same length as 'this'

  if (Array.isArray(a))
  {
    // Return the elementwise product of the two arrays (of same length) this & a
    
    return this.map(function (x, i){return x * a[i];});
  }
  else
  {
    // Multiply each element of 'this' by 'a'
    return this.map(x => x * a);
  }
} // end of Array.times


Array.prototype.which = function(cond=true)
{
  // this: an array with boolean values 
  var which = [];
  
  for (let i = 0; i < this.length; i++)
  {
    if (this[i] == cond) which.push(i);
  }

  return which;
} // end of Array.which


Array.prototype.whichInf = function()
{
  return this.map(z => !isFinite(z)).which();
} // end of Array.whichInf
