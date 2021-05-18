// Functions of general interest
// Author: Patrick Bélisle


// Version 0.3 (May 2021)
//   [distributed]

// Change log
// =============================================================================
//
// Version 0.3 (May 2021)
// ----------------------
//   The definition of the following functions were moved from Array.prototype
//   to classical fct defns
//     - cumsum
//     - diff
//     - divided_by [renamed: ratio]
//     - division_sqSum [renamed: ratio_sqSum]
//     - division_sum   [renamed: ratio_sum]
//     - dot_product
//     - indexed [renamed: elements]
//     - log_sum_exp [renamed: log_psum_exp]
//     - plus [renamed: add]
//     - rep
//     - round
//     - selected
//     - sorted [renamed: is_sorted]
//     - sqSum
//     - sum_of_log
//     - times [renamed: prod]
//     - which
//     - whichInf [renamed: which_Inf]
//
//   The following fcts Array.prototype.* were dropped:
//     - all  => use built-in .every() instead
//     - any  => use built-in .some()  instead
//     - sign => use .map(Math.sign) instead
//     - sq   => use .map(z => z**2) instead
//
//   New fct: unique
//
//
// Version 0.2 (Apr 2021)
// ----------------------
//   Added the following functions:
//     - any_le
//     - division_sum
//     - division_sqSum
//     - order
//     - rank
//     - round
//     - split
//     - sum_of_log
//
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - Original code


add = function(arr, a)
{
  // arr: an array
  // a: either s scalar or an array of same length as 'arr'

  if (Array.isArray(a))
  {
    // Return the elementwise sum of the two arrays (of same length) this & a
    
    return arr.map(function (x, i){return x + a[i];});
  }
  else
  {
    // Add the scalar 'a' to each element of 'this'
    return arr.map(x => x + a);
  }
} // end of add


all_the_same = function(arr)
{
  if (arr.length > 1)
  {
    let last_elem = arr[arr.length-1];
    return arr.map(z => z == last_elem).every(e => e);
  }
  else if (arr.length == 1)
  {
    return true; // !
  }
  else
  {
    return undefined;
  }
} // end of all_the_same


any_NaN = function(arr)
{
  if (arr.length > 0)
  {
    const bool = (element) => element;
    return arr.map(isNaN).some(bool);
  }
  else
  {
    return undefined;
  }
} // end of any_NaN


any_gt = function(arr, k)
{
  // arr: an array
  // k: scalar comparison value
 
  if (arr.length > 0) return arr.some(x => x > k);
  else                return undefined;
} // end of any_gt


any_le = function(arr, k)
{
  // arr: an array
  // k: scalar comparison value
 
  if (arr.length > 0) return arr.some(x => x <= k);
  else                return undefined;
} // end of any_le


any_lt = function(arr, k)
{
  // arr: an array
  // k: scalar comparison value
 
  if (arr.length > 0) return arr.some(x => x < k);
  else                return undefined;
} // end of any_lt


as_array = function(o, default_value, size)
{
  // o is either an array (in which case the function returns that array untouched)
  //             or a number to be repeated 'size' times
  //   (if o is undefined, then the repeated number is 'default_value')
  
  var a = [];
  
  if (typeof o == 'undefined')
  {
    for (let i = 0; i<size; i++) a.push(default_value);
    return a;
  }
  else if (typeof o == 'number')
  {
    for (let i = 0; i<size; i++) a.push(o); 
    return a; 
  }
  else return o;
} // end of as_array


cumsum = function(arr)
{
  var cumsum = [],
      tot = 0;
  
  for (let i=0; i< arr.length; i++) 
  {
    tot += arr[i];
    cumsum.push(tot);
  }
  
  return cumsum;
} // end of cumsum


diff = function(arr)
{
  // arr: an array
  //
  // Returns the equivalent of diff(a) in R
  // -> as an array if length(arr) >= 3, as a number if length(arr) = 2
  
  var diff = [];
  for (let i=1; i<arr.length; i++) diff.push(arr[i]-arr[i-1]);
  
  if (diff.length == 1) diff = diff[0];
  
  return diff;
} // end of diff


dot_product = function(x, a)
{
  // Returns dot-product (or scalar product) of x & a
  
  return x.reduce(function(r, x, i){return r + x*a[i]}, 0);
} // end of dot_product


elements = function(arr, j)
{
  // arr: array
  // j: array of integers (>= 0)
  // -- was named: indexed
  
  // Return arr[j] (as R would do)
  //   [see also the fct 'selected']
  
  var elements = [];
  
  j.forEach(function(j){elements.push(arr[j]);});
  
  return elements;
} // end of elements


is_sorted = function(arr)
{
  var is_sorted = true;
  
  for (let i=1; i<arr.length && is_sorted; i++)
  {
    if (arr[i] < arr[i-1]) is_sorted = false;
  }
  
  return is_sorted;
} // end of is_sorted


log_abs = function(arr)
{
  // Return log(|arr|) -- an array
  
  return arr.map(x => Math.log(Math.abs(x)));
} // end of log_abs


log_psum_exp = function(a, b)
{
  // a & b: two arrays
  // Return log(exp(a) + exp(b)), where a & b are two arrays of same length
  
  var e = pmax(a, b);
  
  var z1 = a.map((x, i) => x - e[i]).map(Math.exp);
  var z2 = b.map((x, i) => x - e[i]).map(Math.exp);
  
  return z1.map((z, i) => z + z2[i]).map(Math.log).map((z,i) => z + e[i]);
} // end of log_psum_exp


log_sum_exp = function(arr)
{
  // arr: an array
  // Return log(sum(exp(arr[i])))
  
  var e = max(arr);
  arr = arr.map(a => a - e).map(Math.exp);
  
  return e + Math.log(sum(arr));
} // end of log_sum_exp


max = function(arr) 
{
  return Math.max.apply(null, arr);
} // end of max


max_abs = function(arr) 
{
  return Math.max.apply(null, arr.map(Math.abs));
} // end of max_abs


mean = function(arr)
{
  // Return the mean of the elements of 'arr'
  
  if (arr.length > 0)
  {
    return sum(arr) / arr.length;
  }
  else
  {
    return undefined;
  }
} // end of mean


median = function(arr, na_rm=true)
{ 
  var median;
  var a = [];
  
  if (na_rm)
  {
    for (let i=0; i<arr.length; i++)
    {
      if (!isNaN(arr[i]) & typeof(arr[i]) != 'undefined') a.push(arr[i]);
    }
  }
  else
  {
    a = arr.slice();
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
} // end of median


min = function(arr) 
{
  return Math.min.apply(null, arr);
} // end of min


order = function(arr)
{
  var rank = rank(arr);
  var o = rank.slice();
  var taken = [];
  
  for (let i=0; i<arr.length; i++) taken.push(false);
  
  for (let i=0; i<arr.length; i++) 
  {
    let r = rank[i];
    
    while (taken[r]) r--;
    if (r > 0) taken[r] = true;
    
    o[r] = i;
  }

  return o;
} // end of order


pmax = function(arr, a)
{
  // arr & a: two arrays of same length
  
  var pmax = [];
  
  for (let i = 0; i < arr.length; i++)
  {
    let x = arr[i] >= a[i] ? arr[i] : a[i];
    pmax.push(x);
  }
  
  return pmax;
} // end of pmax


prod = function(arr, a)
{
  // arr: an array
  // a: either a scalar or an array of same length as 'arr'

  if (Array.isArray(a))
  {
    // Returns the elementwise product of the two arrays (of same length) 'this' & a
    return arr.map(function (x, i){return x * a[i];});
  }
  else
  {
    // Multiply each element of 'this' by 'a'
    return arr.map(x => x * a);
  }
} // end of prod


quadratic_soln = function(A, B, C)
{
  var delta = B**2 - 4*A*C;
  var midpoint = - B / (2*A);
  var m = Math.sqrt(delta) / (2 * Math.abs(A));
  
  return {left: midpoint - m, right: midpoint + m};
} // end of quadratic_soln


rank = function(arr)
{
  // code found on StockExchange (slightly modified)
  return arr.map((x, i) => [x, i]).sort((a, b) => b[0] - a[0]).reduce((a, x, i, s) => (a[x[1]] = i > 0 && s[i - 1][0] ==  x[0] ? a[s[i - 1][1]] : arr.length - i - 1, a), []);
} // end of rank


ratio = function(arr, a)
{
  // arr: an array
  // a: either a scalar or an array of same length as 'arr'

  if (Array.isArray(a))
  {
    // Return the elementwise division of the two arrays (of same length) arr & a  (arr / a)
     
    return arr.map(function (x, i){return x / a[i];});
  }
  else
  {
    // Divide each element of 'arr' by a
    return arr.map(x => x / a);
  }
} // end of ratio


ratio_sqSum = function(num, denom, log_input=false)
{
  // num, denom: two arrays of same length
  // Returns sum (num_i/denom_i)**2
  
  var sqSum = 0;
  
  if (log_input)
  {
    for (let i=0; i<num.length; i++) sqSum += Math.exp(2*(num[i]-denom[i]));
  }
  else
  {
    for (let i=0; i<num.length; i++) sqSum += (num[i]/denom[i])**2;
  }
  
  return sqSum;
} // end of ratio_sqSum


ratio_sum = function(num, denom, log_input=false)
{
  // num, denom: two arrays of same length
  // Returns sum (num_i/denom_i), that is, the sum of element-wise divisions 
  
  var sum = 0;
  
  if (log_input)
  {
    for (let i=0; i<num.length; i++) sum += Math.exp(num[i]-denom[i]);
  }
  else
  {
    for (let i=0; i<num.length; i++) sum += num[i]/denom[i];
  }
  
  return sum;
} // end of ratio_sum


rep = function(a, n_rep)
{
  // a: an array or a number
  // n_rep: a positive integer or an array of positive integers of same length as 'arr'
  
  var rep = [];


  if (typeof a == 'number')
  {
    for (let i=0; i<n_rep; i++) rep.push(a);
    return rep;
  }
  
  
  // 'a' is an array
  
  if (n_rep.length == 1 || typeof n_rep == 'number')
  {
    let m = Array.isArray(n_rep) ? n_rep[0] : n_rep;
  
    for (let i=0; i<a.length; i++)
    {
      for (let j=0; j<m; j++) rep.push(a[i]);
    }
  }
  else if (n_rep.length == a.length)
  {
    for (let i=0; i<a.length; i++)
    {
      for (let j=0; j<n_rep[i]; j++) rep.push(a[i]);
    }
  }
  
  return rep;
} // end of rep


round = function(x, digits=2)
{
  // x: a number or an array
  
  if (typeof x == 'number') return Number((x).toFixed(digits));
  else
  {
    let round = [];
    for (let i=0; i<x.length; i++) round.push(Number((x[i]).toFixed(digits)));
    return round;
  }
} // end of round


sd = function(arr)
{  
  var sd;
  var n = arr.length;
  
  
  if (n > 1)
  {
    let xbar = sum(arr) / n;
    let sum_xi2 = sqSum(arr);
    
    sd = Math.sqrt((sum_xi2 - n * xbar**2) / (n - 1));
  }
  else
  {
    sd = NaN;
  }
  
  return sd;
} // end of sd


selected = function(arr, cond)
{
  // arr: an array
  // cond: array of boolean (of same length as arr)
  
  // Return arr[cond] (as R would do)
  
  var selected = [];
  
  for (let i=0; i<arr.length; i++)
  {
    if (cond[i]) selected.push(arr[i]);
  }
  
  return selected;
} // end of selected


split = function(arr, ref)
{
  // split values in 'arr' into two arrays, le [<=ref] & gt [>ref]
  var le = [], gt = [];
  
  for (let i=0; i<arr.length; i++)
  {
    if (arr[i] <= ref) le.push(arr[i]);
    else gt.push(arr[i]);
  }
  
  return {le: le, gt: gt};
} // end of split


seq = function(stop, from=0, by=1)
{
  // return seq(from=from, to=stop-1), in R terms
  
  var seq = [];
  for (let i=from; i<stop; i += by) seq.push(i);
  
  return seq;
} // end of seq


sqSum = function(arr)
{
  return arr.reduce(function (a, b) { return a + (b * b); }, 0);
} // end of sqSum


substract = function(arr, m)
{
  // arr: array
  // m: either a scalar or an array of same length as 'arr'

  if (Array.isArray(m))
  {
    // Return the elementwise difference between the two arrays (of same length) arr & m (arr - m)
    
    return arr.map(function (x, i){return x - m[i];});
  }
  else
  {
    // Substract the scalar 'm' from each element of 'arr'
    return arr.map(x => x - m);
  }
} // end of substract


sum = function(arr)
{
  // Return the sum of the elements of 'arr'
  
  if (typeof arr == 'number') return arr;
  else                        return arr.reduce(function (a, b) { return a + b; }, 0); // scalar
} // end of sum


sum_of_log = function(arr)
{
  // Returns the sum of the log of the elements of 'arr', that is:
  // sum(log(arr))
  
  return sum(arr.map(Math.log));
} // end of sum_of_log


unique = function(arr)
{
  // we assume that arr was sorted (in ascending order) beforehand
  
  let u = arr[0];
  let unique = [u];
  
  for (i=1; i<arr.length; i++)
  {
    if (arr[i] > u) u = arr[i], unique.push(u);
  }
  
  return unique;
} // end of unique


which = function(bool_arr, cond=true)
{
  // bool_arr: an array with boolean values 
  var which = [];
  
  for (let i=0; i < bool_arr.length; i++)
  {
    if (bool_arr[i] == cond) which.push(i);
  }

  return which;
} // end of which


which_Inf = function(arr)
{
  return which(arr.map(z => !isFinite(z)));
} // end of which_Inf
