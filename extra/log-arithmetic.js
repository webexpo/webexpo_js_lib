// Functions related to calculation with log-objects
// {x: [], s: []}  where x is the log(abs(z)) and s = sign(z)
// PB, Mar 2021


function log_diff(l1, l2)
{
  l1.s = l1.s.times(-1);

  return log_sum(l1, l2);
} // end of log_diff


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


function log_mult(l1, l2)
{
  var l = {x: l1.x.plus(l2.x), s: l1.s.times(l2.s)};
  return l;
} // end of log_mult;


function log_notation(z)
{
  // z: either a number or an array
  
  if (typeof z == 'number')
  {
    out = {x: Math.log(Math.abs(z)), s: Math.sign(z)};
  }
  else
  {
    // z is an array
    out = {x: z.map(Math.abs).map(Math.log), s: z.map(Math.sign)};
  }
  
  return out;
} // end of log_notation


function log_sum(l1, l2)
{
  // l1 & l2 are two objects with dimension (x, s), of same lengths in the two objects
  
  var combo = {l1x: l1.x, l1s: l1.s, 
               l2x: l2.x, l2s: l2.s};

  var same_sign = combo.l1s.equals(combo.l2s);  
  
  var x = [0].rep(l1.x.length), 
      s = [0].rep(l1.x.length);  // Templates for x & s object output dimensions  
  
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
  
  
  var log_sum = {x: x, s: s};
  
  return log_sum;
} // end of log_sum