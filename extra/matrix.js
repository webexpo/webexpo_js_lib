// Matrix Algebra functions
//
// Author: Patrick Bélisle
//
// Version 0.3 (May 2021)
//  [distributed]                          


// Change log
// ======================
// 
// Version 0.3 (May 2021)
// ----------------------
//   Changed a few fct calls, as they are defined through Array.prototype anymore
//     (e.g. dot_product)
//
//   The definition of the following functions were embedded in MyMatrix:
//     - cbind
//     - diag
//     - inverse
//     - matrix_product & matrix_times_vector [combined under the name: times]
//     - rbind
//     - sym_filled
//     - transpose
//
//  New fct embedded in MyMatrix:
//    - reset_diag
//
//
// Version 0.2 (Apr 2021)
// ----------------------
//   - undocumented changes
//
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code


const MyMatrix = 
{
  m: [],
  
  
  cbind: function(o)
  {
    // o: a MyMatrix object
    
    var cbind = Object.create(MyMatrix);
    
    // copy by value (and not by reference)
    M = this.m.map(function(arr){return arr.slice();});
    for (let r=0; r<o.m.length; r++) M[r] = M[r].concat(o.m[r]);
    
    cbind.m = M;
    return cbind;  
  }, // end of cbind
    
    
  diag: function()
  {
    // Return the diagonal elements of the matrix 'this.m'
      
    var v = [];
    
    var R = this.m.length,
        C = this.m[0].length;
        
    var dim = R < C ? R : C;
    
    for (let i=0; i<dim; i++) v.push(this.m[i][i]);
    
    return v;
  }, // end of diag
  
  
  inverse: function()
  {
    // I use Gaussian Elimination to calculate the inverse of matrix 'm'
    
    // copy by value (and not by reference)
    var M = this.m.map(function(arr){return arr.slice();});
    
    var R = this.m.length;
  
    
    // Add Identity matrix to the right of M
    
    for (let i=0; i<R; i++)
    {
      for (let j=0; j<R; j++)
      {
        M[i].push(i == j ? 1 : 0);
      }
    }
    
    // Clean the lower triangular part of M (make its entries = 0)
    
    for (let i=0; i<R-1; i++)
    {
      M[i] = M[i].map(z => z / M[i][i]); // standardize reference row
      for (let j=i+1; j<R; j++) M[j] = substract(M[j], M[i].map(m => m * M[j][i]));
    }
    
    // Standardize last row
    
    var i = R - 1;
    M[i] = M[i].map(z => z / M[i][i]);
    
    // Clean the upper triangular part of M (make its entries = 0)
    
    for (let i=0; i<R-1; i++)
    {
      for (let j=i+1; j<R; j++) M[i] = substract(M[i], M[j].map(m => m * M[i][j]));
    }
  
    // The inverse is sitting in the right-side part of the augmented matrix M
  
    var inverse = new Array(R);
    
    for (let i=0; i<R; i++) 
    {
      inverse[i] = [];
      for (let j=R; j<2*R; j++) inverse[i].push(M[i][j]);
    }
    
    
    var M = Object.create(MyMatrix);
      M.m = inverse;
      
    return M;
  }, // end of inverse
  

  rbind: function(o)
  {
    // o: a MyMatrix object
    
    var rbind = Object.create(MyMatrix);
    
    // copy by value (and not by reference)
    var M = this.m.map(function(arr){return arr.slice();});
    for (let r=0; r<o.m.length; r++) M.push(o.m[r]);
    
    rbind.m = M;
    return rbind;  
  }, // end of rbind
    
  
  set_diag: function(diag)
  {
    // diag: an array to fill the matrix (this.m) diagonal
    
    var M = Object.create(MyMatrix);
      // copy by value (and not by reference)
      M.m = this.m.map(function(arr){return arr.slice();});
    
    for (let i=0; i<diag.length; i++) M.m[i][i] = diag[i];
    
    return M;
  }, // end of set_diag  
  
  
  sym_filled: function(copyLowerTriangle2Upper = true)
  {
    // Return a copy of 'this.m' as a symetric matrix by copying the values in its lower triangular part to the upper triangular part 
    // (when copyLowerTriangle2Upper is true, or vice-versa when it is false)
    
    var M = Object.create(MyMatrix);
      // copy by value (and not by reference)
      m = this.m.map(function(arr){return arr.slice();});
      M.m = m;
    
    for (let i=1; i<M.m.length; i++)
    {
      if (copyLowerTriangle2Upper)
      {
        for (let j=0; j<i; j++) M.m[j][i] = M.m[i][j];
      }
      else
      {
        for (let j=0; j<i; j++) M.m[i][j] = M.m[j][i];
      }
    }
    
    return M;
  }, // end of sym_filled
  
  
  times: function(b)
  {
    // b: a MyMatrix object or an array
    // Returns - the matrix product of 'this.m' %*% b.m (if b is a matrix)
    //         - the product this.m %*% b               (if b is an array/vector) 
    
    var R = this.m.length; // number of Rows of this.m
    
    
    if (Array.isArray(b))
    {
      let times = [];
      for (let i=0; i<R; i++) times.push(dot_product(this.m[i], b));
        
      return times; 
    }
     
     
    // b is a MyMatrix object
    
    var C = b.m[0].length; // number of Columns of b.m
      var Rb = b.m.length; // number of rows in b
    
    var M = create_matrix(R, C);
  
    for (let j=0; j<C; j++)
    {
      let column_j = [];
      for (let i=0; i < Rb; i++) column_j.push(b.m[i][j]);
      
      for (let i=0; i < R; i++) M.m[i][j] = dot_product(this.m[i], column_j);
    }
    
    return M;
  }, // end of times
 
  
  transpose: function()
  {
    // transpose the matrix 'this.m'
    var R = this.m.length;
    var C = this.m[0].length;
    
    var M = create_matrix(C, R);
    
    for (let i=0; i<R; i++)
    {
      for (let j=0; j<C; j++) M.m[j][i] = this.m[i][j];
    }
    
    return M;
  } // end of transpose
} // end of MyMatrix type definition


create_matrix = function(R, C, new_object=true, fill_in_value=0)
{
  var m = Array(R).fill(fill_in_value).map(x => Array(C).fill(fill_in_value));
  
  if (new_object)
  {
    var M = Object.create(MyMatrix); 
      M.m = m;
  
    return M;
  }
  else return m;
} // end of create_matrix

            
diag = function(arr, new_object=true)
{
  // arr: a (1-D) array
  // => Return a matrix with diagonal elements = 'arr'
  
  var dim = arr.length;
  
  if (new_object)
  {
    var M = create_matrix(dim, dim);
    for (let i=0; i<dim; i++) M.m[i][i] = arr[i];

    return M;
  }
  else
  {
    var m = create_matrix(dim, dim, false);
    for (let i=0; i<dim; i++) m[i][i] = arr[i];

    return m;
  }    
} // end of diag


function LinearRegression_BetaHat(Y, Xp)
{   
  return Xp.times(Xp.transpose()).inverse().times(Xp).times(Y);
} // end of LinearRegression_BetaHat


function NewtonRaphsonChange(Jacobian, f, target=[])
{
  // Jacobian: a MyMatrix object, with property.m a 2-D array of dimension m x m
  // f:        an array of length m
  
  // NOTE: Requires the substract function defined in genericFcts.js
    
  if (target.length > 0) f = substract(f, target);
  
  var change = Jacobian.inverse().times(f);
  
  return change;
} // end of NewtonRaphsonChange


function NewtonRaphson_max_change(theta, previous_theta)
{
  let abs_change = theta.map((t,i) => t - previous_theta[i]).map(Math.abs);
  let max_change = abs_change.shift();
  
  for (i=0; i<abs_change.length; i++)
  {
    if (abs_change[i] > max_change) max_change = abs_change[i];
  }
  
  return max_change;
} // end of NewtonRaphsonMaxChange