// Matrix Algebra functions
//
// Author: Patrick Bélisle
//
// Version 0.2 (Apr 2021)


// Change log
// ======================
//
// Version 0.1 (Mar 2021)
//   - original code


Array.prototype.cbind = function(m)
{
  // this & m: 2-d arrays
  
  // copy by value (and not by reference)
  var M = this.map(function(arr){return arr.slice();});
  for (let r=0; r<m.length; r++) M[r] = M[r].concat(m[r]);
  
  return M;  
} // end of Array.cbind


create_matrix = function(R, C, initial_value=0)
{
  var M = Array(R).fill(initial_value).map(x => Array(C).fill(initial_value));
  
  return M;
} // end of create_matrix


Array.prototype.diag = function()
{
  // this: either a
  //   a) 2-D array => Return the diagonal elements of the matrix 'this'
  //   b) 1-D array => Return a matrix with diagonal = 'this'
  
  if (Array.isArray(this[0]))
  {
    // it is a matrix (2-D array)
    
    var v = [];
    
    var R = this.length,
        C = this[0].length;
        
    var dim = R < C ? R : C;
    
    for (let i=0; i<dim; i++) v.push(this[i][i]);
    return v;
  }
  else
  {
    // it is a vector (array)
    var dim = this.length;
    var M = Array(dim).fill(0).map(x => Array(dim).fill(0)); // create a null dim x dim matrix
    for (let i=0; i<dim; i++) M[i][i] = this[i];
    return M;
  }
} // end of Array.diag


Array.prototype.inverse = function()
{
  // I use Gaussian Elimination to calculate the inverse of matrix 'this'
  
  // copy by value (and not by reference)
  var M = this.map(function(arr){return arr.slice();});
  
  var R = this.length;

  
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
    for (let j=i+1; j<R; j++) M[j] = M[j].minus(M[i].times(M[j][i]));
  }
  
  // Standardize last row
  
  var i = R - 1;
  M[i] = M[i].map(z => z / M[i][i]);
  
  // Clean the upper triangular part of M (make its entries = 0)
  
  for (let i=0; i<R-1; i++)
  {
    for (let j=i+1; j<R; j++) M[i] = M[i].minus(M[j].times(M[i][j]));
  }

  // The inverse is sitting in the right-side part of the augmented matrix M

  var inverse = new Array(R);
  
  for (let i=0; i<R; i++) 
  {
    inverse[i] = [];
    for (let j=R; j<2*R; j++) inverse[i].push(M[i][j]);
  }
  
  return inverse;
} // end of Array.inverse


function LinearRegression_BetaHat(Y, Xp)
{   
  return Xp.matrix_product(Xp.transpose()).inverse().matrix_product(Xp).matrix_times_vector(Y);
} // end of LinearRegression_BetaHat


Array.prototype.matrix_product = function(b)
{
  // Returns the matrix product of 'this' x b
  
  var R = this.length; // number of Rows of this
  var C = b[0].length; // number of Columns of b
  
  var m = create_matrix(R, C);

  for (let j=0; j<C; j++)
  {
    let column_j = [];
    for (let i=0; i<b.length; i++) column_j.push(b[i][j]);
    
    for (let i=0; i<R; i++) m[i][j] = column_j.dot_product(this[i]);
  }
  
  return m;
} // end of Array.matrix_product


Array.prototype.matrix_times_vector = function(v)
{
  var res = [];
  for (let i=0; i<this.length; i++) res.push(this[i].dot_product(v));
  
  return res; 
} // end of Array.matrix_times_vector


function NewtonRaphsonChange(Jacobian, f, target=[])
{
  // Jacobian: a 2-D array, of dimension m x m
  // f:        an array of length m
  
  // NOTE: Requires the dot_product & minus functions defined in genericFcts.js
  
  var change = [];
  var J_inv = Jacobian.inverse();
  
  if (target.length > 0) f = f.minus(target);
  
  for (let i=0; i<f.length; i++) change.push(J_inv[i].dot_product(f));
  
  return change;
} // end of NewtonRaphsonChange


Array.prototype.rbind = function(m)
{
  // this & m: 2-d arrays
  
  // copy by value (and not by reference)
  var M = this.map(function(arr){return arr.slice();});
  for (let r=0; r<m.length; r++) M.push(m[r]);
  
  return M;  
} // end of Array.rbind


Array.prototype.set_diag = function(diag)
{
  // this: a 2-D array
  // diag: an array to fill the matrix (this) diagonal
  
  // copy by value (and not by reference)
  var M = this.map(function(arr){return arr.slice();});
  
  for (let i=0; i<diag.length; i++)
  {
    M[i][i] = diag[i]; // set diagonal element to new value
  }
  
  return M;
} // end of Array.set_diag


Array.prototype.sym_filled = function(copyLowerTriangle2Upper = true)
{
  // this: a 2-D array (assumed symetric)

  // Return a copy of 'this' as a symetric matrix by copying the values in its lower triangular part to the upper triangular part 
  // (when copyLowerTriangle2Upper is true, or vice-versa when it is false)
  
  // copy by value (and not by reference)
  var M = this.map(function(arr){return arr.slice();});
  
  for (let i=1; i<M.length; i++)
  {
    if (copyLowerTriangle2Upper)
    {
      for (let j=0; j<i; j++) M[j][i] = M[i][j];
    }
    else
    {
      for (let j=0; j<i; j++) M[i][j] = M[j][i];
    }
  }
  
  return M;
} // end of Array.set_diag


Array.prototype.transpose = function()
{
  // this: a 2-D array
  
  var R = this.length;
  var C = this[0].length;
  
  var M = Array(C).fill(0).map(x => Array(R).fill(0));
  
  for (let i=0; i<R; i++)
  {
    for (let j=0; j<C; j++) M[j][i] = this[i][j];
  }
  
  return M;
} // end of Array.transpose