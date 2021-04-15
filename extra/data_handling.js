// Functions for reading input from 
// and preparing output for
// Html Form
//
// Author: Patrick Bélisle

// Version 0.1 (Apr 2021)


Object.prototype.AnyNegativeOrNullValue = function()
{
  // Fct to be used to check whether WebExpo data includes negative values or not
  // (before taking their log, when a log-Normal distrn is assumed for data)
  
  var AnyNegative = false;
  
  if (                     this.y.length > 0)            AnyNegativeValue = this.y.any_le(0);

  if (!AnyNegativeValue && this.lt.length > 0)           AnyNegativeValue = this.lt.any_le(0);  
  if (!AnyNegativeValue && this.interval.gt.length > 0)  AnyNegativeValue = this.interval.gt.any_le(0); 
  if (!AnyNegativeValue && this.gt.length > 0)           AnyNegativeValue = this.gt.any_le(0);
   
  return AnyNegativeValue; 
} // end of Object.AnyNegativeOrNullValue


Object.prototype.TakeLog = function()
{
  // Fct to be used to log-transform typical WebExpo data
  
  this.y = this.y.map(Math.log);
  this.gt = this.gt.map(Math.log);
  this.lt = this.lt.map(Math.log);
  this.interval.gt = this.interval.gt.map(Math.log);
  this.interval.lt = this.interval.lt.map(Math.log);
} // end of Object.TakeLog