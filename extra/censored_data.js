// Functions for reading input from 
// and preparing output for
// Html Form
//
// Author: Patrick Bélisle

// Version 0.2 (May 2021)
// [distributed]


// Change Log
// =======================
//
// Version 0.2 (May 2021)
// ----------------------
//   In previous version (file was called data_handling.js)
//   the functions AnyNegativeOrNullValue & TakeLog 
//   were defined through Object.prototype declarations, which is not *kosher*...
//
//   Added definition of CensoredData
//
//   Renamed the following functions:
//    - AnyCensored -> any_censored
//    - AnyNegativeOrNullValue -> any_le0


const CensoredData = 
{
  y: [],
  lt: [],
  interval: {gt: [], lt: []},
  gt: [],
  
  
  any_censored: function()
  {
    return this.gt.length > 0 | this.lt.length > 0 | this.interval.gt.length > 0;
  }, // end of any_censored


  any_le0: function()
  {
    // Fct to be used to check whether censored data (typical of WebExpo) includes null or negative values or not
    // (before taking their log, when a log-Normal distrn is assumed for data)
    
    var            any_le0 = this.y.some(z => z <= 0);
  
    if (!any_le0)  any_le0 = this.interval.gt.some(z => z <= 0);
    if (!any_le0)  any_le0 = this.gt.some(z => z <= 0);
    if (!any_le0)  any_le0 = this.lt.some(z => z <= 0); 
     
    return any_le0; 
  }, // end of any_le0
  
  
  size: function()
  {
    return this.y.length + this.gt.length + this.lt.length + this.interval.gt.length;
  }, // end of size
  
  
  TakeLog: function()
  {
    // Fct to be used to log-transform (typical WebExpo) data
    
    this.y = this.y.map(Math.log);
    this.gt = this.gt.map(Math.log);
    this.lt = this.lt.map(Math.log);
    this.interval.gt = this.interval.gt.map(Math.log);
    this.interval.lt = this.interval.lt.map(Math.log);
  } // end of TakeLog
} // end of CensoredData variable type declaration