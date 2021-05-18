// Functions for reading input from 
// and preparing output for
// Html Form
//
// Author: Patrick Bélisle

// Version 0.3 (May 2021)
// [distributed]

// Change log
// ======================
//
// Version 0.3 (May 2021)
// ----------------------
//   Modified function(s) 
//     - ReadDataFromHtmlTextArea  -- now calling Object.create
//
//   The definition of the following functions were moved from Array.prototype.* to classical fct defns:
//     - commasep
//
//
// Version 0.2 (Apr 2021)
// ----------------------
//   - included the fct defined in output.js (deleted)
//   - Added the function DisplayRuntime
//
// Version 0.1 (Mar 2021)
// ----------------------
//   - original code (was called input.jx)


commasep = function(a, m_by_row=5)
{  
  // a: an array
  
  var code = "",
      tmp = [];
  
  for (let i=1; i<=a.length; i++)
  {
    tmp.push(a[i-1]);
    
    if (i%m_by_row == 0 | i == a.length)
    {
      let str = tmp.join(", ");
      let eol = i == a.length ? "" : ",\n";
      str += eol;
      code += str;
      tmp = [];
    }
  }
  
  return code;
} // end of commasep


commasep_thousands = function(a)
{
  // Return a large number in a nice fashion (with commas separating thousands)
  
  var parts = a.toString().split(".");

  if (parts[0].length > 4)
  { 
    while (parts[0].match(/\d{4}/)) parts[0] = parts[0].replace(/(\d+)(\d{3})/, '$1,$2');
  }
  
  return parts.join(".");
} // end of commasep_thousands


DisplayLocalTime = function()
{
  var now = new Date();
  console.log(now.toString());
} // end of DisplayLocalTime


DisplayRuntime = function(t0, t1, n_iter)
{
  if (typeof t1 == 'undefined') t1 = performance.now();
  
  var seconds_exact = (t1-t0) / 1000;
  var seconds = Math.floor(seconds_exact);
  
  var str = "Run time: " + seconds.toString() + " seconds";
  
  if (seconds > 60)
  {
    let min = Math.floor(seconds/60);
    seconds -= min * 60;
    str += " [" + min.toString() + " min " + seconds.toString() + " sec]";
  }
  
  str += ".";
  
  if (typeof n_iter != 'undefined')
  {
    str += '\n\t';
    let speed = Math.round(n_iter / seconds_exact);
    str += '(approx. ' + commasep_thousands(speed) + ' iterations per second)';
  }
  
  console.log(str);
} // end of DisplayRuntime


MCMCParms = function(my_document)
{
  var mcmc = {niter: 0, burnin: 0, monitor_burnin: false};

  mcmc.niter  = Number(my_document.niter.value);
  mcmc.burnin = Number(my_document.burnin.value);
  mcmc.monitor_burnin = my_document.monitor_burnin.checked;
      
  return mcmc;
} // end of MCMCParms

  
ReadData = function(my_document)
{
  // Example of call: var data = ReadData(document.riskban_form);
  
  var data_txt = my_document.data.value;  // Name of text-area must be -data-
  var data = ReadDataFromHtmlTextArea(data_txt);
  
  data.N = data.size();
  data.any_censored = data.any_censored();
  
  return data;  
} // end of ReadData


ReadDataFromHtmlTextArea = function(data_txt)
{
  const regex = /\S/g;
  const regex_lt = /</;
  const regex_gt = />/;
  const regex_bracket = /[\[\]]/g;
  const regex_num = /\d/;
  
  var y = [],
      lt = [],
      gt = [],
      interval = {gt: [], lt: []};
  
  var data_lines = data_txt.split('\n');


  for (let i = 0;i < data_lines.length; i++)
  {
    let not_empty = data_lines[i].match(regex);
    
    if (not_empty)
    {
      if (data_lines[i].match(regex_lt))
      {
        let value = data_lines[i].replace("<", "");
        lt.push(Number(value));
      }
      else if (data_lines[i].match(regex_gt))
      {
        let value = data_lines[i].replace(">", "");
        gt.push(Number(value));
      }
      else if (data_lines[i].match(regex_bracket))
      {
        let values = data_lines[i].replace(regex_bracket, "").split(',');
        interval.gt.push(Number(values[0]));
        interval.lt.push(Number(values[1]));
      }
      else if (data_lines[i].match(regex_num))
      {
        y.push(Number(data_lines[i]));
      }
    }
  }
  
  
  var data = Object.create(CensoredData);
    data.y = y;
    data.gt = gt;
    data.interval = interval;
    data.lt = lt;
  
  return data;
} // end of ReadDataFromHtmlTextArea