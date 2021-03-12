// Functions for reading input from Html Form
//
// Author: Patrick B�lisle

// Version 0.2 (Mar 2021)

// Change log
// ======================
//
// Version 0.2 (Mar 2021)
//   - included the fct defined in output.js
// Version 0.1 (Mar 2021)
//   - original code (was called input.jx)



Array.prototype.commasep = function(m_by_row=5)
{  
  var code = "",
      tmp = [];
  
  for (let i=1; i<=this.length; i++)
  {
    tmp.push(this[i-1]);
    
    if (i%m_by_row == 0 | i == this.length)
    {
      let str = tmp.join(", ");
      let eol = i == this.length ? "" : ",\n";
      str += eol;
      code += str;
      tmp = [];
    }
  }
  
  return code;
}

MCMCParms = function(_)
{
  var mcmc = {iter: 0, burnin: 0, monitor_burnin: false};
  
  mcmc.niter  = Number($('#nIter').val())
  mcmc.burnin = Number($('#nBurnin').val())
  mcmc.monitor_burnin = $('#monitorBurnin').is(':checked')
      
  return mcmc;
}

ReadData = function(my_document)
{
  // Example of call: var data = ReadData(document.riskban_form);
  
  var data_txt = $('#obsValues').val()
  var data = ReadDataFromHtmlTextArea(data_txt);
  
  data.N = data.y.length + data.gt.length + data.lt.length + data.interval.gt.length;
  data.any_censored = data.gt.length > 0 | data.lt.length > 0 | data.interval.gt.length > 0;
  
  return data;  
}


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
        let values = data_lines[i].replace(regex_bracket, "").split(/[,-]/);
        interval.gt.push(Number(values[0]));
        interval.lt.push(Number(values[1]));
      }
      else if (data_lines[i].match(regex_num))
      {
        y.push(Number(data_lines[i]));
      }
    }
  }
  
  
  data = {y: y, gt: gt, lt: lt, interval: interval};
  
  return data;
}

TakeLog = function(dat)
{
  dat.y = dat.y.map(Math.log);
  dat.gt = dat.gt.map(Math.log);
  dat.lt = dat.lt.map(Math.log);
  dat.interval.gt = dat.interval.gt.map(Math.log);
  dat.interval.lt = dat.interval.lt.map(Math.log);
  return dat
}