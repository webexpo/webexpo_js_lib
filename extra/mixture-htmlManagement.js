// Mixture algorithm
// Html Management functions
//
// Author: Patrick Bélisle
//
// Version 0.2 (May 2021)
//  [distributed]


// Change Log
// ==================
//
// Version 0.2 (May 2021)
// ----------------------
//   Changed fct calls to commasep (not defined through an Array.prototype statement anymore)
// 
// Version 0.1 (Apr 2021)
// ----------------------
//   Original code.




function ClearData()
{  
  document.getElementById("data").value = "";
  HideResultsTextbox();
  //HideProgressBar();
} // end of ClearData


function ClearErrorMsg()
{
  document.getElementById("error_msg").innerHTML = "";
} // end of ClearErrorMsg


function CopyRCodeWithMCMCResults()
{
  let textarea = document.getElementById("RCode");
  textarea.select();
  document.execCommand("copy");
} // end of CopyRCodeWithMCMCResults


function DisplayOmegaPriorInfo(show=true)
{

  if (show)
  {
    let btn = "<input type='button' value='Ok' onClick='DisplayOmegaPriorInfo(false)'>";
    
    var lcl_str = document.form.omegaPriorLCL.value;
    var ucl_str = document.form.omegaPriorUCL.value;
    
      var lcl = Number(lcl_str);
      var ucl = Number(ucl_str);
  
    html = "<table border=0><tr><td style='vertical-align: top;'>";
    
    if (lcl_str.length == 0 && ucl_str.length == 0)
    {
      html += "A Beta distribution will be used for &omega;.<br>Enter your lower &amp; upper prior credible interval limits and click the above button again to know its parameters.";
    }
    else if (lcl < 0 || ucl > 1 || lcl >= ucl || isNaN(lcl) || isNaN(ucl))
    {
      html += "A Beta distribution will be used for &omega;.<br>Enter valid lower &amp; upper prior credible interval limits and click the above button again to know its parameters.";
    }
    else
    {
      let level = Number(document.form.level.value); 
      let tmp = BetaParmsFromQuantiles(lcl, ucl, level);
      
      if (tmp.converged)
      {
        html += "A Beta distribution with parameters<br><span class='parm'>&alpha; = " + round(tmp.alpha, 3) +"</span> and<br><span class='parm'>&beta; = " + round(tmp.beta, 3) + "</span> will be used for &omega;."; 
      }
      else
      {
        html += "Sorry. The algorithm calculating the Beta distribution parameters did not converge.";
      }
    }
    
    html += "</td><td style='vertical-align: bottom;'>";
    html += btn + "</td></tr></table>";
  }
  else
  {
    html = "";
  }
  
  document.getElementById("omega_prior_info").innerHTML = html;
} // end of DisplayOmegaPriorInfo


function ErrorMsg(error_code)
{
  var msg;
  
  if      (error_code == 0)  msg = "Please enter some data.";
  else if (error_code == 1)  msg = "Both Lower &amp; Upper Credible Interval Limits must be entered."
  else if (error_code == 2)  msg = "Lower Limit must be lower than Upper Limit.";
  else if (error_code == 3)  msg = "Both Lower \&amp; Upper Limits must be between 0 and 1.";
  else if (error_code == 4)  msg = "Log-Normal distribution assumption is incorrect with datasets including (null or) negative values.";
  else if (error_code == 5)  msg = "Problem in Beta Distrn rnd value generation.";
  else if (error_code == 6)  msg = "The algorithm computing the prior Beta distribution parameters for &omega; did not converge.<br>Sorry.";
  
   
  var btn = "<input type='button' value='Ok' onClick='ClearErrorMsg()'>";
  
  var html = "<table style='background-color: gainsboro;' border=0><tr><td><span style='color: red; font-weight:bold;'>Error</span>&nbsp;&nbsp;" + msg + "</td><td style='text-align: right; padding-left: 20px;'>" + btn + "</td></tr></table>";
  
  document.getElementById("error_msg").innerHTML = html;
} // end of ErrorMsg


function HideOmegaPriorCrITextboxes()
{
  document.getElementById("omegaPriorCrI_level").innerHTML  = "";
  document.getElementById("omegaPriorCrI_limits").innerHTML = "";
  document.getElementById("what_is_omega_prior").innerHTML = "";
  document.getElementById("omega_prior_info").innerHTML = "";
  
  return;
} // end of HideOmegaPriorCrITextboxes


function HideResultsTextbox()
{
  document.getElementById("RCodeWithResults").innerHTML = "";
  return;
} // end of HideResultsTextbox


function LoadDataset(ds)
{
  const regex = /\;/g;
  
  if (ds == 1)
  {
    data = '2.298;1.337;0.309;>0.231;>0.243;[0.212,1.477];[0.1,2.22];<-3;<-4;<-7;<-8';
  }
  else if (ds == 2)
  {
    data = '47;<19.9;116;<19.9;103;202;128;167;137;<19.9;210;86.8;110;<19.9;<19.9;130;320;34.1;<19.9;169;<19.9;28.4;<19.9;<19.9;64.8;216;167;48.9;61.6;50.9';
  }
  
  
  document.getElementById("data").value = data.replace(regex, '\n');
} // end of LoadDataset


function RCodeWithMCMCSample(obj, obj_name='mcmc')
{
  var RCode    = obj_name + " = list(";
  
  var RCode_mu    = "mu = c("    + commasep(obj.mu) + ")";
  var RCode_sigma = "sigma = c(" + commasep(obj.sigma) + ")";
  var RCode_omega = "omega = c(" + commasep(obj.omega) + ")";
  
  RCode += RCode_mu;
  RCode +=  ",\n" + RCode_sigma;
  RCode +=  ",\n" + RCode_omega;
  
  RCode += ")\n";
  
  return RCode;
} // end of RCodeWithMCMCSample


function RCodeWithMCMCSample4Regression(obj)
{
  var RCode =  "reg = list(";
  
  var tmp = "mu_mode = c(" + commasep(obj.mu_mode) + ")";
  RCode += tmp;
  
  var tmp = "sigma_mode = c(" + commasep(obj.sigma_mode) + ")";
  RCode += ",\n" + tmp;
  
  tmp = "omega_mode = c(" + commasep(obj.omega_mode) + ")";
  RCode +=  ",\n" + tmp;
  //RCode +=  ",\n" + RCode_omega;
  
  RCode += ")\n";
  
  return RCode;
} // end of RCodeWithMCMCSample4Regression


function ShowOmegaPriorCrITextboxes()
{
  html  = "<select name='level' id='level'>";
    html += "<option value='0.95'>95% Credible Interval";
    html += "<option value='0.90'>90% Credible Interval";
    html += "<option value='0.80'>80% Credible Interval";
  html += "</select>";
  document.getElementById("omegaPriorCrI_level").innerHTML = html;
  
  html  = "Lower Limit: <input style='text-align:center;' type='text' name='omegaPriorLCL' id='omegaPriorLCL' size=8><br>";
  html += "Upper Limit: <input style='text-align:center;' type='text' name='omegaPriorUCL' id='omegaPriorUCL' size=8>";
  document.getElementById("omegaPriorCrI_limits").innerHTML = html;
  
  html = "<input type='button' id='omega_prior' onClick='DisplayOmegaPriorInfo();' style='omega_prior' value='What prior distribution will be used?'>";
  document.getElementById("what_is_omega_prior").innerHTML = html;
  
} // end of ShowOmegaPriorCrITextboxes


function WorkComplete(sample, burnin, mcmc)
{
  var html = "<br><br>MCMC results:<br><textarea id='RCode' rows=10 cols=40></textarea>";
  html += "<br><input type='button' id='CopyRCodeWithResults' value='Copy code to clipboard' onClick='CopyRCodeWithMCMCResults()'>";
  document.getElementById("RCodeWithResults").innerHTML = html;
  
  RCode = RCodeWithMCMCSample(sample);
  
  if (mcmc.monitor_burnin) RCode += "\n" + RCodeWithMCMCSample(burnin, 'burnin');
  
  // RCode += "\n" + RCodeWithMCMCSample4Regression(regression);  
  
  document.getElementById("RCode").value = RCode; 
} // end of WorkComplete