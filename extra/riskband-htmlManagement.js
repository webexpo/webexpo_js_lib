// Code for Riskband/Banerjee algorithm
// Html-related functions
//
// Author: Patrick B�lisle 
//
//
// Version 0.2 (May 2021)
// [distributed]

// Requires: output.js

// Change Log
// ========================
//
// Version 0.2 (May 2021)
// ----------------------
//  Added an error message when in presence of negative values and assuming log-normal distribution of the data.


function change_R()
{
  // user changes the number of regions
  DisplayRegionProbsTable();
} // end of change_R


function clear_data()
{  
  document.getElementById("data").value = "";
  document.getElementById("RCodeWithResults").innerHTML = "";
} // end of clear_data


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


function DisplayRegionProbsTable(A, region_probs_value)
{  
  var R = document.getElementById("R").value;
  var user_defined_region_probs = region_probs_value == 'user';
  
  // var default_prob = 1/R;
  var default_prob; // leave undefined
  
  var defaultA5 = [0.01, 0.1, 0.5, 1];

  if (typeof A === 'undefined')
  {
    if (R == 5)
    {
      A = defaultA5;
    }
    else
    {
      A = [];
      for (let i=0;i<R-1; i++) A.push("");
    }
  }
  
  var new_table = "<table><tr><td style='text-align: center;'>Cut-offs</td>";
  var udProbs = false
  if (user_defined_region_probs) new_table += "<td style='text-align: center;'>Cut-off-delimited region prior probabilities</td>";
  
  new_table += "</tr>";
  
  for (let i=0; i<R; i++)
  {
    if (user_defined_region_probs) 
    {
      new_table += "<tr><td></td><td style='text-align: center;'><input style='text-align: center;' id='rpp" + i + "' class='ud-probs' type='text'";
      if (typeof default_prob !== 'undefined') new_table += " value='" + default_prob + "'";
      new_table += "></td></tr>";
    }
    
    if (i < R-1)
    {
      udProbs = true
      new_table += "<tr><td style='text-align: center;'><input style='text-align: center;' id='A" + i + "' class='cut-offs' type='text' value='" + A[i] + "'";
      new_table += "></td></tr>";
    }  
  }
    
  new_table += "</table>";
  
  document.getElementById("regions_table").innerHTML = new_table; 
  
  $('#regions_table [data-i18n]').i18n()
  
  if ( zygotine.SEG.dataEntries.riskbandCutOffs != undefined ) {
    zygotine.SEG.dataEntries.riskbandCutOffs.reset()
  }
  
  if ( zygotine.SEG.dataEntries.riskbandUserDefnProbs != undefined && udProbs ) {
    zygotine.SEG.dataEntries.riskbandUserDefnProbs.reset()
  }
} // end of DisplayRegionProbsTable


function ErrorMsg(error_code, A)
{
  var msg;
  
  if (error_code == 0)      msg = "Prior probabilities must sum to 1.";
  else if (error_code == 1) msg = "Prior probabilities must all be non-negative.";
  else if (error_code == 2)
  {
    if (A.length == 1) msg = "The cut-off value "  + A[0]         + " does";
    else               msg = "The cut-off values " + A.join(", ") + " do";
    
    msg += " not cross the (mu, sigma) domain.";
  }
  else if (error_code == 3) msg = "Please fill in data in the opposite text box.";
  else if (error_code == 4) msg = "The cut-off values must be listed in ascending order.";
  else if (error_code == 5) msg = "Cut-off values must not be duplicated.";
  else if (error_code == 6) msg = "Log-Normal distribution assumption is incorrect with datasets including (null or) negative values.";
  
   
  /*var btn = "<input type='button' value='Ok' onClick='ClearErrorMsg()'>";
  
  var html = "<table style='background-color: gainsboro;' border=0><tr><td><span style='color: red; font-weight:bold;'>Error</span>&nbsp;&nbsp;" + msg + "</td><td style='text-align: right; padding-left: 20px;'>" + btn + "</td></tr></table>";
  
  document.getElementById("error_msg").innerHTML = html;*/
  return msg;
} // end of ErrorMsg


function PleaseBePatient()
{
  document.getElementById("RCodeWithResults").innerHTML = "<p style='color: DarkSlateBlue;'>Program running.<br>Please be patient.</p>";
} // end of PleaseBePatient


function Read_A_fromHtml(R)
{
  var R = document.getElementById("R").value;
  var A = [];
  
  for (let i=0; i<R-1; i++)
  {
    let html_varname = "A" + i;
    let A_i = document.getElementById(html_varname).value;
    
    A.push(A_i);
  }
  
  A = A.map(Number); // convert to numeric values
  
  return A;
} // end of Read_A_fromHtml


function regionProbsChange(region_probs_option)
{
  var A = Read_A_fromHtml(); // Preserve Cut-off values if already entered
  
  DisplayRegionProbsTable(A, region_probs_option.value);
} // end of regionProbsChange


function WorkComplete(sample)
{
  var html = "<br><br>MCMC results:<br><textarea id='RCode' rows=10 cols=40></textarea>";
  html += "<br><input type='button' id='CopyRCodeWithResults' value='Copy code to clipboard' onClick='CopyRCodeWithMCMCResults()'>";
  document.getElementById("RCodeWithResults").innerHTML = html;
  
  var RCode_mu    = "mcmc = list(mu = c(" + commasep(sample.mu) + ")";
  var RCode_sigma = "sigma = c(" + commasep(sample.sigma) + ")";
  var RCode = RCode_mu + ",\n" + RCode_sigma + ")\n";
  
  document.getElementById("RCode").value = RCode; 
} // end of WorkComplete