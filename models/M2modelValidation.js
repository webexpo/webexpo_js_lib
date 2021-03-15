/* eslint 
    valid-jsdoc: 1
    no-extra-parens: 0
*/
/// <reference path="A.js" />
/// <reference path="NUM.js" />
/// <reference path="MT.js" />
/// <reference path="O.js" />
/// <reference path="M0.js" />
/// <reference path="M2modelIntro.js" />
/// <reference path="U.js" />
/// <reference path="S.js" />


zygotine.M.rangeValidation = {};

zygotine.M.rangeValidation.PastData = function() {
    var logN = true;
    var normal = false;
    this.label = "PastDataSummary parameters";
    this.mean = { logN: [-6.908, 4.605], normal: [50, 120] }; 
    this.sd = { logN: [0.405, 2.303], normal : [0.5, 10]}; 
    this.n = [1, 1000]; 
};


zygotine.M.validation = { mcmc: true, measureList: true, specificParameters: true, pastData: true };

/**
* @function
* @param {zygotine.M.SEGInformedVarModelResult} result
* @param {zygotine.M.SEGInformedVarModel} model
*/
zygotine.M.validateModelParameters = function (result, model) {

    var className = model.className;
    var sp = model.specificParameters;
    var ml = model.measureList;

    var addError = function (msg) {
        result.addError(msg);
    };

    var measureListValidator0 = function () {
        var ml = model.measureList;
        var lu = ml.measureCountByType.uncensored;

        if ((lu < 3) || ((lu / ml.n) < 0.3)) {
            if (lu < 3) {
                addError("A measure list must contain at least 3 uncensored measures.");
            } else {
                addError("Uncensored measures must represent at least 30% of all the measures.");
            }
        }

        var logN = sp.logN;

        var limInf = logN ? 0.0001 : 40;
        var limSup = logN ? 1000 : 140;

        var failed = [];
        var some = false;
        var type = "uncensored";
        if (ml.measureByType[type].some(function (el) { return (el.a < limInf) || (el.a > limSup); })) {
            failed.push(type);
        }

        type = "lessThan";
        if (ml.measureByType[type].some(function (el) { return (el.a < limInf) || (el.a > limSup); })) {
            failed.push(type);
        }

        type = "greaterThan";
        if (ml.measureByType[type].some(function (el) { return (el.a < limInf) || (el.a > limSup); })) {
            failed.push(type);
        }

        type = "interval";
        if (ml.measureByType[type].some(function (el) { return (el.a < limInf) || (el.a > limSup) || (el.b < limInf) || (el.b > limSup); })) {
            failed.push(type);
        }

        if (failed.length > 0) {
            addError(zU.fmt("All measures must be contained in the interval [{0}, {1}]. Check for the following categor{2}: [{3}].", limInf, limSup, failed.length > 1 ? "ies" : "y", failed.join(", ")));
        }

        some = ml.measureByType.interval.some(function (el) { return (el.a >= el.b); });
        if (some) {
            addError("For an interval [a, b] to be well defined, b must be greater than a.");
        }
    };

    var mcmcValidator0 = function () {

        /** @type {zygotine.M.McmcParameters} */
        var mcmc = model.mcmcParameters;

        var
            nIter = mcmc.nIter,
            nBurnin = mcmc.nBurnin,
            monitorBurnin = mcmc.monitorBurnin;

        var checkRange = getRangeValidator("MCMC");
        checkRange("nIter", mcmc.nIter, [5000, 150000]);
        checkRange("nBurnin", mcmc.nBurnin, [0, 5000]);
        if (typeof mcmc.monitorBurnin === 'undefined') {
            monitorBurnin = false;
            result.addInfo("MCMC parameter monitorBurnin set to false.");
        }
    };

    var getRangeValidator = function (parameterSet) {
        var format = parameterSet + ": parameter {0}, whose value is {1}, must be contained in the interval {2}";
        var fn = function (paramName, param, range) {
            if ((param < range[0]) || (param > range[1])) {
                addError(zU.fmt(format, paramName, param, zU.fmtArray(range)));
                return false;
            } else {
                return true;
            }
        };

        return fn;
    };

    var infVarSpecificParamsValidator = function () {
        /** @type {zygotine.M.SEGInformedVarModelParameters} */
        var logN = sp.logN;
        var range, param, paramName;
        var checkRange = getRangeValidator("SEGInformedVarModel specific parameters");
        if (checkRange("muLower", sp.muLower, logN ? [-100, -0.5] : [20, 85])) {

            var min = sp.logN ? Math.log(ml.min) : ml.min;
            if (sp.muLower >= min) {
                addError(zU.fmt("Parameter muLower ({0}) must be less than {1}the minimum value ({2}) of the measure list.", sp.muLower, sp.logN ? "the logarithm of " : "", min));
            }
        }

        if (checkRange("muUpper", sp.muUpper, logN ? [0.5, 100] : [85, 140])) {
            var max = sp.logN ? Math.log(ml.max) : ml.max;
            if (sp.muUpper <= max) {
                addError(zU.fmt("Parameter muUpper ({0}) must be greather than {1}the maximum value (={2}) of the measure list.", sp.muUpper, sp.logN ? "the logarithm of " : "", max));
            }
        }

        checkRange("logSigmaMu", sp.logSigmaMu, logN ? [-0.90, 0.48] : [-0.69, 2.30]);
        checkRange("logSigmaPrec", sp.logSigmaPrec, logN ? [0.40, 6.0] : [0.40, 6.0]);
        checkRange("initMu", sp.initMu, logN ? [-6.908, 4.605] : [50, 120]);
        checkRange("initSigma", sp.initSigma, logN ? [0.405, 2.303] : [0.5, 10]);
    };

    var infVarPastDataValidator = function () {
        /** @type {zygotine.M.PastDataSummary} */
        var pds = model.pastData;
        if (pds.defined) {
            var logN = sp.logN;
            var checkRange = getRangeValidator("PastDataSummary parameters");
            checkRange("mean", pds.mean, logN ? [-6.908, 4.605] : [50, 120]); // valeurs identiques à celles de initMu.
            checkRange("sd", pds.sd, logN ? [0.405, 2.303] : [0.5, 10]); // valeurs identiques à celles de initSigma.
            checkRange("n", pds.n, [1, 1000]); 
        } else {
            result.addInfo("Model does not use past data.");
        }
    };

    var infVar = function () {
        /** @type {zygotine.M.SEGInformedVarModelParameters} */
        var valid = zygotine.M.validation;
        if (valid.measureList) {
            measureListValidator0();
        } else {
            result.addWarning("Measure list validation skipped!");
        }

        if (valid.mcmc) {
            mcmcValidator0();
        } else {
            result.addWarning("Mcmc parameters validation skipped!");
        }

        if (valid.specificParameters) {
            infVarSpecificParamsValidator();
        } else {
            result.addWarning("InformedVar model specific parameters skipped!");
        }

        if (valid.pastData) {
            infVarPastDataValidator();
        } else {
            result.addWarning("Pastdata parameters validation skipped!");
        }
    };
    
    var riskbandSpecificParamsValidator = function () {
      var R = document.riskband_form.R.value;
      var logNormalDistrn = $('#logN').is(':checked')
      var A = Read_A_fromHtml(R)
      var region_prior_prob = [];
    
      var prior_density = {equal_region_probs: document.getElementById("rp_equalwts").checked, 
                           uniform:            document.getElementById("rp_unif").checked};

      if ( A.filter(isNaN).length > 0 ) {
        addError(ErrorMsg(6))
      }
      if (!A.sorted())
      {
        addError(ErrorMsg(4))
      }
      else if (A.any_duplicated_cutoff())
      {
        addError(ErrorMsg(5))
      }


      if (!prior_density.equal_region_probs && !prior_density.uniform)
      {
        // user-defined region probs

        for (let i=0; i<R; i++)
        {
          let html_varname = "rpp" + i
          let region_prior_prob_i = document.getElementById(html_varname).value
          region_prior_prob.push(region_prior_prob_i)
        }

        region_prior_prob = region_prior_prob.map(Number)

        let tot_prob = region_prior_prob.sum()

        if (tot_prob != 1)
        {
          addError(ErrorMsg(0))                         
        }

        let any_negative_value = region_prior_prob.filter(p => p < 0).length > 0

        if (any_negative_value)
        {
          addError(ErrorMsg(1))
        }
      }
    }

    if (model instanceof zygotine.M.SEGInformedVarModel) {
        infVar();
    } else
    if (model instanceof zygotine.M.SEGRiskbandModel) {
      riskbandSpecificParamsValidator()
    }
      
};



