zygotine.SEG = {};
zygotine.X.common = zygotine.SEG


zygotine.SEG.ready = function () {
    zygotine.SEG.hideNumericalResults();
    zygotine.SEG.setDataEntries();
    zygotine.X.ready()

    $("#resetBtn").click(function () {
        zygotine.SEG.reset();
    });

    $("#calcBtn").click(function () {
      console.log(performance.now(), 'calc')
      zygotine.SEG.hideNumericalResults()
      zygotine.X.showBePatient()
      window.setTimeout(zygotine.SEG.runModel, 20)
    });

    (function () {
        //merci � 'Useless Code'
        //https://stackoverflow.com/questions/21012580/is-it-possible-to-write-data-to-file-using-only-javascript
        var textFile = null,
            makeTextFile = function (text) {
                var data = new Blob([text], { type: 'text/plain' });
                // If we are replacing a previously generated file we need to
                // manually revoke the object URL to avoid memory leaks.
                if (textFile !== null) {
                    window.URL.revokeObjectURL(textFile);
                }

                textFile = window.URL.createObjectURL(data);
                return textFile;
            };

        var
            dnldSampleChains = document.getElementById("dnldSampleChainBtn"),
            dnldBurninChains = document.getElementById("dnldBurninChainBtn"),
            dnldMuTraceplot = document.getElementById("dnldMuTraceplotBtn"),
            dnldSigmaTraceplot = document.getElementById("dnldSigmaTraceplotBtn")


        dnldSampleChains.addEventListener('click', function () {
            var link = document.createElement('a');
            link.setAttribute('download', `chaines_sample_mcmc__${zygotine.SEG.lastModel.result.prngSeed}.csv`);
            link.href = makeTextFile(zygotine.X.concatChains(zygotine.SEG.lastModel.result, 'sample', ';'));
            document.body.appendChild(link);

            // wait for the link to be added to the document
            window.requestAnimationFrame(function () {
                var event = new MouseEvent('click');
                link.dispatchEvent(event);
                document.body.removeChild(link);
            });
        }, false);


        dnldBurninChains.addEventListener('click', function () {
            var link = document.createElement('a');
            link.setAttribute('download', `chaines_burnin_mcmc__${zygotine.SEG.lastModel.result.prngSeed}.csv`);
            link.href = makeTextFile(zygotine.X.concatChains(zygotine.SEG.lastModel.result, 'burnin', ';'));
            document.body.appendChild(link);

            // wait for the link to be added to the document
            window.requestAnimationFrame(function () {
                var event = new MouseEvent('click');
                link.dispatchEvent(event);
                document.body.removeChild(link);
            });
        }, false);
        
        dnldMuTraceplot.addEventListener('click', function() { downloadTraceplot({name: 'mu', symbol: '?'}) }, false)
        dnldSigmaTraceplot.addEventListener('click', function() { downloadTraceplot({name: 'sd', symbol: '?'}) }, false)
    })();
    //$('#select_Chains').change(function () { zygotine.SEG.showChain($('#select_Chains').val()); });
};

zygotine.SEG.lastModel = null;
zygotine.SEG.dataEntries = {};

zygotine.SEG.Model = function () {
    zygotine.M.ErrorLogging.call(this);
    this.result = null;
    this.numericalResult = null;
    var entries = zygotine.SEG.dataEntries;
    this.logN = entries.dstrn.currentValue === 'logN';
    this.modelType = entries.method.currentValue == 'classic' ? (entries["sigmaPrior"].currentValue === 'expostats' ? 'inf' : 'unInf') : 'riskband';
    this.oel = Number(entries.oel.currentValue);
    this.percOfInterest = Number(entries.percOfInterest.currentValue);
    this.fracThreshold = Number(entries.fracThreshold.currentValue);
    this.confidenceLevelForCredibileInterval = Number(entries.confidenceLevelForCredibileInterval.currentValue);
    this.model = entries["sigmaPrior"].currentValue === 'expostats' ? 'inf' : 'unInf';
    this.measureRepartition = '';
    this.elapsedTime = {
        runningModel: NaN,
        numericalResults: NaN
    };
};

zygotine.SEG.Model.prototype = Object.create(zygotine.M.ErrorLogging.prototype);

zygotine.SEG.Model.prototype.className = "SEGModel";

//zygotine.SEG.Model.prototype.showChain = function (muOrSd) {
//    $('#values_Chains').text(this.result.chains[muOrSd].data.join('\r\n'));
//};

zygotine.SEG.Model.prototype.getChainAsString = function (chainId) {
    return this.result.chains[chain].data.join('\r\n');
};

zygotine.SEG.Model.prototype.validate = function () {
    var entries = zygotine.SEG.dataEntries;
    var names = zygotine.SEG.getDataEntryNames(entries["sigmaPrior"].currentValue, entries['withPastData'].currentValue);
    var toutes = names.map(n => entries[n]);
    var fautives = toutes.filter(e => !e.validation.ok);
    if (fautives.length !== 0) {
        this.addError($.i18n('info-manq'));
    }
};

zygotine.SEG.Model.prototype.doCalculation = function () {
    if (this.hasError) {
        //alert('Certaines infos sont incorrectes ou manquantes. Les calculs ne peuvent �tre effectu�s.');
        this.result = null;
    } else {
        var entries = zygotine.SEG.dataEntries;
        var ml = new zygotine.M.MeasureList(entries.obsValues.currentValue);
        var nIter = Number(entries.nIter.currentValue);
        var nBurnin = Number(entries.nBurnin.currentValue);
        var monitorBurnin = entries.monitorBurnin.currentValue;
        var mcmc = new zygotine.M.McmcParameters(nIter, nBurnin, monitorBurnin);
        var initMu = Number(entries.initMu.currentValue);
        var initSigma = Number(entries.initSigma.currentValue);
        var muLower = Number(entries.muLower.currentValue);
        var muUpper = Number(entries.muUpper.currentValue);
        var oel = Number(entries.oel.currentValue);
        
        if (this.modelType === 'inf') {
            let logSigmaMu = Number(entries.logSigmaMu.currentValue);
            let logSigmaPrec = Number(entries.logSigmaPrec.currentValue);
            let modelParameters = zygotine.SEG.createSEGInformedVarModelParameters(this.logN, oel, initMu, initSigma, muLower, muUpper, logSigmaMu, logSigmaPrec);
            var pds = zygotine.M.PastDataSummary.dummyPDS;
            let withPastData = (this.modelType === 'inf') && entries['withPastData'].currentValue;
            if (withPastData) {
                let str = entries.pdMean.currentValue + ";" + entries.pdSd.currentValue + ";" + entries.pdN.currentValue;
                pds = new zygotine.M.PastDataSummary(str);
            }

            mdl = new zygotine.M.SEGInformedVarModel(ml, modelParameters, mcmc, pds);
        } else
        if ( this.modelType === 'unInf' ) {
            let sdRangeInf = Number(entries.sdRangeInf.currentValue);
            let sdRangeSup = Number(entries.sdRangeSup.currentValue);
            let modelParameters = zygotine.SEG.createSEGUninformativeModelParameters(this.logN, oel, initMu, initSigma, muLower, muUpper, sdRangeInf, sdRangeSup);
            mdl = new zygotine.M.SEGUninformativeModel(ml, modelParameters, mcmc); // modif f�vrier 2020: "SEGInformedVarModel" remplac� par "SEGUninformativeModel"
        } else {
          let modelParameters = zygotine.SEG.createSEGRiskbandModelParameters(this.logN)
          mdl = new zygotine.M.SEGRiskbandModel(ml, modelParameters, mcmc)
        }

        let t0 = performance.now();
        this.result = mdl.compute(); // compute de BaseModel accepte un prngSeed
        let t1 = performance.now();
        if (!this.result.hasError) {
            this.elapsedTime.runningModel = Math.round10((t1 - t0) / 1000.0, -4);
        }

        let mtypeUncensored = ml.measureByType.uncensored ? ml.measureByType.uncensored.length : 0;
        let mtypeLessThan = ml.measureByType.uncensored ? ml.measureByType.lessThan.length : 0;
        let mtypeGreaterThan = ml.measureByType.uncensored ? ml.measureByType.greaterThan.length : 0;
        let mtypeInterval = ml.measureByType.uncensored ? ml.measureByType.interval.length : 0;
        this.measureRepartition = zygotine.U.fmt($.i18n('non-cens')+':{0},  <:{1}, >:{2}, ' + $.i18n('intervalle') + ':{3}', mtypeUncensored, mtypeLessThan, mtypeGreaterThan, mtypeInterval);
    }
};

zygotine.SEG.Model.prototype.showResults = function (muSample = undefined, sdSample = undefined) {

    let t0 = performance.now();
    if ( muSample == undefined ) {
      muSample = this.result.chains.muSample.data
    }
    if ( sdSample == undefined ) {
      sdSample = this.result.chains.sdSample.data
    }
    let numRes = zygotine.X.getNumericalResult(this.logN, muSample, sdSample, this.oel, this.confidenceLevelForCredibileInterval, this.fracThreshold, this.percOfInterest);
    zygotine.SEG.displayNumericalResults(numRes);
    let t1 = performance.now();
    this.elapsedTime.numericalResults = Math.round10((t1 - t0) / 1000.0, -4);
    $("#numRes_NumRes").children(".VALUES").children('.numRes_NumRes').text(zygotine.SEG.lastModel.elapsedTime.numericalResults);
};

zygotine.SEG.runModel = function () {
    var model = new zygotine.SEG.Model();
    model.validate();
    if (model.hasError) {
        zygotine.X.hideBePatient();
        window.setTimeout(zygotine.X.alert, 0, model.messages.getErrors(), $.i18n('notez'));
        return;
    }

    model.doCalculation();
    zygotine.X.hideBePatient();

    if (model.result.hasError) {
        zygotine.X.alert(model.result.messages.getErrors(), $.i18n('notez'));
        return;
    }

    zygotine.SEG.lastModel = model;
    var test = model.result.chains.muBurnin.data.length === 0;
    $('option[value="muBurnin"]').prop('disabled', test);
    $('option[value="sdBurnin"]').prop('disabled', test);

    model.showResults();

    //    $('#select_Chains').val('muSample');
    //    model.showChain('muSample');
    zygotine.SEG.lastModel = model;
};

//zygotine.SEG.showChain = function (muOrSd) {
//    if (zygotine.SEG.lastModel === null) {
//        return;
//    } else {
//        zygotine.SEG.lastModel.showChain(muOrSd);
//    }
//};

zygotine.SEG.hideNumericalResults = function () {
    $('#numRes').hide();
    $('#dnldChains').hide();
};

zygotine.SEG.showNumericalResults = function (logN) {
    zygotine.SEG.hideNumericalResults();
    if (logN) {
        $('#numRes').children(".norm").hide();
        $('#numRes').children(".logN").show();
    } else {
        $('#numRes').children(".logN").hide();
        $('#numRes').children(".norm").show();
    }

    $('#numRes').show();
    $('#dnldChains').show();
    $('#dnldBurninChainBtn').hide();
    if (zygotine.SEG.dataEntries.monitorBurnin.currentValue && (parseInt(zygotine.SEG.dataEntries.nBurnin.currentValue) > 0)) {
        $('#dnldBurninChainBtn').show();
    }
};

zygotine.SEG.displayNumericalResults = function (resultObject) {

    var display1Result = zygotine.X.display1NumericalResult;
    zygotine.SEG.hideNumericalResults();
    var t = new Date();
    var element;
    element = $("#" + 'measure' + "_NumRes").children(".VALUES").children('.measure_NumRes').text(zygotine.SEG.lastModel.measureRepartition);
    element = $("#" + 'info' + "_NumRes").children(".VALUES").children('.when_NumRes').text(t.toLocaleDateString() + " " + t.toLocaleTimeString());
    element = $("#" + 'info' + "_NumRes").children(".VALUES").children('.prngRoot_NumRes').text(" (" + zygotine.SEG.lastModel.result.prngSeed.toString() + ")");
    element = $("#" + 'runTime' + "_NumRes").children(".VALUES").children('.runTime_NumRes').text(zygotine.SEG.lastModel.elapsedTime.runningModel);
    var id, result;
    ['.A', '.B', '.C', '.R'].forEach(sel => $(sel).text(''));// on vide tous les champs pertinents ...

    var keys = Object.keys(resultObject);
    for (let iKey = 0; iKey < keys.length; iKey++) {
        id = keys[iKey];
        result = resultObject[id];
        display1Result(id, result);
    }

    zygotine.SEG.showNumericalResults(zygotine.SEG.dataEntries.dstrn.currentValue === 'logN');
};

zygotine.SEG.createSEGInformedVarModelParameters = function (logN, oel, initMu, initSigma, muLower, muUpper, logSigmaMu, logSigmaPrec) {
    var params = new zygotine.M.SEGInformedVarModelParameters(logN, oel);
    params.initMu = initMu;
    params.initSigma = initSigma;
    params.muLower = muLower;
    params.muUpper = muUpper;
    params.logSigmaMu = logSigmaMu;
    params.logSigmaPrec = logSigmaPrec;
    return params;
};

zygotine.SEG.createSEGUninformativeModelParameters = function (logN, oel, initMu, initSigma, muLower, muUpper, sdRangeInf, sdRangeSup) {
    var params = new zygotine.M.SEGUninformativeModelParameters(logN, oel); // modif f�vrier 2020: "SEGInformedVarModelParameters" remplac� par "SEGUninformativeModelParameters"
    params.initMu = initMu;
    params.initSigma = initSigma;
    params.muLower = muLower;
    params.muUpper = muUpper;
    params.sdRange = [sdRangeInf, sdRangeSup];
    return params;
};

zygotine.SEG.createSEGRiskbandModelParameters = function (logN) {
  var params = new zygotine.M.SEGRiskbandModelParameters(logN)
  return params
}

zygotine.SEG.defaultEntryValues = (function () {

    var round = (function (nDigits) {
        var n = -nDigits;
        var f = function (x) {
            if (x === null) {
                return null;
            }

            return Math.round10(x, n);
        };
        return f;
    })(5);


    var _1 = function (val) {
        val = round(val);
        return { logN: { inform: val, unInform: val }, norm: { inform: val, unInform: val } };
    };

    var _4 = function (logN_informVal, logN_unInformVal, norm_informVal, norm_unInformVal) {
        return { logN: { inform: round(logN_informVal), unInform: round(logN_unInformVal) }, norm: { inform: round(norm_informVal), unInform: round(norm_unInformVal) } };
    };

    var logI = new zygotine.M.SEGInformedVarModelParameters(true);
    var logU = new zygotine.M.SEGUninformativeModelParameters(true);
    var normI = new zygotine.M.SEGInformedVarModelParameters(false);
    var normU = new zygotine.M.SEGUninformativeModelParameters(false);
    var dev = {};
    // "obsValues"
    dev.nIter = _1('15000');
    dev.nBurnin = _1('500');
    dev.initMu = _4(logI.initMu, logU.initMu, normI.initMu, normU.initMu);
    dev.initSigma = _4(logI.initSigma, logU.initSd, normI.initSigma, normU.initSd);
    dev.muLower = _4(logI.muLower, logU.muLower, normI.muLower, normU.muLower);
    dev.muUpper = _4(logI.muUpper, logU.muUpper, normI.muUpper, normU.muUpper);
    dev.logSigmaMu = _4(logI.logSigmaMu, null, normI.logSigmaMu, null);
    dev.logSigmaPrec = _4(logI.logSigmaPrec, null, normI.logSigmaPrec, null);
    dev.sdRangeInf = _4(null, logU.sdRange[0], null, normU.sdRange[0]);
    dev.sdRangeSup = _4(null, logU.sdRange[1], null, normU.sdRange[1]);
    return dev;
})();

zygotine.SEG.setDataEntries = function () {

    var entries = zygotine.SEG.dataEntries;
    var valueBased = zygotine.X.ValueBasedDataEntry;
    var integer = true;
    var float = false;
    //zygotine.X.SetDefaultEntryValues();
    var defaults = zygotine.SEG.defaultEntryValues;
    var entry;
    var changeFn;

    entries.obsValues = new valueBased('obsValues', '', "Requis. Voir la documentation quant � la fa�on de pr�senter les observations.");
    entries.obsValues.validate = function () {
        var ml;
        var fmtMeasureListMessages = function () {
            if (!ml.hasError) {
                return "";
            } else {
                let tmp = ml.messages.msgList.filter(z => z.level === 'error').map(y => "<li>" + $('<div>').text(y.msgText).html() + "</li>").join("\r\n");
                return "<div id='obsValuesErr' class='errorMsg'><span>Observations</span>\r\n<ul>\r\n" + tmp + "</ul></div>";
            }
        };

        var tmp = this.currentValue.replace(/\s/g, '');
        if (tmp !== '') {
            var empty = false;

            ml = new zygotine.M.MeasureList(this.currentValue);
            this.hasError = ml.hasError;
            if (this.hasError) {
                this.validation = { ok: false, empty: false, messages: fmtMeasureListMessages(), val: '', ml: '' };
            } else {
                this.validation = { ok: true, empty: false, messages: "(n=" + ml.n + ")", val: this.currentValue, ml: ml.toString() };
            }
        } else {
            let tmp = "<div id='obsValuesErr' class='errorMsg'><span>Observations</span>\r\n<ul><li data-i18n='aucune-obs'></li></ul></div>";
            this.validation = { ok: false, empty: true, messages: tmp, val: '', ml: '' };
        }

        return this.validation;
    };

    changeFn = (function (entry) {
        return function () {
            console.log(performance.now(), 'change');
            entry.currentValue = $('#' + entry.elementId).val();
            entry.validate();
            $('#obsValuesErr').remove();
            $('#obsValueCount').hide();
            if (entry.validation.ok) {
                entry.element.removeClass('invalid');
                $('#obsValueCount').text(entry.validation.messages).show();
            } else {
                entry.element.addClass('invalid');
                $('#inputData').after(entry.validation.messages);
            }
        };
    })(entries.obsValues);

    entries.obsValues.element.change(changeFn);

    //mcmc
    entries.nIter = new valueBased('nIter', defaults.nIter.logN.inform, "Requis. Un nombre entier compris entre 500 et  500000", integer, 500, 500000);
    entries.nBurnin = new valueBased('nBurnin', defaults.nBurnin.logN.inform, "Requis. Un nombre entier compris entre 0 et  50000", integer, 0, 15000);
    entries.initMu = new valueBased("initMu", '', "Requis. Un nombre r�el.", float);

    entries.initSigma = new valueBased("initSigma", '', float);
    //param�tres pour l'interpr�tation des donn�es
    entries.oel = new valueBased("oel", '', "Requis. Un nombre r�el correspondant � la valeur limite d'exposition.", float);
    entries.confidenceLevelForCredibileInterval = new valueBased("confidenceLevelForCredibileInterval", 90, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);
    entries.percOfInterest = new valueBased("percOfInterest", 95, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);
    entries.fracThreshold = new valueBased("fracThreshold", 5, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);

    //prior sur mu
    entries.muLower = new valueBased("muLower", '', 'Requis. Un nombre r�el �tablissant un minimum pour la prior de mu.', float);
    entries.muUpper = new valueBased("muUpper", '', 'Requis. Un nombre r�el �tablissant un maximum pour la prior de mu.', float);
    // sigmaPrior : uniform vs expostats
    entries.logSigmaMu = new valueBased("logSigmaMu", '', "Nombre r�el requis lorsque la prior pour sigma est 'expostats'. Donne la valeur de logSigmaMu.", float);
    entries.logSigmaPrec = new valueBased("logSigmaPrec", '', "Nombre r�el requis lorsque la prior pour sigma est 'expostats'. Donne la valeur de logSigmaPrec.", float);
    entries.sdRangeInf = new valueBased("sdRangeInf", '', float);
    entries.sdRangeSup = new valueBased("sdRangeSup", '', float);
    //past data pour expostats
    entries.pdMean = new valueBased("pdMean", '', "Nombre r�el requis lorsque la prior pour sigma est 'expostats' et que des donn�es externes sont prises en compte.", float);
    entries.pdSd = new valueBased("pdSd", '', "Nombre r�el positif requis lorsque la prior pour sigma est 'expostats' et que des donn�es externes sont prises en compte.", float, Number.MIN_VALUE, Number.MAX_VALUE);
    entries.pdN = new valueBased("pdN", '', "Nombre entier (>1) requis lorsque la prior pour sigma est 'expostats' et que des donn�es externes sont prises en compte.", integer, 2, Number.MAX_SAFE_INTEGER);
    //radio

    zygotine.X.setDataEntries()
    
    var radioBased = zygotine.X.RadioButtonBasedDataEntry;
    entries.sigmaPrior = new radioBased('sigmaPrior', 'spExpostats');
    entry = entries.sigmaPrior;
    changeFn = (function (entry) {
        var fn = function () {
            var oStr = 'input[name=' + entry.name + ']:checked';
            entry.currentValue = $(oStr).val();
            var isExpostats = entry.currentValue === 'expostats';
            if (isExpostats) {
                $('.uniformDep').prop('disabled', true);
                $('.expostatsDep').prop('disabled', false);
                $('#withPastData').trigger('change');
            } else {
                $('.expostatsDep').prop('disabled', true);
                $('.pastDataDep').prop('disabled', true);
                $('.uniformDep').prop('disabled', false);
            }
        };
        return fn;
    })(zygotine.SEG.dataEntries.sigmaPrior);
    entry.element.change(changeFn);

    //    entry.element.change(entry.getChangeFn());    //checkbox

    entries.monitorBurnin = new zygotine.X.checkboxBasedDataEntry('monitorBurnin', false);
    entry = entries.monitorBurnin;
    entry.element.change(entry.getChangeFn());


    entries.withPastData = new zygotine.X.checkboxBasedDataEntry('withPastData', false);
    entry = entries.withPastData;
    changeFn = (function (entry, classe) {
        var fn = function () {
            entry.currentValue = $('#' + entry.elementId).prop('checked');
            $("." + classe).prop('disabled', !entry.currentValue);
        };
        return fn;
    })(zygotine.SEG.dataEntries.withPastData, "pastDataDep");
    entry.element.change(changeFn);

    entries.dstrn = new radioBased('dstrn', 'logN');
    entry = entries.dstrn;
    entry.reset();

    entries.dstrn.element.change(function () {
        var oStr = 'input[name=dstrn]:checked';
        entry.currentValue = $(oStr).val();
        zygotine.SEG.setDefaultsForDistribution(entry.currentValue);
    });

    entries.method = new radioBased('method', 'meth-classic')
    entries.method.element.change(function() {
      zygotine.SEG.disableMethodDiv(this)
      zygotine.SEG.dataEntries.method.currentValue = this.value
    })
    entries.method.reset()
    
    entries.regionsWeightingType = new radioBased('regionProbs', 'rp_equalwts')
    entries.regionsWeightingType.element.change(function() { regionProbsChange(this) })
    entries.regionsWeightingType.reset()
    
    entries.muLowerRiskband = new valueBased("mu_lower", '', 'Requis. Un nombre réel établissant un minimum pour la distrbution a priori de mu.', float)
    entries.muUpperRiskband = new valueBased("mu_upper", '', 'Requis. Un nombre réel établissant un maximum pour la distrbution a priori de mu.', float)
    entries.gsdLowerRiskband = new valueBased("gsd_lower", '', 'Requis. Un nombre réel établissant un minimum pour la distrbution a priori de gsd.', float)
    entries.gsdUpperRiskband = new valueBased("gsd_upper", '', 'Requis. Un nombre réel établissant un maximum pour la distrbution a priori de gsd.', float)
    
  zygotine.SEG.reset()
};

zygotine.SEG.getDataEntryNames = (function () {

    var common = ['obsValues', 'nIter', 'nBurnin', 'monitorBurnin', 'oel', 'confidenceLevelForCredibileInterval', 'fracThreshold', 'percOfInterest', 'initMu', 'initSigma', 'muLower', 'muUpper', 'prngSeed'];
    var unifPrior = ['sdRangeInf', 'sdRangeSup'];
    var expostatsPrior = ['logSigmaMu', 'logSigmaPrec'];
    var pd = ['pdMean', 'pdSd', 'pdN'];

    var fn = function (sigmaPrior, withPastData) {
        var rep = [];

        Array.prototype.push.apply(rep, common);

        if (sigmaPrior === 'uniform') {
            Array.prototype.push.apply(rep, unifPrior);
        } else {
            Array.prototype.push.apply(rep, expostatsPrior);
            if (withPastData) {
                Array.prototype.push.apply(rep, pd);
            }
        }

        return rep;
    };

    return fn;
})();

zygotine.SEG.reset = function () {
    zygotine.X.hideBePatient();
    zygotine.SEG.hideNumericalResults();
    var entries = zygotine.SEG.dataEntries;
    entries.dstrn.reset();
    entries.obsValues.reset();
    entries.nIter.reset();
    entries.nBurnin.reset();
    entries.monitorBurnin.reset();
    entries.withPastData.reset();
    // ++ monitorBurnin ++
    //param�tres pour l'interpr�tation des donn�es
    entries.oel.reset();
    entries.confidenceLevelForCredibileInterval.reset();
    entries.percOfInterest.reset();
    entries.fracThreshold.reset();
    //past data
    entries.pdMean.reset();
    entries.pdSd.reset();
    entries.pdN.reset();
    //radio ... autres que dstrn
    entries.sigmaPrior.reset();
    entries.prngSeed.element.val(zygotine.X.genPseudoRand32Bit())
    $("#RCodeWithResults").html("")
    ClearRiskbandErrorMsg()
    
};

zygotine.SEG.setDefaultsForDistribution = function (loi) {
    var entries = zygotine.SEG.dataEntries;
    var defaults = zygotine.SEG.defaultEntryValues;

    entries.initMu.element.val(defaults.initMu[loi].inform).trigger('change');
    entries.initSigma.element.val(defaults.initSigma[loi].inform).trigger('change');
    //prior sur mu
    entries.muLower.element.val(defaults.muLower[loi].inform).trigger('change');
    entries.muUpper.element.val(defaults.muUpper[loi].inform).trigger('change');
    // sigmaPrior name / uniform vs expostats
    entries.logSigmaMu.element.val(defaults.logSigmaMu[loi].inform).trigger('change');
    entries.logSigmaPrec.element.val(defaults.logSigmaPrec[loi].inform).trigger('change');
    entries.sdRangeInf.element.val(defaults.sdRangeInf[loi].unInform).trigger('change');
    entries.sdRangeSup.element.val(defaults.sdRangeSup[loi].unInform).trigger('change');
    zygotine.X.setDefaultsForDistribution(loi)
};

zygotine.SEG.tests = {
    ecrireLogNExpostatsWoPDS: function () {
        var a = '<0.0923\n0.0595\n<0.222\n[0.138,0.309]\n<0.149\n2.65\n<0.245\n<0.0892\n0.0817\n0.0251\n<0.0639\n[0.815,1.83]\n<0.0748\n0.628\n0.914\n0.949\n0.0377\n<0.153\n0.0996\n0.0499\n0.0471\n0.605\n0.157\n1.28\n[0.153,0.343]\n<0.0706\n>2.48\n0.169\n4.49\n[0.125,0.282]\n<0.206\n<0.109\n0.0869\n<0.123\n0.229\n<0.121\n0.237\n<0.0376\n<0.0773\n0.153\n<0.205\n[0.282,0.633]\n<0.158\n>2.13\n0.347\n3.57\n>1.84\n0.418\n0.442\n<0.22\n<0.236\n>1.77\n3.71\n<0.13\n>2.06\n<0.126\n0.206\n<0.0819\n3.08\n>1.2\n0.138\n<0.177\n0.209\n<0.158\n0.052\n<0.184\n0.0515\n<0.0566\n>3\n0.328\n>1.26\n1.81\n<0.143\n<0.138\n>2.38\n[0.141,0.317]\n<0.239\n<0.0435\n<0.0586\n<0.254\n0.122\n<0.078\n0.105\n<0.0541\n<0.0175\n[0.163,0.367]\n[0.231,0.52]\n0.573\n0.0546\n[0.116,0.26]\n<0.221\n0.164\n<0.144\n[0.485,1.09]\n<0.063\n0.228\n<0.23\n>1.63\n<0.245\n0.42';
        zygotine.SEG.reset();
        $('#logN').prop('checked', true).trigger('change');
        $('#spExpostats').prop('checked', true).trigger('change');
        $('#withPastData').prop('checked', false).trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('.95').trigger('change');
        $('#prngSeed').val('1111').trigger('change')
    },

    ecrireLogNExpostatsWithPDS: function () {
        var a = '<0.0923\n0.0595\n<0.222\n[0.138,0.309]\n<0.149\n2.65\n<0.245\n<0.0892\n0.0817\n0.0251\n<0.0639\n[0.815,1.83]\n<0.0748\n0.628\n0.914\n0.949\n0.0377\n<0.153\n0.0996\n0.0499\n0.0471\n0.605\n0.157\n1.28\n[0.153,0.343]\n<0.0706\n>2.48\n0.169\n4.49\n[0.125,0.282]\n<0.206\n<0.109\n0.0869\n<0.123\n0.229\n<0.121\n0.237\n<0.0376\n<0.0773\n0.153\n<0.205\n[0.282,0.633]\n<0.158\n>2.13\n0.347\n3.57\n>1.84\n0.418\n0.442\n<0.22\n<0.236\n>1.77\n3.71\n<0.13\n>2.06\n<0.126\n0.206\n<0.0819\n3.08\n>1.2\n0.138\n<0.177\n0.209\n<0.158\n0.052\n<0.184\n0.0515\n<0.0566\n>3\n0.328\n>1.26\n1.81\n<0.143\n<0.138\n>2.38\n[0.141,0.317]\n<0.239\n<0.0435\n<0.0586\n<0.254\n0.122\n<0.078\n0.105\n<0.0541\n<0.0175\n[0.163,0.367]\n[0.231,0.52]\n0.573\n0.0546\n[0.116,0.26]\n<0.221\n0.164\n<0.144\n[0.485,1.09]\n<0.063\n0.228\n<0.23\n>1.63\n<0.245\n0.42';
        zygotine.SEG.reset();
        $('#logN').prop('checked', true).trigger('change');
        $('#spExpostats').prop('checked', true).trigger('change');
        $('#withPastData').prop('checked', true).trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('.95').trigger('change');
        $('#pdMean').val('-0.5108256').trigger('change');
        $('#pdSd').val('0.9162907').trigger('change');
        $('#pdN').val('20').trigger('change');
        $('#prngSeed').val('2222').trigger('change')
    },

    ecrireNormExpostatsWoPDS: function () {
        var a = ">93;|[80.5,90.5];|84.2;|<78.7;|[80.7,90.7];|>93.7;|86.6;|87;|85.7;|83.1;|86.9;|84.3;|>93.3;|84.3;|86;|88;|<78.9;|<78.2;|84.4;|83.5;|<78.5;|>92.7;|<76.8;|>94.9;|>92.6;|<76.7;|>93.6;|>92.9;|<79.2;|87.4;|84.3;|<78.5;|<77.1;|88.6;|88.7;|<78.8;|86.8;|<76.8;|87;|<79.4;|84;|[79.9,89.9];|[80.9,90.9];|84.3;|<78.9;|87.9;|<76.8;|88.3;|<75.7;|83.4;|85.8;|88.1;|<77.7;|84;|<79.2;|85.9;|<78.4;|<77.5;|<78.5;|<78.4;|87.2;|[80.9,90.9];|<79.2;|84.2;|<79.4;|<79.3;|[81.8,91.8];|83.5;|<79.3;|85.4;|<79.7;|[81.7,91.7];|<79.3;|<79;|83.9;|[81.2,91.2];|<79.3;|87.9;|[80.5,90.5];|[82.2,92.2];|>93.1;|84.9;|86.2;|<77.4;|87;|<79;|<79;|84.4;|<79.4;|<76.1;|84.7;|<76.7;|<77.2;|83.4;|<79.2;|<78.1;|<79.1;|<78.4;|>92.9;|86.6;";
        a = new zygotine.M.MeasureList(a).toString('\r\n');
        $('#norm').prop('checked', true).trigger('change');
        $('#spExpostats').prop('checked', true).trigger('change');
        $('#withPastData').prop('checked', false);
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('90.0').trigger('change');
        $('#prngSeed').val('3333').trigger('change')
    },

    ecrireNormExpostatsWithPDS: function () {
        var a = ">93;|[80.5,90.5];|84.2;|<78.7;|[80.7,90.7];|>93.7;|86.6;|87;|85.7;|83.1;|86.9;|84.3;|>93.3;|84.3;|86;|88;|<78.9;|<78.2;|84.4;|83.5;|<78.5;|>92.7;|<76.8;|>94.9;|>92.6;|<76.7;|>93.6;|>92.9;|<79.2;|87.4;|84.3;|<78.5;|<77.1;|88.6;|88.7;|<78.8;|86.8;|<76.8;|87;|<79.4;|84;|[79.9,89.9];|[80.9,90.9];|84.3;|<78.9;|87.9;|<76.8;|88.3;|<75.7;|83.4;|85.8;|88.1;|<77.7;|84;|<79.2;|85.9;|<78.4;|<77.5;|<78.5;|<78.4;|87.2;|[80.9,90.9];|<79.2;|84.2;|<79.4;|<79.3;|[81.8,91.8];|83.5;|<79.3;|85.4;|<79.7;|[81.7,91.7];|<79.3;|<79;|83.9;|[81.2,91.2];|<79.3;|87.9;|[80.5,90.5];|[82.2,92.2];|>93.1;|84.9;|86.2;|<77.4;|87;|<79;|<79;|84.4;|<79.4;|<76.1;|84.7;|<76.7;|<77.2;|83.4;|<79.2;|<78.1;|<79.1;|<78.4;|>92.9;|86.6;";
        a = new zygotine.M.MeasureList(a).toString('\r\n');
        zygotine.SEG.reset();
        $('#norm').prop('checked', true).trigger('change');
        $('#spExpostats').prop('checked', true).trigger('change');
        $('#withPastData').prop('checked', true).trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#norm').prop('checked', true).trigger('change');
        $('#oel').val('90.0').trigger('change');
        $('#pdMean').val('90').trigger('change');
        $('#pdSd').val('5').trigger('change');
        $('#pdN').val('20').trigger('change');
        $('#prngSeed').val('4444').trigger('change')
    },

    ecrireLogNUniform: function () {
        var a = '<0.0923\n0.0595\n<0.222\n[0.138,0.309]\n<0.149\n2.65\n<0.245\n<0.0892\n0.0817\n0.0251\n<0.0639\n[0.815,1.83]\n<0.0748\n0.628\n0.914\n0.949\n0.0377\n<0.153\n0.0996\n0.0499\n0.0471\n0.605\n0.157\n1.28\n[0.153,0.343]\n<0.0706\n>2.48\n0.169\n4.49\n[0.125,0.282]\n<0.206\n<0.109\n0.0869\n<0.123\n0.229\n<0.121\n0.237\n<0.0376\n<0.0773\n0.153\n<0.205\n[0.282,0.633]\n<0.158\n>2.13\n0.347\n3.57\n>1.84\n0.418\n0.442\n<0.22\n<0.236\n>1.77\n3.71\n<0.13\n>2.06\n<0.126\n0.206\n<0.0819\n3.08\n>1.2\n0.138\n<0.177\n0.209\n<0.158\n0.052\n<0.184\n0.0515\n<0.0566\n>3\n0.328\n>1.26\n1.81\n<0.143\n<0.138\n>2.38\n[0.141,0.317]\n<0.239\n<0.0435\n<0.0586\n<0.254\n0.122\n<0.078\n0.105\n<0.0541\n<0.0175\n[0.163,0.367]\n[0.231,0.52]\n0.573\n0.0546\n[0.116,0.26]\n<0.221\n0.164\n<0.144\n[0.485,1.09]\n<0.063\n0.228\n<0.23\n>1.63\n<0.245\n0.42';
        zygotine.SEG.reset();
        $('#logN').prop('checked', true).trigger('change');
        $('#spUniform').prop('checked', true).trigger('change');
        $('#withPastData').prop('checked', false);
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('.95').trigger('change');
        $('#prngSeed').val('5555').trigger('change')
    },

    ecrireNormUniform: function () {
        var a = ">93;|[80.5,90.5];|84.2;|<78.7;|[80.7,90.7];|>93.7;|86.6;|87;|85.7;|83.1;|86.9;|84.3;|>93.3;|84.3;|86;|88;|<78.9;|<78.2;|84.4;|83.5;|<78.5;|>92.7;|<76.8;|>94.9;|>92.6;|<76.7;|>93.6;|>92.9;|<79.2;|87.4;|84.3;|<78.5;|<77.1;|88.6;|88.7;|<78.8;|86.8;|<76.8;|87;|<79.4;|84;|[79.9,89.9];|[80.9,90.9];|84.3;|<78.9;|87.9;|<76.8;|88.3;|<75.7;|83.4;|85.8;|88.1;|<77.7;|84;|<79.2;|85.9;|<78.4;|<77.5;|<78.5;|<78.4;|87.2;|[80.9,90.9];|<79.2;|84.2;|<79.4;|<79.3;|[81.8,91.8];|83.5;|<79.3;|85.4;|<79.7;|[81.7,91.7];|<79.3;|<79;|83.9;|[81.2,91.2];|<79.3;|87.9;|[80.5,90.5];|[82.2,92.2];|>93.1;|84.9;|86.2;|<77.4;|87;|<79;|<79;|84.4;|<79.4;|<76.1;|84.7;|<76.7;|<77.2;|83.4;|<79.2;|<78.1;|<79.1;|<78.4;|>92.9;|86.6;";
        a = new zygotine.M.MeasureList(a).toString('\r\n');
        zygotine.SEG.reset();
        $('#norm').prop('checked', true).trigger('change');
        $('#spUniform').prop('checked', true).trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('87.0').trigger('change');
        $('#prngSeed').val('6666').trigger('change')
    },
    
    ecrireLogNRiskband: function () {
      
      var a = "2.298|1.337|0.309|>0.231|>0.243|<3.402|<4.41|[0.212,1.477]|[0.1,2.22]"
      a = new zygotine.M.MeasureList(a).toString('\r\n')
      zygotine.SEG.reset()
      $('#logN').prop('checked', true).trigger('change')
      $('#obsValues').val(a).trigger('change');
      $('#oel').val('0.5').trigger('change')
      $('#prngSeed').val('7777').trigger('change')
      $('#meth-riskband').click()
    }
    
};

zygotine.SEG.disableMethodDiv = function(ev) {
  $currChosen = $('.method-chosen')
  $currChosen.removeClass('method-chosen')
  $currChosen.find('input, select').prop('disabled', true)
  $target = ev instanceof Event ? $(ev.target) : $('input[name="method"]:checked')
  let chosenVal = $target.val()
  $chosen = $(`.method-candidate[data-val="${chosenVal}"]`)
  $chosen.addClass('method-chosen')
  $chosen.find('input, select').prop('disabled', false)
}