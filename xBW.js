zygotine.BW = {};
zygotine.X.common = zygotine.BW

// à des fins de débogage
zygotine.BW.memorize = false;
zygotine.BW.history = [];

zygotine.BW.ready = function () {
    zygotine.BW.hideNumericalResults();
    zygotine.BW.setDataEntries();
    zygotine.X.ready()

    $("#resetBtn").click(function () {
        zygotine.BW.reset();
    });

    $("#calcBtn").click(function () {
        console.log(performance.now(), 'calc');
        zygotine.BW.hideNumericalResults();
        zygotine.X.showBePatient();
        window.setTimeout( zygotine.BW.runModel,20);
    });

    //$('#select_Chains').change(function () { zygotine.BW.lastModel.showChain($('#select_Chains').val()); });
    //resultat.workers.ids === ["w1", "w2", "w3", "w4", "w5"]
    //$('#select_workerChains').change(function () { zygotine.BW.lastModel.showWorkerChain($('#select_workerChains').val()); });

    (function () {
        //merci à 'Useless Code'
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
            dnldMuOverallTraceplot = document.getElementById("dnldMuOverallTraceplotBtn"),
            dnldSigmaBetweenTraceplot = document.getElementById("dnldSigmaBetweenTraceplotBtn"),
            dnldSigmaWithinTraceplot = document.getElementById("dnldSigmaWithinTraceplotBtn")

        dnldSampleChains.addEventListener('click', function () {
            var link = document.createElement('a');
            link.setAttribute('download', zygotine.BW.lastModel.result.prngSeed + '_sampleChains.txt');
            link.href = makeTextFile(zygotine.X.concatChains(zygotine.BW.lastModel.result, 'sample'));
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
            link.setAttribute('download', zygotine.BW.lastModel.result.prngSeed + '_burninChains.txt');
            link.href = makeTextFile(zygotine.X.concatChains(zygotine.BW.lastModel.result, 'burnin'));
            document.body.appendChild(link);

            // wait for the link to be added to the document
            window.requestAnimationFrame(function () {
                var event = new MouseEvent('click');
                link.dispatchEvent(event);
                document.body.removeChild(link);
            });
        }, false);
        
        dnldMuOverallTraceplot.addEventListener('click', function() { downloadTraceplot({name: 'muOverall', symbol: 'μ<sub>y</sub>'}) }, false)
        dnldSigmaBetweenTraceplot.addEventListener('click', function() { downloadTraceplot({name: 'sigmaBetween', symbol: 'σ<sub>b</sub>'}) }, false)
        dnldSigmaWithinTraceplot.addEventListener('click', function() { downloadTraceplot({name: 'sigmaWithin', symbol: 'σ<sub>w</sub>'}) }, false)
    })();
};


zygotine.BW.lastModel = null;
zygotine.BW.dataEntries = {};

zygotine.BW.Model = function () {
    zygotine.M.ErrorLogging.call(this);
    this.result = null;
    this.numericalResult = null;
    var entries = zygotine.BW.dataEntries;
    this.modelType = entries["sigmaPrior"].currentValue;
    this.logN = entries.dstrn.currentValue === 'logN';
    this.oel = Number(entries.oel.currentValue);
    this.confidenceLevelForCredibileInterval = Number(entries.confidenceLevelForCredibileInterval.currentValue);
    this.fracThreshold = Number(entries.fracThreshold.currentValue);
    this.percOfInterest = Number(entries.percOfInterest.currentValue);
    this.wwct = Number(zygotine.BW.dataEntries.wwct.currentValue);
    this.rAppapCoverage = Number(entries.rAppapCoverage.currentValue);
    this.probIndOverXThreshold = Number(entries.probIndOverXThreshold.currentValue);
    this.measureRepartition = '';
    this.measureRepartitionByWorker = '';
    this.elapsedTime = {
        runningModel: NaN,
        numericalResults: NaN
    };
};

zygotine.BW.Model.prototype = Object.create(zygotine.M.ErrorLogging.prototype);

zygotine.BW.Model.prototype.className = "Model";

//zygotine.BW.Model.prototype.showChain = function (chainId) {
//    $('#values_Chains').text(this.result.chains[chainId].data.join('\r\n'));
//};

//zygotine.BW.Model.prototype.showWorkerChain = function (chainId) {
//    $('#values_workerChains').text(this.result.chains[chainId].data.join('\r\n'));
//};

zygotine.BW.Model.prototype.getChainAsString = function (chainId) {
    return this.result.chains[chain].data.join('\r\n');
};

zygotine.BW.Model.prototype.validate = function () {
    var entries = zygotine.BW.dataEntries;
    var names = zygotine.BW.getDataEntryNames(entries["sigmaPrior"].currentValue);
    var toutes = names.map(n => entries[n]);
    var fautives = toutes.filter(e => e.validation == null || !e.validation.ok);
    if (fautives.length !== 0) {
        this.addError($.i18n('info-manq'));
    }
};

zygotine.BW.Model.prototype.doCalculation = function () {
    if (this.hasError) {
        //alert('Certaines infos sont incorrectes ou manquantes. Les calculs ne peuvent être effectués.');
        this.result = null;
    } else {
        var mdl;
        var modelParameters;
        var entries = zygotine.BW.dataEntries;
        var ml = new zygotine.M.MeasureList(entries.obsValues.currentValue);
        var logN = entries.dstrn.currentValue === 'logN';
        var nIter = Number(entries.nIter.currentValue);
        var nBurnin = Number(entries.nBurnin.currentValue);
        var monitorBurnin = entries.monitorBurnin.currentValue;
        var mcmc = new zygotine.M.McmcParameters(nIter, nBurnin, monitorBurnin);
        var initMuOverall = Number(entries.initMuOverall.currentValue);
        var initSigmaWithin = Number(entries.initSigmaWithin.currentValue);
        var muOverallLower = Number(entries.muOverallLower.currentValue);
        var muOverallUpper = Number(entries.muOverallUpper.currentValue);

        if (this.modelType === 'expostats') {
            let logSigmaBetweenMu = Number(entries.logSigmaBetweenMu.currentValue);
            let logSigmaBetweenPrec = Number(entries.logSigmaBetweenPrec.currentValue);
            let logSigmaWithinMu = Number(entries.logSigmaWithinMu.currentValue);
            let logSigmaWithinPrec = Number(entries.logSigmaBetweenPrec.currentValue);
            modelParameters = zygotine.BW.createBWModelParametersWithExpostatsSigmaPriors(this.logN, this.oel, initMuOverall, initSigmaWithin, muOverallLower, muOverallUpper, logSigmaBetweenMu, logSigmaBetweenPrec, logSigmaWithinMu, logSigmaWithinPrec);
        } else {
            let sigmaBetweenRangeInf = Number(entries.sigmaBetweenRangeInf.currentValue); // modif février 2020 : on utilisait entries.sigmaWithinRangeInf.currentValue
            let sigmaBetweenRangeSup = Number(entries.sigmaBetweenRangeSup.currentValue); // modif février 2020 : on utilisait entries.sigmaWithinRangeSup.currentValue
            let sigmaWithinRangeInf = Number(entries.sigmaWithinRangeInf.currentValue);
            let sigmaWithinRangeSup = Number(entries.sigmaWithinRangeSup.currentValue);
            modelParameters = zygotine.BW.createBWModelParametersWithSigmaUniformPriors(this.logN, this.oel, initMuOverall, initSigmaWithin, muOverallLower, muOverallUpper, sigmaBetweenRangeInf, sigmaBetweenRangeSup, sigmaWithinRangeInf, sigmaWithinRangeSup);
        }

        mdl = new zygotine.M.BetweenWorkerModel(ml, modelParameters, mcmc);
        let t0 = performance.now();
        this.result = mdl.compute(zygotine.X.prngSeed);
        let t1 = performance.now();
        if (!this.result.hasError) {
            this.elapsedTime.runningModel = Math.round10((t1 - t0) / 1000.0, -4);
        }

        let mtypeUncensored = ml.measureByType.uncensored ? ml.measureByType.uncensored.length : 0;
        let mtypeLessThan = ml.measureByType.uncensored ? ml.measureByType.lessThan.length : 0;
        let mtypeGreaterThan = ml.measureByType.uncensored ? ml.measureByType.greaterThan.length : 0;
        let mtypeInterval = ml.measureByType.uncensored ? ml.measureByType.interval.length : 0;
        this.measureRepartition = zygotine.U.fmt('<span><b>' + $.i18n('non-cens')+ '</b>:&nbsp;{0},&nbsp; <b>' + $.i18n('plus-pt')+ '</b>:&nbsp;{1},&nbsp; <b>' + $.i18n('plus-gr')+ '</b>:&nbsp;{2},&nbsp; <b>' + $.i18n('intervalle')+ '</b>:&nbsp;{3}</span>', mtypeUncensored, mtypeLessThan, mtypeGreaterThan, mtypeInterval);
        this.measureRepartitionByWorker = "<span>" + this.result.workerIds.map(x => '<b>' + x + "</b>:&nbsp;" + ml.measureByWorker[x].length).join(', ') + "<span>";
    }
};

zygotine.BW.Model.prototype.showResults = function () {
    //on crée un tableau vide pour chacun des travailleurs
    var workerIds = this.result.workerIds;
    zygotine.BW._createWorkerNumResTables(workerIds);
    var muWorkerChainIds = workerIds.map(wid => 'mu_' + wid + 'Sample');
    var swChain = this.result.chains.sigmaWithinSample.data;
    var sbChain = this.result.chains.sigmaBetweenSample.data;
    var muOChain = this.result.chains.muOverallSample.data;
    var result;

    let t0 = performance.now();
    result = zygotine.BW.getGlobalResult(this.logN, swChain, sbChain, muOChain,
        this.oel, this.confidenceLevelForCredibileInterval, this.fracThreshold, this.percOfInterest, this.wwct, this.rAppapCoverage, this.probIndOverXThreshold);
    zygotine.BW.insertGlobaResult(result);

    var insertResult = function (resultObject, workerId) {
        const exp10 = 3;
        var r10 = function (n) {
            if (n === 0) {
                return 0;
            }

            var tmp = Math.log10(Math.abs(n));
            var e = (tmp > exp10 ? exp10 : Math.round(tmp)) - exp10;
            var x = Math.round10(n, e);
            return x.toString();
        };

        var resultId, oneResult;

        var insertOneResult = function () {
            var element;
            var targetElement;
            element = $("#" + resultId + "_NumRes_" + workerId).children(".VALUES");
            targetElement = element.children(".B");
            targetElement.text(r10(oneResult.q[1]));
            targetElement = element.children(".A");
            targetElement.text(r10(oneResult.q[0]));
            targetElement = element.children(".C");
            targetElement.text(r10(oneResult.q[2]));
            if (typeof oneResult.risk !== 'undefined') {
                element = $("#" + resultId + "Risk_NumRes_" + workerId).children(".VALUES");
                targetElement = element.children(".R");
                targetElement.text(r10(oneResult.risk));
            }
        };

        var keys = Object.keys(resultObject);
        for (let iKey = 0; iKey < keys.length; iKey++) {
            resultId = keys[iKey];
            oneResult = resultObject[resultId];
            insertOneResult();
        }
    };

    for (let i = 0; i < workerIds.length; i++) {
        let muWrkrChain = zygotine.BW.lastModel.result.chains[muWorkerChainIds[i]].data;
        result = zygotine.X.getNumericalResult(
            this.logN,
            muWrkrChain, swChain,
            this.oel, this.confidenceLevelForCredibileInterval, this.fracThreshold, this.percOfInterest);
        insertResult(result, workerIds[i]);
    }

    zygotine.BW.showNumericalResults(this.logN);
    let t1 = performance.now();
    this.elapsedTime.numericalResults = Math.round10((t1 - t0) / 1000.0, -4);
    $("#numRes_NumRes").children(".VALUES").children('.numRes_NumRes').text(this.elapsedTime.numericalResults);
};

zygotine.BW._createWorkerNumResTables = function (workerIds) {
    var numResWorkerFmt = '<div class="wrkrNumRes" worker="{0}"><div class="wrkrIdent">' + $.i18n('le-trav')+ ' <span style=\"font-weight: bold; font-size:1.1em;\"> \"{0}\"</span></div>\r\n            <div id="gMean_NumRes_{0}" class="logN">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('moy-geom')+ '</span><span class="VALUES"><span class="B"></span> [<span class="A"></span>,&nbsp;<span class="C"></span>]</span>\r\n            </div>\r\n            <div id="aMean_NumRes_{0}" class="logN  norm">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('moy-arith')+ '</span><span class="VALUES"><span class="B"></span> [<span class="A"></span>,&nbsp;<span class="C"></span>] </span>\r\n            </div>\r\n            <div id="aMeanRisk_NumRes_{0}" class="logN  norm">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('risque-moy-arith')+ '</span><span class="VALUES"><span class="R"></span>%</span>\r\n            </div>\r\n            <div id="gSd_NumRes_{0}" class="logN">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('et-geom')+ '</span><span class="VALUES"><span class="B"></span> [<span class="A"></span>,&nbsp;<span class="C"></span>]</span>\r\n            </div>\r\n            <div id="aSd_NumRes_{0}" class="norm">\r\n                <span class="rightJustLbl3">Écart type :&nbsp;</span><span class="VALUES"><span class="B"></span> [<span class="A"></span>,&nbsp;<span class="C"></span>]</span>\r\n            </div>\r\n            <div id="exceedanceFraction_NumRes_{0}" class="logN  norm">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('frac-depasse')+ '</span><span class="VALUES"><span class="B"></span> [<span class="A"></span>,&nbsp;<span class="C"></span>]</span>\r\n            </div>\r\n            <div id="exceedanceFractionRisk_NumRes_{0}" class="logN  norm">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('risque-frac')+ '</span><span class="VALUES"><span class="R"></span>%</span>\r\n            </div>\r\n            <div id="percOfInterest_NumRes_{0}" class="logN norm">\r\n                <span class="rightJustLbl3">Centile critique :&nbsp;</span><span class="VALUES"><span class="B"></span> [<span class="A"></span>,&nbsp;<span class="C"></span>]</span>\r\n            </div>\r\n            <div id="percOfInterestRisk_NumRes_{0}" class="logN norm">\r\n                <span class="rightJustLbl3 append-colon">' + $.i18n('risque-cc')+ '</span><span class="VALUES"><span class="R"></span>%</span>\r\n            </div>\r\n</div>';
    $(".wrkrNumRes").remove(); // si jamais il existait des résultats...
    var ele = $('#workerNumResContainer'); //on cache le 
    ele.hide();
    for (let i = 0; i < workerIds.length; i++) {
        ele.append(zygotine.U.fmt(numResWorkerFmt, workerIds[i]));
    }
};

zygotine.BW.runModel = function () {
    var model = new zygotine.BW.Model();
    model.validate();
    if (model.hasError) {
        zygotine.X.hideBePatient();
        window.setTimeout(zygotine.X.alert, 0, model.messages.getErrors(), 'Prendre note ...');
        return;
    }

    model.doCalculation();
    zygotine.X.hideBePatient();
    if (model.result.hasError) {
        zygotine.X.alert(model.result.messages.getErrors(), 'Prendre note ...');
        return;
    }

    zygotine.BW.lastModel = model;

    var monitorBurnin = model.result.chains.muOverallBurnin.data.length !== 0;
    zygotine.BW.insertWorkerChainOptions(model.result.workerIds, monitorBurnin);
    zygotine.BW.insertChainOptions(monitorBurnin);
    //$('.optBurnin').prop('disabled', model.result.chains.muOverallBurnin.data.length === 0);

    $('#select_Chains').val('muOverallSample').trigger('change');
    //$('#select_Chains');

    $('#select_workerChains').val('mu_' + model.result.workerIds[0] + 'Sample');
    $('#select_workerChains').trigger('change');
    model.showResults();
    //    zygotine.BW.formatGlobalNumericalResults();
    //    zygotine.BW.showNumericalResults(model.logN);

};

zygotine.BW.insertWorkerChainOptions = function (workerIds, monitorBurnin) {
    $('#select_workerChains option').remove();
    var ele = $('#select_workerChains');
    var optFmtBurnin = '<option value="mu_{0}Burnin" class="optBurnin">{0} mu - burnin</option>\r\n';
    var optFmtSample = '<option value="mu_{0}Sample">{0} mu</option>\r\n';
    for (let i = 0; i < workerIds.length; i++) {
        let wid = workerIds[i];
        if (monitorBurnin) {
            ele.append($(zygotine.U.fmt(optFmtBurnin, wid)));
        }

        ele.append(zygotine.U.fmt(optFmtSample, wid));
    }
};

zygotine.BW.insertChainOptions = function (monitorBurnin) {

    var ids = ['muOverall', 'sigmaBetween', 'sigmaWithin'];
    $('#select_Chains option').remove();
    var ele = $('#select_Chains');
    var optFmtBurnin = '<option value="{0}Burnin" class="optBurnin">{0} - burnin</option>\r\n';
    var optFmtSample = '<option value="{0}Sample">{0}</option>\r\n';
    for (let i = 0; i < ids.length; i++) {
        let id = ids[i];
        if (monitorBurnin) {
            ele.append($(zygotine.U.fmt(optFmtBurnin, id)));
        }

        ele.append(zygotine.U.fmt(optFmtSample, id));
    }
};

zygotine.BW.createBWModelParametersWithExpostatsSigmaPriors = function (logN, oel, initMuOverall, initSigmaWithin, muOverallLower, muOverallUpper, logSigmaBetweenMu, logSigmaBetweenPrec, logSigmaWithinMu, logSigmaWithinPrec) {
    let uupOnSds = false;
    var params = new zygotine.M.BetweenWorkerModelParameters(logN, uupOnSds, oel);
    params.initMuOverall = initMuOverall;
    params.initSigmaWithin = initSigmaWithin;
    params.muOverallLower = muOverallLower;
    params.muOverallUpper = muOverallUpper;
    params.logSigmaBetweenMu = logSigmaBetweenMu;
    params.logSigmaBetweenPrec = logSigmaBetweenPrec;
    params.logSigmaWithinMu = logSigmaWithinMu;
    params.logSigmaWithinPrec = logSigmaWithinPrec;
    return params;
};

zygotine.BW.createBWModelParametersWithSigmaUniformPriors = function (logN, oel, initMuOverall, initSigmaWithin, muOverallLower, muOverallUpper, sigmaBetweenRangeInf, sigmaBetweenRangeSup, sigmaWithinRangeInf, sigmaWithinRangeSup) {
    let uupOnSds = true;
    var params = new zygotine.M.BetweenWorkerModelParameters(logN, uupOnSds, oel);
    params.initMuOverall = initMuOverall;
    params.initSigmaWithin = initSigmaWithin;
    params.muOverallLower = muOverallLower;
    params.muOverallUpper = muOverallUpper;
    params.sigmaBetweenRange[0] = sigmaBetweenRangeInf;
    params.sigmaBetweenRange[1] = sigmaBetweenRangeSup;
    params.sigmaWithinRange[0] = sigmaWithinRangeInf;
    params.sigmaWithinRange[1] = sigmaWithinRangeSup;
    return params;
};

zygotine.BW.formatWorkerNumericalResults = function () {
    var result;
    var rep = [];
    var logN = zygotine.BW.dataEntries.dstrn.currentValue === 'logN';
    var oel = Number(zygotine.BW.dataEntries.oel.currentValue);
    var confidenceLevelForCredibileInterval = Number(zygotine.BW.dataEntries.confidenceLevelForCredibileInterval.currentValue);
    var fracThreshold = Number(zygotine.BW.dataEntries.fracThreshold.currentValue)
    var percOfInterest = Number(zygotine.BW.dataEntries.percOfInterest.currentValue);
    var wwct = Number(zygotine.BW.dataEntries.wwct.currentValue);
    var rAppapCoverage = Number(zygotine.BW.dataEntries.rAppapCoverage.currentValue);
    var probIndOverXThreshold = Number(zygotine.BW.dataEntries.probIndOverXThreshold.currentValue);

    //rAppapCoverage, probIndOverXThreshold, wwct
    //les résultats d'ensemble
    //Les résultats par travailleur
    var swChain = zygotine.BW.lastModel.result.chains.sigmaWithinSample.data;
    var sbChain = zygotine.BW.lastModel.result.chains.sigmaBetweenSample.data;
    var muOChain = zygotine.BW.lastModel.result.chains.muOverallSample.data;
    var gRes = zygotine.BW.getGlobalResult(logN, swChain, sbChain, muOChain, oel, confidenceLevelForCredibileInterval, percOfInterest, wwct, rAppapCoverage, probIndOverXThreshold);

};

zygotine.BW.getGlobalResult = function (
    //oui logN en personne, et de facon plus discrète oel aussi appelé oel ou encore tlv.
    logN,
    // les chaines 
    swChain, sbChain, muOChain,
    //les parmètres spécifiques
    oel, confidenceLevelForCredibileInterval, fracThreshold, percOfInterest, wwct, rAppapCoverage, probIndOverXThreshold) {
    var chaine = [];
    var quantile = new zygotine.S.Quantile([(100 - confidenceLevelForCredibileInterval) / 200, .5, 1 - (100 - confidenceLevelForCredibileInterval) / 200]);
    var pNorm = function (x) {
        return zygotine.S.normal.cdf(x, 0, 1, true, false);
    };
    var qNorm = function (x) {
        return zygotine.S.normal.icdf(x, 0, 1, true, false);
    };

    if (logN) {
        let groupGMeanFn = function () {
            chaine = muOChain.map(b => Math.exp(b));
            var rep = quantile.compute(chaine);
            chaine = [];
            return { src: "gMean", logN: logN, q: rep };
        };

        let gSdFn = function (sigmaChain) {
            chaine = sigmaChain.map(b => Math.exp(b));
            var rep = quantile.compute(chaine);
            chaine = [];
            return { src: "gSd", logN: logN, q: rep };
        };

        let rhoFn = function () {
            //within worker correlation
            var sb2, sw2;
            for (let i = 0; i < swChain.length; i++) {
                sb2 = sbChain[i] * sbChain[i];
                sw2 = swChain[i] * swChain[i];
                chaine[i] = sb2 / (sb2 + sw2);
            }
            //proba that correlation > threshold
            var risk = 100 * chaine.filter(function (b) { return b > wwct; }).length / chaine.length;
            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        //rappaport ratio
        let rRatioFn = function (heteroCrtr) {
            // arithmetic mean
            var counts = [0, 0];
            var c = 2 * qNorm(1 - (100 - rAppapCoverage) / 200);
            for (let i = 0; i < sbChain.length; i++) {
                chaine[i] = Math.exp(c * sbChain[i]);

                for (let iHetero = 0; iHetero < heteroCrtr.length; iHetero++) {
                    if (chaine[i] > heteroCrtr[iHetero]) {
                        counts[iHetero]++;
                    }
                }
            }

            let rep = quantile.compute(chaine);
            var risk = [];
            for (let iHetero = 0; iHetero < heteroCrtr.length; iHetero++) {
                risk.push(100 * counts[iHetero] / chaine.length);
            }

            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        let probIndOverXPercFn = function () {
            //probability of individual overexposure based on the critical percentile 
            var qn = qNorm(percOfInterest / 100);
            var logOel = Math.log(oel);
            for (let i = 0; i < muOChain.length; i++) {
                chaine[i] = 100 * (1 - pNorm((logOel - muOChain[i] - qn * swChain[i]) / sbChain[i]));
            }

            var risk = 100 * chaine.filter(function (b) { return b > probIndOverXThreshold; }).length / chaine.length;
            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        let probIndOverXAMeanFn = function () {
            //probability of individual overexposure based on the arithmetic mean
            var logOel = Math.log(oel);
            for (let i = 0; i < muOChain.length; i++) {
                chaine[i] = 100 * (1 - pNorm((logOel - muOChain[i] - 0.5 * swChain[i] * swChain[i]) / sbChain[i]));
            }

            var risk = 100 * chaine.filter(function (b) { return b > probIndOverXThreshold; }).length / chaine.length;
            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        let rep = {
            groupGMean: groupGMeanFn(),
            gSdB: gSdFn(sbChain),
            gSdW: gSdFn(swChain),
            rho: rhoFn(),
            rRatio: rRatioFn([2, 10]),
            probIndOverXPerc: probIndOverXPercFn(),
            probIndOverXAMean: probIndOverXAMeanFn()


        };

        return rep;
    } else {
        let rhoFn = function () {
            //within worker correlation
            var sb2, sw2;
            for (let i = 0; i < swChain.length; i++) {
                sb2 = sbChain[i] * sbChain[i];
                sw2 = swChain[i] * swChain[i];
                chaine[i] = sb2 / (sb2 + sw2);
            }
            //proba that correlation > threshold
            var risk = 100 * chaine.filter(function (b) { return b > wwct; }).length / chaine.length;
            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        let aSdFn = function (sigmaChain) {
            var rep = quantile.compute(sigmaChain);
            return { logN: logN, q: rep };
        };

        let rDiffFn = function () {
            //R difference : difference between the "worst" and "best" worker relative to group mean
            var groupMean = new zygotine.S.Quantile([0.5]).compute(muOChain);
            var c = 2 * qNorm(1 - (100 - rAppapCoverage) / 200) / groupMean;
            for (let i = 0; i < muOChain.length; i++) {
                chaine[i] = c * sbChain[i];
            }

            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep };
        };

        let probIndOverXPercFn = function () {
            //probability of individual overexposure based on the critical percentile 
            var qn = qNorm(percOfInterest / 100);
            for (let i = 0; i < muOChain.length; i++) {
                chaine[i] = 100 * (1 - pNorm((oel - muOChain[i] - qn * swChain[i]) / sbChain[i]));
            }

            var risk = 100 * chaine.filter(function (b) { return b > probIndOverXThreshold; }).length / chaine.length;
            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        let probIndOverXAMeanFn = function () {
            //probability of individual overexposure based on the arithmetic mean
            for (let i = 0; i < muOChain.length; i++) {
                chaine[i] = 100 * (1 - pNorm((oel - muOChain[i]) / sbChain[i]));
            }

            var risk = 100 * chaine.filter(function (b) { return b > probIndOverXThreshold; }).length / chaine.length;
            var rep = quantile.compute(chaine);
            chaine = [];
            return { logN: logN, q: rep, risk: risk };
        };

        let rep = {
            groupMean: { q: quantile.compute(muOChain) },
            aSdB: aSdFn(sbChain),
            aSdW: aSdFn(swChain),
            rho: rhoFn(),
            rDiff: rDiffFn(),
            probIndOverXPerc: probIndOverXPercFn(),
            probIndOverXAMean: probIndOverXAMeanFn()
        };

        return rep;
    }
};

zygotine.BW.insertGlobaResult = function (resultObject) {
    ['.A', '.B', '.C', '.R'].forEach(sel => $(sel).text(''));
    var display1Result = zygotine.X.display1NumericalResult;
    zygotine.BW.hideNumericalResults();
    var t = new Date();
    var element;
    element = $("#" + 'measure' + "_NumRes").children(".VALUES").children('.measure_NumRes').html(zygotine.BW.lastModel.measureRepartition);
    element = $("#" + 'worker' + "_NumRes").children(".VALUES").children('.worker_NumRes').html(zygotine.BW.lastModel.measureRepartitionByWorker);
    element = $("#" + 'info' + "_NumRes").children(".VALUES").children('.when_NumRes').text(t.toLocaleDateString() + " " + t.toLocaleTimeString());
    element = $("#" + 'info' + "_NumRes").children(".VALUES").children('.prngRoot_NumRes').text(" (" + zygotine.BW.lastModel.result.prngSeed.toString() + ")");
    element = $("#" + 'runTime' + "_NumRes").children(".VALUES").children('.runTime_NumRes').text(zygotine.BW.lastModel.elapsedTime.runningModel);
    var id, result;
    var keys = Object.keys(resultObject);
    for (let iKey = 0; iKey < keys.length; iKey++) {
        id = keys[iKey];
        result = resultObject[id];
        display1Result(id, result);
    }
};

zygotine.BW.hideNumericalResults = function () {
    $('#numRes').hide();
    $('#dnldChains').hide();
    $('#workerNumResContainer').hide();
};

zygotine.BW.showNumericalResults = function (logN) {
    zygotine.BW.hideNumericalResults();
    if (logN) {
        $('.wrkrNumRes').children(".norm").hide();
        $('.wrkrNumRes').children(".logN").show();
        $('#numRes').children(".norm").hide();
        $('#numRes').children(".logN").show();
    } else {
        $('.wrkrNumRes').children(".logN").hide();
        $('.wrkrNumRes').children(".norm").show();
        $('#numRes').children(".logN").hide();
        $('#numRes').children(".norm").show();
    }
    $('#workerNumResContainer').show();
    $('#numRes').show();
    $('#dnldBurninChainBtn').hide();
    if (zygotine.BW.dataEntries.monitorBurnin.currentValue && (parseInt(zygotine.BW.dataEntries.nBurnin.currentValue) > 0)) {
        $('#dnldBurninChainBtn').show();
    }

    $('#dnldChains').show();

};

zygotine.BW.defaultEntryValues = (function () {

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
        return { logN: { expostats: val, uupOnSds: val }, norm: { expostats: val, uupOnSds: val } }; // la même valeur, indépendamment de la loi statistique et de la prior sur sd.
    };

    var _4 = function (logN_expostats, logN_uupOnSds, norm_expostats, norm_uupOnSds) {
        return { logN: { expostats: logN_expostats, uupOnSds: logN_uupOnSds }, norm: { expostats: norm_expostats, uupOnSds: norm_uupOnSds } };
    };

    var logN_expostats = new zygotine.M.BetweenWorkerModelParameters(true, false);

    var logN_uupOnSds = new zygotine.M.BetweenWorkerModelParameters(true, true);

    var norm_expostats = new zygotine.M.BetweenWorkerModelParameters(false, false);

    var norm_uupOnSds = new zygotine.M.BetweenWorkerModelParameters(false, true);

    var dev = {};
    // "obsValues"
    dev.nIter = _1('15000'); // modif février 2020 : 15000 plutôt que 150000
    dev.nBurnin = _1('500');
    dev.muOverallLower = _4(logN_expostats.muOverallLower, logN_uupOnSds.muOverallLower, norm_expostats.muOverallLower, norm_uupOnSds.muOverallLower);
    dev.muOverallUpper = _4(logN_expostats.muOverallUpper, logN_uupOnSds.muOverallUpper, norm_expostats.muOverallUpper, norm_uupOnSds.muOverallUpper);
    dev.initMuOverall = _4(logN_expostats.initMuOverall, logN_uupOnSds.initMuOverall, norm_expostats.initMuOverall, norm_uupOnSds.initMuOverall);
    dev.initSigmaWithin = _4(logN_expostats.initSigmaWithin, logN_uupOnSds.initSigmaWithin, norm_expostats.initSigmaWithin, norm_uupOnSds.initSigmaWithin);
    dev.logSigmaBetweenMu = _4(logN_expostats.logSigmaBetweenMu, null, norm_expostats.logSigmaBetweenMu, null);
    dev.logSigmaBetweenPrec = _4(logN_expostats.logSigmaBetweenPrec, null, norm_expostats.logSigmaBetweenPrec, null);
    dev.logSigmaWithinMu = _4(logN_expostats.logSigmaWithinMu, null, norm_expostats.logSigmaWithinMu, null);
    dev.logSigmaWithinPrec = _4(logN_expostats.logSigmaWithinPrec, null, norm_expostats.logSigmaWithinPrec, null);
    dev.sigmaBetweenRangeInf = _4(null, logN_uupOnSds.sigmaBetweenRange[0], null, norm_uupOnSds.sigmaBetweenRange[0]);
    dev.sigmaBetweenRangeSup = _4(null, logN_uupOnSds.sigmaBetweenRange[1], null, norm_uupOnSds.sigmaBetweenRange[1]);
    dev.sigmaWithinRangeInf = _4(null, logN_uupOnSds.sigmaWithinRange[0], null, norm_uupOnSds.sigmaWithinRange[0]);
    dev.sigmaWithinRangeSup = _4(null, logN_uupOnSds.sigmaWithinRange[1], null, norm_uupOnSds.sigmaWithinRange[1]);
    return dev;
})();

zygotine.BW.setDataEntries = function () {

    var entries = zygotine.BW.dataEntries;
    var valueBased = zygotine.X.ValueBasedDataEntry;
    var integer = true;
    var float = false;
    //zygotine.X.SetDefaultEntryValues();
    var defaults = zygotine.BW.defaultEntryValues;
    var key;
    var entry;
    var changeFn;

    /********* obsValues */
    key = 'obsValues';
    entries[key] = new valueBased(key, '', "Requis. Voir la documentation quant à la façon de présenter les observations.");
    entries[key].validate = function () {
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
            let withWorkerInfo = true;
            ml = new zygotine.M.MeasureList(this.currentValue, withWorkerInfo);
            this.hasError = ml.hasError;
            if (this.hasError) {
                this.validation = { ok: false, empty: false, messages: fmtMeasureListMessages(), val: '', ml: '' };
            } else {
                this.validation = { ok: true, empty: false, messages: "(n=" + ml.n + ")", val: this.currentValue, ml: ml.toString() };
            }
        } else {
            let msg = "<div id='obsValuesErr' class='errorMsg'><span>Observations</span>\r\n<ul><li data-i18n='aucune-obs'></li></ul></div>";
            this.validation = { ok: false, empty: true, messages: msg, val: '', ml: '' };
        }
        return this.validation;
    };



    changeFn = (function (entry) {
        return function () {
            console.log($('#calcBtn').queue(), 'queue');
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
    })(entries[key]);

    entries.obsValues.element.change(changeFn);

    //mcmc

    /********* nIter */
    key = 'nIter';
    entries[key] = new valueBased(key, defaults.nIter.logN.expostats, "Requis. Un nombre entier compris entre 500 et  500000", integer, 500, 500000);

    /********* nIter */
    key = 'nBurnin';
    entries[key] = new valueBased(key, defaults.nBurnin.logN.expostats, "Requis. Un nombre entier compris entre 0 et  50000", integer, 0, 15000);

    //paramètres pour l'interprétation des données

    /********* oel */
    key = 'oel';
    entries[key] = new valueBased(key, '', "Requis. Un nombre réel positif correspondant à la valeur limite d'exposition.", float);

    /********* confidenceLevelForCredibileInterval */
    key = 'confidenceLevelForCredibileInterval';
    entries[key] = new valueBased(key, 90, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);

    /********* fracThreshold */
    key = 'fracThreshold';
    entries[key] = new valueBased(key, 5, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);

    /********* percOfInterest */
    key = 'percOfInterest';
    entries[key] = new valueBased(key, 95, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);

    /********* wwct */
    key = 'wwct';
    entries[key] = new valueBased(key, 0.2, "Requis. Un nombre réel compris entre 0 et 1.", float, 0, 1);

    /********* rAppapCoverage */
    key = 'rAppapCoverage';
    entries[key] = new valueBased(key, 80, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);

    /********* probIndOverXThreshold */
    key = 'probIndOverXThreshold';
    entries[key] = new valueBased(key, 20, "Requis. Un entier compris entre 1 et 99.", integer, 1, 99);

    //prior de mu

    /********* muOverallLower */
    key = 'muOverallLower';
    entries[key] = new valueBased(key, '', 'Requis. Un nombre réel établissant un minimum pour la prior de mu.', float);

    /********* muOverallUpper */
    key = 'muOverallUpper';
    entries[key] = new valueBased(key, '', 'Requis. Un nombre réel établissant un maximum pour la prior de mu.', float);

    //valeurs initiales

    /********* initMuOverall */
    key = 'initMuOverall';
    entries[key] = new valueBased(key, '', 'Requis. Un nombre réel établissant un minimum pour la prior de mu.', float);

    /********* initSigmaWithin */
    key = 'initSigmaWithin';
    entries[key] = new valueBased(key, '', 'Requis. Un nombre réel établissant un maximum pour la prior de mu.', float);

    // sigmaPrior : uniform vs expostats

    /********* logSigmaBetweenMu */
    key = 'logSigmaBetweenMu';
    entries[key] = new valueBased(key, '', "Nombre réel requis lorsque la prior pour sigma est 'expostats'. Donne la valeur de logSigmaMu.", float);


    /********* logSigmaBetweenPrec */
    key = 'logSigmaBetweenPrec';
    entries[key] = new valueBased(key, '', "Nombre réel requis lorsque la prior pour sigma est 'expostats'. Donne la valeur de logSigmaPrec.", float);

    /********* logSigmaWithinMu */
    key = 'logSigmaWithinMu';
    entries[key] = new valueBased(key, '', "Nombre réel requis lorsque la prior pour sigma est 'expostats'. Donne la valeur de logSigmaMu.", float);

    /********* logSigmaWithinPrec */
    key = 'logSigmaWithinPrec';
    entries[key] = new valueBased(key, '', "Nombre réel requis lorsque la prior pour sigma est 'expostats'. Donne la valeur de logSigmaPrec.", float);

    /********* sigmaBetweenRangeInf */
    key = 'sigmaBetweenRangeInf';
    entries[key] = new valueBased(key, '', float);

    /********* sigmaBetweenRangeSup */
    key = 'sigmaBetweenRangeSup';
    entries[key] = new valueBased(key, '', float);

    /********* sigmaWithinRangeInf */
    key = 'sigmaWithinRangeInf';
    entries[key] = new valueBased(key, '', float);

    /********* sigmaWithinRangeSup */
    key = 'sigmaWithinRangeSup';
    entries[key] = new valueBased(key, '', float);

    //radio

    var radioBased = zygotine.X.RadioButtonBasedDataEntry;

    /********* sigmaPrior */
    key = 'sigmaPrior';
    entries[key] = new radioBased(key, 'spExpostats');
    changeFn = (function (entry) {
        var fn = function () {
            var oStr = 'input[name=' + entry.name + ']:checked';
            entry.currentValue = $(oStr).val();
            var isExpostats = entry.currentValue === 'expostats';
            if (isExpostats) {
                $('.uniformDep').prop('disabled', true);
                $('.expostatsDep').prop('disabled', false);
            } else {
                $('.expostatsDep').prop('disabled', true);
                $('.uniformDep').prop('disabled', false);
            }
        };

        return fn;
    })(zygotine.BW.dataEntries[key]);
    entries[key].element.change(changeFn);

    /********* monitorBurnin */
    key = 'monitorBurnin';
    entries[key] = new zygotine.X.checkboxBasedDataEntry(key, false);
    entries[key].element.change(entries[key].getChangeFn());

    /********* dstrn */
    key = 'dstrn';
    entries[key] = new radioBased(key, 'logN');
    entries[key].reset();
    entries[key].element.change(function () {
        var oStr = 'input[name=dstrn]:checked';
        entries[key].currentValue = $(oStr).val();
        zygotine.BW.setDefaultsForDistribution(entries[key].currentValue);
    });

    zygotine.X.setDataEntries()
    zygotine.BW.reset();
};

zygotine.BW.getDataEntryNames = (function () {

    var inCommon = [
        'obsValues',
        'nIter', 'nBurnin', 'monitorBurnin',
        'oel', 'confidenceLevelForCredibileInterval', /*'exceedanceFrac',*/ 'percOfInterest', 'rAppapCoverage', 'probIndOverXThreshold', 'wwct',
        'initMuOverall', 'initSigmaWithin', 'muOverallLower', 'muOverallUpper'];
    var specifics =
        {
            expostats: ['logSigmaBetweenMu', 'logSigmaBetweenPrec', 'logSigmaWithinMu', 'logSigmaWithinPrec'],
            uupOnSds: ['sigmaBetweenRangeInf', 'sigmaBetweenRangeSup', 'sigmaWithinRangeInf', 'sigmaWithinRangeSup']
        };

    var fn = function (sigmaPrior) {
        var rep = [];

        Array.prototype.push.apply(rep, inCommon);
        Array.prototype.push.apply(rep, specifics[sigmaPrior]);
        return rep;
    };

    return fn;
})();

zygotine.BW.reset = function () {
    zygotine.BW.hideNumericalResults();
    zygotine.X.hideBePatient();
    var entries = zygotine.BW.dataEntries;
    entries.dstrn.reset();
    entries.obsValues.reset();
    entries.nIter.reset();
    entries.nBurnin.reset();
    entries.monitorBurnin.reset();
    // ++ monitorBurnin ++
    //paramètres pour l'interprétation des données
    entries.oel.reset();
    entries.confidenceLevelForCredibileInterval.reset();
    entries.fracThreshold.reset();
    entries.percOfInterest.reset();
    entries.wwct.reset();
    //rAppapCoverage, probIndOverXThreshold, wwct
    entries.rAppapCoverage.reset();
    entries.probIndOverXThreshold.reset();

    //radio ... autres que dstrn
    entries.sigmaPrior.reset();
};

zygotine.BW.setDefaultsForDistribution = function (loi) {
    var entries = zygotine.BW.dataEntries;
    var defaults = zygotine.BW.defaultEntryValues;

    //valeurs initiales
    entries.initMuOverall.element.val(defaults.initMuOverall[loi].expostats).trigger('change');
    entries.initSigmaWithin.element.val(defaults.initSigmaWithin[loi].expostats).trigger('change');
    //prior sur mu
    entries.muOverallLower.element.val(defaults.muOverallLower[loi].expostats).trigger('change');
    entries.muOverallUpper.element.val(defaults.muOverallUpper[loi].expostats).trigger('change');
    // sigmaPrior name / uniform vs expostats
    entries.logSigmaBetweenMu.element.val(defaults.logSigmaBetweenMu[loi].expostats).trigger('change');
    entries.logSigmaBetweenPrec.element.val(defaults.logSigmaBetweenPrec[loi].expostats).trigger('change');
    entries.logSigmaWithinMu.element.val(defaults.logSigmaWithinMu[loi].expostats).trigger('change');
    entries.logSigmaWithinPrec.element.val(defaults.logSigmaWithinPrec[loi].expostats).trigger('change');
    entries.sigmaBetweenRangeInf.element.val(defaults.sigmaBetweenRangeInf[loi].uupOnSds).trigger('change');
    entries.sigmaBetweenRangeSup.element.val(defaults.sigmaBetweenRangeSup[loi].uupOnSds).trigger('change');
    entries.sigmaWithinRangeInf.element.val(defaults.sigmaWithinRangeInf[loi].uupOnSds).trigger('change');
    entries.sigmaWithinRangeSup.element.val(defaults.sigmaWithinRangeSup[loi].uupOnSds).trigger('change');
    
    zygotine.X.setDefaultsForDistribution(loi)
};

zygotine.BW.tests = {
    ecrireLogNExpostats: function () {
        var a = '0.303;w1|<0.436;w1|0.31;w1|<0.401;w1|<0.424;w1|<0.422;w1|0.299;w1|<0.43;w1|<0.397;w1|<0.439;w1|0.301;w1|0.331;w1|<0.423;w1|0.303;w1|0.321;w1|<0.365;w1|<0.39;w1|<0.371;w1|<0.431;w1|<0.4;w1|<0.442;w10|<0.415;w10|<0.444;w10|<0.409;w10|0.299;w10|>0.229;w10|<0.425;w10|0.315;w10|0.33;w10|[0.214,0.481];w10|0.341;w10|0.305;w10|>0.23;w10|0.316;w10|0.299;w10|<0.42;w10|0.312;w10|<0.424;w10|<0.434;w10|0.298;w10|0.323;w11|0.303;w11|0.335;w11|[0.208,0.467];w11|0.312;w11|0.317;w11|>0.241;w11|<0.444;w11|0.312;w11|<0.436;w11|>0.23;w10|0.304;w11|0.307;w11|0.299;w11|>0.228;w11|>0.234;w11|0.326;w11|[0.214,0.48];w11|>0.228;w11|0.299;w11|<0.427;w12|>0.241;w11|<0.428;w12|0.325;w12|[0.208,0.467];w12|0.326;w12|[0.199,0.447];w12|0.323;w12|0.315;w12|<0.397;w12|<0.431;w12|<0.395;w12|<0.432;w12|<0.438;w12|<0.429;w12|0.321;w12|0.298;w12|<0.416;w12|<0.443;w12|[0.204,0.459];w12|>0.231;w13|[0.212,0.477];w13|<0.402;w13|>0.243;w13|0.298;w13|0.337;w13|[0.206,0.463];w13|0.326;w12|>0.229;w13|0.305;w13|0.318;w13|0.319;w13|<0.41;w13|<0.41;w13|0.315;w13|0.305;w13|[0.207,0.465];w13|<0.435;w13|0.309;w13|[0.211,0.476];w13|0.314;w14|[0.209,0.471];w14|<0.393;w14|<0.436;w14|<0.406;w14|<0.44;w14|0.301;w14|<0.445;w14|<0.417;w14|0.334;w14|[0.203,0.458];w14|<0.406;w14|0.322;w14|0.316;w14|0.315;w14|0.298;w14|>0.243;w13|<0.435;w14|[0.209,0.47];w14|<0.421;w14|<0.414;w15|<0.443;w15|<0.389;w15|<0.43;w15|0.334;w14|<0.403;w15|0.338;w15|>0.229;w15|0.305;w15|<0.437;w15|<0.441;w15|<0.434;w15|>0.228;w15|0.333;w15|<0.422;w15|0.312;w15|0.302;w15|<0.441;w15|<0.396;w15|0.311;w15|[0.213,0.48];w16|0.331;w16|0.301;w16|<0.363;w16|[0.22,0.495];w16|<0.42;w16|0.328;w16|<0.401;w16|<0.403;w16|0.314;w16|<0.437;w16|0.316;w16|<0.44;w16|<0.421;w16|>0.229;w15|0.31;w16|<0.436;w16|[0.214,0.482];w16|>0.236;w16|<0.408;w16|<0.394;w17|<0.439;w17|0.315;w17|<0.436;w17|<0.445;w17|<0.393;w17|<0.426;w17|0.303;w17|<0.421;w17|<0.401;w17|0.3;w17|[0.202,0.454];w17|0.319;w17|0.317;w17|0.31;w17|[0.222,0.499];w17|<0.43;w17|<0.427;w17|<0.389;w17|<0.42;w17|0.299;w18|[0.213,0.479];w18|<0.413;w18|0.311;w18|<0.421;w18|0.297;w18|<0.416;w18|[0.222,0.499];w17|0.311;w18|<0.407;w18|0.301;w18|0.311;w18|0.303;w18|<0.436;w18|<0.433;w18|<0.406;w18|<0.438;w18|<0.416;w18|0.312;w18|0.315;w18|0.301;w19|>0.251;w19|<0.422;w19|>0.228;w19|0.341;w19|0.336;w19|0.334;w19|[0.214,0.483];w19|[0.209,0.47];w19|[0.213,0.479];w18|[0.217,0.489];w19|>0.243;w19|[0.217,0.489];w19|0.317;w19|0.321;w19|0.325;w19|0.315;w19|0.305;w19|0.314;w19|>0.247;w19|0.298;w2|0.316;w2|<0.437;w2|<0.426;w2|0.317;w2|>0.228;w2|<0.443;w2|>0.252;w2|[0.198,0.447];w2|<0.418;w2|<0.441;w2|<0.43;w2|0.325;w2|<0.434;w2|>0.251;w19|[0.212,0.476];w2|<0.413;w2|<0.42;w2|<0.41;w2|<0.436;w2|[0.207,0.465];w20|[0.2,0.45];w20|<0.408;w20|<0.39;w20|<0.429;w20|0.309;w20|<0.426;w20|>0.252;w2|<0.384;w20|0.3;w20|0.306;w20|<0.439;w20|<0.408;w20|[0.198,0.446];w20|0.318;w20|<0.41;w20|[0.223,0.503];w20|<0.434;w20|<0.405;w20|0.301;w20|>0.228;w3|>0.235;w3|>0.247;w3|>0.25;w3|[0.223,0.503];w20|0.321;w3|0.331;w3|0.321;w3|0.329;w3|<0.407;w3|0.331;w3|[0.213,0.479];w3|0.298;w3|>0.228;w3|0.32;w3|0.313;w3|0.334;w3|[0.216,0.486];w3|0.337;w3|0.327;w3|0.317;w4|0.299;w4|>0.234;w4|<0.42;w4|0.304;w4|<0.384;w4|0.3;w4|<0.379;w4|<0.428;w4|<0.42;w4|>0.25;w3|<0.433;w4|<0.396;w4|<0.404;w4|0.306;w4|<0.414;w4|0.302;w4|0.298;w4|<0.393;w4|<0.384;w4|<0.403;w5|[0.21,0.473];w5|0.331;w5|0.3;w5|>0.236;w5|>0.234;w4|>0.251;w5|<0.436;w5|<0.437;w5|0.298;w5|0.304;w5|0.327;w5|0.304;w5|0.312;w5|>0.229;w5|0.309;w5|0.329;w5|>0.24;w5|<0.425;w5|<0.442;w5|0.315;w6|<0.435;w6|<0.427;w6|<0.42;w6|0.321;w6|<0.425;w6|>0.23;w6|<0.437;w6|<0.414;w6|0.301;w6|<0.443;w6|0.307;w6|0.303;w6|0.309;w6|0.318;w6|0.31;w6|>0.251;w5|0.305;w6|<0.409;w6|0.316;w6|0.302;w7|[0.227,0.51];w7|0.332;w7|>0.23;w6|0.322;w7|>0.241;w7|[0.215,0.483];w7|[0.218,0.49];w7|0.335;w7|0.316;w7|>0.231;w7|0.312;w7|<0.443;w7|>0.238;w7|0.308;w7|>0.232;w7|0.318;w7|>0.242;w7|<0.378;w7|0.302;w7|0.318;w8|0.329;w8|0.323;w8|0.319;w8|>0.231;w8|0.323;w8|<0.437;w8|0.299;w8|0.332;w8|0.315;w8|0.33;w8|>0.238;w8|[0.209,0.47];w8|[0.22,0.495];w8|0.314;w8|0.32;w8|<0.419;w8|>0.25;w8|>0.242;w7|0.333;w8|<0.433;w9|[0.206,0.463];w9|<0.442;w9|0.302;w9|0.303;w9|[0.209,0.47];w9|[0.209,0.469];w9|<0.423;w9|<0.402;w9|<0.395;w9|>0.25;w8|<0.348;w9|0.299;w9|<0.411;w9|<0.445;w9|<0.44;w9|<0.431;w9|<0.423;w9|<0.41;w9|<0.405;w9';
        a = new zygotine.M.MeasureList(a).toString("\r\n", "\t");
        zygotine.BW.reset();
        $('#logN').prop('checked', true).trigger('change');
        $('#spExpostats').prop('checked', true).trigger('change');
        $('#oel').val('1.0').trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('1.0').trigger('change');
        $('#nIter').val('5000').trigger('change');
        $('#nBurnin').val('0').trigger('change');
        $('#prngSeed').val('7777').trigger('change')
    },

    ecrireNormExpostats: function () {
        var a = "<91.5;w01|83.7;w01|<89.4;w01|<91;w01|88.6;w01|85.2;w10|>84.3;w10|88.6;w01|<89.8;w10|87.5;w10|<90.9;w02|<92.3;w02|<86.7;w02|<91.9;w02|<83;w02|>82.6;w03|<92.3;w02|86.8;w03|>80.1;w03|87.6;w03|84.2;w04|88.8;w04|87.8;w04|84.5;w04|>82.6;w03|<90.3;w05|<83.9;w05|<91.4;w05|89.6;w04|<88.8;w05|84.1;w06|<92.7;w05|90.3;w06|89.1;w06|<95.8;w06|<95.4;w07|<94.9;w07|85;w07|[75.9,100];w07|<81.9;w07|[72.8,96.3];w08|[75.9,100];w07|86.5;w08|>80.2;w08|85.3;w08|>85.8;w08|<91.8;w09|83.7;w09|<89.9;w09|[73.5,97.2];w09";
        //new zygotine.M.MeasureList(a).transform(function (x) { return (205.0 / 156)})
        a = new zygotine.M.MeasureList(a).toString("\r\n", "\t");
        zygotine.BW.reset();
        $('#norm').prop('checked', true).trigger('change');
        $('#spExpostats').prop('checked', true).trigger('change');
        $('#oel').val('90.0').trigger('change');
        $('#nIter').val('5000').trigger('change');
        $('#nBurnin').val('0').trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#prngSeed').val('8888').trigger('change')
    },

    ecrireLogNUniform: function () {
        var a = '0.303;w1|<0.436;w1|0.31;w1|<0.401;w1|<0.424;w1|<0.422;w1|0.299;w1|<0.43;w1|<0.397;w1|<0.439;w1|0.301;w1|0.331;w1|<0.423;w1|0.303;w1|0.321;w1|<0.365;w1|<0.39;w1|<0.371;w1|<0.431;w1|<0.4;w1|<0.442;w10|<0.415;w10|<0.444;w10|<0.409;w10|0.299;w10|>0.229;w10|<0.425;w10|0.315;w10|0.33;w10|[0.214,0.481];w10|0.341;w10|0.305;w10|>0.23;w10|0.316;w10|0.299;w10|<0.42;w10|0.312;w10|<0.424;w10|<0.434;w10|0.298;w10|0.323;w11|0.303;w11|0.335;w11|[0.208,0.467];w11|0.312;w11|0.317;w11|>0.241;w11|<0.444;w11|0.312;w11|<0.436;w11|>0.23;w10|0.304;w11|0.307;w11|0.299;w11|>0.228;w11|>0.234;w11|0.326;w11|[0.214,0.48];w11|>0.228;w11|0.299;w11|<0.427;w12|>0.241;w11|<0.428;w12|0.325;w12|[0.208,0.467];w12|0.326;w12|[0.199,0.447];w12|0.323;w12|0.315;w12|<0.397;w12|<0.431;w12|<0.395;w12|<0.432;w12|<0.438;w12|<0.429;w12|0.321;w12|0.298;w12|<0.416;w12|<0.443;w12|[0.204,0.459];w12|>0.231;w13|[0.212,0.477];w13|<0.402;w13|>0.243;w13|0.298;w13|0.337;w13|[0.206,0.463];w13|0.326;w12|>0.229;w13|0.305;w13|0.318;w13|0.319;w13|<0.41;w13|<0.41;w13|0.315;w13|0.305;w13|[0.207,0.465];w13|<0.435;w13|0.309;w13|[0.211,0.476];w13|0.314;w14|[0.209,0.471];w14|<0.393;w14|<0.436;w14|<0.406;w14|<0.44;w14|0.301;w14|<0.445;w14|<0.417;w14|0.334;w14|[0.203,0.458];w14|<0.406;w14|0.322;w14|0.316;w14|0.315;w14|0.298;w14|>0.243;w13|<0.435;w14|[0.209,0.47];w14|<0.421;w14|<0.414;w15|<0.443;w15|<0.389;w15|<0.43;w15|0.334;w14|<0.403;w15|0.338;w15|>0.229;w15|0.305;w15|<0.437;w15|<0.441;w15|<0.434;w15|>0.228;w15|0.333;w15|<0.422;w15|0.312;w15|0.302;w15|<0.441;w15|<0.396;w15|0.311;w15|[0.213,0.48];w16|0.331;w16|0.301;w16|<0.363;w16|[0.22,0.495];w16|<0.42;w16|0.328;w16|<0.401;w16|<0.403;w16|0.314;w16|<0.437;w16|0.316;w16|<0.44;w16|<0.421;w16|>0.229;w15|0.31;w16|<0.436;w16|[0.214,0.482];w16|>0.236;w16|<0.408;w16|<0.394;w17|<0.439;w17|0.315;w17|<0.436;w17|<0.445;w17|<0.393;w17|<0.426;w17|0.303;w17|<0.421;w17|<0.401;w17|0.3;w17|[0.202,0.454];w17|0.319;w17|0.317;w17|0.31;w17|[0.222,0.499];w17|<0.43;w17|<0.427;w17|<0.389;w17|<0.42;w17|0.299;w18|[0.213,0.479];w18|<0.413;w18|0.311;w18|<0.421;w18|0.297;w18|<0.416;w18|[0.222,0.499];w17|0.311;w18|<0.407;w18|0.301;w18|0.311;w18|0.303;w18|<0.436;w18|<0.433;w18|<0.406;w18|<0.438;w18|<0.416;w18|0.312;w18|0.315;w18|0.301;w19|>0.251;w19|<0.422;w19|>0.228;w19|0.341;w19|0.336;w19|0.334;w19|[0.214,0.483];w19|[0.209,0.47];w19|[0.213,0.479];w18|[0.217,0.489];w19|>0.243;w19|[0.217,0.489];w19|0.317;w19|0.321;w19|0.325;w19|0.315;w19|0.305;w19|0.314;w19|>0.247;w19|0.298;w2|0.316;w2|<0.437;w2|<0.426;w2|0.317;w2|>0.228;w2|<0.443;w2|>0.252;w2|[0.198,0.447];w2|<0.418;w2|<0.441;w2|<0.43;w2|0.325;w2|<0.434;w2|>0.251;w19|[0.212,0.476];w2|<0.413;w2|<0.42;w2|<0.41;w2|<0.436;w2|[0.207,0.465];w20|[0.2,0.45];w20|<0.408;w20|<0.39;w20|<0.429;w20|0.309;w20|<0.426;w20|>0.252;w2|<0.384;w20|0.3;w20|0.306;w20|<0.439;w20|<0.408;w20|[0.198,0.446];w20|0.318;w20|<0.41;w20|[0.223,0.503];w20|<0.434;w20|<0.405;w20|0.301;w20|>0.228;w3|>0.235;w3|>0.247;w3|>0.25;w3|[0.223,0.503];w20|0.321;w3|0.331;w3|0.321;w3|0.329;w3|<0.407;w3|0.331;w3|[0.213,0.479];w3|0.298;w3|>0.228;w3|0.32;w3|0.313;w3|0.334;w3|[0.216,0.486];w3|0.337;w3|0.327;w3|0.317;w4|0.299;w4|>0.234;w4|<0.42;w4|0.304;w4|<0.384;w4|0.3;w4|<0.379;w4|<0.428;w4|<0.42;w4|>0.25;w3|<0.433;w4|<0.396;w4|<0.404;w4|0.306;w4|<0.414;w4|0.302;w4|0.298;w4|<0.393;w4|<0.384;w4|<0.403;w5|[0.21,0.473];w5|0.331;w5|0.3;w5|>0.236;w5|>0.234;w4|>0.251;w5|<0.436;w5|<0.437;w5|0.298;w5|0.304;w5|0.327;w5|0.304;w5|0.312;w5|>0.229;w5|0.309;w5|0.329;w5|>0.24;w5|<0.425;w5|<0.442;w5|0.315;w6|<0.435;w6|<0.427;w6|<0.42;w6|0.321;w6|<0.425;w6|>0.23;w6|<0.437;w6|<0.414;w6|0.301;w6|<0.443;w6|0.307;w6|0.303;w6|0.309;w6|0.318;w6|0.31;w6|>0.251;w5|0.305;w6|<0.409;w6|0.316;w6|0.302;w7|[0.227,0.51];w7|0.332;w7|>0.23;w6|0.322;w7|>0.241;w7|[0.215,0.483];w7|[0.218,0.49];w7|0.335;w7|0.316;w7|>0.231;w7|0.312;w7|<0.443;w7|>0.238;w7|0.308;w7|>0.232;w7|0.318;w7|>0.242;w7|<0.378;w7|0.302;w7|0.318;w8|0.329;w8|0.323;w8|0.319;w8|>0.231;w8|0.323;w8|<0.437;w8|0.299;w8|0.332;w8|0.315;w8|0.33;w8|>0.238;w8|[0.209,0.47];w8|[0.22,0.495];w8|0.314;w8|0.32;w8|<0.419;w8|>0.25;w8|>0.242;w7|0.333;w8|<0.433;w9|[0.206,0.463];w9|<0.442;w9|0.302;w9|0.303;w9|[0.209,0.47];w9|[0.209,0.469];w9|<0.423;w9|<0.402;w9|<0.395;w9|>0.25;w8|<0.348;w9|0.299;w9|<0.411;w9|<0.445;w9|<0.44;w9|<0.431;w9|<0.423;w9|<0.41;w9|<0.405;w9';
        a = new zygotine.M.MeasureList(a).toString("\r\n", "\t");
        var ml = new zygotine.M.MeasureList(a);
        ml = ml.transform(function (x) { return Math.round10(100 * (1.474358974 * x - 0.241923077), -4); });
        zygotine.BW.reset();
        $('#logN').prop('checked', true).trigger('change');
        $('#spUniform').prop('checked', true).trigger('change');
        $('#oel').val('.5').trigger('change');
        $('#nIter').val('15000').trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#prngSeed').val('9999').trigger('change')
    },

    ecrireNormUniform: function () {
        var a = "<91.5; w1 | 83.7; w1 |<89.4; w1 |<91; w1 | 88.6; w1 | 85.2; w10 |>84.3; w10 | 88.6; w1 |<89.8; w10 | 87.5; w10 |<90.9; w2 |<92.3; w2 |<86.7; w2 |<91.9; w2 |<83; w2 |>82.6; w3 |<92.3; w2 | 86.8; w3 |>80.1; w3 | 87.6; w3 | 84.2; w4 | 88.8; w4 | 87.8; w4 | 84.5; w4 |>82.6; w3 |<90.3; w5 |<83.9; w5 |<91.4; w5 | 89.6; w4 |<88.8; w5 | 84.1; w6 |<92.7; w5 | 90.3; w6 | 89.1; w6 |<95.8; w6 |<95.4; w7 |<94.9; w7 | 85; w7 | [75.9, 100]; w7 |<81.9; w7 | [72.8, 96.3]; w8 | [75.9, 100]; w7 | 86.5; w8 |>80.2; w8 | 85.3; w8 |>85.8; w8 |<91.8; w9 | 83.7; w9 |<89.9; w9 | [73.5, 97.2]; w9";
        a = new zygotine.M.MeasureList(a).toString("\r\n", "\t");
        zygotine.BW.reset();
        $('#norm').prop('checked', true).trigger('change');
        $('#spUniform').prop('checked', true).trigger('change');
        $('#obsValues').val(a).trigger('change');
        $('#oel').val('90.0').trigger('change');
        $('#nIter').val('15000').trigger('change');
        $('#prngSeed').val('1010').trigger('change')
    },

    // 
    testLogN2020_02_1: function () {
        var a = "0.296;w1|<0.39;w1|0.278;w1|0.296;w1|0.336;w2|>0.235;w2|<0.398;w3|<0.409;w3|[0.222, 0.499];w3|<0.401;w4|0.321;w4|[0.222, 0.499];w3|[0.196, 0.442];w5|<0.403;w5|0.298;w5";
        a = new zygotine.M.MeasureList(a).toString("\r\n", "\t");
        zygotine.BW.reset();
        zygotine.BW.reset();
        $('#logN').prop('checked', true).trigger('change');
        $('#spUniform').prop('checked', true).trigger('change');
        $('#oel').val('1.0').trigger('change');
        $('#nIter').val('20000').trigger('change');
        $('#nBurnin').val('0').trigger('change');
        $('#obsValues').val(a).trigger('change');   
        zygotine.X.prngSeed = 12;
    }
};