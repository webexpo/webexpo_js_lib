/* eslint 
    valid-jsdoc: 0
    no-extra-parens: 0
*/
/// <reference path="A.js" />
/// <reference path="M0.js" />

/***** McmcParameters    *************/

/**
* Represents MCMC parameters to be used when running a model
* Tous les paramètres sont requis et doivent être du type attendu
* @constructor
* @param {number} nIter - sample size!
* @param {number} nBurnin - burnin size
* xparam {number} nThin - thinning parameter (1 no thinning, 2 means that 1 of 2 is retained ...)
* @param {boolean} monitorBurnin
*
* @property {boolean} hasError
* @property {zygotine.M.Messages} messages
* @property {boolean} monitorBurnin
* @property {number} nBurnin
* @property {number} nIter
* @property {number} nThin
*/
zygotine.M.McmcParameters = function (nIter, nBurnin, monitorBurnin) {
    zygotine.M.ErrorLogging.call(this);
    this.monitorBurnin = false;
    // this.nThin = 1;
    this.nBurnin = 500;
    this.nIter = 15000;

    if (typeof (monitorBurnin) === "boolean") {
        this.monitorBurnin = monitorBurnin;
    } else {
        this.monitorBurnin = false;
    }


    if (!this.messages.hasError && (typeof nBurnin !== "undefined")) {
        this.nBurnin = parseInt(nBurnin);
        if (isNaN(this.nBurnin) || this.nBurnin < 0) {
            this.addError("Parameter nBurnin is not well defined!");
        }
    }

    if (!this.messages.hasError && (typeof nIter !== "undefined")) {
        this.nIter = parseInt(nIter);
        if (isNaN(this.nIter) || this.nIter < 0) {
            this.addError("Parameter nIter is not well defined!");
        }
    }
};

zygotine.M.McmcParameters.prototype = Object.create(zygotine.M.ErrorLogging.prototype);
zygotine.M.McmcParameters.prototype.className = "McmcParameters";
zygotine.M.McmcParameters.prototype.toString = function () {
    return zygotine.U.fmt("nBurnin={0}, nIter={1}, nThin={2}, monitorBurnin={3}", this.nBurnin, this.nIter, /* this.nThin,*/this.monitorBurnin);
};

zygotine.M.McmcParameters.getDefaults = function () {
    return new zygotine.M.McmcParameters(15000, 500, false);
};

/*******************************************/

/***** ModelParameters *******/

/**
* Represents parameters for a Model
* @constructor
* @param {boolean} logN - true for lognormal distribution, false for normal distribution
*
* @property {boolean} logN - true for lognormal distribution, false for normal distribution
*/
zygotine.M.ModelParameters = function (logN, oel) {
    this.logN = logN;
    // modif février 2020 : l'interface usager exige un oel, mais pour des tests ce n'est pas le cas, 
    // donc on s'assure d'avoir une valeur.
    if (typeof oel === 'undefined') {
        this.oel = this.logN ? 1 : 90; // modif février 2020
    } else {
        this.oel = oel;
    }
};

zygotine.M.ModelParameters.prototype = {
    toString: function () {
        return zU.toString(this);
    }
};

/*******************************************/



/****** Chain    **************/

/**
* @constructor
* @param {string} label étiquette permettant d'identifier la chaine
* @param {number[]} data un tableau donnant les valeurs de la chaine
* @property {string} label
* @property {number[]} data
*/
zygotine.M.Chain = function (label, data) {
    this.label = label;
    this.data = data;
    this._quantile = [];
    this.quantileProbs = [];
};


/*******************************************/

/***** ModelResult *************************/

/**
* Represents result for the SEGInformedVar model
* @constructor
* @param {zygotine.M.SEGInformedVarModel} model 
*
* @property {zygotine.M.Chain[]} chains;
* @property {string} mcmcParameters
* @property {string} measureList
* @property {zygotine.M.Messages} messages
* @property {zygotine.M.SEGInformedVarModel} model 
* @property {string} specificParameters
* @property {number} totalNumberOfIterations
* @property {number} elapsedRunTime
*/
zygotine.M.ModelResult = function (model) {
    // Un ModelResult doit être instancié dans la méthode compute du Model.
    // Il faut lorsqu'on sort de compute mettre à jour la propriété hasError.
    // 
    zygotine.M.ErrorLogging.call(this);
    this.modelClassName = model.className;
    this.specificParameters = model.specificParameters.toString();
    this.logN = model.specificParameters.logN;
    this.mcmcParameters = model.mcmcParameters.toString();
    this.measureList = model.measureList.toString();
    this.chains = {};
    this.chainByType = { burnin: [], sample: [] };
    this.totalNumberOfIterations = 0;
    this.prngSeed = NaN;
    this.elapsedRunTime = NaN;
    this.quantileProbs = [.025, .05, .1, .25, .5, .75, .9, .95, .975];
};

zygotine.M.ModelResult.prototype = Object.create(zygotine.M.ErrorLogging.prototype);
/** function
* @param {string[]} names - les prefixes des noms de chaines, par ex. mu ou sd.
* @param {number} burninMonitoringLength - doit être 0, si le monitorBurnin est faux.
* @param {number} nIter - la taille de l'échantillon retenue
*/
zygotine.M.ModelResult.prototype.initChains = function (namePrefixes, burninMonitoringLength, nIter) {
    let ln = namePrefixes.length;
    var label, prefix, arr;
    for (let i = 0; i < ln; i++) {
        prefix = namePrefixes[i];
        label = prefix + "Burnin";
        arr = new Array(burninMonitoringLength);
        arr.fill(NaN);
        this.chains[label] = new zygotine.M.Chain(label, arr);
        this.chainByType.burnin.push(this.chains[label]);
        label = prefix + "Sample";
        arr = new Array(nIter);
        arr.fill(NaN);
        this.chains[label] = new zygotine.M.Chain(label, arr);
        this.chainByType.sample.push(this.chains[label]);
    }
};

zygotine.M.ModelResult.prototype.addLogOelToMuChains = function (oel) {
    var wPrefixes = this.workerChainLabelPrefixes;
    var logOel = Math.log(oel);
    var chain;

    chain = this.chains["muBurnin"].data;
    for (let ic = 0; ic < chain.length; ic++) {
        chain[ic] += logOel;
    }

    chain = this.chains["muSample"].data;
    for (let ic = 0; ic < chain.length; ic++) {
        chain[ic] += logOel;
    }
};

zygotine.M.ModelResult.prototype.computeQuantiles = function (probs) {
    var q = new zygotine.S.Quantile(typeof probs !== 'undefined' ? probs : this._quantileProbs);
    var burnin = this.chainByType.burnin;
    for (let i = 0; i < burnin.length; i++) {
        if (burnin[i].data.length === 0) {
            break;
        }

        burnin[i].quantile = q.compute(burnin[i].data);
    }

    var sample = this.chainByType.sample;
    for (let i = 0; i < burnin.length; i++) {
        sample[i].quantile = q.compute(sample[i].data);
    }
};

zygotine.M.ModelResult.prototype.getQuantileArray = function () {
    var rep = [].concat.apply([], this.chainByType.burnin.map(o => o.quantile));
    var rep2 = [].concat.apply([], this.chainByType.sample.map(o => o.quantile));
    rep = rep.concat(rep2);
    return rep;
};

/*******************************************/

/***** BaseModel ***************************/

/**
* Represents a model
* @constructor
* @param {zygotine.M.MeasureList} measureList - la liste de mesures 
* @param {zygotine.M.SEGInformedVarModelParameters} specificParameters - les paramètres spécifiques au modèle
* @param {zygotine.M.McmcParameters} mcmcParameters - MCMC parameters
*
* @property {string} className - at the prototype level
* @property {boolean} hasError - inherits from zygotine.M.ErrorLogging
* @property {zygotine.M.McmcParameters} mcmcParameters - MCMC parameters
* @property {zygotine.M.MeasureList} measureList - list of measures
* @property {zygotine.M.Messages} messages - inherits from zygotine.M.ErrorLogging
* @property {zygotine.M.PastDataSummary} pastData - past data summary
* @property {zygotine.M.SEGInformedVarModelParameters} specificParameters - parameters that are specific to the SEGInformedVarModel
*/
zygotine.M.BaseModel =
    function (
        measureList,
        specificParameters,
        mcmcParameters,
        pds) {

        this.className = "BaseModel";
        zygotine.M.ErrorLogging.call(this);
        if ((typeof measureList === 'undefined') || (typeof specificParameters === 'undefined') || (typeof mcmcParameters === 'undefined')) {
            this.addError(this.paramMissing);
            return;
        }
        this.specificParameters = specificParameters;
        this.mcmcParameters = mcmcParameters;
        this.measureList = measureList;
        if (!this.measureList.hasError) {
            if (this.specificParameters.logN) {
                this.measureList = zygotine.M.MeasureList.divideByOel(this.measureList, this.specificParameters.oel);
            }
        } else {
            this.addError("Syntax error: parameter measureList is not valid.");
            return;
        }

        this.logN = this.specificParameters.logN;
        this.oel = this.specificParameters.oel;
    };

zygotine.M.BaseModel.prototype = Object.create(zygotine.M.ErrorLogging.prototype);

zygotine.M.BaseModel.prototype.validateParameters__ = function (result) {
    // Il se peut que nous ayons des erreurs aient été confirmées lors de l'analyse et la conversion des chaines servant à
    // définir les mesures.
    if (this.hasError) {
        result.addError("Model contains errors.");
        /** @type {zygotine.M.Message[]} msgList} */
        var msgList = this.messages.msgList;
        /** @type {zygotine.M.Message} msgList} */
        var message;
        for (let i = 0; i < msgList.length; i++) {
            /** @type {zygotine.M.Message} msgList} */
            message = msgList[i];
            if (message.level === 'error') {
                result.addError(message.msgText);
            }
        }

        result.addInfo(zU.fmt("{0}: execution aborted. No further validation will be performed.", this.className));
        return;
    }

    zygotine.M.validateModelParameters(result, this);
};

zygotine.M.BaseModel.prototype.doParameterValidation = function () {
    var
        prngSeed = 1,
        validationOnly = true;
    return this.run(prngSeed, validationOnly);
};

zygotine.M.BaseModel.prototype.compute = function (prngSeed) {
    if (typeof prngSeed === 'undefined') {
      if ( zygotine.X.common.dataEntries.prngSeed.currentValue.length > 0 ) {
        prngSeed = parseInt(zygotine.X.common.dataEntries.prngSeed.currentValue)
      } else {
        prngSeed = zygotine.X.genPseudoRand32Bit()
      }
    }

    var validationOnly = false;
    return this.run(prngSeed, validationOnly);
};

