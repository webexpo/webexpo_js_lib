/* eslint 
    valid-jsdoc: 0
    no-extra-parens: 0
*/

/***** SEGInformedVarModelParameters *******/

/**
* Represents parameters for the SEGInformedVarModel
* @constructor
* @param {boolean} logN - true for lognormal distribution, false for normal distribution
*
* @property {number} muLower - lower bound for mu
* @property {number} muUpper - upper bound for mu
* @property {number} logSigmaMu 
* @property {number} logSigmaPrec
* @property {number} initMu
* @property {number} initSigma
*/
zygotine.M.SEGInformedVarModelParameters = function (logN, oel) {
    //    zygotine.M.ErrorLogging.call(this);
    zygotine.M.ModelParameters.call(this, logN, oel);
    this.muLower = -20;
    this.muUpper = 20;
    this.logSigmaMu = -0.1744;
    this.logSigmaPrec = 2.5523;
    this.initMu = -1.2039728043259361; // log(0.3) 
    this.initSigma = 0.91629073187415511; //log(2.5)

    if (!logN) {
        this.muLower = 40;
        this.muUpper = 125;
        this.logSigmaMu = 1.0986122886681098; //#(GM = 3) // corrigé le 31 mai log(3)
        this.logSigmaPrec = 1.191059; //(GSD=2.5)
        this.initMu = 85;
        this.initSigma = 3;
    }
};

zygotine.M.SEGInformedVarModelParameters.prototype = Object.create(zygotine.M.ModelParameters.prototype);

/*******************************************/

/****** SEGInformedVarModelResult **********/

/**
* Represents result for the SEGInformedVar model
* @constructor
* @param {zygotine.M.SEGInformedVarModel} model - le modèle qui produira les résultats seront affichés
*
* @property {string} mcmcParameters
* @property {string} measureList
* @property {zygotine.M.Messages} messages
* @property {zygotine.M.SEGInformedVarModel} model 
* @property {zygotine.M.PastDataSummary} pastData 
* @property {string} specificParameters
* @property {number} totalNumberOfIterations
*/
zygotine.M.SEGInformedVarModelResult = function (model) {
    zygotine.M.ModelResult.call(this, model);
    this.pastDataSummary = model.pastData.toString();
};

zygotine.M.SEGInformedVarModelResult.prototype = Object.create(zygotine.M.ModelResult.prototype);
zygotine.M.SEGInformedVarModelResult.prototype.className = "SEGInformedVarModelResult";

/*******************************************/

/***** SEGInformedVarModel *****************/

/**
* Represents a model
* @constructor
* @param {zygotine.M.MeasureList} measureList - la liste des mesures
* @param {zygotine.M.SEGInformedVarModelParameters} specificParameters  - les paramètres propres au modèle
* @param {zygotine.M.McmcParameters} mcmcParameters - MCMC parameters
* @param {zygotine.M.PastDataSummary} pds - past data summary
*
* @property {boolean} hasError - inherits from zygotine.M.ErrorLogging
* @property {zygotine.M.McmcParameters} mcmcParameters - MCMC parameters
* @property {zygotine.M.MeasureList} measureList - list of measures
* @property {zygotine.M.Messages} messages - inherits from zygotine.M.ErrorLogging
* @property {zygotine.M.PastDataSummary} pastData - past data summary
* @property {zygotine.M.SEGInformedVarModelParameters} specificParameters - parameters that are specific to the SEGInformedVarModel
*/
zygotine.M.SEGInformedVarModel =
    function (
        measureList,
        specificParameters,
        mcmcParameters,
        pds) {

        zygotine.M.BaseModel.call(this, measureList, specificParameters, mcmcParameters);
        this.className = "SEGInformedVarModel";
        if (typeof pds === 'undefined') {
            pds = zygotine.M.PastDataSummary.dummyPDS;
        }

        if (pds.hasError) {
            this.addError("Syntax error: parameter pds is not valid.");
        }

        this.pastData = pds;
    };


zygotine.M.SEGInformedVarModel.prototype = Object.create(zygotine.M.BaseModel.prototype);

zygotine.M.SEGInformedVarModel.prototype.paramMissing = "SEGInformedVar model constructor: the following parameters are required:  measureList, specificParameters, mcmcParameters.";

zygotine.M.SEGInformedVarModel.prototype.MEAny = false;

/**
* méthode privée, ne devrait être appelée que de doParameterValidation ou de compute.
* les deux paramètres doivent être présents et bien définits. Le premier est un entier et le second un booléen.
* @method
*
* @returns {SEGInformedVarModelResult}
*/
zygotine.M.SEGInformedVarModel.prototype.run = function (prngSeed, validationOnly) {
    /** @type {zygotine.M.SEGInformedVarModelResult} */
    var result = new zygotine.M.SEGInformedVarModelResult(this);
    this.validateParameters__(result); // Méthode de BaseModel

    if (validationOnly) {
        result.addInfo("Parameters validation performed.");
    }

    if (result.hasError || validationOnly) {
        return result;
    }

    // L'objectif n'était pas que de valider les paramètres.
    // De plus les paramètres sont valides.
    // On poursuit donc avec le calcul des chaines.



    result.prngSeed = prngSeed;
    zygotine.S.prng.init(result.prngSeed);

    /** @type {zygotine.M.McmcParameters} */
    var mcmc = this.mcmcParameters;
    var // mcmc
        monitorBurnin = mcmc.monitorBurnin,
        nBurnin = mcmc.nBurnin,
        nIter = mcmc.nIter,
        nThin = 1; // mcmc.nThin;

    result.initChains(["mu", "sd"], mcmc.monitorBurnin ? mcmc.nBurnin : 0, nIter);

    var
        burninMu = result.chains.muBurnin.data,
        burninSd = result.chains.sdBurnin.data,
        sampleMu = result.chains.muSample.data,
        sampleSd = result.chains.sdSample.data;

    /** @type {zygotine.M.SEGInformedVarModelParameters} */
    var sp = this.specificParameters;
    var logN = sp.logN;
    var zO = zygotine.O;
    
    var data = new zygotine.M.DataSummary(this.measureList, sp.logN);
    /** @type {zygotine.O.YGen} */
    var genY = new zO.YGen(sp.logN);

    var
        combinedN = data.n,
        mu = sp.initMu,
        muLim = [sp.muLower, sp.muUpper],
        sigma = sp.initSigma,
        logSigmaMu = sp.logSigmaMu,
        logSigmaSd = 1.0 / Math.sqrt(sp.logSigmaPrec); // ok pour 0.16 de PatrickB


    if (this.pastData.defined) {
        combinedN += this.pastData.n;
    }

    var
        muCondMean,
        muCondSd,
        p,
        pLim = [],
        sigmaBeta,
        u,
        yBar;

    var // int
        iter = 0,
        nIterations = nBurnin + (nIter * nThin),
        savedIter = 0;

    for (iter = 0; iter < nIterations; iter++) {

        if (data.anyCensored) {
            genY = zygotine.O.YGen.inits(data, mu, sigma, logN);
        }

        sigmaBeta = (data.uncensoredSum2 + genY.sum2 - 2.0 * mu * (data.uncensoredSum + genY.sum) + data.n * Math.pow(mu, 2.0)) / 2.0;
        if (this.pastData.defined) {
            sigmaBeta += this.pastData.n / 2.0 * Math.pow(this.pastData.mean - mu, 2) + this.pastData.ns2 / 2.0;
        }

        u = zS.uniform.sample(1, 0, 1);
        sigma = zygotine.O.sigmaGenICdf4InformedVar(combinedN, sigmaBeta, logSigmaMu, logSigmaSd, u);

        muCondMean = this.pastData.defined ?
            (data.uncensoredSum + genY.sum + this.pastData.n * this.pastData.mean) / combinedN :
            (data.uncensoredSum + genY.sum) / data.n;

        muCondSd = sigma / Math.sqrt(combinedN); // mu.cond.sd <- sigma/sqrt(data$N)
        pLim = [zygotine.S.normal.cdf((muLim[0] - muCondMean) / muCondSd, 0, 1, true, false),
        zygotine.S.normal.cdf((muLim[1] - muCondMean) / muCondSd, 0, 1, true, false)]; //  pnorm((mu.lim - mu.cond.mean) / mu.cond.sd)
        p = zS.uniform.sample(1, pLim[0], pLim[1]);
        mu = zS.normal.icdf(p, muCondMean, muCondSd, true, false);

        if (iter < nBurnin) {
            if (monitorBurnin) {
                burninMu[iter] = mu;
                burninSd[iter] = sigma;
            }
        }
        else if ((iter - nBurnin) % nThin === 0) {
            sampleMu[savedIter] = mu;
            sampleSd[savedIter] = sigma;
            savedIter++;
        }
    }// for( int iter = 1 ...

    if (this.logN) {
        result.adjustMuChains(this.oel)
    }

    result.totalNumberOfIterations = iter;
    return result;
};

/*******************************************/

/** SEGUninformativeModelParameters ********/

/**
* Represents parameters for the SEGInformedVar model
* @constructor
* @param {boolean} logN : true for the lognormal distribution and false for the normal distribution
* @property {number} muLower
* @property {number} muUpper
* @property {number} initMu
* @property {number} initSd
* @property {number[]} sdRange
*/
zygotine.M.SEGUninformativeModelParameters = function (logN, oel) {
    zygotine.M.ModelParameters.call(this, logN, oel);
    if (logN) {
        // this.logN = logN; // modif février 2020 instruction inutile
        this.muLower = -20;
        this.muUpper = 20;
        this.initMu = -1.2039728043259361; // log(0.3) 
        this.initSd = 0.91629073187415511;  // log(2.5)
        this.sdRange = [0.095, 2.3];
    } else {
        this.muLower = 40;
        this.muUpper = 125;
        this.initMu = 85;
        this.initSd = 3;
        this.sdRange = [0.1, 20];
    }
};


zygotine.M.SEGUninformativeModelParameters.prototype = Object.create(zygotine.M.ModelParameters.prototype);

/*******************************************/

/****** SEGUninformativeModelResult *******/

/**
* Represents result for the SEGUninformativeModel model
* @constructor
* @param {zygotine.M.SEGUninformativeModel} model - le modèle auquel sont rattachés les résultats
*
* @property {string} mcmcParameters
* @property {string} measureList
* @property {zygotine.M.Messages} messages
* @property {zygotine.M.SEGUninformativeModel} model // modif février 2020 : le type du modèle était incorrect.
* @property {zygotine.M.PastDataSummary} pastData 
* @property {string} specificParameters
* @property {number} totalNumberOfIterations
*/
zygotine.M.SEGUninformativeModelResult = function (model) {
    zygotine.M.ModelResult.call(this, model);
};

zygotine.M.SEGUninformativeModelResult.prototype = Object.create(zygotine.M.ModelResult.prototype);
zygotine.M.SEGUninformativeModelResult.prototype.className = "SEGUninformativeModelResult";

/*******************************************/

/***** SEGUninformativeModel ***************/

/**
* Represents a model
* @constructor
* @param {zygotine.M.MeasureList} measureList - une liste de mesures
* @param {zygotine.M.SEGInformedVarModelParameters} specificParameters - des paramètres propres au modèle en balisant l'exécution
* @param {zygotine.M.McmcParameters} mcmcParameters - les paramètres MCMC 
*
* @property {string} className - le nom de la classe, hérité de BaseModel
* @property {boolean} hasError - vrai si une erreur a été détecté, sinon faux, hérité de LoggingError via BaseModel
* @property {zygotine.M.McmcParameters} mcmcParameters = les paramètres MCMC, hérités de BaseModel
* @property {zygotine.M.MeasureList} measureList - la liste de mesures utilisée, héritée de BaseModel
* @property {zygotine.M.Messages} messages - une liste de messages d'erreurs ou de mises en garde, hérités de LoggingError via BaseModel
* @property {zygotine.M.SEGInformedVarModelParameters} specificParameters  - les paramètres propres à l'exécution du modèle, hérités de BaseModel
*/
zygotine.M.SEGUninformativeModel =
    function (
        measureList,
        specificParameters,
        mcmcParameters) {
        zygotine.M.BaseModel.call(this, measureList, specificParameters, mcmcParameters);
    };

zygotine.M.SEGUninformativeModel.prototype = Object.create(zygotine.M.BaseModel.prototype);
zygotine.M.SEGUninformativeModel.prototype.className = "SEGUninformativeModel";
zygotine.M.SEGUninformativeModel.prototype.paramMissing = "SEGUninformativeModel constructor: the following parameters are required:  measureList, specificParameters, mcmcParameters.";
zygotine.M.SEGUninformativeModel.prototype.MEAny = false;

// modif février 2020 ajout du commentaire qui suit
/**
* méthode privée, ne devrait être appelée que de doParameterValidation ou de compute.
* les deux paramètres doivent être présents et bien définis. Le premier est un entier et le second un booléen.
* @method
*
* @returns {SEGUninformativeModelResult}
*/
zygotine.M.SEGUninformativeModel.prototype.run = function (prngSeed, validationOnly) {

    /** @type {zygotine.M.SEGUninformativeModelParameters} */
    var result = new zygotine.M.SEGUninformativeModelResult(this); // modif février 2020 : on a remplacé SEGUninformativeModelParameters par SEGUninformativeModelResult
    this.validateParameters__(result); // Méthode de BaseModel

    if (validationOnly) {
        result.addInfo("Parameters validation performed.");
    }

    if (result.hasError || validationOnly) {
        return result;
    }

    //les paramètres d'entrée sont valides. On sauvegarde la racine du générateur de nbr aléatoires.

    // modif février 2020, 3 lignes suivantes mise en commentaires, la responsabilité des seed relève de la méthode compute.
    /*
    if (typeof prngSeed === 'undefined') {
        prngSeed = (Math.random() * Math.pow(2, 31)) | 0;
    }
    */

    zygotine.S.prng.init(prngSeed);
    result.prngSeed = prngSeed;

    /** @type {zygotine.M.McmcParameters} */
    var mcmc = this.mcmcParameters;
    var // mcmc
        nIter = mcmc.nIter,
        nThin = 1; //mcmc.nThin,
    nBurnin = mcmc.nBurnin,
        monitorBurnin = mcmc.nBurnin;

    result.initChains(["mu", "sd"], mcmc.monitorBurnin ? mcmc.nBurnin : 0, nIter);
    var
        burninMu = result.chains.muBurnin.data,
        burninSd = result.chains.sdBurnin.data,
        sampleMu = result.chains.muSample.data,
        sampleSd = result.chains.sdSample.data;

    var
        pnorm = function (x) {
            return zygotine.S.normal.cdf(x, 0, 1, true, false);
        },
        iter = 0,
        nIterations = nBurnin + (nIter * nThin),
        savedIter = 0;

    /** @type {zygotine.M.SEGUninformativeModelParameters} */
    var sp = this.specificParameters;
    var logN = sp.logN;
    /** @type {zygotine.M.DataSummary} */
    var data = new zygotine.M.DataSummary(this.measureList, logN);
    /** @type {zygotine.O.YGen} */
    var genY = new zygotine.O.YGen(sp.logN);

    var
        mu = sp.initMu,
        muLim = [sp.muLower, sp.muUpper],
        sdRange = [sp.sdRange[0], sp.sdRange[1]],
        tauCondAlpha = (data.n - 1.0) / 2.0, // modif février 2020: on a remplacé "data.n / 2 - 1" par "(data.n - 1.0) / 2.0"
        // tauCondAlpha = (data.n  / 2.0) - 1.0, 
        tauBeta = 0.0,
        tauRange = [1.0 / (sp.sdRange[1] * sp.sdRange[1]), 1.0 / (sp.sdRange[0] * sp.sdRange[0])],
        tau = 1.0 / (sp.initSd * sp.initSd);



    var
        muCondMean,
        muCondSd,
        p,
        pLim = [],
        tauCondBeta;

    for (iter = 0; iter < nIterations; iter++) {

        if (data.anyCensored) {
            genY = zygotine.O.YGen.inits(data, mu, 1 / Math.sqrt(tau), logN);
            muCondMean = (data.uncensoredSum + genY.sum) / data.n;
        } else {
            muCondMean = zU.mean(data.y);
        }

        tauCondBeta = tauBeta + (data.uncensoredSum2 + genY.sum2 - 2 * mu * (data.uncensoredSum + genY.sum) + data.n * mu * mu) / 2.0;
        if (tauCondBeta < 1e-6) {
            tauCondBeta = 0.0; // protection against numeric imprecision
        }

        if (tauCondBeta === 0.0) {
            let tmpRange = this.sdRange;
            if (tmpRange[0] === 0.0) {
                tmpRange[0] = 0.001;
            }

            let sigma = zO.randomPow(data.n, tmpRange);
            tau = 1.0 / (sigma * sigma);
        } else {
            // sample from a truncated gamma distribution  
            tau = zygotine.O.rGammaTruncated(tauCondAlpha, tauCondBeta, tauRange);
        }

        // Sample from f(mu | tau)
        muCondSd = 1 / Math.sqrt(data.n * tau);
        plim = [pnorm((muLim[0] - muCondMean) / muCondSd), pnorm((muLim[1] - muCondMean) / muCondSd)];
        p = zygotine.S.uniform.sample(1, plim[0], plim[1]);
        mu = zygotine.S.normal.icdf(p, muCondMean, muCondSd, true, false);
        if (iter < nBurnin) {
            if (monitorBurnin) {
                burninMu[iter] = mu;
                burninSd[iter] = 1 / Math.sqrt(tau);
            }
        }
        else if ((iter - nBurnin) % nThin === 0) {
            sampleMu[savedIter] = mu;
            sampleSd[savedIter] = 1 / Math.sqrt(tau);
            savedIter++;
        }
    } // for (iter = 0; iter < nIterations; iter++) 

    if (this.logN) {
        result.adjustMuChains(this.oel)
    }

    result.totalNumberOfIterations = iter;
    return result;
}; // compute

/*******************************************/

/***** BetweenWorkerModelParameters ********/

/**
* Represents parameters for a model (BetweenWorkerModel)
* @constructor
* @param {boolean} logNormDstrn - lorsque faux, on considère faire affaire avec la loi normale.
* @param {boolean} uupOnSds  - lorsque vrai, on utilise la loi normale pour les priors des écarts-type
* @param {number} oel - la valeur d'exposition limite  
*
* @property {number} initMuOverall - valeur initiale de muOverall
* @property {number} initSigmaWithin - valeur initiale de sigmaWithin
* @property {number} muOverallLower - borne inférieure pour la prior 'uniforme' de muOverall
* @property {number} muOverallUppper - borne supérieure pour la prior 'uniforme' de muOverall
* @property {number[]} sigmaBetweenRange - lorsque uupOnSds est vrai, définit l'interval balisant la prior 'uniforme' de sigmaBetween
* @property {number[]} sigmaWithinRange - lorsque uupOnSds est vrai, définit l'interval balisant la prior 'uniforme' de sigmaWithin
* @property {number} logSigmaWithinMu - lorsque uupOnSds est faux, paramètre pour la prior de sigmaWithin
* @property {number} logSigmaWithinPrec - lorsque uupOnSds est faux, paramètre pour la prior de sigmaWithin
* @property {number} logsigmaBetweenMu - lorsque uupOnSds est faux, paramètre pour la prior de sigmaBetween
* @property {number} logsigmaBetweenPrec - lorsque uupOnSds est faux, paramètre pour la prior de sigmaBetween

*/
zygotine.M.BetweenWorkerModelParameters = function (logN, uupOnSds, oel) {
    zygotine.M.ModelParameters.call(this, logN, oel);

    this.uupOnSds = uupOnSds;

    //les 2 premiers paramètres, balisant muOverall, sont présents quelle que soit la valeur de uupOnSds

    if (logN) {
        this.initMuOverall = -1.20397; // log(0.3)
        this.initSigmaWithin = 0.91629; // log(2.5)
        this.muOverallLower = -20;
        this.muOverallUpper = 20;
    } else {
        this.initMuOverall = 85;
        this.initSigmaWithin = 3;
        this.muOverallLower = 40;
        this.muOverallUpper = 125;
    }

    // les 6 paramètres suivants sont divisés en 2 groupes selon la valeur de uupOnSds
    if (uupOnSds) {
        if (logN) {
            this.sigmaBetweenRange = [0.0, 2.3]; 
            this.sigmaWithinRange = [0.095, 2.3];
        } else {
            this.sigmaBetweenRange = [0, 20];
            this.sigmaWithinRange = [0.1, 20];
        }
    } else {
        if (logN) {
            this.logSigmaBetweenMu = -0.8786;
            this.logSigmaBetweenPrec = 1.634;
            this.logSigmaWithinMu = -0.4106;
            this.logSigmaWithinPrec = 1.9002;
        } else {
            this.logSigmaBetweenMu = 1.098612;
            this.logSigmaBetweenPrec = 1.191059;
            this.logSigmaWithinMu = 1.098612;
            this.logSigmaWithinPrec = 1.191059;
        }
    }
};

zygotine.M.BetweenWorkerModelParameters.prototype = Object.create(zygotine.M.ModelParameters.prototype);

/*******************************************/

/****** BetweenWorkerModelResult **********/

/**
* Represents result for the BetweenWorkerModel model
* @constructor
* @param {zygotine.M.BetweenWorkerModel} model - an instance of BetweenWorkerModel for which the result is built.
*
* @property {string} mcmcParameters
* @property {string} measureList
* @property {zygotine.M.Messages} messages
* @property {zygotine.M.BetweenWorkerModel} model 
* @property {string} specificParameters
* @property {number} totalNumberOfIterations
*/
zygotine.M.BetweenWorkerModelResult = function (model) {
    zygotine.M.ModelResult.call(this, model);
    this.workerIds = [];
};

zygotine.M.BetweenWorkerModelResult.prototype = Object.create(zygotine.M.ModelResult.prototype);
zygotine.M.BetweenWorkerModelResult.prototype.className = "BetweenWorkerModelResult";

zygotine.M.BetweenWorkerModelResult.prototype.updateChains = function (muOverall, sigmaWithin, sigmaBetween, muWorker, burnin, position) {
    var wPrefixes = this.workerChainLabelPrefixes;
    var burninOrSample = burnin ? "Burnin" : "Sample";
    var burninMuWorker = [], sampleMuWorker = [];
    wPrefixes.map(function (prefix, i) { this.chains[prefix + burninOrSample].data[position] = muWorker[i]; }, this);
    this.chains["muOverall" + burninOrSample].data[position] = muOverall;
    this.chains["sigmaWithin" + burninOrSample].data[position] = sigmaWithin;
    this.chains["sigmaBetween" + burninOrSample].data[position] = sigmaBetween;
};

/*******************************************/

/***** BetweenWorkerModel ***************/

/**
* Represents a model
* @constructor
* @param {zygotine.M.MeasureList} measureList - une liste de mesures
* @param {zygotine.M.BetweenWorkerModelParameters} specificParameters - des paramètres propres au modèle qui en balisent l'exécution
* @param {zygotine.M.McmcParameters} mcmcParameters - les paramètres MCMC 
*
* @property {string} className - le nom de la classe, hérité de BaseModel
* @property {boolean} hasError - vrai si une erreur a été détecté, sinon faux, hérité de LoggingError via BaseModel
* @property {zygotine.M.McmcParameters} mcmcParameters = les paramètres MCMC, hérités de BaseModel
* @property {zygotine.M.MeasureList} measureList - la liste de mesures utilisée, héritée de BaseModel
* @property {zygotine.M.Messages} messages - une liste de messages d'erreurs, de messages d'information ou de mises en garde, héritée de LoggingError via BaseModel
* @property {zygotine.M.BetweenWorkerModelParameters} specificParameters  - les paramètres propres à l'exécution du modèle, hérités de BaseModel
*/
zygotine.M.BetweenWorkerModel =
    function (
        measureList,
        specificParameters,
        mcmcParameters) {

        zygotine.M.BaseModel.call(this, measureList, specificParameters, mcmcParameters);
        this.className = "BetweenWorkerModel";
    };

zygotine.M.BetweenWorkerModel.prototype = Object.create(zygotine.M.BaseModel.prototype);
/* modif février 2020. */
zygotine.M.BetweenWorkerModel.prototype.betaMin = 1.0E-100;
zygotine.M.BetweenWorkerModel.prototype.paramMissing = "BetweenWorkerModel constructor: the following parameters are required:  measureList, specificParameters, mcmcParameters.";
/* modif février 2020. méthode calcSigma : nouvelle version plus bas
zygotine.M.BetweenWorkerModel.prototype.calcSigma = function (n, sigmaRange, u) {
    if (typeof u === 'undefined') {
        u = zygotine.S.uniform.sample(1, 0, 1);
    }

    let a = 1 - n;
    let f1 = [Math.pow(sigmaRange[0], a), Math.pow(sigmaRange[0], a)];
    a = 1.0 / a;
    let f2 = [Math.pow(1 - u, a), Math.pow(u, a)];
    return zygotine.U.sumProduct(f1, f2);
};
*/

// modif février 2020. nouvelle version de la méthode.
zygotine.M.BetweenWorkerModel.prototype.calcSigma = function (n, sigmaRange, u) {
    if (typeof u === 'undefined') {
        u = zygotine.S.uniform.sample(1, 0, 1);
    }

    //sigma.between < - sum((sigma.between.range ^ (1 - a)) * c(1 - u, u)) ^ (1 / (1 - a))
    let expon = 1 - n;
    let f1 = zygotine.U.sum([Math.pow(sigmaRange[0], expon) * (1 - u), Math.pow(sigmaRange[1], expon) * u]);
    expon = 1.0 / expon;
    let f2 = Math.pow(f1, expon);
    return f2;
};




zygotine.M.BetweenWorkerModel.prototype.run = function (prngSeed, validationOnly) {
    var result = new zygotine.M.BetweenWorkerModelResult(this);
    this.validateParameters__(result); // Méthode de BaseModel
    if (validationOnly) {
        result.addInfo("Parameters validation performed.");
    }

    if (result.hasError || validationOnly) {
        return result;
    }


    // L'objectif n'était pas que de valider les paramètres.
    // De plus les paramètres sont valides.
    // On poursuit donc avec le calcul des chaines.

    result.prngSeed = prngSeed;
    zygotine.S.prng.init(result.prngSeed);

    /** @type {zygotine.M.McmcParameters} */
    var mcmc = this.mcmcParameters;
    var // mcmc
        nIter = mcmc.nIter,
        nThin = 1; // mcmc.nThin,
    nBurnin = mcmc.nBurnin,
        monitorBurnin = mcmc.monitorBurnin;

    /** @type {zygotine.M.BetweenWorkerModelParameters} */
    var sp = this.specificParameters;
    var uupOnSds = sp.uupOnSds;
    var logN = sp.logN;
    /** @type {zygotine.M.DataSummary} */
    var data = new zygotine.M.DataSummary(this.measureList, logN);
    /** @type {zygotine.M.WorkerDigest} */
    var workerDigest = new zygotine.M.WorkerDigest(data.measureList);
    result.workerIds = workerDigest.workerIds;

    //On initialise les chainesm, burnin et sample
    result.workerChainLabelPrefixes = workerDigest.workerIds.map(function (w) { return "mu_" + w; }); // workerDigest.workerIds est ordonné alphabétiquement
    var chainLabelPrefixes = ["muOverall", "sigmaWithin", "sigmaBetween"].concat(result.workerChainLabelPrefixes);
    result.initChains(chainLabelPrefixes, mcmc.monitorBurnin ? mcmc.nBurnin : 0, nIter);
    result.workers = { ids: workerDigest.workerIds.slice(0), invertedIds: Object.assign({}, workerDigest.workerIdsInverted) };
    //les valeurs initiales
    var
        tmpObj = workerDigest.getInits();
    var
        muOverall = tmpObj.muOverall,
        muWorker = tmpObj.muWorker,
        sigmaWithin = tmpObj.sigmaWithin,
        sigmaBetween;
    // Suite a une demande, 
    // on utilise pour sigmaWithin et muOverall, plutôt que le résultat d'un calcul dans workerDigest, 
    // les valeurs de initMuOverall et initSigmaWithin incluses dans les paramètres(specificParameters).
    // Les paramètres initMuOverall et initSigmaWithin peuvent être ajustées par l'usager a travers le UI.
    sigmaWithin = this.specificParameters.initSigmaWithin;
    muOverall = this.specificParameters.initMuOverall;


    workerDigest.updateMuValues(muOverall, muWorker);
    // Averages contient une moyenne d'ensemble des mesures (avg) et un tableau de moyennes par travailleur (workerAvg).
    // En autant qu'il y ait des données censurées, on mettra à jour ces infos à chaque itération
    // sinon les valeurs calculées ci-bas demeureront inchangées.
    var
        averages = workerDigest.getAverages();

    /** @type {zygotine.O.YGen} */
    var genY = new zygotine.O.YGen(sp.logN);


    // Rapplelons que logSigmaWithinSd et logSigmaBetweenSd ne seront utilisés que si uupOnSds est faux!
    var
        logSigmaWithinSd = NaN,
        logSigmaBetweenSd = NaN;

    if (!sp.uupOnSds) {
        logSigmaBetweenSd = 1 / Math.sqrt(sp.logSigmaBetweenPrec);
        logSigmaWithinSd = 1 / Math.sqrt(sp.logSigmaWithinPrec);
    }

    var
        beta, // b dans le code R de PB
        iter = 0,
        nIterations = nBurnin + (nIter * nThin),
        savedIter = 0;

    // 2 variables qui ne seront utilisées ssi uupOnSds est vrai, autremenent elles demeureront de type undefined.
    var tauWithinRange, tauBetweenRange;

    if (uupOnSds) {
        tauWithinRange = [1.0 / (sp.sigmaWithinRange[1] * sp.sigmaWithinRange[1]), 1.0 / (sp.sigmaWithinRange[0] * sp.sigmaWithinRange[0])];
        tauBetweenRange = [1.0 / (sp.sigmaBetweenRange[1] * sp.sigmaBetweenRange[1]), 1.0 / (sp.sigmaBetweenRange[0] * sp.sigmaBetweenRange[0])]; // il se peut que  sp.sigmaWihtinRange soit 0, alors la borne supérieure sera Infinity.
    }

    var // des doubles
        tmpMean,
        tmpSd,
        sw2,
        sb2;
    var // des tableaux de doubles
        muA,
        sigma2A,
        tmpArr,
        tmpVariance;

    for (iter = 0; iter < nIterations; iter++) {

        if (data.anyCensored) {
            workerDigest.generateValues(muOverall, muWorker, sigmaWithin, data);
            //On prend les moyennes qui varient en compte les valeurs générées associées aux mesures censurées
            averages = workerDigest.getAverages();
        }

        let tmp = workerDigest.getBeta(muOverall, muWorker);
        beta = tmp.beta;

        // Sample from f(sigma.within | other parms)
        if (uupOnSds) {
            if (beta < this.betaMin) {
                let u = zygotine.S.uniform.sample(1, 0, 1);
                sigmaWithin = this.calcSigma(data.n, sp.sigmaWithinRange, u);
            }
            else {
                let tau = zygotine.O.rGammaTruncated((data.n - 1.0) / 2.0, beta, tauWithinRange); // modif février 2020 : on a remplacé "(data.n / 2) - 1" "(data.n - 1.0) / 2.0"
                //let tau = zygotine.O.rGammaTruncated((data.n / 2.0) - 1.0, beta, tauWithinRange); // modif février 2020: mettre en commentaire
                sigmaWithin = 1.0 / Math.sqrt(tau);
            }
        } else {
            sigmaWithin = zygotine.O.sigmaGenICdf4InformedVar(data.n, beta, sp.logSigmaWithinMu, logSigmaWithinSd);
        } //if (uupOnSds)

        // Sample from f(sigma.between | other parms)
        beta = zygotine.U.sumSqr(muWorker) / 2.0;
        if (uupOnSds) {
            if (beta < this.betaMin) {
                let u = zygotine.S.uniform.sample(1, 0, 1);
                sigmaBetween = this.calcSigma(workerDigest.nWorkers, sp.sigmaBetweenRange, u);
            }
            else {
                let tau = zygotine.O.rGammaTruncated((workerDigest.nWorkers - 1.0) / 2.0, beta, tauBetweenRange); // modif février 2020 : on a remplacé "(data.nWorkers / 2) - 1" "(data.nWorkers - 1.0) / 2.0"
                //let tau = zygotine.O.rGammaTruncated((workerDigest.nWorkers / 2.0) - 1.0, beta, tauBetweenRange) // modif février 2020: mettre en commentaire
                sigmaBetween = 1.0 / Math.sqrt(tau);
            }
        } else {
            sigmaBetween = zygotine.O.sigmaGenICdf4InformedVar(workerDigest.nWorkers, beta, sp.logSigmaBetweenMu, logSigmaBetweenSd);
        }



        // Sample from from f(mu.overall | other parms)
        tmpMean = averages.avg - zU.sumProduct(workerDigest.measureCountByWorker, muWorker) / data.n;
        tmpSd = sigmaWithin / Math.sqrt(data.n);
        muOverall = zygotine.O.rNormCensored(tmpMean, tmpSd, [sp.muOverallLower], [sp.muOverallUpper], false)[0];
        muA = zygotine.U.substract(averages.workerAvg, muOverall);

        // Sample from f(mu.worker's | other parms)
        sw2 = sigmaWithin * sigmaWithin;
        sigma2A = workerDigest.measureCountByWorker.map(function (n) { return sw2 / n; });
        sb2 = sigmaBetween * sigmaBetween;
        tmpMean = muA.map(function (x, i) { return (x * sb2) / (sigma2A[i] + sb2); });
        tmpVar = sigma2A.map(function (x) { return (x * sb2) / (x + sb2); });
        muWorker = tmpMean.map(function (mean, i) { return zygotine.S.normal.sample(1, mean, Math.sqrt(tmpVar[i])); });
        workerDigest.updateMuValues(muOverall, muWorker);

        if (iter < nBurnin) {
            if (monitorBurnin) {
                result.updateChains(muOverall, sigmaWithin, sigmaBetween, muWorker, true, iter);
            }
        } else {
            if ((iter - nBurnin) % nThin === 0) {
                result.updateChains(muOverall, sigmaWithin, sigmaBetween, muWorker, false, savedIter);
                savedIter++;
            }
        }
    } // for (iter = 0; iter < nIterations; iter++) 

    if (this.logN) {
        result.adjustMuChains(this.oel)
    }

    result.totalNumberOfIterations = iter;
    return result;
}; // compute

/*******************************************/

/***** BetweenWorkerModelParameters *******/

/**
* Représente les paramètres de départ du modèle BetweenWorkerModel
* @constructor
* @property {boolean} logN - un booléen, vrai pour la loi log-normale et faux pour la loi normale
* @property {number} logSigmaBetweenMu - un nombre, pertinent ssi uupOnSds est faux
* @property {number} logSigmaBetweenPrec - un nombre, pertinent ssi uupOnSds est faux
* @property {number} logSigmaWithinMu - un nombre, pertinent ssi uupOnSds est faux
* @property {number} logSigmaWithinPrec - un nombre, pertinent ssi uupOnSds est faux
* @property {number} muOverallLower - paramètre commun au 2 types de priors possibles (uupOnSds)
* @property {number} muOverallUpper - paramètre commun au 2 types de priors possibles (uupOnSds)
* @property {number[]} sigmaBetweenRange - un interval, pertinent ssi uupOnSds est vrai
* @property {number[]} sigmaWithinRange - un interval, pertinent ssi uupOnSds est vrai
* @property {boolean} uupOnSds -  un booléen, vrai si l'on utilise une prior uniforme, faux sinon ...
*/


/***** Worker ******************************/

zygotine.M.Worker = function (id) {
    this._mean = NaN; // champ privé, facilite le calcul de mu (dans getInits du WorkerDigest), sera mis à NaN lorsque devenu inutile.
    this.mu = NaN;
};

/*******************************************/

/***** WorkerDigest ************************/

/**
* Represents information on workers organized for the BetweenWorker model.
* @constructor
* @param {zygotine.M.Measure[]} measureArray - un tableau de mesures qui sont dans l'ordre original ...
*
* @property {string} className - le nom de la classe, hérité de BaseModel
* @property {boolean} hasError - vrai si une erreur a été détecté, sinon faux, hérité de LoggingError via BaseModel
* @property {zygotine.M.McmcParameters} mcmcParameters = les paramètres MCMC, hérités de BaseModel
* @property {zygotine.M.MeasureList} measureList - la liste de mesures utilisée, héritée de BaseModel
* @property {zygotine.M.Messages} messages - une liste de messages d'erreurs, de messages d'information ou de mises en garde, héritée de LoggingError via BaseModel
* @property {number} nWorkers - le nombre de travailleurs
* @property {zygotine.M.BetweenWorkerModelParameters} specificParameters  - les paramètres propres à l'exécution du modèle, hérités de BaseModel
*/
zygotine.M.WorkerDigest = function (measureArray) {

    this._measureArray = []; // selon l'ordre de présentation, ie celui de la chaine de caractère ayant servi à définir la MeasureList, ordre que préserve measureArray ...
    this._measureByType = { uncensored: [], greaterThan: [], lessThan: [], interval: [] };
    this._measureByWorker = [];
    //    this.workerInfo = {};
    this._measureOrderedAsR = []; // l'ordre de R, ie "non-censuré" suivi de  ">" suivi de  "<" suivi de [ , ]! Dans chaque groupe l'ordre est donné par celui de la MeasureList
    var m;
    var o;

    for (let i = 0; i < measureArray.length; i++) {
        m = measureArray[i];
        extM = new zygotine.M.ExtendedMeasure(m);

        this._measureArray.push(extM); // l'ordre est celui du paramètre measureArray, donc l'ordre original des mesures.
        this._measureByType[extM.type].push(extM); // pour un type donné, l'ordre est celui du parametre measureArray.
        this._measureByWorker[extM.workerId] = this._measureByWorker[extM.workerId] || [];
        this._measureByWorker[extM.workerId].push(extM); // pour un travailleur, les mesures sont dans l'ordre de measureArray
    }

    this.workerIds = Object.keys(this._measureByWorker).sort(); // pour un travailleur, les mesures sont dans l'ordre de measureArray).sort(); // le tri est alphabétique.
    this.measureCountByWorker = this.workerIds.map(function (id) { return this._measureByWorker[id].length; }, this); // l'ordre sera celui des id triés alphabétiquement

    //this.measureCountByWorker = [];
    //for (let i = 0; i < this.workerIds.length; i ++) {
    //    this.measureCountByWorker.push(this._measureByWorker[this.workerIds[i]].length);
    //}

    this.workerIdsInverted = {};
    for (let i = 0; i < this.workerIds.length; i++) {
        this.workerIdsInverted[this.workerIds[i]] = i;
    }

    this.n = this._measureArray.length;
    this.nWorkers = this.workerIds.length;
    this._measureOrderedAsR = this._measureByType.uncensored.concat(this._measureByType.greaterThan, this._measureByType.lessThan, this._measureByType.interval);
    this.muOverall = NaN;
    this.muWorker = [];
};

zygotine.M.WorkerDigest.prototype = {

    updateMuValues: function (muOverall, muWorker) {
        this.muOverall = muOverall;
        this.muWorker = muWorker.slice(0); // on crée une copie du tableau.
    },

    getAverages: function () {
        var globalSum = 0.0;
        var globalN = 0;
        var wAverages = [];
        for (let i = 0; i < this.workerIds.length; i++) {
            let workerId = this.workerIds[i]; // id d'un travailleur
            let nMeasures = this._measureByWorker[workerId].length; // nombre de mesures associées à un travailleur (workerId)
            let tmpSum = zU.sum(this._measureByWorker[workerId].map(function (o) { return o.generatedValue; })); // somme des mesures pour le travailleur
            wAverages.push(tmpSum / nMeasures);
            globalSum += tmpSum;
            globalN += nMeasures;
        }

        var average = globalSum / globalN;
        return { avg: average, workerAvg: wAverages };
    },

    getInits: function () {
        var wIds = this.workerIds;
        var wId;
        var sum = 0;
        var workerInfo = []; // une entrée pour chaque travailleur
        for (let i = 0; i < wIds.length; i++) {
            wId = wIds[i];
            workerInfo[wId] = new zygotine.M.Worker();
            let currentWInfo = workerInfo[wId];
            currentWInfo._mean = zU.sum(this._measureByWorker[wId].map(function (o) { return o.getInitialValue(); })) / this._measureByWorker[wId].length;
            sum += currentWInfo._mean;
        }

        var inits = {
            muOverall: sum / wIds.length,
            muWorker: [],
            sigmaWithin: NaN
        };

        for (let i = 0; i < wIds.length; i++) {
            wId = wIds[i];
            let currentWInfo = workerInfo[wId];
            let cm = currentWInfo._mean - inits.muOverall;
            inits.muWorker.push(cm);
            currentWInfo.mu = cm;
        }

        var predicted = this._measureOrderedAsR.map(function (m) {
            return m.getInitialValue() - workerInfo[m.workerId]._mean;
        }, this);

        inits.sigmaWithin = Math.sqrt(zygotine.U.variance(predicted));
        return inits;
    },
    /** 
     * @function
     * @param {number} muOverall
     * @param {number[]} muWorker
     * @param {number} sigmaWithin
     * @param {zygotine.M.DataSummary} data
     * @returns {undefined}
     */
    generateValues: function (muOverall, muWorker, sigmaWithin, data) {
        var tmpMuWorker = {};
        for (let i = 0; i < muWorker.length; i++) {
            tmpMuWorker[this.workerIds[i]] = muWorker[i];
        }

        var getMeans = function (corpus) {
            let tmp = corpus.map(function (measure) {
                return tmpMuWorker[measure.workerId] + muOverall;
            }, this);
            return tmp;
        };

        var getGenFunction = function (type) {

            switch (type) {
                case 'lt':
                    return function (corpus, means) {
                        var tmp = [];
                        var len = corpus.length;
                        for (let i = 0; i < len; i++) {
                            let generated = zygotine.O.rNormCensored(means[i], sigmaWithin, [], [data.lt[i]], false)[0];
                            tmp.push(generated);
                            corpus[i].generatedValue = generated;
                        }

                        return tmp;
                    };

                case 'gt':
                    return function (corpus, means) {
                        var tmp = [];
                        var len = corpus.length;
                        if (len > 0) {
                            for (let i = 0; i < len; i++) {
                                let generated = zygotine.O.rNormCensored(means[i], sigmaWithin, [data.gt[i]], [], false)[0];
                                tmp.push(generated);
                                corpus[i].generatedValue = generated;
                            }
                        }

                        return tmp;
                    };

                case 'i':
                    return function (corpus, means) {
                        var tmp = [];
                        var len = corpus.length;
                        if (len > 0) {
                            for (let i = 0; i < len; i++) {
                                let generated = zygotine.O.rNormCensored(means[i], sigmaWithin, [data.i.lower[i]], [data.i.upper[i]], false)[0];
                                tmp.push(generated);
                                corpus[i].generatedValue = generated;
                            }
                        }

                        return tmp;
                    };
            }
        };

        var corpus, tmp;

        corpus = this._measureByType.greaterThan;
        tmp = getGenFunction("gt")(corpus, getMeans.call(this, corpus));

        corpus = this._measureByType.lessThan;
        //zygotine.S.prng.init(12);
        tmp = getGenFunction("lt")(corpus, getMeans.call(this, corpus));

        corpus = this._measureByType.interval;
        tmp = getGenFunction("i")(corpus, getMeans.call(this, corpus));
    },

    getBeta: function (muOverall, muWorker) {
        var residuals = this._measureOrderedAsR.map(
            function (m, i) { return m.generatedValue - muWorker[this.workerIdsInverted[m.workerId]] - muOverall; },
            this);
        var beta = zygotine.U.sumSqr(residuals) / 2.0;
        // le programme appelant n'a pas besoin de residuals, mais fournit un moyen de se comparer à R.
        return { residuals: residuals, beta: beta };
    }
};

/*******************************************/
