/***** SEGInformedMeanModelParameters *******/

/**
* Represents parameters for the SEGInformedMeanModel
* @constructor
* @param {boolean} logN - true for lognormal distribution, false for normal distribution
* @param {number} oel - occupational exposure limit
*
* @property {number} muMean
* @property {number} muSd 
* @property {number} logSigmaMu 
* @property {number} logSigmaPrec
* @property {number} initMu
* @property {number} initSigma
*/

zygotine.M.SEGInformedMeanModelParameters = function (logN, oel) {
    //    zygotine.M.ErrorLogging.call(this);
    zygotine.M.ModelParameters.call(this, logN, oel);
    this.muMean = -20;
    this.muSd = 20;
    this.logSigmaMu = -0.1744;
    this.logSigmaPrec = 2.5523;
    this.initSigma = 0.91629073187415511; //log(2.5)

    if (!logN) {
        this.muMean = 40;
        this.muSd = 125;
        this.logSigmaMu = 1.0986122886681098; //#(GM = 3) // corrigé le 31 mai log(3)
        this.logSigmaPrec = 1.191059; //(GSD=2.5)
        this.initSigma = 3;
    }
};

zygotine.M.SEGInformedMeanModelParameters.prototype = Object.create(zygotine.M.ModelParameters.prototype);
/**
 * @method
 * @param {number} mean
 */
zygotine.M.SEGInformedMeanModelParameters.prototype.setMuMean = function (muMean) {
    this.muMean = muMean;
};

zygotine.M.SEGInformedMeanModelParameters.prototype.setMuMean = function (muSd) {
    this.muSd = muSd;
};

/*******************************************/

/***** SEGInformedVarModel *****************/

/**
* Represents a model
* @constructor
* @param {zygotine.M.MeasureList} measureList - la liste des mesures
* @param {zygotine.M.SEGInformedVarModelParameters_informedMeanPrior} specificParameters  - les paramètres propres au modèle
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
zygotine.M.SEGInformedMeanModel =
    function (
        measureList,
        specificParameters,
        mcmcParameters,
        pds) {

        zygotine.M.BaseModel.call(this, measureList, specificParameters, mcmcParameters);
        this.className = "SEGInformedVarModelParameters_priorMeanSd";
        if (typeof pds === 'undefined') {
            pds = zygotine.M.PastDataSummary.dummyPDS;
        }

        if (pds.hasError) {
            this.addError("Syntax error: parameter pds is not valid.");
        }

        this.pastData = pds;
    };


zygotine.M.SEGInformedMeanModel.prototype = Object.create(zygotine.M.BaseModel.prototype);

zygotine.M.SEGInformedMeanModel.prototype.paramMissing = "SEGInformedMean model (using prior on mu defined by muMean and muSd) constructor: the following parameters are required:  measureList, specificParameters, mcmcParameters.";

zygotine.M.SEGInformedMeanModel.prototype.MEAny = false;

/**
* méthode privée, ne devrait être appelée que de doParameterValidation ou de compute.
* les deux paramètres doivent être présents et bien définits. Le premier est un entier et le second un booléen.
* @method
*
* @returns {SEGInformedVarModelResult}
*/
zygotine.M.SEGInformedMeanModel.prototype.run = function (prngSeed, validationOnly) {
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
    var
        sp = this.specificParameters,
        logN = sp.logN,
        data = new zygotine.M.DataSummary(this.measureList, sp.logN),
        genY = new zygotine.O.YGen(logN),
        combinedN = data.n + (this.pastData.defined ? this.pastData.n : 0),
        informedMeanPrior = { a: 1 / (sp.muSd * sp.muSd), b: sp.muMean / (sp.muSd * sp.muSd) },
        sigma = sp.initSigma,
        logSigmaMu = sp.logSigmaMu,
        logSigmaSd = 1.0 / Math.sqrt(sp.logSigmaPrec), // ok pour 0.16 de PatrickB
        mu = zygotine.S.normal.sample(1, sp.muMean, sp.muSd);

    var
        muCondMean,
        sigmaBeta,
        u;

    var // int
        iter = 0,
        nIterations = nBurnin + (nIter * nThin),
        savedIter = 0,
        tmpA,
        tmpB;


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

        tmpA = combinedN / (sigma * sigma) + informedMeanPrior.a;
        tmpB = (data.uncensoredSum + genY.sum + ( this.pastData.defined ? this.pastData.sum : 0)) / (sigma * sigma) + informedMeanPrior.b;
        muCondMean = tmpB / tmpA;
        mu = zygotine.S.normal.sample(1, muCondMean, 1 / Math.sqrt(tmpA));

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
        result.addLogOelToMuChains(this.oel);
    }

    result.totalNumberOfIterations = iter;
    return result;
};

/*******************************************/