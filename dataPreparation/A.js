zygotine = {};
zygotine.B = {};
zygotine.O = {};

zygotine.udef = "undefined";


zygotine.isUndefined = function (a) {
    return typeof a === zygotine.udef;
};

zygotine.isnan = function (a) {
    if (typeof a === 'number') {
        return isNaN(a);
    }

    return true;
};

zygotine.isABoolean = function (a) {
    return typeof a === 'boolean';
};


zygotine.isANumber = function (a) {
    return typeof a === 'number';
};

/*

zygotine.M

    Measure = function () 
        @constructor
        @property {number} a, b - double
        @property {number} ordinal - int
        @property {string} type - measure type
        @property {string} worker - the worker, if available, to which the measure is associated


    MeasureList = function (s)
        @constructor
        {string} s - une chaine représantant les mesures sans 'erreur de mesure'
        @property {boolean} anyCensored
        @property {boolean} hasError
        @property {boolean} hasME
        @property {zygotine.M.Messages} messages
        @property {object} measureByType
        @property {object} measureCountByType
        @property {zygotine.M.Measure[]} measureList
        @property {string} METext
        @property {RegExp} rgxME
    
    Message = function (msg, error, src)
        @constructor
        {string} msg
        {boolean} error
        {string} src
        @property {string} msgText, level, source

    Messages = function (errSrc) 
        @constructor
        {string} errSrc 
        @property {boolean} hasError - by default it is false
        @property {string} errorSrc; 
        @property {zygotine.M.Message[]} msgList

zygotine.O 

    rNormCensored = function (mu, sd, lower, upper, negativeDisallowed)
        {number} mu, sd
        {number[]} lower, upper
        {boolean} negativeDisallowed
        @returns {number[]}

    rUnifLogP1 = function (logPLim, lowerTail, usample)
        {number[]} logPLim
        {boolean} lowerTail 
        {number} usample
        @returns {number}

    YGen = function()
        @constructor
        @property {number[]} logGt, gt, logLt, lt, logI, i

    YGen.getInstance(data, mu, sigma, logNormalDistrn)
        {zygotine.M.DataSummary} data
        {number} mu, sigma
        {boolean} logNormalDistrn
        @returns {zygotine.O.YGen}

    YGen.inits = function (data, mu, sigma, logNormalDistrn)        /*
        {zygotine.M.DataSummary} data
        {number} mu, sigma
        {boolean} logNormalDistrn
        @returns {zygotine.O.YGen}



*/
