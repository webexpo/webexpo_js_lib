/*eslint 
    valid-jsdoc: 0,
    no-extra-parens: 0
*/
/// <reference path="A.js" />
/// <reference path="M0.js" />
/// <reference path="M0.js" />
/// <reference path="M2modelIntro.js" />
/// <reference path="M2modelPerSe.js" />
/// <reference path="NUM.js" />
/// <reference path="S.js" />
/// <reference path="U.js" />

var zO = zygotine.O;
/**
@function zygotine.O.rUnifLogP1
@param {number[]} logPLim
@param {boolean} lowerTail 
@param {number} usample - uniform sample in (0,1) not required
@returns {number}
*/
zO.rUnifLogP1 = function (logPLim, lowerTail, usample) {
    // Le premier paramètre est requis. Si lowerTail est undefined alors il sera mis à true.
    // Nous avons maintenu le paramètre usample qui dans la pratique, autre que le débogage, ne sera pas utilisé.

    //La version R donne un paramètre "size" qui n'est pas utilisé
    //Ici nous faisons l'hypothèse que logPLim est de longueur 2 ce que laisse entendre le commentaire de PB.

    if (typeof usample === "undefined") {
        usample = zS.uniform.sample(1, 0, 1);
    }

    if (typeof lowerTail === "undefined") {
        lowerTail = true;
    }

    //# To sample 'size' values uniformly in c(exp(logp.lim[1]), exp(logp.lim[2]))
    var w = logPLim[1] - logPLim[0];
    if (!lowerTail) {
        w = Math.abs(w);
    }

    var logP = Math.log(usample) + zygotine.S.exponential.cdf(w, 1, true, true); // var logP = Math.log(u) + Math.log(jStat.exponential.cdf(w, 1.0));
    var x = zygotine.S.exponential.icdf(logP, 1, true, true);     //var x = jStat.exponential.inv(Math.exp(logP), 1); // ExponentialDistribution.QExp(logP, logP: true);
    x = Math.max(logPLim[0], logPLim[1]) - x;
    return x;
};

zygotine.O.rGammaTruncated = function (shape, rate, range) {
    // shape = alpha et rate = beta
    var pgamma = function (x) { return zS.gamma.cdf.usingRateParameter(x, shape, rate, true, true); };
    var logPLim = [pgamma(range[0]), pgamma(range[1])];
    var logP = zO.rUnifLogP1(logPLim);
    return zS.gamma.icdf.usingRateParameter(logP, shape, rate, true, true);
}; // rGammaTruncated

/**
@function zygotine.O.randomPow 
@param {number} a
@param {number[]} range - un tableau de longueur 2 décrivant un interval.
@returns {number} 
*/
zygotine.O.randomPow = function (a, range) {
    // sample a value from f(x) = 1/x^a on the range specified

    var u = zS.uniform.sample(1);
    var z;
    if (a === 1) {
        let z0 = (1 - u) * Math.log(range[0]);
        let z1 = u * Math.log(isFinite(range[1]) ? range[1] : 1.0e8);
        z = Math.exp(z1 - z0);
    } else {
        let tmp = 1 - a;
        let z0 = Math.pow(range[0], tmp);
        let z1 = Math.pow(range[1], tmp);
        z = Math.pow(((1 - u) * z0) + (u * z1), 1.0 / tmp);
    }

    return z;
}; //random.pow 

/**
@function zygotine.O.rNormRightCensored 
@param {number} mu
@param {number} sd
@param {number} upperBound
@param {number} usample - un nombre aléatoire tiré de façon uniforme de l'intervalle (0,1]
@returns {number}
*/
zygotine.O.rNormRightCensored = function (mu, sd, upperBound, usample) {

    if (typeof usample === 'undefined') {
        usample = zS.uniform.sample(1, 0, 1);  //jStat.uniform.sample(0, 1);
    }

    var logP = Math.log(usample) + zS.normal.cdf(upperBound, mu, sd, true, true); // Math.log(UniformDistribution.RUnif()) + NormalDistribution.PNorm(upper, mu: mu, sigma: sd, log_p: true);
    var y = zS.normal.icdf(logP, mu, sd, true, true);
    return y;
}; //# end of rnorm.right.censored

/**
@function zygotine.O.rNormLeftCensored 
@param {number} mu
@param {number} sd
@param {number} lowerBound
@returns {number}
*/
zygotine.O.rNormLeftCensored = function (mu, sd, lowerBound, usample) {
    // tous les paramètres sont des scalaires 
    var y = -zO.rNormRightCensored(-mu, sd, -lowerBound, usample);
    return y;
}; //# end of rnorm.left.censored


/**
@function zygotine.O.rNormIntervalCensored 
@param {number} mu
@param {number} sd
@param {number} lowerBound
@param {number} upperBound
@returns {number}
*/
zygotine.O.rNormIntervalCensored = function (mu, sd, lowerBound, upperBound, usample) {
    // tous les paramètres sont des scalaires 
    if (typeof usample === 'undefined') {
        usample = zS.uniform.sample(1, 0, 1); //jStat.uniform.sample(0, 1);
    }

    var logPLower = zS.normal.cdf(lowerBound, mu, sd, true, true);   // NormalDistribution.PNorm(lower, mu: mu, sigma: sd, log_p: true);
    var logPUpper = zS.normal.cdf(upperBound, mu, sd, true, true);
    var logP = zO.rUnifLogP1([logPLower, logPUpper], true, usample);
    var y = zS.normal.icdf(logP, mu, sd, true, true); // NormalDistribution.QNorm(logP, mu: mu, sigma: sd, log_p: true);
    return y;
};//# end of rnorm.interval.censored


/**
@function zygotine.O.rNormCensored
@param {number} mu
@param {number} sd
@param {number[]} lower
@param {number[]} upper
@param {boolean} negativeDisallowed
@returns {number[]}
*/
zygotine.O.rNormCensored = function (mu, sd, lower, upper, negativeDisallowed) {
    // tous les paramètres doivent être définis a l'exception de negativeDisallowed qui s'il ne l'est pas sera false
    // mu, sd sont des scalaires alors que lower et upper sont des tableaux.
    // Concernant les tableaux lower et upper, un seul peut etre vide. Si aucun des tableaux n'est vide alors ils sont de même longueur.
    // Si les longueurs de lower et de upper sont toutes 2 > 0, alors elles doivent être égales.



    //nouveau 24 avril
    //function rNormIntervalCensored(mu, sd, lower, upper) {
    //    // tous les paramètres sont des scalaires 

    //    var logPLower =   zS.normal.cdf(lower, mu, sd, true, true);   // NormalDistribution.PNorm(lower, mu: mu, sigma: sd, log_p: true);
    //    var logPUpper = zS.normal.cdf(upper, mu, sd, true, true);
    //    var logP = rUnifLogP1([logPLower, logPUpper], true, usample);
    //    var y =  zS.normal.icdf( logP, mu, sd, true, true); // NormalDistribution.QNorm(logP, mu: mu, sigma: sd, log_p: true);
    //    return y;
    //}//# end of rnorm.interval.censored

    //nouveau le 24 avril
    //function rNormRightCensored(mu, sd, upper) {
    //    // tous les paramètres sont des scalaires 
    //    //# This function was incorrectly named rnorm.left.censored in previous versions
    //    //# mu, sd: same length (1, or same as length(upper))

    //    /*
    //        log.p <- log(runif(length(upper))) + pnorm(upper, mean=mu, sd=sd, log.p=T)
    //        y <- qnorm(log.p, mean=mu, sd=sd, log.p=T)
    //    */
    //    var logP = Math.log(usample) + zS.normal.cdf(upper, mu, sd, true, true); // Math.log(UniformDistribution.RUnif()) + NormalDistribution.PNorm(upper, mu: mu, sigma: sd, log_p: true);
    //    var y = zS.normal.icdf(logP, mu, sd, true, true); 
    //    return y;
    //} //# end of rnorm.right.censored

    //nouveau 24 avril
    //function rNormLeftCensored(mu, sd, lower) {
    //    // tous les paramètres sont des scalaires 
    //    var y = -rNormRightCensored(-mu, sd, -lower);
    //    return y;
    //} //# end of rnorm.left.censored

    var maxLength = Math.max(lower.length, upper.length);
    var lowerBound = lower;
    var upperBound = upper;

    var /** {boolean} */ leftCensored = lowerBound.length > 0;
    var /** {boolean} */ rightCensored = upperBound.length > 0;
    var /** {number[]} */ z = new Array(maxLength);
    //les fonctions
    var
        rNormIntervalCensored = zO.rNormIntervalCensored,
        rNormLeftCensored = zO.rNormLeftCensored,
        rNormRightCensored = zO.rNormRightCensored;


    var i;
    if (leftCensored) {
        if (rightCensored) {
            for (i = 0; i < maxLength; i++) {
                //# it is interval-censored
                z[i] = rNormIntervalCensored(mu, sd, lowerBound[i], upperBound[i]);
            }
        }
        else {
            for (i = 0; i < maxLength; i++) {
                //# is is left-censored
                z[i] = rNormLeftCensored(mu, sd, lowerBound[i]);
            }
        }
    }
    else {
        for (i = 0; i < maxLength; i++) {
            //# then it is right-censored
            z[i] = rNormRightCensored(mu, sd, upperBound[i]);
        }
    }

    return z;
};

/**
@constructor
@param {boolean} logNormalDistrn
@property {number[]} gt
@property {number[]} lt
@property {number[]} i
@property {number} sum
@property {number} sum2
@property {boolean} logNormalDistrn
*/
zygotine.O.YGen = function (logNormalDistrn) {
    this.gt = [];
    this.lt = [];
    this.i = [];
    this.sum = 0.0;
    this.sum2 = 0.0;
    this.logNormalDistrn = logNormalDistrn;
};


/**
* @function 
* @param {zygotine.M.DataSummary} data
* @param {number} mu
* @param {number} sigma
* @param {boolean} logNormalDistrn
*
* @returns {zygotine.O.YGen}
*/
zygotine.O.YGen.getInstance =
    function (data, mu, sigma, logNormalDistrn) {
        return zygotine.O.YGen.inits(data, mu, sigma, logNormalDistrn);
    };


/**
    @function
    @param {zygotine.M.DataSummary} data
    @param {number} mu
    @param {number} sigma
    @param {boolean} logNormalDistrn
    @returns {zygotine.O.YGen}
*/
zygotine.O.YGen.inits = function (data, mu, sigma, logNormalDistrn) {
    var newYGen = new zygotine.O.YGen(logNormalDistrn);

    /**
    * @type {zygotine.m.RawMeasures} subset
    */
    var rNormCensored = zygotine.O.rNormCensored;
    newYGen.sum = 0;
    if (data.gt.length > 0) {
        newYGen.gt = rNormCensored(mu, sigma, data.gt, [], false);
        newYGen.sum += zU.sum(newYGen.gt);
        newYGen.sum2 += zU.sumSqr(newYGen.gt);
    }

    if (data.lt.length > 0) {
        newYGen.lt = rNormCensored(mu, sigma, [], data.lt, false);
        newYGen.sum += zU.sum(newYGen.lt);
        newYGen.sum2 += zU.sumSqr(newYGen.lt);
    }

    if (data.i.lower.length > 0) {
        newYGen.i = rNormCensored(mu, sigma, data.i.lower, data.i.upper, false);
        newYGen.sum += zU.sum(newYGen.i);
        newYGen.sum2 += zU.sumSqr(newYGen.i);
    }

    return newYGen;
};


/**
* @constructor
* @param {number} n - int
* @param {number} beta - double
* @param {number} lNormMu - double
* @param {number} lNormSd - double
*
* @property {number} n - int
* @property {number} b  - double
* @property {number} k - double - sert à la définition de f
* @property {number} m  - double
* @property {number} s  - double
* @property {number} s2  - double

*/
zygotine.O.SigmaGenObject = function (n, beta, lNormMu, lNormSd) {
    this.n = n;
    this.b = beta;
    this.m = lNormMu;
    this.s = lNormSd;
    this.s2 = lNormSd * lNormSd;
    this.k = NaN;
};

zygotine.O.SigmaGenObject.prototype = {

    f: function (sigma) {
        // f <- function(sigma, n=N, b=beta, m=lnorm.mu, s=lnorm.sd, k=log.k){exp(log.f(sigma, n, b, m, s)+k)}
        var tmp = Math.exp(this.logF(sigma) + this.k);
        return tmp;
    },

    /** @function
    *
    *@returns {zygotine.Num.NumericIntegration.Results}
    */
    fCum: function (lowerBound, upperBound) {
        if (isNaN(this.k)) {
            return NaN;
        }

        var self = this;
        var F = function (xA) {
            return xA.map(self.f, self);
        };

        return new zNum.NumericIntegration(F, lowerBound, upperBound).compute();
    },

    logF: function (sigma) {
        var A = this;
        return -(A.n + 1.0) * Math.log(sigma) - A.b / Math.pow(sigma, 2) - Math.pow(Math.log(sigma) - A.m, 2) / (2.0 * A.s2);
    },

    logFPrime: function (sigma) {
        var A = this;
        return (-(A.n + 1.0) / sigma) + (2.0 * A.b / Math.pow(sigma, 3)) - ((Math.log(sigma) - A.m) / (sigma * A.s2));
    }, //LogFPrime

    logFSecond: function (sigma) {
        var A = this;
        var sigma2 = Math.pow(sigma, 2);
        var x = (A.n + 1.0) / sigma2 - (6.0 * A.b / Math.pow(sigma, 4)) + ((Math.log(sigma) - A.m - 1.0) / (A.s2 * sigma2));
        return x;
    }, //LogFSecond

    LogFThird: function (sigma) {
        var A = this;
        return -2.0 * (A.n + 1) / Math.pow(sigma, 3.0) + 24 * A.b / Math.pow(sigma, 5) + (3 - 2 * Math.log(sigma) + 2 * A.m) / (A.s2 * Math.pow(sigma, 3));
    },

    show: function () {
        var fmt = zU.fmt;
        return fmt("n={0}, b={1}, m={2}, s={3}, s2={4}, k={5}", this.n, this.b, this.m, this.s, this.s2, this.k);
    }
};


/**
* @constructor 
* @param {number} val - double
* @param {number} change - double
*
* @property {number} value - double
* @property {number} change - double
*/
zygotine.O.HistoryPoint = function (val, change) {
    this.value = val;
    this.change = change;
};


/**
* @constructor
* @param {zygotine.O.SigmaGenObject} funs
* @param {number} xStart
* @param {object} xRange
*
* @property {number} epsilon
* @property {number} fRatioRemote
* @property {number} mode
* @property {number} inflexionPoint
* @property {number} remoteLeft
* @property {number} remoteRight
*/
zygotine.O.ReferencePoints = function (funs, xStart, xRange) {

    // Find x with log.f.prime > 0
    var x = xStart;
    var cntn = funs.logFPrime(x) <= 0;
    while (cntn) {
        x = x / 2.0;
        cntn = funs.logFPrime(x) <= 0;
    }

    var xLeft = x;
    var change;
    var previousX;

    // Find mode
    cntn = true;
    while (cntn) {
        change = funs.logFPrime(x) / funs.logFSecond(x);
        x = x - change;
        cntn = Math.abs(change) > this.epsilon;
    }

    this.mode = x;
    var target = funs.logF(this.mode) + Math.log(this.fRatioRemote);

    // Find inflexion point by Newton-Raphson
    var hstory = [];
    cntn = true;
    while (cntn) {
        hstory.push(x);
        previousX = x;
        change = funs.logFSecond(x) / funs.LogFThird(x);
        hstory.push(change);
        x = x - change;
        if (x < 0) {
            x = previousX / 2.0;
        }

        cntn = Math.abs(change) > this.epsilon;
    }

    hstory.push(x);
    this.inflexionPoint = x;

    // Find right-side remote point
    cntn = true;
    var maxXRange = zU.max(xRange);
    while (cntn) {
        change = (funs.logF(x) - target) / funs.logFPrime(x);
        x = x - change;
        cntn = (Math.abs(change) > this.epsilon) && (x < maxXRange);
    }

    this.remoteRight = x;

    // Find left-side remote point

    x = xLeft;
    cntn = true;
    while (cntn) {
        previousX = x;
        change = (funs.logF(x) - target) / funs.logFPrime(x);
        x = x - change;
        if (x < 0) {
            x = previousX / 2.0;
        }

        change = x - previousX;
        cntn = Math.abs(change) > this.epsilon;
    }

    this.remoteLeft = x;
}; //ref.points


zygotine.O.ReferencePoints.prototype = {
    fRatioRemote: 1e-8,
    epsilon: 1e-6
};


zygotine.O.sigmaGenICdf = function (n, beta, lNormMu, lNormSd, range, newtonRaphson, interval, returnProb, banerjee, u) {
    // a l'appel il doit y avoir un argument pour chacun des 4 premiers paramètres paramètre.

    if (typeof u === "undefined") {
        u = jStat.uniform.sample(0, 1);
    }

    if (typeof banerjee === "undefined") {
        banerjee = false;
    }

    if (typeof returnProb === "undefined") {
        returnProb = false;
    }

    if (typeof interval === "undefined") {
        interval = [];
    }

    if (typeof newtonRaphson === "undefined") {
        newtonRaphson = { epsilon: 1e-8, maxNiter: 100 };
    }

    if (typeof range === "undefined") {
        range = [0, Infinity];
    }

    var funs = new zO.SigmaGenObject(n, beta, lNormMu, lNormSd);

    //if (!Banerjee.null.beta)
    //{
    /** @type {number} */
    var xStart = Math.exp(lNormMu);
    /** @type {zygotine.O.ReferencePoints} */
    var z = new zO.ReferencePoints(funs, xStart, range);
    /** @type {zygotine.Num.NumericIntegration.Result} */
    var areaResult;

    var area; // double
    //}

    if (beta === 0) {
        var fMode = Math.exp(funs.logF(z.mode));
        var rightLim = [z.mode, range[1]];
        area = (z.mode - range[0]) * fMode - (Math.pow(rightLim[1], -n) - Math.pow(rightLim[0], -n)) / n;
        var logK = Math.log(area);
    } else {
        let alpha = (n - 2.0) / 2.0;
        logK = alpha * Math.log(beta) - zNum.logGamma(alpha) - Math.log(2 * Math.PI) / 2 - Math.log(lNormSd);
    }

    // Define f and F (fcum) (somewhat) standardized
    funs.k = logK;

    /*
    La définition de fCum est intégré à zygotine.O.SigmaGenObject (ici l'objet funs).
    L'affectation d'une valeur au paramètre k permet de rendre opérationnel funs.f et funs.fCum.  
    */

    if (range[0] < z.remoteLeft) {
        range[0] = z.remoteLeft;
    }

    if (range[1] > z.remoteRight) {
        range[1] = z.remoteRight;
    }

    areaResult = funs.fCum(range[0], range[1]); // areaResult est un objet, la valeur de l'intégrale est donnée par areaResult.result
    area = (areaResult.ier === 0) ? areaResult.result : NaN;

    // Sample a value !!! 

    var target = area * u;

    // Run Newton-Raphson algorithm to find x such that fcum(x) = target
    /** @type {number} */
    var x = z.mode; // Start at inflexion point of fcum (= mode of f)
    /** @type {boolean} */
    var cntn = true;
    /** @type {number} */
    var iter = 0;
    /** @type {number} */
    var change = NaN;
    /** @type {boolean} */
    var bisectionalSearch = false;

    // we will save values of x & change along the march to solution, in case Newton-Raphson algorithm
    // does not converge: from these we will have a lower & upper bound on which to run a bisectional search

    /** @type {zygotine.O.HistoryPoint[]} */
    var changeHistory = []; // (start with empty vectors, which will be augmented at each passage
    /** @type {zygotine.O.HistoryPoint} */
    var point;
    /** @type {boolean} */
    var converged;
    /** @type {number} */
    var lastX;
    while (cntn) {
        iter = iter + 1;
        change = (funs.fCum(range[0], x).result - target) / funs.f(x);
        changeHistory.push(new zygotine.O.HistoryPoint(x, change));
        x = x - change;
        converged = Math.abs(change) < newtonRaphson.epsilon;
        cntn = !converged && (iter < newtonRaphson.maxNiter);

        // If x steps out of (z$remote.left, z$remote.right),
        // then stop Newton-Raphson algorithm and proceed to Bisectional-search (below)
        if (x < z.remoteLeft || (x > z.remoteRight)) {
            lastX = x + change;
            x = (change < 0) ? z.remoteRight : z.remoteLeft;
            change = lastX - x;
            bisectionalSearch = true;
            cntn = false;
        }
    }

    if (!converged) {
        if (!bisectionalSearch) {
            // Find x lower & upper bound from visited (history) values
            /** @type {zygotine.O.HistoryPoint[]} */
            var w = changeHistory.filter(function (point) { return point.change < 0; });
            var xLower = w.reduce(function (acc, point) { return Math.max(acc, point.value); }, -Infinity);
            var xUpper = w.reduce(function (acc, point) { return Math.min(acc, point.value); }, Infinity);
            if (isFinite(xLower) && isFinite(xUpper)) {
                x = xUpper;                // scalar
                change = xLower - xUpper; // scalar
                bisectionalSearch = true;
            } else if (isFinite(z.remoteLeft) && isFinite(z.remoteRight)) {
                x = z.remoteLeft;
                change = z.remoteRight - z.remoteLeft;
                bisectionalSearch = true;
            }
        }

        if (bisectionalSearch) {
            x = x + change / 2.0;
            change = Math.abs(change / 2.0);
            var dir;
            while (!converged) {
                dir = funs.fCum(range[0], x).result > target ? -1.0 : 1.0; // scalar
                change = change / 2.0;
                x = x + dir * change;
                converged = change < newtonRaphson.epsilon;
            }
        }
        else {
            x = NaN;
        }
    }

    return x;
};

/**
* @function
* @param {number} n - int
* @param {number} beta - double
* @param {number} lNormMu - double
* @param {number} lNormSd - double
* @param {number} u - prn tiré de l'interval (0,1] à des fins de débogage 
* 
* @returns {number} 
*/
zygotine.O.sigmaGenICdf4InformedVar = function (n, beta, lNormMu, lNormSd, u) {
    if (typeof u === "undefined") {
        u = zS.uniform.sample(1, 0, 1);
    }

    var banerjee = false;
    var returnProb = false;
    var interval = [];
    var newtonRaphson = { epsilon: 1e-8, maxNiter: 100 };
    var range = [0, Infinity];
    return zygotine.O.sigmaGenICdf(n, beta, lNormMu, lNormSd, range, newtonRaphson, interval, returnProb, banerjee, u);
};

zygotine.O.sigmaGenICdf4BW = function (n, beta, lNormMu, lNormSd, range, u) {
    if (typeof u === "undefined") {
        u = zS.uniform.sample(1, 0, 1);
    }
    var banerjee = false;
    var returnProb = false;
    var interval = [];
    var newtonRaphson = { epsilon: 1e-8, maxNiter: 100 };
    return zygotine.O.sigmaGenICdf(n, beta, lNormMu, lNormSd, range, newtonRaphson, interval, returnProb, banerjee, u);
};
