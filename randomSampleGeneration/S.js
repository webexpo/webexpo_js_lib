/* eslint 
    xvalid-jsdoc: 0,
    no-extra-parens : 0
*/
/// <reference path="A.js" />
/// <reference path="U.js" />

zygotine.S = {};
zS = zygotine.S;
zygotine.S.normal = {};
zygotine.S.exponential = {};

zS.prng = {
    generator: undefined,
    seed: undefined
};

zS.prng.init = function (seed) {

    if (zygotine.isnan(seed)) {
        seed = new Date().getTime();
    }

    seed = seed >>> 0; // obtenir un uint32
    zS.prng.seed = seed; // on stocke la racine
    zS.prng.generator = new zygotine.MT.MersenneTwister(seed);
};

zS.prng.init(12);

/***> loi uniforme <***/

zS.uniform = {};

zS.uniform.sample = function (n, a, b) {

    var isUndefined = zygotine.isUndefined;
    var rep = [];

    // on copie le comportement de R sauf pour n, qui lorsque non défini, prendra la valeur de 1.
    // pour la fonction runif de R le paramètre n est requis.
    if (isUndefined(b)) {
        b = 1.0;
    }

    if (isUndefined(a)) {
        a = 0.0;
    }

    if (isUndefined(n)) {
        n = 1;
    }

    if (!isFinite(a) || !isFinite(b) || b < a) {
        rep = Array.apply(null, { length: n }).map(xx => NaN);
    } else {
        /**  
        * @type {zygotine.MT.MersenneTwister}
        */
        if (a === b) {
            rep = Array.apply(null, { length: n }).map(xx => a);
        } else {
            let rand = zS.prng.generator.genrand_real3;
            let len = (b - a);
            rep = Array.apply(null, { length: n }).map(xx => a + len * zS.uniform.rand());
        }
    }

    return rep.length === 1 ? rep[0] : rep;
};

zS.uniform.rand = function () {
    return zS.prng.generator.genrand_real3();
}; //END of zS.uniform.sample

zS.uniform.randN = function (n) {
    let runif = zS.prng.generator.genrand_real3();
    return Array.apply(null, { length: n }).map(xx => runif());
};

/***> loi normale <***/

zS.normal.pdf = function (x, mu, sigma, giveLog) {

    var isUndefined = zygotine.isUndefined;

    if (isUndefined(giveLog)) {
        giveLog = false;
    }

    if (isUndefined(sigma)) {
        sigma = 1.0;
    }

    if (isUndefined(mu)) {
        mu = 0.0;
    }

    if (isNaN(x) || isNaN(mu) || isNaN(sigma)) {
        return NaN;
    }

    if (!isFinite(sigma)) {
        return giveLog ? -Infinity : 0.0;
    }

    if (!isFinite(x) && (mu === x)) {
        return Nan;/* x-mu is NaN */
    }

    if (sigma <= 0) {
        if (sigma < 0) {
            return NaN; // ML_ERR_return_NAN;
        }
        /* sigma == 0 */
        return (x === mu) ? Infinity : (giveLog ? -Infinity : 0.0);
    }

    x = (x - mu) / sigma;

    if (!isFinite(x)) {
        return giveLog ? -Infinity : 0.0;
    }

    x = Math.abs(x);
    if (x >= 2 * Math.sqrt(zU.CONST.DBL_MAX)) {
        return giveLog ? -Infinity : 0.0;
    }

    if (giveLog) {
        return -(zU.CONST.M_LN_SQRT_2PI + 0.5 * x * x + Math.log(sigma));
    }

    if (x < 5) {
        return zU.CONST.M_1_SQRT_2PI * Math.exp(-0.5 * x * x) / sigma;
    }

    /* ELSE:
    
     * x*x  may lose upto about two digits accuracy for "large" x
     * Morten Welinder's proposal for PR#15620
     * https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620
    
     * -- 1 --  No hoop jumping when we underflow to zero anyway:
    
     *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
     *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
     * but "thanks" to denormalized numbers, underflow happens a bit later,
     *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
     * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
     * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
     *              =IEEE=  38.58601
     * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
     */
    if (x > Math.sqrt(-2 * zU.CONST.M_LN2 * (zU.CONST.DBL_MIN_EXP + 1 - zU.CONST.DBL_MANT_DIG))) {
        return 0.;
    }

    /* Now, to get full accurary, split x into two parts,
     *  x = x1+x2, such that |x2| <= 2^-16.
     * Assuming that we are using IEEE doubles, that means that
     * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).
    
     * If we do not have IEEE this is still an improvement over the naive formula.
     */
    var x1 = Math.round(x * 65536) / 65536;
    var x2 = x - x1;
    return zU.CONST.M_1_SQRT_2PI / sigma * (Math.exp(-0.5 * x1 * x1) * Math.exp((-0.5 * x2 - x1) * x2));

};

/*
* @function
* @param {number} x
* @param {number} mu - mean - default to 0
* @param {number} sigma - sd - default to 1
* @param {boolean} lowerTail  - default to true
* @param {boolean} logP - default to false
* @returns {number}
*/
zygotine.S.normal.cdf = function (x, mu, sigma, lowerTail, logP) {
    var ud = zygotine.udef;
    var isABoolean = zygotine.isABoolean;
    var isANumber = zygotine.isANumber;

    // on vérifie les paramètres ayant des valeurs implicites
    if (!isABoolean(logP)) {
        logP = false;
    }

    if (!isABoolean(lowerTail)) {
        lowerTail = true;
    }

    if (!isANumber(sigma)) {
        sigma = 1.0;
    }

    if (!isANumber(mu)) {
        mu = 0.0;
    }

    var
        p, cp = 0.0;

    if (isNaN(x) || isNaN(mu) || isNaN(sigma)) {
        return NaN;
    }

    if (!isFinite(x) && (mu === x)) {
        return NaN;
    }


    if (sigma <= 0) {
        if (sigma < 0) {
            return NaN;
        }

        /* sigma = 0 : */
        return (x < mu) ? (lowerTail ? (logP ? -Infinity : 0.0) : (logP ? 0.0 : 1.0)) : (lowerTail ? (logP ? 0.0 : 1.0) : (logP ? -Infinity : 0.0));
    }

    p = (x - mu) / sigma;
    if (!isFinite(p)) {
        return (x < mu) ? (lowerTail ? (logP ? -Infinity : 0.0) : (logP ? 0.0 : 1.0)) : (lowerTail ? (logP ? 0.0 : 1.0) : (logP ? -Infinity : 0.0));
    }

    x = p;
    var o = zygotine.S._pnorm_both(x, (lowerTail ? 0 : 1), logP);
    return lowerTail ? o.cum : o.ccum;

};

/**
    private
*/
zygotine.S.normal.pnorm_bothCONST = {

    a: [
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    ],

    b: [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956
    ],

    c: [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8
    ],

    d: [
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
    ],

    p: [
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
    ],

    q: [
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5
    ]
};

/**
* privé
* @function 
* @param {number} x - double
* @param {number} iTail - int
* @param {boolean} logP
* @returns {object} un objet composé de 2 propriétés : cum et ccum
*/
zygotine.S._pnorm_both = function (x, iTail, logP) {
    /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
       if(lower) return  *cum := P[X <= x]
       if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
    */
    const SIXTEN = 16;
    var // double
        xden, xnum, temp, del, eps = zU.CONST.DBL_EPSILON, xsq, y;
    var
        min = zU.CONST.DBL_MIN;
    var // int
        i, lower, upper;
    var
        rep = { cum: NaN, ccum: NaN };

    if (isNaN(x)) {
        rep.cum = x;
        rep.ccum = x;
        return rep;
    }

    var
        tmp = zygotine.S.normal.pnorm_bothCONST;

    var
        a = tmp.a,
        b = tmp.b,
        c = tmp.c,
        d = tmp.d,
        p = tmp.p,
        q = tmp.q;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = iTail !== 1;
    upper = iTail !== 0;
    y = Math.abs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
        if (y > eps) {
            xsq = x * x;
            xnum = a[4] * xsq;
            xden = xsq;
            for (i = 0; i < 3; ++i) {
                xnum = (xnum + a[i]) * xsq;
                xden = (xden + b[i]) * xsq;
            }
        }
        else {
            xnum = xden = 0.0;
        }

        temp = x * (xnum + a[3]) / (xden + b[3]);
        if (lower) {
            rep.cum = 0.5 + temp;
        }

        if (upper) {
            rep.ccum = 0.5 - temp;
        }

        if (logP) {
            if (lower) {
                rep.cum = Math.log(rep.cum);
            }

            if (upper) {
                rep.ccum = Math.log(rep.ccum);
            }
        }
    }
    else if (y <= zU.CONST.M_SQRT_32) {
        /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
        xnum = c[8] * y;
        xden = y;
        for (i = 0; i < 7; ++i) {
            xnum = (xnum + c[i]) * y;
            xden = (xden + d[i]) * y;
        }

        temp = (xnum + c[7]) / (xden + d[7]);
        xsq = Math.trunc(y * SIXTEN) / SIXTEN;
        // expansion de do_del(y)
        del = (y - xsq) * (y + xsq);
        if (logP) {
            rep.cum = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
            if ((lower && x > 0.0) || (upper && x <= 0.0))
                rep.ccum = Math.log1p(-Math.exp(-xsq * xsq * 0.5) *
                      Math.exp(-del * 0.5) * temp);
        }
        else {
            rep.cum = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
            rep.ccum = 1.0 - rep.cum;
        }

        if (x > 0.0) {/* swap  ccum <--> cum */
            temp = rep.cum;
            if (lower) {
                rep.cum = rep.ccum;
            }

            rep.ccum = temp;
        }
    }
    else if ((logP && y < 1e170) /* avoid underflow below */
        /*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
         * Then, make use of  Abramowitz & Stegun, 26.2.13, something like
    
         xsq = x*x;
    
         if(xsq * DBL_EPSILON < 1.)
            del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
         else
            del = 0.;
         *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
         *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
    
         swap_tail;
    
         [Yes, but xsq might be infinite.]
    
        */
        || (lower && -37.5193 < x && x < 8.2924)
        || (upper && -8.2924 < x && x < 37.5193)) {

        /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
        xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
        xnum = p[5] * xsq;
        xden = xsq;
        for (i = 0; i < 4; ++i) {
            xnum = (xnum + p[i]) * xsq;
            xden = (xden + q[i]) * xsq;
        }

        temp = xsq * (xnum + p[4]) / (xden + q[4]);
        temp = (zU.CONST.M_1_SQRT_2PI - temp) / y;

        // do_del(x);
        xsq = Math.trunc(x * SIXTEN) / SIXTEN;
        del = (x - xsq) * (x + xsq);
        if (logP) {
            rep.cum = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
            if ((lower && x > 0.0) || (upper && x <= 0.0))
                rep.ccum = Math.log1p(-Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp);
        }
        else {
            rep.cum = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
            rep.ccum = 1.0 - rep.cum;
        }

        if (x > 0.0) {/* swap  ccum <--> cum */
            temp = rep.cum;
            if (lower) {
                rep.cum = rep.ccum;
            }

            rep.ccum = temp;
        }
    }
    else { /* large x such that probs are 0 or 1 */
        if (x > 0) {
            rep.cum = (logP ? 0.0 : 1.0);
            rep.ccum = (logP ? -Infinity : 0.0);
        }
        else {
            rep.cum = (logP ? -Infinity : 0.0);
            rep.ccum = (logP ? 0.0 : 1.0);
        }
    }

    if (logP) {
        if (rep.cum > -min) {
            rep.cum = -0.0;
        }

        if (rep.ccum > -min) {
            rep.ccum = -0.0;
        }
    }
    else {
        if (rep.cum < min) {
            rep.cum = 0.0;
        }

        if (rep.ccum < min) {
            rep.ccum = 0.0;
        }
    }

    return rep;

};

/*
* @function
* @param {number} p
* @param {number} mu - mean - default to 0
* @param {number} sigma - sd - default to 1
* @param {boolean} lowerTail  - default to true
* @param {boolean} logP - default to false
* @returns {number}
*/
zygotine.S.normal.icdf = function (p, mu, sigma, lowerTail, logP) {

    var ud = zygotine.udef;
    if (typeof logP === ud) {
        logP = false;
    }

    if (typeof lowerTail === ud) {
        lowerTail = true;
    }

    if (typeof sigma === ud) {
        sigma = 1.0;
    }

    if (typeof mu === ud) {
        mu = 0.0;
    }

    var
        p_, q, r, val;
    var
        R_DT_CIv = zS.R_DT_CIv,
        R_DT_qIv = zS.R_DT_qIv;

    if (isNaN(p) || isNaN(mu) || isNaN(sigma)) {
        return NaN;
    }

    //R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

    if (logP) {
        if (p > 0) {
            return NaN;
        }

        if (p === 0) { /* upper bound*/
            return lowerTail ? Infinity : -Infinity;
        }

        if (p === -Infinity) {
            return lowerTail ? -Infinity : Infinity;
        }
    }
    else { /* !log_p */
        if (p < 0 || p > 1) {
            return NaN; //ML_ERR_return_NAN
        }

        if (p === 0) {
            return lowerTail ? -Infinity : Infinity;
        }

        if (p === 1) {
            return lowerTail ? Infinity : -Infinity;
        }
    }


    if (sigma < 0) {
        return NaN;
    } //ML_ERR_return_NAN;

    if (sigma === 0) {
        return mu;
    }

    p_ = R_DT_qIv(p, lowerTail, logP);/* real lower_tail prob. p */
    q = p_ - 0.5;

    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
    
            Produces the normal deviate Z corresponding to a given lower
            tail area of P; Z is accurate to about 1 part in 10**16.
    
            (original fortran code used PARAMETER(..) for the coefficients
             and provided hash codes for checking them...)
    */
    if (Math.abs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
                q * (((((((r * 2509.0809287301226727 +
                           33430.575583588128105) * r + 67265.770927008700853) * r +
                         45921.953931549871457) * r + 13731.693765509461125) * r +
                       1971.5909503065514427) * r + 133.14166789178437745) * r +
                     3.387132872796366608)
                / (((((((r * 5226.495278852854561 +
                         28729.085735721942674) * r + 39307.89580009271061) * r +
                       21213.794301586595867) * r + 5394.1960214247511077) * r +
                     687.1870074920579083) * r + 42.313330701600911252) * r + 1.0);
    }
    else { /* closer than 0.075 from {0,1} boundary */

        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = R_DT_CIv(p, lowerTail, logP);/* 1-p */
        else
            r = p_;/* = R_DT_Iv(p) ^=  p */

        r = Math.sqrt(-((logP && ((lowerTail && q <= 0) || (!lowerTail && q > 0))) ? p : Math.log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
        if (r <= 5.0) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.0);
        }
        else { /* very close to  0 or 1 */
            r += -5.0;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.0);
        }

        if (q < 0.0)
            val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }

    return mu + sigma * val;
};

zygotine.S.normal.rand = function () {
    // In C, case INVERSION: ...

    var u1;
    const BIG = 134217728; /* 2^27 */
    /* unif_rand() alone is not of high enough precision */
    u1 = zS.uniform.rand();
    u1 = Math.trunc(BIG * u1) + zS.uniform.rand();
    return zS.normal.icdf(u1 / BIG, 0.0, 1.0, true, false);
};

zygotine.S.normal.sample = function (n, mu, sigma) {

    var
        isUndefined = zygotine.isUndefined,
        rep = [];

    if (isUndefined(sigma)) {
        sigma = 1.0;
    }

    if (isUndefined(mu)) {
        mu = 0.0;
    }

    if (isUndefined(n)) {
        n = 1;
    }

    var retValue = false, value;
    if (isNaN(mu) || !isFinite(sigma) || (sigma < 0.0)) {
        value = NaN; // ML_ERR_return_NAN;
        retValue = true;
    }

    if (sigma === 0.0 || !isFinite(mu)) {
        value = mu; /* includes mu = +/- Inf with finite sigma */
        retValue = true;
    }

    if (retValue) {
        rep = Array.apply(null, { length: n }).map(xx => value);
    } else {
        let rand = zS.normal.rand;
        rep = Array.apply(null, { length: n }).map(xx => mu + sigma * rand());
    }

    return (n === 1) ? rep[0] : rep;
};

/***> loi exponentielle <***/

/*
        rate (intensité) === (	λ > 0 ) ~~~  1 / β où β est le paramètre scale (échelle)
donc    scale ==== β > 0 et < Inf

Dans les sources de R écrits en C, les fonctions de la loi Exponentielle utilisent rate plutôt que scale alors que
dexp, qexp, pexp et rexp dans R utilisent rate.

Ici toutes les fonctions utilisent scale, comme le fait le code C de R.

L'usage qu'on fait de la loi exponentielle est limité: tous les appels faits à l'une ou l'autre des fonctions de zS.exponential utilisent un rate (intensité) fixé à 1, ce qui implique un scale (échelle) de 1.

Les tests faits des fonctions de la loi ne vise qu'à couvrir les cas autres que ceux où rate === scale === 1.0.
*/

/*
    R file : dnorm.c
    function name dnorm4
*/

/**
 *@function
* @param {number} x
* @param {number} scale - defaults to 1
* @param {boolean} lowerTail  - default to true
* @param {boolean} logP - default to false
* @returns {number}
* 
*/
zygotine.S.exponential.cdf = function (x, scale, lowerTail, logP) {

    if (!zygotine.isABoolean(logP)) {
        logP = false;
    }

    if (!zygotine.isABoolean(lowerTail)) {
        lowerTail = true;
    }

    if (!zygotine.isANumber(scale)) {
        scale = 1.0;
    }

    //if (zygotine.isUndefined(x)) {
    //    throw new Error("In zygotine.S.exponential.cdf, argument x is required!")
    //}

    if (isNaN(x) || isNaN(scale)) {
        return NaN;
    }

    if (scale < 0) {
        return NaN;
    }

    if (x <= 0.0) {
        return (lowerTail ? (logP ? -Infinity : 0.0) : (logP ? 0.0 : 1.0));
    }

    x = -(x / scale);
    if (lowerTail) {
        return logP ? (x > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(x)) : Math.log1p(-Math.exp(x))) : -Math.expm1(x);
    }

    return logP ? (x) : Math.exp(x);
}; // END OF zygotine.S.exponential.cdf

zygotine.S.exponential.cdf.usingScaleParameter = zygotine.S.exponential.cdf;

zygotine.S.exponential.cdf.usingRateParameter = function (x, rate, lowerTail, logP) {

    if (zygotine.isUndefined(rate)) {
        rate = 1.0;
    }

    if (!zygotine.isANumber(rate)) {
        throw new Error("In exponential.cdf.usingRateParameter, supplied parameter rate must be a number!");
    }

    return zygotine.S.exponential.cdf.usingScaleParameter(x, 1 / rate, lowerTail, logP);
};

/**
* R file: qexp.c
*
* @function
* @param {number} p
* @param {number} scale - defaults to 1
* @param {boolean} lowerTail  - default to true
* @param {boolean} logP - default to false
* @returns {number}
* 
*/
zygotine.S.exponential.icdf = function (p, scale, lowerTail, logP) {

    if (!zygotine.isABoolean(logP)) {
        logP = false;
    }

    if (!zygotine.isABoolean(lowerTail)) {
        lowerTail = true;
    }

    if (!zygotine.isANumber(scale)) {
        scale = 1.0;
    }

    if (!zygotine.isANumber(p) || !zygotine.isANumber(scale)) {
        return NaN;
    }

    if (scale < 0) {
        return NaN;
    }

    if ((logP && p > 0) || (!logP && (p < 0 || p > 1))) {
        return NaN;
    }

    if (p === (lowerTail ? (logP ? -Infinity : 0.0) : (logP ? 0.0 : 1.0))) {
        return 0;
    }

    return -scale * (lowerTail ? (logP ? ((p) > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(p)) : Math.log1p(-Math.exp(p))) : Math.log1p(-p)) : (logP ? (p) : Math.log(p)));
}; // END OF zygotine.S.exponential.icdf

zygotine.S.exponential.icdf.usingScaleParameter = zygotine.S.exponential.icdf;

zygotine.S.exponential.icdf.usingRateParameter = function (p, rate, lowerTail, logP) {
    if (zygotine.isUndefined(rate)) {
        rate = 1.0;
    }

    if (!zygotine.isANumber(rate)) {
        throw new Error("In exponential.icdf.usingRateParameter, supplied parameter rate must be a number!");
    }

    return zygotine.S.exponential.icdf.usingScaleParameter(p, 1 / rate, lowerTail, logP);
};

zygotine.S.exponential.sample = function (n, scale) {

    var isUndefined = zygotine.isUndefined;
    var rep = [];

    // on copie le comportement de R sauf pour n, qui lorsque non défini, prendra la valeur de 1.
    // pour runif le paramètre n est requis.

    if (zygotine.isUndefined(scale)) {
        scale = 1.0;
    }

    if (zygotine.isUndefined(n)) {
        n = 1;
    }


    if (!isFinite(scale) || scale <= 0.0) {

        let tmp = (scale === 0.0) ? 0.0 : NaN;
        rep = Array.apply(null, { length: n }).map(xx => tmp);
    } else {
        rep = Array.apply(null, { length: n }).map(xx => scale * zygotine.S.exponential.rand());
    }

    return rep.length === 1 ? rep[0] : rep;
};

zygotine.S.exponential.rand = function () {
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 16) is determined by q[n-1] = 1.0 */
    /* within standard precision */


    var /*double*/
        a = 0.0,
        u = zS.uniform.sample();
    /* precaution if u = 0 is ever returned */
    // unif_rand based on MersenneTwister returns values in the open interval (0,1).
    // so following line is superflous

    // while (u <= 0.0 || u >= 1.0) u = unif_rand(); 

    var q = zygotine.S.exponential.rand.q;

    for (; ;) {
        u += u;
        if (u > 1.0) {
            break;
        }

        a += q[0];
    }

    u -= 1.0;
    if (u <= q[0]) {
        return a + u;
    }

    var /*int*/ i = 0;

    var /*double*/
        ustar = zS.uniform.sample(),
        umin = ustar;
    do {
        ustar = zS.uniform.sample();
        if (umin > ustar) {
            umin = ustar;
        }

        i++;
    } while (u > q[i]);

    return a + umin * q[0];
};

zygotine.S.exponential.rand.q = [
    0.6931471805599453,
    0.9333736875190459,
    0.9888777961838675,
    0.9984959252914960,
    0.9998292811061389,
    0.9999833164100727,
    0.9999985691438767,
    0.9999998906925558,
    0.9999999924734159,
    0.9999999995283275,
    0.9999999999728814,
    0.9999999999985598,
    0.9999999999999289,
    0.9999999999999968,
    0.9999999999999999,
    1.0000000000000000
]; // END of zygotine.S.exponential.rand.q

/***> loi de Poisson <***/

zS.poisson = {
    cutoff: zU.CONST.M_LN2 * zU.CONST.DBL_MAX_EXP / zU.CONST.DBL_EPSILON,
    pdf: function () {
        throw new Error("Not implemented yet.");
    }
};

zS.poisson.cutoff = zU.CONST.M_LN2 * zU.CONST.DBL_MAX_EXP / zU.CONST.DBL_EPSILON;

zS.poisson.pdf = function () { };

zS.poisson.pdf.dpois_raw = function (x, lambda, logP) {

    /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
        lambda >= 0
    */
    if (lambda === 0) {
        return ((x === 0) ? (logP ? 0.0 : 1.0) : (logP ? -Infinity : 0.0));
    }

    if (!isFinite(lambda)) {
        return (logP ? -Infinity : 0.0);
    }

    if (x < 0) {
        return ((logP ? -Infinity : 0.0));
    }
    if (x <= lambda * zU.CONST.DBL_MIN) return ((logP ? (-lambda) : Math.exp(-lambda)));
    if (lambda < x * zU.CONST.DBL_MIN)
        return (logP
            ? (-lambda + x * Math.log(lambda) - zygotine.Num.logGamma(x + 1))
            : Math.exp(-lambda + x * Math.log(lambda) - zygotine.Num.logGamma(x + 1)));

    // ci-bas réécriture de : return (RVaria.R_D_fexp(RVaria.M_2PI * x, -Functions.StirlingError(x) - RVaria.bd0(x, lambda), logP));
    var vf = zU.CONST.M_2PI * x;
    var vx = -zygotine.Num.stirlingError(x) - zygotine.Num.bd0(x, lambda);
    if (logP) {
        return -0.5 * Math.log(vf) + (vx);
    } else {
        return Math.exp(vx) / Math.sqrt(vf);
    }
};

//internal static double R_D_fexp(double f, double x, bool giveLog)
//{
//    return (giveLog ? -0.5 * Math.Log(f) + (x) : Math.Exp(x) / Math.Sqrt(f));
//}

/*
    R file : pgamma.c
*/
zS.poisson.pdf.dpois_wrap = function (x_plus_1, lambda, giveLog) {
    var x;
    if (!isFinite(lambda)) {
        return (giveLog ? -Infinity : 0.0);
    }

    if (x_plus_1 > 1) {
        return zS.poisson.pdf.dpois_raw(x_plus_1 - 1, lambda, giveLog);
    }

    if (lambda > Math.abs(x_plus_1 - 1) * zS.poisson.cutoff) {
        x = -lambda - zygotine.Num.logGamma(x_plus_1);
        return (giveLog ? x : Math.exp(x));
    }
    else {
        var d = zS.poisson.pdf.dpois_raw(x_plus_1, lambda, giveLog);
        return giveLog
            ? d + Math.log(x_plus_1 / lambda)
            : d * (x_plus_1 / lambda);
    }
};

/* R file: pgamma.c
 * 
 * Asymptotic expansion to calculate the probability that Poisson variate * has value <= x.
 * Various assertions about this are made (without proof) at http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
zS.poisson.cdf = {}; // On a pas la fn cdf comme telle (ppois).
zS.poisson.cdf.ppois_asymp = function (x, lambda, lowerTail, logP) {
    var /*double*/ elfb, elfb_term;
    var /*double*/ res12, res1_term, res1_ig, res2_term, res2_ig;
    var /*double*/ dfm, pt_, s2pt, f, np;
    var /*int*/ i;
    var a = zS.poisson.cdf.ppois_asymp.coefs_a;
    var b = zS.poisson.cdf.ppois_asymp.coefs_b;
    dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    pt_ = -zS.gamma.log1pmx(dfm / x);
    s2pt = Math.sqrt(2 * x * pt_);
    if (dfm < 0) {
        s2pt = -s2pt;
    }

    res12 = 0;
    res1_ig = res1_term = Math.sqrt(x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++) {
        res12 += res1_ig * a[i];
        res12 += res2_ig * b[i];
        res1_term *= pt_ / i;
        res2_term *= 2 * pt_ / (2 * i + 1);
        res1_ig = res1_ig / x + res1_term;
        res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++) {
        elfb += elfb_term * b[i];
        elfb_term /= x;
    }

    if (!lowerTail) {
        elfb = -elfb;
    }

    f = res12 / elfb;

    np = zygotine.S.normal.cdf(s2pt, 0.0, 1.0, !lowerTail, logP);

    if (logP) {
        var n_d_over_p = zygotine.S.poisson.cdf.ppois_asymp.dpnorm(s2pt, !lowerTail, np);
        return np + Math.log1p(f * n_d_over_p);
    } else {
        var nd = zS.normal.pdf(s2pt, 0.0, 1.0, logP); // DNorm
        return np + f * nd;
    }
}; /* ppois_asymp() */

zS.poisson.cdf.ppois_asymp.dpnorm = function (x, lowerTail, lp) {
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

    if (x < 0) {
        x = -x;
        lowerTail = !lowerTail;
    }

    if (x > 10 && !lower_tail) {
        var
            sum = term = 1.0 / x,
            x2 = x * x,
            i = 1.0;

        do {
            term *= -i / x2;
            sum += term;
            i += 2.0;
        } while (Math.abs(term) > zU.CONST.DBL_EPSILON * sum);

        return 1.0 / sum;
    } else {
        var d = zygotine.S.normal.pdf(x, 0.0, 1.0, false);
        return d / Math.exp(lp);
    }
};

zS.poisson.cdf.ppois_asymp.coefs_a = [
    -1e99, /* placeholder used for 1-indexing */
    2.0 / 3.0,
	-4 / 135.0,
    8.0 / 2835.0,
    16.0 / 8505.0,
    -8992.0 / 12629925.0,
    -334144.0 / 492567075.0,
    698752 / 1477701225.0];

zS.poisson.cdf.ppois_asymp.coefs_b = [
    -1e99, /* placeholder */
    1.0 / 12.0,
    1.0 / 288.0,
    -139.0 / 51840.0,
    -571.0 / 2488320.0,
    163879.0 / 209018880.0,
    5246819.0 / 75246796800.0,
    -534703531 / 902961561600.0];

/***> loi Gamma <***/

zS.gamma = {};

/* R file: pgamma.c
 * Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). 
 * 
 * fonction utilisée tant pour pgamma que qgamma.
 */
zS.gamma.lgamma1p = function (a) {


    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */

    var CONST = zS.gamma.lgamma1p.CONST;
    var
        c = CONST.c, // 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
        euler = CONST.eulers_const,
        N = CONST.N,
        tol_logcf = CONST.tol_logcf; //  1e-14;

    var coeffs = zS.gamma.lgamma1p.coeffs;
    var /*double*/ lgam;
    var /*int*/ i;

    if (Math.abs(a) >= 0.5) {
        return zygotine.Num.logGamma(a + 1);
    }

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    lgam = c * zS.gamma.logcf(-a / 2, N + 2, 1, tol_logcf);
    for (i = N - 1; i >= 0; i--) {
        lgam = coeffs[i] - a * lgam;
    }

    return (a * lgam - euler) * a - zS.gamma.log1pmx(a);
}; /* lgamma1p */

/* 
* R file: pgamma.c
*/
zS.gamma.lgamma1p.CONST = {
    c: 0.2273736845824652515226821577978691e-12, /* zeta(N+2)-1 */
    eulers_const: 0.5772156649015328606065120900824024,
    N: 40,
    tol_logcf: 1e-14
};

zS.gamma.lgamma1p.coeffs = [
    0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
    0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
    0.2058080842778454787900092413529198e-1,
    0.7385551028673985266273097291406834e-2,
    0.2890510330741523285752988298486755e-2,
    0.1192753911703260977113935692828109e-2,
    0.5096695247430424223356548135815582e-3,
    0.2231547584535793797614188036013401e-3,
    0.9945751278180853371459589003190170e-4,
    0.4492623673813314170020750240635786e-4,
    0.2050721277567069155316650397830591e-4,
    0.9439488275268395903987425104415055e-5,
    0.4374866789907487804181793223952411e-5,
    0.2039215753801366236781900709670839e-5,
    0.9551412130407419832857179772951265e-6,
    0.4492469198764566043294290331193655e-6,
    0.2120718480555466586923135901077628e-6,
    0.1004322482396809960872083050053344e-6,
    0.4769810169363980565760193417246730e-7,
    0.2271109460894316491031998116062124e-7,
    0.1083865921489695409107491757968159e-7,
    0.5183475041970046655121248647057669e-8,
    0.2483674543802478317185008663991718e-8,
    0.1192140140586091207442548202774640e-8,
    0.5731367241678862013330194857961011e-9,
    0.2759522885124233145178149692816341e-9,
    0.1330476437424448948149715720858008e-9,
    0.6422964563838100022082448087644648e-10,
    0.3104424774732227276239215783404066e-10,
    0.1502138408075414217093301048780668e-10,
    0.7275974480239079662504549924814047e-11,
    0.3527742476575915083615072228655483e-11,
    0.1711991790559617908601084114443031e-11,
    0.8315385841420284819798357793954418e-12,
    0.4042200525289440065536008957032895e-12,
    0.1966475631096616490411045679010286e-12,
    0.9573630387838555763782200936508615e-13,
    0.4664076026428374224576492565974577e-13,
    0.2273736960065972320633279596737272e-13,
    0.1109139947083452201658320007192334e-13
];

/*

With a shape parameter κ and a scale parameter θ.
With a shape parameter α = k and an inverse scale parameter β = 1/θ, called a rate parameter.
With a shape parameter k and a mean parameter μ = k/β.

*/
/*
 * Fichier R : pgamma.c
 */
zS.gamma.cdf = function (x, alpha, scale, lowerTail, logP) {

    var isnan = zygotine.isnan; // exige qu'on ait affaire à un nombre et que pour celui-ci la fonction isNaN ramène true.
    var isUndefined = zygotine.isUndefined;
    var isABoolean = zygotine.isABoolean;
    var isANumber = zygotine.isANumber; // attention isANumber(NaN) est vrai car typeof(NaN) === 'number'!

    if (!isABoolean(logP)) {
        logP = false;
    }

    if (!isABoolean(lowerTail)) {
        lowerTail = true;
    }

    if (!isANumber(scale)) {
        scale = 1.0;
    }

    if (isnan(x) || isnan(alpha) || isnan(scale)) {
        return NaN;
    }

    if ((alpha < 0.0) || (scale <= 0.0)) {
        return NaN;
    }

    x /= scale;
    if (alpha === 0.0) /* limit case; useful e.g. in pnchisq() */ {
        return (x <= 0) ? (lowerTail ? (logP ? -Infinity : 0.0) : (logP ? 0.0 : 1.0)) : (lowerTail ? (logP ? 0.0 : 1.0) : (logP ? -Infinity : 0.0)); /* <= assert  pgamma(0,0) ==> 0 */
    }

    return zS.gamma.cdf.pgamma_raw(x, alpha, lowerTail, logP);
};

zS.gamma.cdf.usingScaleParameter = zS.gamma.cdf;

zS.gamma.cdf.usingRateParameter = function (x, alpha, rate, lowerTail, logP) {

    if (typeof rate === 'undefined') {
        rate = 1.0;
    }

    return zS.gamma.cdf.usingScaleParameter(x, alpha, 1.0 / rate, lowerTail, logP);
};

zS.gamma.cdf.scalefactor = ((4294967296.0 * 4294967296.0) * (4294967296.0 * 4294967296.0)) * ((4294967296.0 * 4294967296.0) * (4294967296.0 * 4294967296.0));

zS.gamma.cdf.pgamma_raw = function (x, alph, lowerTail, logP) {

    /* Here, assume that  (x,alph) are not NA  &  alph > 0 . */
    var res;
    //#define R_P_bounds_01(x, x_min, x_max) 	\
    //if (x <= x_min) return R_DT_0;		\
    //if (x >= x_max) return R_DT_1

    //R_P_bounds_01(x, 0.0, ML_POSINF);
    if (x <= 0.0) {
        return (lowerTail ? (logP ? -Infinity : 0.0) : (logP ? 0.0 : 1.0));
    }

    if (x >= Infinity) {
        return (lowerTail ? (logP ? 0.0 : 1.0) : (logP ? -Infinity : 0.0));
    }

    var sum, d;
    if (x < 1) {
        res = zS.gamma.cdf.pgamma_smallx(x, alph, lowerTail, logP);
    }
    else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
        /* incl. large alph compared to x */
        sum = zS.gamma.cdf.pd_upper_series(x, alph, logP); /* = x/alph + o(x/alph) */
        d = zS.poisson.pdf.dpois_wrap(alph, x, logP);
        if (!lowerTail) {
            res = logP
            ? zS.R_Log1_Exp(d + sum)
            : 1 - d * sum;
        }
        else {
            res = logP ? sum + d : sum * d;
        }
    }
    else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
        /* incl. large x compared to alph */
        d = zS.poisson.pdf.dpois_wrap(alph, x, logP);

        if (alph < 1) {
            if (x * zU.CONST.DBL_EPSILON > 1 - alph) {
                sum = (logP ? 0.0 : 1.0);
            }
            else {
                var f = zS.gamma.cdf.pd_lower_cf(alph, x - (alph - 1)) * x / alph;
                /* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
                sum = logP ? Math.log(f) : f;
            }
        }
        else {
            sum = zS.gamma.cdf.pd_lower_series(x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
            sum = logP ? Math.log1p(sum) : 1 + sum;
        }

        if (!lowerTail) {
            res = logP ? sum + d : sum * d;
        }
        else {
            if (logP) {
                let y = d + sum;
                res = y > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(y)) : Math.log1p(-Math.exp(y));
            } else {
                res = 1 - d * sum;
            }
        }
    }
    else { /* x >= 1 and x fairly near alph. */
        res = zS.poisson.cdf.ppois_asymp(alph - 1, x, !lowerTail, logP);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.	 In those
     * cases, simply redo via log space.
     */
    if (!logP && res < zU.CONST.DBL_MIN / zU.CONST.DBL_EPSILON) {
        /* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
        return Math.exp(zS.gamma.cdf.pgamma_raw(x, alph, lowerTail, true));
    }
    else {
        return res;
    }
}; // END OF pgamma_raw : private static double

/* R file: pgamma.c
 * 
 * Abramowitz and Stegun 6.5.29 [right]
 */
zS.gamma.cdf.pgamma_smallx = function (x, alph, lowerTail, logP) {
    var sum = 0, c = alph, n = 0, term;
    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
        n++;
        c *= -x / n;
        term = c / (alph + n);
        sum += term;
    } while (Math.abs(term) > zU.CONST.DBL_EPSILON * Math.abs(sum));

    if (lowerTail) {
        var f1 = logP ? Math.log1p(sum) : 1 + sum;
        var f2;
        if (alph > 1) {
            f2 = zS.poisson.pdf.dpois_raw(alph, x, logP);
            f2 = logP ? f2 + x : f2 * Math.exp(x);
        }
        else if (logP) {
            f2 = alph * Math.log(x) - zS.gamma.lgamma1p(alph);
        }
        else {
            f2 = Math.pow(x, alph) / Math.exp(zS.gamma.lgamma1p(alph));
        }

        return logP ? f1 + f2 : f1 * f2;
    }
    else {
        var lf2 = alph * Math.log(x) - zS.gamma.lgamma1p(alph);
        if (logP) {
            return zS.R_Log1_Exp(Math.log1p(sum) + lf2);
        }
        else {
            var f1m1 = sum;
            var f2m1 = Math.expm1(lf2);
            return -(f1m1 + f2m1 + f1m1 * f2m1);
        }
    }
}; /*  END OF pgamma_smallx() */

/* R file: pgamma.c
 * 
 */
zS.gamma.cdf.pd_upper_series = function (x, y, logP) {
    var term = x / y;
    var sum = term;

    do {
        y++;
        term *= x / y;
        sum += term;
    } while (term > sum * zU.CONST.DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
     *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return logP ? Math.log(sum) : sum;
};  /*  END OF zS.gamma.cdf.pd_upper_series */

/* R file: pgamma.c
 * 
 * Continued fraction for calculation of scaled upper-tail F_{gamma}
 * ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
zS.gamma.cdf.pd_lower_cf = function (y, d) {

    var scalefactor = zS.gamma.cdf.scalefactor;

    var f = 0.0 /* -Wall */, of, f0;
    var i, c2, c3, c4, a1, b1, a2, b2;
    var max_it = 200000;

    if (y === 0) {
        return 0;
    }

    f0 = y / d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */

    if (Math.abs(y - 1) < Math.abs(d) * zU.CONST.DBL_EPSILON) { /* includes y < d = Inf */
        return (f0);
    }

    if (f0 > 1.0) {
        f0 = 1.0;
    }

    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while (b2 > scalefactor) {
        a1 /= scalefactor;
        b1 /= scalefactor;
        a2 /= scalefactor;
        b2 /= scalefactor;
    }

    i = 0; of = -1.0; /* far away */
    while (i < max_it) {

        i++; c2--; c3 = i * c2; c4 += 2;
        /* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
        a1 = c4 * a2 + c3 * a1;
        b1 = c4 * b2 + c3 * b1;

        i++; c2--; c3 = i * c2; c4 += 2;
        /* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
        a2 = c4 * a1 + c3 * a2;
        b2 = c4 * b1 + c3 * b2;

        if (b2 > scalefactor) {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        }

        if (b2 !== 0) {
            f = a2 / b2;
            /* convergence check: relative; "absolute" for very small f : */
            if (Math.abs(f - of) <= zU.CONST.DBL_EPSILON * Math.max(f0, Math.abs(f))) {
                return f;
            }

            of = f;
        }
    }

    // MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n", f);
    return f;/* should not happen ... */
}; /* pd_lower_cf() */

zS.gamma.logcf = function (x, i, d, eps /* ~ relative tolerance */) {
    var /*double*/
        c1 = 2 * d,
        c2 = i + d,
        c4 = c2 + d,
        a1 = c2,
        b1 = i * (c2 - i * x),
        b2 = d * d * x,
        a2 = c4 * c2 - b2;

    var c3;
    var scalefactor = zS.gamma.cdf.scalefactor;

    b2 = c4 * b1 - i * b2;

    while (Math.abs(a2 * b1 - a1 * b2) > Math.abs(eps * b1 * b2)) {
        c3 = c2 * c2 * x;
        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;

        c3 = c1 * c1 * x;
        c1 += d;
        c4 += d;
        a2 = c4 * a1 - c3 * a2;
        b2 = c4 * b1 - c3 * b2;
        if (Math.abs(b2) > scalefactor) {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        }
        else if (Math.abs(b2) < 1 / scalefactor) {
            a1 *= scalefactor;
            b1 *= scalefactor;
            a2 *= scalefactor;
            b2 *= scalefactor;
        }
    }

    return a2 / b2;
};

/* R file: pgamma.c
 * 
 */
zS.gamma.cdf.pd_lower_series = function (lambda, y) {
    var term = 1, sum = 0;

    while (y >= 1 && term > sum * zU.CONST.DBL_EPSILON) {
        term *= y / lambda;
        sum += term;
        y--;
    }

    /* 
     * sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     * y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
     * y/lambda + o(y/lambda)
     */

    if (y !== Math.floor(y)) {
        /*
        * The series does not converge as the terms start getting
        * bigger (besides flipping sign) for y < -lambda.
        */
        var f;
        /* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
         *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
        f = zS.gamma.cdf.pd_lower_cf(y, lambda + 1 - y);
        sum += term * f;
    }

    return sum;
}; /* pd_lower_series() */

zS.gamma.cdf.minLog1Value = -0.79149064;

zS.gamma.cdf.two = 2.0;

zS.gamma.cdf.tol_logcf = 1e-14;

zS.gamma.icdf = function (p, shape, scale, lowerTail, logP) {

    var isnan = zygotine.isnan; // exige qu'on ait affaire à un nombre et que pour celui-ci la fonction isNaN ramène true.
    var isUndefined = zygotine.isUndefined;
    var isABoolean = zygotine.isABoolean;
    var isANumber = zygotine.isANumber; // attention isANumber(NaN) est vrai car typeof(NaN) === 'number'!
    var cntn = true;

    if (!isABoolean(logP)) {
        logP = false;
    }

    if (!isABoolean(lowerTail)) {
        lowerTail = true;
    }

    if (!isANumber(scale)) {
        scale = 1.0;
    }

    if (isnan(p) || isnan(shape) || isnan(scale)) {
        return NaN;
    }

    if ((shape < 0.0) || (scale <= 0.0)) {
        return NaN;
    }

    var
        i420 = 1.0 / 420.0,
        i2520 = 1.0 / 2520.0,
        i5040 = 1.0 / 5040.0;

    var
        CONST = zS.gamma.icdf.CONST;

    var
        EPS1 = CONST.EPS1;
    EPS2 = CONST.EPS2,
    EPS_N = CONST.EPS_N,
    MAXIT = CONST.MAXIT,
    pMIN = CONST.pMIN,
    pMAX = CONST.pMAX;

    var p_, a, b, c, g, ch, ch0, p1;
    var p2, q, s1, s2, s3, s4, s5, s6, t, x;
    var /*int*/ i, max_it_Newton = 1;

    //R_Q_P01_boundaries(p, 0.0, ML_POSINF);
    if (logP) {
        if (p > 0.0) {
            return NaN;
        }

        if (p === 0) { /* upper bound*/
            return lowerTail ? -Infinity : 0.0;
        }

        if (p === -Infinity) {
            return lowerTail ? 0.0 : -Infinity;
        }
    }
    else { /* !log_p */
        if ((p < 0) || (p > 1)) {
            return NaN;
        }

        if (p === 0) {
            return lowerTail ? 0.0 : -Infinity;
        }

        if (p === 1) {
            return lowerTail ? -Infinity : 0.0;
        }
    }

    if (shape < 0 || scale <= 0) {
        return NaN; //ML_ERR_return_NAN;
    }

    if (shape === 0) { /* all mass at 0 : */
        return 0.0;
    }

    if (shape < 1e-10) {
        /* Warning seems unnecessary now: */
        //MATHLIB_WARNING("value of shape (%g) is extremely small: results may be unreliable", alpha);
        max_it_Newton = 7;/* may still be increased below */
    }

    p_ = zS.R_DT_qIv(p, lowerTail, logP);/* lower_tail prob (in any case) */

    g = zygotine.Num.logGamma(shape);/* log Gamma(v/2) */

    /*----- Phase I : Starting Approximation */
    ch = zS.gamma.icdf.qchisq_appr(p, 2 * shape, g, lowerTail, logP, EPS1);
    if (!isFinite(ch)) {
        /* forget about all iterations! */
        max_it_Newton = 0;
        cntn = false;
    }

    // la boucle ne sera exécuté qu'une fois! Pour gérer le goto END du langage C.
    while (cntn) {
        if (ch < EPS2) {/* Corrected according to AS 91; MM, May 25, 1999 */
            max_it_Newton = 20;
            cntn = false;
            break;/* and do Newton steps */
        }

        /* FIXME: This (cutoff to {0, +Inf}) is far from optimal
         * -----  when log_p or !lower_tail, but NOT doing it can be even worse */
        if (p_ > pMAX || p_ < pMIN) {
            /* did return ML_POSINF or 0.;	much better: */
            max_it_Newton = 20;
            break;/* and do Newton steps */
        }

        /*----- Phase II: Iteration
         *	Call pgamma() [AS 239]	and calculate seven term taylor series
         */
        c = shape - 1;
        s6 = (120 + c * (346 + 127 * c)) * i5040; /* used below, is "const" */

        ch0 = ch;/* save initial approx. */
        for (i = 1; i <= MAXIT; i++) {
            q = ch;
            p1 = 0.5 * ch;
            p2 = p_ - zS.gamma.cdf.pgamma_raw(p1, shape, /*lower_tail*/true, /*log_p*/false);
            if (!isFinite(p2) || ch <= 0) {
                ch = ch0;
                max_it_Newton = 27;
                break;
            }

            t = p2 * Math.exp(shape * zU.CONST.M_LN2 + g + p1 - c * Math.log(ch));
            b = t / ch;
            a = 0.5 * t - b * c;
            s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) * i420;
            s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) * i2520;
            s3 = (210 + a * (462 + a * (707 + 932 * a))) * i2520;
            s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) * i5040;
            s5 = (84 + 2264 * a + c * (1175 + 606 * a)) * i2520;
            ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
            if (Math.abs(q - ch) < EPS2 * ch) {
                break; // on sort du for
            }

            if (Math.abs(q - ch) > 0.1 * ch) { /* diverging? -- also forces ch > 0 */
                if (ch < q) {
                    ch = 0.9 * q;
                }
                else {
                    ch = 1.1 * q;
                }
            }
        } //for (i = 1; i <= MAXIT; i++)

        break; // on sort du while, pour arriver à l'étiquette END:
    }
    /* no convergence in MAXIT iterations -- but we add Newton now... */
    END:
    /* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
       --------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision
    
       * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
       *
       * Improved (MM): - only if rel.Err > EPS_N (= 1e-15);
       *		    - also for lower_tail = FALSE	 or log_p = TRUE
       * 		    - optionally *iterate* Newton
       */
        x = 0.5 * scale * ch;
    if (max_it_Newton > 0) {
        /* always use log scale */
        if (!logP) {
            p = Math.log(p);
            logP = true;
        }
        if (x === 0) {
            x = zU.CONST.DBL_MIN;
            p_ = zS.gamma.cdf.usingScaleParameter(x, shape, scale, lowerTail, logP);
            if ((lowerTail && (p_ > p * (1.0 + 1e-7))) || (!lowerTail && (p_ < p * (1.0 - 1e-7)))) {
                return (0.0);
            }
            /* else:  continue, using x = DBL_MIN instead of  0  */
        }
        else {
            p_ = zS.gamma.cdf.usingScaleParameter(x, shape, scale, lowerTail, logP);
        }

        if (p_ === -Infinity) {
            return 0; /* PR#14710 */
        }

        for (i = 1; i <= max_it_Newton; i++) {
            p1 = p_ - p;
            if (Math.abs(p1) < Math.abs(EPS_N * p)) {
                break;
            }

            /* else */
            // dgamma (l'affectation de g était dans la condition if ci-dessous.
            g = zS.gamma.pdf.usingScaleParameter(x, shape, scale, logP);
            if (g === (logP ? -Infinity : 0.0)) {
                break;
            }

            /* else :
             * delta x = f(x)/f'(x);
             * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
             * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))
             */
            t = logP ? p1 * Math.exp(p_ - g) : p1 / g;/* = "delta x" */
            t = lowerTail ? x - t : x + t;
            p_ = zS.gamma.cdf.usingScaleParameter(t, shape, scale, lowerTail, logP);
            if (
                (Math.abs(p_ - p) > Math.abs(p1)) ||
                (i > 1 && Math.abs(p_ - p) === Math.abs(p1))) {
                /* <- against flip-flop */
                /* no improvement */
                break;
            } /* else : */
            //# ifdef Harmful_notably_if_max_it_Newton_is_1
            /* control step length: this could have started at
               the initial approximation */
            if (t > 1.1 * x) {
                t = 1.1 * x;
            } else if (t < 0.9 * x) {
                t = 0.9 * x;
            }

            x = t;
        }
    } //if (max_it_Newton > 0)

    return x;
};

zS.gamma.icdf.usingScaleParameter = zS.gamma.icdf;

zS.gamma.icdf.usingRateParameter = function (p, shape, rate, lowerTail, logP) {

    if (zygotine.isUndefined(rate)) {
        rate = 1.0;
    }

    return zS.gamma.icdf.usingScaleParameter(p, shape, 1.0 / rate, lowerTail, logP);
};

/* 
* R file: dgamma.c
*/
zS.gamma.pdf = function (x, shape, scale, giveLog) {
    var isnan = zygotine.isnan; // exige qu'on ait affaire à un nombre et que pour celui-ci la fonction isNaN ramène true.
    var isUndefined = zygotine.isUndefined;
    var isABoolean = zygotine.isABoolean;
    var isANumber = zygotine.isANumber; // attention isANumber(NaN) est vrai car typeof(NaN) === 'number'!

    if (!isABoolean(giveLog)) {
        giveLog = false;
    }

    if (!isANumber(scale)) {
        scale = 1.0;
    }

    if (isnan(x) || isnan(shape) || isnan(scale)) {
        return NaN;
    }

    if ((shape < 0.0) || (scale <= 0.0)) {
        return NaN;
    }

    if (x < 0) {
        return (giveLog ? -Infinity : 0.0);
    }

    if (shape === 0) {
        /* point mass at 0 */
        return (x === 0) ? -Infinity : (giveLog ? -Infinity : 0.0);
    }

    if (x === 0) {
        if (shape < 1) {
            return -Infinity;
        }

        if (shape > 1) {
            return (giveLog ? -Infinity : 0.0);
        }

        /* else */
        return giveLog ? -Math.log(scale) : 1 / scale;
    }

    if (shape < 1) {
        pr = zS.poisson.pdf.dpois_raw(shape, x / scale, giveLog);
        return giveLog ? pr + Math.log(shape / x) : pr * shape / x;
    }
    /* else  shape >= 1 */
    pr = zS.poisson.pdf.dpois_raw(shape - 1, x / scale, giveLog);
    return giveLog ? pr - Math.log(scale) : pr / scale;
};

zS.gamma.pdf.usingScaleParameter = zS.gamma.pdf;

zS.gamma.pdf.usingRateParameter = function (x, shape, rate, giveLog) {

    if (zygotine.isUndefined(rate)) {
        rate = 1.0;
    }

    return zS.gamma.pdf.usingScaleParameter(x, shape, 1 / rate, giveLog);
};

zS.gamma.icdf.CONST = {
    EPS1: 1e-2,
    EPS2: 5e-7,    /* final precision of AS 91 */
    EPS_N: 1e-15,  /* precision of Newton step / iterations */
    MAXIT: 1000,   /* was 20 */
    pMIN: 1e-100,  /* was 0.000002 : 2e-6 */
    pMAX: (1 - 1e-14)
};

zS.gamma.icdf.qchisq_appr = function (p, nu, g, lowerTail, logP, tol) {

    const
        C7 = 4.67,
        C8 = 6.66,
        C9 = 6.73,
        C10 = 13.32;


    var /*double*/
        alpha, a, c, ch, p1, p2, q, t, x;

    /* test arguments and initialise */

    if (isNaN(p) || isNaN(nu)) {
        return p + nu;
    }

    if ((logP && (p > 0)) ||
        (!logP && ((p < 0) || (p > 1)))) {
        return NaN; //ML_ERR_return_NAN
    }

    //R_Q_P01_check(p);
    if (nu <= 0) {
        return double.NaN; //ML_ERR_return_NAN;
    }

    alpha = 0.5 * nu;/* = [pq]gamma() shape */
    c = alpha - 1;
    //p1 = (lower_tail ? (log_p ? (p) : Math.Log(p)) : ((p) > -M_LN2 ? Math.Log(-expm1(p)) : log1p(-Math.Exp(p))));
    p1 = (lowerTail ? (logP ? (p) : Math.log(p)) : (logP ? ((p) > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(p)) : Math.log1p(-Math.exp(p))) : Math.log1p(-p)));
    if (nu < (-1.24 * p1)) {   /* for small chi-squared */
        /* Math.Log(alpha) + g = Math.Log(alpha) + Math.Log(gamma(alpha)) =
         *        = Math.Log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
         *  catastrophic cancellation when alpha << 1
         */
        let lgam1pa = (alpha < 0.5) ? zS.gamma.lgamma1p(alpha) : (Math.log(alpha) + g);
        ch = Math.exp((lgam1pa + p1) / alpha + zU.CONST.M_LN2);
    }
    else if (nu > 0.32) {   /*  using Wilson and Hilferty estimate */

        x = zS.normal.icdf(p, 0, 1, lowerTail, logP);
        p1 = 2.0 / (9 * nu);
        ch = nu * Math.pow(x * Math.sqrt(p1) + 1 - p1, 3);

        /* approximation for p tending to 1: */
        if (ch > 2.2 * nu + 6) {
            // tmp === R_DT_Clog(p, lower_tail, log_p)
            let tmp =
                lowerTail ?
                (logP ? (p > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(p)) : Math.log1p(-Math.exp(p))) : Math.log1p(-p)) :
                (logP ? p : Math.log(p));
            ch = -2 * (tmp - c * Math.log(0.5 * ch) + g);
        }
    }
    else { /* "small nu" : 1.24*(-Math.Math.Log(p)) <= nu <= 0.32 */
        ch = 0.4;
        // tmp === R_DT_Clog(p, lower_tail, log_p)
        let tmp =
            lowerTail ?
            (logP ? (p > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(p)) : Math.log1p(-Math.exp(p))) : Math.log1p(-p)) :
            (logP ? p : Math.log(p));
        a = tmp + g + c * zU.CONST.M_LN2;

        do {
            q = ch;
            p1 = 1.0 / (1 + ch * (C7 + ch));
            p2 = ch * (C9 + ch * (C8 + ch));
            t = -0.5 + (C7 + 2 * ch) * p1 - (C9 + ch * (C10 + 3 * ch)) / p2;
            ch -= (1 - Math.exp(a + 0.5 * ch) * p2 * p1) / t;
        } while (Math.abs(q - ch) > tol * Math.abs(ch));
    }

    return ch;
}; /* END of zS.gamma.icdf.qchisq_appr */

/* R file: pgamma.c
 * Accurate calculation of log(1+x)-x, particularly for small x.
 */
zS.gamma.log1pmx = function (x) {
    if (x > 1.0 || x < -0.79149064) {
        return Math.log1p(x) - x;
    }
    else {
        /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
         * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
         * ---------------------------------------------
         * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
         */
        var
            r = x / (2.0 + x),
            y = r * r;

        if (Math.abs(x) < 1.0e-2) {

            return r * ((((2.0 / 9.0 * y + 2.0 / 7.0) * y + 2.0 / 5.0) * y + 2.0 / 3.0) * y - x);
        }
        else {
            return r * (2 * y * zS.gamma.logcf(y, 3, 2, zS.gamma.cdf.tol_logcf) - x);
        }
    }
};


zS.gamma.validateShapeAndScale = function (shape, scale) {

    // On assume, pour l'instant, que shape et scale, si définis, sont des nombres au sens stricts

    //if (ISNAN(a) || ISNAN(scale))
    //    ML_ERR_return_NAN;
    //if (a <= 0.0 || scale <= 0.0) {
    //    if (scale == 0. || a == 0.) return 0.;
    //    ML_ERR_return_NAN;
    //}
    //if (!R_FINITE(a) || !R_FINITE(scale)) return ML_POSINF;

    var reponse = { retval: undefined, throwMsg: null, isOK: false };

    var isnan = zygotine.isnan; // isnan(x) est vrai si x n'es pas un nombre; dans le cas d'un nombre on retourne isNaN(x).


    if (zygotine.isUndefined(scale)) {
        scale = 1.0;
    }

    if (zygotine.isUndefined(shape)) {
        reponse.msg = "Parameter shape is required.";
    } else if (isnan(shape) || isnan(scale)) {
        reponse.retval = NaN; //ML_ERR_return_NAN;
    } else {
        if ((shape <= 0.0) || (scale <= 0.0)) {
            if ((shape === 0.0) || (scale === 0.0)) {
                reponse.retval = 0.0;
            } else {
                reponse.retval = NaN; // ML_ERR_return_NAN
            }
        } else {
            if (!isFinite(shape) || !isFinite(scale)) {
                reponse.retval = Infinity;
            } else {
                reponse.isOK = true;
            }
        }
    }

    return reponse;
};


//publique
zS.gamma.sample = function (n, shape, scale) {

    if (zygotine.isUndefined(n)) {
        // test sans grand intérêt car shape est requis et n'a pas de valeur par défaut.
        n = 1;
    }

    /* 
    
    tous les paramètres fournies pour lesquels
    pour l'instant on assume que si n est défini (non undefined), il s'agit d'un nombre (Number).

    if (!zygotine.isANumber(n)) {
        throw new Error("Parameter n is required.");
    }

    */

    //if (ISNAN(a) || ISNAN(scale))
    //    ML_ERR_return_NAN;
    //if (a <= 0.0 || scale <= 0.0) {
    //    if (scale == 0. || a == 0.) return 0.;
    //    ML_ERR_return_NAN;
    //}
    //if (!R_FINITE(a) || !R_FINITE(scale)) return ML_POSINF;

    // On doit avoir: shape et scale finis, shape >= 0 et scale > 0
    let test = zS.gamma.validateShapeAndScale(shape, scale);
    if (test.throwMsg !== null) {
        throw new Error(test.throwMsg);
    } else {
        var rep = [];
        if (!test.isOK) {
            rep = Array.apply(null, { length: n }).map(xx => test.retval);
        } else {
            let rand = zS.gamma.rand;
            rep = Array.apply(null, { length: n }).map(xx => rand(shape, scale));
        }

        return (n === 1) ? rep[0] : rep;
    }
};

zS.gamma.sample.usingScaleParameter = zS.gamma.sample;

zS.gamma.sample.usingRateParameter = function (n, shape, scale) {
    return zS.gamma.sample(n, shape, 1 / scale);
};

zS.gamma.rand = function (shape, scale) {

    //if (!R_FINITE(a) || !R_FINITE(scale) || a < 0.0 || scale <= 0.0)
    //{
    //    if (scale == 0.0) return 0.0;
    //    return double.NaN; //ML_ERR_return_NAN;
    //}
    /* State variables [FIXME for threading!] :*/
    var /*double*/
        aa = 0.0, aaa = 0.0,
        s = 0, s2 = 0, d = 0,    /* no. 1 (step 1) */
        q0 = 0, b = 0, si = 0, c = 0;/* no. 2 (step 4) */
    var
        e, p, q, r, t, u, v, w, x, ret_val;
    var
        CONST = zS.gamma.rand.CONST;
    var
        Q = CONST.q,
        A = CONST.a;

    if (shape < 1.0) {
        /* GS algorithm for parameters a < 1 */
        if (shape === 0.0) {
            return 0.0;
        }

        e = 1.0 + CONST.exp_m1 * shape;
        for (; ;) {
            p = e * zS.uniform.rand(); // interval (0,1)
            if (p >= 1.0) {
                x = -Math.log((e - p) / shape);
                if (zS.exponential.rand() >= (1.0 - shape) * Math.log(x)) {
                    break;
                }
            }
            else {
                x = Math.exp(Math.log(p) / shape);
                if (zS.exponential.rand() >= x) {
                    break;
                }
            }
        }

        return scale * x;
    }

    /* --- a >= 1 : GD algorithm --- */

    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (shape !== aa) {
        aa = shape;
        s2 = shape - 0.5;
        s = Math.sqrt(s2);
        d = CONST.sqrt32 - s * 12.0;
    }
    /* Step 2: t = standard normal deviate,
               x = (s,1/2) -normal deviate. */

    /* immediate acceptance (i) */
    t = zS.normal.rand();
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0) {
        return scale * ret_val;
    }

    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = zS.uniform.rand();
    if (d * u <= t * t * t) {
        return scale * ret_val;
    }

    /* Step 4: recalculations of q0, b, si, c if necessary */

    if (shape !== aaa) {
        aaa = shape;
        r = 1.0 / shape;
        q0 = ((((((Q[7] * r + Q[6]) * r + Q[5]) * r + Q[4]) * r + Q[3]) * r
               + Q[2]) * r + Q[1]) * r;

        /* Approximation depending on size of parameter a */
        /* The constants in the expressions for b, si and c */
        /* were established by numerical experiments */

        if (shape <= 3.686) {
            b = 0.463 + s + 0.178 * s2;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.16 * s;
        } else if (shape <= 13.022) {
            b = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
            c = 0.062 / s + 0.024;
        } else {
            b = 1.77;
            si = 0.75;
            c = 0.1515 / s;
        }
    } // if (shape != aaa)
    /* Step 5: no quotient test if x not positive */

    if (x > 0.0) {
        /* Step 6: calculation of v and quotient q */
        v = t / (s + s);
        if (Math.abs(v) <= 0.25)
            q = q0 + 0.5 * t * t * ((((((A[7] * v + A[6]) * v + A[5]) * v + A[4]) * v
                          + A[3]) * v + A[2]) * v + A[1]) * v;
        else {
            q = q0 - s * t + 0.25 * t * t + (s2 + s2) * Math.log(1.0 + v);
        }

        /* Step 7: quotient acceptance (q) */
        if (Math.log(1.0 - u) <= q) {
            return scale * ret_val;
        }
    }

    while (true) {
        /* Step 8: e = standard exponential deviate
         *	u =  0,1 -uniform deviate
         *	t = (b,si)-double exponential (laplace) sample */
        e = zS.exponential.rand();
        u = zS.uniform.rand();
        u = u + u - 1.0;
        if (u < 0.0) {
            t = b - si * e;
        }
        else {
            t = b + si * e;
        }

        /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
        if (t >= -0.71874483771719) {
            /* Step 10:	 calculation of v and quotient q */
            v = t / (s + s);
            if (Math.abs(v) <= 0.25) {
                q = q0 + 0.5 * t * t *
                    ((((((A[7] * v + A[6]) * v + A[5]) * v + A[4]) * v + A[3]) * v
                      + A[2]) * v + A[1]) * v;
            } else {
                q = q0 - s * t + 0.25 * t * t + (s2 + s2) * Math.log(1.0 + v);
            }

            /* Step 11:	 hat acceptance (h) */
            /* (if q not positive go to step 8) */
            if (q > 0.0) {
                w = Math.expm1(q);
                /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
                /* if t is rejected sample again at step 8 */
                if (c * Math.abs(u) <= w * Math.exp(e - 0.5 * t * t)) {
                    break;
                }
            }
        } // if (t >= -0.71874483771719)
    } /* repeat .. until  `t' is accepted */

    x = s + 0.5 * t;
    return scale * x * x;
}; // END of zS.gamma.rand

zS.gamma.rand.CONST = {
    sqrt32: 5.656854,
    exp_m1: 0.36787944117144232159,/* exp(-1) : 1/e */
    q: [
        NaN,
        0.04166669,
        0.02083148,
        0.00801191,
        0.00144121,
        -7.388e-5,
        2.4511e-4,
        2.424e-4],

    a: [
        NaN,
        0.3333333,
        -0.250003,
        0.2000062,
        -0.1662921,
        0.1423657,
        -0.1367177,
        0.1233795]
};

/**
* * R file: quantile.R
* subset of quantile.defaultl function
*
* @constructor
*
* @param {number[]} probs - un tableau de nombres pris dans l'interval ]0, 1[;
*
* @property {number[]} probs
*/
zygotine.S.Quantile = function (probs) {
    this.probs = (typeof probs !== "undefined") ? probs : [.025, .05, .1, .25, .5, .75, .9, .95, .975];
    this.probs.filter(x => (isFinite(x) && (x > this.eps) && (x < 1 - this.eps)));
    probsOk = this.probs.length > 0;
    if (!probsOk) {
        throw new Error("Propability vector in zygotine.S.Quantile is not well defined!");
    }
};

zygotine.S.Quantile.prototype = {

    //private const double eps = 2.22045E-14;
    eps: 2.22045E-14,

    /**
    * @function 
    * @param {number[]} x
    */
    compute: function (x) {
        var n = x.length;
        var index = zU.product(this.probs, n - 1);
        var lo = index.map(Math.floor);
        var hi = index.map(Math.ceil);
        var y = new Array(n);
        for (var i = 0; i < n; i++) {
            y[i] = x[i];
        }

        y.sort(function (a, b) { return a - b; });
        var qs = lo.map(function (iLo) {
            return y[iLo];
        });

        var iTmp = [];
        for (let k = 0; k < index.length; k++) {
            if (index[k] > lo[k])
                iTmp.push(k);
        }

        iTmp.forEach(
            function (k) {
                let h = (index[k] - lo[k]);
                qs[k] = ((1 - h) * qs[k]) + h * y[hi[k]];
            });

        return qs;
    }
};

zygotine.S.R_D_Cval = function (p, lowerTail) {
    return (lowerTail ? (0.5 - (p) + 0.5) : (p));	/*  1 - p */
};

zygotine.S.R_DT_CIv = function (p, lowerTail, logP) {
    return (logP ? (lowerTail ? -Math.expm1(p) : Math.exp(p)) : zygotine.S.R_D_Cval(p, lowerTail));
};

zygotine.S.R_DT_qIv = function (p, lowerTail, logP) {
    return (logP ? (lowerTail ? Math.exp(p) : -Math.expm1(p)) : zygotine.S.R_D_Lval(p, lowerTail));
};

zygotine.S.R_D_Lval = function (p, lowerTail) {
    return (lowerTail ? (p) : (0.5 - (p) + 0.5));
};

// R file dpq.h  
// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
// #define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
zygotine.S.R_Log1_Exp = function (x) {
    return ((x) > -zU.CONST.M_LN2 ? Math.log(-Math.expm1(x)) : Math.log1p(-Math.exp(x)));
};

zygotine.S.sampleWithReplacement = function (p, n, replace) {

    var rep = [];
    if (zygotine.isUndefined(n)) {
        n = 1;
    }

    var totalMass = p.reduce(function (a, total) { return a + total; });
    var sortedX = p.map(function (a, i) { return { p: a, index: i }; });
    sortedX.sort(function (a, b) { return b.p - a.p; });
    var weightedU;
    var mass;
    var position;

    for (let i = 0; i < n; i++) {
        mass = 0.0;
        weightedU = totalMass * zS.uniform.sample(1);
        for (position = 0; position < p.length; position++) {
            mass += sortedX[position].p;
            if (weightedU <= mass) {
                break;
            }
        }

        rep.push({
            i: sortedX[position].index,
            p: sortedX[position].p
        });
    }

    return rep;

};

zygotine.S.sample = function (p, n, replace) {

    var rep = [];

    if (zygotine.isUndefined(replace)) {
        replace = true;
    }

    if (zygotine.isUndefined(n)) {
        n = 1;
    }

    if (p.length === 0) {
        return null;
    }

    if (!replace && (n > p.length)) {
        return null;
    }

    var totalMass = p.reduce(function (a, total) { return a + total; });
    var sortedX = p.map(function (a, i) { return { p: a, index: i }; });
    sortedX.sort(function (a, b) { return b.p - a.p; });
    var weightedU;
    var mass;
    var position;

    for (let i = 0; i < n; i++) {
        mass = 0.0;
        weightedU = totalMass * zS.uniform.sample(1);
        for (position = 0; position < sortedX.length; position++) {
            mass += sortedX[position].p;
            if (weightedU <= mass) {
                rep.push({
                    i: sortedX[position].index,
                    p: sortedX[position].p
                });
                if (!replace) {
                    totalMass -= sortedX[position].p;
                    sortedX.splice(position, 1); // on retire l'élément dont l'index est 'position'
                }

                break;
            }
        }
    }

    return rep;

};



var p = [0.154162841732613742, 0.446221745456568897, 0.740049691521562636, 0.873387157334946096, 0.263315018615685403, 0.086434275493957102, 0.533739388803951442, 0.723601803998462856, 0.014574962784536183, 0.350478273001499474];


