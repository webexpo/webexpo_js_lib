/* eslint valid-jsdoc: 0 */
/// <reference path="A.js" />

if (!zygotine) {
    zygotine = {};
}

zygotine.U = {};
zU = zygotine.U;

zU.CONST =
    {
        M_LN2: 0.693147180559945309417, // math.h
        DBL_MAX_EXP: 1024, // float.h ieee-754
        DBL_MAX: 1.7976931348623158e+308, // C float.h ieee-754
        DBL_MIN: 2.2250738585072014E-308, // C float.h  Smallest positive normalized double- voir FPExtreme
        DBL_MIN_EXP: -1021, // C float.h
        DBL_MANT_DIG: 53, // C float.h
        DBL_EPSILON: 2.220446049250313e-16, // float.h.  DBL_EPSILON est le plus petit double positif tel que 1.0 + x > 1.0.
        R_AccuracyInfo_eps: this.DBL_EPSILON,
        M_2PI: 6.283185307179586476925286766559, // Rmath.h
        M_1_SQRT_2PI: 0.398942280401432677939946059934, //Rmath.h ...	/* 1/sqrt(2pi) */
        M_LN_SQRT_2PI: 0.918938533204672741780329736406,	// log(sqrt(2*pi)) :: log(2*pi)/2  // Rmath.h
        M_LN_SQRT_PId2: 0.225791352644727432363097614947,	// log(sqrt(pi/2)) // RMath.h
        M_SQRT2: 1.41421356237309504880,
        M_SQRT1_2: 0.707106781186547524401, // 1/sqrt(2)  
        M_SQRT_32: 5.656854249492380195206754896838, // Rmath.h
        M_PI: 3.141592653589793238462643383279502884197169399375, // Constants.h 
        CHAR_BIT: 8,
        INT_MAX: 2147483647,
        FLT_RADIX: 2,
        M_LOG10_2: 0.301029995663981195213738894724

    };

var CONST = zU.CONST;

zU.hypot = function (a, b) {
    var p, r, s, t, tmp, u; // des nombres (double)

    if (isNaN(a) || isNaN(b)) {
        return NaN;
    }

    if (!isFinite(a) || !isFinite(b)) {
        return Infinity;
    }

    p = Math.max(Math.abs(a), Math.abs(b));
    if (p !== 0.0) {
        /* r = (min(|a|,|b|) / p) ^2 */
        tmp = Math.min(Math.abs(a), Math.abs(b)) / p;
        r = tmp * tmp;
        for (; ;) {
            t = 4.0 + r;
            if (Math.abs(r) < 2 * CONST.DBL_EPSILON) {
                break;
            }
            s = r / t;
            u = 1.0 + 2.0 * s;
            p *= u;

            /* r = (s / u)^2 * r */
            tmp = s / u;
            r *= tmp * tmp;
        }
    }

    return p;
};

zU.fmt = (function () {
    var reg = [], tmp;
    addRE(20);

    return function (fmt, args) {
        tmp = arguments[0];
        for (var i = 1; i < arguments.length; i++) {
            tmp = tmp.replace(reg[i - 1], arguments[i]);
        }

        return tmp;
    };

    function addRE(stop) {
        var i, tmp;
        for (i = reg.length; i < stop; i++) {
            reg[i] = new RegExp("\\{" + (i) + "\\}", "g");
        }
    }
})();

zU.fmtArray = function (arr) {
    var rep = [];
    rep.push("[");
    rep.push(arr.join(", "));
    rep.push("]");
    return rep.join("");
};

/**
@function zygotine.U.toR
@param {number[]} arr 
@returns {string}
*/
zU.toR = function (arr) {
    return "c(" + arr.join(",\n") + ")";
};

zU.sum = function (x) {
    if (!Array.isArray(x)) {
        return typeof x === 'number' ? x : NaN;
    }
    else {
        /** @type {array} */
        return x.reduce(function (a, b) { return a + b; }, 0);
    }
};

zU.mean = function (x) {
    if (!Array.isArray(x)) {
        return typeof x === 'number' ? x : NaN;
    }

    if (x.length === 0) {
        return NaN;
    }

    var sum = x.reduce(function (a, b) { return a + b; }, 0);
    return sum / x.length;
};

zygotine.U.variance = function (source) {
    if (!Array.isArray(source) || source.length < 2) {
        return NaN;
    }

    var
        n = 0,
        mean = 0.0,
        M2 = 0.0;

    for (let i = 0; i < source.length; i++) {
        let x = source[i];
        let delta = x - mean;
        mean += delta / (i + 1);
        M2 += delta * (x - mean);
    }

    return M2 / (source.length - 1);
};


zU.max = function (x) {
    if (!Array.isArray(x)) {
        return typeof x === 'number' ? x : NaN;
    }
    else {
        return x.reduce(function (a, b) {
            return Math.max(a, b);
        }, -Infinity);
    }
};

zU.min = function (x) {
    if (!Array.isArray(x)) {
        return typeof x === 'number' ? x : NaN;
    }
    else {
        return x.reduce(function (a, b) {
            return Math.min(a, b);
        }, Infinity);
    }
};

zU.add = function (x, y) {

    var xIsArray = Array.isArray(x);
    var yIsArray = Array.isArray(y);

    // 2 scalaires
    if (!xIsArray && !yIsArray) {
        return x + y;
    }

    // 1 scalaire et un tableau
    var nx, ny;
    if ((xIsArray && !yIsArray) || (!xIsArray && yIsArray)) {
        if (xIsArray) {
            nx = x;
            ny = y;
        }
        else {
            nx = y;
            ny = x;
        }
        if (nx.length === 0) {
            return [];
        }

        return nx.map(function (b) {
            return b + ny;
        });
    }


    // 2 tableaux
    // on additionne point à point en se basant sur la longueur du long et en bouclant sur les éléments du plus court en prenant un modulo.


    if (x.length >= y.length) {
        nX = x;
        nY = y;
    }
    else {
        nX = y;
        nY = x;
    }

    if ((x.length === 0) || (y.length === 0)) {
        return [];
    }

    var lY = nY.length;
    return nX.map(function (a, i) {
        return a + nY[i % lY];
    });
};

zU.sumSqr = function (x) {
    if (!Array.isArray(x)) {
        return Math.pow(x, 2);
    }
    else {
        /** @type {array} */
        return x.reduce(function (a, b) { return a + (b * b); }, 0);
    }
};

zU.log = function (x) {
    if (!Array.isArray(x)) {
        return Math.log(x);
    }
    else {
        return x.map(function (a) {
            return Math.log(a);
        });
    }
};

zU.substract = function (x, y) {

    var xIsArray = Array.isArray(x);
    var yIsArray = Array.isArray(y);
    var permuted;
    // 2 scalaires
    if (!xIsArray && !yIsArray) {
        return x - y;
    }

    // 1 scalaire et un tableau
    var nx, ny;
    var fn;
    if ((xIsArray && !yIsArray) || (!xIsArray && yIsArray)) {
        if (xIsArray) {
            nx = x;
            ny = y;
            permuted = false;
        }
        else {
            permuted = true;
            nx = y;
            ny = x;
        }

        if (nx.length === 0) {
            return [];
        }

        if (permuted) {
            fn = function (b) {
                return ny - b;
            };

        } else {
            fn = function (b) {
                return b - ny;
            };
        }

        return nx.map(fn);
    }

    // 2 tableaux
    // on soustrait point à point en se basant sur la longueur du + long et en bouclant sur les éléments du plus court en prenant un modulo.
    if ((x.length === 0) || (y.length === 0)) {
        return [];
    }

    if (x.length >= y.length) {
        nX = x;
        nY = y;
        permuted = false;
    }
    else {
        nX = y;
        nY = x;
        permuted = true;
    }

    var lY = nY.length;

    fn = permuted ? function (a, i) { return nY[i % lY] - a; } : function (a, i) { return a - nY[i % lY]; };

    return nX.map(fn);
};

zU.product = function (x, y) {

    var xIsArray = Array.isArray(x);
    var yIsArray = Array.isArray(y);

    // 2 scalaires
    if (!xIsArray && !yIsArray) {
        return x * y;
    }

    // 1 scalaire et un tableau
    var nx, ny;
    if ((xIsArray && !yIsArray) || (!xIsArray && yIsArray)) {
        if (xIsArray) {
            nx = x;
            ny = y;
        }
        else {
            nx = y;
            ny = x;
        }

        if (nx.length === 0) {
            return [];
        }

        return nx.map(function (b) {
            return b * ny;
        });
    }


    // 2 tableaux
    // on additionne point à point en se basant sur la longueur du plus long et en bouclant sur les éléments du plus court en prenant un modulo.


    if (x.length >= y.length) {
        nX = x;
        nY = y;
    }
    else {
        nX = y;
        nY = x;
    }

    if ((x.length === 0) || (y.length === 0)) {
        return [];
    }

    var lY = nY.length;
    return nX.map(function (a, i) {
        return a * nY[i % lY];
    });
};


zygotine.U.sumProduct = function (x, y) {
    var p = zygotine.U.product(x, y);
    if (!Array.isArray(p)) {
        return p;
    }

    return zygotine.U.sum(p);
};

zU.divide = function (x, y) {

    var xIsArray = Array.isArray(x);
    var yIsArray = Array.isArray(y);
    var permuted;
    // 2 scalaires
    if (!xIsArray && !yIsArray) {
        return x / y;
    }

    // 1 scalaire et un tableau
    var nx, ny;
    var fn;
    if ((xIsArray && !yIsArray) || (!xIsArray && yIsArray)) {
        if (xIsArray) {
            nx = x;
            ny = y;
            permuted = false;
        }
        else {
            permuted = true;
            nx = y;
            ny = x;
        }

        if (nx.length === 0) {
            return [];
        }


        if (permuted) {
            fn = function (b) {
                return ny / b;
            };

        } else {
            fn = function (b) {
                return b / ny;
            };
        }

        return nx.map(fn);
    }

    // 2 tableaux
    // on soustrait point à point en se basant sur la longueur du long et en bouclant sur les éléments du plus court en prenant un modulo.
    if ((x.length === 0) || (y.length === 0)) {
        return [];
    }

    if (x.length >= y.length) {
        nX = x;
        nY = y;
        permuted = false;
    }
    else {
        nX = y;
        nY = x;
        permuted = true;
    }

    var lY = nY.length;

    fn = permuted ? function (a, i) { return nY[i % lY] / a; } : function (a, i) { return a / nY[i % lY]; };

    return nX.map(fn);
};

zU.minMax = function (x) {

    var isArray = Array.isArray(x);
    if (!isArray || (x.length === 1)) {
        let y = isArray ? x[0] : x;
        return [y, y];
    }

    if (x.length === 0) {
        return [Infinity, -Infinity];
    }

    // un tableau dont la longueur >= 2
    let fn = function (acc, b) { if (b < acc[0]) acc[0] = b; if (b > acc[1]) acc[1] = b; return acc; };
    return x.reduce(fn, [Infinity, -Infinity]);
};

zU.toString = function (object) {
    // ici object est un objet simple qui contient des nombres, des booléens et des chaines et  des tableaux de ces types.
    // c'est la responsabilité du programme appelant de satisfaire à ces conditions
    // 
    var keys = Object.keys(object).sort();
    var f = function (key) {
        var param = this[key];
        var txt = (Array.isArray(param)) ? zU.fmtArray(param) : param.toString();
        return key + ": " + txt;
    };

    return keys.map(function (key, index) {
        return f.call(object, key);
    }, this).join(", ");
};

zygotine.U.devConsoleHelper = {};

zygotine.U.devConsoleHelper.getIds = function () {

    $('[id]').each(function (index) {
        var o = $(this);
        console.log(index + " :id=" + o.attr('id') + " :tag=" + o.prop('tagName') + " :class=" + o.attr('class'));
    });
};

zygotine.U.permutation = function (n, perm, nans, ans) {
    var
        //double
        rT, mass, totalmass,
        //int
        i, j, k, n1, p = 1.0 / n;

    /* Record element identities */
    for (i = 0; i < n; i++) {
        perm[i] = i + 1;
    }
    var unif_rand = Math.random;
    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    //revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n - 1; i < nans; i++ , n1--) {
        rT = unif_rand();
        mass = 0;
        for (j = 0; j < n1; j++) {
            mass += p;
            if (rT <= mass)
                break;
        }

        ans[i] = perm[j];
        totalmass -= p;
        for (k = j; k < n1; k++) {
            //p[k] = p[k + 1];
            perm[k] = perm[k + 1];
        }
    }

    return { ans: ans, perm: perm };
};

zygotine.U.shuffle = function (n, k) {

    if ((typeof (k) === 'undefined') || (k < 1) || (k > n)) {
        k = n;
    }

    var x = new Array(n);
    var ik = new Array(k);
    for (let i = 0; i < n; i++) {
        x[i] = i;
    }

    for (let i = 0; i < k; i++) {
        let j = Math.floor(Math.random() * n);
        ik[i] = x[j] + 1;
        x[j] = x[--n];
    }

    return ik;
};
