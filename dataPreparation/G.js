
// très très fortement inspiré de https://developer.mozilla.org/fr/docs/Web/JavaScript/Reference/Objets_globaux/Array/filter
if (!Array.prototype.which) {
    Array.prototype.which = function (fun /*, thisArg */) {
        'use strict';

        if (this === undefined || this === null) {
            throw new TypeError();
        }

        var t = Object(this);
        var len = t.length >>> 0;

        // NOTE : fix to avoid very long loop on negative length value
        if (len > t.length || typeof fun !== 'function') {
            throw new TypeError();
        }

        var res = [];
        var thisArg = arguments.length >= 2 ? arguments[1] : void 0;
        for (var i = 0; i < len; i++) {
            if (i in t) {
                var val = t[i];
                if (fun.call(thisArg, val, i, t))
                    res.push(i);
            }
        }

        return res;
    };
}

if (!Math.trunc) {
    Math.trunc = function (v) {
        v = +v;
        if (!isFinite(v)) return v;

        return v - v % 1 || (v < 0 ? -0 : v === 0 ? v : 0);

        // returns:
        //  0        ->  0
        // -0        -> -0
        //  0.2      ->  0
        // -0.2      -> -0
        //  0.7      ->  0
        // -0.7      -> -0
        //  Infinity ->  Infinity
        // -Infinity -> -Infinity
        //  NaN      ->  NaN
        //  null     ->  0
    };
}

if (!Array.prototype.fill) {
    Object.defineProperty(Array.prototype, 'fill', {
        value: function (value) {

            // Steps 1-2.
            if (this === null) {
                throw new TypeError('this is null or not defined');
            }

            var O = Object(this);

            // Steps 3-5.
            var len = O.length >>> 0;

            // Steps 6-7.
            var start = arguments[1];
            var relativeStart = start >> 0;

            // Step 8.
            var k = relativeStart < 0 ?
              Math.max(len + relativeStart, 0) :
              Math.min(relativeStart, len);

            // Steps 9-10.
            var end = arguments[2];
            var relativeEnd = end === undefined ?
                len : end >> 0;

            // Step 11.
            var final = relativeEnd < 0 ?
              Math.max(len + relativeEnd, 0) :
              Math.min(relativeEnd, len);

            // Step 12.
            while (k < final) {
                O[k] = value;
                k++;
            }

            // Step 13.
            return O;
        }
    });
}

(function () {
    /**
     * Ajustement décimal d'un nombre
     *
     * @param {String}  type : Le type d'ajustement souhaité.
     * @param {Number}  value : le nombre à traité The number.
     * @param {Integer} exp  : l'exposant (le logarithme en base 10 de l'ajustement).
     * @returns {Number} la valeur ajustée.
     */
    function decimalAdjust(type, value, exp) {
        // Si la valeur de exp n'est pas définie ou vaut zéro...
        if (typeof exp === 'undefined' || +exp === 0) {
            return Math[type](value);
        }
        value = +value;
        exp = +exp;
        // Si la valeur n'est pas un nombre
        // ou si exp n'est pas un entier...
        if (value === null || isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0)) {
            return NaN;
        }
        // Si la valeur est négative
        if (value < 0) {
            return -decimalAdjust(type, -value, exp);
        }
        // Décalage
        value = value.toString().split('e');
        value = Math[type](+(value[0] + 'e' + (value[1] ? +value[1] - exp : -exp)));
        // Décalage inversé
        value = value.toString().split('e');
        return +(value[0] + 'e' + (value[1] ? +value[1] + exp : exp));
    }

    // Arrondi décimal
    if (!Math.round10) {
        Math.round10 = function (value, exp) {
            return decimalAdjust('round', value, exp);
        };
    }
    // Arrondi décimal inférieur
    if (!Math.floor10) {
        Math.floor10 = function (value, exp) {
            return decimalAdjust('floor', value, exp);
        };
    }
    // Arrondi décimal supérieur
    if (!Math.ceil10) {
        Math.ceil10 = function (value, exp) {
            return decimalAdjust('ceil', value, exp);
        };
    }
})();


