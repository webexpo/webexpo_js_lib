/// <reference path="A.js" />
/// <reference path="M0.js" />
/// <reference path="M.js" />
/// <reference path="U.js" />

zygotine.X = {};
zygotine.X.lastModel = null;
zygotine.X.common = null

zygotine.X.genPseudoRand32Bit = function() {
  return (Math.random() * Math.pow(2, 31)) | 0
}

zygotine.X.setDefaultsForDistribution = function(loi) {
  zygotine.X.common.dataEntries.prngSeed.reset()
}

zygotine.X.ready = function() {
  $('.toggle-show').click(function(e) {
    $targ = $(this).find('svg[data-icon]')
    $targ.toggleClass('fa-plus-square')
    $targ.toggleClass('fa-minus-square')
    $span = $(this).find('span')
    let show = $span.attr('data-show')
    show = 1-parseInt(show)
    $span.attr('data-show', show)
    $span.text($.i18n(`show-examples-${show}`))
    $('#demoBtns').toggle()
  })
}
zygotine.X.setDataEntries = function() {
  let entries = zygotine.X.common.dataEntries
  entries.prngSeed = new zygotine.X.ValueBasedDataEntry("prngSeed", zygotine.X.genPseudoRand32Bit(), zygotine.X.i18n('algo-seed-expl', 'prngSeed'), true, 1, Math.pow(2,31)-1)
}

zygotine.X.getNumericalResult = function (
    logN,
    muChain,
    sigmaChain,
    oel,
    confidenceLevelForCredibileInterval,
    fracThreshold,
    percOfInterest
) {

    var reponse;
    var muC = muChain;
    var sigmaC = sigmaChain;
    var chaine = [];
    var quantile = new zygotine.S.Quantile([(100 - confidenceLevelForCredibileInterval) / 200, .5, 1 - (100 - confidenceLevelForCredibileInterval) / 200]);
    var t0 = performance.now(), t;

    if (logN) {

        reponse = {
            //gm pour J
            gMean: (function () {
                chaine = muC.map(b => Math.exp(b));
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "gMean", logN: logN, q: rep };
            })(),
            //gsd pour J
            gSd: (function () {
                chaine = sigmaChain.map(b => Math.exp(b));
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "gSd", logN: logN, q: rep };
            })(),

            //frac et frac.risk pour J
            exceedanceFraction: (function () {
                var logOel = Math.log(oel);
                var pNorm = zygotine.S.normal.cdf;
                var lowerTail = true;
                var logP = false;
                for (let i = 0; i < muC.length; i++) {
                    chaine[i] = 100 * (1 - pNorm((logOel - muC[i]) / sigmaC[i], 0, 1, lowerTail, logP));
                }

                //Overexposure risk based on exceedance fraction
                var risk = 100 * chaine.filter(function (b) { return b >= fracThreshold; }).length / chaine.length;
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "exceedanceFraction", logN: logN, q: rep, risk: risk };
            })(),

            //perc et perc.risk pour J
            percOfInterest: (function () { //percentile of interest

                var qNorm = zygotine.S.normal.icdf;
                var lowerTail = true;
                var logP = false;
                var icdfOfTargetQuantile = qNorm(percOfInterest / 100, 0, 1, lowerTail, logP);
                for (let i = 0; i < muC.length; i++) {
                    chaine[i] = Math.exp(muC[i] + (icdfOfTargetQuantile * sigmaC[i]));
                }

                var risk = 100 * chaine.filter(function (b) { return b > oel; }).length / chaine.length;
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "percOfInterest", logN: logN, q: rep, risk: risk };
            })(),

            // am et am.risk pour J
            aMean: (function () {  // arithmetic mean

                for (let i = 0; i < muC.length; i++) {
                    chaine[i] = Math.exp(muC[i] + (0.5 * sigmaC[i] * sigmaC[i]));
                }

                var risk = 100 * chaine.filter(function (b) { return b > oel; }).length / chaine.length;
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "aMean", logN: logN, q: rep, risk: risk };
            })()
        };

    } else {

        reponse = {
            //am pour J
            aMean: (function () {
                chaine = muC;
                var rep = quantile.compute(chaine);
                var risk = Math.round10(100.0 * (chaine.filter(function (b) { return b > oel; }).length / chaine.length), -1);
                return { src: "aMean", logN: logN, q: rep, risk: risk };
            })(),

            //asd pour J
            aSd: (function () {
                chaine = sigmaC;
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "aSd", logN: logN, q: rep };
            })(),

            //frac **** et frac.risk **** pour J
            exceedanceFraction: (function () {
                var pNorm = zygotine.S.normal.cdf;
                var lowerTail = true;
                var logP = false;
                for (let i = 0; i < muC.length; i++) {
                    chaine[i] = 100 * (1 - pNorm((oel - muC[i]) / sigmaC[i], 0, 1, lowerTail, logP));
                }

                var risk = 100 * chaine.filter(function (b) { return b >= fracThreshold; }).length / chaine.length;
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "exceedanceFraction", logN: logN, q: rep, risk: risk };
            })(),

            //perc et perc.risk pour J
            percOfInterest: (function () { //percentile of interest

                var qNorm = zygotine.S.normal.icdf;
                var lowerTail = true;
                var logP = false;
                var qNormOfPercOfInterest = qNorm(percOfInterest / 100, 0, 1, lowerTail, logP);
                for (let i = 0; i < muC.length; i++) {
                    chaine[i] = muC[i] + (qNormOfPercOfInterest * sigmaC[i]);
                }

                var risk = 100 * chaine.filter(function (b) { return b > oel; }).length / chaine.length;
                var rep = quantile.compute(chaine);
                chaine = [];
                return { src: "percOfInterest", logN: logN, q: rep, risk: risk };
            })()
        };
    }
    
    let calcAihaRiskBand = function(chainType = "p95") {
      var qNorm = zygotine.S.normal.icdf
      var lowerTail = true
      var logP = false
      var qNormOfPercOfInterest = qNorm(percOfInterest / 100, 0, 1, lowerTail, logP)
      for (let i = 0; i < muC.length; i++) {
        if ( chainType == "p95" ) {
          chaine[i] = muC[i] + (qNormOfPercOfInterest * sigmaC[i])
          if ( logN ) { chaine[i] = Math.exp(chaine[i]) }
        } else if ( chainType == "am" ) {
          chaine[i] = logN ? Math.exp(muC[i] + (0.5 * sigmaC[i] * sigmaC[i])) : muC[i];
        }
      }

      let overexpo = function(scaleFactor = 1) { return 100 * chaine.filter(function (b) { return b > scaleFactor*oel; }).length / chaine.length; }
      let riskVals = [
        (100 - overexpo(.01)).toFixed(1),
        (overexpo(.01) - overexpo(.1)).toFixed(1),
        (overexpo(.1) - overexpo(.5)).toFixed(1),
        (overexpo(.5) - overexpo()).toFixed(1),
        overexpo().toFixed(1)
      ]
      return riskVals
    }
    reponse.aihaBandP95 = (function () {
        let vals = calcAihaRiskBand("p95")
        return { src: "aihaBandP95", logN: logN, q: vals, risk: vals.join(' / ') };
    })()
    reponse.aihaBandAM = (function () {
        let vals = calcAihaRiskBand("am")
        return { src: "aihaBandPAM", logN: logN, q: vals, risk: vals.join(' / ') };
    })()
    
    t = performance.now() - t0;
    return reponse;
}; //zygotine.X.getNumericalResultFunctions


zygotine.X.ValueBasedDataEntry = function (elementId, initialValue, title, isInt, min, max) {
    this.elementId = elementId;
    this.initialValue = initialValue;
    this.element = $("#" + this.elementId);
    this.element.attr("title", title);
    this.min = (typeof min === 'number') ? min : (this.isInteger ? -500000 : Number.NEGATIVE_INFINITY);
    this.max = (typeof max === 'number') ? max : (this.isInteger ? 500000 : Number.POSITIVE_INFINITY);
    this.currentValue = '';
    this.element.change(this.getChangeFn());
    this.isInteger = (typeof isInt === 'boolean') ? isInt : false;
    this.validation = null;
};

zygotine.X.ValueBasedDataEntry.prototype = {

    type: 'value',

    reset: function () {
        this.element.val(this.initialValue).trigger('change');
    },

    getChangeFn() {
        var entry = this;
        return function () {
            entry.currentValue = $('#' + entry.elementId).val();
            entry.validate();
            if (entry.validation.ok) {
                entry.element.removeClass('invalid');
            } else {
                entry.element.addClass('invalid');
            }
        };
    },

    validate: function () {
        return this.isInteger ? this.__isInt() : this.__isFloat();
    },

    __isInt: function () {
        var x = this.currentValue.replace(/\s/g, '');
        var y = Number(x);
        var ok = Number.isInteger(y);
        var empty = (x === '');
        var withinBounds = (y >= this.min) && (y <= this.max);
        this.validation = { ok: ok && !empty && withinBounds, empty: empty, withinBounds: withinBounds, sign: Math.sign(y), val: y };
        return this.validation;
    },


    __isFloat: function () {
        var x = this.currentValue.replace(/\s/g, '');
        var y = Number(x);
        var ok = isFinite(y);
        var empty = (x === '');
        var withinBounds = (y >= this.min) && (y <= this.max);
        ok = ok && withinBounds;
        this.validation = { ok: ok && !empty && withinBounds, empty: empty, withinBounds: withinBounds, sign: Math.sign(y), val: y };
        return this.validation;
    }
};



zygotine.X.RadioButtonBasedDataEntry = function (name, initialObjectId) {
    this.initialElement = $('#' + initialObjectId);
    this.name = name;
    this.element = $('input[name=' + name + ']');
    this.currentValue = '';
};




zygotine.X.RadioButtonBasedDataEntry.prototype = {

    type: 'radio',

    reset: function () {
        this.initialElement.prop('checked', true).trigger('change');

    },

    getChangeFn: function () {
        var entry = this;
        var oStr = 'input[name=' + entry.name + ']:checked'
        return function () {
            entry.currentValue = $(oStr).val();
        }
    }
};

zygotine.X.checkboxBasedDataEntry = function (elementId, initialValue) {
    this.elementId = elementId;
    this.initialValue = initialValue;
    this.element = $("#" + this.elementId);
    this.currentValue = '';
    this.validation = { ok: true, val: this.currentValue };
    this.element.change(this.getChangeFn());
    //this.element.val(initialValue).trigger("change");
};

zygotine.X.checkboxBasedDataEntry.prototype = {

    reset: function () {
        this.element.prop('checked', this.initialValue).trigger('change');
    },

    getChangeFn() {
        var entry = this;
        return function () {
            entry.currentValue = $('#' + entry.elementId).prop('checked');
            entry.validation.val = entry.currentValue;
        };
    }
};

zygotine.X.r10 = function (n) {
    const exp10 = 3;
    var sign = Math.sign(n);
    if (sign !== 0) {
        if (sign === -1) {
            n = -n;
        }

        var tmp = Math.log10(n);
        var e = (tmp > exp10 ? exp10 : Math.round(tmp)) - exp10;
        var x = Math.round10(n, e);
        if (sign === -1) {
            x = -x;
        }

        return x.toString();
    } else {
        return '0.0';
    }
};

zygotine.X.display1NumericalResult = (function () {
    var r10 = zygotine.X.r10;

    var f = function (id, result) {
        var element;
        var targetElement;
        element = $("#" + id + "_NumRes").children(".VALUES");
        if (typeof result.singleValue === 'undefined') {
            targetElement = element.children(".B");
            targetElement.text(r10(result.q[1]));
            targetElement = element.children(".A");
            targetElement.text(r10(result.q[0]));
            targetElement = element.children(".C");
            targetElement.text(r10(result.q[2]));
            if (typeof result.risk !== 'undefined') {
                element = $("#" + id + "Risk_NumRes").children(".VALUES");
                if (typeof result.risk !== 'object') {
                    targetElement = element.children(".R");
                    targetElement.text(r10(result.risk));
                } else {
                    if (!Array.isArray(result.risk)) {
                        return;
                    } else {
                        for (let iR = 0; iR < result.risk.length; iR++) {
                            element = $("#" + id + "Risk" + iR + "_NumRes").children(".VALUES");
                            targetElement = element.children(".R");
                            targetElement.text(r10(result.risk[iR]));
                        }
                    }
                }
            }
        } else {
            targetElement = element.children(".V");
            targetElement.text(r10(result.singleValue));
        }
    };

    return f;
})();

zygotine.X.concatChains = function (result, flavor, chainSep = '\t', eol = '\r\n', decimalSeparator) {
    var replaceDecimalPoint = (typeof decimalSeparator === 'string');
    if (replaceDecimalPoint) {
        replaceDecimalPoint = decimalSeparator[0] !== '.';
    }

    var shrKut = result.chainByType[flavor];
    var headers = shrKut.map(x => x.label).join('\t');
    var lRes = shrKut[0].data.length + 1;
    var res = new Array(lRes);
    res[0] = headers;
    var nChains = shrKut.length;
    var oneLineData = new Array(nChains);

    for (iL = 0; iL < lRes - 1; iL++) {
        for (let iC = 0; iC < nChains; iC++) {
            oneLineData[iC] = shrKut[iC].data[iL];
        }

        res[iL + 1] = oneLineData.join(chainSep);
    }

    var rep = res.join(eol);
    return replaceDecimalPoint ? rep.replace(/\./g, decimalSeparator[0]) : rep;
}; 

zygotine.X.hideBePatient = function () {
    $('#bePatient').hide();
};

zygotine.X.showBePatient = function () {
    $('#bePatient').show();
};

zygotine.X.alert = function (message, title) {
    var jqobj = $('#dialog-message');
    jqobj.attr('title', title);
    jqobj.children('p').html(message);
    zygotine.X.a = $("#dialog-message").dialog({
        width: 500,
        modal: true,
        //open: function (event, ui) {
        //    $(".ui-dialog-titlebar-close").hide();
        //},
        buttons: {
            Ok: function () {
                $(this).dialog("close");
            }
        }
    });
}

zygotine.X.i18n = function(msgid, elemid) {
  (waitForTranslation = function() {
    if ( !$('body').data('trans-done') ) {
      setTimeout(waitForTranslation, 500)
    } else {
      $(`#${elemid}`).attr('title', $.i18n(msgid))
    }
  })()
}

zygotine.X.RiskbandFloatSetDataEntry = function (className, title, min, max) {
  this.className = className
  this.currentValue = null
  this.elements = $(`.${this.className}`)
  this.elements.attr("title", title)
  this.min = min == undefined ? Number.NEGATIVE_INFINITY : min
  this.max = max == undefined ? Number.POSITIVE_INFINITY : max
  this.validation = null
  this.reset()
}

zygotine.X.RiskbandFloatSetDataEntry.prototype = {
  type: 'riskband-float-set',
  reset: function () {
    this.elements = $(`.${this.className}`)
    if ( this.elements.length > 0 ) {
      this.elements.change(this.getChangeFn())
      this.elements.first().trigger('change')
    }
  },
  getChangeFn() {
    var dis = this
    return function () {
      dis.currentValue = dis.elements.map(function() { return $(this).val() }).get().map(parseFloat)
      dis.validate()
      dis.elements.removeClass('invalid')
      if (!dis.validation.ok) {
        dis.elements.filter((i, _) => !dis.validation.elems[i]).addClass('invalid')
      }
    }
  },
  validate: function () {
    let validElems = this.currentValue.map(val => !isNaN(val) && val >= this.min && val <= this.max)
    let ok = validElems.filter(v => v).length == validElems.length
    this.validation = { ok, elems: validElems }
  }
}