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
  entries.prngSeed = new zygotine.X.ValueBasedDataEntry("prngSeed", zygotine.X.genPseudoRand32Bit(), zygotine.X.i18n('algo-seed-expl', 'prngSeed'), true, 1, Math.pow(2,31)-1);
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

zygotine.X.concatChains = function (result, flavor, decimalSeparator) {
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

        res[iL + 1] = oneLineData.join('\t');
    }

    var rep = res.join('\r\n');
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

function initLocale()
{
  $.i18n.debug = true;
        
  var url = new URL(window.location.href);
  var lang = url.searchParams.get('lang');
  lang = lang == null ? 'fr' : lang.toLowerCase();
  $.i18n.locale = lang;
}

function translatePage()
{
  var i18n = $.i18n({locale: $.i18n.locale});
  $('body').attr('data-lang', $.i18n.locale);
  $('body').data('trans-done', 0)
  $.i18n().load( `i18n/trans-${i18n.locale}.json`, i18n.locale ).done(function(x) {
    $('html').i18n();
      $('body').data('trans-done', 1)
  });
  $('.lang-switcher .lang').each(function() {
    var this_locale = $(this)[0].classList[1];
    if ( $.i18n.locale != this_locale ) {
      $(this).attr('href', window.location.origin + window.location.pathname + "?lang=" + this_locale);
    } else {
      $(this).replaceWith("<span class='lang selected'>" + $(this).text() + "</span>");
    }
  });
}

function downloadTraceplot(mcmcParam) {
  let modelType = typeof(zygotine.SEG) !== "undefined" ? "SEG" : "BW"
  let burninChain = zygotine[modelType].lastModel.result.chains[`${mcmcParam.name}Burnin`].data
  let mainChain = zygotine[modelType].lastModel.result.chains[`${mcmcParam.name}Sample`].data
  if ( mainChain.length > 0 ) {
    let plotElem = document.createElement('div')
    let data = [
      {
        x: [...Array(burninChain.length).keys()].map(x => x+1),
        y: burninChain,
        line: {
          color: 'red'
        },
        name: 'Burnin'
      },
      {
        x: [...Array(mainChain.length).keys()].map(x => burninChain.length+x+1),
        y: mainChain,
        line: {
          color: "#1f77b4"
        },
        name: $.i18n('Posterior sample')
      }
    ]
    
    var layout = {
      title: {
        text: `${$.i18n('traceplot-title')} <i>${mcmcParam.symbol}</i>`,
        font: {
          family: 'Arial Black', 
          color: 'black',
          size: 16
        },
        y: 0.98,
        yanchor: 'top'
      },
      xaxis: {
        title: $.i18n('Iteration'),
        range: [0, burninChain.length + mainChain.length + 500],
        gridcolor: 'darkgrey',
        dtick: 1000
      },
      yaxis: {
        title: mcmcParam.symbol,
        showgrid: false
      },
      margin: {
        t: 0
      },
      showlegend: true,
      legend: {
        orientation: 'h',
        borderwidth: 1
      }
    }
    
    Plotly.newPlot(
      plotElem,
      data,
      layout
    ).then(function(gd) {
      let now = new Date()
      let dateTimeFormat = new Intl.DateTimeFormat('en', { year: 'numeric', month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit', second: '2-digit', hour12: false }) 
      let [{ value: month },,{ value: day },,{ value: year },,{ value: hour },,{ value: minute },,{ value: sec },,] = dateTimeFormat .formatToParts(now) 
      let plotFilename = `${$.i18n('traceplot-filename')}_${mcmcParam.symbol.replace(/<[^>]+>(?=.)/g, '_').replace(/<[^>]+>$/, '')}_${year}${month}${day}_${hour}${minute}${sec}`
      Plotly.downloadImage(gd, {
        format: 'png',
        height: 600,
        width: 1000,
        filename: plotFilename
      })
      /* gd.remove() */
    })
  }
}