zygotine.M = {};

/*****  ErrorLogging  **********/

zygotine.M.ErrorLogging = function () {
    this.messages = new zygotine.M.Messages(this.className); // no message, either warning or error
    this.hasError = false;
};

zygotine.M.ErrorLogging.prototype = {

    /**
    * used to add an error Message to the property messages. The method returns nothing
    * @method 
    * @param {string} msgText - string representing the message
    * @returns {undefined}
    */
    addError: function (msgText) {
        this.messages.addError(msgText);
        this.hasError = true;
    },

    /**
     * used to add a warning (Message) to the property messages. The method returns nothing
     * @method
     * @param {string} msgText - string representing the message
     * @returns {undefined}
     */
    addWarning: function (msgText) {
        this.messages.addWarning(msgText);
    },

    addInfo: function (msgText) {
        this.messages.addInfo(msgText);
    },

    /**
 * add each message contained in msgs using method addMessage 
 * @method
 * @param {zygotine.M.Messages} msgs - the Messages that contains the individual messages to be merged
 * @param {boolean} takeCopy - doit être un véritable booléen true ou false!
 *
 * @returns {undefined}
 */
    mergeMessagesFrom: function (msgs) {
        var i;

        var /** @type {zygotine.M.Message} */ msg;
        for (i = 0; i < msgs.msgList.length; i++) {
            msg = msgs.msgList[i];
            this.messages.addMessage(new zygotine.M.Message(msg.msgText, msg.level, msg.source));
            if (msg.level === 'error') {
                this.hasError = true;
            }
        }
    }
};

/*******************************************/

/**
* Represents a message (error message or warning)
* @constructor
* @param {string} msg - the message
* @param {string} level - "info" | "warning" | "error"
* @param {string} src - the message src
*
* @property {string} level - "info" | "warning" | "error"
* @property {string} msgText - the message
* @property {string} src - object class name from which originates the message
*/
zygotine.M.Message = function (msg, level, src) {
    // msg and src are string, error is a boolean( true for error, false for warning)
    this.msgText = msg;
    this.level = level;
    this.source = src;
};

zygotine.M.Message.prototype = {
    toString: function () {
        return zU.fmt("{0}: {1}", this.level, this.msgText);
    }
};

/*****  Messages  **********/

/**
* Represents a list of messages
* @constructor
* @param {string} errSrc - default error source
*
* @property {boolean} hasError - by default it is false
* @property {zygotine.M.Message[]} msgList - an array of Message 
* @property {string} errorSrc; 
*/
zygotine.M.Messages = function (errSrc) {
    this.hasError = false;
    this.msgList = [];
    this.errorSrc = errSrc;
};

zygotine.M.Messages.prototype = {


    /**
     * add a warning Message. The method returns nothing
     * @method
     * @param {zygotine.O.Message} msgObject - a message object
     * @returns {undefined}
     */
    addMessage: function (msgObject) {
        this.msgList[this.msgList.length] = msgObject;
        if (msgObject.level === "error") {
            this.hasError = true;
        }
    },

    /**
    * create an error Message and add it. The method returns nothing.
    * @method
    * @param {string} msgText - string representing the message
    * @returns {undefined}
    */
    addError: function (msgText) {
        this.addMessage(new zygotine.M.Message(msgText, "error", this.errorSrc));
    },

    /**
    * create a warning Message and add it. The method returns nothing.
    * @method
    * @param {string} msgText - string representing the message
    * @returns {undefined}
    */
    addWarning: function (msgText) {
        this.addMessage(new zygotine.M.Message(msgText, "warning", this.errorSrc));
    },

    /**
    * create an informative Message and add it. The method returns nothing.
    * @method
    * @param {string} msgText - string representing the message
    * @returns {undefined}
    */
    addInfo: function (msgText) {
        this.addMessage(new zygotine.M.Message(msgText, "info", this.errorSrc));
    },

    toString: function () {
        let rep = [];
        rep.push(zU.fmt("hasError: {0}", this.hasError));
        rep.push(this.msgList.map(z => z.toString()).join("\n"));
        rep = rep.join("\n");
        return rep;
    },

    getErrors: function () {
        var rep = '';
        if (!this.hasError) {
            return rep;
        }

        rep = this.msgList.filter(z => z.level === 'error').map(m => m.msgText + "<br/>");
        return rep;
    }

};

/*******************************************/



/*****  Measure  **********/

/**
* Represents a Measure.
* @constructor
* @property {number} a - double
* @property {number} b - double
* @property {number} ordinal - int
* @property {string} type - measure type
* @property {string} workerId - the id of the worker, if available, to which the measure is associated
*/
zygotine.M.Measure = function () {
    // s a string
    this.a = Number.NaN;
    this.b = Number.NaN;
    this.ordinal = 0;
    this.type = undefined;
    this.workerId = "";
    this.intervalSeparator = '-'
};

/**
 * gives a representation of a Measure.
 * @method toString 
 * @returns {string} - a string that represents the Measure
 */
zygotine.M.Measure.prototype.toString = function (fieldSep) {

    if (typeof fieldSep === "undefined") {
        fieldSep = ';';
    }

    var rep = [];
    switch (this.type) {
        case 'uncensored':
            rep[0] = this.a;
            break;

        case 'greaterThan':
        case 'lessThan':
            rep[0] = (this.type === 'lessThan') ? "<" : ">";
            rep[1] = this.a;
            break;

        case 'interval':
            rep[0] = "[";
            rep[1] = this.a;
            rep[2] = this.intervalSeparator;
            rep[3] = this.b;
            rep[4] = "]";
            break;
    }

    var ret = rep.join("");
    if (this.workerId !== "") {
        ret += fieldSep + this.workerId;
    }

    return ret;
};

zygotine.M.Measure.prototype.clone = function () {
    var m = new zygotine.M.Measure();
    m.a = this.a;
    m.b = this.b;
    m.ordinal = this.ordinal;
    m.type = this.type;
    m.workerId = this.workerId;
    return m;
};

/*******************************************/


/*****  ExtendedMeasure  **********/

/**
* Represents a Measure.
* @constructor
* @param {zygotine.M.Measure} m - all properties of m will be duplicated
*
* @property {number} a - double
* @property {number} b - double
* @property {number} currentValue - a generated value
* @property {number} ordinal - int
* @property {string} type - measure type
* @property {string} workerId - the id of the worker, if available, to which the measure is associated

*/
zygotine.M.ExtendedMeasure = function (m) {
    // s a string
    this.a = m.a;
    this.b = m.b;
    this.ordinal = m.ordinal;
    this.type = m.type;
    this.workerId = m.workerId;
    this.intervalSeparator = '-';
    this.generatedValue = this.getInitialValue();
};

zygotine.M.ExtendedMeasure.prototype.getInitialValue = function () {
    return this.type === 'interval' ? (this.a + this.b) / 2.0 : this.a;
};

/**
 * gives a representation of a ExtendedMeasure.
 * @method toString 
 * @returns {string} - a string that represents the Measure
 */
zygotine.M.ExtendedMeasure.prototype.toString = function () {

    var rep = [];
    switch (this.type) {
        case 'uncensored':
            rep[0] = this.a;
            break;

        case 'greaterThan':
        case 'lessThan':
            rep[0] = (this.type === 'lessThan') ? "<" : ">";
            rep[1] = this.a;
            break;

        case 'interval':
            rep[0] = "[";
            rep[1] = this.a;
            rep[2] = this.intervalSeparator;
            rep[3] = this.b;
            rep[4] = "]";
            break;
    }
    rep.push(" => ");
    rep.push(this.generatedValue);
    rep.push("(");
    rep.push(this.workerId);
    rep.push(")");
    return rep.join("; ");
};


/*******************************************/


/*****  MeasureList  **********/

/**
* @description Represents a list of measures.
* @constructor
* @param {string} s - une chaine représantant les mesures sans 'erreur de mesure'
* @param {string} withWorkerInfo - un booléen - par défaut 'false'
* @property {boolean} anyCensored
* @property {boolean} hasError
* @property {boolean} hasME
* @property {zygotine.M.Messages} messages
* @property {object} measureByType
* @property {object} measureCountByType
* @property {zygotine.M.Measure[]} measureList
* @property {string} METext
* @property {number} n
* @property {RegExp} rgxME
*/
zygotine.M.MeasureList = function (s, withWorkerInfo) {

    zygotine.M.ErrorLogging.call(this);
    if (zygotine.isUndefined(withWorkerInfo)) {
        withWorkerInfo = false;
    }

    this.withWorkerInfo = withWorkerInfo;
    this.anyCensored = false;
    this.hasME = false;
    this.measureByType = { 'interval': [], 'lessThan': [], 'greaterThan': [], 'uncensored': [] };
    this.measureCountByType = { 'interval': 0, 'lessThan': 0, 'greaterThan': 0, 'uncensored': 0 };
    this.measureByWorker = {};
    this.measureList = [];
    this.METext = "";
    this.n = 0; // aucune mesure reconnue.
    this.nWAssigned = 0; // aucune mesure n'a un travailleur associé
    this.min = Infinity;
    this.max = -Infinity;
    this.intervalSeparator = '-'

    var fmt = zygotine.U.fmt;

    if ( this.nullList == undefined ) {
        this.nullList = false;
    }

    if (!s) {
        if (  !this.nullList ) {
            this.addError("Parameter \"s\" is required.");
        }
        return;
    }

    s = s.replace(/ /g, "");
    var records = this.splitRecords(s);
    if (records.length === 0) {
        this.addError(fmt("No measure found in string '{0}'!", s));
        return;
    }

    var i;
    /*@type {zygotine.M.Measure} */
    var m;

    this.hasME = this.rgxME.test(records[0]);
    if (this.hasME) {
        this.METext = records[0];
        this.addWarming("Measurement error not supported, it will be ignored. First record skipped!");
    }

    for (i = this.hasME ? 1 : 0; i < records.length; i++) {
        //parseMeasure retourne une mesure
        let fields = this.splitFields(records[i]);
        if (fields.length === 0) {
            // on passe à l'enregistrement suivant
            continue;
        }

        m = this.parseMeasure(fields[0]);

        if (m.invalid) {
            this.addError(fmt("Invalid measure '{0}' : record {1}.", records[i], i + 1));
            //break //*** on ne s'occupe pas de la suite en sortant de la boucle for.
            continue;
        } else if ((fields.length > 1) && (fields[1] !== "")) {
            m.workerId = fields[1];
        }

        this.addMeasure(m);
    } // for

    //On vérifie n
    if (!this.hasError) {
        this.n = this.measureCountByType.uncensored + this.measureCountByType.interval + this.measureCountByType.lessThan + this.measureCountByType.greaterThan;

        if ((this.n === 0)) {
            this.addError("No measure recognized.");
        } else {
            var aOrB;
            for (let i = 0; i < this.measureList.length; i++) {
                m = this.measureList[i];
                if (m.a < this.min) {
                    this.min = m.a;
                }

                let aOrB = (m.type === "interval") ? m.b : m.a;
                if (aOrB > this.max) {
                    this.max = aOrB;
                }
            }

            if (this.min <= 0) {
                this.addError("All values must be stricly positive!");
            }

            if (this.withWorkerInfo && (this.nWAssigned !== this.n)) {
                this.addError('All measures must be assigned to a worker.');
            }

        }
    }
}; // zygotine.M.MeasureList

zygotine.M.MeasureList.prototype = Object.create(zygotine.M.ErrorLogging.prototype);

zygotine.M.MeasureList.prototype.className = "MeasureList";

zygotine.M.MeasureList.prototype.rgxME = new RegExp("(sd|cv)\(.*\)", 'i');

// s is a string
// return an array of strings
zygotine.M.MeasureList.prototype.splitRecords = function (s) {
    var tmp = s.split(/[|\r\n]+/);
    tmp = tmp.filter(e => e !== "");
    return tmp;
};

/**
* @function - split string based on ';', '\t'
* @param {string} s - string to be splitted into fields.
* @returns {string[]} - the input string splitted
*/
zygotine.M.MeasureList.prototype.splitFields = function (s) {
    return s.split(/[;\t]+/);
};

/**
* Add a Measure. No test is made to check if for the measure to be added,  m.invalid is true.
* @method
* @param {Measure} m - a measure
* @returns {undefined}
*/
zygotine.M.MeasureList.prototype.addMeasure = function (/** @type  {zygotine.M.Measure}*/ m) {
    m.ordinal = this.measureList.length;
    this.measureList[this.measureList.length] = m;
    if (m.type !== 'uncensored') {
        this.anyCensored = true;
    }
    //Update collections
    this.measureCountByType[m.type] = this.measureCountByType[m.type] + 1;
    /** @type {number[]} */
    this.n++;
    if (m.workerId !== "") {
        this.nWAssigned++;
        this.measureByWorker[m.workerId] = this.measureByWorker[m.workerId] || [];
        this.measureByWorker[m.workerId].push(m);
    }

    this.measureByType[m.type].push(m);
};

/**
* parse a string and try to convert it to a measure. If the parsing fail, the invalid property of the returned measure is  set to true.
* @method
* @param {string} str
* @returns {zygotine.M.Measure} - a measure. 
*/
zygotine.M.MeasureList.prototype.parseMeasure = function (str) {
    // str is a string
    var pm = new zygotine.M.Measure();
    switch (str[0]) {
        case '[':
            pm.invalid = str.slice(-1) !== ']';
            if (!pm.invalid) {
                //tout ce qu'il y a entre les 2 crochets
                str = str.substring(1, str.length - 1);
                var numbers = str.split(/[,-]/);
                if (numbers.length === 2) {
                    pm.type = "interval";
                    pm.a = Number(numbers[0]);
                    pm.invalid = !isFinite(pm.a);
                    if (!pm.invalid) {
                        pm.b = Number(numbers[1]);
                        pm.invalid = !isFinite(pm.b);
                        if (!pm.invalid) {
                            pm.invalid = pm.a >= pm.b;
                        }
                    }
                }
            }

            break;

        case '<':
            pm.type = "lessThan";
            pm.a = Number(str.substring(1, str.length));
            pm.invalid = !isFinite(pm.a);
            break;

        case '>':
            pm.type = "greaterThan";
            pm.a = Number(str.substring(1, str.length));
            pm.invalid = !isFinite(pm.a);
            break;

        default:
            pm.type = "uncensored";
            pm.a = Number(str);
            pm.invalid = !isFinite(pm.a);
            break;
    }

    return pm;
};

zygotine.M.MeasureList.prototype.transform = function (xFunction) {

    var ml = new zygotine.M.MeasureList(this.toString());
    for (let iM = 0; iM < ml.measureList.length; iM++) {
        /** type {zygotine.M.Measure} */
        m = ml.measureList[iM];
        m.a = xFunction(m.a);
        if (!isNaN(m.b)) {
            m.b = xFunction(m.b);
        }
    }
    // la chaine ml n'est pas cohérente; entre autres, le min et le max n'ont pas été ajustés.
    ml = new zygotine.M.MeasureList(ml.toString());
    return ml;
};

/**
* 
* @method toString - gives a representation of the MeasureList.
* @returns {string}
*/
zygotine.M.MeasureList.prototype.toString = function (recSep, fieldSep) {
    /** @type {zygotine.M.MeasureList} */
    if (zygotine.isUndefined(recSep)) {
        recSep = '|';
    }

    if (zygotine.isUndefined(fieldSep)) {
        fieldSep = ';';
    }

    var self = this;
    if (self.messages.invalid) {
        return "Invalid measure list!";
    }

    /** @type {Measure[]} */
    var x = self.measureList;
    var a = x.map(function (x) { return x.toString(fieldSep); });
    var ret = a.join(recSep);
    if (self.METext !== '') {
        ret = self.METext + " | " + ret;
    }

    return ret;
};



zygotine.M.MeasureList.divideByOel = function (measureList, oel) {
    //permet de créer des données compatibles avec ce qu'attend les modèles de McGill.
    var ml = new zygotine.M.MeasureList(measureList.toString());
    if (oel !== 1.0) {
        for (let iM = 0; iM < ml.measureList.length; iM++) {
            /** type {zygotine.M.Measure} */
            m = ml.measureList[iM];
            m.a = m.a / oel;
            if (!isNaN(m.b)) {
                m.b = m.b / oel;
            }
        }
    }
    // la chaine ml n'est pas cohérente; entre autres, le min et le max n'ont pas été ajustés.
    ml = new zygotine.M.MeasureList(ml.toString());
    return ml;
};

/*******************************************/

/*****  RawMeasures  **********/

/**
* Represents the measures in a simplified manner
* @constructor
*
* @property {number[]} gt
* @property {number[]} intervalGt
* @property {number[]} intervalLt
* @property {number[]} lt
* @property {number} uncensoredSum
* @property {number} uncensoredSum2
* @property {number[]} y
*/
zygotine.M.RawMeasures =
    function () {
        this.y = [];
        this.uncensoredSum = NaN;
        this.uncensoredSum2 = NaN;
        this.lt = [];
        this.gt = [];
        this.intervalGt = [];
        this.intervalLt = [];
        this.uncensoredSum = NaN;
        this.uncensoredSum2 = NaN;
    };

/*******************************************/

/*****  Range  **********/

/**
* un couple de nombres définissant un interval
* @param {number} l - la valeur à donner à la borne inférieure
* @param {number} u - la borne supérieure
* @constructor
* @property {number} lower - la borne inférieure
* @property {number} upper - la borne supérieure
*/
zygotine.M.Range = function (l, u) {
    this.lower = l;
    this.upper = u;
};

/*******************************************/



/***** DataSummary  ***********/

/**
* Represents the data provided by the user
* @constructor
* @param {zygotine.M.MeasureList} ml - a MeasureList
* @param {boolean} logNormalDistrn - either true, log normal or false, normal.
*
* @property {boolean} anyCensored - true if there is at least one censored measure
* @property {number[]} gt - number array for the "greater than" measures (>)
* @property {boolean} hasError - at least one error messages is an error.
* @property {zygotine.M.Range} i - un couple de nombres composé de lower et upper définissant les bornes d'un interval 
* @property {zygotine.M.RawMeasures} logTakenData - the log taken data as opposed to the origninal data.
* @property {boolean} logNormalDistrn - either true, log normal or false, normal. Same as the parameter of the same name
* @property {number[]} lt - an array of numbers for the "less than" measures (<)
* @property {zygotine.M.Measure[]} measureList - array of measures, these are duplicata of the measureList fiel of parameter ml.
* @property {zygotine.M.Messages} messages -  a container for warning or error messages
* @property {number} n - the number of measures
* @property {zygotine.M.RawMeasures} originalData - the original data as opposed to log taken data.
*/
zygotine.M.DataSummary = function (ml, logNormalDistrn) {

    zygotine.M.ErrorLogging.call(this);
    this.anyCensored = ml.anyCensored;
    this.logNormalDistrn = logNormalDistrn;
    this.n = ml.measureList.length; //nombre de mesures

    var f = logNormalDistrn ?
        function (m) {
            /** @type {zygotine.M.Measure} */
            var mc = m.clone();
            mc.a = Math.log(m.a);
            if (mc.type === 'interval') {
                mc.b = Math.log(m.b);
            }

            return mc;
        } :
        function (m) {
            /** @type {zygotine.M.Measure} */
            return m.clone();
        };

    this.measureList = ml.measureList.map(function (m) {
        /** @type {zygotine.M.Measure} */
        return f(m);
    });

    //DATA
    this.y = ml.measureByType["uncensored"].map(function (mes) { return mes.a; });
    this.gt = ml.measureByType["greaterThan"].map(function (mes) { return mes.a; });
    this.lt = ml.measureByType["lessThan"].map(function (mes) { return mes.a; });
    this.intervalGt = ml.measureByType["interval"].map(function (mes) { return mes.a; });
    this.intervalLt = ml.measureByType["interval"].map(function (mes) { return mes.b; });
    if (logNormalDistrn) {
        let f = function (a) {
            let l = a.length;
            for (let i = 0; i < l; i++) {
                a[i] = Math.log(a[i]);
            }
        };

        f(this.y);
        f(this.gt);
        f(this.lt);
        f(this.intervalGt);
        f(this.intervalLt);
    }

    var sumFun = zU.sum;
    var sumSqrFun = zU.sumSqr;
    this.i = new zygotine.M.Range(this.intervalGt, this.intervalLt);
    this.uncensoredSum = sumFun(this.y); //.reduce(function (accu, currentValue) { return accu + currentValue.a; }, 0.0);
    this.uncensoredSum2 = sumSqrFun(this.y); // .reduce(function (accu, currentValue) { return accu + (currentValue.a * currentValue.a); }, 0.0);
};


zygotine.M.DataSummary.prototype = Object.create(zygotine.M.ErrorLogging.prototype);
zygotine.M.DataSummary.prototype.className = "DataSummary";

/*******************************************/

/***** PastDataSummary *************/

/**
* Represents past data (from the litterature ...)
* @constructor
* @param {string} def chaine d'entrée
*
* @property {boolean} hasError
* @property {boolean} defined
* @property {number} logSigmaMu
* @property {number} logSigmaSd
* @property {number} mean
* @property {zygotine.M.Messages} messages
* @property {number} n
* @property {number} ns2
* @property {number} sd
* @property {number} shape
*/
zygotine.M.PastDataSummary =
    function (/** string */ def, isLogNDist = false, oel ) {
        zygotine.M.ErrorLogging.call(this);
        this.defined = false; // a priori l'objet est bidon
        if (typeof def === "undefined") {
            return this;
        }

        //valeurs fournies à l'entrée
        this.sd = this.mean = this.n = NaN;
        //valeurs calculés
        this.ns2 = NaN;

        def = def.replace(/ /g, "");
        var parameters = def.split(/[;\t]+/);
        if (parameters.length < 3) {
            this.addError("Invalid data. 3 parameters are expected: 'mean, sd and n'.");
        }
        else {
            this.mean = Number(parameters[0]);
            let standardizeMean = false
            if (!isFinite(this.mean)) {
                this.addError("Invalid data. Mean can not be parsed to a finite float.");
            }
            /* STANDARDISATION */
            if (zygotine.SEG != undefined) {
                standardizeMean = zygotine.SEG.dataEntries.dstrn.currentValue == "logN"
                oel = parseFloat(zygotine.SEG.dataEntries.oel.currentValue)
            } else
            if ( isLogNDist && oel != undefined ) {
                standardizeMean = true
            }

            if ( standardizeMean ) {
              this.mean -= Math.log(oel)
            }

            this.sd = Number(parameters[1]);
            if (!isFinite(this.sd) || this.sd <= 0) {
                this.addError("Invalid data. sd can not be parsed to a finite positive float.");
            }

            this.n = Number(parameters[2]);
            if (!isFinite(this.n) || !Number.isInteger(this.n)) {
                this.addError("Invalid data. n can not be parsed to an integer.");
            } else if (this.n < 2) {
                this.addError("Invalid data. n must be greater or equal to 2.");
            }

            if (!this.hasError) {
                this.defined = true;
                this.ns2 = (this.n - 1) * (this.sd * this.sd);
                this.sum = this.n * this.mean
            }
        }
    };

zygotine.M.PastDataSummary.prototype = Object.create(zygotine.M.ErrorLogging.prototype);

zygotine.M.PastDataSummary.prototype.className = "PastDataSummary";
zygotine.M.PastDataSummary.prototype.level = 0.95;

zygotine.M.PastDataSummary.prototype.toString = function () {

    return zygotine.U.fmt(
        "mean={0}, sd={1}, n={2}, ns2={3}, defined={4}",
        this.mean,
        this.sd,
        this.n,
        this.ns2,
        this.defined);
};

zygotine.M.PastDataSummary.dummyPDS = new zygotine.M.PastDataSummary();

/*******************************************/