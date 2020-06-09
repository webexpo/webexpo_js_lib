/* eslint 
    valid-jsdoc: 0,
    no-extra-parens : 0
*/
/// <reference path="A.js" /> 
/// <reference path="U.js" /> 

zygotine.Num = {};

var zNum = zygotine.Num;

/**
* @constructor
*
* @property {number[]} alist
*
*/
zNum.NumericIntegration = function (f, lowerBound, upperBound) {
    this.a = lowerBound;
    this.b = upperBound;
    this.fn = f;
    // results
    this.result = 0.0;
    this.neval = 0;
    this.ier = 6;
    this.last = 0;
    this.abserr = 0.0;
    // listes de travail
    //this.alist = []; // double
    //this.blist = []; // double
    //this.rlist = []; // double
    //this.elist = []; // double
    //this.iord = []; // int

};

zNum.NumericIntegration.convertFn4Integration = function (fn) {
    var Fn = function (xArray) {
        return xArray.map(function (a) { return fn(a); });
    };
    return Fn;
};

zNum.NumericIntegration.stopi = 8888;

zNum.NumericIntegration.stopd = 8888.7;

zNum.NumericIntegration.setArray = function (arr, _stop) {
    if (typeof _stop === 'undefined') {
        _stop = zNum.NumericIntegration.stopd;
    }

    for (let i = 0; i < arr.length; i++) {
        arr[i] = _stop;
    }
};

zNum.NumericIntegration.displayArray = function (name, arr, n) {
    var fmt = zygotine.U.fmt;
    var tmp = [];
    var last = NaN;

    tmp.push(fmt("{0}=", name));
    const stopd = zNum.NumericIntegration.stopd;
    const stopi = zNum.NumericIntegration.stopi;

    tmp.push(((arr[0] === stopd) || (arr[0] === stopi)) ? "Ok" : fmt("*{0}*", arr[0]));
    tmp.push(", ");

    for (let i = 1; i < n; i++) {
        if ((arr[i] === stopd) || (arr[i] === stopi)) {
            tmp.push(fmt("{0}({1})!", last, i - 1));
            break;
        }
        else {
            last = arr[i];
        }
    }

    return tmp.join("");
};

zNum.NumericIntegration.Workspace = function (limit) {

    this.alist = new Array(limit + 1);
    this.alist.fill(zNum.NumericIntegration.stopd);

    this.blist = new Array(limit + 1);
    this.blist.fill(zNum.NumericIntegration.stopd);

    this.rlist = new Array(limit + 1);
    this.rlist.fill(zNum.NumericIntegration.stopd);

    this.elist = new Array(limit + 1);
    this.elist.fill(zNum.NumericIntegration.stopd);

    this.iord = new Array(limit + 1); //int
    this.iord.fill(zNum.NumericIntegration.stopi);

};

zNum.NumericIntegration.RO_rdqk21 = function (result, abserr, resabs, resasc, rdqk21_count) {
    /*
    rdqk21 : function(f, a, b, out double result, out double abserr, out double resabs, out double resasc, ref int rdqk21_count)
    */
    this.result__1 = result;
    this.abserr__2 = abserr;
    this.resabs__3 = resabs;
    this.resasc__4 = resasc;
    this.rdqk21_count__5 = rdqk21_count;
};

zNum.NumericIntegration.RO_rdqelg = function (n, result, abserr, nres) {
    /*
    private static void rdqelg(ref int n, double[] epstab, out double result, out double abserr, double[] res3la, ref int nres)
    */
    this.n__1 = n;
    this.result__2 = result;
    this.abserr__3 = abserr;
    this.nres__4 = nres;
};

zNum.NumericIntegration.RO_rdqpsrt = function (maxerr, ermax, nrmax) {
    /*
    private static void rdqpsrt(int limit, int last, ref int maxerr, out double ermax, double[] elist, int[] iord, ref int nrmax)
    */
    this.maxerr__1 = maxerr;
    this.ermax__2 = ermax;
    this.nrmax__3 = nrmax;
};

/**
* @constructor
* @property {number} result - numerical value 
* @property {number} ier - int indication the result - 0 === Ok
*/
zNum.NumericIntegration.Result = function (result, abserr, neval, ier, limit, last) {

    this.result = result;
    this.abser = abserr;
    this.neval = neval;
    this.ier = ier;
    this.limit = limit;
    this.last = last;
    /*
    result - double, approximation to the integral (value)
    abserr - double, estimate of the modulus of the absolute error, which should equal or exceed abs(i-result)
    neval  - int, number of integrand evaluations
    ier     int, 
            ier = 0 normal and reliable termination of the routine...
            ier = 1 maximum number of subdivisions allowed has been achieved. 
            ier = 2 the occurrence of roundoff error is detected ...
            ier = 3 extremely bad integrand behaviour occurs at some points of the integration interval.
            ier = 4 the algorithm does not converge...
            ier = 5 the integral is probably divergent, or slowly convergent. ...
            ier = 6 the input is invalid
    limit - int, limit determines the maximum number of subintervals
    lenw  - int, >= limit * 4
    last  - int, on return, last equals the number of subintervals ...
    */
};

zNum.NumericIntegration.Result.prototype = {
    messages: [
                "normal and reliable termination of the routine ...",
                "maximum number of subdivisions allowed has been achieved.",
                "the occurrence of roundoff error is detected ...",
                "extremely bad integrand behaviour occurs at some points of the integration interval.",
                "the algorithm does not converge...",
                "the integral is probably divergent, or slowly convergent ...",
                "the input is invalid"],
    getMessage: function () {
        return this.messages[this.ier];
    }
};

zNum.NumericIntegration.prototype = {

    accuracy: 0.0001220703125, // 2.220446049250313e-16 ^ .25 == .Machine$double.eps ^ .25
    limit: 100,
    epsabs: 0.0001220703125,
    epsrel: 0.0001220703125,

    wg: [
            .066671344308688137593568809893332,
            .149451349150580593145776339657697,
            .219086362515982043995534934228163,
            .269266719309996355091226921569469,
            .295524224714752870173892994651338
    ],

    xgk: [
    .995657163025808080735527280689003,
    .973906528517171720077964012084452,
    .930157491355708226001207180059508,
    .865063366688984510732096688423493,
    .780817726586416897063717578345042,
    .679409568299024406234327365114874,
    .562757134668604683339000099272694,
    .433395394129247190799265943165784,
    .294392862701460198131126603103866,
    .14887433898163121088482600112972,
    0.0
    ],

    wgk: [
    .011694638867371874278064396062192,
    .03255816230796472747881897245939,
    .05475589657435199603138130024458,
    .07503967481091995276704314091619,
    .093125454583697605535065465083366,
    .109387158802297641899210590325805,
    .123491976262065851077958109831074,
    .134709217311473325928054001771707,
    .142775938577060080797094273138717,
    .147739104901338491374841515972068,
    .149445554002916905664936468389821
    ],


    compute: function () {

        /*
            f : double, function subprogram defining the integrand
            epsabs - double, absolute accuracy requested, > 0.0
            epsrel - double precision, relative accuracy requested
    
    
            on return
            result - double, approximation to the integral (value)
            abserr - double, estimate of the modulus of the absolute error, which should equal or exceed abs(i-result)
            neval  - int, number of integrand evaluations
            ier     int, 
                    ier = 0 normal and reliable termination of the routine...
                    ier = 1 maximum number of subdivisions allowed has been achieved. 
                    ier = 2 the occurrence of roundoff error is detected ...
                    ier = 3 extremely bad integrand behaviour occurs at some points of the integration interval.
                    ier = 4 the algorithm does not converge...
                    ier = 5 the integral is probably divergent, or slowly convergent. ...
                    ier = 6 the input is invalid
            limit - int, limit determines the maximum number of subintervals
            lenw  - int, >= limit * 4
            last  - int, on return, last equals the number of subintervals ...
    
            work arrays
                iwork - int
                        vector of dimension at least limit, the first k
                        elements of which contain pointers
                        to the error estimates over the subintervals
                        such that work(limit*3+iwork(1)),... ,
                        work(limit*3+iwork(k)) form a decreasing
                        sequence, with k = last if last <= (limit/2+2),
                        and k = limit+1-last otherwise
    
                work  - double precision
                        vector of dimension at least lenw
                        on return
                        work(1), ..., work(last) contain the left
                        end-points of the subintervals in the
                        partition of (a,b),
                        work(limit+1), ..., work(limit+last) contain
                        the right end-points,
                        work(limit*2+1), ..., work(limit*2+last) contain
                        the integral approximations over the subintervals,
                        work(limit*3+1), ..., work(limit*3+last)
                        contain the error estimates.
        */

        var a = this.a;
        var b = this.b;

        if (this.limit < 1) {
            return;
        }


        // Workspace 


        //this.alist = new Array(this.limit + 1);
        //this.blist = new Array(this.limit + 1);
        //this.rlist = new Array(this.limit + 1);
        //ws.elist = new Array(this.limit + 1);
        //this.iord = new Array(this.limit + 1);
        //this.alist.fill(0.0);
        //this.blist.fill(0.0);
        //this.rlist.fill(0.0);
        //ws.elist.fill(0.0);
        //this.iord.fill(0.0);

        return (this.rdqagse());

    },


    rdqagse: function () {
        /* Local variables */

        var ws = new zNum.NumericIntegration.Workspace(this.limit);


        var jump = {
            continue: true,
            label: -1
        };

        var noext = false, extrap = false;

        var  //int
            k, ksgn = 0, nres = 0,
            ierro,
            ktmin = 0, nrmax = 0,
            iroff1 = 0, iroff2 = 0, iroff3 = 0,
            id,
            numrl2 = 0,
            jupbnd,
            maxerr = 0,
            last,
            ier,
            neval,
            rdqk21_count = 0;

        var // number[]
            res3la = new Array(4),
            rlist2 = new Array(53);

        var // double
            result = 0.0,
            abserr = 0.0,
            abseps,
            area = 0,
            area1,
            area2,
            area12,
            dres,
            epmach;

        var // double
            a1,
            a2,
            b1,
            b2,
            defabs,
            defab1,
            defab2,
            oflow,
            uflow,
            resabs,
            reseps;

        var // double
            error1,
            error2,
            erro12,
            errbnd,
            erlast,
            errmax = 0,
            errsum = 0;

        var // double
            correc = 0.0,
            erlarg = 0.0,
            ertest = 0.0,
            small = 0.0;

        //Let's go
        epmach = zU.CONST.DBL_EPSILON;
        ier = 0;
        neval = 0;
        last = 0;



        ws.alist[1] = this.a;
        ws.blist[1] = this.b;
        ws.rlist[1] = 0.0;
        ws.elist[1] = 0.0;

        //TODO déplacer(?) dans compute
        if ((this.epsabs <= 0.0) && (this.epsrel < Math.max(epmach * 50.0, 5e-29))) {
            ier = 6;
            return;
        }
        /*           first approximation to the integral */
        /*           ----------------------------------- */

        uflow = zU.CONST.DBL_MIN;
        oflow = zU.CONST.DBL_MAX;
        ierro = 0;

        // result, abserr sont des propriétés
        // out defabs, out resabs, ref rdqk21_count sont des variables locales


        //boxing
        var ro_rdqk21 = new zNum.NumericIntegration.RO_rdqk21(result, abserr, defabs, resabs, rdqk21_count);

        //rdqk21(parameters.f, parameters.a, parameters.b, out result, out abserr, out defabs, out resabs, ref rdqk21_count);
        this.rdqk21(this.fn, this.a, this.b, ro_rdqk21);

        //unboxing
        result = ro_rdqk21.result__1;
        abserr = ro_rdqk21.abserr__2;
        defabs = ro_rdqk21.resabs__3;
        resabs = ro_rdqk21.resasc__4;
        rdqk21_count = ro_rdqk21.rdqk21_count__5;



        /*           test on accuracy. */

        dres = Math.abs(result);
        errbnd = Math.max(this.epsabs, this.epsrel * dres);
        last = 1;
        ws.rlist[1] = result;
        ws.elist[1] = abserr;
        ws.iord[1] = 1;
        if ((abserr <= (epmach * 100.0 * defabs)) && (abserr > errbnd))
            ier = 2;
        if (this.limit === 1)
            ier = 1;
        if ((ier !== 0) || (abserr <= errbnd && abserr !== resabs) || (abserr === 0.0)) {
            jump.continue = false;
            jump.label = 140;
        }

        /*           initialization */
        /*           -------------- */

        if (jump.continue) {
            rlist2[0] = result;
            errmax = abserr;
            maxerr = 1;
            area = result;
            errsum = abserr;
            abserr = oflow;
            nrmax = 1;
            nres = 0;
            numrl2 = 2;
            ktmin = 0;
            extrap = false;
            noext = false;
            iroff1 = 0;
            iroff2 = 0;
            iroff3 = 0;
            ksgn = -1;
            if (dres >= (1.0 - epmach * 50.0) * defabs) {
                ksgn = 1;
            }
        }
        /*------------------------*/


        for (last = 2; last <= this.limit; ++(last)) {
            if (!jump.continue) {
                // on sort immmédiatement
                break;
            }

            /*           bisect the subinterval with the nrmax-th largest error estimate. */

            a1 = ws.alist[maxerr];
            b1 = (ws.alist[maxerr] + ws.blist[maxerr]) * .5;
            a2 = b1;
            b2 = ws.blist[maxerr];
            erlast = errmax;
            //boxing 
            ro_rdqk21 = new zNum.NumericIntegration.RO_rdqk21(area1, error1, resabs, defab1, rdqk21_count);
            this.rdqk21(this.fn, a1, b1, ro_rdqk21);
            //unboxing
            area1 = ro_rdqk21.result__1;
            error1 = ro_rdqk21.abserr__2;
            resabs = ro_rdqk21.resabs__3;
            defab1 = ro_rdqk21.resasc__4;
            rdqk21_count = ro_rdqk21.rdqk21_count__5;

            //boxing
            ro_rdqk21 = new zNum.NumericIntegration.RO_rdqk21(area2, error2, resabs, defab2, rdqk21_count);
            this.rdqk21(this.fn, a2, b2, ro_rdqk21);
            //unboxing
            area2 = ro_rdqk21.result__1;
            error2 = ro_rdqk21.abserr__2;
            resabs = ro_rdqk21.resabs__3;
            defab2 = ro_rdqk21.resasc__4;
            rdqk21_count = ro_rdqk21.rdqk21_count__5;

            /*           improve previous approximations to integral
                     and error and test for accuracy. */

            area12 = area1 + area2;
            erro12 = error1 + error2;
            errsum = errsum + erro12 - errmax;
            area = area + area12 - ws.rlist[maxerr];
            if (!(defab1 === error1 || defab2 === error2)) {
                if (Math.abs(ws.rlist[maxerr] - area12) <= Math.abs(area12) * 1e-5 &&
                erro12 >= errmax * .99) {
                    if (extrap)
                        ++iroff2;
                    else /* if(! extrap) */
                        ++iroff1;
                }
                if (last > 10 && erro12 > errmax)
                    ++iroff3;
            }
            ws.rlist[maxerr] = area1;
            ws.rlist[last] = area2;
            errbnd = Math.max(this.epsabs, this.epsrel * Math.abs(area));

            /*           test for roundoff error and eventually set error flag. */

            if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
                ier = 2;

            if (iroff2 >= 5)
                ierro = 3;

            /* set error flag in the case that the number of subintervals equals limit. */
            if (last === this.limit)
                ier = 1;

            /*           set error flag in the case of bad integrand behaviour
                     at a point of the integration range. */

            if (Math.max(Math.abs(a1), Math.abs(b2)) <=
                (epmach * 100.0 + 1.0) * (Math.abs(a2) + uflow * 1e3)) {
                ier = 4;
            }

            /*           append the newly-created intervals to the list. */

            //console.log(zU.fmt("a1={0}, b1={1}, area1={2}, error1={3}", a1, b1, area1, error1));
            //console.log(zU.fmt("a2={0}, b2={1}, area2={2}, error2={3}", a2, b2, area2, error2));

            if (error2 > error1) {
                ws.alist[maxerr] = a2;
                ws.alist[last] = a1;
                ws.blist[last] = b1;
                ws.rlist[maxerr] = area2;
                ws.rlist[last] = area1;
                ws.elist[maxerr] = error2;
                ws.elist[last] = error1;
            }
            else {
                ws.alist[last] = a2;
                ws.blist[maxerr] = b1;
                ws.blist[last] = b2;
                ws.elist[maxerr] = error1;
                ws.elist[last] = error2;
            }




            /*  call subroutine dqpsrt to maintain the descending ordering
                in the list of error estimates and select the subinterval
                with nrmax-th largest error estimate (to be bisected next). */

            /*L30:*/

            /*

               this.maxerr__1 = maxerr;
    this.ermax__2 = ermax;
    this.nrmax__3 = nrmax;
            */
            //boxing
            var ro_rdqpsrt = new zNum.NumericIntegration.RO_rdqpsrt(maxerr, errmax, nrmax);
            this.rdqpsrt(this.limit, last, ws.elist, ws.iord, ro_rdqpsrt);
            //unboxing
            maxerr = ro_rdqpsrt.maxerr__1;
            errmax = ro_rdqpsrt.ermax__2;
            nrmax = ro_rdqpsrt.nrmax__3;

            displayArray = zNum.NumericIntegration.displayArray;
            //console.log(displayArray("alist", ws.alist, 100));
            //console.log(displayArray("blist", ws.blist, 100));
            //console.log(displayArray("elist", ws.elist, 100));
            //console.log(displayArray("rlist", ws.rlist, 100));
            //console.log(displayArray("iord", ws.iord, 100));




            if (errsum <= errbnd) {
                // goto L115 
                jump.continue = false;
                jump.label = 115;
                /* ***jump out of do-loop */
                break; // on sort de for (last = 2; last <= this.limit; ++(last)) {
            }

            if (ier !== 0) break; // on sort de for (last = 2; last <= this.limit; ++(last)) {
            if (last === 2) { /* L80: */
                small = Math.abs(this.b - this.a) * .375;
                erlarg = errsum;
                ertest = errbnd;
                rlist2[1] = area;
                continue;  // continue du for (last = 2; last <= this.limit; ++(last)) {
            }

            if (noext)
                continue; // continue du for (last = 2; last <= this.limit; ++(last)) {

            erlarg -= erlast;
            if (Math.abs(b1 - a1) > small) {
                erlarg += erro12;
            }
            if (!extrap) {

                /*          test whether the interval to be bisected next is the
                        smallest interval. */

                if (Math.abs(ws.blist[maxerr] - ws.alist[maxerr]) > small) {
                    continue;
                }
                extrap = true;
                nrmax = 2;
            }

            if (ierro !== 3 && erlarg > ertest) {

                /*           the smallest interval has the largest error.
                         before bisecting decrease the sum of the errors over the
                         larger intervals (erlarg) and perform extrapolation. */

                id = nrmax;
                jupbnd = last;
                if (last > Math.floor(this.limit / 2) + 2) {
                    jupbnd = this.limit + 3 - last;
                }
                for (k = id; k <= jupbnd; ++k) {
                    maxerr = ws.iord[nrmax];
                    errmax = ws.elist[maxerr];
                    if (Math.abs(ws.blist[maxerr] - ws.alist[maxerr]) > small) {
                        // goto L90;
                        jump.continue = false;
                        jump.label = 90;
                        break; // on sort de 'for (k = id; k <= jupbnd; ++k) ...'
                    }
                    ++nrmax;
                    /* L50: */
                }
            }
            /*           perform extrapolation.  L60: */

            // si on a un goto 90 actif on n'exécutera pas le code du bloc suivant
            if (jump.continue) {
                ++numrl2;
                rlist2[numrl2 - 1] = area;

                //boxing
                var ro_rdqelg = new zNum.NumericIntegration.RO_rdqelg(numrl2, reseps, abseps, nres);
                this.rdqelg(rlist2, res3la, ro_rdqelg);
                // unboxing
                numrl2 = ro_rdqelg.n__1;
                reseps = ro_rdqelg.result__2;
                abseps = ro_rdqelg.abserr__3;
                nres = ro_rdqelg.nres__4;

                ++ktmin;
                if (ktmin > 5 && abserr < errsum * .001) {
                    ier = 5;
                }
                if (abseps < abserr) {
                    ktmin = 0;
                    abserr = abseps;
                    result = reseps;
                    correc = erlarg;
                    ertest = Math.max(this.epsabs, this.epsrel * Math.abs(reseps));
                    if (abserr <= ertest) {
                        break;
                    }
                }

                /*           prepare bisection of the smallest interval.  L70: */

                if (numrl2 === 1) {
                    noext = true;
                }

                if (ier === 5) {
                    break;
                }

                maxerr = ws.iord[1];
                errmax = ws.elist[maxerr];
                nrmax = 1;
                extrap = false;
                small *= 0.5;
                erlarg = errsum;
            }

            // L90
            if (!jump.continue && jump.label === 90) {
                jump.continue = true;
            }
            L90:;
        } //for( last = 2 ...)

        /* L100:	set final result and error estimate. */
        /*		------------------------------------ */
        while (jump.continue) {
            if (abserr === oflow) {
                // goto L115;
                jump.continue = false;
                jump.label = 115;
                break;
            }

            if (ier + ierro === 0) {
                // goto L110;
                jump.continue = false;
                jump.label = 110;
                break;
            }

            if (ierro === 3)
                abserr += correc;

            if (ier === 0)
                ier = 3;

            if (result === 0.0 || area === 0.0) {
                if (abserr > errsum) {
                    // goto L115;
                    jump.continue = false;
                    jump.label = 115;
                    break;
                }

                if (area === 0.0) {
                    //goto L130;
                    jump.continue = false;
                    jump.label = 130;
                    break;
                }
            }
            else { /* L105:*/
                if (abserr / Math.abs(result) > errsum / Math.abs(area)) {
                    // goto L115;
                    jump.continue = false;
                    jump.label = 115;
                    break;
                }
            }

            jump.continue = true;
            break;
        } // while (jump.Continue)


        while (jump.continue || (jump.label === 110)) {
            /* L110: */
            /*		test on divergence. */
            if (ksgn === -1 && Math.max(Math.abs(result), Math.abs(area)) <= defabs * .01) {
                //goto L130;
                jump.continue = false;
                jump.label = 130;
                break;
            }

            if (.01 > result / area || result / area > 100.0 || errsum > Math.abs(area)) {
                ier = 5;
            }

            //goto L130;
            jump.continue = false;
            jump.label = 130;
            break;
        }

        if (jump.continue || (jump.label === 115)) {
            /*L115:*/
            /*		compute global integral sum. */
            result = 0.0;
            for (k = 1; k <= last; ++k)
                result += ws.rlist[k];
            abserr = errsum;
            jump.continue = true;
        }

        if (jump.continue || (jump.label === 130) || (jump.label === 139) || (jump.label === 140)) {
            /* L130:
               L139:
               L140: */

            if ((ier > 2) || jump.label === 140) {
                /*L140:*/
                neval = last * 42 - 21;
            }
            this.abserr = abserr;
            this.ier = ier;
            this.last = last;
            this.neval = neval;
            this.result = result;

        }

        // encapsuler les résultats

        return new zNum.NumericIntegration.Result(this.result, this.abserr, this.neval, this.ier, this.limit, this.last);

        // ---------------------------*/
    }, // rdqagse

    /* R file: integrate.c
     */

    /*
    rdqk21(this.fn, this.a, this.b, out result, out abserr, out defabs, out resabs, ref rdqk21_count);
    rdqk21(this.f, a1, b1, out area1, out error1, out resabs, out defab1, ref rdqk21_count);
    rdqk21(this.f, a2, b2, out area2, out error2, out resabs, out defab2, ref rdqk21_count);
    */

    /**
    * @function
    * @param {zygotine.NumericIntegration.RO_rdqk21} ro
    */
    rdqk21: function (f, a, b, ro /* out double result, out double abserr, out double resabs, out double resasc, ref int rdqk21_count*/) {

        //static void  rdqk21(integr_fn f, void *ex, double *a, double *b, double *result,
        // double *abserr, double *resabs, double *resasc)
        var // {number[]} double
            fv1 = new Array(10),
            fv2 = new Array(10),
            vec = new Array(21);

        var // double 
            absc, resg, resk, fsum, fval1, fval2,
            hlgth, centr, reskh, uflow,
            fc, epmach, dhlgth;

        var
            j, jtw, jtwm1;


        ro.rdqk21_count__5++;
        epmach = zU.CONST.DBL_EPSILON;
        uflow = zU.CONST.DBL_MIN;
        centr = (a + b) * 0.5;
        hlgth = (b - a) * 0.5;
        dhlgth = Math.abs(hlgth);

        /* compute the 21-point kronrod approximation to the integral, and estimate the absolute error. */

        resg = 0.0;
        vec[0] = centr;
        for (j = 1; j <= 5; ++j) {
            jtw = j << 1;
            absc = hlgth * this.xgk[jtw - 1];
            vec[(j << 1) - 1] = centr - absc;

            L5:
                vec[j * 2] = centr + absc;
        }

        for (j = 1; j <= 5; ++j) {
            jtwm1 = (j << 1) - 1;
            absc = hlgth * this.xgk[jtwm1 - 1];
            vec[(j << 1) + 9] = centr - absc;
            vec[(j << 1) + 10] = centr + absc;
        }

        var vec_y = f(vec);//, 21); // un tableau de longueur 21

        fc = vec_y[0];
        resk = this.wgk[10] * fc;
        ro.resabs__3 = Math.abs(resk);

        for (j = 1; j <= 5; ++j) {
            jtw = j << 1;
            absc = hlgth * this.xgk[jtw - 1];
            fval1 = vec_y[(j << 1) - 1];
            fval2 = vec_y[j * 2];
            fv1[jtw - 1] = fval1;
            fv2[jtw - 1] = fval2;
            fsum = fval1 + fval2;
            resg += this.wg[j - 1] * fsum;
            resk += this.wgk[jtw - 1] * fsum;
            ro.resabs__3 += this.wgk[jtw - 1] * (Math.abs(fval1) + Math.abs(fval2));
            //L10: 
        }


        for (j = 1; j <= 5; ++j) {
            jtwm1 = (j << 1) - 1;
            absc = hlgth * this.xgk[jtwm1 - 1];
            fval1 = vec_y[(j << 1) + 9];
            fval2 = vec_y[(j << 1) + 10];
            fv1[jtwm1 - 1] = fval1;
            fv2[jtwm1 - 1] = fval2;
            fsum = fval1 + fval2;
            resk += this.wgk[jtwm1 - 1] * fsum;
            ro.resabs__3 += this.wgk[jtw
                - 1] * (Math.abs(fval1) + Math.abs(fval2));
            // L15:
        }


        reskh = resk * .5;
        ro.resasc__4 = this.wgk[10] * Math.abs(fc - reskh);
        for (j = 1; j <= 10; ++j) {
            ro.resasc__4 += this.wgk[j - 1] * (Math.abs(fv1[j - 1] - reskh) +
                         Math.abs(fv2[j - 1] - reskh));
            /* L20: */
        }
        //vec[0] = 0;

        ro.result__1 = resk * hlgth;
        ro.resabs__3 *= dhlgth;
        ro.resasc__4 *= dhlgth;
        ro.abserr__2 = Math.abs((resk - resg) * hlgth);
        if ((ro.resasc__4 !== 0.0) && (ro.abserr !== 0.0)) {
            ro.abserr__2 = ro.resasc__4 * Math.min(1.0, Math.pow(ro.abserr__2 * 200.0 / ro.resasc__4, 1.5));
        }

        if (ro.resabs__3 > uflow / (epmach * 50.0)) {
            ro.abserr__2 = Math.max(epmach * 50.0 * ro.resabs__3, ro.abserr__2);
        }


        return;
        // throw new System.NotImplementedException();
    }, // rdqk21


    rdqelg: function (epstab, res3la, ro) {
        /*  private static void rdqelg(ref int n, double[] epstab, out double result, out double abserr, double[] res3la, ref int nres)
                n <-> ro.n__1;
                result <-> ro.result__2;
                abserr <->  ro.abserr__3;
                nres <->  ro.nres__4;
         */
        var jump = {
            continue: true,
            label: -1
        };

        /* Local variables */
        var  // int
            i__, indx, ib, ib2, ie, k1, k2, k3, num = 0, newelm = 0, limexp = 0;
        var // double
            delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach, epsinf;
        var // double 
             oflow, ss = 0, res;
        var // double 
            errA, err1, err2, err3, tol1, tol2, tol3;

        /*  ***begin prologue  dqelg
            ***refer to  dqagie,dqagoe,dqagpe,dqagse
            ***revision date  830518   (yymmdd)
            ***keywords  epsilon algorithm, convergence acceleration,
            extrapolation
            ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
            de doncker,elise,appl. math & progr. div. - k.u.leuven
            ***purpose  the routine determines the limit of a given sequence of
            approximations, by means of the epsilon algorithm of
            p.wynn. an estimate of the absolute error is also given.
            the condensed epsilon table is computed. only those
            elements needed for the computation of the next diagonal
            are preserved.
            ***description
        
            epsilon algorithm
            standard fortran subroutine
            double precision version
        
            parameters
            n      - int
            epstab(n) contains the new element in the
            first column of the epsilon table.
        
            epstab - double precision
            vector of dimension 52 containing the elements
            of the two lower diagonals of the triangular
            epsilon table. the elements are numbered
            starting at the right-hand corner of the
            triangle.
        
            result - double precision
            resulting approximation to the integral
        
            abserr - double precision
            estimate of the absolute error computed from
            result and the 3 previous results
        
            res3la - double precision
            vector of dimension 3 containing the last 3
            results
        
            nres   - int
            number of calls to the routine
            (should be zero at first call)
        
            ***end prologue  dqelg
        
        
            list of major variables
            -----------------------
        
            e0     - the 4 elements on which the computation of a new
            e1       element in the epsilon table is based
            e2
            e3                 e0
            e3    e1    new
            e2
        
            newelm - number of elements to be computed in the new diagonal
            errA   - errA = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
            ro.result__2 - the element in the new diagonal with least value of errA
        
            machine dependent constants
            ---------------------------
        
            epmach is the largest relative spacing.
            oflow is the largest positive magnitude.
            limexp is the maximum number of elements the epsilon
            table can contain. if this number is reached, the upper
            diagonal of the epsilon table is deleted. */
        // -res3la;
        //--epstab;

        /* Function Body */
        epmach = zU.CONST.DBL_EPSILON;
        oflow = zU.CONST.DBL_MAX;
        ++ro.nres__4;
        ro.abserr__3 = oflow;
        ro.result__2 = epstab[ro.n__1 - 1];
        if (ro.n__1 < 3) {
            // goto L100;
            jump.continue = false;
            jump.label = 100;
        }

        if (jump.continue) {
            limexp = 50;
            epstab[ro.n__1 + 2 - 1] = epstab[ro.n__1 - 1];
            newelm = Math.floor((ro.n__1 - 1) / 2);
            epstab[ro.n__1 - 1] = oflow;
            num = ro.n__1;
            k1 = ro.n__1;
            for (i__ = 1; i__ <= newelm; ++i__) {
                k2 = k1 - 1;
                k3 = k1 - 2;
                res = epstab[k1 + 2 - 1];
                e0 = epstab[k3 - 1];
                e1 = epstab[k2 - 1];
                e2 = res;
                e1abs = Math.abs(e1);
                delta2 = e2 - e1;
                err2 = Math.abs(delta2);
                tol2 = Math.max(Math.abs(e2), e1abs) * epmach;
                delta3 = e1 - e0;
                err3 = Math.abs(delta3);
                tol3 = Math.max(e1abs, Math.abs(e0)) * epmach;
                if (err2 <= tol2 && err3 <= tol3) {
                    /*           if e0, e1 and e2 are equal to within machine
                    accuracy, convergence is assumed. */
                    ro.result__2 = res;/*		result = e2 */
                    ro.abserr__3 = err2 + err3;/*	abserr = fabs(e1-e0)+fabs(e2-e1) */

                    jump.continue = false;
                    jump.label = 100;
                    break;
                    // goto L100;  /* ***jump out of do-loop */
                }


                e3 = epstab[k1 - 1];
                epstab[k1 - 1] = e1;
                delta1 = e1 - e3;
                err1 = Math.abs(delta1);
                tol1 = Math.max(e1abs, Math.abs(e3)) * epmach;

                /*           if two elements are very close to each other, omit
                a part of the table by adjusting the value of n */

                if (err1 > tol1 && err2 > tol2 && err3 > tol3) {
                    ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
                    epsinf = Math.abs(ss * e1);

                    /*           test to detect irregular behaviour in the table, and
                    eventually omit a part of the table adjusting the value of n. */

                    if (epsinf > 1e-4) {
                        // goto L30;
                        jump.continue = false;
                        jump.label = 30;
                    }
                } // if (err1 > tol1 && err2 > tol2 && err3 > tol3)

                //il faut sauter le bloc suivant lorsqu'un goto L30 est actif.
                if (jump.continue) {
                    ro.n__1 = i__ + i__ - 1;
                    // goto L50
                    jump.continue = false;
                    jump.label = 50;
                    break; // on sort de la boucle
                    // goto L50; /* ***jump out of do-loop */
                }

                if (jump.continue || (jump.label === 30)) {
                    jump.continue = true; // on remet à true pour satisfaire le goto L30
                    /*L30:*/
                    /* compute a new element and eventually adjust the value of ro.result__2. */

                    res = e1 + 1.0 / ss;
                    epstab[k1 - 1] = res;
                    k1 += -2;
                    errA = err2 + Math.abs(res - e2) + err3;
                    if (errA <= ro.abserr__3) {
                        ro.abserr__3 = errA;
                        ro.result__2 = res;
                    }
                }

            } // for (i__ = 1; i__ <= newelm; ++i__)

            /*           shift the table. */
        } // if( jump.continue) 

        if (jump.continue || (jump.label === 50)) {
            /* L50:*/

            if (ro.n__1 === limexp) {
                ro.n__1 = (limexp / 2 << 1) - 1;
            }

            if (num / 2 << 1 === num) ib = 2; else ib = 1;
            ie = newelm + 1;
            for (i__ = 1; i__ <= ie; ++i__) {
                ib2 = ib + 2;
                epstab[ib - 1] = epstab[ib2 - 1];
                ib = ib2;
            }

            if (num !== ro.n__1) {
                indx = num - ro.n__1 + 1;
                for (i__ = 1; i__ <= ro.n__1; ++i__) {
                    epstab[i__ - 1] = epstab[indx - 1];
                    ++indx;
                }
            }

            /*L80:*/
            if (ro.nres__4 >= 4) {
                /* L90: */
                ro.abserr__3 = Math.abs(ro.result__2 - res3la[3]) +
                    Math.abs(ro.result__2 - res3la[2]) +
                    Math.abs(ro.result__2 - res3la[1]);
                res3la[1] = res3la[2];
                res3la[2] = res3la[3];
                res3la[3] = ro.result__2;
            }
            else {
                res3la[ro.nres__4] = ro.result__2;
                ro.abserr__3 = oflow;
            }
        }

        L100:/* compute error estimate */
            ro.abserr__3 = Math.max(ro.abserr__3, epmach * 5.0 * Math.abs(ro.result__2));
        return;
    }, // rdqelg

    rdqpsrt: function (limit, last, elist, iord, ro) {
        /*
        ro.maxerr__1 = maxerr;
        ro.ermax__2 = ermax;
        ro.nrmax__3 = nrmax;
        */

        var
            gotoLast = false;
        /* Local variables */
        var //int 
            i, j, k, ido, jbnd, isucc, jupbn;
        var //double 
            errmin, errmax;
        /* ***begin prologue  dqpsrt
         ***refer to  dqage,dqagie,dqagpe,dqawse
         ***routines called  (none)
         ***revision date  810101   (yymmdd)
         ***keywords  sequential sorting
         ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
                   de doncker,elise,appl. math. & progr. div. - k.u.leuven
         ***purpose  this routine maintains the descending ordering in the
                    list of the local error estimated resulting from the
                    interval subdivision process. at each call two error
                    estimates are inserted using the sequential search
                    method, top-down for the largest error estimate and
                    bottom-up for the smallest error estimate.
         ***description
        
                   ordering routine
                   standard fortran subroutine
                   double precision version
        
                   parameters (meaning at output)
                      limit  - int
                               maximum number of error estimates the list
                               can contain
        
                      last   - int
                               number of error estimates currently in the list
        
                      maxerr - int
                               maxerr points to the nrmax-th largest error
                               estimate currently in the list
        
                      ermax  - double precision
                               nrmax-th largest error estimate
                               ermax = elist(maxerr)
        
                      elist  - double precision
                               vector of dimension last containing
                               the error estimates
        
                      iord   - int
                               vector of dimension last, the first k elements
                               of which contain pointers to the error
                               estimates, such that
                               elist(iord(1)),...,  elist(iord(k))
                               form a decreasing sequence, with
                               k = last if last <= (limit/2+2), and
                               k = limit+1-last otherwise
        
                      nrmax  - int
                               maxerr = iord(nrmax)
        
        ***end prologue  dqpsrt
        */

        /* Parameter adjustments */
        // mis en commentaires par zygotine
        //--iord;
        //--elist;

        /* Function Body */

        /*           check whether the list contains more than
                 two error estimates. */
        if (last <= 2) {
            iord[1] = 1;
            iord[2] = 2;
            gotoLast = true;
        }
        /*           this part of the routine is only executed if, due to a
                 difficult integrand, subdivision increased the error
                 estimate. in the normal case the insert procedure should
                 start after the nrmax-th largest error estimate. */
        while (!gotoLast) {
            errmax = elist[ro.maxerr__1];
            if (ro.nrmax__3 > 1) {
                ido = ro.nrmax__3 - 1;
                for (i = 1; i <= ido; ++i) {
                    isucc = iord[ro.nrmax__3 - 1];
                    if (errmax <= elist[isucc])
                        break; /* out of for-loop */
                    iord[ro.nrmax__3] = isucc;
                    --ro.nrmax__3;
                    /* L20: */
                }
            }

            /*L30:       compute the number of elements in the list to be maintained
                     in descending order. this number depends on the number of
                     subdivisions still allowed. */
            if (last > Math.floor(limit / 2) + 2)
                jupbn = limit + 3 - last;
            else
                jupbn = last;

            errmin = elist[last];

            /*           insert errmax by traversing the list top-down,
                     starting comparison from the element elist(iord(nrmax+1)). */

            jbnd = jupbn - 1;
            for (i = ro.nrmax__3 + 1; i <= jbnd; ++i) {
                isucc = iord[i];
                if (errmax >= elist[isucc]) {/* ***jump out of do-loop */
                    /* L60: insert errmin by traversing the list bottom-up. */
                    iord[i - 1] = ro.maxerr__1;
                    for (j = i, k = jbnd; j <= jbnd; j++, k--) {
                        isucc = iord[k];
                        if (errmin < elist[isucc]) {
                            /* goto L80; ***jump out of do-loop */
                            iord[k + 1] = last;

                            gotoLast = true;
                            break; // on sort de for (j = i, k = jbnd; j <= jbnd; j++, k--) {
                        }
                        iord[k + 1] = isucc;
                    }

                    if (gotoLast) {
                        break; // on sort de  for (i = ro.nrmax__3 + 1; i <= jbnd; ++i) {
                    }

                    iord[i] = last;
                    gotoLast = true; // // on sort de  for (i = ro.nrmax__3 + 1; i <= jbnd; ++i) {
                    break;
                }

                iord[i - 1] = isucc;
            } // for (i = ro.nrmax__3 + 1; i <= jbnd; ++i) {

            if (gotoLast) {
                break; // on sort de while (!gotoLast) {
            }
            iord[jbnd] = ro.maxerr__1;
            iord[jupbn] = last;
            break; // il faut sortir du while
        } //while (!gotoLast)

        /* Last: */
        /* set maxerr and ermax. */
        ro.maxerr__1 = iord[ro.nrmax__3];
        ro.ermax__2 = elist[ro.maxerr__1];
        return;
    } /* rdqpsrt_ */


};

/**
* @function 
* @param {number} x
* @param {number[]} a
* @param {number} n
* 
* @returns {number}
* Fichier R : chebyshev.c 
*/
zNum.chebyshev_eval = function (x, a, n) {
    var //double
        b0, b1, b2, twox;
    var
        i;

    if (n < 1 || n > 1000) {
        return NaN;
    }

    if (x < -1.1 || x > 1.1) {
        return NaN;
    }

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }

    return (b0 - b2) * 0.5;
};

/* R file: stirlerr.c
 * name used in R: stirlerr
 */
zNum.stirlingError = function (n) {
    /*
      error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
    */
    var nn;

    if (n <= 15.0) {
        nn = n + n;
        var nnTruncate = Math.trunc(nn);
        if (nn === nnTruncate) {
            return (zNum.stirlingError.sferr_halves[nnTruncate]);
        }

        return (zNum.logGamma(n + 1.0) - (n + 0.5) * Math.log(n) + n - zU.CONST.M_LN_SQRT_2PI);
    }

    var S = zNum.stirlingError.S;
    nn = n * n;
    if (n > 500) return ((S[0] - S[1] / nn) / n);
    if (n > 80) return ((S[0] - (S[1] - S[2] / nn) / nn) / n);
    if (n > 35) return ((S[0] - (S[1] - (S[2] - S[3] / nn) / nn) / nn) / n);
    /* 15 < n <= 35 : */
    return ((S[0] - (S[1] - (S[2] - (S[3] - S[4] / nn) / nn) / nn) / nn) / n);
};

zNum.stirlingError.S = [
/*S0 = */ 0.083333333333333333333,        /* 1/12 */
/*S1 = */ 0.00277777777777777777778,      /* 1/360 */
/*S2 = */ 0.00079365079365079365079365,   /* 1/1260 */
/*S3 = */ 0.000595238095238095238095238,  /* 1/1680 */
/*S4 = */ 0.0008417508417508417508417508  /* 1/1188 */
];

// R file: stirlerr.c
zNum.stirlingError.sferr_halves = [
    0.0, /* n=0 - wrong, place holder only */
    0.1534264097200273452913848,  /* 0.5 */
    0.0810614667953272582196702,  /* 1.0 */
    0.0548141210519176538961390,  /* 1.5 */
    0.0413406959554092940938221,  /* 2.0 */
    0.03316287351993628748511048, /* 2.5 */
    0.02767792568499833914878929, /* 3.0 */
    0.02374616365629749597132920, /* 3.5 */
    0.02079067210376509311152277, /* 4.0 */
    0.01848845053267318523077934, /* 4.5 */
    0.01664469118982119216319487, /* 5.0 */
    0.01513497322191737887351255, /* 5.5 */
    0.01387612882307074799874573, /* 6.0 */
    0.01281046524292022692424986, /* 6.5 */
    0.01189670994589177009505572, /* 7.0 */
    0.01110455975820691732662991, /* 7.5 */
    0.010411265261972096497478567, /* 8.0 */
    0.009799416126158803298389475, /* 8.5 */
    0.009255462182712732917728637, /* 9.0 */
    0.008768700134139385462952823, /* 9.5 */
    0.008330563433362871256469318, /* 10.0 */
    0.007934114564314020547248100, /* 10.5 */
    0.007573675487951840794972024, /* 11.0 */
    0.007244554301320383179543912, /* 11.5 */
    0.006942840107209529865664152, /* 12.0 */
    0.006665247032707682442354394, /* 12.5 */
    0.006408994188004207068439631, /* 13.0 */
    0.006171712263039457647532867, /* 13.5 */
    0.005951370112758847735624416, /* 14.0 */
    0.005746216513010115682023589, /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
];

/*
 * R file: bd0.c
 * 
 *  AUTHOR
 *	Catherine Loader, catherine@research.bell-labs.com.
 *	October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000, The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *	Evaluates the "deviance part"
 *	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
 *		  =  x * log(x/M) + M - x
 *	where M = E[X] = n*p (or = lambda), for	  x, M > 0
 *
 *	in a manner that should be stable (with small relative error)
 *	for all x and M=np. In particular for x/np close to 1, direct
 *	evaluation fails, and evaluation is based on the Taylor series
 *	of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).
 */
zNum.bd0 = function (x, np) {
    var ej, s, s1, v;
    var /*int*/ j;

    if (!isFinite(x) || !isFinite(np) || np === 0.0) {
        return NaN;
    }

    if (Math.abs(x - np) < 0.1 * (x + np)) {
        v = (x - np) / (x + np);  // might underflow to 0
        s = (x - np) * v; /* s using v -- change by MM */
        if (Math.abs(s) < zU.CONST.DBL_MIN) {
            return s;
        }

        ej = 2 * x * v;
        v = v * v;
        for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop as |v| < .1,  v^2000 is "zero" */
            ej *= v;// = v^(2j+1)
            s1 = s + ej / ((j << 1) + 1);
            if (s1 === s) { /* last term was effectively 0 */
                return s1;
            }

            s = s1;
        }
    }

    /* else:  | x - np |  is not too small */
    return (x * Math.log(x / np) + np - x);
};

/* R file: lgamma.c
 * name used in R: lgammafn
 */
zNum.logGamma = function (x) {
    var sign = null;
    // le paramètre sign de lgamma_fn ne devrait pas être modifié suite à l'appel de lgammafn_sign.
    // Don on passe null et on testera pour null (au sens de javascript). 
    //
    // Cette solution de contournement ne serait pas suffisante sous d'autres conditions.
    // Par exemple, la fn choose (fichier choose.c des sources de R) utilise la fn lfastchoose2 qui elle utilise réellement le paramètre sign. 
    // Dans ce cas, la fn choose appelle lfastchoose2 qui, elle, appelle lgammafn_sign.
    // 
    return zNum.logGamma.lgammafn_sign(x, sign);
};

/* R file: gamma.c
 * name used in R: gammafn
 */
zNum.gamma = function (x) {

    var /*int*/ i, n;
    var /*double*/
        y,
        sinpiy,
        value;

    var /*double*/
        xmin = -170.5674972726612,
        xmax = 171.61447887182298,
        xsml = 2.2474362225598545e-308,
        dxrel = 1.490116119384765696e-8;

    var /*int*/
        ngam = 22;


    if (isNaN(x)) {
        return x;
    }


    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x === 0.0 || (x < 0.0 && (x === Math.round(x)))) {
        return NaN;
    }

    y = Math.abs(x);

    if (y <= 10) {
        /* Compute gamma(x) for -10 <= x <= 10
         * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
         * first of all. 
         */

        n = Math.trunc(x); // on avait un typecast n = (int) x;

        if (x < 0) {
            --n;
        }

        y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
        --n;
        value = zNum.chebyshev_eval(y * 2 - 1, zNum.gamma.gamcs, ngam) + .9375;
        if (n === 0) {
            return value;/* x = 1.dddd = 1+y */
        }

        if (n < 0) {
            /* compute gamma(x) for -10 <= x < 1 */
            /* exact 0 or "-n" checked already above */
            /* The answer is less than half precision */
            /* because x too near a negative integer. */
            if ((x < -0.5) && (Math.abs(x - Math.trunc(x - 0.5) / x) < dxrel)) // if ((x < -0.5) && (Math.abs(x - (int)(x - 0.5) / x) < dxrel))
            {
                // ML_ERROR(ME_PRECISION, "gammafn");
            }

            /* The argument is so close to 0 that the result would overflow. */
            if (y < xsml) {
                // ML_ERROR(ME_RANGE, "gammafn");
                if (x > 0) {
                    return Infinity;
                }

                else {
                    return -Infinity;
                }
            }

            n = -n;
            for (i = 0; i < n; i++) {
                value /= (x + i);
            }

            return value;
        }
        else {
            /* gamma(x) for 2 <= x <= 10 */

            for (i = 1; i <= n; i++) {
                value *= (y + i);
            }

            return value;
        }
    } else {
        /* gamma(x) for	 y = |x| > 10. */

        if (x > xmax) {			/* Overflow */
            // ML_ERROR(ME_RANGE, "gammafn");
            return Infinity;
        }

        if (x < xmin) {			/* Underflow */
            // ML_ERROR(ME_UNDERFLOW, "gammafn");
            return 0.0;
        }

        if ((y <= 50.0) && (y === Math.trunc(y))) { /* compute (n - 1)! */
            value = 1.0;
            for (i = 2; i < y; i++) {
                value *= i;
            }
        }
        else {
            /* normal case */
            // pas de sens en Javascript ... ((2 * y == (int)2 * y)
            //value = Math.exp((y - 0.5) * Math.log(y) - y + zU.CONST.M_LN_SQRT_2PI + ((2 * y == (int)2 * y) ? StirlingError(y) : lgammacor(y)));
            value = Math.exp((y - 0.5) * Math.log(y) - y + zU.CONST.M_LN_SQRT_2PI + zNum.lgammacor(y));
        }

        if (x > 0) {
            return value;
        }

        if (Math.abs((x - Math.trunc(x - 0.5)) / x) < dxrel) {

            /* The answer is less than half precision because */
            /* the argument is too near a negative integer. */
            // ML_ERROR(ME_PRECISION, "gammafn");
        }

        sinpiy = zNum.sinpi(y);
        if (sinpiy === 0.0) {		/* Negative integer arg - overflow */
            // ML_ERROR(ME_RANGE, "gammafn");
            return Infinity;
        }

        return -zU.CONST.M_PI / (y * sinpiy * value);
    }
};

// R file: gamma.c
zNum.gamma.gamcs = [
    +.8571195590989331421920062399942e-2,
    +.4415381324841006757191315771652e-2,
    +.5685043681599363378632664588789e-1,
    -.4219835396418560501012500186624e-2,
    +.1326808181212460220584006796352e-2,
    -.1893024529798880432523947023886e-3,
    +.3606925327441245256578082217225e-4,
    -.6056761904460864218485548290365e-5,
    +.1055829546302283344731823509093e-5,
    -.1811967365542384048291855891166e-6,
    +.3117724964715322277790254593169e-7,
    -.5354219639019687140874081024347e-8,
    +.9193275519859588946887786825940e-9,
    -.1577941280288339761767423273953e-9,
    +.2707980622934954543266540433089e-10,
    -.4646818653825730144081661058933e-11,
    +.7973350192007419656460767175359e-12,
    -.1368078209830916025799499172309e-12,
    +.2347319486563800657233471771688e-13,
    -.4027432614949066932766570534699e-14,
    +.6910051747372100912138336975257e-15,
    -.1185584500221992907052387126192e-15,
    +.2034148542496373955201026051932e-16,
    -.3490054341717405849274012949108e-17,
    +.5987993856485305567135051066026e-18,
    -.1027378057872228074490069778431e-18,
    +.1762702816060529824942759660748e-19,
    -.3024320653735306260958772112042e-20,
    +.5188914660218397839717833550506e-21,
    -.8902770842456576692449251601066e-22,
    +.1527474068493342602274596891306e-22,
    -.2620731256187362900257328332799e-23,
    +.4496464047830538670331046570666e-24,
    -.7714712731336877911703901525333e-25,
    +.1323635453126044036486572714666e-25,
    -.2270999412942928816702313813333e-26,
    +.3896418998003991449320816639999e-27,
    -.6685198115125953327792127999999e-28,
    +.1146998663140024384347613866666e-28,
    -.1967938586345134677295103999999e-29,
    +.3376448816585338090334890666666e-30,
    -.5793070335782135784625493333333e-31];

// R file: cospi.c
zNum.sinpi = function (x) {
    if (!isFinite(x)) {
        return NaN;
    }


    // x = Math.mod(x, 2.0); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
    // en C#
    x = x % 2.0;
    // map (-2,2) --> (-1,1] :
    if (x <= -1) {
        x += 2.0;
    } else if (x > 1.0) {
        x -= 2.0;
    }

    if ((x === 0.0) || (x === 1.0)) {
        return 0.0;
    }

    if (x === 0.5) {
        return 1.0;
    }

    if (x === -0.5) {
        return -1.0;
    }
    // otherwise
    return Math.sin(zU.CONST.M_PI * x);
};

/* R file: lgamma.c
 */
zNum.logGamma.xmax = 2.5327372760800758e+305;

zNum.logGamma.dxrel = 1.490116119384765625e-8;

zNum.logGamma.lgammafn_sign = function (x, sgn) {
    var /*double*/ ans, y, sinpiy;
    var xmax = zNum.logGamma.xmax; // 2.5327372760800758e+305,
    var dxrel = zNum.logGamma.dxrel; // 1.490116119384765625e-8;

    // attention modif. voulue du code original
    if (sgn !== null) sgn = 1;

    if (isNaN(x)) {
        return x;
    }

    if ((sgn !== null) && (x < 0) && ((Math.floor(-x) % 2.0) === 0)) { // fmod était utilisé : fmod(floor(-x), 2.) == 0
        sgn = -1;
    }

    if ((x <= 0) && (x === Math.trunc(x))) { /* Negative integer argument */
        return Infinity;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = Math.abs(x);

    if (y < 1e-306) {
        return -Math.log(y); // denormalized range, R change
    }

    if (y <= 10) {
        return Math.log(Math.abs(zNum.gamma(x)));
    }
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
        // ML_ERROR(ME_RANGE, "lgamma");
        return Infinity;
    }

    if (x > 0) { /* i.e. y = x > 10 */
        if (x > 1e17) {
            return (x * (Math.log(x) - 1.0));
        }
        else if (x > 4934720.0) {
            return (zU.CONST.M_LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x);
        }
        else {
            return zU.CONST.M_LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x + zNum.lgammacor(x);
        }
    }
    /* else: x < -10; y = -x */
    sinpiy = Math.abs(zNum.sinpi(y));

    if (sinpiy === 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
        // MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
        return NaN;
        // ML_ERR_return_NAN;
    }

    ans = zU.CONST.M_LN_SQRT_PId2 + (x - 0.5) * Math.log(y) - x - Math.log(sinpiy) - zNum.lgammacor(y);

    if (Math.abs((x - Math.trunc(x - 0.5)) * ans / x) < dxrel) {

        /* The answer is less than half precision because
         * the argument is too near a negative integer. */

        // ML_ERROR(ME_PRECISION, "lgamma");
    }

    return ans;
}; //zNum.logGamma.lgammafn_sign

// R file: lgammacor.c
zNum.lgammacor = function (x) {

    var tmp;

    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     *   xbig = 2 ^ 26.5
     *   xmax = DBL_MAX / 48 =  2^1020 / 3 
     */

    var /*int*/ nalgm = 5;
    var /*double*/
        xbig = 94906265.62425156,
        xmax = 3.745194030963158e306;

    if (x < 10) {
        return NaN;
    } else if (x >= xmax) {
        // ML_ERROR(ME_UNDERFLOW, "lgammacor");
        /* allow to underflow below */
    } else if (x < xbig) {
        tmp = 10 / x;
        return zNum.chebyshev_eval(tmp * tmp * 2 - 1, zNum.lgammacor.algmcs, nalgm) / x;
    }

    return 1 / (x * 12);
};

// R file: lgammacor.c
zNum.lgammacor.algmcs = [
    +.1666389480451863247205729650822e+0,
    -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13,
    -.3399615005417721944303330599666e-15,
    +.2683181998482698748957538846666e-17,
    -.2868042435334643284144622399999e-19,
    +.3962837061046434803679306666666e-21,
    -.6831888753985766870111999999999e-23,
    +.1429227355942498147573333333333e-24,
    -.3547598158101070547199999999999e-26,
    +.1025680058010470912000000000000e-27,
    -.3401102254316748799999999999999e-29,
    +.1276642195630062933333333333333e-30
];

zygotine.Num.digamma = (function () {

    //    pour digamma_imp voir
    //    https://github.com/compute-io/digamma/blob/master/lib/number.js 
    //    pour polyval voir
    //    https://github.com/compute-io/polynomial/blob/master/lib/index.js

    //  Adapted from the Boost library implementation:
    //  (C) Copyright John Maddock 2006.
    //  Use, modification and distribution are subject to the
    //  Boost Software License, Version 1.0. (See accompanying file
    //  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)



    var
        floor = Math.floor,
        log = Math.log,
        tan = Math.tan;


    function polyval(c, x) {
        var len = c.length,
            p = 0,
            i = 0;
        for (; i < len; i++) {
            p = p * x + c[i];
        }

        return p;
    }

    /**
     * FUNCTION: digamma_imp_large( x )
     *	Evaluate digamma function via asymptotic expansion.
     *	This gives 17-digit precision for x >= 10:
     *
     * @param {Number} x - input value
     * @returns {Number} function value
     */
    function digamma_imp_large(x) {
        var P,
            result,
            z;
        P = [
            -0.44325980392156862745098039215686274509803921568627,
            0.083333333333333333333333333333333333333333333333333,
            -0.021092796092796092796092796092796092796092796092796,
            0.0075757575757575757575757575757575757575757575757576,
            -0.0041666666666666666666666666666666666666666666666667,
            0.003968253968253968253968253968253968253968253968254,
            -0.0083333333333333333333333333333333333333333333333333,
            0.083333333333333333333333333333333333333333333333333
        ];
        x -= 1;
        result = log(x);
        result += 1 / (2 * x);
        z = 1 / (x * x);
        result -= z * polyval(P, z);
        return result;
    } // end FUNCTION digamma_imp_large()

    /**
* FUNCTION: digamma_imp_1_2( x )
*	Evaluates digamma function over interval [1,2]. This gives 17-digit precision.
*
* @param {Number} x - input value
* @returns {Number} function value
*/
    function digamma_imp_1_2(x) {
        /*
            Now the approximation, we use the form:
    
            digamma(x) = (x - root) * (Y + R(x-1))
    
            Where root is the location of the positive root of digamma,
            Y is a constant, and R is optimised for low absolute error
            compared to Y.
    
            Maximum Deviation Found:               1.466e-18
            At double precision, max error found:  2.452e-17
        */

        var Y = 0.99558162689208984,
            root1 = 1569415565 / 1073741824,
            root2 = (381566830 / 1073741824) / 1073741824,
            root3 = 0.9016312093258695918615325266959189453125e-19,
            P,
            Q,
            g,
            r,
            result;

        P = [
            -0.0020713321167745952,
            -0.045251321448739056,
            -0.28919126444774784,
            -0.65031853770896507,
            -0.32555031186804491,
            0.25479851061131551
        ];
        Q = [
            -0.55789841321675513e-6,
            0.0021284987017821144,
            0.054151797245674225,
            0.43593529692665969,
            1.4606242909763515,
            2.0767117023730469,
            1.0
        ];
        g = x - root1;
        g -= root2;
        g -= root3;
        r = polyval(P, x - 1) / polyval(Q, x - 1);
        result = g * Y + g * r;
        return result;
    } // end FUNCTION digamma_imp_1_2()


    return function (x) {


        var
            result = 0,
            remainder;

        // CONSTANTS //
        var
            MIN_SAFE_ASYMPTOTIC = 10,
            PI = Math.PI;

        // Check for negative arguments and use reflection:
        if (x <= -1) {
            // Reflect:
            x = 1 - x;
            // Argument reduction for tan:
            remainder = x - floor(x);
            // Shift to negative if > 0.5:
            if (remainder > 0.5) {
                remainder -= 1;
            }
            // check for evaluation at a negative pole:
            if (remainder === 0) {
                return NaN;
            }

            result = PI / tan(PI * remainder);
        }

        if (x === 0) {
            return NaN;
        }
        //
        // If we're above the lower-limit for the
        // asymptotic expansion then use it:
        //
        if (x >= MIN_SAFE_ASYMPTOTIC) {
            result += digamma_imp_large(x);
        } else {
            // If x > 2 reduce to the interval [1,2]:
            while (x > 2) {
                x -= 1;
                result += 1 / x;
            }

            // If x < 1 use recurrance to shift to > 1:
            while (x < 1) {
                result -= 1 / x;
                x += 1;
            }

            result += digamma_imp_1_2(x);
        }

        return result;


    };

})();


zNum.expm1 = function (x) {
    // R file: expm1.c
    var
        y, a = Math.abs(x);
    if (a < zU.CONST.DBL_EPSILON) {
        return x;
    }


    if (a > 0.697) {
        return Math.exp(x) - 1;  /* negligible cancellation */
    }

    if (a > 1e-8) {
        y = Math.exp(x) - 1;
    }
    else { /* Taylor expansion, more accurate in this range */
        y = (x / 2 + 1) * x;
    }

    /* Newton step for solving   log(1 + y) = x   for y : */
    /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
    y -= (1 + y) * (zNum.log1p(y) - x);
    return y;
};

zNum.log1p = (function (x) {
    // R file: log1p.c
    /* series for log1p on the interval -.375 to .375
     *				     with weighted error   6.35e-32
     *				      log weighted error  31.20
     *			    significant figures required  30.93
     *				 decimal places required  32.01
     */
    var alnrcs = [
           +0.10378693562743769800686267719098e+1,
            -0.13364301504908918098766041553133e+0,
            +0.19408249135520563357926199374750e-1,
            -0.30107551127535777690376537776592e-2,
            +0.48694614797154850090456366509137e-3,
            -0.81054881893175356066809943008622e-4,
            +0.13778847799559524782938251496059e-4,
            -0.23802210894358970251369992914935e-5,
            +0.41640416213865183476391859901989e-6,
            -0.73595828378075994984266837031998e-7,
            +0.13117611876241674949152294345011e-7,
            -0.23546709317742425136696092330175e-8,
            +0.42522773276034997775638052962567e-9,
            -0.77190894134840796826108107493300e-10,
            +0.14075746481359069909215356472191e-10,
            -0.25769072058024680627537078627584e-11,
            +0.47342406666294421849154395005938e-12,
            -0.87249012674742641745301263292675e-13,
            +0.16124614902740551465739833119115e-13,
            -0.29875652015665773006710792416815e-14,
            +0.55480701209082887983041321697279e-15,
            -0.10324619158271569595141333961932e-15,
            +0.19250239203049851177878503244868e-16,
            -0.35955073465265150011189707844266e-17,
            +0.67264542537876857892194574226773e-18,
            -0.12602624168735219252082425637546e-18,
            +0.23644884408606210044916158955519e-19,
            -0.44419377050807936898878389179733e-20,
            +0.83546594464034259016241293994666e-21,
            -0.15731559416479562574899253521066e-21,
            +0.29653128740247422686154369706666e-22,
            -0.55949583481815947292156013226666e-23,
            +0.10566354268835681048187284138666e-23,
            -0.19972483680670204548314999466666e-24,
            +0.37782977818839361421049855999999e-25,
            -0.71531586889081740345038165333333e-26,
            +0.13552488463674213646502024533333e-26,
            -0.25694673048487567430079829333333e-27,
            +0.48747756066216949076459519999999e-28,
            -0.92542112530849715321132373333333e-29,
            +0.17578597841760239233269760000000e-29,
            -0.33410026677731010351377066666666e-30,
            +0.63533936180236187354180266666666e-31
    ];

    var fn = function (x) {
        var nlnrel = 22;
        var xmin = -0.999999985;

        //if (x == 0.0) return 0.0;/* speed */
        //if (x == -1) return(ML_NEGINF);
        //if (x  < -1) ML_ERR_return_NAN;

        if (Math.abs(x) <= .375) {
            /* Improve on speed (only);
           again give result accurate to IEEE double precision: */
            if (Math.abs(x) < .5 * zU.CONST.DBL_EPSILON) {
                return x;
            }

            if ((0 < x && x < 1e-8) || (-1e-9 < x && x < 0)) {
                return x * (1 - .5 * x);
            }
            /* else */
            return x * (1 - x * zNum.chebyshev_eval(x / .375, alnrcs, nlnrel));
        }
        /* else */
        if (x < xmin) {
            /* answer less than half precision because x too near -1 */
            throw new Exception("log1p: ME_PRECISION");
        }

        return Math.log(1 + x);
    };

    return fn;
})();




