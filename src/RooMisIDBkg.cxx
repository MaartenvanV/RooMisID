/*****************************************************************************
 * Project: invariant mass pdf for mis ID invariant mass                     *
 *          correction without resolution                                    *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 08-2018                                                            *
 *****************************************************************************/

#include "RooFit.h"
#include "Riostream.h" 

#include <cassert>
#include <cmath>
#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsProxy.h"
#include "RooRecursiveFraction.h"
#include "RooMath.h"
#include "RooRealVar.h"
#include "RooMsgService.h"

#include "RooPowerLawApproxVar.h"
#include "RooMomentumFractionPdf.h"
#include "RooMisIDBkg.h" 

// constructors
RooMisIDBkg::RooMisIDBkg(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _dm2,
            RooAbsReal& _power,
            RooAbsReal& _xmin,
            RooAbsReal& _xmax,
            RooAbsReal& _dx) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        dm2("dm2","dm2",this,_dm2),
        power("power","power",this,_power),
        xmin("xmin","xmin",this,_xmin),
        xmax("xmax","xmax",this,_xmax),
        dx("dx","dx",this,_dx),
        slopelist("slopelist","slopelist",this,kTRUE,kFALSE),
        coefslist("coefslist","coefslist",this,kTRUE,kFALSE),
        slopeSize(5)
{
    // approximation of a power law with exponentials
    Double_t beta_approx = 3.0;
    RooPowerLawApproxVar::Flags flag_slope = RooPowerLawApproxVar::Slope;
    RooPowerLawApproxVar::Flags flag_coeff = RooPowerLawApproxVar::Coeff;
    for (Int_t i = 0; i < slopeSize; i++) {
        RooAbsReal* coeff = new RooPowerLawApproxVar(Form("%s_power_law_approx_coeff_%i",GetName(),i),"Approx Power Law Coeff",
                                    _power, i, beta_approx, slopeSize, flag_coeff);
        RooAbsReal* slope = new RooPowerLawApproxVar(Form("%s_power_law_approx_slope_%i",GetName(),i),"Approx Power Law Slope", 
                                    _power, i, beta_approx, slopeSize, flag_slope);
        coefslist.add(*coeff);
        slopelist.add(*slope);
    }
}

RooMisIDBkg::RooMisIDBkg(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _dm2,
            RooAbsReal& _xmin,
            RooAbsReal& _xmax,
            RooAbsReal& _dx,
            RooArgList& _slopelist,
            RooArgList& _coefslist) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        dm2("dm2","dm2",this,_dm2),
        xmin("xmin","xmin",this,_xmin),
        dx("dx","dx",this,_dx),
        slopelist("slopelist","slopelist",this,kTRUE,kFALSE),
        xmax("xmax","xmax",this,_xmax),
        coefslist("coefslist","coefslist",this,kTRUE,kFALSE),
        slopeSize(_slopelist.getSize())
{
    // assert the (relative) sizes of the lists
    if (!((slopeSize == _coefslist.getSize()) || (slopeSize == _coefslist.getSize() + 1))) {
        coutE(InputArguments) << "RooMomentumFractionPdf::RooMomentumFractionPdf(" << GetName() << ") number of slopes and coefficients inconsistent" << std::endl;
        coutE(InputArguments) << "... use N_slopes = N_coefs or use N_slopes = N_coefs + 1 for recursive fraction (like RooAddPdf)" << std::endl;
        assert(0);
    }

    // iterate over the lists
    TIterator* slopeIter = _slopelist.createIterator();
    TIterator* coefsIter = _coefslist.createIterator();
    RooAbsReal* slope;
    RooAbsReal* coeff;
    RooArgList partinCoefList;
    Bool_t first(kTRUE);

    if (_coefslist.getSize() == 0) {
        RooRealVar* onecoeff = new RooRealVar(Form("coeff_one_%s",GetName()),"coeff one",1.0);
        addOwnedComponents(*onecoeff);
        coefslist.add(*onecoeff);
        slopelist.add(*((RooAbsReal*)slopeIter->Next()));
    }
    else {
        while ((slope = (RooAbsReal*)slopeIter->Next())) {
            coeff = (RooAbsReal*)coefsIter->Next();
            slopelist.add(*slope);
            if ( slopeSize == _coefslist.getSize() + 1 ) {
                partinCoefList.add(*coeff);
                if ( first ) {
                    coefslist.add(*coeff);
                    first = kFALSE;
                } 
                else {
                    RooAbsReal* rfrac = new RooRecursiveFraction(Form("%s_recursive_fraction_%s",
                                GetName(),slope->GetName()),"Recursive Fraction",partinCoefList);
                    addOwnedComponents(*rfrac);
                    coefslist.add(*rfrac);
                }
            }
            else {
                coefslist.add(*coeff);
            }
        }
        if (_coefslist.getSize() == 0) {
            RooRealVar* onecoeff = new RooRealVar(Form("coeff_one_%s",GetName()),"coeff one",1.0);
            addOwnedComponents(*onecoeff);
            coefslist.add(*onecoeff);
        }
    }
}

RooMisIDBkg::RooMisIDBkg(const RooMisIDBkg& other, const char* name) :  
    RooAbsPdf(other,name),
    m("m",this,other.m),
    mean("mean",this,other.mean),
    dm2("dm2",this,other.dm2),
    x("x",this,other.x),
    xmin("xmin",this,other.xmin),
    xmax("xmax",this,other.xmax),
    dx("dx",this,other.dx),
    power("power",this,other.power),
    slopelist("slopelist",this,other.slopelist),
    coefslist("coefslist",this,other.coefslist),
    slopeSize(other.slopeSize)
{
} 

// functions
Double_t RooMisIDBkg::getValExponentials(const double xval) const 
{
    Double_t ret = 0.0;
    Double_t c[slopeSize],s[slopeSize];
    for (int i = 0; i < slopeSize; ++i) {
        c[i] = ((RooAbsReal*)coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)slopelist.at(i))->getVal();
        ret += c[i] * exp(s[i] * xval);
    }
    return ret;
}

Double_t RooMisIDBkg::getValSpline(const double xval) const 
{
    // retrieve exponent function and build spline
    Double_t c[slopeSize],s[slopeSize];
    Double_t y1 = 0.0;
    Double_t y2 = 0.0;
    Double_t  z = 0.0;
    Double_t k2 = 0.0;
    Double_t val;
    for (int i = 0; i < slopeSize; ++i) {
        c[i] = ((RooAbsReal*)coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)slopelist.at(i))->getVal();
        val = c[i] * exp(s[i]*(xmin+dx));
        y2 += val;
        k2 += s[i] * val;
        z  += s[i] * s[i] * val;
    }
    Double_t t  = (xval - xmin) / dx;
    Double_t b  = -k2 * dx + (y2-y1);
    Double_t a  = 2.0 * b + 0.5 * z * dx * dx ;
    return (1.0 - t) * y1 + t * y2 + t * (1.0 - t) * (a * (1.0 - t) + b * t);
}

Double_t RooMisIDBkg::evaluate() const 
{
    // start with basic return value
    Double_t ret = 0.0;

    // transform mass to x
    Double_t xval = ( pow(m,2) - pow(mean,2) ) / dm2;

    // split into outer/turnon/exp
    if ( xval <= xmin ) {
        // outer region where pdf = 0
        ret = 0.0;
        return ret;
    } 
    else if ( xval >= xmin + dx ) {
        // sum over exponent terms
        ret = getValExponentials(xval);
        return ret;
    }
    else if ( xval < xmin + dx && xval > xmin ) {
        // spline part
        ret = getValSpline(xval);        
        return ret;
    }
    else {
        // must not go here
        assert(1==0);
    }   
    // multiply jacobian
    ret *= abs(m/dm2);
    return ret;
}

Int_t RooMisIDBkg::getAnalyticalIntegral(RooArgSet& allVars, 
        RooArgSet& analVars, const char* /*rangeName*/) const
{
    return 0;
}

Double_t RooMisIDBkg::analyticalIntegral(Int_t code, const char* rangeName) const
{
    assert(1 == code);
    return 1.0;
}


