/*****************************************************************************
 * Project: Momentum fraction pdf for mis ID invariant mass correction       *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 01-2016                                                            *
 *****************************************************************************/

#include "RooFit.h"
#include "Riostream.h" 

#include <cassert>
#include <cmath>
#include "TIterator.h"
#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooArgList.h"
#include "RooRecursiveFraction.h"
#include "RooRealVar.h"
#include "RooMsgService.h"

#include "RooPowerLawApproxVar.h"
#include "RooMomentumFractionPdf.h" 

// constructors
RooMomentumFractionPdf::RooMomentumFractionPdf(const char *name, const char *title, 
            RooAbsReal& _x,
            RooAbsReal& _xmin,
            RooAbsReal& _dx,
            RooAbsReal& _power) :
        RooAbsPdf(name,title),
        m_x("x","x",this,_x),
        m_xmin("xmin","xmin",this,_xmin),
        m_dx("dx","dx",this,_dx),
        m_power("power","power",this,_power),
        m_slopelist("slopelist","slopelist",this,kTRUE,kFALSE),
        m_coefslist("coefslist","coefslist",this,kTRUE,kFALSE),
        m_slopeSize(5)
{
    // approximation of a power law with exponentials
    Double_t beta_approx = 3.0;
    RooPowerLawApproxVar::Flags flag_slope = RooPowerLawApproxVar::Slope;
    RooPowerLawApproxVar::Flags flag_coeff = RooPowerLawApproxVar::Coeff;
    for (Int_t i = 0; i < m_slopeSize; i++) {
        RooAbsReal* coeff = new RooPowerLawApproxVar(Form("%s_power_law_approx_coeff_%i",GetName(),i),"Approx Power Law Coeff", 
                _power, i, beta_approx, m_slopeSize, flag_coeff);
        RooAbsReal* slope = new RooPowerLawApproxVar(Form("%s_power_law_approx_slope_%i",GetName(),i),"Approx Power Law Slope", 
                _power, i, beta_approx, m_slopeSize, flag_slope);
        m_coefslist.add(*coeff);
        m_slopelist.add(*slope);
    }
}

RooMomentumFractionPdf::RooMomentumFractionPdf(const char *name, const char *title, 
            RooAbsReal& _x,
            RooAbsReal& _xmin,
            RooAbsReal& _dx,
            RooArgList& _slopelist,
            RooArgList& _coefslist) :
        RooAbsPdf(name,title),
        m_x("x","x",this,_x),
        m_xmin("xmin","xmin",this,_xmin),
        m_dx("dx","dx",this,_dx),
        m_slopelist("slopelist","slopelist",this,kTRUE,kFALSE),
        m_coefslist("coefslist","coefslist",this,kTRUE,kFALSE),
        m_slopeSize(_slopelist.getSize())
{ 
    // assert the (relative) sizes of the lists
    if (!((m_slopeSize == _coefslist.getSize()) || (m_slopeSize == _coefslist.getSize() + 1))) {
        coutE(InputArguments) << "RooMomentumFractionPdf::RooMomentumFractionPdf(" << GetName() << ") number of slopes and coefficients inconsistent" << std::endl;
        coutE(InputArguments) << "... use N_slopes = N_coefs or use N_slopes = N_coefs + 1 for recursive fraction (like RooAddPdf)" << std::endl;
        assert(0);
    }

    // iterate over the lists
    TIterator* m_slopeIter = _slopelist.createIterator();
    TIterator* m_coefsIter = _coefslist.createIterator();
    RooAbsReal* slope;
    RooAbsReal* coeff;
    RooArgList partinCoefList;
    Bool_t first(kTRUE);

    if (_coefslist.getSize() == 0) {
        RooRealVar* onecoeff = new RooRealVar(Form("coeff_one_%s",GetName()),"coeff one",1.0);
        addOwnedComponents(*onecoeff);
        m_coefslist.add(*onecoeff);
        m_slopelist.add(*((RooAbsReal*)m_slopeIter->Next()));
    }
    else {
        while ((slope = (RooAbsReal*)m_slopeIter->Next())) {
            coeff = (RooAbsReal*)m_coefsIter->Next();
            m_slopelist.add(*slope);
            if ( m_slopeSize == _coefslist.getSize() + 1 ) {
                partinCoefList.add(*coeff);
                if ( first ) {
                    m_coefslist.add(*coeff);
                    first = kFALSE;
                } 
                else {
                    RooAbsReal* rfrac = new RooRecursiveFraction(Form("%s_recursive_fraction_%s",
                                GetName(),slope->GetName()),"Recursive Fraction",partinCoefList);
                    addOwnedComponents(*rfrac);
                    m_coefslist.add(*rfrac);
                }
            }
            else {
                m_coefslist.add(*coeff);
            }
        }
        if (_coefslist.getSize() == 0) {
            RooRealVar* onecoeff = new RooRealVar(Form("coeff_one_%s",GetName()),"coeff one",1.0);
            addOwnedComponents(*onecoeff);
            m_coefslist.add(*onecoeff);
        }
    }
}

RooMomentumFractionPdf::RooMomentumFractionPdf(const RooMomentumFractionPdf& other, const char* name) :  
    RooAbsPdf(other,name), 
    m_x("x",this,other.m_x),
    m_xmin("xmin",this,other.m_xmin),
    m_dx("dx",this,other.m_dx),
    m_power("power",this,other.m_power),
    m_slopelist("slopelist",this,other.m_slopelist),
    m_coefslist("coefslist",this,other.m_coefslist),
    m_slopeSize(other.m_slopeSize)
{ 
} 

// private functions
Double_t RooMomentumFractionPdf::getValExponentials(const double xval) const 
{
    Double_t ret = 0.0;
    Double_t c[m_slopeSize],s[m_slopeSize];
    for (int i = 0; i < m_slopeSize; ++i) {
        c[i] = ((RooAbsReal*)m_coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)m_slopelist.at(i))->getVal();
        ret += c[i] * exp(s[i] * xval);
    }
    return ret;
}

Double_t RooMomentumFractionPdf::getValSpline(const double xval) const 
{
    // get range vars
    Double_t dxval = m_dx;
    Double_t xminval = m_xmin;

    // retrieve exponent function and build spline
    Double_t c[m_slopeSize],s[m_slopeSize];
    Double_t y1 = 0.0;
    Double_t y2 = 0.0;
    Double_t  z = 0.0;
    Double_t k2 = 0.0;
    Double_t val;
    for (int i = 0; i < m_slopeSize; ++i) {
        c[i] = ((RooAbsReal*)m_coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)m_slopelist.at(i))->getVal();
        val = c[i] * exp(s[i]*(xminval+dxval));
        y2 += val;
        k2 += s[i] * val;
        z  += s[i] * s[i] * val;
    }
    Double_t t  = (xval - xminval) / dxval;
    Double_t b  = -k2 * dxval + (y2-y1);
    Double_t a  = 2.0 * b + 0.5 * z * dxval * dxval ;
    return (1.0 - t) * y1 + t * y2 + t * (1.0 - t) * (a * (1.0 - t) + b * t);
}

Double_t RooMomentumFractionPdf::getIntExponentials(const double x0, const double x1) const 
{
    // integrate over exponential part
    // sum over exponent terms
    Double_t ret = 0.0;
    Double_t c[m_slopeSize],s[m_slopeSize];
    for (int i = 0; i < m_slopeSize; ++i) {
        c[i] = ((RooAbsReal*)m_coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)m_slopelist.at(i))->getVal();
        if ( fabs(s[i]) > 0.0 ) { 
            ret += c[i] * (exp(s[i] * x1) - exp(s[i] * x0)) / s[i];
        } else {
            ret += c[i] * (x1 - x0);
        }
    }
    return ret;
}

Double_t RooMomentumFractionPdf::getIntSpline(const double x0, const double x1) const {
    // get range vars
    Double_t dxval = m_dx;
    Double_t xminval = m_xmin;

    // retrieve exponent function and build spline
    Double_t c[m_slopeSize],s[m_slopeSize];
    Double_t y1 = 0.0;
    Double_t y2 = 0.0;
    Double_t  z = 0.0;
    Double_t k2 = 0.0;
    Double_t val;
    for (int i = 0; i < m_slopeSize; ++i) {
        c[i] = ((RooAbsReal*)m_coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)m_slopelist.at(i))->getVal();
        val = c[i] * exp(s[i]*(xminval+dxval));
        y2 += val;
        k2 += s[i] * val;
        z  += s[i] * s[i] * val;
    }
    Double_t tmin = (x0 - xminval) / dxval;
    Double_t tmax = (x1 - xminval) / dxval;
    Double_t b  = -k2 * dxval + (y2-y1);
    Double_t a  = 2.0 * b + 0.5 * z * dxval * dxval ;
    Double_t t1 = tmax - tmin;
    Double_t t2 = pow(tmax,2) - pow(tmin,2);
    Double_t t3 = pow(tmax,3) - pow(tmin,3);
    Double_t t4 = pow(tmax,4) - pow(tmin,4);
    return dxval * ( t4 * (a-b)/4.0 + t3 * (b-2.0*a)/3.0 + t2 * (a-y1+y2)/2.0 + t1 * y1 );
}

// public functions

Double_t RooMomentumFractionPdf::evaluate() const 
{
    Double_t ret = 0.0;
    Double_t xval = m_x;
    Double_t dxval = m_dx;
    Double_t xminval = m_xmin;
    // split into: outer / spline like turnon / exponentials 
    if ( xval < xminval ) {
        // outer region where pdf = 0
        return ret;
    } 
    else if ( xval >= xminval + dxval ) {
        // sum over exponent terms
        ret = getValExponentials(xval);
        return ret;
    }
    else {
        ret = getValSpline(xval);        
        return ret;
    }
} 

Int_t RooMomentumFractionPdf::getAnalyticalIntegral(RooArgSet& allVars,
        RooArgSet& analVars, const char* /*rangeName*/) const
{
    if (matchArgs(allVars, analVars, m_x)) return 1;
    return 0;
}

Double_t RooMomentumFractionPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
    // assert one has to integrate over x
    assert(1 == code);
    Double_t ret = 0.0;

    // integration ranges
    Double_t x0 = m_x.min(rangeName);
    Double_t x1 = m_x.max(rangeName);
    
    Double_t xminval = m_xmin;
    Double_t dxval   = m_dx;
    Double_t xmidval = xminval + dxval;
    
    // split in different region options
    if ( x0 <= xminval ) {
        if ( x1 <= xminval ) { 
            // integrate over lower region
            return ret; 
        }
        else if ( x1 <= xmidval ) {
            // integrate over spline like function and lower region
            ret = getIntSpline(xminval,x1);
            return ret;   
        }
        else {
            // integrate over entire range
            if ( fabs(dxval) > 0.0 ) {
                // integrate over spline and exponentials
                ret = getIntSpline(xminval,xmidval) + getIntExponentials(xmidval,x1);
                return ret;
            }
            else {
                // integrate over exponentials (no spline, since dx = 0)
                ret = getIntExponentials(xminval,x1);
                return ret;
            }
        }
    } 
    else if ( dxval > 0.0 && x0 <= xmidval ) { 
        if ( x1 <= xmidval && x0 > xminval ) {
            // integrate over spline like function only
            ret = getIntSpline(x0,x1);
            return ret;
        }
        else {
            // integrate over spline and exponential
            ret = getIntSpline(x0,xmidval) + getIntExponentials(xmidval,x1);
            return ret;
        }
    }
    else {
        // integrate over exponential part
        ret = getIntExponentials(x0,x1);
        return ret;
    }
}

RooArgList* RooMomentumFractionPdf::parseDependents() const 
{
    // parse a list of dependents for Gaussian MisID function
    RooArgList *parseList = new RooArgList("dependents");
    
    // add parameters (in order)
    parseList -> add(m_x.arg());
    parseList -> add(m_xmin.arg());
    parseList -> add(m_dx.arg());
    for (int i = 0; i < m_coefslist.getSize(); ++i) {
        parseList -> add(m_coefslist[i]);    
    }
    for (int i = 0; i < m_slopelist.getSize(); ++i) {
        parseList -> add(m_slopelist[i]);
    }
    return parseList;
}
