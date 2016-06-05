/*****************************************************************************
 * Project: invariant mass pdf for mis ID invariant mass correction          *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 01-2016                                                            *
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
#include "TVectorD.h"
#include "RooRecursiveFraction.h"
#include "RooConstVar.h"
#include "RooMath.h"
#include "RooRealVar.h"
#include "RooMsgService.h"

#include "RooPowerLawApproxVar.h"
#include "RooMomentumFractionPdf.h"
#include "RooMisIDGaussian.h" 

// much used constants
static const std::complex<double> iComplex(0,1);
static const double sqrtPiOver2 = 1.2533141373155001;
static const double sqrt2 = 1.4142135623730951;
static const double sqrt2pi = 2.5066282746310002;

// constructors
RooMisIDGaussian::RooMisIDGaussian(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _sigma,
            RooAbsReal& _dm2,
            RooMomentumFractionPdf& _xpdf) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        sigma("sigma","sigma",this,_sigma),
        dm2("dm2","dm2",this,_dm2),
        x("x","x",this,kTRUE,kFALSE),
        xmin("xmin","xmin",this,kTRUE,kFALSE),
        dx("dx","dx",this,kTRUE,kFALSE),
        slopelist("slopelist","slopelist",this,kTRUE,kFALSE),
        coefslist("coefslist","coefslist",this,kTRUE,kFALSE),
        slopeSize(0)
{
    // add coefficients and slopes from referenced pdf 
    RooArgList *deps = (RooArgList*)_xpdf.parseDependents();
    x.setArg((RooAbsReal&)*(deps->at(0)));
    xmin.setArg((RooAbsReal&)*(deps->at(1)));
    dx.setArg((RooAbsReal&)*(deps->at(2)));
    int nSlopes = ((deps -> getSize()) - 3)/2;
    for (int i = 3; i < nSlopes+3; ++i) {
        coefslist.add(*(deps->at(i)));
    }
    for (int i = nSlopes+3; i < 2*nSlopes+3; ++i) {
        slopelist.add(*(deps->at(i)));
    }
    slopeSize = nSlopes;
    // ensure that functions know we deal with 2D function
    flag_2D = true;
    // use analytical integrals as standard
    flag_num_int = false;
    // other flags
    flag_abspdf = false;
    flag_conditional_x = false;
}

RooMisIDGaussian::RooMisIDGaussian(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _sigma,
            RooAbsReal& _dm2,
            RooAbsReal& _x,
            RooAbsPdf& _xpdf) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        sigma("sigma","sigma",this,_sigma),
        dm2("dm2","dm2",this,_dm2),
        x("x","x",this,_x),
        abspdf("abspdf","abspdf",this,(RooAbsReal&)_xpdf)
{
    // ensure that functions know we deal with 2D function
    flag_2D = true;
    // use numerical integrals 
    flag_num_int = false;
    // other flags
    flag_abspdf = true;
    flag_conditional_x = false;
}

RooMisIDGaussian::RooMisIDGaussian(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _sigma,
            RooAbsReal& _dm2,
            RooAbsReal& _x) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        sigma("sigma","sigma",this,_sigma),
        dm2("dm2","dm2",this,_dm2),
        x("x","x",this,_x)
{
    // ensure that functions know we deal with 2D function
    flag_2D = true;
    // use numerical integrals 
    flag_num_int = false;
    // other flags
    flag_abspdf = false;
    flag_conditional_x = true;
}

RooMisIDGaussian::RooMisIDGaussian(const char *name, const char *title, 
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _sigma,
            RooAbsReal& _dm2,
            RooAbsReal& _x,
            RooAbsReal& _xmin,
            RooAbsReal& _dx,
            RooArgList& _slopelist,
            RooArgList& _coefslist) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        sigma("sigma","sigma",this,_sigma),
        dm2("dm2","dm2",this,_dm2),
        x("x","x",this,_x),
        xmin("xmin","xmin",this,_xmin),
        dx("dx","dx",this,_dx),
        slopelist("slopelist","slopelist",this,kTRUE,kFALSE),
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
    
    // ensure functions that they know to use 2D functions
    flag_2D = true;
    // use analytical integrals as standard
    flag_num_int = false;
    // other flags
    flag_abspdf = false;
    flag_conditional_x = false;
} 

RooMisIDGaussian::RooMisIDGaussian(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _sigma,
            RooAbsReal& _dm2,
            RooAbsReal& _power,
            RooAbsReal& _xmin,
            RooAbsReal& _xmax,
            RooAbsReal& _dx) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        sigma("sigma","sigma",this,_sigma),
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

    // ensure that functions know we deal with a 1D function
    flag_2D = false;
    // use analytical integrals as standard
    flag_num_int = false;
    // other flags
    flag_abspdf = false;
    flag_conditional_x = false;
}

RooMisIDGaussian::RooMisIDGaussian(const char *name, const char *title,
            RooAbsReal& _m,
            RooAbsReal& _mean,
            RooAbsReal& _sigma,
            RooAbsReal& _dm2,
            RooAbsReal& _xmin,
            RooAbsReal& _dx,
            RooArgList& _slopelist,
            RooArgList& _coefslist,
            RooAbsReal& _xmax) :
        RooAbsPdf(name,title),
        m("m","m",this,_m),
        mean("mean","mean",this,_mean),
        sigma("sigma","sigma",this,_sigma),
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

    // ensure that functions know we deal with a 1D function
    flag_2D = false;
    // use analytical integrals as standard
    flag_num_int = false;
    // other flags
    flag_abspdf = false;
    flag_conditional_x = false;
}

RooMisIDGaussian::RooMisIDGaussian(const RooMisIDGaussian& other, const char* name) :  
    RooAbsPdf(other,name),
    m("m",this,other.m),
    mean("mean",this,other.mean),
    sigma("sigma",this,other.sigma),
    dm2("dm2",this,other.dm2),
    x("x",this,other.x),
    xmin("xmin",this,other.xmin),
    xmax("xmax",this,other.xmax),
    dx("dx",this,other.dx),
    power("power",this,other.power),
    abspdf("abspdf",this,other.abspdf),
    slopelist("slopelist",this,other.slopelist),
    coefslist("coefslist",this,other.coefslist),
    slopeSize(other.slopeSize)
{
    flag_2D = other.flag_2D; 
    flag_num_int = other.flag_num_int;
    flag_abspdf = other.flag_abspdf;
    flag_conditional_x = other.flag_conditional_x;
} 

// private functions 
Double_t RooMisIDGaussian::ApproxErf(const double arg) 
{
    static const double erflim = 5.0;
    if( arg > erflim )
        return 1.0;
    if( arg < -erflim )
        return 1.0;
    return RooMath::erf(arg);
}

Double_t RooMisIDGaussian::primitiveExpoGauss2D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double cval, const double aval) const 
{
    // primitive function that integrates over x and m for the exponential part
    Double_t retval = 0.0;
    Double_t sqrtdm2 = sqrt(dm2val);
    Double_t alpha   = dm2val - 2 * aval * sigval * sigval;
    Double_t mucorr = sqrt( meanval * meanval + dm2val * xval );
    Double_t arg = ( mval - mucorr ) /  ( sqrt2 * sigval );
    std::complex<double> sqrtAlpha = ( alpha > 0. ) ? 
        std::complex<double>(sqrt(fabs(alpha)),0.) : std::complex<double>(0.,sqrt(fabs(alpha)));
    // ensure that function works with a = 0
    if ( fabs(aval) > 0. ) {
        retval = sqrtPiOver2 * cval * sigval / aval * ( exp(aval * xval) * RooMath::erf(arg) 
                + ( sqrtdm2 * exp ( - aval * ( pow(meanval,2) / dm2val - pow(mval,2) / alpha ) ) / sqrtAlpha * 
                RooMath::erf( ( - mval * dm2val + mucorr * alpha ) / ( sqrt2 * sqrtdm2 * sigval * sqrtAlpha ) ) ).real() );
        return retval;
    }
    else { 
        retval = - cval * sigval / dm2val * ( (mval + mucorr) * sigval * exp(-arg*arg) 
               + sqrtPiOver2 * ( pow(mval,2) - dm2val * xval - pow(mucorr,2) + pow(sigval,2) ) * RooMath::erf(arg) );
        return retval;
    }    
}

Double_t RooMisIDGaussian::primitiveExpoGauss1D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double cval, const double aval) const 
{
    // primitive function that integrates over x for the exponential part
    Double_t retval = 0.0;
    Double_t sqrtdm2 = sqrt(dm2val);
    Double_t alpha   = dm2val - 2 * aval * sigval * sigval;
    Double_t mucorr = sqrt( meanval * meanval + dm2val * xval );
    Double_t arg = ( mval - mucorr ) /  ( sqrt2 * sigval );
    std::complex<double> sqrtAlpha = ( alpha > 0. ) ? 
        std::complex<double>(sqrt(fabs(alpha)),0.) : std::complex<double>(0.,sqrt(fabs(alpha)));
    // ensure that function works with a = 0
    if ( fabs(aval) > 0. ) {
        retval = cval * sigval * ( - 2 * exp( aval * xval ) * exp(-arg*arg) * sigval / alpha 
                + sqrtdm2 * mval * sqrt2pi * exp( aval * ( pow(mval,2) / alpha - pow(meanval,2) / dm2val ) ) * 
                ( RooMath::erf( ( alpha * mucorr - mval * dm2val ) / ( sqrt2 * sqrtdm2 * sigval * sqrtAlpha ) ) 
                  / ( sqrtAlpha * sqrtAlpha * sqrtAlpha ) ).real() );
        
        return retval;
    }
    else { 
        retval = - 2 * cval * sigval / dm2val * ( sigval * exp( - arg * arg ) + mval * sqrtPiOver2 * RooMath::erf( arg ) ) ;
        return retval;
    }
}

Double_t RooMisIDGaussian::primitiveSplineGauss2D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double dxval, const double xminval) const
{
    // primitive function that integrates over x and m for the spline part
    Double_t retval;
    // retrieve exponent function
    Double_t c[slopeSize],s[slopeSize];
    Double_t y1 = 0.0;
    Double_t y2 = 0.0;
    Double_t  z = 0.0;
    Double_t k2 = 0.0;
    Double_t val;
    Double_t xmidval = xminval + dxval;
    for (int i = 0; i < slopeSize; ++i) {
        c[i] = ((RooAbsReal*)coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)slopelist.at(i))->getVal();
        val = c[i] * exp(s[i]*xmidval);
        y2 += val;
        k2 += s[i] * val;
        z  += s[i] * s[i] * val;
    }
    Double_t t  = (xval - xminval) / dxval;
    Double_t b  = -k2 * dxval + (y2-y1);
    Double_t a  = 2.0 * b + 0.5 * z * dxval * dxval ;

    // convert to standard polynomial constants: s1 * x^3 + s2 * x^2 + s1 * x + s0
    Double_t s0, s1, s2, s3;
    Double_t dxval3 = pow(dxval,3);
    s0 = ( xmidval * (b * pow(xminval,2) - a * xminval * xmidval + pow(dxval,2) * y1 ) - pow(dxval,2) * xminval * y2 ) / dxval3;
    s1 = ( a * xmidval * (dxval + 3 * xminval) - b * xminval * ( 2 * dxval + 3 * xminval ) + pow(dxval,2) * (y2 - y1) ) / dxval3;
    s2 = ( - 2 * a * dxval + b * dxval - 3 * a * xminval + 3 * b * xminval ) / dxval3;
    s3 = ( a - b ) / dxval3;

    // most obvious repeating factors
    Double_t dm2xval = dm2val * xval;
    Double_t mucorr = sqrt( meanval * meanval + dm2xval );
    Double_t arg = ( mval - mucorr ) / ( sqrt2 * sigval ) ;
    Double_t erf_mx = RooMath::erf(arg);
    Double_t mvals[9], muvals[9], sigvals[9];
    for (int i = 0; i < 9; ++i) {
        mvals[i] = pow(mval,i);
        muvals[i] = pow(meanval,i);
        sigvals[i] = pow(sigval,i);
    }

    // combined polynomial primitive
    retval = (sigval*(-((sigval*(pow(dm2val,3)*(12*s0 + xval*(6*s1 + 4*s2*xval + 3*s3*pow(xval,2)))*(mval + mucorr) + 
        pow(dm2val,2)*(mvals[3]*(6*s1 + 4*s2*xval + 3*s3*pow(xval,2)) - mval*(6*s1 + 4*s2*xval + 3*s3*pow(xval,2))*muvals[2] + 
        mvals[2]*(6*s1 + 4*s2*xval + 3*s3*pow(xval,2))*mucorr + 3*mval*(10*s1 + 12*s2*xval + 13*s3*pow(xval,2))*sigvals[2] + 
        mucorr*((-6*s1 - 4*s2*xval - 3*s3*pow(xval,2))*muvals[2] + (18*s1 + 20*s2*xval + 21*s3*pow(xval,2))*sigvals[2])) + 
        dm2val*(mvals[5]*(4*s2 + 3*s3*xval) + mvals[4]*(4*s2 + 3*s3*xval)*mucorr + 
        2*mvals[2]*mucorr*((-4*s2 - 3*s3*xval)*muvals[2] + 3*(8*s2 + 9*s3*xval)*sigvals[2]) + 
        mvals[3]*(-2*(4*s2 + 3*s3*xval)*muvals[2] + 2*(28*s2 + 33*s3*xval)*sigvals[2]) + 
        mucorr*((4*s2 + 3*s3*xval)*muvals[4] - 2*(8*s2 + 9*s3*xval)*muvals[2]*sigvals[2] + 15*(4*s2 + 7*s3*xval)*sigvals[4]) + 
        mval*((4*s2 + 3*s3*xval)*muvals[4] - 6*(4*s2 + 5*s3*xval)*muvals[2]*sigvals[2] + 3*(44*s2 + 87*s3*xval)*sigvals[4])) + 
        3*s3*(mvals[7] + mvals[6]*mucorr - 3*mvals[5]*(muvals[2] - 9*sigvals[2]) + mvals[4]*mucorr*(-3*muvals[2] + 25*sigvals[2]) + 
        3*mvals[2]*mucorr*(muvals[4] - 10*muvals[2]*sigvals[2] + 47*sigvals[4]) + mvals[3]*(3*muvals[4] - 34*muvals[2]*sigvals[2] + 185*sigvals[4]) - 
        mval*(muvals[6] - 7*muvals[4]*sigvals[2] + 45*muvals[2]*sigvals[4] - 279*sigvals[6]) - 
        mucorr*(muvals[6] - 5*muvals[4]*sigvals[2] + 25*muvals[2]*sigvals[4] - 105*sigvals[6]))))/exp(arg*arg)) + 
        sqrtPiOver2*((pow(dm2,4)*xval*(12*s0 + xval*(6*s1 + 4*s2*xval + 3*s3*pow(xval,2))) + 12*pow(dm2val,3)*s0*muvals[2] - 6*pow(dm2val,2)*s1*muvals[4] + 4*dm2*s2*muvals[6] - 
        3*s3*muvals[8])*erf_mx + (mvals[2]*(12*pow(dm2val,3)*s0 + 6*pow(dm2val,2)*s1*(mvals[2] - 2*muvals[2]) + 4*dm2val*s2*(mvals[4] - 3*mvals[2]*muvals[2] + 3*muvals[4]) + 
        3*s3*(mvals[6] - 4*mvals[4]*muvals[2] + 6*mvals[2]*muvals[4] - 4*muvals[6])) + 
        12*(pow(dm2val,3)*s0 + pow(dm2val,2)*s1*(3*mvals[2] - muvals[2]) + s3*pow(mvals[2] - muvals[2],2)*(7*mvals[2] - muvals[2]) + 
        dm2val*s2*(5*mvals[4] - 6*mvals[2]*muvals[2] + muvals[4]))*sigvals[2] + 
        18*(pow(dm2val,2)*s1 + 10*dm2val*mvals[2]*s2 + 35*mvals[4]*s3 - 2*(dm2val*s2 + 15*mvals[2]*s3)*muvals[2] + 3*s3*muvals[4])*sigvals[4] + 
        60*(dm2*s2 + 21*mvals[2]*s3 - 3*s3*muvals[2])*sigvals[6] + 315*s3*sigvals[8])*-1.0*erf_mx)))/(12.*pow(dm2val,4));

    return retval;
}

Double_t RooMisIDGaussian::primitiveSplineGauss1D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double dxval, const double xminval) const
{
    // primitive function that integrates over x and m for the spline part
    Double_t retval;
    // retrieve exponent function
    Double_t c[slopeSize],s[slopeSize];
    Double_t y1 = 0.0;
    Double_t y2 = 0.0;
    Double_t  z = 0.0;
    Double_t k2 = 0.0;
    Double_t val;
    Double_t xmidval = xminval + dxval;
    for (int i = 0; i < slopeSize; ++i) {
        c[i] = ((RooAbsReal*)coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)slopelist.at(i))->getVal();
        val = c[i] * exp(s[i]*xmidval);
        y2 += val;
        k2 += s[i] * val;
        z  += s[i] * s[i] * val;
    }
    Double_t t  = (xval - xminval) / dxval;
    Double_t b  = -k2 * dxval + (y2-y1);
    Double_t a  = 2.0 * b + 0.5 * z * dxval * dxval ;

    // convert to standard polynomial constants: s1 * x^3 + s2 * x^2 + s1 * x + s0
    Double_t s0, s1, s2, s3;
    Double_t dxval3 = pow(dxval,3);
    s0 = ( xmidval * (b * pow(xminval,2) - a * xminval * xmidval + pow(dxval,2) * y1 ) - pow(dxval,2) * xminval * y2 ) / dxval3;
    s1 = ( a * xmidval * (dxval + 3 * xminval) - b * xminval * ( 2 * dxval + 3 * xminval ) + pow(dxval,2) * (y2 - y1) ) / dxval3;
    s2 = ( - 2 * a * dxval + b * dxval - 3 * a * xminval + 3 * b * xminval ) / dxval3;
    s3 = ( a - b ) / dxval3;

    // most obvious repeating factors
    Double_t dm2xval = dm2val * xval;
    Double_t mucorr = sqrt( meanval * meanval + dm2xval );
    Double_t arg = ( mval - mucorr ) / ( sqrt2 * sigval ) ;
    Double_t erf_x = RooMath::erf(arg);
    Double_t mvals[9], muvals[9], sigvals[9];
    for (int i = 0; i < 9; ++i) {
        mvals[i] = pow(mval,i);
        muvals[i] = pow(meanval,i);
        sigvals[i] = pow(sigval,i);
    }

    retval = -((sigvals[1]*(2*exp(-arg*arg)*sigvals[1]*(s3*(mvals[6] + 48*sigvals[6] - 2*mvals[4]*(muvals[2] 
             - 10*sigvals[2]) + mvals[2]*(muvals[4] - 12*muvals[2]*sigvals[2] + 87*sigvals[4]) + mvals[5]*sqrt(muvals[2] + dm2xval) 
             - 2*mvals[3]*(muvals[2] - 9*sigvals[2])*sqrt(muvals[2] + dm2xval) + mvals[1]*(muvals[4] - 10*muvals[2]*sigvals[2] 
             + 57*sigvals[4])*sqrt(muvals[2] + dm2xval)) + pow(dm2val,3)*(s0 + xval*(s1 + xval*(s2 + s3*xval))) + pow(dm2val,2)*(mvals[2]*(s1 + xval*(s2 + s3*xval)) 
             + m*sqrt(muvals[2] + dm2xval)*(s1 + xval*(s2 + s3*xval)) + 2*sigvals[2]*(s1 + xval*(2*s2 + 3*s3*xval))) + dm2val*(mvals[4]*(s2 + s3*xval) + mvals[3]*sqrt(muvals[2] 
             + dm2xval)*(s2 + s3*xval) + 8*sigvals[4]*(s2 + 3*s3*xval) - mvals[2]*(muvals[2]*(s2 + s3*xval) - 3*sigvals[2]*(3*s2 + 5*s3*xval)) 
             - m*sqrt(muvals[2] + dm2xval)*(muvals[2]*(s2 + s3*xval) - sigvals[2]*(7*s2 + 11*s3*xval)))) 
             + mvals[1]*sqrt2pi*(pow(dm2val,3)*s0 + pow(dm2val,2)*s1*(mvals[2] - muvals[2] + 3*sigvals[2]) 
             + dm2val*s2*(mvals[4] + muvals[4] - 6*muvals[2]*sigvals[2] + 15*sigvals[4] - 2*mvals[2]*(muvals[2] - 5*sigvals[2])) 
             + s3*(mvals[6] - muvals[6] + 9*muvals[4]*sigvals[2] - 45*muvals[2]*sigvals[4] + 105*sigvals[6] - 3*mvals[4]*(muvals[2] 
             - 7*sigvals[2]) + 3*mvals[2]*(muvals[4] - 10*muvals[2]*sigvals[2] + 35*sigvals[4])))*RooMath::erf(arg)))/
            pow(dm2val,4));
    
    return retval;
}

Double_t RooMisIDGaussian::getValExponentials(const double xval) const 
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

Double_t RooMisIDGaussian::getValSpline(const double xval) const 
{
    // get range vars
    Double_t dxval = dx;
    Double_t xminval = xmin;

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

// integration over x
Double_t RooMisIDGaussian::getIntExponentialsGaussOverX(const double x0, const double x1, const double mval) const
{
    Double_t retval = 0.;
    Double_t c[slopeSize],s[slopeSize];
    for (int i = 0; i < slopeSize; ++i) {
        c[i] = ((RooAbsReal*)coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)slopelist.at(i))->getVal();
        retval += + primitiveExpoGauss1D(mval,x1,mean,sigma,dm2,c[i],s[i]);
        retval += - primitiveExpoGauss1D(mval,x0,mean,sigma,dm2,c[i],s[i]);
    }
    return retval;  
}

Double_t RooMisIDGaussian::getIntSplineGaussOverX(const double x0, const double x1, const double mval) const
{
    Double_t retval = 0.;
    retval += + primitiveSplineGauss1D(mval,x1,mean,sigma,dm2,dx,xmin);
    retval += - primitiveSplineGauss1D(mval,x0,mean,sigma,dm2,dx,xmin);
    return retval;
}

// integration over m 
Double_t RooMisIDGaussian::getIntGaussNotXintegrated(const double m0, const double m1, const double xval) const
{
    // integral JUST over gaussian term
    Double_t sig = sigma;
    Double_t mucorr = sqrt( mean * mean + dm2 * xval );
    Double_t ret = sqrtPiOver2 * sig * ( RooMath::erf( ( m1 - mucorr ) / ( sqrt2 * sig ) ) 
            - RooMath::erf( ( m0 - mucorr ) / ( sqrt2 * sig ) ) );
    return ret;
}

// multidimensional integration
Double_t RooMisIDGaussian::getMultiIntExponentialsGauss(const double x0, const double x1, 
        const double m0, const double m1) const
{
    Double_t retval = 0.;
    Double_t c[slopeSize],s[slopeSize];
    for (int i = 0; i < slopeSize; ++i) {
        c[i] = ((RooAbsReal*)coefslist.at(i))->getVal();
        s[i] = ((RooAbsReal*)slopelist.at(i))->getVal();
        retval += + primitiveExpoGauss2D(m0,x0,mean,sigma,dm2,c[i],s[i]);
        retval += + primitiveExpoGauss2D(m1,x1,mean,sigma,dm2,c[i],s[i]);
        retval += - primitiveExpoGauss2D(m1,x0,mean,sigma,dm2,c[i],s[i]);
        retval += - primitiveExpoGauss2D(m0,x1,mean,sigma,dm2,c[i],s[i]);
    }
    return retval;   
}

Double_t RooMisIDGaussian::getMultiIntSplineGauss(const double x0, const double x1, 
        const double m0, const double m1) const
{
    Double_t retval = 0.;
    retval += + primitiveSplineGauss2D(m0,x0,mean,sigma,dm2,dx,xmin);
    retval += + primitiveSplineGauss2D(m1,x1,mean,sigma,dm2,dx,xmin);
    retval += - primitiveSplineGauss2D(m0,x1,mean,sigma,dm2,dx,xmin);
    retval += - primitiveSplineGauss2D(m1,x0,mean,sigma,dm2,dx,xmin);
    return retval;    
}

// public functions

Double_t RooMisIDGaussian::evaluate() const 
{
    // start with basic return value
    Double_t ret = 0.0;
    // split into two dimensional (x,m) or one dimensional (m) pdf
    if ( flag_2D ) {
        // return value depending on m and x
        Double_t mval = m;
        Double_t meanval = mean;
        Double_t sval = sigma;
        Double_t xval = x;
        Double_t dm2val = dm2; 
        
        Double_t arg = mval - sqrt(meanval * meanval + dm2val * xval);

        // base gaussian
        ret = exp(-0.5*arg*arg/(sval*sval));

        // x as conditional observable
        if ( flag_conditional_x ) return ret;

        // pdf for x
        if ( flag_abspdf ) {
            ret *= ((RooAbsPdf&)abspdf.arg()).getValV();
            return ret;
        }

        // x pdf specific variables
        Double_t dxval = dx;
        Double_t xminval = xmin;
        
        // split into outer/turnon/exp
        if ( xval <= xminval ) {
            // outer region where pdf = 0
            ret = 0.0;
            return ret;
        } 
        else if ( xval >= xminval + dxval ) {
            // sum over exponent terms
            ret *= getValExponentials(xval);
            return ret;
        }
        else if ( xval < xminval + dxval && xval > xminval ) {
            // spline part
            ret *= getValSpline(xval);        
            return ret;
        }
        else {
            // must not go here
            assert(1==0);
        }
    } else {
        // return value depending only on m
        Double_t mval = m;
        Double_t xmaxval = xmax;  
        Double_t xminval = xmin;
        Double_t dxval = dx;
        if ( dxval > 0. ) {
            if ( xminval + dxval >= xmaxval ) {
                ret = getIntSplineGaussOverX(xminval,xmaxval,mval);
            } else {
                ret = getIntSplineGaussOverX(xminval,xminval+dxval,mval) + getIntExponentialsGaussOverX(xminval+dxval,xmaxval,mval);
            }
        }
        else {
            ret = getIntExponentialsGaussOverX(xminval,xmaxval,mval);
        }
        // ensure that -0.000... is seen as +0.000... (with fabs())
        return fabs(ret);
    }
} 

Int_t RooMisIDGaussian::getAnalyticalIntegral(RooArgSet& allVars,
        RooArgSet& analVars, const char* /*rangeName*/) const
{
    if ( flag_num_int ) return 0;
    if ( flag_2D ) {
        if ( !flag_conditional_x && !flag_abspdf ) {
            if (matchArgs(allVars, analVars, x, m)) return 3;
            if (matchArgs(allVars, analVars, x)) return 1;
        }
    }
    if (matchArgs(allVars, analVars, m)) return 2;
    return 0;
}

Double_t RooMisIDGaussian::analyticalIntegral(Int_t code, const char* rangeName) const
{
    // return value
    Double_t ret = 0.0;
    
    // bool
    Bool_t multi = ( 3 == code ) || ( !flag_2D && 2 == code );

    if ( 1 == code || multi ) {
        // Ranges
        Double_t x0,x1;
        if ( flag_2D ) {
            x0 = x.min(rangeName);
            x1 = x.max(rangeName);
        }
        else {
            x0 = 0.;
            x1 = xmax; // should be limit to infinity
        }

        Double_t xminval = xmin;
        Double_t dxval   = dx;
        Double_t xmidval = xminval + dxval;

        // mass range 
        Double_t m0 = m.min(rangeName);
        Double_t m1 = m.max(rangeName);
        
        Double_t mval = m;

        // split in three regions
        if ( x0 <= xminval ) {
            if ( x1 <= xminval ) { 
                // integrate over lower region
                return ret; 
            }
            else if ( x1 <= xmidval ) {
                // integrate over spline like function and lower region
                if ( multi ) {
                    ret = getMultiIntSplineGauss(xminval,x1,m0,m1);
                }
                else {
                    ret = getIntSplineGaussOverX(xminval,x1,mval);
                }
                return ret;   
            }
            else {
                // integrate over entire range
                if ( fabs(dxval) > 0.0 ) {
                    // integrate over spline and exponentials
                    if ( multi ) {
                        ret = getMultiIntSplineGauss(xminval,xmidval,m0,m1)
                            + getMultiIntExponentialsGauss(xmidval,x1,m0,m1);
                    } 
                    else {
                        ret = getIntSplineGaussOverX(xminval,xmidval,mval) 
                            + getIntExponentialsGaussOverX(xmidval,x1,mval);
                    }
                    return ret;
                }
                else {
                    // integrate over exponentials (no spline)
                    if ( multi ) {
                        ret = getMultiIntExponentialsGauss(xminval,x1,m0,m1);
                    }
                    else {
                        ret = getIntExponentialsGaussOverX(xminval,x1,mval);
                    }
                    return ret;
                }
            }
        } 
        else if ( dxval > 0.0 && x0 <= xmidval ) { 
            if ( x1 <= xmidval && x0 > xminval ) {
                // integrate over spline like function only
                if ( multi ) {
                    ret = getMultiIntSplineGauss(x0,x1,m0,m1);
                }
                else {
                    ret = getIntSplineGaussOverX(x0,x1,mval);
                }
                return ret;
            }
            else {
                // integrate over spline and exponential
                if ( multi ) {
                    ret = getMultiIntSplineGauss(x0,xmidval,m0,m1) 
                        + getMultiIntExponentialsGauss(xmidval,x1,m0,m1);
                }
                else {
                    ret = getIntSplineGaussOverX(x0,xmidval,mval) 
                        + getIntExponentialsGaussOverX(xmidval,x1,mval);
                }
                return ret;
            }
        }
        else {
            // integrate over exponential part
            if ( multi ) {
                ret = getMultiIntExponentialsGauss(x0,x1,m0,m1);
            }
            else {
                ret = getIntExponentialsGaussOverX(x0,x1,mval);
            }
            return ret;
        }

    }
    else if ( 2 == code && flag_2D ) {
        // integrate over m
        // Ranges
        Double_t xval = x;
        Double_t dm2val = dm2; 

        // ensure that mass range is compatible with x values
        Double_t m0 = m.min(rangeName);
        Double_t m1 = m.max(rangeName);
        
        // base gaussian
        ret = getIntGaussNotXintegrated(m0,m1,xval);

        // in case of abspdf or no pdf for x
        if ( flag_conditional_x ) {
           return ret;
        }
        if ( flag_abspdf ) {
           ret *= ((RooAbsPdf&)abspdf.arg()).getValV();
           return ret;
        }

        // pdf for x
        Double_t dxval = dx;
        Double_t xminval = xmin;

        // split into outer/turnon/exp
        if ( xval <= xminval ) {
            // outer region where pdf = 0
            ret = 0.0;
            return ret;
        } 
        else if ( xval >= xminval + dxval ) {
            // sum over exponent terms
            ret *= getValExponentials(xval);
            return ret;
        }
        else {
            // spline part
            ret *= getValSpline(xval);        
            return ret;
        }
        return ret;
    } 
    else {
        // must not go here
        assert(1==0);
    }
}
