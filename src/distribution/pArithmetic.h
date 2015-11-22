#ifndef PRECISECALC
#define PRECISECALC
#include <string>
extern const bool printWarning;
//
//
/// Handle precision of log and exp operations when the value is tiny.
class pArithmetic
{
    public :
    static const long double logInf_ld ;
    static const long double Inf_ld    ;

    static const double logInf ;
    static const double Inf    ;

    ///  This gives a trancation warning message if the value is too small for expl(logInf).
    ///  The return value is double, because the smallest value in exp space is e-2500, which 
    ///  can be represented in double in log space.
    inline static double myLog(long double val, std::string who)
    {
        if(val < Inf_ld)
        {
            if(printWarning)
                Rprintf("[warning] trancating value when log(%Le) performing %s.\n", val, who.c_str());
            return logInf_ld;
        }
        return logl(val);
    }
    inline static long double myExp(double val, std::string who)
    {
        if(val < logInf_ld)
        {
            if(printWarning)
                Rprintf("[warning] exp(%e) generates a 0 performing %s.\n", val, who.c_str());
            return Inf_ld;
        }
        return expl(val);
    }
};
const long double pArithmetic::logInf_ld = logl(1.000e-4500L) ;
const long double pArithmetic::Inf_ld    = expl(logl(1.000e-4500L)) ;
const double pArithmetic::logInf = log(1e-310) ;
const double pArithmetic::Inf    = exp(log(1e-310)) ;



#endif
