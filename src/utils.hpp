//
//  utils.h Created by nrclark on 03/05/2012.
//  BioAid
//
// Misc helper functions.
//
// Copyright (C) 2012 Nick Clark
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of this
// software and associated documentation files (the "Software"), to deal in the Software
// without restriction, including without limitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following
// conditions:
// 
// The above copyright notice and this permission notice shall be included in all copies
// or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef __UTILS_H_B7FA53D3__
#define __UTILS_H_B7FA53D3__

#include <math.h>
#include <iostream>
#include <ostream>


#define LIN_OFFSET (1e-20f) // = -400 dB to protect from -Inf if log of zero attempted


// The DBGM macro version here allows you to use nice syntax like this ...
// DBGM("first: " << first << ", second: " << second);
#if defined(DEBUG) || defined(MATLAB_VERBOSE)
#define DBGM(x) std::cout << x << std::endl;
#else
#define DBGM(x)
#endif



class utils
{
public:    
    
    //---- STATIC METHODS -- (These should all be templatized)
    inline static float lin2db(const float linVal)
    {
        return 20.0f*log10f(fabs(linVal) + LIN_OFFSET) ;
    };
    
    inline static float db2lin(const float dbVal)
    {
        return powf(10.0, dbVal/20.0f) - LIN_OFFSET;
    };
    
    inline static float pa2dbspl(const float paVal)
    {
        return 20.0f*log10f( LIN_OFFSET  +  fabs(paVal)/20e-6f );
    };
    
    inline static float dbspl2pa(const float dbsplVal)
    {
        return 20e-6f * powf(10.0f, dbsplVal / 20.0f) - LIN_OFFSET;
    };       
    
    //inline static void DBG(const std::string x)
    //Problem with this version is that you need to lexical_cast to show values (see macro version above)
    //{
    //#ifdef DEBUG
    //    std::cout << x << std::endl; //VS2008 does not like this
    //#endif
    //}
    
    // This code maps a value within one range to a vlue in another scale
    template<typename T>
    static T mapVal(T ipVal, 
                    T ipMin, T ipMax, 
                    T opMin, T opMax)
    {                
        //ip params
        const T ipRange = ipMax-ipMin;
        const T ipFract = (ipVal - ipMin) / ipRange;
        
        //op params
        const T opRange = opMax-opMin;    
        return opMin + ipFract*opRange;        
    };
};

#include <boost/thread/mutex.hpp>

// Nice little optional ScopedLock (the pointer can be NULL)
class NullCheckingScopedLock
{
public:
    NullCheckingScopedLock( void* _pMutex)
//    : pMutex( _pMutex)
    {
        //A reinterpret_cast from void is used to that this is the only place in the code that directly
        //refers to boost::mutex. Allows non-boost, but boost compatible mutex useage.
        pMutex = reinterpret_cast < boost::mutex* > ( _pMutex );
        if( pMutex) 
        {
            pMutex->lock();
            DBGM("Mutex Lock");
        }
    }
    
    ~NullCheckingScopedLock()
    {
        if( pMutex)
        {
            pMutex->unlock();
            DBGM("Mutex unlock")
        }
    }
private:
    boost::mutex* pMutex;
};




#ifdef MATLAB_VERBOSE
#include "mex.h"
//Following allows us to get debug messages in Matlab console
//http://stackoverflow.com/questions/243696/correctly-over-loading-a-stringbuf-to-replace-cout-in-a-matlab-mex-file
class mystream : public std::streambuf
{
public:   
protected:
    virtual int overflow(int c = EOF) 
    {
        if (c != EOF) {
            mexPrintf("%.1s",&c);
        }
        return 1;
    }

    virtual std::streamsize xsputn(const char *s, std::streamsize n) 
    {
        mexPrintf("%.*s",n,s);
        return n;
    }    
};
#endif //MATLAB_VERBOSE


#endif  // __UTILS_H_B7FA53D3__
