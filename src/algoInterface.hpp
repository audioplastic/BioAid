//
//  algoInterface.h Created by nrclark on 03/05/2012.
//  BioAid
//
// Provides an interface to the underlying data structure.
// Calls should be made through this interface.
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

#ifndef BioAid_algoInterface_h
#define BioAid_algoInterface_h

#include <boost/thread/mutex.hpp>
#include <boost/scoped_ptr.hpp>

#include "params.hpp"
#include "algoComponents.hpp"



class cAidAlgo
{
private:        
    boost::scoped_ptr<cMOCsimContainer > pMOCsimContainer;    //Abstract this away, the user need not be confused by this implementation detail
    boost::scoped_ptr<cAidStereoChannelManager > pManagerL;
    boost::scoped_ptr<cAidStereoChannelManager > pManagerR;
    
    //Mutex now defined as void ptr and cast within the function to cut down on preprocessor required if MT is disabled
    void* pMyMutex;
    
    
public:
    //We have dual constructors here depending on whether you want a mono or a stereo processor
    //Params cannot be const because they have subscribers added
    cAidAlgo(cUniqueStereoParams& _lPars, 
             cUniqueStereoParams& _rPars,
             cSharedStereoParams& _sharedPars, 
             void* _pMyMutex = NULL          
             ) :
    pMOCsimContainer(new cMOCsimContainer(_sharedPars)),
    pManagerL(new cAidStereoChannelManager(_lPars, _sharedPars, *pMOCsimContainer)),
    pManagerR(new cAidStereoChannelManager(_rPars, _sharedPars, *pMOCsimContainer)),
    pMyMutex(_pMyMutex)
    {
     
    };
    
    cAidAlgo(cUniqueStereoParams& _lPars,
             cSharedStereoParams& _sharedPars, 
             void* _pMyMutex = NULL
             ) :
    pMOCsimContainer(new cMOCsimContainer(_sharedPars)),
    pManagerL(new cAidStereoChannelManager(_lPars, _sharedPars, *pMOCsimContainer)),
    pManagerR(NULL),
    pMyMutex(_pMyMutex)
    {
        
    };
    
    //cMOCsimContainer MOCsimContainer(sharedPars); 
    
    void processSampleBlock (const float ** inputChannelData,
                             int 	numInputChannels,
                             float ** 	outputChannelData,
                             int 	numOutputChannels,
                             int 	numSamples )
    {        
                
        bool isStereo = true;        
        // No point processing stereo if there is only one output channel
        // ..or if there is only a MONO processor
        if ( (numInputChannels > numOutputChannels) || (pManagerR==NULL) ) {
            isStereo = false; 
        }
        
        //This code sets the index for the right input channel
        //If both outputs are sharing a common input then the following loop needs to know
        int inRightIdx = 0;
        if(numInputChannels==2)
            inRightIdx = 1;
                
        {//New scope just for lock
            const NullCheckingScopedLock l(pMyMutex); //Lock this as short as possible                
            for (int nn=0; nn< numSamples; ++nn)
            {
                outputChannelData[0][nn] = pManagerL->process(inputChannelData[0][nn]);
                if (isStereo)
                    outputChannelData[1][nn] = pManagerR->process(inputChannelData[inRightIdx][nn]);
            }
        }//Mutex unlocked 
        
        //Copy the data into both channels if we are mono
        if(!isStereo && (numOutputChannels==2) )
        {
            for (int nn=0; nn< numSamples; ++nn)
            {  
                outputChannelData[1][nn] = outputChannelData[0][nn];
            }
        }
        
    }
    
    
};


#endif
