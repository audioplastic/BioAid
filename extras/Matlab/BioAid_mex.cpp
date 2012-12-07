//
//  BioAid_mex.h Created by nrclark on 03/05/2012.
//  BioAid
//
// This code is the gateway function to make a Matlab funtion from the BioAid C++ code.
// The MEX version generated can be used in an identical way to the interpreted equivalent.
// To build ensure your compiler is set up correctly and type
// >> mex BioAid_m 
// ... followed by a number of compiler flags ...
// Matlab appears to be oblivious of your system path, so you'll need to tell it where 
// the boost headers and libs are. For example, if boost is installed on a mac using 
// homebrew (recommended), the following compiler flags need to be issued ...
//
// -I/usr/local/include/ 
// -L/usr/local/lib
//
// Users of other OS will need to determine their boost installation path.
// The mex file needs to be linked against boost's multithreaded system library
//
// -lboost_system-mt
//
// Windows users may also need to link against the thread library and the date-time library
// Windows prebuilt boost libs can be obtained at http://www.boostpro.com/download/
//
// Finally, the MEX file needs to be renamed. This can be done manually afterward,
// or at compile time using the 'output' flag.
//
// TEST WITH:
// >> clear mex; [x,sr] = wavread('twister_44kHz'); tic; y = bioaid(x); toc; soundsc([x;y],sr)
//
// -output bioaid
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


#include "mex.h"

#include "../../src/utils.h"
#include "../../src/params.h"
#include "../../src/algoInterface.h"

#define AUDIO_DATA_IN prhs[0] 
#define AUDIO_DATA_OUT plhs[0]

#define UNIQUE_PARS_IN prhs[1]
#define SHARED_PARS_IN prhs[2]
         
void checkParamValidity(const char* tmpName, mxArray* tmpData);

/*************************************************************************/
/* Gateway function                                                      */
/*************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
  
  #ifdef MATLAB_VERBOSE
    //Allows my DBGM macro to redirect output to matlab console
    mystream mout; //defined in utils.h
    std::streambuf *outbuf = std::cout.rdbuf(&mout); //buffer redirection
  #endif//MATLAB_VERBOSE  

  //=======================================================================
  // *** User error handling and initialization ***      
  // Check the number of input and output parameters
  // This is more important than in regular m-files, as an error at this
  // stage may cause the Matlab instance to crash 
  //=======================================================================           
  if(nrhs>3)
      mexErrMsgTxt("Too many inputs: Should be an in signal vector followed by (optional) parameter structures");
  if(nrhs<1)
      mexErrMsgTxt("Too few inputs: At minimum, a signal vector must be supplied.");
  if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments:\n"
                   "You are attempting to run the mex version of the algorithm.\n"
                   "Inspection of parameters is disabled.");
  
  bool isUniqueParsSupplied = (nrhs>1); //alias for readability
  bool isSharedParsSupplied = (nrhs>2); //alias for readability
  
  if (isUniqueParsSupplied)
      if(!mxIsStruct(UNIQUE_PARS_IN))
          mexErrMsgTxt("(optional) second input argument must be a parameter structure");
  
  if (isSharedParsSupplied)
      if(!mxIsStruct(SHARED_PARS_IN))
          mexErrMsgTxt("(optional) third input argument must be a parameter structure");
          
  int xM = mxGetM(AUDIO_DATA_IN); //columns (NxM)=Matlab
  int xN = mxGetN(AUDIO_DATA_IN); //rows (NxM)=Matlab
  
  if ((xN>2) || (xN<1))
      mexErrMsgTxt("Input must be a single or double channel column matrix i.e. (n,1) or (n,2)");
  bool isStereo = (xN==2); //alias for readability
  
  if(mxIsChar(AUDIO_DATA_IN) ||
     !mxIsNumeric(AUDIO_DATA_IN) ||
     mxIsSparse(AUDIO_DATA_IN) ||
     mxIsComplex(AUDIO_DATA_IN)   )
        mexErrMsgTxt("Audio data must be non-complex, non-sparse, and numeric");
  
  int numUniqueParSets =0; 
  int numUniqueParFields =0;
  if (isUniqueParsSupplied)
  {
    numUniqueParSets   = mxGetNumberOfElements(UNIQUE_PARS_IN);
    numUniqueParFields = mxGetNumberOfFields(UNIQUE_PARS_IN);
  }
  
  int numSharedParFields =0;
  if (isSharedParsSupplied)
  {
    if(  mxGetNumberOfElements(SHARED_PARS_IN)   >   1  )
        mexErrMsgTxt("Parameter structure for pars shared between stereo channels should only have one element");
    numSharedParFields = mxGetNumberOfFields(SHARED_PARS_IN);
  }
  
  bool isStereoPars = (numUniqueParSets == 2); //alias for readability
  
  if (numUniqueParSets > 2)
      mexErrMsgTxt("Too many parameter sets supplied (parameter struct has more than 2 elements).\n"
                   "Algorithm only functions in mono or stereo - please supply 1x1 or 1x2 element struct.");
  
  if (isStereoPars && !isStereo)
      mexErrMsgTxt("Ambiguous parameters: Stereo parameters supplied for a mono signal.\n"
                   "Please supply a 1x1 element struct for a mono input, or convert input to stereo" );
  
  
  int totalNumElements = xM*xN; /* this is elements for all channels combined */  
  
  
  //=======================================================================
  // *** Instantiate and modify parameter class(es) ***      
  // Set parameter object values according to the supplied structures
  //=======================================================================
  
  //Create 2: We might waste one, but they are cheap and it saves dynamic allocation
  cUniqueStereoParams uniquePars[2]; 
  cSharedStereoParams sharedPars; // We only need one of these
  
  mxArray* tmpData;
  const char* tmpName;  
  
  //============= UNIQUE ================
  for (int structIdx=0; structIdx<numUniqueParSets; ++structIdx)
  {
      for (int fieldIdx=0; fieldIdx<numUniqueParFields; ++fieldIdx)
      {
          tmpData = mxGetFieldByNumber(UNIQUE_PARS_IN, structIdx, fieldIdx);
          tmpName = mxGetFieldNameByNumber(UNIQUE_PARS_IN, fieldIdx);
                    
          checkParamValidity(tmpName, tmpData);
                                 
          if(  uniquePars[structIdx].setParam(tmpName, (float)mxGetScalar(tmpData))  ) //The actual setting occurs here (returns TRUE on error)
          {
              std::string errMsg = std::string("Field named \"") + std::string(tmpName) 
              + std::string("\" is not an existing field name in \"UNIQUE\" parameters");
              mexErrMsgTxt(errMsg.c_str());
          }
      }
  }
  
  //============= SHARED ================
  for (int fieldIdx=0; fieldIdx<numSharedParFields; ++fieldIdx)
  {
      tmpData = mxGetFieldByNumber(SHARED_PARS_IN, 0, fieldIdx);
      tmpName = mxGetFieldNameByNumber(SHARED_PARS_IN, fieldIdx);
      
      checkParamValidity(tmpName, tmpData);
      
      if(  sharedPars.setParam(tmpName, (float)mxGetScalar(tmpData))  ) //The actual setting occurs here (returns TRUE on error)
      {
          std::string errMsg = std::string("Field named \"") + std::string(tmpName) 
          + std::string("\" is not an existing field name in \"SHARED\" parameters");
          mexErrMsgTxt(errMsg.c_str());
      }
  }    
  if (isStereo) //User blocked from setting this in the "checkValidity" function
  {
      sharedPars.setParam("IsStereo", 1.f); //So we automatically set it here.
  }
  
  //=======================================================================
  // *** Signal processing ***
  //=======================================================================
  AUDIO_DATA_OUT = mxCreateDoubleMatrix(xM,xN,mxREAL);
  
  double* out  = mxGetPr(AUDIO_DATA_OUT); 
  double* in   = mxGetPr(AUDIO_DATA_IN);    
  
  float*  in32  =  (float*)mxCalloc(totalNumElements, sizeof(float));
  float*  out32 =  (float*)mxCalloc(totalNumElements, sizeof(float));
    
  
  // ---------> DOWNCAST <---------
  for (int mm = 0; mm < totalNumElements; ++mm)
      in32[mm] = (float)in[mm];
  // ---------- DOWNCAST ----------
       
  // Make the data look like 2D C-style arrays
  //(this fits with other audio APIs [like VST] nicely)
  float* in2D[] =  {in32,  &in32[xM]};
  float* out2D[] = {out32, &out32[xM]};
        
  if(!isStereo){
      cAidAlgo myAlgo(uniquePars[0], sharedPars);
      myAlgo.processSampleBlock ((const float**) in2D,  1,  (float**) out2D,   1,   xM);
      
  } else {
      if (isStereoPars){
          cAidAlgo myAlgo(uniquePars[0], uniquePars[1], sharedPars);
          myAlgo.processSampleBlock ((const float**) in2D,  2,  (float**) out2D,   2,   xM);
      } else {
          cAidAlgo myAlgo(uniquePars[0], uniquePars[0], sharedPars);
          myAlgo.processSampleBlock ((const float**) in2D,  2,  (float**) out2D,   2,   xM);
      }
  }
  
  // ---------> UPCAST <---------
  for (int mm = 0; mm < totalNumElements; ++mm)
      out[mm] = (double)out32[mm];
  // ---------- UPCAST ----------
 
  
  mxFree(in32); mxFree(out32);
  
  #ifdef MATLAB_VERBOSE
  // Play nice and restore the std stream buffer (not strictly necessary)
    std::cout.rdbuf(outbuf);  
  #endif//MATLAB_VERBOSE
}


//=======================================================================
// *** Input parameter validity checking ***
//=======================================================================

void checkParamValidity(const char* tmpName, mxArray* tmpData)
{
    std::string errMsg;
    
    //Note, the tmpData == NULL glitch was a PAIN to find. See:
    //http://stackoverflow.com/questions/10706739/strange-goings-on-when-reading-structure-data-in-mex-file
    if(  tmpData == NULL  )
    {
        errMsg = std::string("Field named \"") + std::string(tmpName) + std::string("\" has a NULL value\nHave you greated a field in one element of a struct and not set it in the other stereo channel?");
        mexErrMsgTxt(errMsg.c_str());
    }
    
    //CHECK IF THE FIELD IS NOT EMPTY    
    if(  tmpData == NULL || mxIsEmpty(tmpData)  )
    {
        errMsg = std::string("Field named \"") + std::string(tmpName) + std::string("\" is empty");
        mexErrMsgTxt(errMsg.c_str());
    }
    
    //DONT ALLOW DIRECT ACCESS TO IsStereo PARAM
    if(  !strcmp(tmpName, "IsStereo")  ) //strcmp returns zero if they are equal
    {
        mexErrMsgTxt("Programatic access to \"IsStereo\" is restricted in the mex function.\n"
                "It is set automatically depending on the dimensions of the input array.\n"
                "Please remove it from your parameter structure."); //Just a warning incase we increase the channel count.
    }        
        
    
    //CHECK THE FIELD ONLY CONTAINS A NUMERIC SCALAR
    if(mxIsChar(tmpData) ||
            !mxIsNumeric(tmpData) ||
            mxIsSparse(tmpData) ||
            mxIsComplex(tmpData) ||
            mxGetNumberOfElements(tmpData)!=1   )
    {
        errMsg = std::string("Field named \"") + std::string(tmpName) + std::string("\" must be of type: numeric, non-sparse, non complex, scalar");
        mexErrMsgTxt(errMsg.c_str());
    }       
}


//----END
