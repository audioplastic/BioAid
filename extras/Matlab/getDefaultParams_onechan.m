function [ UNIQUEpars, SHAREDpars ] = getDefaultParams_onechan(  )
%GETDEFAULTPARAMS Returns a pair of default parameter structures
%   This function returns some default parameter sets for the gain model.
%   The parameter sets can be modified once they have been created, or
%   alternatively they can be generated from scratch if desired. The
%   main use of this helper function is to reduce code bloat in the demo
%   files.


%%
SHAREDpars.SampleRate = 44100;
SHAREDpars.NumBands = 1;

for nn = 0:SHAREDpars.NumBands-1 %Channels are indexed from zero!
    cf = 1000;
    bw = 1/2;
    loEdge = cf * (2^(-bw/2)); %#ok<NASGU> Warnings suppressed as value is used in eval statement below
    hiEdge = cf * (2^ (bw/2)); %#ok<NASGU> Warnings suppressed as value is used in eval statement below
    
    eval(['SHAREDpars.Band_' num2str(nn) '_LowBandEdge  = loEdge;']);
    eval(['SHAREDpars.Band_' num2str(nn) '_HighBandEdge = hiEdge;']);
    
    eval(['SHAREDpars.Band_' num2str(nn) '_MOCtc = 0.07;']);
    eval(['SHAREDpars.Band_' num2str(nn) '_MOCfactor = 0.8;']);
    eval(['SHAREDpars.Band_' num2str(nn) '_MOClatency = 0.010;']);
    eval(['SHAREDpars.Band_' num2str(nn) '_MOCthreshold_dBspl = 25;']);
end

%%
UNIQUEpars.InputGain_dB = 0;
UNIQUEpars.OutputGain_dB = 0;
UNIQUEpars.ARthreshold_dBSPL = 110;
UNIQUEpars.ARtc = 0.006;
UNIQUEpars.ARlatency = 0.005;

for nn = 0:SHAREDpars.NumBands-1 %Channels are indexed from zero!    
    eval(['UNIQUEpars.Band_' num2str(nn) '_InstantaneousCmpThreshold_dBspl  = 75;']);
    eval(['UNIQUEpars.Band_' num2str(nn) '_DRNLc = 0.2;']);       
    eval(['UNIQUEpars.Band_' num2str(nn) '_Gain_dB  =  0.0;']);
end

% UNIQUEpars = [UNIQUEpars UNIQUEpars];

end

