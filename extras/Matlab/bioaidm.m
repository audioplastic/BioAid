function [y] = bioaidm(x, UNIQUEpars, SHAREDpars)
%BIOAIDM Pure Matlab code version of the biologically inspired gain model
%   This matlab function mimics the behavior of the MEX equivalent. This
%   script is coded to help interested persons understand how the algorithm
%   functions. For this reason, the script is written using Matlab alone 
%   without object orientation. Importantly, it is not written
%   with computational efficiency as a primary goal. The MEX version of
%   this script built upon the C++ BioAid library operates between 3 and 4
%   orders of magnitude faster for longer stimuli. Therefore, the MEX version
%   should be used for stimulus generation. Full documantation is
%   included in the 'doc' folder in the root folder of this package.

%=================================================================
% Check input
%=================================================================
assert(~(size(x,2) > 2), 'Input must be a single or double channel column matrix i.e. (n,1) or (n,2)');
nChans = size(x,2);
isStereo = (nChans == 2);

%NOTE: important distinction . . .
% nChans = number of stereo channels 
% nBands = number of frequency bands

if ~isStereo
    assert(numel(UNIQUEpars)~=2,...
    ['Ambiguous parameters: Stereo parameters supplied for a mono signal.\n'...
     'Please supply a 1x1 element struct for a mono input, or convert input to stereo'] );
end

% If the stimulus is stereo but the user only supplys mono unique pars,
% assume they intended to use the same settings in each stereo channel.
if isStereo && (numel(UNIQUEpars)==1)
    UNIQUEpars = [UNIQUEpars UNIQUEpars];
end

%=================================================================
% Nick Clark's mystical utility anonymous-functions
%=================================================================
LIN_OFFSET = 1e-50; %NaN/Inf removal offset
lin2db  =  @(linVal)   20 * log10(abs(linVal) + LIN_OFFSET); %#ok<NASGU> %not used, but inclused for completeness
db2lin  =  @(dbVal)    10.^(dbVal./20) - LIN_OFFSET;
pa2dbspl=  @(paVal)    20*log10( LIN_OFFSET  +  abs(paVal)/20e-6 );
dbspl2pa=  @(dbsplVal) 20e-6 * 10 .^ (dbsplVal ./ 20) - LIN_OFFSET;
opf = @(dt,tc)deal(dt/tc,[1 dt/tc-1]); %One-pole filter-> [b,a]=opf(dt,tc)   


%=================================================================
% Set internal parameters from parameter structure
% Convert units to something more useful if necessary
%=================================================================

% Set SHAREDpars values
sr = SHAREDpars.SampleRate;
dt = 1/sr;
nBands = SHAREDpars.NumBands;

% Memory prealloc
bp_b = zeros(5,nBands); %bandpass filter coefficients
bp_a = zeros(5,nBands); %bandpass filter coefficients
MOCb = zeros(1,nBands);
MOCa = zeros(2,nBands);
MOCfactor = zeros(1,nBands);
MOCthresh = zeros(1,nBands);
MOCbuffElements = zeros(1,nBands);

for nn = 0:nBands-1
    eval(['loEdge = SHAREDpars.Band_' num2str(nn) '_LowBandEdge;']);
    eval(['hiEdge = SHAREDpars.Band_' num2str(nn) '_HighBandEdge;']);
    [bp_b(:,1+nn) bp_a(:,1+nn)] = butter(2, [loEdge hiEdge]/(sr/2), 'bandpass');
    
    eval(['tmp_MOCtc = SHAREDpars.Band_' num2str(nn) '_MOCtc;']);
    [MOCb(:,1+nn), MOCa(:,1+nn)] = opf(dt,tmp_MOCtc);
    
    eval(['MOCfactor(1+nn) = SHAREDpars.Band_' num2str(nn) '_MOCfactor;']);
    eval(['MOCthresh(1+nn) = SHAREDpars.Band_' num2str(nn) '_MOCthreshold_dBspl;']);
    
    eval(['tmp_MOClat = SHAREDpars.Band_' num2str(nn) '_MOClatency;']);
    MOCbuffElements(1+nn) = 1 + floor(tmp_MOClat * sr);
end

% Set UNIQUEpars values
% Memory prealloc
DRNLc = zeros(1,nBands,nChans);
BSthreshIN = zeros(1,nBands,nChans);
chanGain = zeros(1,nBands,nChans);

for kk = 1:nChans
    inGain(kk)  = db2lin(UNIQUEpars(kk).InputGain_dB); %#ok<AGROW>
    outGain(kk) = db2lin(UNIQUEpars(kk).OutputGain_dB);%#ok<AGROW>
    
    ARthresh(kk) = dbspl2pa(UNIQUEpars(kk).ARthreshold_dBSPL);%#ok<AGROW>
    [ARb(:,kk), ARa(:,kk)] = opf(dt,UNIQUEpars(kk).ARtc);%#ok<AGROW>
    ARbuffElements(kk) = 1 + floor(UNIQUEpars(kk).ARlatency * sr);%#ok<AGROW>
    
    for nn = 0:nBands-1
        eval(['BSthreshIN(1+nn,kk) = UNIQUEpars(kk).Band_' num2str(nn) '_InstantaneousCmpThreshold_dBspl;']);
        eval(['DRNLc(1+nn,kk) = UNIQUEpars(kk).Band_' num2str(nn) '_DRNLc;']);
        eval(['chanGain(1+nn,kk) = UNIQUEpars(kk).Band_' num2str(nn) '_Gain_dB;']);
    end
end

BSthreshIN = dbspl2pa(BSthreshIN);
chanGain = db2lin(chanGain);

DRNLb = BSthreshIN.^(1-DRNLc);
BSthreshOUT = 10.^  (  (1./(1-DRNLc)) .*  log10(DRNLb)  );

%=================================================================
% Initialize IIR filter memories
%=================================================================
ARz  = zeros( max( size(ARb, 1) , size(ARa, 1) ) -1, nChans);
MOCz = zeros(max(size(MOCb,1), size(MOCa,1))-1,nBands); %MONO
bp_zIN = zeros(max(size(bp_b,1), size(bp_a,1))-1,nBands, nChans);
bp_zOUT = zeros(size(bp_zIN));

%=================================================================
% Initialize ring-buffers for MOC and AR feedback compressors
%=================================================================
ARbuff = cell(1,nChans);
for kk = 1:nChans
    ARbuff{kk} = ones(ARbuffElements(kk), 1);
end

MOCbuff = cell(1,nBands); %MONO
for nn = 1:nBands
    MOCbuff{nn} = ones(MOCbuffElements(nn),1);
end

%=================================================================
% Internal storage parameters for post-analysis
%=================================================================
store.MOC = zeros(size(x,1), nBands);


%=================================================================
% Signal processing
%=================================================================
y = zeros(size(x));
isSample2of2 = false;

for nn=1:size(x,1)
    
    for kk = 1:nChans
        
        s = x(nn,kk)*inGain(kk); %'s' is just a single sample temporary value
        
        %--------------------------------
        % Apply acoustic reflex compression
        %--------------------------------
        s = s/ARbuff{kk}(end);
        %------------ end ---------------
        
        
        %--------------------------------
        % Commence within-band processing
        %--------------------------------
        acc = 0;
        for ff=1:nBands
            %--------------------------------
            % FIRST BP FILTERING
            %--------------------------------
            [sb, bp_zIN(:,ff,kk)] = filter(bp_b(:,ff),bp_a(:,ff),s,bp_zIN(:,ff,kk));
            %------------ end ---------------
            
            
            
            %--------------------------------
            % Apply MOC compression
            %--------------------------------
            sb = sb * MOCbuff{ff}(end);
            %------------ end ---------------
            
            
            
            %--------------------------------
            % Apply instantaneous compression
            %--------------------------------
            abs_x = abs(sb);
            if (  abs_x  >  BSthreshOUT(ff,kk)  )
                sb =  sign(sb) * DRNLb(ff,kk) * abs_x^DRNLc(ff,kk);
            end
            %else do nothing
            %------------ end ---------------
            
                                                            
            
            %--------------------------------
            % SECOND BP FILTERING
            %--------------------------------
            [sb, bp_zOUT(:,ff,kk)] = filter(bp_b(:,ff),bp_a(:,ff),sb,bp_zOUT(:,ff,kk));
            %------------ end ---------------
                        
            
            %--------------------------------
            % Update MOC ring-buffer
            %--------------------------------
            % Before calculating attenuation for MOC loop, add a tiny DC
            % offset to the incoming sample. This is so that Pa values of
            % zero do not cause a value of -1000 dB SPL to be fed to the
            % attenuation calculator. Attenuation is thresholded after
            % filtering, so any preceding zeros in a lab-based simulation
            % would cause the internal MOC level to drop incredibly low. At
            % the onset of a stimulus, it takes the internal MOC vale a
            % long time to recover. The effects of the MOC are only visible
            % once the internal value has exceeded threshold, so zeros
            % preceding the stimulus would increase the MOC latency. The
            % small DC offset rectifies this problem. However, one should
            % note that the MOC latency is intrinsically linked to the MOC
            % time constant because of this architecture.
            sb=sb+2e-05;
            %------------------
            if isStereo
                if isSample2of2
                    stereoAccumulator = stereoAccumulator + pa2dbspl(sb);
                    dBsb = stereoAccumulator/2;
                    [MOCbuff{ff}, MOCz(:,ff)] = updateMOCbuffer(dBsb, MOCbuff{ff}, MOCthresh(ff), MOCfactor(ff), MOCb(:,ff), MOCa(:,ff), MOCz(:,ff), db2lin);
                else 
                    stereoAccumulator = pa2dbspl(sb);
                end
                isSample2of2 = ~isSample2of2; %Flip bit for subsequent loop
            else % for a MONO input signal
                dBsb = pa2dbspl(sb);
                [MOCbuff{ff}, MOCz(:,ff)] = updateMOCbuffer(dBsb, MOCbuff{ff}, MOCthresh(ff), MOCfactor(ff), MOCb(:,ff), MOCa(:,ff), MOCz(:,ff), db2lin);
            end      
            store.MOC(nn,ff) = MOCbuff{ff}(end);
            
            % The small DC offset can be removed after calculation of the
            % MOC strength
            sb=sb-2e-05;
            %------------ end ---------------
            
            
            
            %--------------------------------
            % Apply gain for channel
            %--------------------------------
            acc = acc + sb*chanGain(ff,kk);
            %------------ end ---------------
        end
        s = acc;
        %-- end within-band processing --
        
        
        %--------------------------------
        % Update acoustic reflex ring-buffer
        %--------------------------------
        tmp = s * s; % power
        [tmp,ARz(:,kk)] = filter(ARb(:,kk),ARa(:,kk),tmp,ARz(:,kk)); %smooth power
        tmp = sqrt(tmp) / ARthresh(kk); % RMS relative to threshold
        tmp = max(tmp, 1); % Stop AR giving gain
        ARbuff{kk} = [tmp; ARbuff{kk}(1:end-1)];
        %------------ end ---------------
        
        y(nn,kk) = s*outGain(kk); %pop 's' into the output vector
        
    end
end

end %END OF BIOAIDM FUNCTION



function [buff, MOCz] = updateMOCbuffer(dBsample, buff, MOCthresh, MOCfactor, MOCb, MOCa, MOCz, db2lin_func)   
    dBsample = max(dBsample-MOCthresh,0) * MOCfactor; % dB    
    [dBsample,MOCz] = filter(MOCb,MOCa,dBsample,MOCz);    
    dBsample = db2lin_func(-dBsample); %incase you're confused, this is a function handle passed to this function (D.R.Y.)     
    buff = [dBsample; buff(1:end-1)];
end

