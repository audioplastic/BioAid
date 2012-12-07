% DEMO 1 - processes a wav and play its output back to the user
%
% If you've complied the mex version of bioaid, use that instead by
% replacing 'bioaidm' with 'bioaid'. On my computer, bioaidm takes 31s to
% process the sound as compared with 0.05 seconds for bioaid.
%
% Tweak the input level of the stimulus.

close all; clear all; clc;

ipLVL = 94; %input level in dB SPL

[x,sr] = wavread('twister_44kHz');
x = x/sqrt(mean(x.^2)); %norm RMS to 1 ~ 94 dB SPL
x = x*10^((ipLVL-94)/20);


[ UNIQUEpars, SHAREDpars ] = getRealParams();
SHAREDpars.SampleRate = sr;

tic
y = bioaid( x, UNIQUEpars, SHAREDpars);
toc

% Plotting from here down
dt=1/sr; 
tAxis = dt:dt:(dt*size(x,1));
figure;  plot(tAxis, x(:,1)); hold on; plot(tAxis, y(:,1),':r')
ylabel('Amplitude'); xlabel ('time'); xlim([0 max(tAxis)])
legend('input', 'output')


soundsc([x; y], sr);