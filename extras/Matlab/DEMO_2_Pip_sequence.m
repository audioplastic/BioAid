% DEMO 2 - Simple Pip sequence
%
% This demo shows the envelope of the waveform of a sequence of tones that
% increase in intensity before and after aid processing.

close all; clear all; clc;

[ UNIQUEpars, SHAREDpars ] = getDefaultParams_onechan();

sr = SHAREDpars.SampleRate;
freq = 1000;
dBlevs = 20:10:90;
pulseDur = 0.25;
silDur = 0.15;
x = pipSequence(sr, freq, dBlevs, pulseDur, silDur);

y = bioaid( x, UNIQUEpars, SHAREDpars);


% Plotting from here down
sig2dBSPL = @(sig)20*log10(abs(sig/20e-6)+(1/(2^32)));
dt = 1/sr;
tAxis = dt:dt:dt*numel(x);

plot(tAxis, sig2dBSPL(x), 'k')
hold on
plot(tAxis, sig2dBSPL(y), 'r')

ylim([0 100]); xlim([0 max(tAxis)])
title('Level response')
xlabel('Time in seconds')
ylabel('Level in dB SPL')



