function pipOut = pipSequence(sampleRate, freq, dBlevs, pulseDur, silDur)
%PIPSEQUENCE Makes tone pips for testing
%   This is a simple utility function that generates a sequence of tone
%   pips interleaved with silence. The user gan specify the durations of
%   pips and silence, as well as the frequency and level of each of the
%   tone pips.



if nargin < 5
    silDur = 0.3;
end
if nargin < 4
    pulseDur = 0.1;
end
if nargin < 3
    dBlevs = 20:20:100;
end
if nargin < 2
    freq = 500;
end
if nargin < 1
    sampleRate = 48e3;
end

dt = 1/sampleRate;
tAxis = (dt:dt:pulseDur)';
sPulse = sin(2*pi*freq*tAxis);
sPulse = sPulse./sqrt(mean(sPulse.^2));
rms2dBspl = @(dBspl)20e-6*10^(dBspl/20); %sneaky short-hand function by (ab)using function handles
zPad = zeros(ceil(sampleRate*silDur),1);

pipOut = [];
for nn = 1:numel(dBlevs)
    pipOut = [ pipOut; sPulse*rms2dBspl(dBlevs(nn)); zPad;]; %#ok<AGROW>
end

end% ------ OF pipSequence


