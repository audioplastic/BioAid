close all; clear all; clc;

[x,sr] = wavread('twister_44kHz'); 
% x=[zeros(5000,1); x(1:5000,:)];
x=x(1:5000,:);


[ UNIQUEpars, SHAREDpars ] = getDefaultParams();
SHAREDpars.SampleRate = sr;

tic; yMEX = bioaid( x, UNIQUEpars, SHAREDpars); toc
tic; yMAT = bioaidm( x, UNIQUEpars, SHAREDpars); toc

%% Plotting from here down
dt=1/sr; 
tAxis = dt:dt:(dt*size(x,1));
figure;  subplot(2,1,1); plot(tAxis, yMAT(:,1)); hold on; plot(tAxis, yMEX(:,1),':r')
ylabel('Amplitude'); xlabel ('time'); xlim([0 max(tAxis)])
legend('pure-Matlab', 'mex')

subplot(2,1,2); plot(abs(yMAT(:,1) - yMEX(:,1)));

%%
soundsc([x; yMAT; yMEX], sr);