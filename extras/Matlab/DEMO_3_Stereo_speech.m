close all; clear all; clc;

[x,sr] = wavread('twister_44kHz'); 

[ UNIQUEpars, SHAREDpars ] = getDefaultParams();
SHAREDpars.SampleRate = sr;

UNIQUEpars = [UNIQUEpars UNIQUEpars]; %Make stereo parameters
UNIQUEpars(2).Band_3_Gain_dB = 10; %Modify a parameter in the right channel

x = [x x]; %make a stereo input signal

y = bioaid( x, UNIQUEpars, SHAREDpars);


% Plotting from here down
dt=1/sr; 
tAxis = dt:dt:(dt*size(x,1));
figure;  

subplot(2,1,1); plot(tAxis, x(:,1)); hold on; plot(tAxis, y(:,1),'r')
ylabel('Amplitude'); xlabel ('time'); xlim([0 max(tAxis)])
legend('input LEFT', 'output LEFT')

subplot(2,1,2); plot(tAxis, x(:,2)); hold on; plot(tAxis, y(:,2),'r')
ylabel('Amplitude'); xlabel ('time'); xlim([0 max(tAxis)])
legend('input RIGHT', 'output RIGHT')

soundsc([x; y], sr);


