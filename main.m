% Simulation of a communication system
%https://github.com/flandrade/communication-system-simulator
%
% Copyright 2014. Fernanda Andrade
% Universidad de las Fuerzas Armadas - ESPE
%
% Last modified 16-Dec-2015

clear all
clc

%--------------------------------------------------------------------------
%----------------- COMMUNICATION OPTIONS ----------------------------------
%--------------------------------------------------------------------------

%----------------- Quantization -------------------------------------------
%Number of levels
level=32;

%TYPE OF QUANTIZATION
% Choose: 
% 1 = Univorm
% 2 = Mu-law
% 3 = A-Law
option_quantization=1; 

%----------------- Modulation -------------------------------------------
%TYPE OF MODULATION
% Choose: 
% 1 = QPSK
% 2 = BPSK
% 3 = BPSK Y QPSK
option_modulation=1; 

%----------------- AWGN Chanel -------------------------------------------
% Eb/N0 The energy per bit to noise power spectral density ratio
% Option: 1 to 10 where 10 is the least noisy
ebno=10; 

%--------------------------------------------------------------------------
%---------------------- TRANSMITTER ---------------------------------------
%--------------------------------------------------------------------------

%--------------------- LOADING VOICE --------------------------------------
%Loading voice
[x,fm]=audioread('voz.wav');

%Fundamental frequency
N=floor(0.02*fm);
C=xcorr(x,N,'coeff');
N1=floor(0.002*fm);
[x0,vmax]=max(C(N+N1:2*N+1));
t0=(vmax+N1)/fm;
f0=1/t0;
fundamental_frequency=strcat(num2str(f0),' Hz')

%PLOT
%Plotting input signals (voices)
figure(1)
plot(x)
axis([ 0 4500 min(x) max(x) ])
title('Input signal 1');

% Playing voices
disp('Playing input signals');
soundsc(x,fm);
pause(3);


%------------------------- QUANTIZATION ----------------------------------- 
%Quantization
[y1, x2, errorcuantizacion] = quantize(x,option_quantization,level);

%Quantization error
quantization_error = strcat(num2str(errorcuantizacion),' %')

%Variables to plot 
xg=x2; yq=y1;
xq=x; fmq=fm;

%PLOT
%Plotting input signal with levels of quantization
figure(2)
subplot(2,1,1)
plot(x);
axis([ 0 4500 min(x) max(x) ])
hold on
for vv=1:level
    hold on
    plot(yq(vv,:));
end
grid on
xlabel('samples')
ylabel('x(t)')
title('Input signal with levels')

%Ploting input signal quantized 
subplot(2,1,2)
plot(xq)
axis([ 0 4500 min(xq) max(xq)])
grid on
xlabel('samples')
title('Input signal quantized')


%------------------------- MODULATION ------------------------------------- 

