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
%Number of level
level=32;

%TYPE OF QUANTIZATION
% Choose to graph:
% 1 = Univorm
% 2 = Mu-law
% 3 = A-Law
option_quantization=1;

%----------------- Modulation ---------------------------------------------
%TYPE OF MODULATION
% Choose to graph:
% 1 = BPSK
% 2 = QPSK
% 3 = BPSK and QPSK
option_modulation=3;

%%
%--------------------------------------------------------------------------
%---------------------- LOADING VOIDE -------------------------------------
%--------------------------------------------------------------------------
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

%%
%--------------------------------------------------------------------------
%---------------------- QUANTIZATION --------------------------------------
%--------------------------------------------------------------------------
%Quantization
[y1, x2, errorcuantizacion] = quantize(x,option_quantization,level);

%Quantization error
quantization_error = strcat(num2str(errorcuantizacion),' %')

%Variables to plot
xg=x2; yq=y1;
xq=x; fmq=fm;

%PLOT
%Plotting input signal with level of quantization
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
title('Input signal with level')

%Ploting input signal quantized
subplot(2,1,2)
plot(xq)
axis([ 0 4500 min(xq) max(xq)])
grid on
xlabel('samples')
title('Input signal quantized')

%%
%--------------------------------------------------------------------------
%----------------------- CODIFICATION -------------------------------------
%--------------------------------------------------------------------------
%--------------------- NO CODIFICATION ------------------------------------
%Transform decimal to bin
x3=x2-1;
bits=dec_bin(x3,log2(level));

%Matrix to vector
tem=[];
for i=1:size(bits,1)
    for j=1:size(bits,2)
        tem=[tem bits(i,j)];
    end
end
bitsc1=tem;

%--------------------- HAMMING 7,4 ----------------------------------------
% Dimensiones y matrices
n=7;
k=4;

%P matrix
P=[1 1 0; 0 1 1; 1 1 1; 1 0 1];

%Matrix generator
identity=eye(k);
G=[P identity];

%Variables to divide messages in packages
tamanio=size(bitsc1,2); %Size of message
div=1;
matrizc=[];

while(div<tamanio)
    %Divide message in packets
    m=bitsc1(div:div+3);
    div=div+4;

    %Matrix c=mG
    c=mod(m*G,2);
    matrizc=[matrizc c];
end

%--------------------- CONVOLUTIONAL CODE ---------------------------------
%Generator matrix
g=[1 0 1;1 1 1];

%Codification
bits_conv=Cnv_encd(g,1,bitsc1);

%%
%--------------------------------------------------------------------------
%----------------------- MODULATION ---------------------------------------
%--------------------------------------------------------------------------
% Array to plot constellation
switch option_modulation
    case 1
        constellationArray=[-1,0;1,0];
        modulation_name = 'BPSK';
    case 2
        constellationArray=[-1,1;-1,-1;1,1;1,-1];
        modulation_name = 'QPSK';
    case 3
        constellationArray=[-1,1;-1,-1;1,1;1,-1];
        constellationArray1=[-1,0;1,0];
        modulation_name = 'BPSK and QPSK';
end

%PLOT
%Plot axis
figure(3)
plot([-2 2],[0 0],'b');
hold on
plot([0 0 ],[-2 2],'b');

%Plot constellation
for i=1:size(constellationArray,1)
    hold on
    p=plot(constellationArray(i,1),constellationArray(i,2),'*');
    set(p,'Color','red','LineWidth',2);
    if option_modulation==3 && i<3
        hold on
    	p1=plot(constellationArray1(i,1),constellationArray1(i,2),'*');
        set(p1,'Color','black','LineWidth',2);
    end
end

title(strcat('Constellation Diagram: ', modulation_name))
axis([-2 2 -2 2])

%--------------------- MODULATION WITH NO CODIFICATION --------------------
%Modulation BPSK
bitsm1=bitsc1*2-1;

%Modulation QPSK
bitsmqpsk1=mod_qpsk(bitsc1);

%--------------------- MODULATION WITH HAMMING CODE -----------------------
%Modulation BPSK
bitsm2=2*matrizc-1;

%Modulation QPSK
bitsmqpsk2=mod_qpsk(matrizc);

%--------------------- MODULATION WITH CONVOLUTIONAL CODES ----------------
%Modulation BPSK
bitsm3=bits_conv*2-1;

%Modulation QPSK
bitsmqpsk3=mod_qpsk(bits_conv);

%%
%--------------------------------------------------------------------------
%----------------------- BER CURVES ---------------------------------------
%--------------------------------------------------------------------------
%--------------------- Variables to plot BER ------------------------------
%Probability of error (BPSK)
pet1=[]; %No codification
pet2=[]; %Hamming 7,4
pet3=[]; %Convolutional Hard
pet4=[]; %Convolutional Soft

%Probability of error (QPSK)
pet1q=[]; %No codification
pet2q=[]; %Hamming 7,4
pet3q=[]; %Convolutional Hard

%Data for Hamming codes
n=7;
k=4;

%Matriz x (Hamming)
P=[1 1 0; 0 1 1; 1 1 1; 1 0 1];

%Generator matrix (Hamming)
Identidad=eye(k);
G=[P Identidad];

%Matrix H (Hamming)
Identidad=eye(n-k);
H=[Identidad P'];

%Syndrome decoding table (Hamming)
Identidad=eye(n);
t1=[zeros(1,n); Identidad];
tabla=t1*H';

%Generator matrix for Convolutional codes
g=[1 0 1;1 1 1];

% Eb/N0 The energy per bit to noise power spectral density ratio
% 1 to 6 where 6 is the least noisy
ebn0db=0:1:6;

%--------------------- LOOP TO PLOT BER -----------------------------------
for eb=ebn0db
    %AWGN Channel
    ebn0=10^(eb/10);
    sigma=1/sqrt(2*ebn0);

    %Variables to caculate error (no codification)
    %BPSK
    error1=0;
    pe1=[];
    %QPSK
    error1q=0;
    pe1q=[];

    %Variables to caculate error (Hamming 7, 4)
    %BPSK
    error2=0;
    pe2=[];
    %QPSK
    error2q=0;
    pe2q=[];

    %Variables to caculate error (Convolutional codes, hard decision)
    %BPSK
    error3=0;
    pe3=[];
    %QPSK
    error3q=0;
    pe3q=[];

    %Variables to caculate error (Convolutional codes, soft decision)
    %BPSK
    error4=0;
    pe4=[];

        for numiter=1:5

            %--------------- NO CODIFICATION ------------------------------
            N=size(bitsm1,2);

            % CALCULATE FOR BPSK MODULATION
            if option_modulation==1 || option_modulation==3
                %AWGN Channel (BPSK)
                ruido=normrnd(0,sigma,1,N);
                y1=bitsm1+ruido;

                %BPSK demodulation
                bitsmr1=sign(y1);
                bitsr1=(bitsmr1+1)/2;

                %BER - no codification (BPSK)
                error1=error1+sum(xor(bitsc1,bitsr1));
                pe1=[pe1 sum(xor(bitsc1,bitsr1))/N];
            end

            % CALCULATE FOR QPSK MODULATION
            if option_modulation==2 || option_modulation==3
                %AWGN Channel (QPSK)
                num=size(bitsmqpsk1,1);
                ruidoReal=normrnd(0,sigma,1,num);
                ruidoImaginario=normrnd(0,sigma,1,num);

                %Noise + data
                datosRX(:,1)=bitsmqpsk1(:,1)+ruidoReal'; %real value + noise
                datosRX(:,2)=bitsmqpsk1(:,2)+ruidoImaginario'; %imaginary value + value

                %QPSK demodulation
                bitsr1q=demod_qpsk(datosRX);

                %BER - codification (QPSK)
                error1q=error1q+sum(xor(bitsc1,bitsr1q));
                pe1q=[pe1q sum(xor(bitsc1,bitsr1q))/num];
            end

            %--------------- HAMMING 7,4 ----------------------------------
            N=size(bitsm2,2);

            % CALCULATE FOR BPSK MODULATION
            if option_modulation==1 || option_modulation==3

                %AWGN Channel (BPSK)
                ruido=normrnd(0,sigma,1,N);
                y2=bitsm2+ruido;

                %BPSK demodulation
                bitsmr2=sign(y2);
                bitsr2=(bitsmr2+1)/2;

                %Variables to divide message in packets
                % message + parity bit
                tamanio=size(bitsr2,2);
                div=1; %message + parity bit

                bitsr=[]; %decoded matrix - only message

                while(div<tamanio)
                    %Divide message into packets
                    bits2=bitsr2(div:div+6);

                    %Syndrom
                    sindrome=mod(bits2*H',2);

                    %Correcction
                    for i=1:(n+1)
                        if tabla(i,:)==sindrome
                            if i~=1
                                bits2(i-1)=~bits2(i-1);
                            end
                        end
                    end

                    %Message decoded
                    b=bits2(4:7);
                    bitsr=[bitsr b];

                   %Increase div
                    div=div+7;
                end

                %BER - Hamming 7,4 (BPSK)
                error2=error2+sum(xor(bitsc1,bitsr));
                pe2=[pe2 sum(xor(bitsc1,bitsr))/N];

            end

            % CALCULATE FOR QPSK MODULATION
            if option_modulation==2 || option_modulation==3
                %Noise generation (QPSK)
                num2=size(bitsmqpsk2,1);
                ruidoReal=normrnd(0,sigma,1,num2);
                ruidoImaginario=normrnd(0,sigma,1,num2);

                %Noise + data
                datosRX2(:,1)=bitsmqpsk2(:,1)+ruidoReal'; %real value + noise
                datosRX2(:,2)=bitsmqpsk2(:,2)+ruidoImaginario'; %imaginary value + noise

                %QPSK demodulation
                bitsr2q=demod_qpsk(datosRX2);

                %Variables to divide message in packets
                % message + parity bit
                tamanio=size(bitsr2q,2);
                div=1; %message + parity bit

                bitsrq=[]; %decoded matrix - only message (QPSK)

                while(div<tamanio)

                    %Division of message into packets
                    bits2q=bitsr2q(div:div+6);
                    sindromeq=mod(bits2q*H',2);

                    for i=1:(n+1)
                        if tabla(i,:)==sindromeq
                            if i~=1
                                bits2q(i-1)=~bits2q(i-1);
                            end
                        end
                    end

                    %Message decoded
                    bq=bits2q(4:7);
                    bitsrq=[bitsrq bq];

                    %Increase div
                    div=div+7;
                end

                %BER - Hamming 7,4 (QPSK)
                error2q=error2q+sum(xor(bitsc1,bitsrq));
                pe2q=[pe2q sum(xor(bitsc1,bitsrq))/N];
            end


            %--------------- CONVOLUTIONAL CODES --------------------------
            N=size(bitsm3,2);

            % CALCULATE FOR BPSK MODULATION
            if option_modulation==1 || option_modulation==3
                %AWGN Channel (BPSK)
                ruidoReal=normrnd(0,sigma,1,N);
                ruidoImaginario=normrnd(0,sigma,1,N);
                ruido=ruidoReal+1i*ruidoImaginario;

                y4=bitsm3+ruido;

                %BPSK demodulation (hard decision)
                bitsmr3=sign(real(y4));
                bitsr3=(bitsmr3+1)/2;

                %Hard decision decodification (BPSK)
                bits3=viterbi(g,1,bitsr3);

                %BER - Hard decision (BPSK)
                error3=error3+sum(xor(bitsc1,bits3));
                pe3=[pe3 sum(xor(bitsc1,bits3))/N];

                %Soft decision decodification
                bits4=viterbi_s(g,1,y4);

                %BER - soft decision (BPSK)
                error4=error4+sum(xor(bitsc1,bits4));
                pe4=[pe4 sum(xor(bitsc1,bits4))/N];
            end

            % CALCULATE FOR QPSK MODULATION
            if option_modulation==2 || option_modulation==3
                %AWGN Channel (QPSK)
                num3=size(bitsmqpsk3,1);
                ruidoReal=normrnd(0,sigma,1,num3);
                ruidoImaginario=normrnd(0,sigma,1,num3);

                %Noise + data
                datosRX3(:,1)=bitsmqpsk3(:,1)+ruidoReal'; %real value + noise
                datosRX3(:,2)=bitsmqpsk3(:,2)+ruidoImaginario'; %imaginary value + value

                %QPSK demodulation
                bitsr3q=demod_qpsk(datosRX3);

                %Hard decision decodification
                bits3q=viterbi(g,1,bitsr3q);

                %BER of hard decision (QPSK)
                error3q=error3q+sum(xor(bitsc1,bits3q));
                pe3q=[pe3q sum(xor(bitsc1,bits3q))/N];
                pe1q=[pe1q sum(xor(bitsc1,bitsr1q))/num];

            end

        end

        %BPSK
        if option_modulation==1 || option_modulation==3
            pet1=[pet1 mean(pe1)];
            pet2=[pet2 mean(pe2)];
            pet3=[pet3 mean(pe3)];
            pet4=[pet4 mean(pe4)];
        end

        %QPSK
        if option_modulation==2 || option_modulation==3
        pet1q=[pet1q mean(pe1q)];
        pet2q=[pet2q mean(pe2q)];
        pet3q=[pet3q mean(pe3q)];
        end
end

%Pe error (Probability of error)
%It gives the average rate of occurrence of decoding errors.
%BPSK
if option_modulation==1 || option_modulation==3
    errorpe_bpsk_nocod=mean(pet1)
    errorpe_bpsk_hamming=mean(pet2)
    errorpe_bpsk_hard=mean(pet3)
    errorpe_bpsk_soft=mean(pet4)
end
%QPSK
if option_modulation==2 || option_modulation==3
    errorpe_qpsk_nocod=mean(pet1q)
    errorpe_qpsk_hamming=mean(pet2q)
    errorpe_qpsk_hard=mean(pet3q)
end

%PLOT
%BER according modulation
figure(4)
switch option_modulation
    case 1
        semilogy(ebn0db,pet1,'ko-');
        hold on;
        semilogy(ebn0db,pet2,'rx-');
        hold on;
        semilogy(ebn0db,pet3,'b+-');
        hold on;
        semilogy(ebn0db,pet4,'g*-');
        xlabel('Eb/N0, dB')
        ylabel('Bit Error Rate')
        title('BER Curves (BPSK)')
        ylim([10^(-4) 10^(-1)]);
        hleg = legend('No codification','Hamming','Hard','Soft');
        set(hleg,'Location','EastOutside')
        grid on;
    case 2
        semilogy(ebn0db,pet1q,'ko-');
        hold on;
        semilogy(ebn0db,pet2q,'rx-');
        hold on;
        semilogy(ebn0db,pet3q,'b+-');
        xlabel('Eb/N0, dB')
        ylabel('Bit Error Rate')
        title('BER Curves (QPSK)')
        ylim([10^(-4) 10^(-1)]);
        hleg = legend('No codification','Hamming','Hard');
        set(hleg,'Location','EastOutside')
        grid on;
    case 3
        semilogy(ebn0db,pet1,'r');
        hold on;
        semilogy(ebn0db,pet2,'b');
        hold on;
        semilogy(ebn0db,pet3,'k');
        hold on;
        semilogy(ebn0db,pet4,'m');
        hold on
        semilogy(ebn0db,pet1q,'--rx');
        hold on;
        semilogy(ebn0db,pet2q,'--bx');
        hold on;
        semilogy(ebn0db,pet3q,'--kx');
        xlabel('Eb/N0, dB')
        ylabel('Bit Error Rate')
        title('BER Curves (QPSK and BPSK)')
        ylim([10^(-4) 10^(-1)]);
        hleg = legend('No codification (BPSK)','Hamming (BPSK)','Hard (BPSK)', 'Soft (BPSK)','No codification (QPSK)','Hamming (QPSK)','Hard (QPSK)');
        set(hleg,'Location','EastOutside')
        grid on;
end
