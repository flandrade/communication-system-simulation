function bitsdemod = demod_qpsk(datosRX)

    %QPSK reference
    constelacionArrayRef=[-1-1i;-1+1i;1-1i;1+1i];

    %Minmum distance signal
    senalesRX=datosRX(:,1)+1i*datosRX(:,2); %signal with noise
    senalesRX=senalesRX.';

    %Matrix to determinate minimum distance
    temp1=(repmat(constelacionArrayRef.',length(senalesRX),1)).';
    temp2=repmat(senalesRX,4,1);

    %Symbol according distance
    x1=abs(temp2-temp1);
    [distancia simbolsFinal]=min(x1);

    %Symbols to bits
    simbolsFinal=simbolsFinal-1;
    bits=dec_bin(simbolsFinal, 2);

    %Matrix to vector
    bitsdemod=[];
    for i=1:size(bits,1)
        for j=1:size(bits,2)
            bitsdemod=[bitsdemod bits(i,j)];
        end
    end 
end