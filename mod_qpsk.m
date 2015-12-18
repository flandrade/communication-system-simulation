function constelacionArray2 = mod_qpsk(bits)

    %Variables to divide message into packet of 2 bits
    tamanio=size(bits,2); %Total size of packet
    div=1;
    simbols=[]; %symbol matrix

    while(div<tamanio)
        %Message division into packets
        m=bits(div:div+1);
        div=div+2;

        %Transform to symbols
        s=bin2dec(num2str(m));
        simbols=[simbols s];
    end


    for j=1:size(simbols,2)
       if simbols(j)==0
           constelacionArray2(j,1)=-1; %real value
           constelacionArray2(j,2)=-1; %imaginary value 
       elseif simbols(j)==1
           constelacionArray2(j,1)=-1; %real value
           constelacionArray2(j,2)=1;  %imaginary value
       elseif simbols(j)==3
           constelacionArray2(j,1)=1; %real value
           constelacionArray2(j,2)=1; %imaginary value
       else
           constelacionArray2(j,1)=1; %real value
           constelacionArray2(j,2)=-1;%imaginary value
       end    
    end   

end