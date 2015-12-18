function output=cnv_encd(g,k0,input)
%  		cnv_encd(g,k0,input)
%  		determines the output sequence of a binary convolutional encoder
%  		g is the generator matrix of the convolutional code
%        	with n0 rows and l*k0 columns. Its rows are g1,g2,...,gn.
%  		k0 is the number of bits entering the encoder at each clock cycle.
%  		input The binary input seq.

    %  check to see if extra zero padding is necessary
    if rem(length(input),k0) > 0
      input=[input,zeros(size(1:k0-rem(length(input),k0)))];
    end
    n=length(input)/k0;
    %  check the size of matrix g
    if rem(size(g,2),k0) > 0
      error('Error, g is not of the right size.')
    end
    %  determine l and n0
    l=size(g,2)/k0;
    n0=size(g,1);
    %  add extra zeros
    u=[zeros(size(1:(l-1)*k0)),input,zeros(size(1:(l-1)*k0))];
    %  generate uu, a matrix whose columns are the contents of 
    %  conv. encoder at various clock cycles.
    u1=u(l*k0:-1:1);
    for i=1:n+l-2
      u1=[u1,u((i+l)*k0:-1:i*k0+1)];
    end
    uu=reshape(u1,l*k0,n+l-1);
    %  determine the output
    output=reshape(rem(g*uu,2),1,n0*(l+n-1));
end
  
