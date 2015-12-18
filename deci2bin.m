function y=deci2bin(x,l)
y = zeros(1,l);	
i = 1;
while x>=0 & i<=l
	y(i)=rem(x,2);
	x=(x-y(i))/2;
	i=i+1;
end
y=y(l:-1:1);