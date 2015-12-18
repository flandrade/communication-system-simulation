function y=bin2deci(x)
    l=length(x);
    y=(l-1:-1:0);
    y=2.^y;
    y=x*y';
end