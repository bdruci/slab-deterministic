function [ a ] = legendre( l, x )
%LEGENDRE Summary of this function goes here
%   Detailed explanation goes here
if(l == 1)
    a = 1;
elseif(l == 2)
    a = x;
elseif(l == 3)
    a = .5*(3*x^2-1);
elseif(l == 4)
    a = .5*(5*x^3-3*x);
elseif(l == 5) 
    a = .125*(35*x^4-30*x^2+3);
else 
    a = 0.125*(63*x^5-70*x^3+15*x);
end


end

