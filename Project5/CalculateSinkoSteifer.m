function u = CalculateSinkoSteifer(t,x,muBar)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
G=9*log(9./(10-x));
u=(t<=G).*sin(pi*((x-10).*exp(t/9)+9)/9).*exp((1/9-muBar).*t);
%u=(t<=G).*sin(pi*(x-1)/9).*exp((1/9-muBar).*t);
end

