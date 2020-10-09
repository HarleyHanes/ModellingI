function [T] = UninsulatedRodEquil(x,coeff,dataType)
%UNINSULATEDROD Summary of this function goes here
%   Detailed explanation goes here
switch dataType
    case 'Aluminum'
        Tambient=21.289686;
    case 'Copper'
        Tambient=22.284733;
end
T=coeff(1)*exp(-coeff(3)*x)+coeff(2)*exp(ceoff(3)*x)+Tambient;
end

