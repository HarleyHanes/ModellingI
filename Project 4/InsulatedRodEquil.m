function T = InsulatedRodEquil(x,params,dataType)
%INSULATEDROD Summary of this function goes here
%   param:[phi, h]
switch dataType
    case 'Aluminum'
        Tambient=21.289686;
        k=2.37*100; %Move from cm to m
    case 'Copper'
        Tambient=22.284733;
        k=4.01*100; %Move from cm to m
end
a=.0095;b=a;L=.7;


phi=params(1);

T=phi/k*(x-L)+Tambient;

% beta=2*(a+b)*h/(a*b*k);
% FractionDenominator=(k*sqrt(beta)*(exp(-sqrt(beta)*L)-exp(sqrt(beta)*L)));
% EndTempDifference=((phi*FractionDenominator/(k*sqrt(beta))-h*Tambient+phi*exp(sqrt(beta)*L))*exp(sqrt(beta)*L)+...
%     (-h*Tambient+phi*exp(sqrt(beta)*L))*exp(-sqrt(beta)*L)+...
%     Tambient*FractionDenominator)/...
%     FractionDenominator-h*(exp(sqrt(beta)*L)+exp(-sqrt(beta)*L));
% 
% FractionNumerator=(h*EndTempDifference+phi*exp(sqrt(beta)*L));
% T=(phi/(k*sqrt(beta))+FractionNumerator/FractionDenominator)*exp(sqrt(beta)*x)+...
%     FractionNumerator/FractionDenominator*exp(-sqrt(beta)*x)+...
%     Tambient;
end

