function T = UninsulatedRodEquil(x,params,dataType)
%UNINSULATEDROD Summary of this function goes here
%   params:phi
switch dataType
    case 'Aluminum'
        Tambient=21.289686;
        k=2.37*100; %Move from cm to m
    case 'Copper'
        Tambient=22.284733;
        k=4.01*100; %Move from cm to m
end
a=.0095;b=a;
L=.7;

phi=params(1);
h=params(2);

gamma=sqrt(2*(a+b)*h/(a*b*k));
kg=k*gamma;
c1=-phi/kg*((exp(gamma*L)*(h+kg))/(exp(-gamma*L)*(h-kg)+exp(gamma*L)*(h+kg)));
c2=phi/kg+c1;
T=c1*exp(-gamma*x)+c2*exp(gamma*x)+Tambient;
end

