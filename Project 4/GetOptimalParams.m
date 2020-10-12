function [CoeffOptimal,RSSoptimal,ResidualsOptimal] = GetOptimalParams(fModelQOI,Data,ModelType)
%GETOPTIMALPARAMS Solves ordinary least squares problem
%   Detailed explanation goes here
switch ModelType
    case 'Insulated'
        Coeff0=1;
    case 'Uninsulated'
        Coeff0=[-1.84e+05,19];
    case 'Combined'
        Coeff0=[-99000,-180000,16];
end

%Make Cost function
    fRSS=@(Coeff)sum(sum((fModelQOI(Data.X,Coeff)-Data.T).^2));
    
%Find Params with fminsearch
    CoeffOptimal=fminsearch(fRSS,Coeff0);
    RSSoptimal=fRSS(CoeffOptimal);
    ResidualsOptimal=Data.T-fModelQOI(Data.X,CoeffOptimal);
end

