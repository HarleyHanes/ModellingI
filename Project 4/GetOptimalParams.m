function [CoeffOptimal,RSSoptimal] = GetOptimalParams(fModelQOI,Data,Coeff0)
%GETOPTIMALPARAMS Solves ordinary least squares problem
%   Detailed explanation goes here

%Make Cost function
    fRSS=@(Coeff)sum((fModelQOI(Coeff)-Data).^2);
%Find Params with fminsearch
    CoeffOptimal=fminsearch(fRSS,Coeff0);
    RSSoptimal=fRSS(CoeffOptimal);
end

