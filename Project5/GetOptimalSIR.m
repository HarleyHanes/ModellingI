function [paramsOptimal,RSS] = GetOptimalSIR(Data,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Define RSS
RSSfunc=@(params)SIRrss(params,Data);

%Get Optimal Params
paramsOptimal=fminsearch(RSSfunc,params);

%Get Optimal RSS
RSS=RSSfunc(paramsOptimal);
end

