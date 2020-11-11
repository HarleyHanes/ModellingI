function [RSS,Residual] = SIRrss(params,Data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[t,y]=ode45(@(t,y)SIRodeFunc(y,params),Data.t,Data.inits);
Residual=Data.y-y(:,2);
RSS=sum(Residual.^2);
end

