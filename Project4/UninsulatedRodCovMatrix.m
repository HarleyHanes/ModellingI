function [CovMat,yS] = UninsulatedRodCovMatrix(x,fModel,residualsOptimal,coeffOptimal,modelType)
%UNINSULATEDRODCOVMATRIX Summary of this function goes here
%   Detailed explanation goes here

%Redfine inputs so column by column
if strcmpi(modelType,'Combined')
    fCov=@(Coeff)reshape(fModel(x,Coeff),2*length(x),1);
    residualsOptimal=reshape(residualsOptimal,2*length(x),1);
else
    fCov=@(Coeff)fModel(x,Coeff);
end

n=length(residualsOptimal);
p=length(coeffOptimal);
yS=sqrt(1/(n-p)*(residualsOptimal'*residualsOptimal));
%
JacMat=getJacobian(fCov,coeffOptimal);
CovMat=yS^2*(JacMat'*JacMat)^(-1);
end

