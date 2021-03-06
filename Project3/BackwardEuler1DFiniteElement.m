function [uSol] = BackwardEuler1DFiniteElement(x,t,uInit,D,bounds,varargin)
%ForwardTime1DFiniteElement Numerical PDE integrator for 1D 
%diffusion with forcing and variable boundary conditions
%   Uses backward euler in time and centered difference in space
%   D- diffusion term
if nargin>6
    error('Too many input arguments for BackwardEuler1DCenteredSpace')
elseif nargin==6
    sourceFunc=varargin{1};
end


%Check Inputs

%Define Method Varaibles
    %x- Recannotate to remove boundary values
    x=x(2:end-1);
    %uSol- Solution matrix where each column is a solution at ti
    uSol=NaN(length(x),1);
        uSol(:,1)=uInit(2:end-1);
    %deltaX
    deltaX=x(2)-x(1);
        %Confirm constant spacing in x
        if sum(abs((x(2:end)-x(1:end-1))-deltaX))>10^(-12)
            error('Spacing in x is not constant')
        end
    %deltaT
    deltaT=t(2)-t(1);
        %Confirm constant spacing in t
        if max(abs((t(2:end)-t(1:end-1))-deltaT))>10^(-12)
            error('Spacing in t is not constant')
        end
    %M
        mMainDiag=2/3*ones(1,length(x));
        mOffDiag=ones(1,length(x)-1)/6;
        M=deltaX*(diag(mMainDiag,0)+diag(mOffDiag,-1)+diag(mOffDiag,1));
    %K
        kMainDiag=2*ones(1,length(x));
        kOffDiag=-ones(1,length(x)-1);
        K=D/deltaX*(diag(kMainDiag,0)+diag(kOffDiag,-1)+diag(kOffDiag,1));
    %A-Linear Finite Element operator
    	A=eye(length(x))+deltaT*(M\K);
        Ainv=inv(A);
   % keyboard
%Iterate over time
for it=1:length(t)-1
    uSol=Ainv*uSol;
end
uSol=[bounds(1);uSol;bounds(2)];
end



