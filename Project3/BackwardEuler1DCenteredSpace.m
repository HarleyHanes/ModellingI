function [uSol,uX1] = BackwardEuler1DCenteredSpace(x,t,uInit,D,bounds,varargin)
%ForwardTime1DCenteredSpace Numerical PDE integrator for 1D 
%diffusion with out forcing and with 0 boundary conditions
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
        %Make sure x is a column vec
        if size(x,1)==1
            x=x';
            sprintf('Transitioned x from row to column vector')
        end
    %uSol- Solution matrix where each column is a solution at ti
    uSol=NaN(length(x),1);
        uSol(:,1)=uInit(2:end-1);
    %deltaX
    deltaX=x(2)-x(1);
        %Confirm constant spacing in x
        if max(abs((x(2:end)-x(1:end-1))-deltaX))>10^(-12)
            error('Spacing in x is not constant')
        end
    %deltaT
    deltaT=t(2)-t(1);
        %Confirm constant spacing in t
        if max(abs((t(2:end)-t(1:end-1))-deltaT))>10^(-14)
            error('Spacing in x is not constant')
        end
    %A- Centered difference linear operator
        lambda=D*deltaT/deltaX^2;
        d1=1+lambda*2*ones(1,length(x));
        d2=-lambda*ones(1,length(x)-1);
    A=diag(d1,0)+diag(d2,1)+diag(d2,-1);
    Ainv=A^(-1);
    
    %SPECIFIC TO P5 OF PROJECT 3-Store value at x=1
    uX1=NaN(1,length(t));
    uX1(1)=uSol(x==1);
%Iterate over time
if nargin==6
    for it=1:length(t)-1
        uSol=Ainv*(uSol+deltaT*sourceFunc(x,t(it+1)));
        uX1(it+1)=uSol(x==1); %SPECIFIC TO P5 OF PROJECT 3- Store value at x=1
    end
else
    for it=1:length(t)-1
        uSol=Ainv*uSol;
    end
end
uSol=[bounds(1);uSol;bounds(2)];
end

