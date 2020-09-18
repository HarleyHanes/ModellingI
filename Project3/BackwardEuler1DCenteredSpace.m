function [uSol] = BackwardEuler1DCenteredSpace(x,t,uInit,D)
%ForwardTime1DCenteredSpace Numerical PDE integrator for 1D 
%diffusion with out forcing and with 0 boundary conditions
%   Uses backward euler in time and centered difference in space
%   D- diffusion term


%Check Inputs

%Define Method Varaibles
    %uSol- Solution matrix where each column is a solution at ti
    uSol=NaN(length(x),length(t));
        uSol(:,1)=uInit;
    %deltaX
    deltaX=x(2)-x(1);
        %Confirm constant spacing in x
        if sum(abs((x(2:end)-x(1:end-1))-deltaX))>10^(-12)
            error('Spacing in x is not constant')
        end
    %deltaT
    deltaT=t(2)-t(1);
        %Confirm constant spacing in t
        if sum(abs((t(2:end)-t(1:end-1))-deltaT))>10^(-12)
            error('Spacing in x is not constant')
        end
    %A- Centered difference linear operator
        lambda=D*deltaT/deltaX^2;
        d1=1+lambda*2*ones(1,length(x));
        d2=-lambda*ones(1,length(x)-1);
    A=diag(d1,0)+diag(d2,1)+diag(d2,-1);
    
%Iterate over time
for it=1:length(t)
    uSol(:,it+1)=A\uSol(:,it);
end

end

