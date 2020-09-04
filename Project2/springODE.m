function Position = springODE(Coeff,tspan)
%springODE ODE for spring model for input into ode45
%   Coef=[C,K]
    y0=[2,10];
    dydt = @(t,y)[y(2); -Coeff(1).*y(2) - Coeff(2).*y(1)];
    [x,y]=ode45(dydt,tspan,y0);
    Position=y(:,1);
end

