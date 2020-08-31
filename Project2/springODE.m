function dydt = springODE(y,Coeff)
%springODE ODE for spring model for input into ode45
%   Coef=[C,K]
dydt=zeros(2,1);
dydt(1)=y(2);
dydt(2)=-Coeff(1)*y(2)-Coeff(2)*y(1);
end

