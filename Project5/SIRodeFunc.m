function dydt = SIRodeFunc(y,params)
%SIRodeFunc Summary of this function goes here
%   Detailed explanation goes here
beta=params(1);
gamma=params(2);
dydt=NaN(3,1);
dydt(1)=-beta*y(2)/(sum(y))*y(1);
dydt(2)=beta*y(2)/(sum(y))*y(1)-gamma*y(2);
dydt(3)=gamma*y(2);

end

