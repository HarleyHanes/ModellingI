%This is the function file for the spring equation
%In first order form

function dzdt = spring_eq(t,z,m,c,k)
%This function returns a vector representing the derivative at that point
%Parameters m,c,k passed into function


dzdt = [z(2); (-k/m)*z(1) + (-c/m)*z(2)]; %assuming no applied force

