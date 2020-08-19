%This program plots the solutions for various numerical techniques against
%the analytical solution for the spring equation
%It also plots the convergence rate for the various methods

%Define parameters
clear
m = 4; %mass
k = 16; %stiffness
c = 2; %damping

dt = input('dt = ');  %0.5 is a good one to see  
T = 0:dt:20;
N = length(T);

y0 = 2; %initial displacement
y1 = 30; %initial velocity

%Define analytical solution parameters
nu = sqrt(4*k*m - c^2)/(2*m);
a1 = y0;
a2 = (y1 + c*y0/(2*m))/nu;

%Create system matrices and vectors
A = [0 1; -k/m -c/m];
AA = inv(eye(2,2) - dt*A); %this is [I-kA]^(-1)
AA1 = (eye(2,2) + 0.5*dt*A); %this is [I+0.5kA]
AA2 = inv(eye(2,2) - 0.5*dt*A); %this is [I-0.5kA]^(-1)
z0 = [y0;y1];
r0 = [y0;y1];
n0 = [y0;y1];

%Solve analytically, using implicit euler, trapezoid method
z(:,1) = z0; %z((pos,vel),time) I think
r(:,1) = r0; %same
n(:,1) = n0;
y(1) = y0; %initial position

for j=2:N %number of time points
   t = T(j); %current time
   y(j) = exp(-c*t/(2*m))*(a1*cos(nu*t) + a2*sin(nu*t));   % True solution
   z(:,j) = AA*z(:,j-1);                                   % Implicit Euler
   r(:,j) = AA2*AA1*r(:,j-1);                              % Trapezoidal
end

%Solve using ode23
[t,n] = ode23(@(t,n) spring_eq(t,n,m,c,k),T,n0);

%Plot the solutions 
figure()
DisplacementZ = z(1,:); %vector of positions at all times
DisplacementR = r(1,:);
DisplacementO = n(:,1);

plot(T,DisplacementZ,'-g',T,DisplacementR,'-b',T,DisplacementO,'-y',T,y,'--r',T,0*DisplacementZ,'linewidth',2);
h = gca;
set(h,'FontSize',[18]);
xlabel('Time (s)')
ylabel('Displacement (m)')
legend('Implicit Euler','Trapezoidal','Ode23','True')



