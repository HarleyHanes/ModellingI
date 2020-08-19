%      PROGRAM TRAP_METHOD.M
%
  clear all

% 
% Define variables - everything is the same as implicit Euler until AA
%

  m = 4; %mass
  k = 16; %stiffness
  c = 2; %damping

  dt = input('dt = ');    %smallest step size for reasonable runtime:  10^-6
  T = 0:dt:20;
  N = length(T);

  y0 = 2; %initial displacement
  y1 = 30; %initial velocity

  %Define analytical solution parameters
  nu = sqrt(4*k*m - c^2)/(2*m);
  a1 = y0;
  a2 = (y1 + c*y0/(2*m))/nu;
  
%
% Create system matrices and vectors
%

  A = [0 1; -k/m -c/m];
  %New matrices:
  AA1 = (eye(2,2) + 0.5*dt*A); %this is [I+0.5kA]
  AA2 = inv(eye(2,2) - 0.5*dt*A); %this is [I-0.5kA]^(-1)
  z0 = [y0;y1];
  
%
% Iterate to approximate the solution using the trapezoid rule
%

  z(:,1) = z0; %z((pos,vel),time) I think
  y(1) = y0; %initial position
  for j=2:N %number of time points
     t = T(j); %current time
     y(j) = exp(-c*t/(2*m))*(a1*cos(nu*t) + a2*sin(nu*t));   % True solution
     z(:,j) = AA2*AA1*z(:,j-1);                              % Approximate solution - NEW
  end
  
%
% Plot the results
%

  figure(1)
  Displacement = z(1,:); %vector of positions at all times
  plot(T,Displacement,'-g',T,y,'--r',T,0*Displacement,'linewidth',2);
  h = gca;
  set(h,'FontSize',[18]);
  xlabel('Time (s)')
  ylabel('Displacement (m)')
  legend('Trapezoidal','True') %removed ,4 from end
  
  