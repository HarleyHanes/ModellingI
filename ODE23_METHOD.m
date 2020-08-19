%      PROGRAM ODE23_METHOD.M
%
  clear all
  set(0,'DefaultAxesFontSize',18,'defaultlinelinewidth',2);set(gca,'FontSize',18);close(gcf);

% 
% Define variables
%

  m = 4; %mass
  k = 16; %stiffness
  c = 2; %damping

  %Do we need this?
  dt = input('dt = ');    
  T = 0:dt:20;
  N = length(T);

  y0 = 2; %initial displacement
  y1 = 30; %initial velocity
  z0 = [y0;y1];

  %Define analytical solution parameters
  nu = sqrt(4*k*m - c^2)/(2*m);
  a1 = y0;
  a2 = (y1 + c*y0/(2*m))/nu;
  
  
%
% Solve analytically
%

  z(:,1) = z0; %z((pos,vel),time) I think
  y(1) = y0; %initial position
  for j=2:N %number of time points
     t = T(j); %current time
     y(j) = exp(-c*t/(2*m))*(a1*cos(nu*t) + a2*sin(nu*t));   % True solution
  end
  
%
% Solve using ode23
%

  [t,z] = ode23(@(t,z) spring_eq(t,z,m,c,k),T,z0);

%
% Plot the results
%

  figure(1)
  Displacement = z(:,1); %vector of positions at all times (swapped this from previous implemenation)
  plot(T,Displacement,'-g',T,y,'--r',T,0*z,'linewidth',2);
  h = gca;
  set(h,'FontSize',[18]);
  xlabel('Time (s)')
  ylabel('Displacement (m)')
  legend('ode23','True') %removed ,4 from end
  
  