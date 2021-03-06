%      PROGRAM TRAP_METHOD.M
%
  clear all
% 
%% Define variables - everything is the same as implicit Euler until AA
    m = 4; %mass
    k = 16; %stiffness
    c = 2; %damping

    y0 = 2; %initial displacement
    y1 = 30; %initial velocity

    %Define analytical solution parameters
    nu = sqrt(4*k*m - c^2)/(2*m);
    a1 = y0;
    a2 = (y1 + c*y0/(2*m))/nu;
    
    tStart=0;
    tEnd=20;

%% Get True Solutions
    tTrue=linspace(tStart,tEnd,200);
    yTrue= exp(-c*tTrue./(2*m)).*(a1.*cos(nu.*tTrue) + a2.*sin(nu.*tTrue));   % True solution

%% Get Numerical Solutions
  tStep = [1; .5; .25];    %smallest step size for reasonable runtime:  10^-6
  for itStep=1:length(tStep)
      dt=tStep(itStep);
      T = tStart:dt:tEnd;
      N = length(T);
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
      yEst(:,1) = z0; %z((pos,vel),time) I think
      for j=2:N %number of time points
         t = T(j); %current time
         yEst(:,j) = AA2*AA1*yEst(:,j-1);                              % Approximate solution - NEW
      end
      %Save Results to structure
      resultStruct.y=yEst;
      resultStruct.T=T;
      resultStruct.dt=dt;
      %Store structure in cell array
      Results{itStep}=resultStruct;
      %Clear simulation results
      clear(['yEst','T','dt'])
  end
  
%
% Plot the results
%   
    set(0,'DefaultAxesFontSize',18,'defaultlinelinewidth',2);set(gca,'FontSize',18);close(gcf);
    LineStyle=["-","--",":","-."];
    LineColor=['k','b','r','g'];

  figure()
  %Plot True Solution
  p=plot(tTrue,yTrue);
    set(p,'LineStyle',LineStyle(1))
    set(p,'Color',LineColor(1))
  legendList={'True'};
  hold on
  %Plot numerical solutions
  for isol=1:length(Results)
    %Load Results
    sol=Results{isol};
    T=sol.T;
    Displacement = sol.y(1,:); %vector of positions at all times
    p=plot(T,Displacement);
    set(p,'LineStyle',LineStyle(isol+1))
    set(p,'Color',LineColor(isol+1))
    legendList{end+1}=sprintf('dt=%.2f',sol.dt);
  end
    h = gca;
    set(h,'FontSize',[18]);
    xlabel('Time (s)')
  ylabel('Displacement (m)')
  title('Spring Solutions with Trapezoid Method')
  legend(legendList) %removed ,4 from end
  
  