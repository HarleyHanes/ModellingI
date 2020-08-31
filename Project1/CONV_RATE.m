%This program determines convergence rate (i.e. error as a function of step
%size) for the implicit Euler and trapezoidal methods

set(0,'DefaultAxesFontSize',18,'defaultlinelinewidth',2);set(gca,'FontSize',18);close(gcf);
%Define parameters
clear
m = 4; %mass
k = 16; %stiffness
c = 2; %damping

numberofstepsizes = 5; %how many step sizes to test
a = 0; %start at 10^a
b = -6; %end at 10^b (doesn't run in a reasonable time if b=-7)
Ks = (1/2).^(8+(0:4)); %create a vector of step sizes, log-spaced

y0 = 2; %initial displacement
y1 = 30; %initial velocity

%Define analytical solution parameters
nu = sqrt(4*k*m - c^2)/(2*m);
a1 = y0;
a2 = (y1 + c*y0/(2*m))/nu;

%Create system matrices and vectors -- NOT step-size dependent
A = [0 1; -k/m -c/m];
z0 = [y0;y1];
r0 = [y0;y1];

z(:,1) = z0; %initial conditions (position and velocity)
r(:,1) = r0; %same
y(1) = y0; %initial position

EulerError = [];
TrapError = [];

%Iterate over step sizes
for stepnumber = 1:numberofstepsizes
    
    dt = Ks(stepnumber);
    
    %Define the time grid
    T = 0:dt:20; %always testing 20 seconds
    N = length(T);

    %Define system matrices and vectors -- step-size dependent
    AA = inv(eye(2,2) - dt*A); %this is [I-kA]^(-1)
    AA1 = (eye(2,2) + 0.5*dt*A); %this is [I+0.5kA]
    AA2 = inv(eye(2,2) - 0.5*dt*A); %this is [I-0.5kA]^(-1)
    
    %Solve analytically, implicit Euler, and trapezoid
    for j=2:N %number of time points
        %Do I need to re-set y,z,r??
        
        
        t = T(j); %current time
        y(j) = exp(-c*t/(2*m))*(a1*cos(nu*t) + a2*sin(nu*t));   % True solution
        z(:,j) = AA*z(:,j-1);                                   % Implicit Euler
        r(:,j) = AA2*AA1*r(:,j-1);                              % Trapezoidal
    end
    
    %Calculate total error IN POSITION across timespan for both methods
    EulerError = [EulerError max(abs((z(1,:)-transpose(y(:)))))]; %maximum error at any point in timespan
    TrapError = [TrapError max(abs((r(1,:)-transpose(y(:)))))]; %same for trapezoidal method
    
end

%Plot error for both methods as a function of step size
figure()
loglog(Ks,EulerError,'s-g',Ks,TrapError,'o-b');%Ks,0*EulerError,
h = gca;
set(h, 'Xdir', 'reverse')
set(h,'FontSize',[18]);
xlabel('Step Size')
ylabel('Max Error')

%Add lines making clear the slope of the convergences
hold on
d=[1,.1,.01,.001, .0001,.00001, .000001];
baseError=mean([max(TrapError),max(EulerError)]);
o1Convergence=baseError*d;
o2Convergence=baseError*d.^2;
hold on
plot(d,o1Convergence,'s:g')
plot(d,o2Convergence,'o:b')
legend('Implicit Euler','Trapezoidal','O(h) Convergence','O(h^2) Convergence','Location','northeast')

