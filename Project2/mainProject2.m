%Project2-main
clear; close all;
set(0,'DefaultAxesFontSize',18,'defaultlinelinewidth',2);set(gca,'FontSize',18);close(gcf);
%% Master Control
%Determine Problem Set
ProblemSet='1a'; %1a, 1b1, lb2, or 2
Confidence=.95;

%Build structure

%% Read Data
if strcmpi(ProblemSet,'2')
    D=load('Data.txt');
    Data.x=D(:,1);
    Data.y=D(:,2);
    clear('D')
else
    Data.y=[78.5 74.3 104.3 87.6 95.9 109.2 102.7 72.5 93.1 115.9 83.8 113.3 109.4]';
    Data.x=[1 7 26 6 60;
            1 1 29 15 52;
            1 11 56 8 20;
            1 11 31 8 47;
            1 7 52 6 33;
            1 11 55 9 22;
            1 3 71 17 6;
            1 1 31 22 44;
            1 2 54 18 22;
            1 21 47 4 26;
            1 1 40 23 34;
            1 11 66 9 12;
            1 10 68 8 12];
    if strcmpi(ProblemSet,'1a')
    elseif strcmpi(ProblemSet,'1b1')
        Data.x=Data.x(:,1:2);
    elseif strcmpi(ProblemSet,'1b2')
        Data.x=Data.x(:,1:3);
    else 
        error('Problem Set not recognized')
    end
end
Data.n=size(Data.x,1); 

%% Estimate Parameters
%Determine Model Functions
switch ProblemSet
    case '1a'
        fModel=@(x,Coeff)sum(Coeff.*x,2);
    case '1b1'
        fModel=@(x,Coeff)sum(Coeff.*x,2);
    case '1b2'
        fModel=@(x,Coeff)sum(Coeff.*x,2);
    case '2'
        fModel=@(tspan,Coeff)springODE(Coeff,tspan);
end
%Estimate Parameters
if strcmpi(ProblemSet,'2')
   Coeff0=[1,1];
   residualSumSquaresFnc=@(Coeff)sum((fModel(Data.x,Coeff)-Data.y).^2);
   Coeff=fminsearch(residualSumSquaresFnc,Coeff0);
else
    X=Data.x;
    Coeff=((X'*X)\(X'*Data.y))';
end
fprintf('\nModel Coeffecients Found!')
switch ProblemSet
    case '1a'
        fprintf('\nbeta0=%.2f, beta1=%.2f,\nbeta2=%.2f, beta3=%.2f, beta4=%.2f\n',Coeff)
    case '1b2'
        fprintf('\nbeta0=%.2f, beta1=%.2f, beta2=%.2f',Coeff)
    case '1b1'
        fprintf('\nbeta0=%.2f, beta1=%.2f',Coeff)
    case '2'
        fprintf('\nC=%.2f, K=%.2f',Coeff)
end
yEst=fModel(Data.x,Coeff);
%Coeff=Data.y\X

%% Calculate Residuals
Residual=Data.y-yEst;


%% Get Confidence Intervals
n=Data.n; p = length(Coeff);
yBar=mean(yEst);
yS=sqrt(1/(n-p)*(Residual'*Residual));

alpha = 1 - Confidence;
pLo = alpha/2;
pUp = 1 - alpha/2;
crit = tinv([pLo pUp], n-p);

%Get t-values


if strcmpi(ProblemSet,'2') %Nonlin Formula
    Jac=getJacobian(@(Coeff)fModel(Data.x,Coeff),Coeff);
    chisq = (Jac'*Jac)^(-1);
    qCov=yS^2*(Jac'*Jac)^(-1);
    del = sqrt([chisq(1,1) chisq(2,2)])';
else %Lin Formula
    qCov=yS^2*(X'*X)^(-1);
end
    ci95 = Coeff' + crit.*diag(qCov)
    ci2sig = Coeff' + [-2 2].*diag(qCov)
  

% q95Confidence=Coeff'+(tInt(1,:).*sqrt(diag(qCov)))
% q2stdConfidence=Coeff'+(tInt(2,:).*sqrt(diag(qCov)))

%% Plot Results
%Plot Residuals
figure(1)
if strcmpi(ProblemSet,'2')
    plot(Data.x,Residual,'*')
    hold on
    plot(Data.x,2*(yS).*ones(1,length(Data.x)),'--b')
    plot(Data.x,mean(Residual)*ones(size(Data.x)),'-')
    plot(Data.x,2*(-yS).*ones(1,length(Data.x)),'--b')
    xlabel('Time')
else
    plot(1:length(Residual),Residual,'*')
    hold on
    plot(1:length(Residual),2*(-yS).*ones(1,length(Residual)),'--b')
    plot(1:length(Residual),mean(Residual)*ones(size(Data.x)),'-r')
    plot(1:length(Residual),2*yS.*ones(1,length(Residual)),'--b')
    xlabel('Observation')
    axis([1 length(Residual), -2.1*yS 2.1*yS])
end
legend('Residuals',...
        sprintf('.95 Confidence, s=%.2f',yS),sprintf('Mean Residual=%.2e',mean(Residual)))

ylabel('$\hat{y}-y$','Interpreter','Latex')


%Plot Model Result
if strcmpi(ProblemSet,'2')
    figure(2)
    plot(Data.x,Data.y,'or')
    hold on
    plot(Data.x,fModel(Data.x,Coeff),'-b')
    legend('Data',sprintf('Spring Model (C=%.2f,K=%.2f)',Coeff(1),Coeff(2)))
    xlabel('Time')
    ylabel('Displacement')
end

%% Analytic Jacobian
if strcmpi(ProblemSet,'2') %only calculate when doing prob 2

C = Coeff(1); K = Coeff(2); t = [0:0.05:5]';

u = C^2-4*K; a = sqrt(u)/2;

dydc_maple = zeros(n,1); dydk_maple = zeros(n,1);
for i = 1:n
   dydc_maple(i) = -1/((C^2-4*K)^1.5)*(5*((-t(i)*(C+2*K/5)*sqrt(C^2-4*K)+(C^2-4*K)*t(i) ...
    +2*C+4*K/5)*exp(-(C-sqrt(C^2-4*K))*t(i)/2)-exp(-(C+sqrt(C^2-4*K))*t(i)/2) ...
    * (t(i)*(C+2*K/5)*sqrt(C^2-4*K)+(C^2-4*K)*t(i) + 2*C + 4*K/5)));

    dydk_maple(i) = -1/(C^2-4*K)^1.5*((t(i)*(C+10)*sqrt(C^2-4*K)+(C^2-4*K)*t(i) ...
        -2*C-20)*exp(-(C-sqrt(C^2-4*K))*t(i)/2)-(-t(i)*(C+10)*sqrt(C^2-4*K) + ...
        (C^2-4*K)*t(i)-2*C-20)*exp(-(C+sqrt(C^2-4*K))*t(i)/2));
end

Jac_analytic = [dydc_maple dydk_maple];
approxJacError=max(abs(Jac- Jac_analytic))

chisq_a = (Jac_analytic'*Jac_analytic)^(-1); 

qCov_analytic=yS^2*chisq_a;


    del_a = sqrt([chisq_a(1,1) chisq_a(2,2)])';
    ci95_a = [Coeff' + crit(1)*yS*del_a Coeff' + crit(2)*yS*del_a]
    ci2sig_a = Coeff' + [-2 2]*yS
end