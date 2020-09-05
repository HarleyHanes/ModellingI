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
    Data.x=[7 26 6 60;
            1 29 15 52;
            11 56 8 20;
            11 31 8 47;
            7 52 6 33;
            11 55 9 22;
            3 71 17 6;
            1 31 22 44;
            2 54 18 22;
            21 47 4 26;
            1 40 23 34;
            11 66 9 12;
            10 68 8 12];
    if strcmpi(ProblemSet,'1a')
    elseif strcmpi(ProblemSet,'1b1')
        Data.x=Data.x(:,1);
    elseif strcmpi(ProblemSet,'1b2')
        Data.x=Data.x(:,1:2);
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
        fprintf('\nbeta0=%.2f, beta1=%.2f,\nbeta2=%.2f, beta3=%.2f\n',Coeff)
    case '1b2'
        fprintf('\nbeta0=%.2f, beta1=%.2f',Coeff)
    case '1b1'
        fprintf('\nbeta0=%.2f',Coeff)
    case '2'
        fprintf('\nC=%.2f, K=%.2f',Coeff)
end
yEst=fModel(Data.x,Coeff);
%Coeff=Data.y\X

%% Calculate Residuals
Residual=Data.y-yEst;


%% Get Confidence Intervals
n=Data.n;
yBar=mean(yEst);
yS=sqrt(1/(n-1)*Residual'*Residual);

%Get t-values
    % 95% CI
    tUpper=tinv(1-.05/2,n-1);
    tLower=tinv(.05/2,n-1);
    tInt=[tLower,tUpper];
    %2std CI
    tUpper=tinv(1-.02275,n-1); %NOT SURE THIS RIGHT. .97725 is half CI for 2std normal
    tLower=tinv(.02275,n-1);
    tInt(2,:)=[tLower,tUpper];

if strcmpi(ProblemSet,'2') %Nonlin Formula
    Jac=getJacobian(@(Coeff)fModel(Data.x,Coeff),Coeff);
    qCov=yS^2*(Jac'*Jac)^(-1);
else %Lin Formula
    qCov=yS^2*(X'*X)^(-1);
end
q95Confidence=Coeff'+(tInt(1,:).*sqrt(diag(qCov)))
q2stdConfidence=Coeff'+(tInt(2,:).*sqrt(diag(qCov)))

%% Plot Results
%Plot Residuals
figure(1)
if strcmpi(ProblemSet,'2')
    plot(Data.x,Residual,'*')
    hold on
    plot(Data.x,mean(Residual)*ones(size(Data.x)),'-')
    plot(Data.x,2*[yS; -yS].*ones(2,length(Data.x)),'--b')
    xlabel('Time')
else
    plot(1:length(Residual),Residual,'*')
    hold on
    plot(1:length(Residual),mean(Residual)*ones(size(Data.x)),'-')
    plot(1:length(Residual),2*[yS; -yS].*ones(2,length(Residual)),'--b')
    xlabel('Observation')
    axis([1 length(Residual), -1.1*max(abs(Residual)) 1.1*max(abs(Residual))])
end
    legend('Residuals',sprintf('Mean Residual=%.2e',mean(Residual)),...
        sprintf('.95 Confidence, s=%.2f',yS))
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