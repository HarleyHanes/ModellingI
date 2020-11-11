%Project5Main
clc; clear; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);
%% Problem 1
muBar=.2;

[timeMesh, xMesh]=meshgrid(0:.1:20,0:.1:10);
B=figure('Renderer', 'painters', 'Position', [200 200 850 550]);
hold on
x=1:.1:9;
plot3(x,9*log(9./(10-x)),zeros(size(x)),LineSpec(1,'line'),'LineWidth',8)
mesh(xMesh,timeMesh,CalculateSinkoSteifer(timeMesh,xMesh,muBar));
set(gca,'CameraPosition',[33.11908518116462,125.9366027167617,6.287048544982447])
set(gca,'View',[1.541233001077586e+02,41.930308711329431])
xlabel('size (cm)','Interpreter','LaTex')
ylabel('time','Interpreter','LaTex')
zlabel('u(t,x)','Interpreter','LaTex')
legend('G(x)','Interpreter','LaTex')
title(sprintf('Density plot of fish length for $t\\leq G(x)$, $\\mu=%.2f$',muBar),'Interpreter','LaTex')

    savefig(B,"../../Figures/Project5/Problem1Mesh")                                       %Save as .fig
    saveas(B,"../../Figures/Project5/Problem2Mesh.png")

%% Problem 2
clear
%Without Vaccination
    %Get params
    [paramsUnVac,DataUnVac]=loadDataParams('Base');
    %Get Simulation
    [t,yUnVac]=ode45(@(t,y)SIRodeFunc(y,paramsUnVac),DataUnVac.t,DataUnVac.inits);
    %Plot Model
    H=figure('Renderer', 'painters', 'Position', [200 200 750 550]);
    hold on
    %Plot S, I, Culumlative Infections
    for i=1:3
        plot(t,yUnVac(:,i),LineSpec(i,'line'))
    end
    legend('S','I','R','Location','East')
    ylabel('Individuals')
    xlabel('Time (Days)')
    title(sprintf('Vaccinated Epidemic ($\\beta=8$, $\\gamma=1.5$, $y_{Total}=%.3g$)',...
        floor(yUnVac(1,2)+(yUnVac(1,1)-yUnVac(end,1)))),'Interpreter','Latex')

    savefig(H,"../../Figures/Project5/Problem2UnVaccinated")                                       %Save as .fig
    saveas(H,"../../Figures/Project5/Problem2UnVaccinated.png")

%With Vaccination
    %Get params
    [paramsUnVac,DataVac]=loadDataParams('Vaccination');
    %Get Simulation
    [t,yVac]=ode45(@(t,y)SIRodeFunc(y,paramsUnVac),DataVac.t,DataVac.inits);
    %Plot Model
    B=figure('Renderer', 'painters', 'Position', [200 200 750 550]);
    hold on
    %Plot S, I, Culumlative Infections
     for i=1:2
         semilogy(t,yVac(:,i),LineSpec(i,'line'))
     end
    legend('S','I','R','Location','East')
    ylabel('Individuals')
    xlabel('Time (Days)')
    title(sprintf('Vaccinated Epidemic ($\\beta=8$, $\\gamma=1.5$, $y_{Total}=%.3g$)',...
        floor(yVac(1,2)+(yVac(1,1)-yVac(end,1)))),'Interpreter','Latex')

    savefig(B,"../../Figures/Project5/Problem2Vaccinated")                                       %Save as .fig
    saveas(B,"../../Figures/Project5/Problem2Vaccinated.png")
    
% R0
    paramSample=rand(150,2).*paramsUnVac;
    R0=paramSample(:,1)./paramSample(:,2);
    yTotal=NaN(100,1);
    for i=1:size(paramSample,1)
        [~,ySamp]=ode45(@(t,y)SIRodeFunc(y,paramSample(i,:)),linspace(0,30,40),[999.9 .1 0]);
        yTotal(i)=sum(ySamp(end,2:3))-ySamp(1,3)-1;
    end
    %Plot Model
    C=figure('Renderer', 'painters', 'Position', [200 200 750 550]);
    
    semilogy(R0,yTotal,LineSpec(1,'marker'),'MarkerSize',4)
    axis([0 8 0 1000])
    ylabel('Total Secondary Infections','Interpreter','Latex')
    xlabel('$R_0$','Interpreter','Latex')
    savefig(C,"../../Figures/Project5/Problem2R0")                                       %Save as .fig
    saveas(C,"../../Figures/Project5/Problem2R0.png")

%% Problem 3
clear
    %Load Data
    [params,Data]=loadDataParams('Boarding');
    %Get Optimal Model
    [paramsOptimal,~]=GetOptimalSIR(Data,params);
    %Get Model SImulation
    [t,yOptimal]=ode45(@(t,y)SIRodeFunc(y,paramsOptimal),linspace(0,24,60),Data.inits);
    %Plot Model
    H=figure('Renderer', 'painters', 'Position', [200 200 750 550]);
    hold on
    %Plot S, I, R
    for i=1:3
        plot(t,yOptimal(:,i),LineSpec(i,'line'))
    end
    %Plot Data
    plot(Data.t,Data.y,LineSpec(4,'marker'))
    legend('S','I','R','Data','Location','East')
    axis([0 18 0 800])
    ylabel('Individuals')
    xlabel('Time (Days)')
    title(sprintf('Boarding School Results ($\\beta=%.3g$, $\\gamma=%.3g$, $y_{Total}=%.3g$)',...
        paramsOptimal,floor(sum(yOptimal(end,2:3))-yOptimal(1,3))),'Interpreter','Latex')

    savefig(H,"../../Figures/Project5/Problem3InfluenzaFit")              %Save as .fig
    saveas(H,"../../Figures/Project5/Problem3InfluenzaFit.png")           %Save as .png
    
    %Plot Residuals
    [RSSoptimal,residualsOptimal]=SIRrss(paramsOptimal,Data);
    D=figure('Renderer', 'painters', 'Position', [200 200 750 550]);
    hold on
    plot(Data.t,residualsOptimal,LineSpec(1,'marker'))
    xlabel('Time (days)')
    ylabel('I_{true}-I_{expected}')
    title(sprintf('Residuals of Boarding School Data Fit ($RSS=%.3g$)',RSSoptimal),'Interpreter','Latex')
    
    
    %Calculate Confidence Intervals
%     residualsOptimal=Data.y-
%     fCov=@(params)InfluenzaCOvFunc(params,Data);
%     n=length(residualsOptimal);
%     p=length(coeffOptimal);
%     yS=sqrt(1/(n-p)*(residualsOptimal'*residualsOptimal));
%     JacMat=getJacobian(fCov,coeffOptimal);
%     CovMat=yS^2*(JacMat'*JacMat)^(-1);
%     
%     
%     function Residual=InfluenzaCovFunc(params,Data)
%         [RSS,Residual] = SIRrss(params,Data);
%     end