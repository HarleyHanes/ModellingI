%Project3 main
clc; clear; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1
    %Assign Data
        Data.rho=[.021 .037 .058 .067 .103];
        Data.u=[14.3 10.3 7.2 4.9 2.7];
    %Plot Data
        plot(Data.rho,Data.u)
        xlabel('$\rho$ (cars/minute)','Interpreter','LaTex')
        ylabel('$u$ (m/s)','Interpreter','LaTex')
    
%% Problem 3
    %Set up Simulation Parameters
            numSpaceSteps=[11 51 501];
            numTimeSteps=20;
    %Initialize Plot
            figure
            hold on
    %Loop over Step Values
    for ixStep=1:length(numSpaceSteps)
        %Calculate Simulation Variables
            xVec=linspace(0,2,numSpaceSteps(ixStep))';
            tVec=linspace(0,1,numTimeSteps)';
        %Get Numerical Solutions
            uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.07);
        %Get Error
            %uSol_Analytic=
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric(:,end))
            clear('xVec','tVec','uSol_numeric')
    end
    %Plot Analytic Solutions
    %Format Plot
        xlabel('x','Interpreter','LaTex')
        ylabel('u(x,1)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(numSpaceSteps)); %Stores legend strings
        for  iLegend=1:length(numSpaceSteps)
            h=2/(numSpaceSteps(iLegend)-1);
            legendCell{iLegend}=sprintf('$\\hat{u}$ for $h=%.2g$',h);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off

%% Problem 4

%% Problem 5