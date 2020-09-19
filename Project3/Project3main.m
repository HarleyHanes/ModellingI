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
    
%% Problem 3- Simulation Plot
    %Set up Simulation Parameters
            numSpacePoints=[11 21 31];
            numTimePoints=[6 11 21];
    %Initialize Plot
            figure
            hold on
    %Loop over Step Values for plotting
    for iStep=1:length(numSpacePoints)
        %Calculate Simulation Variables
            xVec=linspace(0,2,numSpacePoints(iStep))';
            tVec=linspace(0,1,numTimePoints(iStep))';
        %Get Numerical Solutions
            uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.07,[0 0]);
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric(:,end))
            clear('xVec','tVec','uSol_numeric')
    end
    %Plot Analytic Solutions
    %Format Plot
        xlabel('x','Interpreter','LaTex')
        ylabel('u(x,1)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(numSpacePoints)); %Stores legend strings
        for  iLegend=1:length(numSpacePoints)
            h=2/(numSpacePoints(iLegend)-1);
            k=1/(numTimePoints(iLegend)-1);
            legendCell{iLegend}=sprintf('$\\hat{u}$ for $h=%.2g$, $k=%.2g$',h,k);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off
        
%% Problem 3- Convergence Tables
clear;
%Spatial Convergence
dT=10^(-8);
numSpacePoints=2*10^(-2).*[1 1/2 1/4 1/8 1/16];
errorSpace=NaN(length(numSpacePoints),3);
errorSpace(:,1)=numSpacePoints;
for iStep=1:length(numSpacePoints)
        %Calculate Simulation Variables
            xVec=linspace(0,2,numSpacePoints(iStep))';
            tVec=linspace(0,1,dT)';
        %Get Analytic Solution
        %Get Numerical Solutions
            uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.07,[0 0]);
        %Get Error
            errorSpace(2,iStep)=max(abs(uSol_Numeric-uSol_Analytic));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(3,iStep)=errorSpace(2,iStep-1)/errorSpace(2,iStep);
        end
end
%Temporal Convergence
    %Get Spatial and Time Steps
        dX=10^(-8);
            xVec=linspace(0,1,dT)';
        numTimePoints=10^(-2).*[1 1/2 1/4 1/8 1/16];
    %Setup Error Table
        errorTime=NaN(length(numTimePoints),3);
        errorTime(:,1)=numTimePoints;
    %Get Analytic Solution
    %Get Numeric Solutions
    for iStep=1:length(numSpacePoints)
            %Calculate Simulation Variables
                tVec=linspace(0,1,numTimePoints(iStep))';
            %Get Numerical Solutions
                uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.07,[0 0]);
            %Get Error
                errorTime(2,iStep)=max(abs(uSol_Numeric-uSol_Analytic));
            %Caculate Relative Change in error
            if iStep~=1
                errorTime(3,iStep)=errorTime(2,iStep-1)/errorTime(2,iStep);
            end
    end
%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace,SpaceRowNames,SpaceColNames,'filename','P3_SpaceErrorTable')
    %Time Table
        matrixToTexTable(errorTime,SpaceRowNames,SpaceColNames,'filename','P3_TimeErrorTable')
%% Problem 4

%% Problem 5