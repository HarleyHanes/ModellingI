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
            numTimeSteps=[6 11 21];
    %Initialize Plot
            figure
            hold on
    %Loop over Step Values for plotting
    for iStep=1:length(numSpacePoints)
        %Calculate Simulation Variables
            xVec=linspace(0,2,numSpacePoints(iStep))';
            tVec=linspace(0,1,numTimeSteps(iStep))';
        %Get Numerical Solutions
            uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric(:,end))
            clear('xVec','tVec','uSol_numeric')
    end
    %Plot Analytic Solutions
    %Format Plot
        title('$\frac{\partial u}{\partial t}=\alpha\frac{\partial^2 u}{\partial x^2}$','Interpreter','LaTex')
        xlabel('x','Interpreter','LaTex')
        ylabel('u(x,1)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(numSpacePoints)); %Stores legend strings
        for  iLegend=1:length(numSpacePoints)
            h=2/(numSpacePoints(iLegend)-1);
            k=1/(numTimeSteps(iLegend)-1);
            legendCell{iLegend}=sprintf('$\\hat{u}$ for $h=%.2g$, $k=%.2g$',h,k);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off
        
%% Problem 3- Convergence Tables
clear;
%Spatial Convergence
dT=10^(-5);
dXsteps=(10^(-2).*[1 1/2 1/4 1/8 1/16]);
errorSpace=NaN(length(dXsteps),3);
errorSpace(:,1)=dXsteps;
for iStep=1:length(dXsteps)
        %Calculate Simulation Variables
            xVec=(0:dXsteps(iStep):2)';
            tVec=(0:dT:1)';
        %Get Analytic Solution
        %Get Numerical Solutions
            uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Get Error
            errorSpace(2,iStep)=max(abs(uSol_Numeric-uSol_Analytic));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(3,iStep)=errorSpace(2,iStep-1)/errorSpace(2,iStep);
        end
end
%Temporal Convergence
    %Get Spatial and Time Steps
        dX=10^(-5);
            xVec=linspace(0,1,dT)';
        dTsteps=10^(-2).*[1 1/2 1/4 1/8 1/16];
    %Setup Error Table
        errorTime=NaN(length(dTsteps),3);
        errorTime(:,1)=dTsteps;
    %Get Analytic Solution
    %Get Numeric Solutions
    for iStep=1:length(numSpacePoints)
            %Calculate Simulation Variables
                tVec=linspace(0,1,dTsteps(iStep))';
            %Get Numerical Solutions
                uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
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

%% Problem 5- Simulation Plot
    %Set up Simulation Parameters
            numSpacePoints=[11 21 41];
            numTimeSteps=[6 11 21 ];
    %Initialize Plot
            figure
            hold on
    %Loop over Step Values for plotting
    for iStep=1:length(numSpacePoints)
        %Calculate Simulation Variables
            xVec=linspace(0,2,numSpacePoints(iStep))';
            tVec=linspace(0,1,numTimeSteps(iStep))';
        %Get Numerical Solutions
            uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0],@(x,t)(x.*(x-2).*sin(3*pi*t)));
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric(:,end))
            clear('xVec','tVec','uSol_numeric')
    end
    %Format Plot
        title('$\frac{\partial u}{\partial t}=\alpha\frac{\partial^2 u}{\partial x^2}+x(x-2)sin(3\pi t)$','Interpreter','LaTex')
        xlabel('x','Interpreter','LaTex')
        ylabel('u(x,1)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(numSpacePoints)); %Stores legend strings
        for  iLegend=1:length(numSpacePoints)
            h=2/(numSpacePoints(iLegend)-1);
            k=1/(numTimeSteps(iLegend)-1);
            legendCell{iLegend}=sprintf('$\\hat{u}$ for $h=%.2g$, $k=%.2g$',h,k);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off
        
%% Problem 5- Convergence Tables
clear;
%Spatial Convergence
dT=10^(-5);
dXsteps=(10^(-2).*[1 1/2 1/4 1/8 1/16]);
errorSpace=NaN(length(dXsteps),3);
errorSpace(:,1)=dXsteps;
for iStep=1:length(dXsteps)
        %Calculate Simulation Variables
            xVec=(0:dXsteps(iStep):2)';
            tVec=(0:dT:1)';
        %Get Analytic Solution
        %Get Numerical Solutions
            uSol_Numeric(:,1)=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Get Error
            errorSpace(2,iStep)=max(abs(uSol_Numeric-uSol_Analytic));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(3,iStep)=errorSpace(2,iStep-1)/errorSpace(2,iStep);
        end
end
%Temporal Convergence
    %Get Spatial and Time Steps
        dX=10^(-5);
            xVec=linspace(0,1,dT)';
        dTsteps=10^(-2).*[1 1/2 1/4 1/8 1/16];
    %Setup Error Table
        errorTime=NaN(length(dTsteps),3);
        errorTime(:,1)=dTsteps;
    %Get Analytic Solution
    %Get Numeric Solutions
    for iStep=1:length(numSpacePoints)
            %Calculate Simulation Variables
                tVec=linspace(0,1,dTsteps(iStep))';
            %Get Numerical Solutions
                uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
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
