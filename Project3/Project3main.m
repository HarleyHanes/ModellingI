%Project3 main
clc; clear; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Problem 1
    clear; close all
    %Assign Data
        Data.rho=[.021 .037 .058 .067 .103]';
        Data.u=[14.3 10.3 7.2 4.9 2.7]';
    %Define function for u(rho)
        uFuncLinear=@(rho,Coeff)(Coeff(1)*(1-rho/Coeff(2)));
        %uFuncQuad=@(rho,Coeff)(Coeff(1)*(1-rho.^2/Coeff(2)));
        %uFuncRecip=@(rho,Coeff)(Coeff(1)*(1+1./(Coeff(2)*rho)));
        uFuncLog=@(rho,Coeff)(Coeff(1)*(1+log(rho)/(Coeff(2))));
    %Determine Best fit with fminsearch
        %Define Residual Func
        residualSumSquaresLinear=@(Coeff)sum((uFuncLinear(Data.rho,Coeff)-Data.u).^2);
        %residualSumSquaresQuad=@(Coeff)sum((uFuncQuad(Data.rho,Coeff)-Data.u).^2);
        %residualSumSquaresRecip=@(Coeff)sum((uFuncRecip(Data.rho,Coeff)-Data.u).^2);
        residualSumSquaresLog=@(Coeff)sum((uFuncLog(Data.rho,Coeff)-Data.u).^2);
        %Solve
        CoeffLinear=fminsearch(residualSumSquaresLinear,[20 1]);
        %CoeffQuad=fminsearch(residualSumSquaresQuad,[20 1]);
        %CoeffRecip=fminsearch(residualSumSquaresRecip,[20 1]);
        CoeffLog=fminsearch(residualSumSquaresLog,[20 1]);
    %Plot Data and model
        plot(Data.rho,Data.u,'-k')
        hold on
            rhoIter=linspace(.02,.11,50);
        plot(rhoIter,uFuncLinear(rhoIter,CoeffLinear),'--b');
        %plot(rhoIter,uFuncQuad(rhoIter,CoeffQuad));
        %plot(rhoIter,uFuncRecip(rhoIter,CoeffRecip));
        plot(rhoIter,uFuncLog(rhoIter,CoeffLog),'-.g');
        legend('Data',sprintf('$u(\\rho)=%.2g(1-\\frac{\\rho}{%.2g})$, RSS=%.2g',CoeffLinear,residualSumSquaresLinear(CoeffLinear)),...
            ...%sprintf('$%.2g(1-\\frac{\\rho\\textsuperscript{2}}{%.2g})$',CoeffQuad),...%sprintf('$%.2g(1-\\frac{1}{%.2g\\rho})$',CoeffRecip),...
            sprintf('$u(\\rho)=%.2g(1-\\frac{log(\\rho)}{%.2g})$, RSS=%.2g',CoeffLog,residualSumSquaresLog(CoeffLog)),...
            'Interpreter','LaTex')
        xlabel('$\rho$ (cars/minute)','Interpreter','LaTex')
        ylabel('$u$ (m/s)','Interpreter','LaTex')
%% Problem 3
%Make Plots and Get Errors
    [errorTime,errorSpace] = SolveSimpleDiffusion('Finite Differences');
%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace,SpaceRowNames,SpaceColNames,'filename','P3_SpaceErrorTable')
    %Time Table
        matrixToTexTable(errorTime,SpaceRowNames,SpaceColNames,'filename','P3_TimeErrorTable')
        
        
        
        

%% Problem 4
clear;
%Make Plots and Get Errors
    [errorTime,errorSpace] = SolveSimpleDiffusion('Finite Element');
%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace,SpaceRowNames,SpaceColNames,'filename','P4_SpaceErrorTable')
    %Time Table
        matrixToTexTable(errorTime,SpaceRowNames,SpaceColNames,'filename','P4_TimeErrorTable')
        
        

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
%Get Accurate Answer
            xVec=(0:10^(-3)/16:2)';
            tVec=(0:dT:1)';
            uSol_Accurate(:,1)=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
errorSpace=NaN(length(dXsteps),3);
errorSpace(:,1)=dXsteps;
for iStep=1:length(dXsteps)
        %Calculate Simulation Variables
            xVec=(0:dXsteps(iStep):2)';
            tVec=(0:dT:1)';
        %Get Numerical Solutions
            uSol_Numeric(:,1)=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Get Error
            %Extract Indices from Accurate answer
            stepConversion=(length(uSol_Accurate)-1)/(length(xVec)-1);
            if stepConversion~=round(stepConversion)
                error('Accurate and Test stepsizes dont match')
            end
            accurateConverted=uSol_Accurate(1:125:end);
            %Calculate Error
            errorSpace(2,iStep)=max(abs(uSol_Numeric-accurateConverted));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(3,iStep)=errorSpace(2,iStep-1)/errorSpace(2,iStep);
        end
end

%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace',SpaceRowNames,SpaceColNames,'filename','P3_SpaceErrorTable')
%% Temporal Convergence
clear
    %Get Spatial and Time Steps
        dX=10^(-3);
            xVec=linspace(0,1,dX)';
        dTsteps=10^(-2).*[1 1/2 1/4 1/8 1/16];
    %Setup Error Table
        errorTime=NaN(length(dTsteps),3);
        errorTime(:,1)=dTsteps;
%Get Accurate Answer
            tVec=(0:10^(-4)/16:2)';
            xVec=(0:dX:1)';
            uSol_Accurate(:,1)=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
    %Get Test Solutions
    for iStep=1:length(dTsteps)
            %Calculate Simulation Variables
                tVec=0:dTsteps(iStep):1;
            %Get Numerical Solutions
                uSol_Numeric=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
            %Get Error
                errorTime(iStep,2)=max(abs(uSol_Numeric-uSol_Accurate));
            %Caculate Relative Change in error
            if iStep~=1
                errorTime(iStep,3)=errorTime(iStep-1,2)/errorTime(iStep,2);
            end
    end
%Write Tables to LaTex
    %Make Column and table names
        TimeColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        TimeRowNames=cell(0,0);
    %Time Table
        matrixToTexTable(errorTime,TimeRowNames,TimeColNames,'filename','P3_TimeErrorTable')
