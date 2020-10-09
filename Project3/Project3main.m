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
        uFuncLog=@(rho,Coeff)(Coeff(1)*(1-log(rho/(Coeff(2)))));
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
        plot(Data.rho,Data.u,'*k','MarkerSize',10)
        hold on
            rhoIter=linspace(.02,.11,50);
        plot(rhoIter,uFuncLinear(rhoIter,CoeffLinear),'--b');
        %plot(rhoIter,uFuncQuad(rhoIter,CoeffQuad));
        %plot(rhoIter,uFuncRecip(rhoIter,CoeffRecip));
        plot(rhoIter,uFuncLog(rhoIter,CoeffLog),'-.g');
        legend('Data',sprintf('$u(\\rho)=%.2f(1-\\frac{\\rho}{%.2f})$, RSS=%.2g',CoeffLinear,residualSumSquaresLinear(CoeffLinear)),...
            ...%sprintf('$%.2g(1-\\frac{\\rho\\textsuperscript{2}}{%.2g})$',CoeffQuad),...%sprintf('$%.2g(1-\\frac{1}{%.2g\\rho})$',CoeffRecip),...
            sprintf('$u(\\rho)=%.2f\\left(1-ln\\left(\\frac{\\rho}{%.2f}\\right)\\right)$, RSS=%.2g',CoeffLog,residualSumSquaresLog(CoeffLog)),...
            'Interpreter','LaTex')
        xlabel('$\rho$ (cars/m)','Interpreter','LaTex')
        ylabel('$u$ (m/s)','Interpreter','LaTex')
%% Problem 3
%Make Plots and Get Errors
    [errorTime,errorSpace] = SolveSimpleDiffusion('Finite Differences');
%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        TimeColNames={'k','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace,SpaceRowNames,SpaceColNames,'filename','P3_SpaceErrorTable')
    %Time Table
        matrixToTexTable(errorTime,SpaceRowNames,TimeColNames,'filename','P3_TimeErrorTable')
        
        
        
        

%% Problem 4
clear;
%Make Plots and Get Errors
    [errorTime,errorSpace] = SolveSimpleDiffusion('Finite Element');
%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        TimeColNames={'k','Error','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace,SpaceRowNames,SpaceColNames,'filename','P4_SpaceErrorTable')
    %Time Table
        matrixToTexTable(errorTime,SpaceRowNames,TimeColNames,'filename','P4_TimeErrorTable')
        
        
keyboard
%% Problem 5- Steady State
clear;
    %Set up Simulation Parameters
            dXSteps=[.1 .05 .01];
            dTSteps=[.1 .05 .01];
            tmax=6;
    %Initialize Plot
            figure
            hold on
    %Loop over Step Values for plotting
    for iStep=1:length(dXSteps)
        %Calculate Simulation Variables
            xVec=(0:dXSteps(iStep):2)';
            tVec=(0:dTSteps(iStep):tmax)';
        %Get Numerical Solutions
            [uSol_Numeric,uSol_X1]=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0],@(x,t)(x.*(x-2).*sin(3*pi*t)));
        %Plot Numeric Solutions
            plot(tVec,uSol_X1,LineSpec(iStep))
            clear('xVec','tVec','DummyuSol_Numeric','uSol_X1')
    end
    %Format Plot
        title('$\frac{\partial u}{\partial t}=\alpha\frac{\partial^2 u}{\partial x^2}+x(x-2)sin(3\pi t)$','Interpreter','LaTex')
        xlabel('t','Interpreter','LaTex')
        ylabel('u(1,x)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(dXSteps)); %Stores legend strings
        for  iLegend=1:length(dXSteps)
            h=dXSteps(iLegend);
            k=dTSteps(iLegend);
            legendCell{iLegend}=sprintf('$\\hat{u}$ for $h=%.2g$, $k=%.2g$',h,k);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off
%% Problem 5- Solution at t=1
    clear;
    %Set up Simulation Parameters
            dXSteps=[.1 .01 .001 .0005];
            dTSteps=[.1 .01 .001 .0005];
    %Initialize Plot
            figure
            hold on
    %Loop over Step Values for plotting
    for iStep=1:length(dXSteps)
        %Calculate Simulation Variables
            xVec=(0:dXSteps(iStep):2)';
            tVec=(0:dTSteps(iStep):1)';
        %Get Numerical Solutions
            [uSol_Numeric,uSol_X1]=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0],@(x,t)(x.*(x-2).*sin(3*pi*t)));
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric,LineSpec(iStep))
            clear('xVec','tVec','DummyuSol_Numeric','uSol_X1')
    end
    %Format Plot
        title('$\frac{\partial u}{\partial t}=\alpha\frac{\partial^2 u}{\partial x^2}+x(x-2)sin(3\pi t)$','Interpreter','LaTex')
        xlabel('x','Interpreter','LaTex')
        ylabel('u(1,x)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(dXSteps)); %Stores legend strings
        for  iLegend=1:length(dXSteps)
            h=dXSteps(iLegend);
            k=dTSteps(iLegend);
            legendCell{iLegend}=sprintf('$\\hat{u}$ for $h=%.2g$, $k=%.2g$',h,k);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off
        
        
%% Problem 5- Convergence Tables for the max error at T=1
clear;
%Get "True Answer"
            dT=.00001;
            dX=.0005;
            xVec=(0:dX:2)';
            tVec=(0:dT:1)';
            uSol_Accurate(:,1)=BackwardEuler1DCenteredSpace(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
            clear('xVec','tVec')
        %Get Spacing Vectors
            dXsteps=[20 10 5]*2*dX;
            dTsteps=[20 10 5]*2*dT;
%Spatial Convergence
%Set up Error Matrix
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
            accurateConverted=uSol_Accurate(1:stepConversion:end);
            %Calculate Error
            errorSpace(iStep,2)=max(abs(uSol_Numeric-accurateConverted));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(iStep,3)=errorSpace(iStep,2)/errorSpace(iStep-1,2);
        end
        clear('uSol_Numeric')
end

%Write Tables to LaTex
    %Make Column and table names
        SpaceColNames={'h','$E_i=\max\left(\left|x_{\delta x=h}-x_{\delta x=.0005}\right|\right)$','$\frac{|E_i|}{|E_{i-1}|}$'};
        SpaceRowNames=cell(0,0);
    %Space Table
        matrixToTexTable(errorSpace,SpaceRowNames,SpaceColNames,'filename','P5_SpaceErrorTable')
% Temporal Convergence
    %Get Spatial Steps
        xVec=(0:dX:2)';
    %Setup Error Table
        errorTime=NaN(length(dTsteps),3);
        errorTime(:,1)=dTsteps;
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
                errorTime(iStep,3)=errorTime(iStep,2)/errorTime(iStep-1,2);
            end
        clear('uSol_Numeric')
    end
%Write Tables to LaTex
    %Make Column and table names
        TimeColNames={'k','$E_i=\max\left(\left|x_{\delta t=k}-x_{\delta t=.0005}\right|\right)$','$\frac{|E_i|}{|E_{i-1}|}$'};
        TimeRowNames=cell(0,0);
    %Time Table
        matrixToTexTable(errorTime,TimeRowNames,TimeColNames,'filename','P5_TimeErrorTable')
