function [errorTime,errorSpace] = SolveSimpleDiffusion(solverType)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%% Select Solver
    if strcmpi(solverType,'Finite Differences')
        SystemSolver=@BackwardEuler1DCenteredSpace;
    elseif strcmpi(solverType,'Finite Element')
        SystemSolver=@BackwardEuler1DFiniteElement;
    else
        error('Unrecongized solver type')
    end
        uSol_Analytic=@(x,t)exp(-.7*(pi/2)^2*t).*sin(pi*x/2);
%% Simulation Plot
    %Set up Simulation Parameters
            numSpacePoints=[11 41 161];
            numTimeSteps=[11 41 161]*4;
    %Initialize Plot
            figure
            hold on
    %Plot Analytic Solutions
        plot(0:.05:2,uSol_Analytic(0:.05:2,1),LineSpec(1))
    %Loop over Step Values for plotting
    for iStep=1:length(numSpacePoints)
        %Calculate Simulation Variables
            xVec=linspace(0,2,numSpacePoints(iStep))';
            tVec=linspace(0,1,numTimeSteps(iStep))';
        %Get Numerical Solutions
            uSol_Numeric=SystemSolver(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric(:,end),LineSpec(iStep+1))
            clear('xVec','tVec','uSol_numeric')
    end
    %Format Plot
        title('$\frac{\partial u}{\partial t}=.7\frac{\partial^2 u}{\partial x^2}$','Interpreter','LaTex')
        xlabel('x','Interpreter','LaTex')
        ylabel('u(x,1)','Interpreter','Latex')
        %Make Legend,
        legendCell=cell(1,length(numSpacePoints)+1); %Stores legend strings
        legendCell{1}='True Solution';
        for  iLegend=1:length(numSpacePoints)
            h=2/(numSpacePoints(iLegend)-1);
            k=1/(numTimeSteps(iLegend)-1);
            legendCell{iLegend+1}=sprintf('$\\hat{u}$ for $h=%.2g$, $k=%.2g$',h,k);
        end
        legend(legendCell,'Interpreter','Latex')
        hold off
        
%% Convergence Tables
%Spatial Convergence
dT=10^(-7);
dXsteps=[.02 .01 .005];
%dXsteps=(10^(-2).*[1 1/2 1/4 1/8 1/16]);
errorSpace=NaN(length(dXsteps),3);
errorSpace(:,1)=dXsteps;
for iStep=1:length(dXsteps)
        %Calculate Simulation Variables
            xVec=(0:dXsteps(iStep):2)';
            tVec=(0:dT:1)';
        %Get Analytic Solution
         %Get Numerical Solutions
            uSol_Numeric=SystemSolver(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Get Error
            errorSpace(iStep,2)=max(abs(uSol_Numeric-uSol_Analytic(xVec,1)));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(iStep,3)=errorSpace(iStep,2)/errorSpace(iStep-1,2);
        end
end
%Temporal Convergence
    %Get Spatial and Time Steps
        dX=5*10^(-3);
            xVec=(0:dX:2)';
        dTsteps=[.0005 .00025 .000125];
    %Setup Error Table
        errorTime=NaN(length(dTsteps),3);
        errorTime(:,1)=dTsteps;
    %Get Analytic Solution
    %Get Numeric Solutions
    for iStep=1:length(numSpacePoints)
            %Calculate Simulation Variables
                tVec=(0:dTsteps(iStep):1)';
            %Get Numerical Solutions
                uSol_Numeric=SystemSolver(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
            %Get Error
                errorTime(iStep,2)=max(abs(uSol_Numeric-uSol_Analytic(xVec,1)));
            %Caculate Relative Change in error
            if iStep~=1
                errorTime(iStep,3)=errorTime(iStep,2)/errorTime(iStep-1,2);
            end
    end
   
end

