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
%% Simulation Plot
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
            uSol_Numeric=SystemSolver(xVec,tVec,sin(pi*xVec/2),.7,[0 0]);
        %Plot Numeric Solutions
            plot(xVec,uSol_Numeric(:,end))
            clear('xVec','tVec','uSol_numeric')
    end
    %Plot Analytic Solutions
    %Format Plot
        title('$\frac{\partial u}{\partial t}=.7\frac{\partial^2 u}{\partial x^2}$','Interpreter','LaTex')
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
%Spatial Convergence
dT=10^(-4);
dXsteps=(10^(-2).*[1 1/2 1/4 1/8 1/16]);
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
            errorSpace(2,iStep)=max(abs(uSol_Numeric-uSol_Analytic));
        %Caculate Relative Change in error
        if iStep~=1
            errorSpace(3,iStep)=errorSpace(2,iStep-1)/errorSpace(2,iStep);
        end
end
%Temporal Convergence
    %Get Spatial and Time Steps
        dX=10^(-4);
            xVec=linspace(0,1,dX)';
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
                uSol_Numeric=SystemSolver(xVec,tVec,sin(pi*xVec/2),.7,[0 0],0);
            %Get Error
                errorTime(2,iStep)=max(abs(uSol_Numeric-uSol_Analytic));
            %Caculate Relative Change in error
            if iStep~=1
                errorTime(3,iStep)=errorTime(2,iStep-1)/errorTime(2,iStep);
            end
    end
   
end

