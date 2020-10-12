%Project4 main
clc; clear; close all;
set(0,'defaultLineLineWidth',4,'defaultAxesFontSize',20);

%% Master control
    modelType='Uninsulated'; %Insulated vs. Uninsulated vs. Combined
    dataType='Copper'; %Copper vs. Aluminum vs. Combined
    %1) Select Model
        switch modelType
            case 'Insulated'
                fModel=@(x,coeff)InsulatedRodEquil(x,coeff,dataType);
                fCov=@(c1,c2,c3,c4)InsulatedRodCovMatrix(c1,c2,c3,c4);
            case 'Uninsulated'
                fModel=@(x,coeff)UninsulatedRodEquil(x,coeff,dataType);
                fCov=@(c1,c2,c3,c4,c5)UninsulatedRodCovMatrix(c1,c2,c3,c4,c5);
            case 'Combined'
                fModel=@(x,coeff)CombinedRodEquil(x,coeff);
                fCov=@(c1,c2,c3,c4,c5)UninsulatedRodCovMatrix(c1,c2,c3,c4,c5);
            otherwise
                error('Unrecognized model type')
        end
    %2) Select Data
        switch dataType
            case 'Aluminum'
            Data.T=[96.139218	80.122101	67.655241	57.960937	50.900923,...
                44.843707	39.750191	36.159770	33.307622	31.150610,...
                29.279187	27.884737	27.180961	26.395638	25.860308]'; %Aluminum
                Data.X=(.1:.04:.66)';
            case 'Copper'
            Data.T=[66.035828	60.036251	54.807377	50.415894	46.743714,...
                43.663094	40.760390	38.492416	36.420054	34.771424,...
                33.184676	32.355127	31.564928	30.907695	30.561850]';%Copper
            Data.X=(.1:.04:.66)';
            case 'Combined'
            Data.T=[96.139218	80.122101	67.655241	57.960937	50.900923,...
                44.843707	39.750191	36.159770	33.307622	31.150610,...
                29.279187	27.884737	27.180961	26.395638	25.860308;...
                66.035828	60.036251	54.807377	50.415894	46.743714,...
                43.663094	40.760390	38.492416	36.420054	34.771424,...
                33.184676	32.355127	31.564928	30.907695	30.561850]'; %Aluminum
            Data.X=(.1:.04:.66)';
            otherwise 
                error('Unrecognized data type')
        end

%% Fit Model parameters
    [coeffOptimal,RSSoptimal,residualsOptimal]=GetOptimalParams(fModel,Data,modelType);

%% Get confidence intervals
    %Get Covariance and S
         [covMat,yS]=fCov(Data.X,fModel,residualsOptimal,coeffOptimal,modelType);
     %Caclulate Confidence intervals
         coeff95=coeffOptimal' + [-1.97 1.97].*diag(covMat);
         
%% Print Results
    switch modelType
        case 'Insulated'
            fprintf('\nOptimal Parameters: phi=%.3g\n',coeffOptimal)
        case 'Uninsulated'
            fprintf('\nOptimal Parameters: phi=%.3g, h=%.3g\n',coeffOptimal)
        case 'Combined'
            fprintf(['\nOptimal Parameters: phiAluminum=%.3g, phiCopper=%.3g h=%.3g\n'...
                'Confidence Intervals: [%.3g %.3g], [%.3g %.3g], [%.3g %.3g]\n'...
                'S=%.3g\n'],...
                coeffOptimal,coeff95',yS)
    end
    
        
%% Plots
    % Solution vs. Data
    figure('Renderer', 'painters', 'Position', [200 200 750 550])
        hold on
        xVec=(0:.01:.7)';
    if strcmpi(modelType,'Combined')
        modelPlots=fModel(0:.01:.7,coeffOptimal);
        plot(xVec,modelPlots(:,1),'-b')
        plot(xVec,modelPlots(:,2),'--r')
    else
        plot(xVec,fModel(xVec,coeffOptimal))
    end
    plot(Data.X,Data.T,'o')
    switch modelType
        case 'Uninsulated'
            legend(sprintf('$T_s(x)$, $\\Phi=%.3g$, $h=%.3g$',coeffOptimal),'Temperature Data','Interpreter','LaTex')
        case 'Insulated'
            legend(sprintf('$T_s(x)$, $\\Phi=%.3g$',coeffOptimal),'Temperature Data','Interpreter','LaTex')
        case 'Combined'
            legend(sprintf('$T_s(x)$, $\\Phi_{aluminum}=%.3g$, $h=%.3g$',coeffOptimal([1,3])),...
                   sprintf('$T_s(x)$, $\\Phi_{copper}=%.3g$, $h=%.3g$',coeffOptimal([2,3])),...
                   'Aluminum Data', 'Copper Data','Interpreter','LaTex');
    end
    xlabel('Distance (m)')
    ylabel('Temperature (°C)')
    title(sprintf('T_s(x) for %s rod of %s',modelType,dataType))
    
    %%Residual
    figure('Renderer', 'painters', 'Position', [200 200 750 550])
    if strcmpi(modelType,'Combined')
        plot(Data.X,residualsOptimal(:,1),'bo')
        hold on
        plot(Data.X,residualsOptimal(:,2),'rs')
        legend('Aluminum Data','Copper Data')
    else
        plot(Data.X,residualsOptimal,'o')
    end
    
    xlabel('Distance (m)')
    ylabel('Temperature (°C)')
    title(sprintf('Residuals for %s rod of %s, RSS=%.3g',modelType,dataType,RSSoptimal))
    
    
    


