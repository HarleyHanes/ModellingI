%% Master control
    modelType='Insulated'; %Insulated vs. Uninsulated
    dataType='Copper'; %Copper vs. Aluminum
    %1) Select Model
        switch modelType
            case 'Insulated'
                fModel=@(x,coeff)InsulatedRodEquil(x,coeff,dataType);
                fCov=@(x,coeff)InsulatedRodCovMatrix(x,coeff);
            case 'Uninsulated'
                fModel=@(t,params)UnInsulatedRodEquil(x,coeff,dataType);
                fCov=@(x,coeff)UnInsulatedRodCovMatrix(x,coeff);
            otherwise
                error('Unrecognized model type')
        end
    %2) Select Data
        switch dataType
            case 'Aluminum'
            Data.T=[96.139218	80.122101	67.655241	57.960937	50.900923,...
                44.843707	39.750191	36.159770	33.307622	31.150610,...
                29.279187	27.884737	27.180961	26.395638	25.860308]; %Aluminum
            case 'Copper'
            Data.T=[66.035828	60.036251	54.807377	50.415894	46.743714,...
                43.663094	40.760390	38.492416	36.420054	34.771424,...
                33.184676	32.355127	31.564928	30.907695	30.561850];%Copper
            otherwise 
                error('Unrecognized data type')
        end
        Data.X=[];

%% Fit Model parameters
    [coeffOptimal,residulOptimal]=GetOptimalParams(fModel,Data);

%% Get confidence intervals
    %Estimate S
        n=length(residuialOptimal);
        p=length(coeffOptimal);
        yS=sqrt(1/(n-p)*residualOptimal'*residualOptimal);
    %Get Covariance
        covMat=fCov(Data.x,coeffOptimal);
    %Caclulate Confidence intervals
        coeff95=Coeff' + 1.97*diag(qCov);
    
        
%% Plots
    %%R
    
    


