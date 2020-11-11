function [params,Data] = loadDataParams(paramset)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%data

%params
    params=[8, 1.5];
switch paramset
    case 'Base'
        Data.inits=[960, 40, 0];
        Data.t=linspace(0,6,50);
    case 'Vaccination'
        Data.inits=[96,4,900];
        Data.t=linspace(0,6,50);
    case 'Boarding'
        Data.inits=[760 3 0];
        
        DataRaw=[0.0000000e+00   3.0000000e+00;
        1.0000000e+00   6.0000000e+00;
        2.0000000e+00   2.5000000e+01;
        3.0000000e+00   7.3000000e+01;
        4.0000000e+00   2.2200000e+02;
        5.0000000e+00   2.9400000e+02;
        6.0000000e+00   2.5800000e+02;
        7.0000000e+00   2.3700000e+02;
        8.0000000e+00   1.9100000e+02;
        9.0000000e+00   1.2500000e+02;
        1.0000000e+01   6.9000000e+01;
        1.1000000e+01   2.7000000e+01;
        1.2000000e+01   1.1000000e+01;
        1.3000000e+01   4.0000000e+00];
        Data.t=DataRaw(:,1);
        Data.y=DataRaw(:,2);
        
end

