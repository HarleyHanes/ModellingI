%Project2-main

%% Master Control
%Determine Problem Set
ProblemSet='1a'; %1a, 1b1, lb2, or 2

%Build structure

%% Read Data
if strcmpi(ProblemSet,'2')
    D=load('Data.txt');
    Data.x=D(:,1);
    Data.y=D(:,2);
    clear('D')
else
    Data.y=[78.5 74.3 104.3 87.6 95.9 109.2 102.7 72.5 93.1 115.9 83.8 113.3 109.4]';
    Data.x=[7 26 6 60;
            1 29 15 52;
            11 56 8 20;
            11 31 8 47;
            7 52 6 33;
            11 55 9 22;
            3 71 17 6;
            1 31 22 44;
            2 54 18 22;
            21 47 4 26;
            1 40 23 34;
            11 66 9 12;
            10 68 8 12];
    if strcmpi(ProblemSet,'1a')
    elseif strcmpi(ProblemSet,'1b1')
        Data.x=Data.x(:,1);
    elseif strcmpi(ProblemSet,'1b2')
        Data.x=Data.x(:,1:2);
    else 
        error('Problem Set not recognized')
    end
end

%% Estimate Parameters
%Determine Model
switch ProblemSet
    case '1a'
        fModel=@(x,Coeff)(Coeff.*x);
    case '1b1'
        fModel=@(x,Coeff)(Coeff.*x);
    case '1b2'
        fModel=@(x,Coeff)(Coeff.*x);
    case '2'
        fmodel=@(Coeff)(ode45(@(y)(springODE(y,Coeff)),Data.x,[0,2]));
end

%% Calculate Residuals


%% Get Confidence Intervals

%% Plot Results
%Plot Residuals

%Plot Model Result
if strcmpi(ProblemSet,'2')
    
end