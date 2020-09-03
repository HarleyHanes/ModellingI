%Project2-main

%% Master Control
%Determine Problem Set
ProblemSet='1a'; %1a, 1b1, lb2, or 2
Confidence=.95;

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
Data.n=size(Data.x,1);

%% Estimate Parameters
%Determine Model Functions
switch ProblemSet
    case '1a'
        fModel=@(x,Coeff)sum(Coeff.*x,2);
    case '1b1'
        fModel=@(x,Coeff)sum(Coeff.*x,2);
    case '1b2'
        fModel=@(x,Coeff)sum(Coeff.*x,2);
    case '2'
        fModel=@(x,Coeff)(ode45(@(y)(springODE(y,Coeff)),x,[0,2]));
end
%Estimate Parameters
X=Data.x;
Coeff=((X'*X)\(X'*Data.y))';
yEst=fModel(Data.x,Coeff);
%Coeff=Data.y\X

%% Calculate Residuals
Residual=Data.y-yEst;


%% Get Confidence Intervals
n=Data.n;
yBar=mean(yEst);
yS=1/(n-1)*sum(Residual.^2);

%Get t-values
    % 95% CI
    tUpper=tinv(1-.05/2,n-1);
    tLower=tinv(.05/2,n-1);
    tInt=[tLower,tUpper];
    %2std CI- 
    tUpper=tinv(1-.02275,n-1); %NOT SURE THIS RIGHT. .97725 is half CI for 2std normal
    tLower=tinv(.02275,n-1);
    tInt(2,:)=[tLower,tUpper];

%yBar Confidence Intervals
yBarInterval=yBar+tInt*yS/sqrt(n);

%yS Confidence Inrevals



%% Plot Results
%Plot Residuals
plot(1:length(Residual),Residual,'*')
ylabel('$|\hat{y}-y|$','Interpreter','Latex')
xlabel('Observation')


%Plot Model Result
if strcmpi(ProblemSet,'2')
    
end