function matrixToTexTable(matrix, RowNames, ColNames, varargin)
%
%input:
% matrix = the matrix to be writting as a LaTeX table
% 'ColNames', string of column names
% 'RowNames', string of row names

% outputDir/filename  is written as LaTeX matrix
%
% variable input is of the form: 'OptionName',value
% 'outputDir', name of output folder 
% 'filename', name of output file
% example:
%  

%% Add latex table to path
addpath('/Users/hhanes/Documents/Modelling/ModellingICode/latextable')
%% Input Management
p = inputParser;

addRequired(p,'matrix',@isnumeric);
addOptional(p,'outputDir','.',@ischar); % output folder
addOptional(p,'filename','matrix',@ischar); % output filename

% filename={'%s/sensitivity_table_transpose.tex'};

parse(p, matrix, varargin{:});

matrix = p.Results.matrix;
outputDir = p.Results.outputDir;
filename = p.Results.filename;

%%
[nrow,ncols] = size(matrix);

nColNames = length(ColNames);
nRowNames = length(RowNames);
% Init basic param strings, if needed
if nColNames ~= ncols
    if nColNames > 0
        keyboard
        warning('fewer column names than matrix columns')
    end
    ColNames = basicNames(ncols, 'col ');
end
% Init basic quant strings, if needed
if nRowNames ~= nrow
    if nRowNames > 0
        warning('fewer row names than matrix rows')
    end
    RowNames = basicNames(nrow, 'row ');
end

% fname = sprintf('%s/sensitivity_table.tex',outputDir);
fname = sprintf(filename,outputDir);

latextable(matrix, 'Horiz', ColNames, 'Vert', RowNames,...
    'Hline', [0:nRowNames,NaN], 'Vline', [0:nColNames,NaN],...
    'name', fname, 'format', '%.2g');

end
function Names=basicNames(n,head)
Names=cell(n,1);
for i=1:n
    Names{i}=sprintf('%s %i',head,i);
end
end