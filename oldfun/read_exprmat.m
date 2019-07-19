function [X,genelist,sampleid]=read_exprmat(filename,genecolnum,verbose)
% Validate input args
narginchk(1,Inf);

% Get Filename
if ~ischar(filename) && ~(isstring(filename) && isscalar(filename))
    error(message('MATLAB:csvread:FileNameMustBeString')); 
end
filename = char(filename);

% Make sure file exists
if exist(filename,'file') ~= 2 
    error(message('MATLAB:csvread:FileNotFound'));
end


if nargin<2, genecolnum=1; end
if nargin<3, verbose=true; end
if verbose, fprintf('Reading %s ...... ',filename); end
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
T=readtable(filename);
warning('on','MATLAB:table:ModifiedAndSavedVarnames');
X=table2array(T(:,1+genecolnum:end));
genelist=string(table2array(T(:,1:genecolnum)));
if verbose, fprintf('done.\n'); end
if nargin>2
    sampleid=string(T.Properties.VariableNames(1+genecolnum:end)');
end
