function [idx]=strfindcell(X,a,extractmatch)
%[idx]=strfindcell(X,a)

if nargin < 1 || isempty(X)
    error('strfindcell:TooFewInputs', ...
          'Requires a cell X.');
end
if nargin < 2 || isempty(a)
    error('strfindcell:TooFewInputs', ...
          'Requires a string a.');
end

%pnames = {'extractmatch','type'  'rows' 'tail'};
%dflts  = {fales, 'p'     'a'    'both'};

if nargin<3
    extractmatch=true;
end

if extractmatch
    %idx=cellfun(@(u) ~isempty( strmatch (u, a) ), X);
    idx=ismember(X,a);
    % find(ismember(rsidlist,Irsid{ky-1}))
    %({'rs1271','rs12718444'},'rs12718444')
else
    s=strfind(X,a);
    idx=~cellfun(@isempty,s);
    %warning('strfind returns partial match');
end

