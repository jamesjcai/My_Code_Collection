function xlsreadmyown_ver2007(filename,sheetname,nonemptyidx,conv2logical)
if nargin<4
    conv2logical=false;
end

if isempty(findstr(filename,'.'))
    filename=strcat(filename, '.xlsx');
end
if nargin<2,
    [a,descr] = xlsfinfo(filename);
    descr
    k=input('Which sheet? ');    
    sheetname=descr{k};
end
if nargin<3
    nonemptyidx=[];
end


[A,b]=xlsread(filename,sheetname);
varname=b(1,:);
b=b(2:end,:);

if ~isempty(nonemptyidx)
    b=b(~isnan(A(:,nonemptyidx)),:);
    A=A(~isnan(A(:,nonemptyidx)),:);
end
    

idx=find(sum(strcmp(b,''))>0.5*size(b,1)); 
%idx=find(sum(strcmp(b,''))==size(b,1));    numeric fields
idx2=find(sum(isnan(A))==size(A,1));        % full nan column

%idx2
%idx
A(:,idx2)=[];
%if size(A,2)~=length(idx)
%    error('xxx')
%end


if isempty(idx)
for k=1:length(varname)
    try
        %varname{idx(k)}
        %exist(varname{idx(k)},'var')
 %       if ~ismember('chrid',existingvar)
            if isvarname(varname{k})
                assignin('base',varname{k},A(:,k));
            else
                fprintf('Invalid name: %s\n',varname{idx(k)});
                warning('Invalid name %s\n',varname{idx(k)});
            end                
                catch ME2
        rethrow(ME2);
        disp(varname)
    end    
end
return;
end


%existingvar=who;
for k=1:length(idx)    
    try
        %varname{idx(k)}
        %exist(varname{idx(k)},'var')
 %       if ~ismember('chrid',existingvar)
            if isvarname(varname{idx(k)})
                varA=A(:,k);                
                    if conv2logical&&isnumeric(varA)                    
                        ub=unique(varA);
                        if length(ub)==2
                            if ub(1)==0&&ub(2)==1
                                varA=logical(varA);                            
                            end
                        end
                    end
                assignin('base',varname{idx(k)},varA);
            else
                fprintf('Invalid name: %s\n',varname{idx(k)});
                warning('Invalid name %s\n',varname{idx(k)});
            end
  %      else
  %          warning(sprintf('Overwritten %s\n',varname{idx(k)}));
  %      end        
    catch ME0
        rethrow(ME0);
        disp(varname{idx(k)})
    end
end


idx2=setdiff(1:size(b,2),idx);
for k=1:length(idx2)    
    try
   %     if  ~exist(varname{idx(k)},'var')
            if isvarname(varname{idx2(k)})
            assignin('base',varname{idx2(k)},b(:,idx2(k)));
            else
                fprintf('2 Invalid name: %s\n',varname{idx2(k)});
                %warning('2 Invalid name %s\n',varname{idx2(k)});
            end
    %    else
    %        warning(sprintf('Overwritten %s\n',varname{idx(k)}));
    %    end        
    catch ME1
        rethrow(ME1);
        disp(varname{idx2(k)})
    end
end


    