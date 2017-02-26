function xlsreadmyown(filename,sheetname)

if isempty(findstr(filename,'.'))
    filename=strcat(filename, '.xls');
end

if nargin<2,
    [a,descr] = xlsfinfo(filename);
    descr
    k=input('Which sheet? ');    
    sheetname=descr{k};
end



[A,b]=xlsread(filename,sheetname);
varname=b(1,:);
b=b(2:end,:);

if size(A,1)~=0
    idx=find(sum(cellfun('isempty',b))==size(A,1));
    idx2=find(sum(cellfun('isempty',b))<size(A,1));
    A(:,find(sum(isnan(A))==size(A,1)))=[];
else
    idx=[];
    idx2=find(sum(cellfun('isempty',b))<=size(b,1));
end

%A=[zeros(size(A,1),size(b,2)-size(A,2)),A];

%idx=find(sum(strcmp(b,''))==size(b,1));
%idx2=find(sum(isnan(A))==size(A,1));
%idx2

%A(:,idx2)=[];

%if size(A,2)~=length(idx)
%    error('xxx')
%end

%existingvar=who;


for k=1:length(idx)    
    try
        %varname{idx(k)}
        %exist(varname{idx(k)},'var')
 %       if ~ismember('chrid',existingvar)
            if isvarname(varname{idx(k)})
                if size(A,1)~=0
                    assignin('base',varname{idx(k)},A(:,find(idx==idx(k))));
                else
                    assignin('base',varname{idx(k)},A(:,idx(k)));
                end
            else
            disp(sprintf('Invalid name: %s\n',varname{idx(k)}));
            warning(sprintf('Invalid name %s\n',varname{idx(k)}));                
            end
  %      else
  %          warning(sprintf('Overwritten %s\n',varname{idx(k)}));
  %      end        
    catch
        disp(varname{idx(k)})
    end
end


%idx2=setdiff(1:size(b,2),idx);
for k=1:length(idx2)
    try
   %     if  ~exist(varname{idx(k)},'var')
            if isvarname(varname{idx2(k)})
            assignin('base',varname{idx2(k)},b(:,idx2(k)));
            else
                disp(sprintf('2 Invalid name: %s\n',varname{idx2(k)}));
                warning(sprintf('2 Invalid name %s\n',varname{idx2(k)}));                                
            end
    %    else
    %        warning(sprintf('Overwritten %s\n',varname{idx(k)}));
    %    end        
    catch ME1
        rethrow(ME1);
        disp(varname{idx2(k)})
    end
end


    