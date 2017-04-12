function [OUT]=readmsoutput(filename)
%READMSOUTPUT - reads Hudson's MS output
%
% Usage: [OUT]=readmsoutput(filename)
%
% OUT.commandline - e.g., 'ms 120 1000 -t 10 '
% OUT.segsites - e.g., {1x1000 cell}
% OUT.positions - e.g., {1x1000 cell}
% OUT.gametes - e.g., {1x1000 cell}

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


%disp(['Reading MS output file ',filename,' ...']);
%txt = textread(filename,'%s','delimiter','\n','whitespace','');

fid=fopen(filename,'r');
txt=textscan(fid,'%s','delimiter','\n','whitespace','');
txt=txt{1};
fclose(fid);


% eliminate empty lines
txt(cellfun('isempty',txt)) = [];

OUT.commandline=txt{1};
%a=sscanf(txt{1},'ms%d%d');
%nsam=a(1);
%OUT.nsam=nsam;
%OUT.nreps=a(2);

idx1=find(cellfun(@isempty,strfind(txt,'//'))==0);
n=length(idx1);

x=idx1+3;
y=idx1(2:end)-1;
y=[y; length(txt)];

    idx1=[idx1;length(txt)];
    for kk=1:n
        thistxt=txt(idx1(kk):idx1(kk+1));
        OUT.segsites{kk}=i_getfield1('segsites:',thistxt);
        OUT.positions{kk}=i_getfield1('positions:',thistxt);
        %OUT.probs{kk}=i_getfield1('probs:',thistxt);
        %OUT.times{kk}=i_getfield1('times:',thistxt);
        
        a=x(kk); b=y(kk);
        Gx=logical([]);
        for j=a:b        
            c=strtok(txt{j});
            g=false(1,length(c));
            for i=1:length(c)
                if c(i)=='1'
                    g(i)=true;
                end
            end
            Gx=[Gx;g];
        end
        OUT.gametes{kk}=Gx;        
    end

%{    
for k=1:n-1
        a=idx1(k)+3;
        b=a+nsam-1;
        Gx=[];
        
    if (idx1(k+1)-idx1(k))>4
    for j=a:b        
        c=strtok(txt{j});
        g=false(1,length(c));
        for i=1:length(c)
            if c(i)=='1'
                g(i)=true;
            end
            %g(i)=str2double(c(i));
        end
        Gx=[Gx;g];
    end
    end
    %gametes{k}=logical(Gx);
    gametes{k}=Gx;
end
OUT.gametes=gametes;
%}



function [X]=i_getfield1(word,txt)
    lword=length(word);
    idx2=find(cellfun(@isempty,strfind(txt,word))==0);
    if isempty(idx2)
        X=[];    
    else
        X=zeros(1,length(idx2));
        for k=1:length(idx2)
            linetxt=txt{idx2(k)};                        
            X=str2num(linetxt(lword+1:end));
            %segsites(k)=strread(linetxt(strfind(linetxt,':')+2:end),'%d');
        end
    end

function [X]=i_getfield2(word,txt)
    lword=length(word);
    idx2=find(cellfun(@isempty,strfind(txt,word))==0);
    if isempty(idx2)
        X={};    
    else        
        for k=1:length(idx2)
            linetxt=txt{idx2(k)};
            X{k}=str2double(linetxt(lword+1:end));
            %segsites(k)=strread(linetxt(strfind(linetxt,':')+2:end),'%d');
        end
    end


%{
$segsites
[1] 5 5 5 5

$gametes
$gametes[[1]]
     [,1] [,2] [,3] [,4] [,5]
[1,]    1    0    1    0    1
[2,]    1    0    0    0    1
[3,]    0    1    0    1    0
[4,]    0    1    0    1    0
[5,]    0    1    0    1    0

$gametes[[2]]
     [,1] [,2] [,3] [,4] [,5]
[1,]    1    1    0    0    0
[2,]    0    1    0    0    0
[3,]    0    1    0    0    0
[4,]    1    1    0    0    0
[5,]    0    0    1    1    1

$gametes[[3]]
     [,1] [,2] [,3] [,4] [,5]
[1,]    1    0    0    0    0
[2,]    0    0    0    1    0
[3,]    0    0    0    1    0
[4,]    0    0    0    1    0
[5,]    1    1    1    0    1

$gametes[[4]]
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    1    1    0    0
[2,]    1    0    0    0    1
[3,]    0    0    0    1    0
[4,]    0    1    1    0    0
[5,]    0    1    1    0    0


$probs
list()

$times

[1,]

$positions
$positions[[1]]
[1] 0.0480 0.1846 0.2939 0.5769 0.5924

$positions[[2]]
[1] 0.2031 0.2734 0.6270 0.7904 0.8603

$positions[[3]]
[1] 0.3335 0.3380 0.3796 0.4202 0.6175

$positions[[4]]
[1] 0.0750 0.0886 0.1544 0.8918 0.9935


$nsam
[1] 5

$nreps
[1] 4
%}