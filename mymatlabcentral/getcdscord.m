function [transid,transcord]=getcdscord(id)
%[transid,transcord,DX]=getcdscord('ENSG00000183431')

%{
load('y:/humangenes/humangenelist','D');
[ok,idx]=ismember(id,D.geneid);
if ~ok
    transid=[];transcord=[];DX=[];
    return;
end
DX.strand=D.strand(idx);
DX.chrid=D.chrid(idx);
%}


%id='ENSG00000183431';
filename=sprintf('Y:/GENES/humangenes/res_v51/%s.tab',id);
txt = textread(filename,'%s','delimiter','\n','whitespace','','bufsize',4095*2);
mt = find(cellfun('isempty',txt));
txt(mt) = [];
totalline=length(txt);
if (mod(totalline,4)>0),
    error('xx')
end

counter=1;
transid={};
for k=1:4:totalline
    transid{counter}=txt{k};
        s1=str2num(txt{k+1});
        e1=str2num(txt{k+3});
        a=txt{k+2}(1:end-1);
    transcord{counter}=i_makecoord(a,s1,e1);
    counter=counter+1;
end



function [data]=i_makecoord(a,s1,e1)
%%%%%%
[s] = strread(a,'%s','delimiter',',');
data=zeros(length(s),2);
for k=1:length(s)
    [data(k,1),data(k,2)] = strread(s{k},'%d%d','delimiter','-');
end
data=sortrows(data);
id1=inladder2(data,s1);
id2=inladder2(data,e1);
%id1=sum(cdspos1>=e(:,1));
%id2=sum(cdspos2>=e(:,1));
data=data([id1:id2],:);
data(1)=s1;
data(end)=e1;


%if (mod(sum(data(:,2)-data(:,1)+1),3)>0)
%	data(end)=e1-mod(sum(data(:,2)-data(:,1)+1),3);
%end



