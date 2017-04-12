function [genodata2,markinfo2]=snp_pickmarker(genodata,markinfo,s)
%[genodata2,markinfo2]=snp_pickmarker(genodata,markinfo,idx)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 18:40:46 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 526 $
% $LastChangedBy: jcai $

genodata2=genodata;
markinfo2=[];
if ~isempty(markinfo)
    markinfo2=markinfo;
end

if nargin<3    
[s,v] = choosebox('Name','Pick including marker(s)','PromptString',...
    'Markers available:','SelectString','Selected markers:',...
    'ListString',markinfo.rsid'); 
else
    v=1;
end
if (v==1),
    if islogical(s)
        s=find(s);
    end
    
    s2=s*2-1;
    sx=zeros(1,length(s)*2);
    for k=1:length(s)
        sx(k*2-1)=s2(k);
        sx(k*2)=s2(k)+1;
    end
    %sx
    genodata2 = genodata(:,sx); 
if ~isempty(markinfo)    
    a=fieldnames(markinfo);
    for k=1:length(a)
        if strcmp(a{k},'popid')
            markinfo2.popid=markinfo.popid;
        else
            try
                markinfo2=i_saftmapping(markinfo2,markinfo,a{k},s);
            catch ME
                %warning('Not all MARKINFO fields are extracted.')
            end
        end
    end
end

end

function markinfo2=i_saftmapping(markinfo2,markinfo1,fieldtxt,s)
    if isfield(markinfo2,fieldtxt)&&isfield(markinfo1,fieldtxt)
        v=getfield(markinfo1,fieldtxt,{s});        
        %markinfo2=setfield(markinfo2,fieldtxt,v);
        markinfo2.(fieldtxt)=v;
    end

