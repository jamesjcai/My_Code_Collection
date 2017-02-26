function [standardname]=yeastgenename(name)
%view-source:http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=ACE1
standardname='';
%name='ACE1';
urlFetch=sprintf('http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=%s',...
    name);
try
    pagecontent=urlread(urlFetch);
catch E01
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(E01);
end

fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);
k=1;
c='';
while isempty(c)&&k<=length(fetchResults);
    a=fetchResults{k};
    c=i_getit(a);
    k=k+1;
end
standardname=c;
    



function c=i_getit(a)
if  strfind(a,'<span class="gene_name">')==1
    b=strfind(a,'</');
    c=a(25:b(1)-1);
else
    c='';
end
    