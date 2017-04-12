function [status] =  snp_writelinkage(geno,mark,filename,forplink)
%SNP_WRITELINKAGE - saves as linkage format
%snp_writelinkage(geno,mark,filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 18:40:46 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 526 $
% $LastChangedBy: jcai $

if nargin<4
    forplink=true;
end
if (isempty(geno)), status=0; return; end
if nargin<2
    mark=[];
end
%if (isempty(mark)), status=0; return; end

if nargin < 3
    [filename, pathname,filterindex] = uiputfile( ...
        {'*.pedigree;*.ped', 'Linkage Format Files (*.pedigree, *.ped)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
    
	if filterindex==1
        if isempty(find(filename=='.', 1))
            filenameped=[filename,'.ped'];
            filenamemap=[filename,'.map'];
        else
            filenameped=filename;
            filenamemap=[filename,'.map'];
        end
    end           
else
       if isempty(find(filename=='.', 1))
			filenameped=[filename,'.ped'];
            filenamemap=[filename,'.map'];
       else
            filenameped=filename;
            filenamemap=[filename,'.map'];
       end
end

fid = fopen(filenameped,'wt');
if fid == -1
   status=0;
   warning('Unable to open file.');
   return;
end

[samplen,marklen]=snp_samplen(geno);
indvlen=samplen/2;

if ~isempty(mark)
if ~isfield(mark,'rsid')
    mark.rsid=num2cellstr(1:marklen);
end
end


%ACGT='12340';
ACGT='ACGT0';

for k=1:indvlen
    
%fprintf(fid,['%d\n'],indvlen);
%fprintf(fid,['%d\n'],marklen);
%fprintf(fid,['P %s\n'],sprintf('%d ',mark.pos));
%fprintf(fid,[char(ones(1,marklen)*['S']),'\n']);
      fprintf(fid,'%d %d 0 0 1 1 ',k,1);      
      %for j=1:marklen*2
	  %    fprintf(fid,'%s\t',ACGT(geno(k,j)));
      %end
      fprintf(fid,'%c ',ACGT(geno(k,1:end-1)));
      fprintf(fid,'%c\n',ACGT(geno(k,end)));
end
fclose(fid);

if ~isempty(mark)
fid = fopen(filenamemap,'wt');
    %fprintf(fid,'MARKER_ID\tSNP_rs#\tbp_POSITION\n');
for k=1:marklen
    if mark.rsid{k}=='.'
        mark.rsid{k}=sprintf('%d_%d',mark.chrid(k),mark.pos(k));
    end
    if forplink
        try
        fprintf(fid,'%d\t',mark.chrid(k));   % chromosome (1-22, X, Y or 0 if unplaced)
        catch
            if isfield(mark,'chrid')
                fprintf(fid,'%d\t',mark.chrid);   % chromosome (1-22, X, Y or 0 if unplaced)    
            else
                fprintf(fid,'%d\t',0);   % chromosome (1-22, X, Y or 0 if unplaced)    
            end
        end
        fprintf(fid,'%s\t',mark.rsid{k}); % rs# or snp identifier
        fprintf(fid,'%d\t',0);   % Genetic distance (morgans)
        fprintf(fid,'%d\n',mark.pos(k));  % Base-pair position (bp units)
    else
        fprintf(fid,'%s\t',mark.rsid{k});
        fprintf(fid,'%d\n',mark.pos(k));
    end
    
end
fclose(fid);
end
status=1;


