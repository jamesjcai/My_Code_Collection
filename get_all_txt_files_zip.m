
% files=ReadFileNames('CIBERSORT_COMMONMIND',{'txt','csv','tsv'});
d = uigetdir(pwd, 'Select a folder');
files=ReadFileNames(d,{'txt','csv','tsv','fa','mtx'});
for k=1:length(files)
    f=files{k};
    s=dir(f);
    if s.bytes >= 2000000
        try
            fz=gzip(f);
            s2=dir(fz{1});
            fprintf('%s......%.2f\n',f,s2.bytes/s.bytes);
            delete(f);
        catch
            fprintf('%s......ERROR\n',f);
        end
    end
end


function [ FList ] = ReadFileNames(DataFolder,extList)
% Author: Thokare Nitin D.
% 
% This function reads all file names contained in Datafolder and it's subfolders
% with extension given in extList variable in this code...
% Note: Keep each extension in extension list with length 3
% i.e. last 3 characters of the filename with extension
% if extension is 2 character length (e.g. MA for mathematica ascii file), use '.'
% (i.e. '.MA' for given example)
% Example:
% extList={'jpg','peg','bmp','tif','iff','png','gif','ppm','pgm','pbm','pmn','xcf'};
% Gives the list of all image files in DataFolder and it's subfolder
% 

if nargin<2
    extList={'txt','csv'};
else
    % extList
end

DirContents=dir(DataFolder);
FList=[];

if(strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64'))
    NameSeperator='\';
elseif(strcmpi(computer,'GLNX86') || strcmpi(computer,'GLNXA86'))
    NameSeperator='/';
end

% extList={'jpg','peg','bmp','tif','iff','png','gif','ppm','pgm','pbm','pmn','xcf'};


% Here 'peg' is written for .jpeg and 'iff' is written for .tiff
for i=1:numel(DirContents)
    if(~(strcmpi(DirContents(i).name,'.') || strcmpi(DirContents(i).name,'..')))
        if(~DirContents(i).isdir)            
            if ~isempty(strfind(DirContents(i).name,'.'))
                idx=strfind(DirContents(i).name,'.');
                extension=DirContents(i).name(idx(end)+1:end);
            % extension=DirContents(i).name(end-2:end);
            if(numel(find(strcmpi(extension,extList)))~=0)
                FList=cat(1,FList,{[DataFolder,NameSeperator,DirContents(i).name]});
            end
            end
        else
            getlist=ReadFileNames([DataFolder,NameSeperator,DirContents(i).name],extList);
            FList=cat(1,FList,getlist);
        end
    end
end

end


