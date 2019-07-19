a=dir('../*.m');

for k=1:length(a)
   f=fullfile(a(k).folder, a(k).name);
[~, toolboxList] = matlab.codetools.requiredFilesAndProducts(f);
    for kk = 1 : length(toolboxList)        
        fprintf('%s   %s\n', a(k).name, toolboxList(kk).Name);
    end   
end
