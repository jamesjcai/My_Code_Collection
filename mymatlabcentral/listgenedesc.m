function listgenedesc(genenames)
olddir=pwd;
cd 'C:/biodata/GENES/HGNC'
load HGNCGeneMap;
fprintf('\n====================================\n');
for k=1:length(genenames)
    a=genenames{k};
    if isKey(HGNCMap,a)
        fprintf('%s\t%s\n',a,HGNCMap(a));
    else
        fprintf('%s\t___\n',a);
    end
end
fprintf('====================================\n\n');
cd(olddir);