function [status]=saveSeuratRdsOLD(sce,filename)
    [status]=0;
    isdebug=true;
    if nargin<2, error('run.r_saveSeuratRds(sce,filename)'); end
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_SeuratSaveRds');
    if ~isok, error(msg); return; end
    if exist('output.Rds','file'), delete('output.Rds'); end
    sc_writefile('input.txt',sce.X,sce.g);
    %if isdebug, return; end
    Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

    [status]=copyfile('output.Rds',filename,'f');
    %if exist('input.txt','file'), delete('input.txt'); end
    if exist('output.Rds','file'), delete('output.Rds'); end
    cd(oldpth);
end