function [ss]=waitforbuttonpressFS
%waitforbuttonpressFS produces a warning message when user closes the figure he is interacting with 
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit
%
%

%% Beginning of code
try
    ss=1;
    ss=waitforbuttonpress;
catch ME
    if strcmp(ME.message,'waitforbuttonpress exit because all figures have been deleted')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        disp('The figure you are interacting with has been deleted: PROCESS ENDED');
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        ss=1;
    end
end
end