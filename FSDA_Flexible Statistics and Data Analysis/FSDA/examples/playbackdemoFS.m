function playbackdemoFS(demo_name)
%playbackdemoFS opens FSDA web site and launches playback device
%
%  PLAYBACKDEMOFS(DEMO_NAME) launches a playback demo in "www.riani.it/MATLAB/".
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit
%
%   Example:
%      playbackdemoFS('fishery')

if nargin <1
    url = 'www.riani.it/MATLAB/index.html';
else
    
    mlroot = 'www.riani.it/MATLAB/';
    
    % Assemble file paths
    url = [mlroot,demo_name,'/index.html'];
end
% Launch browser
web(url,'-browser');

end

