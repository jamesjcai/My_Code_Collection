%function cll

x=questdlg('Close Clear all?');
if strcmp(x,'Yes')
    close all force
    % clear all
    % clearvars
    clear
end