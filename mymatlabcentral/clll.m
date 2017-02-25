%function cll

x=questdlg('Close Clear all?');
if strcmp(x,'Yes')
    close all
    clear all
    clc
end