:MENU
@echo off
COLOR 17
Echo(
SET HOMEDIR=%USERPROFILE%\OneDrive\k2pdfopt

%HOMEDIR%\MailAlert\MailAlert.exe -a %1 -r "jamesjcai@kindle.com"
GOTO EXITEND

:EXITEND
ECHO 
TIMEOUT 50 > NUL
