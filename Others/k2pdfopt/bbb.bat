:MENU
@echo off
COLOR 17
Echo(
Echo 1. Send to Kindle Keyboard
Echo 2. Send to Kindle Oasis
Echo 3. Exit
Echo(
SET HOMEDIR=%USERPROFILE%\OneDrive\k2pdfopt
REM SET INFILE="%~dpnx1"
SET TFILE1="%TMP%\%~nx1"
SET TFILE2="%TMP%\chopped_%~nx1"
SET /P Choice= Please make a choice from the menu above: 


IF NOT EXIST %TFILE2% JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d %TFILE2% -s %1

IF %Choice%==1 GOTO K1
IF %Choice%==2 GOTO K2
IF %Choice%==3 GOTO EXITEND

:K1
%HOMEDIR%\k2pdfopt.exe -fs 12 -fc- -o %TFILE1% %TFILE2%
%HOMEDIR%\MailAlert\MailAlert.exe -a %TFILE1% -r "jamesjcai@kindle.com"
GOTO EXITEND

:K2
%HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o %TFILE1% %TFILE2%
%HOMEDIR%\MailAlert\MailAlert.exe -a %TFILE1% -r "jamesjjcai@kindle.com"
GOTO EXITEND

:EXITEND
ECHO 
TIMEOUT 5 > NUL
