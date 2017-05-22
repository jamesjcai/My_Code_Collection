@ECHO OFF
setlocal enabledelayedexpansion
COLOR 17
ECHO(
IF %1=="" GOTO EXITEND

SET argC=0
FOR %%x in (%*) DO SET /A argC+=1

SET HOMEDIR=%USERPROFILE%\OneDrive\k2pdfopt
SET TEMAIL="jamesjjcai@kindle.com"

SET LongFile=%~n1
SET ShortFile=%LongFile:~0,80%%~x1
SET TFILE="%TMP%\%ShortFile%"
SET CFILE="%TMP%\c_%ShortFile%"
IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d %CFILE% -s %1
%HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o %TFILE% %CFILE%
%HOMEDIR%\MailAlert\MailAlert.exe -a %TFILE% -r %TEMAIL%
MOVE %TFILE% "%HOMEDIR%\pdfpool\processed\"
MOVE %1 "%HOMEDIR%\pdfpool\original\%ShortFile%"


IF %argC% GTR 1 (
    SET LongFile=%~n2
    SET ShortFile=!LongFile:~0,80!%~x2
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %2
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %2 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)


IF %argC% GTR 2 (
    SET LongFile=%~n3
    SET ShortFile=!LongFile:~0,80!%~x3
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %3
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %3 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)


IF %argC% GTR 3 (
    SET LongFile=%~n4
    SET ShortFile=!LongFile:~0,80!%~x4
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %4
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %4 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)

IF %argC% GTR 4 (
    SET LongFile=%~n5
    SET ShortFile=!LongFile:~0,80!%~x5
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %5
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %5 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)

IF %argC% GTR 5 (
    SET LongFile=%~n6
    SET ShortFile=!LongFile:~0,80!%~x6
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %6
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %6 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)

IF %argC% GTR 6 (
    SET LongFile=%~n7
    SET ShortFile=!LongFile:~0,80!%~x7
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %7
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %7 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)

IF %argC% GTR 7 (
    SET LongFile=%~n8
    SET ShortFile=!LongFile:~0,80!%~x8
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %8
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %8 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)

IF %argC% GTR 8 (
    SET LongFile=%~n9
    SET ShortFile=!LongFile:~0,80!%~x9
    SET TFILE="%TMP%\!ShortFile!"
    SET CFILE="%TMP%\c_!ShortFile!"
    IF NOT EXIST !CFILE! JAVA -jar %HOMEDIR%\briss\briss-0.9.jar -d !CFILE! -s %9
    %HOMEDIR%\k2pdfopt.exe -dev kv -fs 12 -fc- -o !TFILE! !CFILE!
    %HOMEDIR%\MailAlert\MailAlert.exe -a !TFILE! -r %TEMAIL%
    MOVE !TFILE! "%HOMEDIR%\pdfpool\processed\"
    MOVE %9 "%HOMEDIR%\pdfpool\original\!ShortFile!"
)

IF %argC% GTR 9 (
    ECHO Only up to 9 files can be processed at a time.
)

endlocal
GOTO EXITEND

:EXITEND
ECHO 
TIMEOUT 5 > NUL
