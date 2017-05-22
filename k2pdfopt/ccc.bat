@ECHO OFF
setlocal enableextensions enabledelayedexpansion
IF %1=="" GOTO EXITEND

for %%x in (%*) do (
    SET "LongFile=%%~x"
    ECHO !LongFile!
)

pause
