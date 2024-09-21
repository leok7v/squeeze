::: amalgamate.bat
@echo off
if exist ..\shl (
    if not exist ..\shl\squeeze (
        mkdir ..\shl\squeeze
    )
    (
        type ..\inc\squeeze\squeeze.h
        ::: LF
        echo.
        echo #ifdef squeeze_implementation
        type ..\src\squeeze.c
        echo.
        echo #endif // squeeze_implementation
        echo.
    ) > ..\shl\squeeze\squeeze.h
)
