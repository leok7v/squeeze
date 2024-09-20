::: amalgamate.bat
@echo off
if exist ..\shl (
    if not exist ..\shl\squeeze (
        mkdir ..\shl\squeeze
    )
    (
        findstr /v /c:"#endif // squeeze_h" ..\inc\squeeze\squeeze.h
        echo #ifdef squeeze_implementation
        type ..\src\squeeze.c
        echo.
        echo #endif // squeeze_implementation
        echo.
        echo #endif // squeeze_h
    ) > ..\shl\squeeze\squeeze.h
)
