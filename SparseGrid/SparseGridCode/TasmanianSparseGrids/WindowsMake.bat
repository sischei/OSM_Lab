if /i "%1" == "clean" GOTO Clean

ECHO Compiling

cl -c *.cpp /Ox /EHsc /openmp

lib tsg*.obj Tas*.obj /OUT:libtasmaniansparsegrid.lib

cl tasgrid*.obj libtasmaniansparsegrid.lib /Fe:tasgrid.exe

GOTO End

:Clean

ECHO Cleaning

del *.obj *.dll *.lib *.exe

GOTO End

:End