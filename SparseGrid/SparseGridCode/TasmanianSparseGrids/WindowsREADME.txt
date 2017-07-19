In order to compile TASMANIAN using VC++

1. Install Microsoft Visual C++

2. Extract the .zip files, for example into C:\User\Bill\Downloads\Tasmanian

3. Run "Development Command Prompt VS..." from the windows menu

4. Navigate to the folder with the source code using the cd command, in the example above: cd C:\User\Bill\Downloads\Tasmanian

5. Run the WindowsMake.bat script, i.e., type: WindowsMake.bat

The WindowsMake scrip will compile the static library and the tasgrid.exe executable.

Note: to delete all files created by the compile use: "WindowsMake.bat clean"

Note: to if you are experiencing problems with OpenMP, you can disable it by: Right-Click on WindowsMake.bat, select Edit, remove the /openmp directive from the third command, call WindowsMake.bat clean, followed by WindowsMake.bat