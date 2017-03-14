#!/bin/bash
# Compilation script, to be launched from the source folder
# Render the file executable if needed : chmod +x compilation.sh
cmakeFile="CMakeLists.txt"
makefile="Makefile"
execFile="sph"
if [ -f "$cmakeFile" ]
then
	if [ -f $execFile ]
	then
		rm $execFile
	fi
	echo "Building Cmake."
	cd build
	cmake ..
	if [ -f "$makefile" ]
	then
		echo "Compiling."
		make
		if [ -f "$execFile" ]
		then
			echo "Compilation succeeded."
			mv $execFile ../
		else
			echo "Compilation failed"
		fi
	else
		echo "Cmake Makefile generation failed."
	fi
		
else
	echo "$cmakeFile not found."
fi


