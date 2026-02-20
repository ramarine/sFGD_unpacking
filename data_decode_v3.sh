#!/bin/bash

input_file="$1"
outputpath="$2"

# Use .cpp extension — unambiguous for clang++ on macOS arm64
g++ -std=c++14 -o data_decode_v3.exe data_decode_v3.cpp `root-config --cflags --libs`

# Stop here if compilation failed
if [ $? -ne 0 ]; then
    echo "Compilation failed, aborting."
    exit 1
fi

if [ -z "$1" ]; then
    ./data_decode_v3.exe
elif [ -n "$2" ]; then
    ./data_decode_v3.exe "${input_file}" "${outputpath}"
else
    ./data_decode_v3.exe "${input_file}"
fi
