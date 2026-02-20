#!/bin/bash


# Get the input filename from the first argument
input_file="$1"
outputpath="$2"

#Compile data_decode.c
g++ -o data_decode_v2.exe data_decode_v2.c `root-config --cflags --libs`

#Run data_decode.exe based on the number of arguments passed
if [ -z "$1" ]; then
    #no file given
    ./data_decode_v2.exe
fi

if [ -n "$2" ]; then
    #path not empty
    ./data_decode_v2.exe "${input_file}" "${outputpath}"
else
    #path empty
    ./data_decode_v2.exe "${input_file}"
fi
