#!/bin/bash

# check if $1 is a file
if [ -f $1 ]; then
    # check if $1 is a .c file
    if [[ $1 == *.c ]]; then
        # compile $1 to assembly
        gcc -S $1

        # get the name of the file without the extension
        name=$(echo $1 | cut -d'.' -f1)

        # check if $name.s exists
        if [ -f $name.s ]; then
            cat $name.s
        else
            echo "Error: $name.s does not exist"
        fi
    else
        echo "Error: $1 is not a .c file"
    fi
else
    echo "Error: $1 is not a file"
fi