#!/bin/bash

# check if $1 is a file
if [ -f $1 ]; then
    # check if $1 is a .s file
    if [[ $1 == *.s ]]; then

        # get the name of the file without the extension
        name=$(echo $1 | cut -d'.' -f1)

        # compile $name.s to executable
        gcc -o $name.out $name.s
        # check if $name.out exists
        if [ -f $name.out ]; then
            # run $name.out
            ./$name.out
        else
            echo "Error: $name.out does not exist"
        fi
    else
        echo "Error: $1 is not a .s file"
    fi
else
    echo "Error: $1 is not a file"
fi