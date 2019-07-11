#!/bin/bash

while IFS= read -r line
do
 # echo $line
 ./main.exe $line 1

done < $1
