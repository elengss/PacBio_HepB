#!/bin/bash
echo "Find full path lengths for use in shell script"
echo ""
 
full_path=$(realpath $0)
home=$(dirname $full_path)
echo "If this is the top level directory, then the entry for home="
echo "in the shell scripts is :"
echo $home
echo "with the following subfolders :"
echo $home/data
echo $home/reference
