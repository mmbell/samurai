#!/bin/bash
path='/home/quinting/Diplomarbeit/Daten/Dropsondendaten_C130/'
filename='D20080920_054223_PQC.eol.Wwind'
#filename='D20080929_233718_PQC.eol.Wwind'
here=$PWD
cd $path
files=`ls`
#echo $files
cd $here

#for file in $files
for file in $filename
do
  echo EXECUTING:  ncl 'dir="'$path'"' 'fname="'$file'"' skewt_tparc_samurai.ncl
  ncl 'dir="'$path'"' 'fname="'$file'"' skewt_tparc_samurai.ncl
done