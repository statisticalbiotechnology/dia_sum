#!/bin/bash

echo "Generating .tsv from .txt spectral library files."
for filename in ecoli**.txt
		do
				cp $filename "${filename::-4}.tsv"
				echo "Generated ${filename::-4}.tsv"
		done;



