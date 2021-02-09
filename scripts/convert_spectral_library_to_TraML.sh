#!/bin/bash

echo "Converting spectral libraries to .TraML format."
for filename in ecoli**.tsv
		do
				TargetedFileConverter -in $filename -out "${filename::-4}.TraML"
				echo "Generated ${filename::-4}.TraML"
		done;


