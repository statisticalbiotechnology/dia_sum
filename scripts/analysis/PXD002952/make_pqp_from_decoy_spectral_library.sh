#!/bin/bash

echo "Generating spectral libraries with decoy."
for filename in ecoli**decoy.TraML
		do
			    #echo $filename
				TargetedFileConverter -in ${filename} -out ${filename::-6}.pqp
				echo "Generated ${filename::-6}_target_decoy.pqp"
		done;


