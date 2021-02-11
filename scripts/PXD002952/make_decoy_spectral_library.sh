#!/bin/bash

echo "Generating spectral libraries with decoy."
for filename in ecoli**.TraML
		do
				OpenSwathDecoyGenerator -in $filename -out "${filename::-6}_target_decoy.TraML" -method pseudo-reverse
				echo "Generated ${filename::-6}_target_decoy.TraML"
		done;


