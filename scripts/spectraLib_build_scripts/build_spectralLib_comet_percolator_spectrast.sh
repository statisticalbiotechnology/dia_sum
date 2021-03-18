#!/bin/bash

# Run comet library search 
~/tools/crux/crux-3.2.Linux.x86_64/bin/crux comet --decoy_search 1 HYE124_TTOF6600_DDA_1ug_Ecoli_lgillet_I150212_002-Pedro_-_Ecoli_Library_-_IDA20_-_1ug_-_Repl2.mzML ecoli_UP000000625_83333.fasta

# Run Percolator 
~/tools/crux/crux-3.2.Linux.x86_64/bin/crux percolator --pepxml-output T crux-output/comet.target.pep.xml

# Copy .mzML to the filename format that spectraST wants
cp HYE124_TTOF6600_DDA_1ug_Ecoli_lgillet_I150212_002-Pedro_-_Ecoli_Library_-_IDA20_-_1ug_-_Repl2.mzML crux-output/comet.target.pep.mzML

# msconvert mzML to mzXML
msconvert --mzXML crux-output/comet.target.pep.mzML

cp 'HYE124_TTOF6600_DDA_1ug_Ecoli_lgillet_I150212_002-Pedro - Ecoli Library - IDA20 - 1ug - Repl2.mzXML' crux-output/comet.target.pep.mzXML


# Run spectraST
spectrast -cNraw -cP0.9 crux-output/percolator.target.pep.xml
