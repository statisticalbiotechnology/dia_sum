#!/bin/bash

echo "PyProhet merge..."

pyprophet score --in=merged.osw --level=ms1 score --in=merged.osw --level=ms2 score --in=merged.osw --level=transition --out merged.scored.osw



