#!/bin/bash

../bin/preprocess 2REX/2REX_r_u.pdb -c AB -u AMSE:BMSE -o 2REX/2REXR.pdb -a ../data/ATM
../bin/preprocess 2REX/2REX_l_u.pdb -c D -s @GTP -o 2REX/2REXL.pdb -a ../data/ATM
../bin/watmap 2REX/2REXR.pdb 2REX/2REXRwat.pdb -a ../data/ATM
../bin/watmap 2REX/2REXL.pdb 2REX/2REXLwat.pdb -a ../data/ATM
../bin/sdock 2REX/2REXR.pdb 2REX/2REXL.pdb -o 2REX/2REX_record -r ../data/so3layer.qua -p 2.0
../bin/build 2REX/2REXR.pdb 2REX/2REXL.pdb -o 2REX/2REX_record -r ../data/so3layer.qua -c 1 -m 1 -n 2REX -d 2REX/
../bin/sdock 2REX/2REXR.pdb 2REX/2REXL.pdb 2REX/2REXRwat.pdb 2REX/2REXLwat.pdb -o 2REX/2REX_watrecord -r ../data/so3layer.qua -p 0
../bin/build 2REX/2REXR.pdb 2REX/2REXL.pdb -o 2REX/2REX_watrecord -r ../data/so3layer.qua -c 1 -m 1 -n 2rex -d 2REX/
