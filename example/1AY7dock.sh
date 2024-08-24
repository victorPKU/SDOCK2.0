#!/bin/bash

../bin/preprocess 1AY7/1AY7_r_u.pdb -c A -o 1AY7/1AY7R.pdb -a ../data/ATM
../bin/preprocess 1AY7/1AY7_l_u.pdb -c B -o 1AY7/1AY7L.pdb -a ../data/ATM
../bin/watmap 1AY7/1AY7R.pdb 1AY7/1AY7Rwat.pdb -a ../data/ATM
../bin/watmap 1AY7/1AY7L.pdb 1AY7/1AY7Lwat.pdb -a ../data/ATM
../bin/sdock 1AY7/1AY7R.pdb 1AY7/1AY7L.pdb -o 1AY7/1AY7_record -r ../data/so3layer.qua
../bin/build 1AY7/1AY7R.pdb 1AY7/1AY7L.pdb -o 1AY7/1AY7_record -r ../data/so3layer.qua -c 66 -m 1 -n 1AY7 -d 1AY7/model/
../bin/sdock 1AY7/1AY7R.pdb 1AY7/1AY7L.pdb 1AY7/1AY7Rwat.pdb 1AY7/1AY7Lwat.pdb -o 1AY7/1AY7_watrecord -r ../data/so3layer.qua
../bin/build 1AY7/1AY7R.pdb 1AY7/1AY7L.pdb -o 1AY7/1AY7_watrecord -r ../data/so3layer.qua -c 66 -m 1 -n 1ay7 -d 1AY7/model/

