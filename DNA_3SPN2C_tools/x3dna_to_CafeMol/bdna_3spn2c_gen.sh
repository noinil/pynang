#!/usr/bin/env bash

PYNANG_BIN_PATH=/home/ctan/Workspace/pynang/DNA_3SPN2C_tools/x3dna_to_CafeMol
export PATH=$PATH:/home/ctan/Workspace/x3dna-v2.4/bin
export X3DNA=/home/ctan/Workspace/x3dna-v2.4

# This tool is originally built by de Pablo's group and modified to generate
# input files for CafeMol.

if [ $# -ne 1 ]; then
    echo "Usage: $0 <sequence file>"
    exit 1
fi

echo "================================================================================"
echo "Making DNA curvature parameter file..."
$PYNANG_BIN_PATH/seq2curv_DNA2c.py $1

echo "--------------------------------------------------------------------------------"
echo "Running X3DNA..."
x3dna_utils cp_std BDNA
rebuild -atomic dna2c.curv bdna_aa.pdb
rm -f Atomic*
rm -f ref_frames.dat
rm -f dna2c.curv

echo "--------------------------------------------------------------------------------"
echo "Making CafeMol ninfo files..."
$PYNANG_BIN_PATH/pdb2ninfo_DNA2c.py bdna_aa.pdb
echo " DONE!"
echo "================================================================================"
