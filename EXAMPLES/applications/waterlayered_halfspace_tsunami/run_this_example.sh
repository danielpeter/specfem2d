#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

# SPECFEM2D root directory
ROOT=../../../

# links executables
mkdir -p bin
cd bin/
rm -f x*
for f in ../$ROOT/bin/x*; do
  #echo "link executable: $f"
  ln -s $f
done
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
echo
echo "running mesher on $NPROC processors..."
echo
mpirun -np $NPROC ./bin/xmeshfem2D | tee OUTPUT_FILES/output_mesher.txt
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
echo
echo "running solver on $NPROC processors..."
echo
mpirun -np $NPROC ./bin/xspecfem2D | tee OUTPUT_FILES/output_solver.txt
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

