#!/bin/bash


# Set the directories
# Directories removed prior to uploading to GitHub

home=
inputs=
outputs=
mat_files=${inputs}matrix_files/

# Initialize
echo "Extraction pipeline started at " $(date)
echo "Getting a list of all of the matrix files to run"
files=`cat ${inputs}allen_files_1.txt`
echo ${files}


echo "Getting a list of all queries to run"
cd ${inputs}
queries=$(ls search_split*)
cd ${home}
echo ${queries}



# Extract each cell and corresponding expression data
# corresponding to distinct rows in the df defined by unique cell ID. 
# Due to memory limitations, divided the search queries and single-cell df into smaller chunks
# using slice, then iterated through each before cat'ing the files back together to create the 
# final df

echo "Extracting cells"
head -n 1 ${mat_files}counts_split00 > ${outputs}19850212-2500_all_subclass_exprs_1.txt

for i in ${files[@]}
  do
    echo "Extracting: " ${i}
    for x in ${queries}
      do
        echo "Grep using " ${x}
        grep -f ${inputs}${x} ${mat_files}${i} >> ${outputs}19850212-2500_all_subclass_exprs_1.txt
        echo "Done"
      done
      echo "Extraction for " ${i} " complete at" $(date)
  done
echo "Extraction pipeline complete at " $(date)


# To run with logs
# nohup subset_all_expanded_1.sh >> subset_all_expanded_1.sh.log 2>> subset_all_expanded_1.sh.error &
