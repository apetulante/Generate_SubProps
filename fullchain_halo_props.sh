#!/bin/bash

echo "Usage: create a table of halo properties for trees with a range of numbers, from tree_a_b_c.dat (input as 'abc') to tree_d_e_f.dat (input as 'def')"

tree_start1=$(echo "$1" | awk -F '_' '{print $2}' )
tree_start2=$(echo "$1" | awk -F '_' '{print $3}' )
tree_start3=$(echo "$1" | awk -F '_' '{print $4}' )
tree_end1=$(echo "$2" | awk -F '_' '{print $2}' )
tree_end2=$(echo "$2" | awk -F '_' '{print $3}' )
tree_end3=$(echo "$2" | awk -F '_' '{print $4}' )
mass_cut=1000

python get_tree_range.py $tree_start1 $tree_start2 $tree_start3 $tree_end1 $tree_end2 $tree_end3

while IFS='' read -r tree || [[ -n "$tree" ]]; do
    echo " Working on: $tree"

    cp /fs2/shared/petulaa/sorted_trees/${tree}.dat .
    
    treenum1=$(echo "$tree" | awk -F '_' '{print $2}' )
    treenum2=$(echo "$tree" | awk -F '_' '{print $3}' )
    treenum3=$(echo "$tree" | awk -F '_' '{print $4}' )    

    echo "Removing the header..."
    tail -n +51 ${tree}.dat > ${tree}.txt
    rm ${tree}.dat

    echo "Getting rid of comment lines..."
    sed '/^#/ d' < ${tree}.txt > temp_${tree}.txt    
    mv temp_${tree}.txt ${tree}.txt

    echo "Counting lines in file..."
    lines=`cat ${tree}.txt | wc -l`
   
    echo "Starting the splitter code"
    ./full_split_halos ${tree}.txt $lines $treenum1 $treenum2 $treenum3
    rm ${tree}.txt  

    python halo_properties_fullChain.py $tree

    cat halo_props_${tree}.txt >> halo_properties_${1}_to_${2}.txt
    cat halo_props_${tree}_survivors.txt >> halo_properties_${1}_to_${2}.txt    

    rm ${tree}_halo_1*
    rm ${tree}_halo_2*
    rm ${tree}_halo_3*
    rm ${tree}_halo_4*
    rm ${tree}_halo_5*
    rm ${tree}_halo_*
    rm ${tree}_subs_today.txt
    rm halo_props_${tree}*
done < "tree_range_tree_${tree_start1}_${tree_start2}_${tree_start3}_to_tree_${tree_end1}_${tree_end2}_${tree_end3}.txt" 
