#!/bin/sh

Windname="Halo_Trees"

Trees_vals=("0_0_0" "0_6_0" "1_0_0" "1_6_0" "2_0_0" "2_6_0" "3_0_0" "3_6_0" "4_0_0" "4_6_0" "5_0_0" "5_6_0" "6_0_0" "6_6_0" "7_0_0" "7_6_0" "8_0_0" "8_6_0" "9_0_0" "9_6_0" "10_0_0" "10_6_0")
Trees_end=("0_5_10" "0_10_10" "1_5_10" "1_10_10" "2_5_10" "2_10_10" "3_5_10" "3_10_10" "4_5_10" "4_10_10" "5_5_10" "5_10_10" "6_5_10" "6_10_10" "7_5_10" "7_10_10" "8_5_10" "8_10_10" "9_5_10" "9_10_10" "10_5_10" "10_10_10" ) 

# Initialize screen session
screen -mdS ${Windname}

## Runnning Code
for (( val = 0; val < ${#Trees_vals[@]}; val++)); do
	screen -S ${Windname} -X screen -t ${Windname}_${val}
	screen -S ${Windname} -p ${Windname}_${val} -X stuff $"/fs1/Abbie/Get_HaloProps/fullchain_halo_props.sh tree_${Trees_vals[val]} tree_${Trees_end[val]}"
	screen -S ${Windname} -p ${Windname}_${val} -X stuff $'\n'
done

