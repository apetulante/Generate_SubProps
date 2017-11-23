#!/bin/sh

Windname="Environments"

Trees_vals=("0.63880" "0.288500" "0.370270" "0.440330" "0.50320" "0.558640" "0.612730" "0.662810" "0.708740" "0.754360" "0.795540" "0.833160" "0.872560" "0.913830" "0.957050")
Trees_end=("0.288500" "0.370270" "0.440330" "0.50320" "0.558640" "0.612730" "0.662810" "0.708740" "0.754360" "0.795540" "0.833160" "0.872560" "0.913830" "0.957050" "1.007690")

# Initialize screen session
screen -mdS ${Windname}

## Runnning Code
for (( val = 0; val < ${#Trees_vals[@]}; val++)); do
	screen -S ${Windname} -X screen -t ${Windname}_${val}
	screen -S ${Windname} -p ${Windname}_${val} -X stuff $"python get_environments_outto4Mpc.py halo_properties_scale${Trees_vals[val]}_to_scale${Trees_end[val]}.txt"
	screen -S ${Windname} -p ${Windname}_${val} -X stuff $'\n'
done

