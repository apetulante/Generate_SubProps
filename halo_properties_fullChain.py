# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:21:36 2017

@author: petulaa
"""

import numpy as np
import pandas as pd
import scipy.integrate as integ
import sys
import os

# constants
omega_m = .250000
omega_lam = .750000
w = 1
h0 = .7
t_H = 10e9*(1./h0)
Grav_const = 4.483e-24 #(kpc^2/Msun*yr^2)

tree = sys.argv[1]

subs_today_filename = '%s_subs_today.txt' %tree
    
file = open('halo_props_%s_survivors.txt' %tree, 'w')
                
counter = 0
if os.stat(subs_today_filename).st_size != 0:
    subs_today = pd.read_csv(subs_today_filename, delim_whitespace=True, header=None)
    surviving_sub = subs_today[0].values
    surviving_sub_host = subs_today[1].values
    for q in range(0, len(surviving_sub)):
        if surviving_sub_host[q] != 0:
            if counter%100 == 0:
                print('Working on %s subhalo today #%d...' %(tree, surviving_sub[q]))
            sub_filename = '%s_halo_%d_MainHalo.txt' %(tree, surviving_sub[q])
            main_filename = '%s_halo_%d_MainHalo.txt' %(tree, surviving_sub_host[q])
            subtree_ID = 0    
            
            sub_data = pd.read_csv(sub_filename, delim_whitespace=True, header=None)
            halo_num = surviving_sub[q]
            
            a_sub = sub_data[0].values
            ID_sub = sub_data[1].values
            Upid_sub = sub_data[6].values
            mvir_sub = sub_data[10].values
            rvir_sub = sub_data[11].values
            Vmax_sub = sub_data[16].values
            x_coord_sub = sub_data[17].values
            y_coord_sub = sub_data[18].values
            z_coord_sub = sub_data[19].values
            v_x_sub = sub_data[20].values
            v_y_sub = sub_data[21].values
            v_z_sub = sub_data[22].values
            r_s_klyp_sub = sub_data[37].values
            spin_sub = sub_data[26].values
            spin_bull_sub = sub_data[45].values
            b_to_a_sub = sub_data[46].values
            c_to_a_sub = sub_data[47].values
            scale_of_last_MM_sub = sub_data[15].values
            redshift_sub = (1./a_sub) - 1
            conc_sub = rvir_sub/r_s_klyp_sub
            num_subs = sub_data[59].values
            max_sub_mass = sub_data[60].values
            
            main_props = pd.read_csv(main_filename, delim_whitespace=True, header=None)
            
            a_main = main_props[0].values
            ID_main = main_props[1].values
            scale_of_last_MM_main = main_props[15].values
            x_coord_main = main_props[17].values
            y_coord_main = main_props[18].values
            z_coord_main = main_props[19].values
            mvir_main = main_props[10].values
            rvir_main = main_props[11].values
            v_x_main = main_props[20].values
            v_y_main = main_props[21].values
            v_z_main = main_props[22].values
            r_s_klyp_main = main_props[37].values
            spin_main = main_props[26].values
            spin_bull_main = main_props[45].values
            b_to_a_main = main_props[46].values    
            c_to_a_main = main_props[47].values
            conc_main = np.array(rvir_main)/np.array(r_s_klyp_main)
            num_subs_main = main_props[59].values    
            max_sub_mass_main = main_props[60].values
                        
            #isolate section where the guy is in the thing
            entry_index = 0
            for k in range(len(Upid_sub)):
                if Upid_sub[k]!=ID_main[k]:
                    entry_index = k
                    break 
                if k == len(ID_main)-1:
                    entry_index = k
                    break
                
            
            if entry_index > 0:
                if mvir_sub[entry_index] >= 1000*3.3e7:
                    sub_properties = [[] for x in range(49)]
                    sub_properties[0] = tree
                    sub_properties[1] = halo_num
                    sub_properties[2] = entry_index
                    sub_properties[3] = subtree_ID
                    initial_sub_mass = mvir_sub[entry_index]
                    initial_host_mass = mvir_main[entry_index]
                    initial_sub_rvir = rvir_sub[entry_index]
                    initial_host_rvir = rvir_main[entry_index]
                    
                    initial_sub_conc = conc_sub[entry_index]
                    initial_host_conc = conc_main[entry_index]
                    
                    initial_sub_spin = spin_sub[entry_index]
                    initial_host_spin = spin_main[entry_index]
                    initial_sub_spin_bull = spin_bull_sub[entry_index]
                    initial_host_spin_bull = spin_bull_main[entry_index]
            
                    initial_sub_rel_vx = v_x_sub[entry_index] - v_x_main[entry_index] 
                    initial_sub_rel_vy = v_y_sub[entry_index] - v_y_main[entry_index] 
                    initial_sub_rel_vz = v_z_sub[entry_index] - v_z_main[entry_index]            
                    initial_sub_x = x_coord_sub[entry_index]
                    initial_host_x = x_coord_main[entry_index] 
                    initial_sub_y = y_coord_sub[entry_index]
                    initial_host_y = y_coord_main[entry_index]
                    initial_sub_z = z_coord_sub[entry_index]
                    initial_host_z = z_coord_main[entry_index] 
                    
                    initial_scale = a_sub[entry_index] 
                    initial_sub_Vmax = Vmax_sub[entry_index]                    
                    
                    initial_sub_b_to_a = b_to_a_sub[entry_index]
                    initial_sub_c_to_a = c_to_a_sub[entry_index]
                    initial_host_b_to_a = b_to_a_main[entry_index]
                    initial_host_c_to_a = c_to_a_main[entry_index]
                    
                    initial_sub_scale_last_MM = scale_of_last_MM_sub[entry_index]
                    initial_host_scale_last_MM = scale_of_last_MM_main[entry_index]
                    
                    initial_num_subs = num_subs[entry_index]
                    initial_max_sub_mass = max_sub_mass[entry_index]
                    initial_num_host_subs = num_subs_main[entry_index]
                    initial_max_host_sub_mass = max_sub_mass_main[entry_index]
             
         
                    sub_properties[4] = ID_sub[0]
                    sub_properties[5] = initial_sub_mass
                    sub_properties[6] = initial_host_mass
                    sub_properties[7] = initial_sub_rvir
                    sub_properties[8] = initial_host_rvir
                    sub_properties[9] = initial_sub_conc
                    sub_properties[10] = initial_host_conc
                    sub_properties[11] = initial_sub_spin
                    sub_properties[12] = initial_host_spin
                    sub_properties[13] = initial_sub_spin_bull
                    sub_properties[14] = initial_host_spin_bull
                    sub_properties[15] = initial_sub_rel_vx
                    sub_properties[16] = initial_sub_rel_vy
                    sub_properties[17] = initial_sub_rel_vz
                    sub_properties[18] = initial_sub_x
                    sub_properties[19] = initial_host_x
                    sub_properties[20] = initial_sub_y
                    sub_properties[21] = initial_host_y
                    sub_properties[22] = initial_sub_z
                    sub_properties[23] = initial_host_z
                    sub_properties[24] = initial_scale
                    sub_properties[25] = initial_sub_Vmax
                    sub_properties[26] = initial_sub_b_to_a
                    sub_properties[27] = initial_sub_c_to_a
                    sub_properties[28] = initial_host_b_to_a
                    sub_properties[29] = initial_host_c_to_a
                    sub_properties[30] = initial_sub_scale_last_MM
                    sub_properties[31] = initial_host_scale_last_MM
                    sub_properties[32] = initial_num_subs
                    sub_properties[33] = initial_max_sub_mass
                    sub_properties[34] = initial_num_host_subs
                    sub_properties[35] = initial_max_host_sub_mass
                    
                    survives = 1
                    final_scale = a_sub[0]
                    final_sub_mass = mvir_sub[0]
                    final_host_mass = mvir_main[0]
                    final_host_rvir = rvir_main[0]
                    final_sub_rel_x = x_coord_sub[0] - x_coord_main[0]
                    final_sub_rel_y = y_coord_sub[0] - y_coord_main[0]
                    final_sub_rel_z = z_coord_sub[0] - z_coord_main[0]
                    final_sub_conc = conc_sub[0]
                    final_sub_spin = spin_sub[0]
                    final_sub_spin_bull = spin_bull_sub[0]
                    final_sub_scale_last_MM = scale_of_last_MM_sub[0]
                    final_host_scale_last_MM = scale_of_last_MM_main[0]
                    
                    sub_properties[36] = survives
                    sub_properties[37] = final_scale
                    sub_properties[38] = final_sub_mass 
                    sub_properties[39] = final_host_mass 
                    sub_properties[40] = final_host_rvir 
                    sub_properties[41] = final_sub_rel_x
                    sub_properties[42] = final_sub_rel_y
                    sub_properties[43] = final_sub_rel_z
                    sub_properties[44] = final_sub_conc
                    sub_properties[45] = final_sub_spin
                    sub_properties[46] = final_sub_spin_bull
                    sub_properties[47] = final_sub_scale_last_MM
                    sub_properties[48] = final_host_scale_last_MM
                    
                    file.write(sub_properties[0])
                    file.write(' ')
                    for j in range(1,len(sub_properties)):
                        file.write('%lf' %sub_properties[j])
                        file.write(' ')
                    
                    counter += 1
                    file.write('\n')
    
file.close()