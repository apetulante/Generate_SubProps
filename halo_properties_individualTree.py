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
num_trees = int(sys.argv[2])

print("For this tree there are %d halos" %num_trees)

subs_today_filename = '%s_subs_today.txt' %tree
    
file = open('%s_halo_properties.txt' %tree, 'w')

for q in range(1,num_trees):
    print('Working on %s halo %d...' % (tree, q))
    tree_filename = '%s_halo_%d.txt' % (tree, q)
    main_filename = '%s_halo_%d_MainHalo.txt' % (tree, q)
    sub_filename = '%s_halo_%d_submerger.txt' % (tree, q)
    
    if os.stat(sub_filename).st_size != 0:
    
        tree_data = pd.read_csv(tree_filename, delim_whitespace=True, header=None)
        halo_num = q
        
        a = tree_data[0].values
        ID = tree_data[1].values
        descID = tree_data[3].values
        num_prog = tree_data[4].values
        pid = tree_data[5].values
        Upid = tree_data[6].values
        desc_pid = tree_data[7].values
        mvir = tree_data[10].values
        rvir = tree_data[11].values
        r_s = tree_data[12].values
        mmp = tree_data[14].values
        Vmax = tree_data[16].values
        x_coord = tree_data[17].values
        y_coord = tree_data[18].values
        z_coord = tree_data[19].values
        v_x = tree_data[20].values
        v_y = tree_data[21].values
        v_z = tree_data[22].values
        r_s_klyp = tree_data[37].values
        spin = tree_data[26].values
        spin_bull = tree_data[45].values
        b_to_a = tree_data[46].values
        c_to_a = tree_data[47].values
        redshift = (1./a) - 1
        conc = rvir/r_s_klyp
            
        main_props = pd.read_csv(main_filename, delim_whitespace=True, header=None)
        
        a_main = main_props[0].values
        ID_main = main_props[1].values
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
    
        sub_merger = np.loadtxt(sub_filename, unpack=True)        
        
        subs_init = sub_merger[0]
        subs_final = sub_merger[1]
        main_init = sub_merger[2]
        main_final = sub_merger[3]
        
        for n in range(0, subs_init.size):
            print('There are subs. Working on %d.....' %n)
            if subs_init.size == 1:
                merge_index = int(subs_init)
                index_sub_initial = int(subs_init)     
                index_sub_final = int(subs_final)
                index_main_initial = int(main_init)
                index_main_final = int(main_final)
            else:
                merge_index = int(subs_init[n])
                index_sub_initial = int(subs_init[n])  
                index_sub_final = int(subs_final[n])
                index_main_initial = int(main_init[n])
                index_main_final = int(main_final[n])
            sub_properties = [[] for x in range(36)]
            sub_properties[0] = tree
            sub_properties[1] = halo_num
            sub_properties[2] = merge_index
            
            if index_main_initial != index_main_final:
                initial_sub_mass = mvir[index_sub_initial]
                initial_host_mass = mvir_main[index_main_initial]
                initial_sub_rvir = rvir[index_sub_initial]
                initial_host_rvir = rvir_main[index_main_initial]
                
                initial_sub_conc = conc[index_sub_initial]
                initial_host_conc = conc_main[index_main_initial]
                
                initial_sub_spin = spin[index_sub_initial]
                initial_host_spin = spin_main[index_main_initial]
                initial_sub_spin_bull = spin_bull[index_sub_initial]
                initial_host_spin_bull = spin_bull_main[index_main_initial]
        
                initial_sub_rel_vx = v_x[index_sub_initial] - v_x_main[index_main_initial] 
                initial_sub_rel_vy = v_y[index_sub_initial] - v_y_main[index_main_initial] 
                initial_sub_rel_vz = v_z[index_sub_initial] - v_z_main[index_main_initial]            
                initial_sub_rel_x = x_coord[index_sub_initial] - x_coord_main[index_main_initial] 
                initial_sub_rel_y = y_coord[index_sub_initial] - y_coord_main[index_main_initial]
                initial_sub_rel_z = z_coord[index_sub_initial] - z_coord_main[index_main_initial] 
                
                initial_scale = a[index_sub_initial] 
                initial_sub_Vmax = Vmax[index_sub_initial]                   
                
                initial_sub_b_to_a = b_to_a[index_sub_initial]
                initial_sub_c_to_a = c_to_a[index_sub_initial]
                initial_host_b_to_a = b_to_a_main[index_main_initial]
                initial_host_c_to_a = c_to_a_main[index_main_initial]
         
                sub_properties[3] = ID[index_sub_final]
                sub_properties[4] = initial_sub_mass
                sub_properties[5] = initial_host_mass
                sub_properties[6] = initial_sub_rvir
                sub_properties[7] = initial_host_rvir
                sub_properties[8] = initial_sub_conc
                sub_properties[9] = initial_host_conc
                sub_properties[10] = initial_sub_spin
                sub_properties[11] = initial_host_spin
                sub_properties[12] = initial_sub_spin_bull
                sub_properties[13] = initial_host_spin_bull
                sub_properties[14] = initial_sub_rel_vx
                sub_properties[15] = initial_sub_rel_vy
                sub_properties[16] = initial_sub_rel_vz
                sub_properties[17] = initial_sub_rel_x
                sub_properties[18] = initial_sub_rel_y
                sub_properties[19] = initial_sub_rel_z
                sub_properties[20] = initial_scale
                sub_properties[21] = initial_sub_Vmax
                sub_properties[22] = initial_sub_b_to_a
                sub_properties[23] = initial_sub_c_to_a
                sub_properties[24] = initial_host_b_to_a
                sub_properties[25] = initial_host_c_to_a
                
                
                survives = 0
                final_scale = a[index_sub_final]
                final_sub_mass = mvir[index_sub_final]
                final_host_mass = mvir_main[index_main_final]
                final_sub_rel_x = x_coord[index_sub_final] - x_coord_main[index_main_final]
                final_sub_rel_y = y_coord[index_sub_final] - y_coord_main[index_main_final]
                final_sub_rel_z = z_coord[index_sub_final] - z_coord_main[index_main_final]
                final_sub_conc = conc[index_sub_final]
                final_sub_spin = spin[index_sub_final]
                final_sub_spin_bull = spin_bull[index_sub_final]
                
                sub_properties[26] = survives
                sub_properties[27] = final_scale
                sub_properties[28] = final_sub_mass   
                sub_properties[29] = final_host_mass
                sub_properties[30] = final_sub_rel_x
                sub_properties[31] = final_sub_rel_y
                sub_properties[32] = final_sub_rel_z
                sub_properties[33] = final_sub_conc
                sub_properties[34] = final_sub_spin
                sub_properties[35] = final_sub_spin_bull      
                
                file.write(sub_properties[0])
                file.write(' ')
                for j in range(1,len(sub_properties)):
                    file.write('%lf' %sub_properties[j])
                    file.write(' ')
                
                file.write('\n')
                
counter = 0
if os.stat(subs_today_filename).st_size != 0:
    subs_today = pd.read_csv(subs_today_filename, delim_whitespace=True, header=None)
    surviving_sub = subs_today[0].values
    surviving_sub_host = subs_today[1].values
    for q in range(0, len(surviving_sub)):
        if surviving_sub_host[q] != 0:
            print('Working on %s file %d...' %(tree, surviving_sub[q]))
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
            redshift_sub = (1./a_sub) - 1
            conc_sub = rvir_sub/r_s_klyp_sub
            
            main_props = pd.read_csv(main_filename, delim_whitespace=True, header=None)
            
            a_main = main_props[0].values
            ID_main = main_props[1].values
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
                        
            #isolate section where the guy is in the thing
            entry_index = 0
            for k in range(len(Upid_sub)):
                if Upid_sub[k]!=ID_main[k]:
                    entry_index = k-1
                    break 
                else:
                    entry_index = k
                    
                if k == len(ID_main)-1:
                    entry_index = k
                    break
                
            
            if entry_index > 0:
                if mvir_sub[entry_index] >= 1000*3.3e7:
                    sub_properties = [[] for x in range(36)]
                    sub_properties[0] = tree
                    sub_properties[1] = halo_num
                    sub_properties[2] = subtree_ID
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
                    initial_sub_rel_x = x_coord_sub[entry_index] - x_coord_main[entry_index] 
                    initial_sub_rel_y = y_coord_sub[entry_index] - y_coord_main[entry_index]
                    initial_sub_rel_z = z_coord_sub[entry_index] - z_coord_main[entry_index] 
                    
                    initial_scale = a_sub[entry_index] 
                    initial_sub_Vmax = Vmax_sub[entry_index]                    
                    
                    initial_sub_b_to_a = b_to_a_sub[entry_index]
                    initial_sub_c_to_a = c_to_a_sub[entry_index]
                    initial_host_b_to_a = b_to_a_main[entry_index]
                    initial_host_c_to_a = c_to_a_main[entry_index]
             
         
                    sub_properties[3] = ID_sub[0]
                    sub_properties[4] = initial_sub_mass
                    sub_properties[5] = initial_host_mass
                    sub_properties[6] = initial_sub_rvir
                    sub_properties[7] = initial_host_rvir
                    sub_properties[8] = initial_sub_conc
                    sub_properties[9] = initial_host_conc
                    sub_properties[10] = initial_sub_spin
                    sub_properties[11] = initial_host_spin
                    sub_properties[12] = initial_sub_spin_bull
                    sub_properties[13] = initial_host_spin_bull
                    sub_properties[14] = initial_sub_rel_vx
                    sub_properties[15] = initial_sub_rel_vy
                    sub_properties[16] = initial_sub_rel_vz
                    sub_properties[17] = initial_sub_rel_x
                    sub_properties[18] = initial_sub_rel_y
                    sub_properties[19] = initial_sub_rel_z
                    sub_properties[20] = initial_scale
                    sub_properties[21] = initial_sub_Vmax
                    sub_properties[22] = initial_sub_b_to_a
                    sub_properties[23] = initial_sub_c_to_a
                    sub_properties[24] = initial_host_b_to_a
                    sub_properties[25] = initial_host_c_to_a
                    
                    survives = 1
                    final_scale = a_sub[0]
                    final_sub_mass = mvir_sub[0]
                    final_host_mass = mvir_main[0]
                    final_sub_rel_x = z_coord_sub[0] - x_coord_main[0]
                    final_sub_rel_y = y_coord_sub[0] - y_coord_main[0]
                    final_sub_rel_z = x_coord_sub[0] - x_coord_main[0]
                    final_sub_conc = conc_sub[0]
                    final_sub_spin = spin_sub[0]
                    final_sub_spin_bull = spin_bull_sub[0]
                    
                    sub_properties[26] = survives
                    sub_properties[27] = final_scale
                    sub_properties[28] = final_sub_mass 
                    sub_properties[29] = final_host_mass 
                    sub_properties[30] = final_sub_rel_x
                    sub_properties[31] = final_sub_rel_y
                    sub_properties[32] = final_sub_rel_z
                    sub_properties[33] = final_sub_conc
                    sub_properties[34] = final_sub_spin
                    sub_properties[35] = final_sub_spin_bull
                    
                    file.write(sub_properties[0])
                    file.write(' ')
                    for j in range(1,len(sub_properties)):
                        file.write('%lf' %sub_properties[j])
                        file.write(' ')
                    
                    counter += 1
                    file.write('\n')
    
file.close()