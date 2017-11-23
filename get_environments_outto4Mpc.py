#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 11:32:02 2017

@author: petulaa
"""

import pandas as pd
import numpy as np
import periodic_kdtree as pdtree
import sys
from multiprocessing import Pool


filename = sys.argv[1]

halo_props = pd.read_csv(filename, header=None, delim_whitespace=True,
                         names = ['tree', 'halo_num', 'entry_index', 'merge_index', 'initial_ID',
                         'initial_sub_mass', 'initial_host_mass', 'initial_sub_rvir', 
                         'initial_host_rvir', 'initial_sub_conc', 'initial_host_conc' ,
                         'initial_sub_spin','initial_host_spin','initial_sub_spin_bull',
                         'initial_host_spin_bull', 'initial_sub_rel_vx','initial_sub_rel_vy',
                         'initial_sub_rel_vz','initial_sub_x','initial_host_x','initial_sub_y',
                         'initial_host_y','initial_sub_z','initial_host_z','initial_scale',
                         'initial_sub_Vmax','initial_sub_b_to_a','initial_sub_c_to_a',
                         'initial_host_b_to_a','initial_host_c_to_a','initial_sub_scale_last_MM',
                         'initial_host_scale_last_MM','initial_num_subs','initial_max_sub_mass',
                         'initial_num_host_subs','initial_max_host_sub_mass','survives',
                         'final_scale','final_sub_mass','final_host_mass','final_host_rvir',
                         'final_sub_rel_x','final_sub_rel_y','final_sub_rel_z','final_sub_conc',
                         'final_sub_spin','final_sub_spin_bull','final_sub_scale_last_MM',
                         'final_host_scale_last_MM'])

scale_list = halo_props['initial_scale'].unique()
box_bounds = np.array([130.0,130.0,130.0])

for scale in scale_list:
    print("Checking at scale = %lf" %scale)
    if len(np.where(halo_props['initial_scale'].values==scale)[0]) > 0:
        print("Get the hlist file")
        hlist_filepath = '/fs2/shared/new_130_mpc_box/hires/5045/rockstar_halos/so_m200b/full_res/hlists/hlist_%.5lf.list.gz' %scale
        hlist_all = pd.read_csv(hlist_filepath, comment = '#', index_col = False, 
                                compression='gzip', header = 0, usecols = [1,5,10,17,18,19], engine = 'c',
                                delim_whitespace=True, names = ['ID','pid','mvir', 'x', 'y', 'z'])
        #print(hlist_all.head())
        # PERIODIC BOUNDARY CONDITIONS
        hlist_x = hlist_all['x'].values
        hlist_y = hlist_all['y'].values
        hlist_z = hlist_all['z'].values
        hlist_positions = np.transpose([hlist_x, hlist_y, hlist_z])
        hlist_KDtree = pdtree.PeriodicCKDTree(box_bounds, hlist_positions)
        for i in halo_props.index[halo_props['initial_scale'] == scale].tolist():
            print("Searching.....")
            sub_x = halo_props.iloc[i].initial_sub_x
            sub_y = halo_props.iloc[i].initial_sub_y
            sub_z = halo_props.iloc[i].initial_sub_z
            sub_position = [sub_x,sub_y,sub_z]
            indicies_4Mpc = hlist_KDtree.query_ball_point(sub_position,r=4)
            
            print("Recording.....")
            print("We have %d neighbors to get through." %len(indicies_4Mpc))
            density_8Mpc = 0
            tidal_8Mpc = 0
            vol_8Mpc = (4/3.0)*np.pi*(8**3)
            density_4Mpc = 0
            tidal_4Mpc = 0
            vol_4Mpc = (4/3.0)*np.pi*(4**3)
            density_2Mpc = 0
            tidal_2Mpc = 0
            vol_2Mpc = (4/3.0)*np.pi*(2**3)
            density_1Mpc = 0
            tidal_1Mpc = 0
            vol_1Mpc = (4/3.0)*np.pi*(1**3)
            for neighbor in indicies_4Mpc:
           #     print('Looking at neighbor %d',%neighbor)
                if hlist_all.iloc[neighbor].pid == -1.0:
                    rel_x = hlist_all.iloc[neighbor].x-sub_x
                    rel_y = hlist_all.iloc[neighbor].y-sub_y
                    rel_z = hlist_all.iloc[neighbor].z-sub_z
                    if rel_x >= 130.0/2:
                        rel_x = 130.0 - rel_x
                    if rel_y >= 130.0/2:
                        rel_y = 130.0 - rel_y
                    if rel_z >= 130.0/2:
                        rel_z = 130.0 - rel_z
                    neighbor_dist = np.sqrt(rel_x**2 + rel_y**2 + rel_z**2)
                  #  print("the neighbor distance is %lf" %neighbor_dist)
                    if neighbor_dist == 0:
                        continue
                    if hlist_all.iloc[neighbor].ID == halo_props.iloc[i].initial_ID:
                        continue
                    neighbor_vector = (rel_x/neighbor_dist + rel_y/neighbor_dist + rel_z/neighbor_dist)
                    if neighbor_vector == 0:
                        continue
                    #density_8Mpc += hlist_all.iloc[neighbor].mvir/vol_8Mpc
                    #tidal_8Mpc += (hlist_all.iloc[neighbor].mvir/(neighbor_dist**3))*neighbor_vector
                    if (neighbor_dist < 4):
                        density_4Mpc += hlist_all.iloc[neighbor].mvir/vol_4Mpc
                        tidal_4Mpc += (hlist_all.iloc[neighbor].mvir/(neighbor_dist**3))/neighbor_vector
                        if (neighbor_dist < 2):
                            density_2Mpc += hlist_all.iloc[neighbor].mvir/vol_2Mpc
                            tidal_2Mpc += (hlist_all.iloc[neighbor].mvir/(neighbor_dist**3))/neighbor_vector
                            if (neighbor_dist < 1):
                                density_1Mpc += hlist_all.iloc[neighbor].mvir/vol_1Mpc
                                tidal_1Mpc += (hlist_all.iloc[neighbor].mvir/(neighbor_dist**3))/neighbor_vector
                
            halo_props.set_value(i,'density_1Mpc',density_1Mpc)
            halo_props.set_value(i,'density_2Mpc',density_2Mpc)
            halo_props.set_value(i,'density_4Mpc',density_4Mpc)
           # halo_props.set_value(i,'density_8Mpc',density_8Mpc)
            
            halo_props.set_value(i,'tidal_1Mpc',tidal_1Mpc)
            halo_props.set_value(i,'tidal_2Mpc',tidal_2Mpc)
            halo_props.set_value(i,'tidal_4Mpc',tidal_4Mpc)
           # halo_props.set_value(i,'tidal_8Mpc',tidal_8Mpc)
    
            
out_filename = "COMPLETE_%s" %filename
output_file = open(out_filename,'w')
halo_props.to_csv(output_file, sep=' ', index=False, header=False)
output_file.close()
                
            
                
            
    
    
    