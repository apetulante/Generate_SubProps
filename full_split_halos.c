#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void write_submerger(FILE *out_file, int *sub_merger_start, int *sub_merger_end, int *main_halo_start, int *main_halo_end, int num_subs) {

     int l;
 //    printf("Writing the subs merger file...\n");
     for (l=0; l<num_subs; l++) {
        fprintf(out_file, "%d ", sub_merger_start[l]);
        fprintf(out_file, "%d ", sub_merger_end[l]);
        fprintf(out_file, "%d ", main_halo_start[l]);
        fprintf(out_file, "%d ", main_halo_end[l]);
        fprintf(out_file, "\n"); 
     }
}

void write_MainHalo(FILE *out_file2, double *a_main, double *ID_main, double *desc_scale_main, double *descID_main,
  double *num_prog_main, double *pid_main, double *Upid_main, double *desc_pid_main, double *phantom_main, double *SAM_mvir_main, 
  double *mvir_main, double *rvir_main, double *r_s_main, double *Vrms_main, double *mmp_main, double *last_MM_scale_main, 
  double *Vmax_main, double *x_coord_main, double *y_coord_main, double *z_coord_main, double *x_vel_main, double *y_vel_main,
  double *z_vel_main, double *x_ang_main, double *y_ang_main, double *z_ang_main, double *spin_main, double *breadth_first_ID_main,
  double *depth_first_ID_main, double *tree_root_ID_main, double *orig_halo_ID_main, double *snap_num_main, double *next_coprog_dep_ID_main,
  double *last_mainleaf_dep_ID_main, double *tidal_F_main, double *tidal_ID_main, double *r_s_klyp_main, double *M200b_all_main,
  double *M200b_1_main, double *M200b_2_main, double *M200b_3_main, double *M200b_4_main, double *Xoff_main, double *Voff_main, 
  double *spin_bull_main, double *b_to_a_main, double *c_to_a_main, double *A_x_main, double *A_y_main, double *A_z_main, 
  double *T_U_rat_main, double *M_pe_main, double *halfmass_r_main, double *misc_54_main, double *misc_55_main, double *misc_56_main,
  double *misc_57_main, double *misc_58_main, double *misc_59_main, double *misc_60_main, int len_main) {

     int k;
  //   printf("Writing the main halo file...\n");
     //printf("len_main inside the function is %d\n", len_main);
     for (k=0; k<len_main; k++) {
      //  printf("the loop started\n");
      //  printf("k is %d\n", k);
        fprintf(out_file2, "%lf ", a_main[k]); 
      //  printf("We're printing the main scale and it is: %lf\n", a_main[k]);
        fprintf(out_file2, "%lf ", ID_main[k]);
        fprintf(out_file2, "%lf ", desc_scale_main[k]);
        fprintf(out_file2, "%lf ", descID_main[k]);
        fprintf(out_file2, "%lf ", num_prog_main[k]);
        fprintf(out_file2, "%lf ", pid_main[k]);
        fprintf(out_file2, "%lf ", Upid_main[k]);
        fprintf(out_file2, "%lf ", desc_pid_main[k]);
        fprintf(out_file2, "%lf ", phantom_main[k]);
        fprintf(out_file2, "%lf ", SAM_mvir_main[k]);
        fprintf(out_file2, "%lf ", mvir_main[k]);
        fprintf(out_file2, "%lf ", rvir_main[k]);
        fprintf(out_file2, "%lf ", r_s_main[k]);
        fprintf(out_file2, "%lf ", Vrms_main[k]);
        fprintf(out_file2, "%lf ", mmp_main[k]);
        fprintf(out_file2, "%lf ", last_MM_scale_main[k]);
        fprintf(out_file2, "%lf ", Vmax_main[k]);
        fprintf(out_file2, "%lf ", x_coord_main[k]);
        fprintf(out_file2, "%lf ", y_coord_main[k]);
        fprintf(out_file2, "%lf ", z_coord_main[k]);
        fprintf(out_file2, "%lf ", x_vel_main[k]);
        fprintf(out_file2, "%lf ", y_vel_main[k]);
        fprintf(out_file2, "%lf ", z_vel_main[k]);
        fprintf(out_file2, "%lf ", x_ang_main[k]);
        fprintf(out_file2, "%lf ", y_ang_main[k]);
        fprintf(out_file2, "%lf ", z_ang_main[k]);
        fprintf(out_file2, "%lf ", spin_main[k]);
        fprintf(out_file2, "%lf ", breadth_first_ID_main[k]);
        fprintf(out_file2, "%lf ", depth_first_ID_main[k]);
        fprintf(out_file2, "%lf ", tree_root_ID_main[k]);
        fprintf(out_file2, "%lf ", orig_halo_ID_main[k]);
        fprintf(out_file2, "%lf ", snap_num_main[k]);
        fprintf(out_file2, "%lf ", next_coprog_dep_ID_main[k]);
        fprintf(out_file2, "%lf ", last_mainleaf_dep_ID_main[k]);
        fprintf(out_file2, "%lf ", tidal_F_main[k]);
        fprintf(out_file2, "%lf ", tidal_ID_main[k]);
        fprintf(out_file2, "%lf ", r_s_klyp_main[k]);
        fprintf(out_file2, "%lf ", M200b_all_main[k]);
        fprintf(out_file2, "%lf ", M200b_1_main[k]);
        fprintf(out_file2, "%lf ", M200b_2_main[k]);
        fprintf(out_file2, "%lf ", M200b_3_main[k]);
        fprintf(out_file2, "%lf ", M200b_4_main[k]);
        fprintf(out_file2, "%lf ", Xoff_main[k]);
        fprintf(out_file2, "%lf ", Voff_main[k]);
        fprintf(out_file2, "%lf ", spin_bull_main[k]);
        fprintf(out_file2, "%lf ", b_to_a_main[k]);
        fprintf(out_file2, "%lf ", c_to_a_main[k]);
        fprintf(out_file2, "%lf ", A_x_main[k]);
        fprintf(out_file2, "%lf ", A_y_main[k]);
        fprintf(out_file2, "%lf ", A_z_main[k]);
        fprintf(out_file2, "%lf ", T_U_rat_main[k]);
        fprintf(out_file2, "%lf ", M_pe_main[k]);
        fprintf(out_file2, "%lf ", misc_54_main[k]);
        fprintf(out_file2, "%lf ", misc_55_main[k]);
        fprintf(out_file2, "%lf ", misc_56_main[k]);
        fprintf(out_file2, "%lf ", misc_57_main[k]);
        fprintf(out_file2, "%lf ", misc_58_main[k]);
        fprintf(out_file2, "%lf ", misc_59_main[k]);
        fprintf(out_file2, "%lf ", misc_60_main[k]);
        fprintf(out_file2, "\n");
      }
  }

void calc_subs(int *len_main, int *num_subs, double *a, double *ID, double *desc_scale, double *descID, 
  double *num_prog, double *pid, double *Upid, double *desc_pid, double *phantom, double *SAM_mvir, double *mvir, 
  double *rvir, double *r_s, double *Vrms, double *mmp, double *last_MM_scale, double *Vmax, double *x_coord, double *y_coord, 
  double *z_coord, double *x_vel, double *y_vel, double *z_vel, double *x_ang, double *y_ang, double *z_ang, double *spin, 
  double *breadth_first_ID, double *depth_first_ID, double *tree_root_ID, double *orig_halo_ID, double *snap_num, 
  double *next_coprog_dep_ID, double *last_mainleaf_dep_ID, double *tidal_F, double *tidal_ID, double *r_s_klyp, 
  double *M200b_all, double *M200b_1, double *M200b_2, double *M200b_3, double *M200b_4, double *Xoff, double *Voff, 
  double *spin_bull, double *b_to_a, double *c_to_a, double *A_x, double *A_y, double *A_z, double *T_U_rat, double *M_pe, 
  double *halfmass_r, double *misc_54, double *misc_55, double *misc_56, double *misc_57, double *misc_58, double *misc_59, 
  double *misc_60, int len_tree, double *a_main, double *ID_main, double *desc_scale_main, double *descID_main,
  double *num_prog_main, double *pid_main, double *Upid_main, double *desc_pid_main, double *phantom_main, double *SAM_mvir_main, 
  double *mvir_main, double *rvir_main, double *r_s_main, double *Vrms_main, double *mmp_main, double *last_MM_scale_main, 
  double *Vmax_main, double *x_coord_main, double *y_coord_main, double *z_coord_main, double *x_vel_main, double *y_vel_main,
  double *z_vel_main, double *x_ang_main, double *y_ang_main, double *z_ang_main, double *spin_main, double *breadth_first_ID_main,
  double *depth_first_ID_main, double *tree_root_ID_main, double *orig_halo_ID_main, double *snap_num_main, double *next_coprog_dep_ID_main,
  double *last_mainleaf_dep_ID_main, double *tidal_F_main, double *tidal_ID_main, double *r_s_klyp_main, double *M200b_all_main,
  double *M200b_1_main, double *M200b_2_main, double *M200b_3_main, double *M200b_4_main, double *Xoff_main, double *Voff_main, 
  double *spin_bull_main, double *b_to_a_main, double *c_to_a_main, double *A_x_main, double *A_y_main, double *A_z_main, 
  double *T_U_rat_main, double *M_pe_main, double *halfmass_r_main, double *misc_54_main, double *misc_55_main, double *misc_56_main,
  double *misc_57_main, double *misc_58_main, double *misc_59_main, double *misc_60_main, int *sub_merger_start, int *sub_merger_end,
  int *main_halo_start, int *main_halo_end) {
    
    a_main[0] = a[0];
    ID_main[0] = ID[0];
    desc_scale_main[0] = desc_scale[0];
    descID_main[0] = descID[0];
    num_prog_main[0] = num_prog[0];
    pid_main[0] = pid[0];
    Upid_main[0] = Upid[0];
    desc_pid_main[0] = desc_pid[0];
    phantom_main[0] = phantom[0];
    SAM_mvir_main[0] = SAM_mvir[0];
    mvir_main[0] = mvir[0];
    rvir_main[0] = rvir[0];
    r_s_main[0] = r_s[0];
    Vrms_main[0] = Vrms[0];
    mmp_main[0] = mmp[0];
    last_MM_scale_main[0] = last_MM_scale[0];
    Vmax_main[0] = Vmax[0];
    x_coord_main[0] = x_coord[0];
    y_coord_main[0] = y_coord[0];
    z_coord_main[0] = z_coord[0];
    x_vel_main[0] = x_vel[0];
    y_vel_main[0] = y_vel[0];
    z_vel_main[0] = z_vel[0];
    x_ang_main[0] = x_ang[0];
    y_ang_main[0] = y_ang[0];
    z_ang_main[0] = z_ang[0];
    spin_main[0] = spin[0];
    breadth_first_ID_main[0] = breadth_first_ID[0];
    depth_first_ID_main[0] = depth_first_ID[0];
    tree_root_ID_main[0] = tree_root_ID[0];
    orig_halo_ID_main[0] = orig_halo_ID[0];
    snap_num_main[0] = snap_num[0];
    next_coprog_dep_ID_main[0] = next_coprog_dep_ID[0];
    last_mainleaf_dep_ID_main[0] = last_mainleaf_dep_ID[0];
    tidal_F_main[0] = tidal_F[0];
    tidal_ID_main[0] = tidal_ID[0];
    r_s_klyp_main[0] = r_s_klyp[0];
    M200b_all_main[0] = M200b_all[0];
    M200b_1_main[0] = M200b_1[0];
    M200b_2_main[0] = M200b_2[0];
    M200b_3_main[0] = M200b_3[0];
    M200b_4_main[0] = M200b_4[0];
    Xoff_main[0] = Xoff[0];
    Voff_main[0] = Voff[0];
    spin_bull_main[0] = spin_bull[0];
    b_to_a_main[0] = b_to_a[0];
    c_to_a_main[0] = c_to_a[0];
    A_x_main[0] = A_x[0];
    A_y_main[0] = A_y[0];
    A_z_main[0] = A_z[0];
    T_U_rat_main[0] = T_U_rat[0];
    M_pe_main[0] = M_pe[0];
    halfmass_r_main[0] = halfmass_r[0];
    misc_54_main[0] = misc_54[0];
    misc_55_main[0] = misc_55[0];
    misc_56_main[0] = misc_56[0];
    misc_57_main[0] = misc_57[0];
    misc_58_main[0] = misc_58[0];
    misc_59_main[0] = misc_59[0];
    misc_60_main[0] = misc_60[0];


    int i;
    int j = 1;
    for (i=1; i<len_tree; i++) {
        //printf("i is: %d\n", i);
        //printf("descID is %lf and ID_main is %lf\n", descID[i], ID_main[j-1]);
        if (descID[i] == ID_main[j-1]) {
//            if (pid[i] == -1.0) {
//                printf("Is it the mmp? %lf\n", mmp[i]);
                if (mmp[i] == 1.0) {
//                    printf("Writing main halo arrays for scale: %lf\n", a[j]);
                    a_main[j] = a[i];
//                    printf("just wrote a_main, it is now %lf\n", a_main[j]);
                    ID_main[j] = ID[i];
                    desc_scale_main[j] = desc_scale[i];
                    descID_main[j] = descID[i];
                    num_prog_main[j] = num_prog[i];
                    pid_main[j] = pid[i];
                    Upid_main[j] = Upid[i];
                    desc_pid_main[j] = desc_pid[i];
                    phantom_main[j] = phantom[i];
                    SAM_mvir_main[j] = SAM_mvir[i];
                    mvir_main[j] = mvir[i];
                    rvir_main[j] = rvir[i];
                    r_s_main[j] = r_s[i];
                    Vrms_main[j] = Vrms[i];
                    mmp_main[j] = mmp[i];
                    last_MM_scale_main[j] = last_MM_scale[i];
                    Vmax_main[j] = Vmax[i];
                    x_coord_main[j] = x_coord[i];
                    y_coord_main[j] = y_coord[i];
                    z_coord_main[j] = z_coord[i];
                    x_vel_main[j] = x_vel[i];
                    y_vel_main[j] = y_vel[i];
                    z_vel_main[j] = z_vel[i];
                    x_ang_main[j] = x_ang[i];
                    y_ang_main[j] = y_ang[i];
                    z_ang_main[j] = z_ang[i];
                    spin_main[j] = spin[i];
                    breadth_first_ID_main[j] = breadth_first_ID[i];
                    depth_first_ID_main[j] = depth_first_ID[i];
                    tree_root_ID_main[j] = tree_root_ID[i];
                    orig_halo_ID_main[j] = orig_halo_ID[i];
                    snap_num_main[j] = snap_num[i];
                    next_coprog_dep_ID_main[j] = next_coprog_dep_ID[i];
                    last_mainleaf_dep_ID_main[j] = last_mainleaf_dep_ID[i];
                    tidal_F_main[j] = tidal_F[i];
                    tidal_ID_main[j] = tidal_ID[i];
                    r_s_klyp_main[j] = r_s_klyp[i];
                    M200b_all_main[j] = M200b_all[i];
                    M200b_1_main[j] = M200b_1[i];
                    M200b_2_main[j] = M200b_2[i];
                    M200b_3_main[j] = M200b_3[i];
                    M200b_4_main[j] = M200b_4[i];
                    Xoff_main[j] = Xoff[i];
                    Voff_main[j] = Voff[i];
                    spin_bull_main[j] = spin_bull[i];
                    b_to_a_main[j] = b_to_a[i];
                    c_to_a_main[j] = c_to_a[i];
                    A_x_main[j] = A_x[i];
                    A_y_main[j] = A_y[i];
                    A_z_main[j] = A_z[i];
                    T_U_rat_main[j] = T_U_rat[i];
                    M_pe_main[j] = M_pe[i];
                    halfmass_r_main[j] = halfmass_r[i];
                    misc_54_main[j] = misc_54[i];
                    misc_55_main[j] = misc_55[i];
                    misc_56_main[j] = misc_56[i];
                    misc_57_main[j] = misc_57[i];
                    misc_58_main[j] = misc_58[i];
                    misc_59_main[j] = misc_59[i];
                    misc_60_main[j] = misc_60[i];
                    j ++;
                 }
             //}
         }
     }

     int l = 0;
     double sub_endID;
     double entry_mass = 0;
     int main_ind;
     int n,m,p;
     int entry_index = 0;
     int entry_host_index = 0;
     for (n=0; n < j; n++) {
         for (m=0; m < len_tree; m++) {
            if (descID[m] == ID_main[n]) {
                if (Upid[m] == ID_main[n+1]) {
             //       printf("Found a sub! At scale: %lf\n", a[m]);
        sub_endID = ID[m];
                    main_ind = n+2;
                    for(p=m+1; p<len_tree; p++) {
                        if (descID[p] == sub_endID) {
              sub_endID = ID[p];
              if (Upid[p] == ID_main[main_ind]) {
              entry_mass = mvir[p];
          entry_index = p;
              entry_host_index = main_ind;
        main_ind ++;
          //    printf("It's still living inside\n");
          }
          else 
                                break;
                        }
                    }
                    if (entry_mass > 1000*3.3e7) {
                        sub_merger_end[l] = m;
                        sub_merger_start[l] = entry_index;
                        main_halo_end[l] = n+1;
                        main_halo_start[l] = entry_host_index;
                        l ++;
                    }
                }
            }
         }
     }
     *len_main = j;
     *num_subs = l;
}

void find_main_subs(double *a, double *ID, double *desc_scale, double *descID, 
  double *num_prog, double *pid, double *Upid, double *desc_pid, double *phantom, double *SAM_mvir, 
  double *mvir, double *rvir, double *r_s, double *Vrms, double *mmp, double *last_MM_scale, 
  double *Vmax, double *x_coord, double *y_coord, double *z_coord, double *x_vel, double *y_vel, 
  double *z_vel, double *x_ang, double *y_ang, double *z_ang, double *spin, double *breadth_first_ID, 
  double *depth_first_ID, double *tree_root_ID, double *orig_halo_ID, double *snap_num, 
  double *next_coprog_dep_ID, double *last_mainleaf_dep_ID, double *tidal_F, double *tidal_ID, 
  double *r_s_klyp, double *M200b_all, double *M200b_1, double *M200b_2, double *M200b_3, 
  double *M200b_4, double *Xoff, double *Voff, double *spin_bull, double *b_to_a, double *c_to_a, 
  double *A_x, double *A_y, double *A_z, double *T_U_rat, double *M_pe, double *halfmass_r, 
  double *misc_54, double *misc_55, double *misc_56, double *misc_57, double *misc_58, double *misc_59, 
  double *misc_60, int tree_size, int *num_subs, int tree_num, int *subs_today, 
  int *hosts_today, double *newtree_main_ids, int file_ident_1, int file_ident_2, int file_ident_3) {
     
  int halo_num = tree_num;   
  int len_tree = tree_size;

     int b;
    //  int current_num_subs = *num_subs;
    // printf("The current number of subs is: %d\n", current_num_subs);

     newtree_main_ids[tree_num] = ID[0];
    // printf("We just added a new tree main ID, it's %lf", newtree_main_ids[tree_num]);

     if (Upid[0] != -1.0) {
    //     printf("We found a sub in %d\n", tree_num);
         subs_today[*num_subs] = tree_num;
         for (b=0; b<tree_num; b++) {
             if (newtree_main_ids[b] == Upid[0]) {
                 hosts_today[*num_subs] = b;
 //                printf("This tree is for a sub! Tree number %d is a sub of tree %d\n", subs_today[*num_subs], hosts_today[*num_subs]);
             }
         }
         *num_subs = *num_subs + 1;
  //       printf("The number of subs is: %d\n", *num_subs);
     }

        //create arrays to calculate sub start indicies
   double *a_main = malloc(1000*sizeof(double));
   double *ID_main = malloc(1000*sizeof(double));
   double *desc_scale_main = malloc(1000*sizeof(double));
   double *descID_main = malloc(1000*sizeof(double));
   double *num_prog_main = malloc(1000*sizeof(double));
   double *pid_main = malloc(1000*sizeof(double));
   double *Upid_main = malloc(1000*sizeof(double));
   double *desc_pid_main = malloc(1000*sizeof(double));
   double *phantom_main = malloc(1000*sizeof(double));
   double *SAM_mvir_main = malloc(1000*sizeof(double));
   double *mvir_main = malloc(1000*sizeof(double));
   double *rvir_main = malloc(1000*sizeof(double));
   double *r_s_main = malloc(1000*sizeof(double));
   double *Vrms_main = malloc(1000*sizeof(double));
   double *mmp_main = malloc(1000*sizeof(double));
   double *last_MM_scale_main = malloc(1000*sizeof(double));
   double *Vmax_main = malloc(1000*sizeof(double));
   double *x_coord_main = malloc(1000*sizeof(double));
   double *y_coord_main = malloc(1000*sizeof(double));
   double *z_coord_main = malloc(1000*sizeof(double)); 
   double *x_vel_main = malloc(1000*sizeof(double));
   double *y_vel_main = malloc(1000*sizeof(double));
   double *z_vel_main = malloc(1000*sizeof(double));
   double *x_ang_main = malloc(1000*sizeof(double));
   double *y_ang_main = malloc(1000*sizeof(double));
   double *z_ang_main = malloc(1000*sizeof(double));
   double *spin_main = malloc(1000*sizeof(double));
   double *breadth_first_ID_main = malloc(1000*sizeof(double));
   double *depth_first_ID_main = malloc(1000*sizeof(double));
   double *tree_root_ID_main = malloc(1000*sizeof(double));
   double *orig_halo_ID_main = malloc(1000*sizeof(double));
   double *snap_num_main = malloc(1000*sizeof(double));
   double *next_coprog_dep_ID_main = malloc(1000*sizeof(double));
   double *last_mainleaf_dep_ID_main = malloc(1000*sizeof(double));
   double *tidal_F_main = malloc(1000*sizeof(double));
   double *tidal_ID_main = malloc(1000*sizeof(double));
   double *r_s_klyp_main = malloc(1000*sizeof(double));
   double *M200b_all_main = malloc(1000*sizeof(double));
   double *M200b_1_main = malloc(1000*sizeof(double));
   double *M200b_2_main = malloc(1000*sizeof(double));
   double *M200b_3_main = malloc(1000*sizeof(double));
   double *M200b_4_main = malloc(1000*sizeof(double));
   double *Xoff_main = malloc(1000*sizeof(double));
   double *Voff_main = malloc(1000*sizeof(double));
   double *spin_bull_main = malloc(1000*sizeof(double));
   double *b_to_a_main = malloc(1000*sizeof(double));
   double *c_to_a_main = malloc(1000*sizeof(double));
   double *A_x_main = malloc(1000*sizeof(double));
   double *A_y_main = malloc(1000*sizeof(double));
   double *A_z_main = malloc(1000*sizeof(double));
   double *T_U_rat_main = malloc(1000*sizeof(double));
   double *M_pe_main = malloc(1000*sizeof(double));
   double *halfmass_r_main = malloc(1000*sizeof(double));
   double *misc_54_main = malloc(1000*sizeof(double));
   double *misc_55_main = malloc(1000*sizeof(double));
   double *misc_56_main = malloc(1000*sizeof(double));
   double *misc_57_main = malloc(1000*sizeof(double));
   double *misc_58_main = malloc(1000*sizeof(double));
   double *misc_59_main = malloc(1000*sizeof(double));
   double *misc_60_main = malloc(1000*sizeof(double));
   
   //and an array to keep track of sub-mergers
   int *sub_merger_start = malloc(len_tree*sizeof(int));
   int *sub_merger_end = malloc(len_tree*sizeof(int));
   int *main_halo_start = malloc(len_tree*sizeof(int));
   int *main_halo_end = malloc(len_tree*sizeof(int));

   //fill out main halo info, and find the sub-merger indicies
   int len_main, number_of_subs;
   calc_subs(&len_main, &number_of_subs, a, ID, desc_scale, descID, num_prog, pid, Upid, desc_pid, phantom, SAM_mvir, mvir, rvir, r_s, 
   Vrms, mmp, last_MM_scale, Vmax, x_coord, y_coord, z_coord, x_vel, y_vel, z_vel, x_ang, y_ang, z_ang, spin, breadth_first_ID, 
   depth_first_ID, tree_root_ID, orig_halo_ID, snap_num, next_coprog_dep_ID, last_mainleaf_dep_ID, tidal_F, tidal_ID, r_s_klyp, 
   M200b_all, M200b_1, M200b_2, M200b_3, M200b_4, Xoff, Voff, spin_bull, b_to_a, c_to_a, A_x, A_y, A_z, T_U_rat, M_pe, halfmass_r, 
   misc_54, misc_55, misc_56, misc_57, misc_58, misc_59, misc_60, len_tree,
   a_main, ID_main, desc_scale_main, descID_main, num_prog_main, pid_main, Upid_main, desc_pid_main, phantom_main, SAM_mvir_main,
   mvir_main, rvir_main, r_s_main, Vrms_main, mmp_main, last_MM_scale_main, Vmax_main, x_coord_main, y_coord_main, z_coord_main,
   x_vel_main, y_vel_main, z_vel_main, x_ang_main, y_ang_main, z_ang_main, spin_main, breadth_first_ID_main, depth_first_ID_main,
   tree_root_ID_main, orig_halo_ID_main, snap_num_main, next_coprog_dep_ID_main, last_mainleaf_dep_ID_main, tidal_F_main, tidal_ID_main,
   r_s_klyp_main, M200b_all_main, M200b_1_main, M200b_2_main, M200b_3_main, M200b_4_main, Xoff_main, Voff_main, spin_bull_main, b_to_a_main,
   c_to_a_main, A_x_main, A_y_main, A_z_main, T_U_rat_main, M_pe_main, halfmass_r_main, misc_54_main, misc_55_main, misc_56_main, 
   misc_57_main, misc_58_main, misc_59_main, misc_60_main, sub_merger_start, sub_merger_end, main_halo_start, main_halo_end);

 //  printf("len_main is: %d\n", len_main);
 //  printf("num_subs is: %d\n", num_subs);

   //reallocate arrays to get rid of the extra elements
   a_main = realloc(a_main, len_main*sizeof(double));
   ID_main = realloc(ID_main, len_main*sizeof(double));
   desc_scale_main = realloc(desc_scale_main, len_main*sizeof(double));
   descID_main = realloc(descID_main, len_main*sizeof(double));
   num_prog_main = realloc(num_prog_main, len_main*sizeof(double));
   pid_main = realloc(pid_main, len_main*sizeof(double));
   Upid_main = realloc(Upid_main, len_main*sizeof(double));
   desc_pid_main = realloc(desc_pid_main, len_main*sizeof(double));
   phantom_main = realloc(phantom_main, len_main*sizeof(double));
   SAM_mvir_main = realloc(SAM_mvir_main, len_main*sizeof(double));
   mvir_main = realloc(mvir_main, len_main*sizeof(double));
   rvir_main = realloc(rvir_main, len_main*sizeof(double));
   r_s_main = realloc(r_s_main, len_main*sizeof(double));
   Vrms_main = realloc(Vrms_main, len_main*sizeof(double));
   mmp_main = realloc(mmp_main, len_main*sizeof(double));
   last_MM_scale_main = realloc(last_MM_scale_main, len_main*sizeof(double));
   Vmax_main = realloc(Vmax_main, len_main*sizeof(double));
   x_coord_main = realloc(x_coord_main, len_main*sizeof(double));
   y_coord_main = realloc(y_coord_main, len_main*sizeof(double));
   z_coord_main = realloc(z_coord_main, len_main*sizeof(double));
   x_vel_main = realloc(x_vel_main, len_main*sizeof(double));
   y_vel_main = realloc(y_vel_main, len_main*sizeof(double));
   z_vel_main = realloc(z_vel_main, len_main*sizeof(double));
   x_ang_main = realloc(x_ang_main, len_main*sizeof(double));
   y_ang_main = realloc(y_ang_main, len_main*sizeof(double));
   z_ang_main = realloc(z_ang_main, len_main*sizeof(double));
   spin_main = realloc(spin_main, len_main*sizeof(double));
   breadth_first_ID_main = realloc(breadth_first_ID_main, len_main*sizeof(double));
   depth_first_ID_main = realloc(depth_first_ID_main, len_main*sizeof(double));
   tree_root_ID_main = realloc(tree_root_ID_main, len_main*sizeof(double));
   orig_halo_ID_main = realloc(orig_halo_ID_main, len_main*sizeof(double));
   snap_num_main = realloc(snap_num_main, len_main*sizeof(double));
   next_coprog_dep_ID_main = realloc(next_coprog_dep_ID_main, len_main*sizeof(double));
   last_mainleaf_dep_ID_main = realloc(last_mainleaf_dep_ID_main, len_main*sizeof(double));
   tidal_F_main = realloc(tidal_F_main, len_main*sizeof(double));
   tidal_ID_main = realloc(tidal_ID_main, len_main*sizeof(double));
   r_s_klyp_main = realloc(r_s_klyp_main, len_main*sizeof(double));
   M200b_all_main = realloc(M200b_all_main, len_main*sizeof(double));
   M200b_1_main = realloc(M200b_1_main, len_main*sizeof(double));
   M200b_2_main = realloc(M200b_2_main, len_main*sizeof(double));
   M200b_3_main = realloc(M200b_3_main, len_main*sizeof(double));
   M200b_4_main = realloc(M200b_4_main, len_main*sizeof(double));
   Xoff_main = realloc(Xoff_main, len_main*sizeof(double));
   Voff_main = realloc(Voff_main, len_main*sizeof(double));
   spin_bull_main = realloc(spin_bull_main, len_main*sizeof(double));
   b_to_a_main = realloc(b_to_a_main, len_main*sizeof(double));
   c_to_a_main = realloc(c_to_a_main, len_main*sizeof(double));
   A_x_main = realloc(A_x_main, len_main*sizeof(double));
   A_y_main = realloc(A_y_main, len_main*sizeof(double));
   A_z_main = realloc(A_z_main, len_main*sizeof(double));
   T_U_rat_main = realloc(T_U_rat_main, len_main*sizeof(double));
   M_pe_main = realloc(M_pe_main, len_main*sizeof(double));
   halfmass_r_main = realloc(halfmass_r_main, len_main*sizeof(double));
   misc_54_main = realloc(misc_54_main, len_main*sizeof(double));
   misc_55_main = realloc(misc_55_main, len_main*sizeof(double));
   misc_56_main = realloc(misc_56_main, len_main*sizeof(double));
   misc_57_main = realloc(misc_57_main, len_main*sizeof(double));
   misc_58_main = realloc(misc_58_main, len_main*sizeof(double));
   misc_59_main = realloc(misc_59_main, len_main*sizeof(double));
   misc_60_main = realloc(misc_60_main, len_main*sizeof(double));
   sub_merger_start = realloc(sub_merger_start, number_of_subs*sizeof(int));
   sub_merger_end = realloc(sub_merger_end, number_of_subs*sizeof(int));
   main_halo_start = realloc(main_halo_start, number_of_subs*sizeof(int));
   main_halo_end = realloc(main_halo_end, number_of_subs*sizeof(int));
 
   //write the files
   if (halo_num % 100 == 0)
       printf("Writing sub/main files for halo %d\n", halo_num);
   FILE *out_file;
   char *filename = malloc(40*sizeof(char));
   snprintf(filename, sizeof(char) * 40, "tree_%d_%d_%d_halo_%d_submerger.txt",file_ident_1,file_ident_2,file_ident_3,halo_num);
   out_file = fopen(filename, "w");
   write_submerger(out_file, sub_merger_start, sub_merger_end, main_halo_start, main_halo_end, number_of_subs);
   fclose(out_file);
   free(filename);

   FILE *out_file2;
   char *filename2 = malloc(40*sizeof(char));
   snprintf(filename2, sizeof(char) * 40, "tree_%d_%d_%d_halo_%d_MainHalo.txt",file_ident_1,file_ident_2,file_ident_3,halo_num);
   out_file2 = fopen(filename2, "w");
   write_MainHalo(out_file2, a_main, ID_main, desc_scale_main, descID_main, num_prog_main, pid_main, Upid_main, desc_pid_main, 
              phantom_main, SAM_mvir_main, mvir_main, rvir_main, r_s_main, Vrms_main, mmp_main, last_MM_scale_main, Vmax_main, 
              x_coord_main, y_coord_main, z_coord_main, x_vel_main, y_vel_main, z_vel_main, x_ang_main, y_ang_main, z_ang_main, 
              spin_main, breadth_first_ID_main, depth_first_ID_main, tree_root_ID_main, orig_halo_ID_main, snap_num_main, 
              next_coprog_dep_ID_main, last_mainleaf_dep_ID_main, tidal_F_main, tidal_ID_main, r_s_klyp_main, 
              M200b_all_main, M200b_1_main, M200b_2_main, M200b_3_main, M200b_4_main, Xoff_main, Voff_main, spin_bull_main,
              b_to_a_main, c_to_a_main, A_x_main, A_y_main, A_z_main, T_U_rat_main, M_pe_main, halfmass_r_main, 
              misc_54_main, misc_55_main, misc_56_main, misc_57_main, misc_58_main, misc_59_main, misc_60_main, len_main);
   fclose(out_file2);
   free(filename2);

   
   free(sub_merger_start);
   free(sub_merger_end);
   free(main_halo_start);
   free(main_halo_end);

   //free the arrays
   free(a_main);
   free(ID_main);
   free(desc_scale_main);
   free(descID_main);
   free(num_prog_main);
   free(pid_main);
   free(Upid_main);
   free(desc_pid_main);
   free(phantom_main);
   free(SAM_mvir_main);
   free(mvir_main);
   free(rvir_main);
   free(r_s_main);
   free(Vrms_main);
   free(mmp_main);
   free(last_MM_scale_main);
   free(Vmax_main);
   free(x_coord_main);
   free(y_coord_main);
   free(z_coord_main);
   free(x_vel_main);
   free(y_vel_main);
   free(z_vel_main);
   free(x_ang_main);
   free(y_ang_main);
   free(z_ang_main);
   free(spin_main);
   free(breadth_first_ID_main);
   free(depth_first_ID_main);
   free(tree_root_ID_main);
   free(orig_halo_ID_main);
   free(snap_num_main);
   free(next_coprog_dep_ID_main);
   free(last_mainleaf_dep_ID_main);
   free(tidal_F_main);
   free(tidal_ID_main);
   free(r_s_klyp_main);
   free(M200b_all_main);
   free(M200b_1_main);
   free(M200b_2_main);
   free(M200b_3_main);
   free(M200b_4_main);
   free(Xoff_main);
   free(Voff_main);
   free(spin_bull_main);
   free(b_to_a_main);
   free(c_to_a_main);
   free(A_x_main);
   free(A_y_main);
   free(A_z_main);
   free(T_U_rat_main);
   free(M_pe_main);
   free(halfmass_r_main);
   free(misc_54_main);
   free(misc_55_main);
   free(misc_56_main);
   free(misc_57_main);
   free(misc_58_main);
   free(misc_59_main);
   free(misc_60_main);

    
}

int main(int argc, char ** argv)
{

   // check command line arguments
   if ( argc != 6 ) {
      printf("This program splits tree files into individuals, for n lines of halo data for one host\n");
      printf("Usage: ./tree_filesplitter filename (tree_treenum1_treenum2_treenum3.dat), # of lines, treenum1, treenum2, treenum3\n");
      exit(-1);
   }

   //parse input arguments
   FILE *tree_data = fopen(argv[1],"r");
   int num_lines = atoi(argv[2]);
   int file_ident_1 = atoi(argv[3]);
   int file_ident_2 = atoi(argv[4]);
   int file_ident_3 = atoi(argv[5]);
   
   int num_halos = 150000;
   //allocate arrays for all of the different parameters
   double *a = malloc(num_halos*sizeof(double));
   double *ID = malloc(num_halos*sizeof(double));
   double *desc_scale = malloc(num_halos*sizeof(double));
   double *descID = malloc(num_halos*sizeof(double));
   double *num_prog = malloc(num_halos*sizeof(double));
   double *pid = malloc(num_halos*sizeof(double));
   double *Upid = malloc(num_halos*sizeof(double));
   double *desc_pid = malloc(num_halos*sizeof(double));
   double *phantom = malloc(num_halos*sizeof(double));
   double *SAM_mvir = malloc(num_halos*sizeof(double));
   double *mvir = malloc(num_halos*sizeof(double));
   double *rvir = malloc(num_halos*sizeof(double));
   double *r_s = malloc(num_halos*sizeof(double));
   double *Vrms = malloc(num_halos*sizeof(double));
   double *mmp = malloc(num_halos*sizeof(double));
   double *last_MM_scale = malloc(num_halos*sizeof(double));
   double *Vmax = malloc(num_halos*sizeof(double));
   double *x_coord = malloc(num_halos*sizeof(double));
   double *y_coord = malloc(num_halos*sizeof(double));
   double *z_coord = malloc(num_halos*sizeof(double));
   double *x_vel = malloc(num_halos*sizeof(double));
   double *y_vel = malloc(num_halos*sizeof(double));
   double *z_vel = malloc(num_halos*sizeof(double));
   double *x_ang = malloc(num_halos*sizeof(double));
   double *y_ang = malloc(num_halos*sizeof(double));
   double *z_ang = malloc(num_halos*sizeof(double));
   double *spin = malloc(num_halos*sizeof(double));
   double *breadth_first_ID = malloc(num_halos*sizeof(double));
   double *depth_first_ID = malloc(num_halos*sizeof(double));
   double *tree_root_ID = malloc(num_halos*sizeof(double));
   double *orig_halo_ID = malloc(num_halos*sizeof(double));
   double *snap_num = malloc(num_halos*sizeof(double));
   double *next_coprog_dep_ID = malloc(num_halos*sizeof(double));
   double *last_mainleaf_dep_ID = malloc(num_halos*sizeof(double));
   double *tidal_F = malloc(num_halos*sizeof(double));
   double *tidal_ID = malloc(num_halos*sizeof(double));
   double *r_s_klyp = malloc(num_halos*sizeof(double));
   double *M200b_all = malloc(num_halos*sizeof(double));
   double *M200b_1= malloc(num_halos*sizeof(double));
   double *M200b_2 = malloc(num_halos*sizeof(double));
   double *M200b_3= malloc(num_halos*sizeof(double));
   double *M200b_4 = malloc(num_halos*sizeof(double));
   double *Xoff = malloc(num_halos*sizeof(double));
   double *Voff = malloc(num_halos*sizeof(double));
   double *spin_bull = malloc(num_halos*sizeof(double));
   double *b_to_a = malloc(num_halos*sizeof(double));
   double *c_to_a = malloc(num_halos*sizeof(double));
   double *A_x = malloc(num_halos*sizeof(double));
   double *A_y = malloc(num_halos*sizeof(double));
   double *A_z = malloc(num_halos*sizeof(double));
   double *T_U_rat = malloc(num_halos*sizeof(double));
   double *M_pe = malloc(num_halos*sizeof(double));
   double *halfmass_r = malloc(num_halos*sizeof(double));
   double *misc_54 = malloc(num_halos*sizeof(double));
   double *misc_55 = malloc(num_halos*sizeof(double));
   double *misc_56 = malloc(num_halos*sizeof(double));
   double *misc_57 = malloc(num_halos*sizeof(double));
   double *misc_58 = malloc(num_halos*sizeof(double));
   double *misc_59 = malloc(num_halos*sizeof(double));
   double *misc_60 = malloc(num_halos*sizeof(double));

   double *newtree_main_ids = malloc(num_halos*sizeof(double));
   int *subs_today = malloc(num_halos*sizeof(int));
   int *hosts_today = malloc(num_halos*sizeof(int));

  int i = 0;
  int j;
  double current_scale;
  int tree_num = 0;
  int tree_size = 0;
  int num_subs = 0;
  for (j=0; j<num_lines; j++) {
        if (i > num_halos) {
            i++;
            a = realloc(a, i*sizeof(double));
            ID = realloc(ID, i*sizeof(double));
            desc_scale = realloc(desc_scale, i*sizeof(double));
            descID = realloc(descID, i*sizeof(double));
            num_prog = realloc(num_prog, i*sizeof(double));
            pid = realloc(pid, i*sizeof(double));
            Upid = realloc(Upid, i*sizeof(double));
            desc_pid = realloc(desc_pid, i*sizeof(double));
            phantom = realloc(phantom, i*sizeof(double));
            SAM_mvir = realloc(SAM_mvir, i*sizeof(double));
            mvir = realloc(mvir, i*sizeof(double));
            rvir = realloc(rvir, i*sizeof(double));
            r_s = realloc(r_s, i*sizeof(double));
            Vrms = realloc(Vrms, i*sizeof(double));
            mmp = realloc(mmp, i*sizeof(double));
            last_MM_scale = realloc(last_MM_scale, i*sizeof(double));
            Vmax = realloc(Vmax, i*sizeof(double));
            x_coord = realloc(x_coord, i*sizeof(double));
            y_coord = realloc(y_coord, i*sizeof(double));
            z_coord = realloc(z_coord, i*sizeof(double));
            x_vel = realloc(x_vel, i*sizeof(double));
            y_vel = realloc(y_vel, i*sizeof(double));
            z_vel = realloc(z_vel, i*sizeof(double));
            x_ang = realloc(x_ang, i*sizeof(double));
            y_ang = realloc(y_ang, i*sizeof(double));
            z_ang = realloc(z_ang, i*sizeof(double));
            spin = realloc(spin, i*sizeof(double));
            breadth_first_ID = realloc(breadth_first_ID, i*sizeof(double));
            depth_first_ID = realloc(depth_first_ID, i*sizeof(double));
            tree_root_ID = realloc(tree_root_ID, i*sizeof(double));
            orig_halo_ID = realloc(orig_halo_ID, i*sizeof(double));
            snap_num = realloc(snap_num, i*sizeof(double));
            next_coprog_dep_ID = realloc(next_coprog_dep_ID, i*sizeof(double));
            last_mainleaf_dep_ID = realloc(last_mainleaf_dep_ID, i*sizeof(double));
            tidal_F = realloc(tidal_F, i*sizeof(double));
            tidal_ID = realloc(tidal_ID, i*sizeof(double));
            r_s_klyp = realloc(r_s_klyp, i*sizeof(double));
            M200b_all = realloc(M200b_all, i*sizeof(double));
            M200b_1= realloc(M200b_1, i*sizeof(double));
            M200b_2 = realloc(M200b_2, i*sizeof(double));
            M200b_3= realloc(M200b_3, i*sizeof(double));
            M200b_4 = realloc(M200b_4, i*sizeof(double));
            Xoff = realloc(Xoff, i*sizeof(double));
            Voff = realloc(Voff, i*sizeof(double));
            spin_bull = realloc(spin_bull, i*sizeof(double));
            b_to_a = realloc(b_to_a, i*sizeof(double));
            c_to_a = realloc(c_to_a, i*sizeof(double));
            A_x = realloc(A_x, i*sizeof(double));
            A_y = realloc(A_y, i*sizeof(double));
            A_z = realloc(A_z, i*sizeof(double));
            T_U_rat = realloc(T_U_rat, i*sizeof(double));
            M_pe = realloc(M_pe, i*sizeof(double));
            halfmass_r = realloc(halfmass_r, i*sizeof(double));
            misc_54 = realloc(misc_54, i*sizeof(double));
            misc_55 = realloc(misc_55, i*sizeof(double));
            misc_56 = realloc(misc_56, i*sizeof(double));
            misc_57 = realloc(misc_57, i*sizeof(double));
            misc_58 = realloc(misc_58, i*sizeof(double));
            misc_59 = realloc(misc_59, i*sizeof(double));
            misc_60 = realloc(misc_60, i*sizeof(double));
            i--;
        }
        fscanf(tree_data, "%lf", &a[i]);
        current_scale = a[i];
        // printf("Working on line: %d, scale is: %lf\n", j, current_scale);
        if (current_scale == 1.0000) {
           if (j > 0) {
              tree_num++;
              if (tree_num % 100 == 0) {
                  printf("Working on tree %d.....\n", tree_num);
              }
              find_main_subs(a, ID, desc_scale, descID, num_prog, pid, Upid, desc_pid, phantom, SAM_mvir, 
                mvir, rvir, r_s, Vrms, mmp, last_MM_scale, Vmax, x_coord, y_coord, z_coord, x_vel, y_vel, z_vel, 
                x_ang, y_ang, z_ang, spin, breadth_first_ID, depth_first_ID, tree_root_ID, orig_halo_ID, snap_num, 
                next_coprog_dep_ID, last_mainleaf_dep_ID, tidal_F, tidal_ID, r_s_klyp, M200b_all, M200b_1, M200b_2, 
                M200b_3, M200b_4, Xoff, Voff, spin_bull, b_to_a, c_to_a, A_x, A_y, A_z, T_U_rat, M_pe, halfmass_r, 
                misc_54, misc_55, misc_56, misc_57, misc_58, misc_59, misc_60, tree_size, &num_subs,
                tree_num, subs_today, hosts_today, newtree_main_ids, file_ident_1, file_ident_2, file_ident_3);
              tree_size = 0;
              i = 0;
              }
        }
        fscanf(tree_data, "%lf", &ID[i]);
        fscanf(tree_data, "%lf", &desc_scale[i]);
        fscanf(tree_data, "%lf", &descID[i]);
        fscanf(tree_data, "%lf", &num_prog[i]);
        fscanf(tree_data, "%lf", &pid[i]);
        fscanf(tree_data, "%lf", &Upid[i]);
        fscanf(tree_data, "%lf", &desc_pid[i]);
        fscanf(tree_data, "%lf", &phantom[i]);
        fscanf(tree_data, "%lf", &SAM_mvir[i]);
        fscanf(tree_data, "%lf", &mvir[i]);
        fscanf(tree_data, "%lf", &rvir[i]);
        fscanf(tree_data, "%lf", &r_s[i]);
        fscanf(tree_data, "%lf", &Vrms[i]);
        fscanf(tree_data, "%lf", &mmp[i]);
        fscanf(tree_data, "%lf", &last_MM_scale[i]);
        fscanf(tree_data, "%lf", &Vmax[i]);
        fscanf(tree_data, "%lf", &x_coord[i]);
        fscanf(tree_data, "%lf", &y_coord[i]);
        fscanf(tree_data, "%lf", &z_coord[i]);
        fscanf(tree_data, "%lf", &x_vel[i]);
        fscanf(tree_data, "%lf", &y_vel[i]);
        fscanf(tree_data, "%lf", &z_vel[i]);
        fscanf(tree_data, "%lf", &x_ang[i]);
        fscanf(tree_data, "%lf", &y_ang[i]);
        fscanf(tree_data, "%lf", &z_ang[i]);
        fscanf(tree_data, "%lf", &spin[i]);
        fscanf(tree_data, "%lf", &breadth_first_ID[i]);
        fscanf(tree_data, "%lf", &depth_first_ID[i]);
        fscanf(tree_data, "%lf", &tree_root_ID[i]);
        fscanf(tree_data, "%lf", &orig_halo_ID[i]);
        fscanf(tree_data, "%lf", &snap_num[i]);
        fscanf(tree_data, "%lf", &next_coprog_dep_ID[i]);
        fscanf(tree_data, "%lf", &last_mainleaf_dep_ID[i]);
        fscanf(tree_data, "%lf", &tidal_F[i]);
        fscanf(tree_data, "%lf", &tidal_ID[i]);
        fscanf(tree_data, "%lf", &r_s_klyp[i]);
        fscanf(tree_data, "%lf", &M200b_all[i]);
        fscanf(tree_data, "%lf", &M200b_1[i]);
        fscanf(tree_data, "%lf", &M200b_2[i]);
        fscanf(tree_data, "%lf", &M200b_3[i]);
        fscanf(tree_data, "%lf", &M200b_4[i]);
        fscanf(tree_data, "%lf", &Xoff[i]);
        fscanf(tree_data, "%lf", &Voff[i]);
        fscanf(tree_data, "%lf", &spin_bull[i]);
        fscanf(tree_data, "%lf", &b_to_a[i]);
        fscanf(tree_data, "%lf", &c_to_a[i]);
        fscanf(tree_data, "%lf", &A_x[i]);
        fscanf(tree_data, "%lf", &A_y[i]);
        fscanf(tree_data, "%lf", &A_z[i]);
        fscanf(tree_data, "%lf", &T_U_rat[i]);
        fscanf(tree_data, "%lf", &M_pe[i]);
        fscanf(tree_data, "%lf", &halfmass_r[i]);
        fscanf(tree_data, "%lf", &misc_54[i]);
        fscanf(tree_data, "%lf", &misc_55[i]);
        fscanf(tree_data, "%lf", &misc_56[i]);
        fscanf(tree_data, "%lf", &misc_57[i]);
        fscanf(tree_data, "%lf", &misc_58[i]);
        fscanf(tree_data, "%lf", &misc_59[i]);
        fscanf(tree_data, "%lf", &misc_60[i]);
        i++;
        tree_size++;
   }

   tree_num++;

//   printf("Writing file...\n");
   FILE *out_file2;
   char *filename2 = malloc(30*sizeof(char));
   snprintf(filename2, sizeof(char) * 30, "tree_%d_%d_%d_subs_today.txt",file_ident_1,file_ident_2,file_ident_3);
   out_file2 = fopen(filename2, "w");
   int q;
//   printf("Writing the subs today file...\n");
   for (q=0; q<num_subs; q++) {
        fprintf(out_file2, "%d ", subs_today[q]);
        fprintf(out_file2, "%d ", hosts_today[q]);
        fprintf(out_file2, "\n"); 
     }
   fclose(out_file2);
   free(filename2);

   fclose(tree_data); 
   free(subs_today);
   free(hosts_today);

   free(a);
   free(ID);
   free(desc_scale);
   free(descID);
   free(num_prog);
   free(pid);
   free(Upid);
   free(desc_pid);
   free(phantom);
   free(SAM_mvir);
   free(mvir);
   free(rvir);
   free(r_s);
   free(Vrms);
   free(mmp);
   free(last_MM_scale);
   free(Vmax);
   free(x_coord);
   free(y_coord);
   free(z_coord);
   free(x_vel);
   free(y_vel);
   free(z_vel);
   free(x_ang);
   free(y_ang);
   free(z_ang);
   free(spin);
   free(breadth_first_ID);
   free(depth_first_ID);
   free(tree_root_ID);
   free(orig_halo_ID);
   free(snap_num);
   free(next_coprog_dep_ID);
   free(last_mainleaf_dep_ID);
   free(tidal_F);
   free(tidal_ID);
   free(r_s_klyp);
   free(M200b_all);
   free(M200b_1);
   free(M200b_2);
   free(M200b_3);
   free(M200b_4);
   free(Xoff);
   free(Voff);
   free(spin_bull);
   free(b_to_a);
   free(c_to_a);
   free(A_x);
   free(A_y);
   free(A_z);
   free(T_U_rat);
   free(M_pe);
   free(halfmass_r);
   free(misc_54);
   free(misc_55);
   free(misc_56);
   free(misc_57);
   free(misc_58);
   free(misc_59);
   free(misc_60);

   return 0;     
}
