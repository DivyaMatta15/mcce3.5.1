#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"

PROT monte2_load_conflist(char *fname);
PROT monte2_load_pk1out(char *pk1out, char *dconf2);
void monte2_get_biglist(PROT prot);
void monte2_set_toggle(PROT prot);
PROT monte2_reduce(PROT prot);
int  monte2_load_pairwise(PROT prot);
void monte2_mc(PROT *prot_p);
int  monte2_check_toggle(PROT prot);
void zero_counters(PROT *prot_p);
void do_free_energy(PROT *prot_p);
double free_unf(PROT prot);
int curve_fitting(PROT prot, int n_titra);
int monte_out(PROT prot, int n_titra);
int   load1_pairwise();
void ENTROPY_CORRECTION_FUNCTION(PROT *prot);
double  **my_pairwise;
RES     **flip_res;
GROUP     **flip_group;
long    idum;
double  beta;
ISTATES *TheArray;
char **shead2;

int monte2() {
    // printf("monte2 I added the entropy correction function\n");
    PROT prot, prot_w, prot_red;
    int i_titra;
    float ph, eh, titr_pnt, temp, free_eng;
    float rmsd;
    int i_res, j_res, i_conf, i, i_cycle, n_cycle, n_cycle_min, n_cycle_max,n_cycle_chk,i_step;
    int i_red, j_red;
    int ic, jc, i_saved;
    FILE *fp;
	int r1=0;
	int r2=0;
	// int r=0; // Commented out by Salah on Dec. 2015 (not used)
    time_t   timer_start, timer_end;

    timer_start = time(NULL);
    flip_res = (RES **) malloc(env.monte_flips * sizeof(RES *));
    
    /* Load conflist */
    if (!env.monte_old_input) {
       // printf("before load conflist\n"); 
        prot = monte2_load_conflist(FN_CONFLIST3);
       // printf("after load conflist\n"); 
    }
    else {
        prot = monte2_load_pk1out("pk1.out","d_conf2.dat");
    }
    if (!prot.n_res) return USERERR;

    
    // Load my_pairwise
    printf("   Load pairwise interactions ...\n"); fflush(stdout);
    if (monte2_load_pairwise(prot)) {
        printf("   FATAL: pairwise interaction not loaded\n");
        return USERERR;
    }
    printf("   Done\n\n");
    fflush(stdout);
	
    monte2_check_toggle(prot);
    
    if (env.monte_seed < 0) idum = time(NULL);
    else idum = env.monte_seed;
    timer_end = time(NULL);
    
    fp = fopen(MC_OUT,"w"); 
    fprintf(fp,"random number seed = %ld\n",idum);
    fprintf(fp,"Monte Carlo set up time: %ld seconds.\n", timer_end-timer_start); fflush(stdout);
    fprintf(fp,"Starting number of residue = %d, number of conformer = %d\n",prot.n_res,prot.nc);
    fclose(fp);

    remove(DETAIL);

    for (i=0;i<env.monte_nstart;i++) ran2(&idum);
    for (ic=0;ic<prot.nc;ic++) {
        prot.conf[ic]->occ_table = (float *) malloc(env.titr_steps * sizeof(float));
    }
    
    /* pH/Eh titration */
    for (i_titra=0;i_titra<env.titr_steps;i_titra++) {
        timer_start = time(NULL);
        ph = env.titr_ph0;
        eh = env.titr_eh0;
        if (env.titr_type == 'p') ph = env.titr_ph0 + ((float)i_titra)*env.titr_phd;
        else eh = env.titr_eh0 + ((float)i_titra)*env.titr_ehd;
        fp = fopen(MC_OUT,"a"); fprintf(fp,"\npH = %6.2f  Eh = %6.2f\n", ph, eh); fclose(fp);
        
        for (ic=0; ic<prot.nc; ic++) {
            prot.conf[ic]->E_ph =
            env.monte_temp/ROOMT * prot.conf[ic]->H * (ph-prot.conf[ic]->pKa) * PH2KCAL;
            prot.conf[ic]->E_eh =
            env.monte_temp/ROOMT * prot.conf[ic]->e * (eh-prot.conf[ic]->Em) * PH2KCAL/58.0;
            
            prot.conf[ic]->E_self = prot.conf[ic]->E_vdw0
            + prot.conf[ic]->E_vdw1
            + prot.conf[ic]->E_epol
            + prot.conf[ic]->E_tors
            + prot.conf[ic]->E_dsolv
            + prot.conf[ic]->E_extra
            + prot.conf[ic]->E_ph
            + prot.conf[ic]->E_eh
            + prot.conf[ic]->E_TS
            - prot.conf[ic]->PE_TS;
        }
        
        prot.E_free_unfold = free_unf(prot);
        
        prot_w = monte2_reduce(prot);
        for (ic=0;ic<prot_w.nc;ic++) prot_w.conf[ic]->occ_table = (float *) malloc(sizeof(float));
        
        i_red=0;
        while(1) {
            if (env.monte_n_red > 0) {
                if (i_red >= env.monte_n_red) break;
            }
            if (!prot_w.nc) break;
            
            if (i_red>=2) {
                rmsd = 0.;
                for (ic=0;ic<prot_w.nc;ic++) {
                    float ave_occ = 0.;
                    for (j_red=0;j_red<i_red;j_red++) {
                        ave_occ += prot_w.conf[ic]->occ_table[j_red];
                    }
                    ave_occ = ave_occ/(float) i_red;
                    rmsd += (prot_w.conf[ic]->occ_table[i_red-1] - ave_occ)/(float) (i_red-1.)
                    * (prot_w.conf[ic]->occ_table[i_red-1] - ave_occ)/(float) (i_red-1.);
                }
                rmsd = sqrt(rmsd/(float)prot_w.nc);
                fp = fopen(MC_OUT,"a"); fprintf(fp,"Reduce = %3d, rmsd = %.5f\n",i_red,rmsd); fclose(fp);
                if (rmsd < env.monte_converge) break;
            }
            
            if (i_red) monte2_set_toggle(prot_w);
            prot_red = monte2_reduce(prot_w);
            if (!prot_red.nc) break;
            monte2_get_biglist(prot_red);
            i_red++;
            for (ic=0;ic<prot_w.nc;ic++) prot_w.conf[ic]->occ_table= (float *) realloc(prot_w.conf[ic]->occ_table,i_red*sizeof(float));
            
            fp = fopen(MC_OUT,"a"); fprintf(fp,"\nReduce = %3d, number of active residue: %4d, active conformer %5d\n",i_red,prot_red.n_res,prot_red.nc); fclose(fp);
            
            /* Initialize */
            for (i_res=0;i_res<prot_red.n_res;i_res++) {
                i_conf = (int) (ran2(&idum) * (prot_red.res[i_res].n_conf - 1.)) + 1;
                prot_red.res[i_res].conf_w = &prot_red.res[i_res].conf[i_conf];
            }
            zero_counters(&prot_red);
            
            prot_red.E_state = prot_red.E_base;
            for (i_res=0;i_res<prot_red.n_res;i_res++) {
                prot_red.E_state += prot_red.res[i_res].conf_w->E_self;
                ic = prot_red.res[i_res].conf_w->i_conf_prot;
                for (j_res=i_res+1;j_res<prot_red.n_res;j_res++) {
                    jc = prot_red.res[j_res].conf_w->i_conf_prot;
                    prot_red.E_state += my_pairwise[ic][jc];
                }
            }
            prot_red.E_min = prot_red.E_state;

            /* Equilibrating */
            n_cycle = (prot_red.nc*env.monte_neq - 1)/env.monte_niter_cycle + 1;
            if (env.anneal_nstep) beta = KCAL2KT/(env.anneal_temp_start/ROOMT);
            else beta = KCAL2KT/(env.monte_temp/ROOMT);
            fp = fopen(MC_OUT,"a"); fprintf(fp,"Equilibrating:   %5d cycle(s) at temperature = %10.2fK\n",n_cycle,env.anneal_temp_start); fclose(fp);
            // printf("%d",n_cycle);
            for (i_cycle=0;i_cycle<n_cycle;i_cycle++) {
                monte2_mc(&prot_red);
                //if(i_red==0)//my code start here
            }
            /* Set zero */
            zero_counters(&prot_red);
            
            /* Annealling */
            n_cycle = (prot_red.nc*env.anneal_niter_step - 1)/env.monte_niter_cycle + 1;
            for (i_step=0;i_step<env.anneal_nstep;i_step++) {
                temp = env.anneal_temp_start + ((float)(i_step+1)/(float) env.anneal_nstep) * (env.monte_temp - env.anneal_temp_start);
                beta = KCAL2KT/(temp/ROOMT);
                fp = fopen(MC_OUT,"a"); fprintf(fp,"Annealling:      %5d cycle(s) at temperature = %10.2fK\n",n_cycle,temp); fclose(fp);
                for (i_cycle=0;i_cycle<n_cycle;i_cycle++)                 {
                    monte2_mc(&prot_red);
                    ENTROPY_CORRECTION_FUNCTION(&prot_red);
                }
            }
            zero_counters(&prot_red);
            
            /* Data */
            n_cycle_min = (prot_red.nc*env.monte_niter_min - 1)/env.monte_niter_cycle + 1;
            if (prot_red.nc*env.monte_niter_max >= 0)
                n_cycle_max = (prot_red.nc*env.monte_niter_max - 1)/env.monte_niter_cycle + 1;
            else n_cycle_max = -1;
            n_cycle_chk = (prot_red.nc*env.monte_niter_chk - 1)/env.monte_niter_cycle + 1;
            beta = KCAL2KT/(env.monte_temp/ROOMT);
            fp = fopen(MC_OUT,"a");
            fprintf(fp,"Collecting data:                at temperature = %10.2fK\n",env.monte_temp);
            fprintf(fp,"   Minimum/maximum number of cycles: %5d/%5d\n",n_cycle_min,n_cycle_max);
            fclose(fp);
            i_cycle = 0;
            while (1) {
                if (i_cycle >= n_cycle_min) 
                {
                    if (n_cycle_max != -1) {
                        if (i_cycle >= n_cycle_max) break;
                    }
                    if (!fmod(i_cycle,n_cycle_chk)) {

                    for(r1=0;r1<prot_red.n_res;r1++)
                    for(r2=0;r2<prot_red.res[r1].num_of_states;r2++)
                    prot_red.res[r1].IS[r2].totalOCC=0;
                    for (ic=0;ic<prot_red.nc;ic++) 
                    {
                        prot_red.conf[ic]->occ_old = prot_red.conf[ic]->occ;
                        prot_red.conf[ic]->occ = (float) prot_red.conf[ic]->counter_accept / (float) (i_cycle * env.monte_niter_cycle);
                    } 
                        

                        rmsd = 0.;
                        for (ic=0;ic<prot_red.nc;ic++) {
                            rmsd += (prot_red.conf[ic]->occ - prot_red.conf[ic]->occ_old) * (prot_red.conf[ic]->occ - prot_red.conf[ic]->occ_old);
                        }
                        rmsd = sqrt(rmsd/(float) prot_red.nc);
                        //printf("%4d, %8.3f\n",i_cycle,rmsd);
                        if (rmsd < env.monte_converge) break;
                    }
               
                }   
                i_cycle++;
                monte2_mc(&prot_red);   
                ENTROPY_CORRECTION_FUNCTION(&prot_red);
                for(r1=0;r1<prot_red.n_res;r1++)
                for(r2=0;r2<prot_red.res[r1].num_of_states;r2++)
                prot_red.res[r1].IS[r2].totalOCC=0;
            }

            fp = fopen(MC_OUT,"a"); fprintf(fp,"   Actual  number of cycles to converge: %5d\n",i_cycle); fclose(fp);
            
            for (ic=0;ic<prot_w.nc;ic++) {
                prot_w.conf[ic]->counter_trial = 0;
            }
            
            jc = 0;
            for (ic=0;ic<prot_red.nc;ic++) {
                prot_red.conf[ic]->occ = (float) prot_red.conf[ic]->counter_accept / (float) (i_cycle * env.monte_niter_cycle);
                while (strcmp(prot_red.conf[ic]->uniqID, prot_w.conf[jc]->uniqID)) jc++;
                prot_w.conf[jc]->occ = prot_red.conf[ic]->occ;
                prot_w.conf[jc]->counter_trial = prot_red.conf[ic]->counter_trial;
            }
            //ENTROPY_CORRECTION_FUNCTION(&prot_red);
            for (ic=0;ic<prot_w.nc;ic++) {
                prot_w.conf[ic]->occ_table[i_red-1] = prot_w.conf[ic]->occ;
            }
            
            
            fp = fopen(DETAIL,"a");
            fprintf(fp,"\n");
            fprintf(fp,"pH = %6.2f Eh = %6.2f, Unf. Energy = %8.2f Kcal/mol\n", ph, eh, prot.E_free_unfold);
            fprintf(fp,"pH = %6.2f Eh = %6.2f, Ave. Energy = %8.2f Kcal/mol\n", ph, eh, prot_red.E_accum / (float) (i_cycle * env.monte_niter_cycle));
            fprintf(fp,"pH = %6.2f Eh = %6.2f, Min. Energy = %8.2f Kcal/mol\n", ph, eh, prot_red.E_min);
            if (env.monte_do_energy) {
                free_eng = 0;
                for (i_saved=0;i_saved<prot_red.n_saved;i_saved++) {
                    free_eng += exp(-beta*(prot_red.saved_states[i_saved].E - prot_red.E_min));
                }
                free_eng = - (1./beta) * log(free_eng) + prot_red.E_min;
                fprintf(fp,"%d unique states collected\n", prot_red.n_saved);
                fprintf(fp,"Free Energy = %8.2f Kcal/mol\n", free_eng);
            }
            if (env.titr_type == 'p') titr_pnt = ph;
            else titr_pnt = eh;
            for (ic=0;ic<prot_w.nc;ic++) {
                fprintf(fp,"%6.2f %s occ=%5.3f trial=%8d/%8d %c\n",
                titr_pnt,
                prot_w.conf[ic]->uniqID,
                prot_w.conf[ic]->occ,
                prot_w.conf[ic]->counter_trial,
                i_cycle * env.monte_niter_cycle,
                prot_w.conf[ic]->toggle);
            }
            fclose(fp);
            
            //printf("free prot_red, biglist\n");
            for (i_res=0;i_res<prot_red.n_res;i_res++) {
                if (prot_red.res[i_res].n_ngh) free(prot_red.res[i_res].ngh);
            }
            //printf("free prot_red, saved_states\n");
            free(prot_red.saved_states);
            //printf("free prot_red, conf\n");
            if (prot_red.nc) free(prot_red.conf);
            //printf("free prot_red\n");
            del_prot(&prot_red);
            }
        
            jc = 0;
            for (ic=0;ic<prot_w.nc;ic++) {
                while (strcmp(prot_w.conf[ic]->uniqID, prot.conf[jc]->uniqID)) jc++;
                
                prot.conf[jc]->occ = 0.;
                for (j_red=0;j_red<i_red;j_red++) {
                    prot.conf[jc]->occ += prot_w.conf[ic]->occ_table[j_red];
                }
                prot.conf[jc]->occ = prot.conf[jc]->occ/(float) i_red;
            }
            for (ic=0;ic<prot.nc;ic++) {
                prot.conf[ic]->occ_table[i_titra] = prot.conf[ic]->occ;
            }
            
            //printf("free prot_w, occ_table\n");
            for (ic=0;ic<prot_w.nc;ic++) free(prot_w.conf[ic]->occ_table);
            //printf("free prot_w, conf\n");
            if (prot_w.nc) free(prot_w.conf);
            del_prot(&prot_w);
            timer_end = time(NULL);
            fp = fopen(MC_OUT,"a"); fprintf(fp,"Monte Carlo running time: %ld seconds.\n", timer_end-timer_start); fclose(fp);
            monte_out(prot, i_titra+1);
    }
    
    monte_out(prot, env.titr_steps);
    
	// curve fitting 
    //printf("   Fit titration curves to get pKa/Em ...\n");
    //fflush(stdout);
    //curve_fitting(prot,env.titr_steps);
    
    
    //printf("free prot, occ_table\n");
    for (ic=0;ic<prot.nc;ic++) free(prot.conf[ic]->occ_table);
    //printf("free prot, conf\n");
    if (prot.nc) free(prot.conf);
    del_prot(&prot);
    
    //printf("free flip_res\n");
    free(flip_res);
    return 0;
}

PROT monte2_load_conflist(char *fname) {
    PROT prot;
    FILE *fp;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    CONF conf;
    int  k_res, k_conf;
    
    memset(&prot,0,sizeof(PROT));
    if (!(fp=fopen(fname, "r"))) {
        printf("   FATAL: Can't open file %s\n", fname);fflush(stdout);
        return prot;
    }
    
    fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
    
    while(fgets(sbuff, sizeof(sbuff), fp)) {
        if (strlen(sbuff) < 20) continue;
        
        /*
        iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    self
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        00001 NTR01_0001_001 f 0.00  0.000     0  0.00  0  0  -0.002   1.335   0.000   0.005   0.000   0.000 01O000M000t
        */
        
        /* load this line to a conf */
        memset(&conf,0,sizeof(CONF));
        strncpy(stemp, sbuff,    5); stemp[5] = '\0'; conf.i_conf_prot = atoi(stemp);
        strncpy(conf.uniqID, sbuff+6, 14); conf.uniqID[14]  = '\0';
        conf.toggle = sbuff[21];
        strncpy(stemp, sbuff+23, 4); stemp[4] = '\0'; conf.occ     = atof(stemp);
        strncpy(stemp, sbuff+28, 6); stemp[6] = '\0'; conf.netcrg  = atof(stemp);
        strncpy(stemp, sbuff+35, 5); stemp[5] = '\0'; conf.Em      = atof(stemp);
        strncpy(stemp, sbuff+41, 5); stemp[5] = '\0'; conf.pKa     = atof(stemp);
        strncpy(stemp, sbuff+47, 2); stemp[2] = '\0'; conf.e       = atof(stemp);
        strncpy(stemp, sbuff+50, 2); stemp[2] = '\0'; conf.H       = atof(stemp);
        strncpy(stemp, sbuff+53, 7); stemp[7] = '\0'; conf.E_vdw0  = atof(stemp) * env.scale_vdw0;
        strncpy(stemp, sbuff+61, 7); stemp[7] = '\0'; conf.E_vdw1  = atof(stemp) * env.scale_vdw1;
        strncpy(stemp, sbuff+69, 7); stemp[7] = '\0'; conf.E_tors  = atof(stemp) * env.scale_tor;
        strncpy(stemp, sbuff+77, 7); stemp[7] = '\0'; conf.E_epol  = atof(stemp) * env.scale_ele;
        strncpy(stemp, sbuff+85, 7); stemp[7] = '\0'; conf.E_dsolv  = atof(stemp)* env.scale_dsolv;
        strncpy(stemp, sbuff+93, 7); stemp[7] = '\0'; conf.E_extra  = atof(stemp);
        strncpy(conf.history, sbuff+101, 11); conf.history[11]  = '\0';

        strncpy(conf.resName,  conf.uniqID, 3); conf.resName[3] = '\0';
        strncpy(conf.confName, conf.uniqID, 5); conf.confName[5] = '\0';
        conf.chainID = conf.uniqID[5];
        strncpy(stemp, conf.uniqID+6, 4); stemp[4] = '\0'; conf.resSeq = atoi(stemp);
        conf.iCode = conf.uniqID[10];

        /* search for the residue: using same procedure as in load_pdb() */
        for (k_res = prot.n_res - 1; k_res >= 0; k_res--) {
            if (!strcmp(conf.resName, prot.res[k_res].resName) &&
                conf.chainID == prot.res[k_res].chainID &&
            conf.resSeq  == prot.res[k_res].resSeq  &&
            conf.iCode   == prot.res[k_res].iCode  )
            {
                break;
            }
        }
        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(&prot, prot.n_res);
            strcpy(prot.res[k_res].resName, conf.resName);
            prot.res[k_res].chainID = conf.chainID;
            prot.res[k_res].resSeq  = conf.resSeq;
            prot.res[k_res].iCode   = conf.iCode;
            
            /* Insert a backbone conformer */
            ins_conf(&prot.res[k_res], 0, 0);
        }
        
        k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, 0);
        prot.res[k_res].conf[k_conf] = conf;
        if (conf.toggle == 't') prot.res[k_res].conf[k_conf].toggle = 'f';
        else if (conf.toggle == 'f') prot.res[k_res].conf[k_conf].toggle = 't';
        else if (conf.toggle == 'a') prot.res[k_res].conf[k_conf].toggle = 'a';
        else prot.res[k_res].conf[k_conf].toggle = 't';
    }
    fclose(fp);
    
    for (k_res=0;k_res<prot.n_res;k_res++) {
        for (k_conf=1;k_conf<prot.res[k_res].n_conf;k_conf++) {
            prot.nc++;
            prot.conf = (CONF**) realloc(prot.conf, prot.nc*sizeof(void *));
            prot.res[k_res].conf[k_conf].i_conf_prot = prot.nc-1;
            prot.conf[prot.res[k_res].conf[k_conf].i_conf_prot] = &prot.res[k_res].conf[k_conf];
        }
    }
    return prot;
}

int monte2_load_pairwise(PROT prot) {
    int   ic, jc;
    int   i_res,i_conf,j_res,j_conf;
    int   n_miss_jc = 0,i_miss_ic,j_miss_ic,i_miss_jc,j_miss_jc;
    int   *n_miss_ic = NULL, **miss_ic = NULL, *miss_jc = NULL;
    
	
	//int i, j, kc; // Salah added on Dec. 2015
    EMATRIX ematrix; // Salah added on Dec. 2015
	// Salah added on Dec. 2015 =====
	ematrix.n = 0;  
    if (load_energies(&ematrix, ".", 0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }
	
    /* declare memory */
    if (!(my_pairwise = (double **) malloc(prot.nc * sizeof(double *)))) {
        printf("   FATAL: memory error in load_pairwise()\n");fflush(stdout);
        return USERERR;
    }
    for (ic=0; ic<prot.nc; ic++) {
        if (!(my_pairwise[ic] = (double *) malloc(prot.nc * sizeof(double)))) {
            printf("   FATAL: memory error in load_pairwise()\n");fflush(stdout);
            return USERERR;
        }
        memset(my_pairwise[ic],0,prot.nc * sizeof(double));
    }
    
    for (ic=0; ic<ematrix.n; ic++) {
       for (jc=0; jc<ematrix.n; jc++) {
		    if (ematrix.pw[ic][jc].vdw < 500)
                my_pairwise[ic][jc] = ematrix.pw[ic][jc].ele*env.scale_ele + ematrix.pw[ic][jc].vdw*env.scale_vdw;
            else
                my_pairwise[ic][jc] = ematrix.pw[ic][jc].ele*env.scale_ele + ematrix.pw[ic][jc].vdw;
			
		}
    }
	// free memory
    free_ematrix(&ematrix);
        
        
    
    if (n_miss_jc) {
        for (i_miss_jc=0;i_miss_jc<n_miss_jc;i_miss_jc++) {
            jc = miss_jc[i_miss_jc];
            for (i_miss_ic=0;i_miss_ic<n_miss_ic[i_miss_jc];i_miss_ic++) {
                ic = miss_ic[i_miss_jc][i_miss_ic];
                
                for (j_miss_jc=0;j_miss_jc<n_miss_jc;j_miss_jc++) {
                    if (ic == miss_jc[j_miss_jc]) {
                        for (j_miss_ic=0;j_miss_ic<n_miss_ic[j_miss_jc];j_miss_ic++) {
                            if (jc == miss_ic[j_miss_jc][j_miss_ic]) break;
                        }
                        if (j_miss_ic<n_miss_ic[j_miss_jc]) {
                            printf("   FATAL: Interaction between %s and %s does not exist both opp file\n", prot.conf[ic]->uniqID, prot.conf[jc]->uniqID);fflush(stdout);
                            return USERERR;
                        }
                    }
                }
                my_pairwise[ic][jc] = my_pairwise[jc][ic];
            }
            free(miss_ic[i_miss_jc]);
            n_miss_ic[i_miss_jc] = 0;
        }
        free(n_miss_ic);
        free(miss_jc);
    }
    
    /* average or pick smaller symetric elements, self to sel set to 0 */
    printf("   WARNING: Big difference (> %.1f%%) is reported\n", env.warn_pairwise);
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            ic = prot.res[i_res].conf[i_conf].i_conf_prot;
            for (j_res=0; j_res<prot.n_res; j_res++) {
                for (j_conf=1; j_conf<prot.res[j_res].n_conf; j_conf++) {
                    jc = prot.res[j_res].conf[j_conf].i_conf_prot;

                    if (i_res != j_res) {
                        if (200.*fabs(my_pairwise[ic][jc]-my_pairwise[jc][ic])/(my_pairwise[ic][jc]+my_pairwise[jc][ic]) > env.warn_pairwise &&
                            fabs(my_pairwise[ic][jc] - my_pairwise[jc][ic]) > 0.5) {
                            if ( (my_pairwise[ic][jc]+my_pairwise[jc][ic]) > 0.2 )
                                printf("%s -> %s = %10.3f, <- = %10.3f\n", prot.conf[ic]->uniqID, prot.conf[jc]->uniqID, my_pairwise[ic][jc], my_pairwise[jc][ic]);
                        }
                        if (env.average_pairwise) {   /* average */
                            my_pairwise[ic][jc] = my_pairwise[jc][ic] = 0.5*(my_pairwise[ic][jc]+my_pairwise[jc][ic]);
                        }
                        else {                        /* get the smaller */
                            my_pairwise[ic][jc] = my_pairwise[jc][ic] = fabs(my_pairwise[ic][jc]) > fabs(my_pairwise[jc][ic]) ? my_pairwise[jc][ic] : my_pairwise[ic][jc];
                        }
                    }
                    else {
                        my_pairwise[ic][jc] = my_pairwise[jc][ic] = 0.0; /* Necessary for fast E updating */
                    }
                }
            }
        }
    }
    return 0;
}

PROT monte2_load_pk1out(char *pk1out, char *dconf2) {
    PROT prot;
    FILE *fp1,*fp2;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    CONF conf;
    int  k_res, k_conf;
    int  i,ic,ic_start;
    
    memset(&prot,0,sizeof(PROT));
    if (!(fp1=fopen(pk1out, "r"))) {
        printf("   FATAL: Can't open file %s\n", pk1out);fflush(stdout);
        return prot;
    }

    if (!(fp2=fopen(dconf2, "r"))) {
        printf("   FATAL: Can't open file %s\n", dconf2);fflush(stdout);
        return prot;
    }

    for (i=0;i<11;i++) {
        fgets(sbuff, sizeof(sbuff), fp1); /* skip the first 11 lines */
        if (!strncmp(sbuff,"                  pk0   d_solv  pol    pk1  hC eC   nonel   Em    dg_i    dg_0",78)) break;
    }

    while(fgets(sbuff, sizeof(sbuff), fp1)) {
        if (strlen(sbuff) < 20) continue;

        /*
        .iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    self
        .01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        . 972 DCU X 003A   0.00   4.11  -2.87   0.00  0  0 2207.72    0. 2208.97 2208.97
        . 973 DC+ X 003B   0.00   3.76  -0.09   0.00  0  1 2207.72  260. 2215.88 2215.88
        .00000 NTG01_0001_001 f 0.00  0.000     0  0.00  0  0  -0.027  -0.418   0.000   0.004  -0.000   0.000  -0.441
        .                  pk0   d_solv  pol    pk1  hC eC   nonel   Em    dg_i    dg_0
        .   4 ARG A 007A   0.00   0.84  -0.57   0.00  0  0   -0.60    0.   -0.33   -0.33
        .PHE A 018A   1 f 0.00  1  f
        */

        /* load this line to a conf */
        memset(&conf,0,sizeof(CONF));
        //strncpy(stemp, sbuff,    4); stemp[4] = '\0'; conf.i_conf_prot = atoi(stemp);
        strncpy(conf.uniqID, sbuff+5, 10); conf.uniqID[10]  = '\0';
        strncpy(stemp, sbuff+16, 6); stemp[6] = '\0'; conf.pKa     = atof(stemp);
        strncpy(stemp, sbuff+23, 6); stemp[6] = '\0'; conf.E_dsolv  = atof(stemp) * PH2KCAL;
        strncpy(stemp, sbuff+30, 6); stemp[6] = '\0'; conf.E_epol  = atof(stemp) * PH2KCAL * env.scale_ele;
        strncpy(stemp, sbuff+44, 2); stemp[2] = '\0'; conf.H       = atof(stemp);
        strncpy(stemp, sbuff+47, 2); stemp[2] = '\0'; conf.e       = -atof(stemp);
        strncpy(stemp, sbuff+50, 7); stemp[7] = '\0'; conf.E_vdw1  = atof(stemp) * PH2KCAL * env.scale_vdw;
        strncpy(stemp, sbuff+58, 5); stemp[5] = '\0'; conf.Em      = atof(stemp);
        strncpy(stemp, sbuff+72, 7); stemp[7] = '\0'; conf.E_self0 = atof(stemp) * PH2KCAL;

        strncpy(conf.resName,  conf.uniqID, 3); conf.resName[3] = '\0';
        strncpy(conf.confName, conf.uniqID, 3); conf.confName[3] = '\0';
        conf.chainID = conf.uniqID[4];
        strncpy(stemp, conf.uniqID+6, 3); stemp[3] = '\0'; conf.resSeq = atoi(stemp);
        conf.iCode = ' ';
        conf.netcrg = conf.H - conf.e;
        
        /* search for the residue: using same procedure as in load_pdb() */
        for (k_res = prot.n_res - 1; k_res >= 0; k_res--) {
            if (!strncmp(conf.resName, prot.res[k_res].resName,2) &&
                conf.chainID == prot.res[k_res].chainID &&
            conf.resSeq  == prot.res[k_res].resSeq  &&
            conf.iCode   == prot.res[k_res].iCode  )
            {
                break;
            }
        }
        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(&prot, prot.n_res);
            strcpy(prot.res[k_res].resName, conf.resName);
            prot.res[k_res].chainID = conf.chainID;
            prot.res[k_res].resSeq  = conf.resSeq;
            prot.res[k_res].iCode   = conf.iCode;
            
            /* Insert a backbone conformer */
            ins_conf(&prot.res[k_res], 0, 0);
        }
        
        k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, 0);
        prot.res[k_res].conf[k_conf] = conf;
    }
    fclose(fp1);
    
    /* Add DOD */
    for (k_res=0;k_res<prot.n_res;k_res++) {
        if (strcmp(prot.res[k_res].resName,"HOH")
            && strcmp(prot.res[k_res].resName,"_CL")) continue;
        
        k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, 0);
        
        strcpy(prot.res[k_res].conf[k_conf].resName, prot.res[k_res].resName);
        strcpy(prot.res[k_res].conf[k_conf].confName, prot.res[k_res].resName);
        prot.res[k_res].conf[k_conf].confName[2] = 'D';
        prot.res[k_res].conf[k_conf].confName[3] = '\0';
        prot.res[k_res].conf[k_conf].chainID = prot.res[k_res].chainID;
        prot.res[k_res].conf[k_conf].resSeq  = prot.res[k_res].resSeq;
        prot.res[k_res].conf[k_conf].iCode = ' ';
        prot.res[k_res].conf[k_conf].toggle = 'a';
        
        /* uniqID */
        strncpy(prot.res[k_res].conf[k_conf].uniqID, prot.res[k_res].conf[k_conf].confName, 3);
        prot.res[k_res].conf[k_conf].uniqID[3] = ' ';
        prot.res[k_res].conf[k_conf].uniqID[4] = prot.res[k_res].conf[k_conf].chainID;
        prot.res[k_res].conf[k_conf].uniqID[5] = ' ';
        sprintf(prot.res[k_res].conf[k_conf].uniqID+6,"%03d",prot.res[k_res].conf[k_conf].resSeq);
        prot.res[k_res].conf[k_conf].uniqID[9] = 'z';
        prot.res[k_res].conf[k_conf].uniqID[10]  = '\0';
    }
    
    for (k_res=0;k_res<prot.n_res;k_res++) {
        for (k_conf=1;k_conf<prot.res[k_res].n_conf;k_conf++) {
            prot.nc++;
            prot.conf = (CONF**) realloc(prot.conf, prot.nc*sizeof(void *));
            prot.res[k_res].conf[k_conf].i_conf_prot = prot.nc-1;
            prot.conf[prot.nc-1] = &prot.res[k_res].conf[k_conf];
        }
    }
    
    for (ic=0;ic<prot.nc;ic++) prot.conf[ic]->on = 0;
    ic = 0;
    while(fgets(sbuff, sizeof(sbuff), fp2)) {
        if (strlen(sbuff) < 20) continue;

        /*
        .01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        .PHE A 018A   1 f 0.00  1  f
        */
        strncpy(conf.uniqID, sbuff, 10); conf.uniqID[10]  = '\0';
        conf.toggle = sbuff[15];
        strncpy(stemp, sbuff+17, 4); stemp[4] = '\0'; conf.occ     = atof(stemp);

        ic_start = ic;
        while (strcmp(conf.uniqID, prot.conf[ic]->uniqID)) {
            ic++;
            if (ic >= prot.nc) ic = 0;
            if (ic == ic_start) break;
        }
        if (!strcmp(conf.uniqID, prot.conf[ic]->uniqID)) {
            if (conf.toggle == 'f') prot.conf[ic]->toggle = 't';
            else if (conf.toggle == 't') prot.conf[ic]->toggle = 'f';
            else prot.conf[ic]->toggle = conf.toggle;
            
            prot.conf[ic]->occ = conf.occ;
            prot.conf[ic]->on = 1;
        }
    }
    fclose(fp2);
    
    return prot;
}

void monte2_get_biglist(PROT prot) {
    int i_res, j_res, i_conf, j_conf, ic, jc;
    int added;
    
    for (i_res=0;i_res<prot.n_res;i_res++) {
        prot.res[i_res].n_ngh=0;
        
        for (j_res=0;j_res<prot.n_res;j_res++) {
            if (i_res == j_res) continue;
            
            added = 0;
            for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                ic = prot.res[i_res].conf[i_conf].i_conf_prot;
                for (j_conf=1;j_conf<prot.res[j_res].n_conf;j_conf++) {
                    jc = prot.res[j_res].conf[j_conf].i_conf_prot;
                    
                    if (fabs(my_pairwise[ic][jc])>env.big_pairwise) {
                        prot.res[i_res].n_ngh++;
                        prot.res[i_res].ngh = (RES_STRUCT**) realloc(prot.res[i_res].ngh, prot.res[i_res].n_ngh*sizeof(void *));
                        prot.res[i_res].ngh[prot.res[i_res].n_ngh-1] = &prot.res[j_res];
                        added = 1;
                        break;
                    }
                }
                if (added) break;
            }
        }
        
        prot.res[i_res].n_flip_max = (prot.res[i_res].n_ngh+1) < env.monte_flips ? (prot.res[i_res].n_ngh+1) : env.monte_flips;
    }
}

void monte2_set_toggle(PROT prot) {
    int i_res, i_conf;
    int n_fixed;
    
    for (i_res=0;i_res<prot.n_res;i_res++) {
        n_fixed = 0;
        for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
            if (prot.res[i_res].conf[i_conf].toggle == 'a') continue;
            if (prot.res[i_res].conf[i_conf].occ < env.monte_reduce) {
                prot.res[i_res].conf[i_conf].toggle = 'f';
                prot.res[i_res].conf[i_conf].occ = 0.;
                n_fixed++;
            }
        }
        /* if only one conformer is free, fix it too */
        if ((prot.res[i_res].n_conf - n_fixed) == 2) {
            for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                if (prot.res[i_res].conf[i_conf].toggle != 'f') {
                    prot.res[i_res].conf[i_conf].toggle = 'f';
                    prot.res[i_res].conf[i_conf].occ = 1. - prot.res[i_res].fixed_occ;
                    prot.res[i_res].fixed_occ = 1.;
                    n_fixed++;
                    break;
                }
            }
        }
    }
}

PROT monte2_reduce(PROT prot) {
    /* Uses toggle in prot to creat prot_red, reduced version of prot */
    PROT prot_red;
    int i_res, i_conf, j_res,j_conf;
    int ic, jc;
    
    memset(&prot_red,0,sizeof(PROT));
    cpy_prot(&prot_red, &prot);
    
    for (i_res=prot_red.n_res-1;i_res>=0;i_res--) {
        for (i_conf=prot_red.res[i_res].n_conf-1;i_conf>=1;i_conf--) {
            if (prot_red.res[i_res].conf[i_conf].toggle == 'f') {       /* turn off */
                ic = prot_red.res[i_res].conf[i_conf].i_conf_prot;
                
                prot_red.res[i_res].fixed_occ += prot_red.res[i_res].conf[i_conf].occ;
                
                /* Fix energy by mfe */
                prot_red.E_base += prot_red.res[i_res].conf[i_conf].E_self * prot_red.res[i_res].conf[i_conf].occ;
                
                for (j_res=0;j_res<prot_red.n_res;j_res++) {
                    if (i_res==j_res) continue;
                    for (j_conf=1;j_conf<prot_red.res[j_res].n_conf;j_conf++) {
                        jc = prot_red.res[j_res].conf[j_conf].i_conf_prot;
                        prot_red.res[j_res].conf[j_conf].E_self += my_pairwise[ic][jc] * prot_red.res[i_res].conf[i_conf].occ;
                       // f->E_self=f->E_self+f->E_TS-f->PE_TS;
                   //  f->PE_TS=f->E_TS;
                      //  prot_red.res[j_res].conf[j_conf].PE_TS= prot_red.res[j_res].conf[j_conf].E_TS;
                    }
                }
                
                /* delete this conformer */
                del_conf(&prot_red.res[i_res], i_conf);
                prot_red.nc--;
            }
        }
        /* checking number of free conformers, total occupancy,  */
        
        /* delete the residue if no free conformer left */
        if (prot_red.res[i_res].n_conf == 1) {
            del_res(&prot_red, i_res);
        }
    }
    
    prot_red.conf = (CONF**) malloc(prot_red.nc*sizeof(void *));
    ic = 0;
    for (i_res=0;i_res<prot_red.n_res;i_res++) {
        for (i_conf=1;i_conf<prot_red.res[i_res].n_conf;i_conf++) {
            prot_red.conf[ic] = &prot_red.res[i_res].conf[i_conf];
            ic++;
        }
    }
    return prot_red;
}

void monte2_mc(PROT *prot_p) {
    int k_res, k_conf, k_ngh;
    int n_flip, i_flip, j_flip, i_iter;
    int ic_old, ic_new, jc_w;
    int j_res;
    double E_delta;

    i_iter = env.monte_niter_cycle;
    while (i_iter) {
        
        /* select one residue */
        k_res = (int) (ran2(&idum) * prot_p->n_res);
        /* decide how many flips */
        n_flip = (int) (ran2(&idum) * prot_p->res[k_res].n_flip_max) + 1;
        /* get the list of residue to flip */
        flip_res[0] = &prot_p->res[k_res];
        for (i_flip=1;i_flip<n_flip;i_flip++) {
            j_flip = 0;
            while (j_flip < i_flip) {
                k_ngh = (int) (ran2(&idum) * prot_p->res[k_res].n_ngh);
                flip_res[i_flip] = prot_p->res[k_res].ngh[k_ngh];
                
                for (j_flip=1;j_flip<i_flip;j_flip++) {
                    if (flip_res[i_flip] == flip_res[j_flip]) break;
                }
            }
        }
        
        /* each residue pick new conformer */
        for (i_flip=0;i_flip<n_flip;i_flip++) {
            flip_res[i_flip]->conf_old = flip_res[i_flip]->conf_w;
            
            //flip_res[i_flip]->conf_new = NULL;
            //while (flip_res[i_flip]->conf_new == flip_res[i_flip]->conf_old) {
            //}
            k_conf = (int) (ran2(&idum) * (flip_res[i_flip]->n_conf - 1.)) + 1;
            flip_res[i_flip]->conf_new = &flip_res[i_flip]->conf[k_conf];
            flip_res[i_flip]->counter_trial++;
            flip_res[i_flip]->conf_new->counter_trial++;
        }
        
        /* calc. deltaE */
        E_delta = 0.;
        for (i_flip=0;i_flip<n_flip;i_flip++) {
            ic_old = flip_res[i_flip]->conf_old->i_conf_prot;
            ic_new = flip_res[i_flip]->conf_new->i_conf_prot;
            
            E_delta += (flip_res[i_flip]->conf_new->E_self - flip_res[i_flip]->conf_old->E_self);
            for (j_res=0; j_res<prot_p->n_res; j_res++) {
                jc_w = prot_p->res[j_res].conf_w->i_conf_prot;
                E_delta += (my_pairwise[ic_new][jc_w] - my_pairwise[ic_old][jc_w]);
            }
            
            /* real flip */
            flip_res[i_flip]->conf_w = flip_res[i_flip]->conf_new;
        }
        
        /* Metropolis */
        if (exp(-beta*E_delta) > ran2(&idum)) {
            prot_p->E_state += E_delta;
        }
        else {
            for (i_flip=0;i_flip<n_flip;i_flip++) {
                flip_res[i_flip]->conf_w = flip_res[i_flip]->conf_old;
            }
        }
        /*
        FILE *fp;
        fp = fopen(DETAIL,"a");
        for (k_res=0;k_res<prot_p->n_res;k_res++) {
            if (k_res==59||k_res==57)
            fprintf(fp,"%4d %3d %s",k_res,prot_p->res[k_res].conf_w->i_conf_prot,prot_p->res[k_res].conf_w->uniqID);
        }
        fprintf(fp,"%10.3f\n",prot_p->E_state);
        fclose(fp);
        */
        for (k_res=0;k_res<prot_p->n_res;k_res++) {
            prot_p->res[k_res].conf_w->counter_accept++;
        }
        if (prot_p->E_state < prot_p->E_min) prot_p->E_min = prot_p->E_state;
        prot_p->E_accum += prot_p->E_state;
        
        if (env.monte_do_energy) do_free_energy(prot_p);
        i_iter--;
    }
    
    /*
    for (k_res=0;k_res<prot_p->n_res;k_res++) {
        for (k_conf=1;k_conf<prot_p->res[k_res].n_conf;k_conf++) {
            printf("%s %d\n",prot_p->res[k_res].conf[k_conf].uniqID,prot_p->res[k_res].conf[k_conf].counter_accept);
        }
    }
    */
}

void zero_counters(PROT *prot_p) {
    int i_res,i_conf;
    prot_p->E_accum = 0.;
    for (i_res=0;i_res<prot_p->n_res;i_res++) {
        prot_p->res[i_res].counter_trial = 0;
        for (i_conf=1;i_conf<prot_p->res[i_res].n_conf;i_conf++) {
            prot_p->res[i_res].conf[i_conf].counter_trial = 0;
            prot_p->res[i_res].conf[i_conf].counter_accept = 0;
        }
    }
}

int monte2_check_toggle(PROT prot) {
    int i_res, i_conf, n_fixed;
    float fixed_occ;
    
    for (i_res=0;i_res<prot.n_res;i_res++) {
        fixed_occ = 0.;
        n_fixed = 0;
        for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
            if (prot.res[i_res].conf[i_conf].toggle == 'f') {
                n_fixed++;
                fixed_occ += prot.res[i_res].conf[i_conf].occ;
            }
        }
        if (fixed_occ > 1.) printf(" Warning! Total fixed occupancy on %s %c%04d is over 1.\n",
            prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq);fflush(stdout);
        
        if (prot.res[i_res].n_conf - n_fixed == 2) {
            for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                if (prot.res[i_res].conf[i_conf].toggle != 'f') {
                    prot.res[i_res].conf[i_conf].toggle = 'f';
                    prot.res[i_res].conf[i_conf].occ = 1.-fixed_occ;
                    if (prot.res[i_res].conf[i_conf].occ < 0.) prot.res[i_res].conf[i_conf].occ = 0.;
                    n_fixed++;
                    fixed_occ += prot.res[i_res].conf[i_conf].occ;
                    break;
                }
            }
        }
        
        /* if total fixed occupancy is 1 fixed all the other free conformers */
        prot.res[i_res].fixed_occ = fixed_occ;
        if (prot.res[i_res].fixed_occ > (1.- 1e-6)) {
            if ((prot.res[i_res].n_conf - n_fixed) > 1) {
                for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                    if (prot.res[i_res].conf[i_conf].toggle != 'f') {
                        prot.res[i_res].conf[i_conf].toggle = 'f';
                        prot.res[i_res].conf[i_conf].occ = 0.;
                        n_fixed++;
                    }
                }
            }
        }
        
        /* if all confomers are fixed check if total fixed occupany is 1. */
        if (prot.res[i_res].n_conf - n_fixed == 1) {
            if (fixed_occ < 1.) printf(" Warning! Total fixed occupancy on %s %c%04d is smaller than 1.\n",
                prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq);fflush(stdout);
        }
    }
    return 0;
}

#define SAME_STATE_THR 0.0001
void do_free_energy(PROT *prot_p) {
    int i_res;
    int sum_kc, sum_krkc;
    int search1, search2, search_point;
    int state_found;
    
    sum_kc=0;
    for (i_res=0;i_res<prot_p->n_res;i_res++) sum_kc += prot_p->res[i_res].conf_w->i_conf_prot;
    sum_krkc = 0;
    for (i_res=0;i_res<prot_p->n_res;i_res++) sum_krkc += prot_p->res[i_res].conf_w->i_conf_prot * (i_res+1);
    
    search1 = -1;
    search2 = prot_p->n_saved;
    state_found = 0;

    while (1) {
        if (search1 >= search2-1) break;
        
        search_point = search1 + (search2-search1)/2.;
        
        if (prot_p->saved_states[search_point].sum_kc < sum_kc) {
            search1 = search_point;
            continue;
        }
        else if (prot_p->saved_states[search_point].sum_kc > sum_kc) {
            search2 = search_point;
            continue;
        }
        
        if (prot_p->saved_states[search_point].sum_kc == sum_kc) {
            if (prot_p->saved_states[search_point].sum_krkc < sum_krkc) {
                search1 = search_point;
                continue;
            }
            else if (prot_p->saved_states[search_point].sum_krkc > sum_krkc) {
                search2 = search_point;
                continue;
            }
        }
        
        if (prot_p->saved_states[search_point].sum_krkc == sum_krkc) {
            if ( fabs(prot_p->saved_states[search_point].E - prot_p->E_state) > SAME_STATE_THR ) {
                if (prot_p->saved_states[search_point].E < prot_p->E_state) {
                    search1 = search_point;
                    continue;
                }
                if (prot_p->saved_states[search_point].E > prot_p->E_state) {
                    search2 = search_point;
                    continue;
                }
            }
            else {
                state_found = 1;
                break;
            }
        }
    }
    if (!state_found) {
        prot_p->n_saved++;
        prot_p->saved_states = (UNIQ_STATE*) realloc(prot_p->saved_states, prot_p->n_saved*sizeof(UNIQ_STATE));
        memmove(prot_p->saved_states+search2+1, prot_p->saved_states+search2, (prot_p->n_saved-search2-1)*sizeof(UNIQ_STATE));
        prot_p->saved_states[search2].sum_kc = sum_kc;
        prot_p->saved_states[search2].sum_krkc = sum_krkc;
        prot_p->saved_states[search2].E = prot_p->E_state;
    }
}

double free_unf(PROT prot) {
    /*  (Adopted from Marilyn's free_unf.f)
    .           CALCULATION OF REFERENCE FREE ENERGY
    .   G(unfolded)= -kT Sum(i)(# of ionizable) ln(1 +
    .               exp(-2.3gamma(i)(pH-pKsol(i))))
    .
    . An-Suei Yang and Barry Honig, J.Mol.Biol.(1993), 231, 459-474 
    .
    . Expand this equation so it applies to 
    . residues with multiple protonation/electron states
    .   G(unfolded) = -kT * 
    .     Sum(i_res)  { ln( Sum(k_state) { 
    .     exp(-2.3*gamma_proton(i_res,k_state)*(pH-pKa(i_res)) -(2.3/58.)*gamma_electron(i_res,k_state)*(Eh-Em(i_res))) }
    */
    
    float free_eng = 0., sum_expE;
    int   i_res,i_conf,j_conf;
    
    beta = KCAL2KT/(env.monte_temp/ROOMT);
    for (i_res=0;i_res<prot.n_res;i_res++) {
        sum_expE = 0.;
        for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
            for (j_conf=1;j_conf<i_conf;j_conf++) {
                if (fabs(prot.res[i_res].conf[i_conf].H - prot.res[i_res].conf[j_conf].H) < 1e-4 &&
                    fabs(prot.res[i_res].conf[i_conf].e - prot.res[i_res].conf[j_conf].e) < 1e-4) break;
            }
            if (j_conf < i_conf) continue;
            sum_expE += exp(-beta*(prot.res[i_res].conf[i_conf].E_ph+prot.res[i_res].conf[i_conf].E_eh));
            //printf("%s,%8.3f\n",prot.res[i_res].conf[i_conf].uniqID,prot.res[i_res].conf[i_conf].E_ph);
        }
        if (sum_expE < 1.e-8) continue;
        free_eng += - log(sum_expE) / beta;
    }
    return free_eng;
}

int my_Nx;       /* number of titration points */
float *xp2, *yp2;   /* titration points */

int monte_out(PROT prot, int n_titra) {
    FILE *fp;
    int  N_crg, N_res, j, Counter, i, ic, i_titra, i_res;
    int n_protons, n_electrons, n_crg;
    int n_protons_grnd, n_electrons_grnd, n_crg_grnd;
    //float H, e;
    char line[20];
    char **head2, **mhead;  // add mhead for mfe header
    float **ypp2, **ysp2;
    //int my_Nx;       /* number of titration points */
    float *xp2, *yp2;   /* titration points */
    //float    **occ_table;   /* occ of conformers at various pH/Eh */
    float *netcrg;
    float *new_netcrg;
    float *crg, *protons, *electrons;    /* total net charge at pHs */
    char sbuff[21];

    /* writing occ table */
    if (!(fp = fopen(OCC_TABLE, "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", OCC_TABLE);fflush(stdout);
        return USERERR;
    }
    if (env.titr_type == 'p') {
        fprintf(fp, " ph           ");
        for (i=0; i<n_titra; i++) fprintf(fp, " %5.1f", env.titr_ph0+i*env.titr_phd);
    }
    else {
        fprintf(fp, " eh           ");
        for (i=0; i<n_titra; i++) fprintf(fp, " %5.0f", env.titr_eh0+i*env.titr_ehd);
    }
    fprintf(fp, "\n");
    for (ic=0; ic<prot.nc; ic++) {
        fprintf(fp, "%s", prot.conf[ic]->uniqID);
        for (i_titra=0; i_titra<n_titra; i_titra++) {
            fprintf(fp, " %5.3f", prot.conf[ic]->occ_table[i_titra]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);


    
    my_Nx = env.titr_steps;
    

    head2 = (char **) malloc(prot.nc * sizeof(char *));
    for (i=0; i<prot.nc; i++) head2[i] = (char *) malloc(20 * sizeof(char));
    shead2 = (char **) malloc(prot.nc * sizeof(char *));
    for (i=0; i<prot.nc; i++) shead2[i] = (char *) malloc(20 * sizeof(char));
    mhead = (char **) malloc(prot.nc * sizeof(char *));
    for (i=0; i<prot.nc; i++) mhead[i] = (char *) malloc(20 * sizeof(char));
    netcrg = (float *) malloc(prot.nc*sizeof(float));
    new_netcrg = (float *) malloc(prot.nc*sizeof(float));

    /*--- Assign titration points to an array ---*/
    xp2 = (float *) malloc(my_Nx * sizeof(float));      /* store x points */
    yp2 = (float *) malloc(my_Nx * sizeof(float));      /* store y points */
    ypp2 = (float **) malloc(prot.nc * sizeof(float *));
    for (i=0; i<prot.nc; i++) ypp2[i] = (float *) malloc(my_Nx * sizeof(float));
    ysp2 = (float **) malloc(prot.nc * sizeof(float *));
    for (i=0; i<prot.nc; i++) ysp2[i] = (float *) malloc(my_Nx * sizeof(float));

    
    //<<< Group into residues >>>
    // keep only charged conformers 
    Counter = 0;
    
    // Loop through all the conf in prot
    for (i=0; i <prot.nc; i++) {
        if (strchr(prot.conf[i]->uniqID, '+') || strchr(prot.conf[i]->uniqID, '-')) {
            strncpy(head2[Counter], prot.conf[i]->uniqID, 5); // copies up to 5 characters from the string pointed to, by prot.conf[i]->uniqID to dest head[Counter]
            head2[Counter][5] = '\0';
            strncat(head2[Counter], prot.conf[i]->uniqID+5, 8); // appends the string pointed to by prot.conf[i]->uniqID+5 to the end of the string pointed to by head[Counter] up to 8 characters long
            head2[Counter][11] = '\0';
            //printf("uniqID %s has charge of %5.2f\n",head[Counter],prot.conf[i]->netcrg);
            new_netcrg[Counter] = prot.conf[i]->netcrg; // This is the same as of column 5 in head3.lst
            for (i_titra=0; i_titra<n_titra; i_titra++) {
                ypp2[Counter][i_titra] = prot.conf[i]->occ_table[i_titra];
            } 
            Counter++;
        }
    }
    N_crg = Counter;

    
    // group residues 
    // my_Nx is number of titration points
    
    Counter = 0;
    strncpy(line, head2[0], 11);
    netcrg[0] = new_netcrg[0];
    
    for (i_titra=0; i_titra<n_titra; i_titra++) {
        ysp2[0][i_titra] = ypp2[0][i_titra];
    }

    for (i=1; i<N_crg; i++) {
        if (strncmp(head2[i], line, 11)) {      // not equal, a new residue // compares at most the first 11 bytes of head[i] and line
            strncpy(shead2[Counter], line, 11); 
            shead2[Counter][11] = '\0';
            Counter++;
            strncpy(line, head2[i], 11);
            netcrg[Counter] = new_netcrg[i];
            for (i_titra=0; i_titra<n_titra; i_titra++) ysp2[Counter][i_titra] = ypp2[i][i_titra];
        }
        else {                                 // same residue
            for (i_titra=0; i_titra<n_titra; i_titra++) ysp2[Counter][i_titra] += ypp2[i][i_titra];
        }
    }

    strncpy(shead2[Counter], line, 11); 
    shead2[Counter][11] = '\0';
    N_res = Counter + 1;
    
    crg       = (float *) malloc(my_Nx * sizeof(float));
    protons   = (float *) malloc(my_Nx * sizeof(float));
    electrons = (float *) malloc(my_Nx * sizeof(float));

    memset(crg, 0, my_Nx * sizeof(float));
    memset(protons, 0, my_Nx * sizeof(float));
    memset(electrons, 0, my_Nx * sizeof(float));

    // writing sum_crg 
    if (!(fp = fopen(TOT_CRG, "w"))) {
        printf("   FATAL: Can not write sumcrg to file \"%s\"\n", OCC_TABLE);fflush(stdout);
        return USERERR;
    }
    if (env.titr_type == 'p') {
        fprintf(fp, " pH       ");
        for (i=0; i<n_titra; i++) fprintf(fp, " %5.1f", env.titr_ph0+i*env.titr_phd);
    }
    else {
        fprintf(fp, " Eh       ");
        for (i=0; i<n_titra; i++) fprintf(fp, " %5.1f", env.titr_eh0+i*env.titr_ehd);
    }
    fprintf(fp, "\n");
    
    for (i_res=0;i_res<N_res; i_res++) {
        fprintf(fp, "%s", shead2[i_res]);
        strncpy(sbuff, shead2[i_res], 4); sbuff[4] = '1'; sbuff[5] = '\0';
        if (param_get( "PROTON", sbuff, "", &n_protons)) n_protons = 0;
        if (param_get( "ELECTRON", sbuff, "", &n_electrons)) n_electrons = 0;
        /* n_crg = n_protons-n_electrons; */
        n_crg = n_protons-n_electrons; 
        
        /* number of protons on ground conformer type */
        strncpy(sbuff, shead2[i_res], 3); sbuff[3] = '0'; sbuff[4] = '1'; sbuff[5] = '\0';
        if (param_get( "PROTON", sbuff, "", &n_protons_grnd)) n_protons_grnd = 0;
        if (param_get( "ELECTRON", sbuff, "", &n_electrons_grnd)) n_electrons_grnd = 0;
        n_crg_grnd = n_protons_grnd - n_electrons_grnd;

        for (i_titra=0; i_titra<n_titra; i_titra++) {
            //fprintf(fp, " %5.2f", prot.res[i_res].sum_crg[i_titra]);
            //printf("three: %s %5.2f\n", shead2[i_res],new_netcrg[i_res]);
            fprintf(fp, " %5.2f", netcrg[i_res]*ysp2[i_res][i_titra]);
            crg[i_titra] += netcrg[i_res]*ysp2[i_res][i_titra];
            protons[i_titra] += n_protons*ysp2[i_res][i_titra] + n_protons_grnd*(1.0-ysp2[i_res][i_titra]);
            electrons[i_titra] += n_electrons*ysp2[i_res][i_titra] + n_electrons_grnd*(1.0-ysp2[i_res][i_titra]);
        }
        fprintf(fp, "\n");
    }
    
    
    fprintf(fp, "----------\n");
    fprintf(fp, "Tot Protn");
    for(j=0; j<my_Nx; j++) {
        fprintf(fp, " %5.2f", protons[j]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Tot Elec ");
    for(j=0; j<my_Nx; j++) {
        fprintf(fp, " %5.2f", electrons[j]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Net crg  ");
    for(j=0; j<my_Nx; j++) {
        fprintf(fp, " %5.2f", crg[j]);
    }
    fprintf(fp, "\n");
    
    fclose(fp);

    free(crg);
    free(protons);
    free(electrons);
    free(xp2);
    free(yp2);
    for (i=0; i<prot.nc; i++) free(ypp2[i]);
    free(ypp2);
    for (i=0; i<prot.nc; i++) free(ysp2[i]);
    free(ysp2);
    for (i=0; i<prot.nc; i++) free(head2[i]);
    free(head2);
    for (i=0; i<prot.nc; i++) free(shead2[i]);
    free(shead2);
    return 0;
}

typedef struct {
    float pK;
    float n;
    float chi2;
} PKA;

PKA fitting(float guessed_pK, float guessed_n);
float get_chi2(float v[]);
extern void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk);
float *pH, *sumcrg;
float sumcrg_start,sumcrg_end;

/*
int curve_fitting(PROT prot, int n_titra) {
    int i_res, i_conf;
    FILE *fp;
    PKA  pka;
    float guessed_pK, guessed_n;
    int i_titra;
    float H, e;
    pH = malloc(env.titr_steps * sizeof(float));

    /////////
    for (i_res=0;i_res<prot.n_res;i_res++) {
        // Check if there's more than one charge state in this residue
        if (env.monte_print_nonzero) {
            for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                if (fabs(prot.res[i_res].conf[i_conf].H) > 1e-4) break;
                if (fabs(prot.res[i_res].conf[i_conf].e) > 1e-4) break;
            }
            if (i_conf >= prot.res[i_res].n_conf) continue;
        }
        else {
            H = prot.res[i_res].conf[1].H;
            e = prot.res[i_res].conf[1].e;
            for (i_conf=2;i_conf<prot.res[i_res].n_conf;i_conf++) {
                if (fabs(prot.res[i_res].conf[i_conf].H - H) > 1e-4) break;
                if (fabs(prot.res[i_res].conf[i_conf].e - e) > 1e-4) break;
            }
            if (i_conf >= prot.res[i_res].n_conf) continue;
        }
        
        prot.res[i_res].sum_crg = malloc(n_titra*sizeof(float));
        memset(prot.res[i_res].sum_crg,0,n_titra*sizeof(float));
        for (i_titra=0; i_titra<n_titra; i_titra++) {
            for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                prot.res[i_res].sum_crg[i_titra] += prot.res[i_res].conf[i_conf].occ_table[i_titra] * prot.res[i_res].conf[i_conf].netcrg;
            }
        }
    }
    prot.H = malloc(n_titra*sizeof(float));
    memset(prot.H,0,n_titra*sizeof(float));
    prot.e = malloc(n_titra*sizeof(float));
    memset(prot.e,0,n_titra*sizeof(float));
    prot.sum_crg = malloc(n_titra*sizeof(float));
    memset(prot.sum_crg,0,n_titra*sizeof(float));
    for (i_titra=0; i_titra<n_titra; i_titra++) {
        for (i_res=0;i_res<prot.n_res;i_res++) {
            for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
                prot.H[i_titra] += prot.res[i_res].conf[i_conf].occ_table[i_titra] * prot.res[i_res].conf[i_conf].H;
                prot.e[i_titra] += prot.res[i_res].conf[i_conf].occ_table[i_titra] * prot.res[i_res].conf[i_conf].e;
                prot.sum_crg[i_titra] += prot.res[i_res].conf[i_conf].occ_table[i_titra] * prot.res[i_res].conf[i_conf].netcrg;
            }
        }
    }
    //////////
    if (!(fp = fopen(CURVE_FITTING, "w"))) {
        printf("   FATAL: can not write to file \"%s\".", CURVE_FITTING);
        return USERERR;
    }
    if (env.titr_type == 'p') {   // pH titration
        fprintf(fp, "Residue            pK   ");
    }
    else {      // Eh titration assumed 
        fprintf(fp, "Residue            Em   ");
    }
    fprintf(fp, "   n(slope)   1000*chi2\n");
    
    for (i_titra=0; i_titra<env.titr_steps; i_titra++) {
        if (env.titr_type == 'p') pH[i_titra] = env.titr_ph0 + ((float)i_titra)*env.titr_phd;
        else pH[i_titra] = env.titr_eh0 + ((float)i_titra)*env.titr_ehd;
    }
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        //printf("%s\n", prot.res[i_res].resName);
        //printf("%5.2f\n", prot.res[i_res].sum_crg);
        if (!prot.res[i_res].sum_crg) continue;
        sumcrg = prot.res[i_res].sum_crg;
        
        // a reasonable guess makes optimization easier
        guessed_n = 0;
        guessed_pK = (pH[0]+pH[env.titr_steps-1])/2.;
        
        if (prot.res[i_res].sum_crg[0] >= 0.) sumcrg_start = (int) (prot.res[i_res].sum_crg[0] + 0.5);
        else sumcrg_start = (int) (prot.res[i_res].sum_crg[0] - 0.5);
        if (prot.res[i_res].sum_crg[env.titr_steps-1] >= 0.) sumcrg_end   = (int) (prot.res[i_res].sum_crg[env.titr_steps-1] + 0.5);
        else sumcrg_end   = (int) (prot.res[i_res].sum_crg[env.titr_steps-1] - 0.5);

        if (fabs(sumcrg_start - sumcrg_end)<1e-6) {
            fprintf(fp, "%s %c%04d%c         pK or Em out of range\n", prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].iCode);
            continue;
        }
        if (env.titr_type == 'p') {
            if (fabs(env.titr_phd) < 1e-6) {
                fprintf(fp, "%s %c%04d%c         pK or Em out of range\n", prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].iCode);
                continue;
            }
        }
        else if (env.titr_type == 'e') {
            if (fabs(env.titr_ehd) < 1e-6) {
                fprintf(fp, "%s %c%04d%c         pK or Em out of range\n", prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].iCode);
                continue;
            }
        }
        //printf("%f,%f,%f,%f\n",sumcrg_start,sumcrg_end,guessed_pK, guessed_n);
        pka = fitting(guessed_pK, guessed_n);
        
        pka.n = pka.n * (env.monte_temp/ROOMT);
        fprintf(fp, "%s %c%04d%c    %9.3f %9.3f %9.3f\n",
        prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].iCode,pka.pK,pka.n,pka.chi2*1000.);
    }
    free(pH);
    
    return 0;
}*/

PKA fitting(float guessed_pK, float guessed_n) {
    /* initialize the simplex and set up optimization */
    float **param;
    float chi2[3];
    int   i, n_iter;
    PKA   pka;
    
    param = (float**) malloc(3 * sizeof(float *));
    for (i=0;i<3;i++) param[i] = (float*) malloc(2 * sizeof(float));

    param[0][0] = guessed_pK;      param[0][1] = guessed_n;      chi2[0] = get_chi2(param[0]);
    param[1][0] = guessed_pK+0.01; param[1][1] = guessed_n;      chi2[1] = get_chi2(param[1]);
    param[2][0] = guessed_pK;      param[2][1] = guessed_n+1.0;  chi2[2] = get_chi2(param[2]);

    //printf("%f,%f,%f\n",chi2[0],chi2[1],chi2[2]);
    dhill(param, chi2, 2, 0.00001, &get_chi2, &n_iter);
    
    /*  fit this eq. sumcrg = sumstart + (sumend - sumstart) * exp(n(pH-pK))/(1+exp(n(pH-pK)))  */
    
    pka.pK =   param[0][0];
    pka.n  =   param[0][1];
    pka.chi2 = chi2[0];
    
    for (i=0;i<3;i++)
        free(param[i]);
    free(param);
    
    return pka;
}

float get_chi2(float param[]){
    float chi2 = 0.0;
    float y_fit;
    float T;
    int i;
    
    for (i=0; i<env.titr_steps; i++) {
        if (env.titr_type == 'p') T = 2.303*param[1]*(pH[i]-param[0]);
        else T = 2.303*param[1]*(pH[i]-param[0])/58.;

        if (T<50.) y_fit = sumcrg_start + (sumcrg_end - sumcrg_start) * exp(T)/(1.0+exp(T));
        else y_fit = sumcrg_end;

        chi2 += (y_fit-sumcrg[i])*(y_fit-sumcrg[i]);
        //printf("%f,%f,%f,%f,%f,%f\n",pH[i],param[0],param[1],T,y_fit,chi2);
    }
    
    return chi2;
}

int gmonte_group(PROT *prot_p) {
    int n_conf;
    int i_group,i_res,i_subres,i_conf,i_state,n_state;
    
    /* Group */
    prot_p->n_group = 0;
    for (i_res=0;i_res<prot_p->n_res;i_res++) {
        if (prot_p->res[i_res].groupID == 0) {
            prot_p->n_group++;
            prot_p->group = (GROUP *) realloc(prot_p->group,prot_p->n_group*sizeof(GROUP));
            i_group = prot_p->n_group-1;
            memset(&prot_p->group[i_group],0,sizeof(GROUP));
        }
        else {
            for (i_group=prot_p->n_group-1;i_group>=0;i_group--) {
                if (prot_p->group[i_group].res[0]->groupID == prot_p->res[i_res].groupID) break;
            }
            if (i_group == -1) {
                prot_p->n_group++;
                prot_p->group = (GROUP *) realloc(prot_p->group,prot_p->n_group*sizeof(GROUP));
                i_group = prot_p->n_group-1;
                memset(&prot_p->group[i_group],0,sizeof(GROUP));
            }
        }
        
        prot_p->group[i_group].n_res++;
        prot_p->group[i_group].res = (RES **) realloc(prot_p->group[i_group].res,prot_p->group[i_group].n_res*sizeof(void *));
        prot_p->group[i_group].res[prot_p->group[i_group].n_res-1] = &prot_p->res[i_res];
    }
    
    /* Subres */
    for (i_res=0;i_res<prot_p->n_res;i_res++) {
        prot_p->res[i_res].n_subres = 0;
        for (i_conf=1;i_conf<prot_p->res[i_res].n_conf;i_conf++) {
            /* If subres ID is 0, group them with # of H and e */
            if (prot_p->res[i_res].conf[i_conf].subresID == 0) {
                for (i_subres = prot_p->res[i_res].n_subres-1;i_subres>=0;i_subres--) {
                    if (prot_p->res[i_res].subres[i_subres].conf[0]->subresID == 0) {
                        if (fabs(prot_p->res[i_res].conf[i_conf].H - prot_p->res[i_res].subres[i_subres].conf[0]->H) < 1e-4 &&
                            fabs(prot_p->res[i_res].conf[i_conf].e - prot_p->res[i_res].subres[i_subres].conf[0]->e) < 1e-4) break;
                    }
                }
                if (i_subres == -1) {
                    prot_p->res[i_res].n_subres++;
                    prot_p->res[i_res].subres = (SUBRES *) realloc(prot_p->res[i_res].subres,prot_p->res[i_res].n_subres*sizeof(SUBRES));
                    i_subres = prot_p->res[i_res].n_subres-1;
                    memset(&prot_p->res[i_res].subres[i_subres],0,sizeof(SUBRES));
                }
            }
            else {
                for (i_subres = prot_p->res[i_res].n_subres-1;i_subres>=0;i_subres--) {
                    if (prot_p->res[i_res].subres[i_subres].conf[0]->subresID 
                        == prot_p->res[i_res].conf[i_conf].subresID) break;
                }
                if (i_subres == -1) {
                    prot_p->res[i_res].n_subres++;
                    prot_p->res[i_res].subres = (SUBRES *) realloc(prot_p->res[i_res].subres,prot_p->res[i_res].n_subres*sizeof(SUBRES));
                    i_subres = prot_p->res[i_res].n_subres-1;
                    memset(&prot_p->res[i_res].subres[i_subres],0,sizeof(SUBRES));
                }
            }
            prot_p->res[i_res].subres[i_subres].n_conf++;
            n_conf = prot_p->res[i_res].subres[i_subres].n_conf;
            prot_p->res[i_res].subres[i_subres].conf = (CONF **) realloc(prot_p->res[i_res].subres[i_subres].conf,n_conf*sizeof(void *));
            prot_p->res[i_res].subres[i_subres].conf[n_conf-1] = &prot_p->res[i_res].conf[i_conf];
        }
    }
    
    /* State */
    for (i_group=0;i_group<prot_p->n_group;i_group++) {
        prot_p->group[i_group].n_state = 1;
        for (i_res=0;i_res<prot_p->group[i_group].n_res;i_res++) {
            n_state = n_state * prot_p->group[i_group].res[i_res]->n_subres;
        }
        
        for (i_state=0;i_state<n_state;i_state++) {
            
        }
    }
    return 0;
}

void gmonte_mc(PROT *prot_p) {
    int k_group, k_state, k_conf, k_ngh;
    int n_flip, i_flip, j_flip, i_iter;
    float E_delta;
    int ic_old, ic_new, jc;
    int i_res, j_res;
    
    i_iter = env.monte_niter_cycle;
    while (i_iter) {
        /* select one group */
        k_group = (int) (ran2(&idum) * prot_p->n_group);
        /* decide how many flips */
        n_flip = (int) (ran2(&idum) * prot_p->group[k_group].n_flip_max) + 1;
        /* get the list of groups to flip */
        flip_group[0] = &prot_p->group[k_group];
        for (i_flip=1;i_flip<n_flip;i_flip++) {
            j_flip = 0;
            while (j_flip < i_flip) {
                k_ngh = (int) (ran2(&idum) * prot_p->group[k_group].n_ngh);
                flip_group[i_flip] = prot_p->group[k_group].ngh[k_ngh];
                
                for (j_flip=1;j_flip<i_flip;j_flip++) {
                    if (flip_group[i_flip] == flip_group[j_flip]) break;
                }
            }
        }
        
        /* each group pick new state */
        for (i_flip=0;i_flip<n_flip;i_flip++) {
            flip_group[i_flip]->state_old = flip_group[i_flip]->state_w;
            k_state = (int) (ran2(&idum) * flip_group[i_flip]->n_state);
            flip_group[i_flip]->state_new = &flip_group[i_flip]->state[k_state];
            
            /* each residue pick new subres */
            for (i_res=0;i_res<flip_group[i_flip]->n_res;i_res++) {
                flip_group[i_flip]->res[i_res]->subres_new = flip_group[i_flip]->state[k_state].subres[i_res];
                flip_group[i_flip]->res[i_res]->conf_old = flip_group[i_flip]->res[i_res]->conf_w;
                
                /* pick new conf */
                k_conf = (int) (ran2(&idum) * flip_group[i_flip]->res[i_res]->subres_new->n_conf);
                flip_group[i_flip]->res[i_res]->conf_new = flip_group[i_flip]->res[i_res]->subres_new->conf[k_conf];
            }
        }
        
        /* calc. deltaE */
        E_delta = 0.;
        for (i_flip=0;i_flip<n_flip;i_flip++) {
            for (i_res=0;i_res<flip_group[i_flip]->n_res;i_res++) {
                ic_old = flip_group[i_flip]->res[i_res]->conf_old->i_conf_prot;
                ic_new = flip_group[i_flip]->res[i_res]->conf_new->i_conf_prot;
                for (j_res=0; j_res<prot_p->n_res; j_res++) {
                    jc = prot_p->res[j_res].conf_w->i_conf_prot;
                    E_delta = E_delta + (my_pairwise[ic_new][jc] - my_pairwise[ic_old][jc]);
                }
                /* real flip */
                flip_group[i_flip]->res[i_res]->conf_w = flip_group[i_flip]->res[i_res]->conf_new;
            }
        }
        
        /* Metropolis */
        if (exp(-beta*E_delta) > ran2(&idum)) {
            prot_p->E_state += E_delta;
        }
        else {
            for (i_flip=0;i_flip<n_flip;i_flip++) {
                for (i_res=0;i_res<flip_group[i_flip]->n_res;i_res++) {
                    flip_group[i_flip]->res[i_res]->conf_w = flip_group[i_flip]->res[i_res]->conf_old;
                }
            }
        }

        for (i_res=0;i_res<prot_p->n_res;i_res++) {
            prot_p->res[i_res].conf_w->counter_accept++;
        }
        if (prot_p->E_state < prot_p->E_min) prot_p->E_min = prot_p->E_state;
        prot_p->E_accum += prot_p->E_state;
        
        if (env.monte_do_energy) do_free_energy(prot_p);
        i_iter--;
    }
}

/*
int load1_pairwise()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, "",0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    
    
    
    if (!(my_pairwise = (double **) malloc(ematrix.n * sizeof(double *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    for (kc=0; kc<ematrix.n; kc++) {
        if (!(my_pairwise[kc] = (double *) malloc(ematrix.n * sizeof(double)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           
           my_pairwise[i][j] = my_pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele
                                             +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           
        }
    }
    free_ematrix(&ematrix);
    return 0;
}*/


/////////////////////////////////////////////////////////////////////////////////////
int AddToArray (ISTATES item,RES *s)
{
	if(s->num_of_states == s->memory_for_states) { // need more refs?
		if (s->num_of_states == 0)
			s->memory_for_states = 3; // start off with 3 refs
		else
			s->memory_for_states *= 2; // double the refs allocated

		s->IS = (ISTATES*)realloc(s->IS, (s->memory_for_states * sizeof(ISTATES)));
		if (s->IS == NULL)
		{
			fprintf(stdout, "ERROR: Couldn't realloc memory!");
			exit(2);
		}
	}
	//s->IS=TheArray;
	s->IS[s->num_of_states] = item; //adding the element to the array
	s->num_of_states++;

	return 0;
}

int exist(CONF con,RES *s)
{
int i;
	if (s->num_of_states==0) // no state exist yet so we need to start adding the states
		return 0;
if (con.uniqID[3] == 'D' && con.uniqID[4] == 'M') 
return 2;
	for(i=0;i<s->num_of_states;i++)
	{
		if(s->IS[i].Nelectrons==con.e && s->IS[i].Nprotons==con.H)//check if this confermer has the same electron
		{                                                               //and protons if it does so this state already exist
			s->IS[i].Nconfermers++;     // so here we increment the number of confermer for this state
			s->IS[i].totalOCC=s->IS[i].totalOCC+con.occ;
              			return 1;

		}
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void load_states(RES *r)// to load the ionization states 
{

	int i;
    int j=0;
    //float x=0.00;

    for(j=0;j<r->num_of_states;j++)
    {
    r->IS[j].totalOCC=0;
    }
    int ret;
    	for(i=1;i< r->n_conf;i++)//loop over all the confermer of this residue
    	{
    ret=exist(r->conf[i],r);
    		if(ret==0)// if this is a new ionizatio state add it to the array of the ionization state
    		{
    			ISTATES s;
    			s.totalOCC=0.0;
    			s.Nconfermers=0;
    			s.Nelectrons=r->conf[i].e;//load the number of electron
    			s.Nprotons=r->conf[i].H; //load the number of protons
    			s.Nconfermers++;//change the number of states for this confermer to one
    			s.totalOCC=r->conf[i].occ;
    			//Console::WriteLine("here");
    			AddToArray(s,r);//add the state to the array

    		}
                  if(ret==2)// if this is a new ionizatio state add it to the array of the ionization state
    		{
                      ISTATES s;
    			s.totalOCC=0.0;
                         s.Nelectrons=r->conf[i].e=-1;//load the number of electron
    			s.Nprotons=r->conf[i].H=-1; //load the number of protons
    			s.Nconfermers++;//change the number of states for this confermer to one
    			s.totalOCC=r->conf[i].occ;
    			//Console::WriteLine("here");
    			AddToArray(s,r);//add the state to the array
    //printf("%d",r->conf[i].H);
                  }
    		
    	}
    //for(j=0;j<r->num_of_states;j++)
    //{
    //r->IS[j].totalOCC=0;
    //printf("%f",r->IS[j].totalOCC);
    //printf("%d",r->IS[j].totalOCC);
    //}	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

void calculate_correstion(RES *s){
	double occ=0;
    int j;
    int i;
    for( j=0;j<s->num_of_states;j++)
    {
    	occ=0;
    	for(i=1;i<s->n_conf;i++){
        	if(s->conf[i].e==s->IS[j].Nelectrons&&s->conf[i].H==s->IS[j].Nprotons){
                if(s->conf[i].real_occ>.00000000001)
                occ-=(s->conf[i].real_occ)*log(s->conf[i].real_occ);// use the relation TS=-1.36*sum(Pi)log(Pi))and at this part we only calculate Pi log Pi
                //printf("%f",s->conf[i].real_occ);
            }
        }
    	s->IS[j].Cfactor=(occ/1.693);//multiplaying by 1.36 in order to complete the calculations
    	//Console::WriteLine("the correction factor for this state is    "+s.IS[j].Cfactor.ToString());
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_ETS(CONF *f,RES *s)//modifying the E_TS value for each conformer by adding the correction factor
{
    int i;
	//float correction=0;
	for( i=0;i<s->num_of_states;i++)//loop allover the states
	{
		if(s->IS[i].Nelectrons==f->e&&s->IS[i].Nprotons==f->H)//check if this confermer belong to this state
        {
        //  f->E_TS=0;
        //printf("%f",f->E_self);
			f->E_TS=s->IS[i].Cfactor;//if so add the correction factor for the energy term
        //printf("%f",f->E_TS);
            f->E_self=f->E_self+f->E_TS-f->PE_TS;
            f->PE_TS=f->E_TS;
        //printf("%f",f->PE_TS);
        }
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void load_real_occ(CONF *f,RES s)
{
int i;
	for( i=0;i<s.num_of_states;i++)//loop allover the states
	{
		if(s.IS[i].Nelectrons==f->e&&s.IS[i].Nprotons==f->H)
		{
			f->real_occ=f->occ/(s.IS[i].totalOCC);
//printf("%f",s.IS[i].totalOCC);
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void loop_residues(RES s)
{
int i;
	for( i=1;i<s.n_conf;i++)
	{
		load_real_occ(&s.conf[i],s);
	}
	//Console::WriteLine("here");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void load_correction_factor(RES *s)//add the correction factor to all the confermers in this residue
{
int i;
	for( i=1;i<s->n_conf;i++)//loop allover the confermers
	{

		get_ETS(&s->conf[i],s);//this function will check the ionization state for this confermers and add the appropiate correction to the energy

			}
//s->IS=NULL;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ENTROPY_CORRECTION_FUNCTION(PROT *p)
{
    //printf("start");
    int i;
    //int j;
    	for( i=0;i<p->n_res;i++)//loop all over the whole protein
    	{
    		//p->res[i].num_of_states=0;
    		//p->res[i].memory_for_states=0;
    		load_states(&p->res[i]);//load the ionization state for each residue
    		loop_residues(p->res[i]);
    		calculate_correstion(&p->res[i]);
    		load_correction_factor(&p->res[i]);//add the correction factor to each confermer
    		free(TheArray);
            //         for(j=0;j<p->res[i].n_conf;j++)
            //{
            //printf(p->res[i].conf[j].uniqID);
            //printf("  ");
            //printf("%f",p->res[i].conf[j].E_TS);
            //}
            TheArray=NULL;
    	}

    //	free(TheArray);	//free the memory
}

