#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <zlib.h>
#include "mcce.h"

int refresh_prot(PROT prot);

int make_matrices(PROT prot, char *dir)
{
	int i, kr, kc, counter, verbose;
	int n_conf, n_dummies;
	char fname[MAXCHAR_LINE];
	FILE *fp;
	char sbuff[MAXCHAR_LINE];
	char sbuff2[MAXCHAR_LINE];
	char confName[6];
	int natom;
	int serial;
	char uniqID[15];
	char dummy;
	EMATRIX ematrix;


	if (verbose != 1) {
		verbose = 0;
	}
	n_conf = 0;
	for (kr=0; kr<prot.n_res; kr++) {
		n_conf+=(prot.res[kr].n_conf-1);
	}

	/* Try to load an existing energy table */
	ematrix.n = 0; /* the load_energies program relies on n to determine the validation of pointers */
	if (load_energies(&ematrix, dir, verbose) == n_conf) {   /* partial run exits */
		printf("   File %s/%s exists, energy calculation will update this file\n", dir, ENERGY_TABLE);
	}
	else {  /* initialize a new one with current run */
		free_ematrix(&ematrix); /* if ematrix.n is 0, nothing would happen */
		ematrix.n = n_conf;
		if (!(ematrix.conf = (CONF_HEAD *) calloc(ematrix.n, sizeof(CONF_HEAD)))) {
			printf("   Allocate memory error\n");
			return USERERR;
		}
		if (!(ematrix.pw = (PAIRWISE **) calloc(ematrix.n, sizeof(PAIRWISE *)))) {
			printf("   Allocate memory error\n");
			return USERERR;
		}
		for (i=0; i<ematrix.n; i++) {
			if (!(ematrix.pw[i] = (PAIRWISE *) calloc(ematrix.n, sizeof(PAIRWISE)))) {
				printf("   Allocate memory error\n");
				return USERERR;
			}
		}
		for (i=0; i<ematrix.n; i++) {
			ematrix.conf[i].on = 'f';
		}
	}


	counter = 0; n_dummies = 0;
	for (kr=0; kr<prot.n_res; kr++) {
		for (kc=1; kc<prot.res[kr].n_conf; kc++) {
			strcpy(ematrix.conf[counter].uniqID, prot.res[kr].conf[kc].uniqID);
			strcpy(ematrix.conf[counter].history, prot.res[kr].conf[kc].history);
			strncpy(confName, ematrix.conf[counter].uniqID, 5); confName[5] = '\0';
			if (param_get("NATOM", confName, "", &natom)) {
				printf("   WARNING: no NATOM for %s, 0 assumed\n", confName);
				natom = 0;
				param_sav("NATOM", confName, "", &natom, sizeof(int));
			}
			if (natom == 0) { /* dummy */
				n_dummies ++;
				ematrix.conf[counter].on = 't';
				counter ++;
				continue;
			}
			ematrix.conf[counter].netcrg = prot.res[kr].conf[kc].netcrg;
			if (env.recalc_tors)
				ematrix.conf[counter].E_tors = torsion_conf(&prot.res[kr].conf[kc]);
			/* other parameters will be updated by head3lst_update once uniqID is defined*/
			counter++;
		} 
	}

	/* updating self energy terms for conformers in current delphi */
	refresh_prot(prot);
	counter = 0; n_dummies = 0;
	for (kr=0; kr<prot.n_res;kr++) {
		for (kc=1; kc<prot.res[kr].n_conf; kc++) {
			strncpy(confName, ematrix.conf[counter].uniqID, 5); confName[5] = '\0';
			if (param_get("NATOM", confName, "", &natom)) {
				printf("   WARNING: no NATOM for %s, 0 assumed\n", confName);
				natom = 0;
				param_sav("NATOM", confName, "", &natom, sizeof(int));
			}
			if (natom == 0) { /* dummy */
				n_dummies ++;
				ematrix.conf[counter].on = 't';
				counter ++;
				continue;
			}

			if (((counter-n_dummies+1) >= env.pbe_start && (counter-n_dummies+1) <= env.pbe_end)) {
				ematrix.conf[counter].E_vdw0 = prot.res[kr].conf[kc].E_vdw0;
				ematrix.conf[counter].E_vdw1 = prot.res[kr].conf[kc].E_vdw1;
				ematrix.conf[counter].E_tors = prot.res[kr].conf[kc].E_tors;
				ematrix.conf[counter].E_epol = prot.res[kr].conf[kc].E_epol;
				ematrix.conf[counter].E_rxn  = prot.res[kr].conf[kc].E_rxn;
				if (!strcmp(env.pbe_solver, "apbs") || !strcmp(env.rxn_method, "self")) ematrix.conf[counter].E_rxn0 = 0.0;
				ematrix.conf[counter].E_dsolv = ematrix.conf[counter].E_rxn - ematrix.conf[counter].E_rxn0;
			}
			counter++;
		}
	}
	head3lst_param(ematrix); /* fine tune and complete the rest param  */

	/* updating the pw of matrix */
	n_dummies = 0;
	for (kc=0; kc<ematrix.n; kc++) {
		/* dummy conformer doesn't have opp file */
		strncpy(confName, ematrix.conf[kc].uniqID, 5); confName[5] = '\0';
		if (param_get("NATOM", confName, "", &natom)) {
			printf("   WARNING: no NATOM for %s, 0 assumed\n", confName);
			natom = 0;
			param_sav("NATOM", confName, "", &natom, sizeof(int));
		}
		if (natom == 0) { /* dummy */
			n_dummies ++;
			ematrix.conf[kc].on = 't';
			for (i=0; i<ematrix.n; i++) {
				ematrix.pw[kc][i].ele = ematrix.pw[kc][i].vdw = ematrix.pw[kc][i].crt = ematrix.pw[kc][i].ori = 0.0;
				ematrix.pw[kc][i].mark[0] = '\0';
			}
			continue;
		}

		if (((kc-n_dummies+1) >= env.pbe_start && (kc-n_dummies+1) <= env.pbe_end)) { /* in current calculation */
			ematrix.conf[kc].on = 't';
			sprintf(fname, "%s.opp", ematrix.conf[kc].uniqID);
			if (!(fp = fopen(fname, "r"))) {
				printf("   FATAL: read file error %s\n", fname);
				return USERERR;
			}

			/* load pairwise energy */
			for (i=0; i<ematrix.n; i++) {
				/* Skip dummies */
				while (1) { /* a mismatch, dummy */
					strncpy(confName, ematrix.conf[i].uniqID, 5); confName[5] = '\0';
					if (param_get("NATOM", confName, "", &natom)) {
						printf("   FATAL: no NATOM for %s\n", confName);
						return USERERR;
					}

					if (natom == 0) {
						ematrix.pw[kc][i].ele = ematrix.pw[kc][i].vdw = ematrix.pw[kc][i].crt = ematrix.pw[kc][i].ori = 0.0;
						ematrix.pw[kc][i].mark[0] = '\0';
						dummy = 1;
						i++;  /* skip to the next valid */
					}
					else {
						dummy = 0;
						break;          /* break as non dummy */
					}

					if (i>=ematrix.n) break; /* break as done, dummy would continue */
				}
				if (dummy) break;

				/* everything reaches here is non dummy */
				if (!fgets(sbuff, sizeof(sbuff), fp)) {
					printf("   FATAL: Unexpected end of file %s.\n", fname);
					return USERERR;
				}


				strncpy(sbuff2, sbuff, 5); sbuff2[5] = '\0';
				serial = atoi(sbuff2);

				strncpy(uniqID, sbuff+6, 14); uniqID[14] = '\0';

				strncpy(sbuff2, sbuff+20, 10); sbuff2[10] = '\0';
				ematrix.pw[kc][i].crt = atof(sbuff2);

				strncpy(sbuff2, sbuff+30, 10); sbuff2[10] = '\0';
				ematrix.pw[kc][i].vdw = atof(sbuff2);

				strncpy(sbuff2, sbuff+40, 10); sbuff2[10] = '\0';
				ematrix.pw[kc][i].ori = atof(sbuff2);

				strncpy(ematrix.pw[kc][i].mark, sbuff+50, 3);
				*(strchr(ematrix.pw[kc][i].mark,'\n')) = '\0';

				//sscanf(sbuff, "%d %s %f %f %f %s", &serial, uniqID, &pairwise_raw[kc][i].ele, &pairwise_raw[kc][i].vdw, &pairwise_raw[kc][i].ori, pairwise_raw[kc][i].mark);
				if (strcmp(uniqID, ematrix.conf[i].uniqID)) {
					printf("   FATAL: Mismatch %s in protein structure and %s in %s.\n", ematrix.conf[i].uniqID,
							uniqID,
							fname);
					return USERERR;
				}
			}
			fclose(fp);
		}
	}


	if (write_energies(&ematrix, dir, verbose)) {
		printf("   Error in writing energy lookuptable %s/%s\n", dir, ENERGY_TABLE);
		return USERERR;
	}


	free_ematrix(&ematrix);

	return 0;
}

int load_energies(EMATRIX *ematrix, char *dir, int verbose)
/* this program returns number of conformers loaded, or -1 if no exsiting energy table */
{
	int i, n_conf;
    char sbuff[128];
    FILE *fp, *fp2;
    CONF_HEAD *conf;
    PAIRWISE *pw;
    char version[128];

    // Obtain the first line of the energy lookup table, which is the number of conformers 
    sprintf(sbuff, "%s/%s", dir, ENERGY_TABLE);
    if (!(fp = fopen(sbuff, "r"))) {
        printf("energies.opp not found\n");
        return -1;
    }

    // unzip the energies.opp
    fp2 = tmpfile();
    inf(fp, fp2);
    rewind(fp2);
    fclose(fp);


    fgets(sbuff, sizeof(sbuff), fp2);
    if (sscanf(sbuff, "%d %s", &n_conf, version) != 2) {

       printf("   Version mismatch: Opp file was made by pre-MCCE2.3 and this program is for %s\n", VERSION);
       printf("                     Use oppconvert to fix the mismatch\n");
       return USERERR;
    }

    if (strncmp(VERSION, version, 7)) {
    	// Do not check the version any more, xz 07-06-2013.
       //printf("   Version mismatch: Opp file was made by %s and this program is for %s\n", version, VERSION);
       //printf("                     Use oppconvert to fix the mismatch\n");
       return USERERR;
         
    }

    // allocate memeory 
    if (ematrix->n > 0) {  // existing table 
       if (ematrix->n != n_conf) {
          printf("error in loading Energy lookup table size %d on to %d\n", n_conf, ematrix->n);
          return USERERR;
       }
    }
    else {
       if (!(ematrix->conf = (CONF_HEAD *) calloc(n_conf, sizeof(CONF_HEAD)))) {
          printf("   Memory error in A load_energies\n");
          return 0; // none loaded and memory cleared 
       }
       if (!(ematrix->pw = (PAIRWISE **) calloc(n_conf, sizeof(PAIRWISE *)))) {
          printf("   Memory error in B load_energies\n");
          return 0; // none loaded and memory cleared 
       }
       for (i=0; i<n_conf; i++) {
          if (!(ematrix->pw[i] = (PAIRWISE *) calloc(n_conf, sizeof(PAIRWISE)))) {
             printf("   Memory error in C load_energies\n");
             return 0; // none loaded and memory cleared 
          }
       }
       ematrix->n = n_conf;
    }

    // Buffer of head part 
    if (!(conf = (CONF_HEAD *) calloc(n_conf, sizeof(CONF_HEAD)))){
          printf("   Memory error loading head in load_energies\n");
          return USERERR;
    }
    if (!(pw = (PAIRWISE *) calloc(n_conf, sizeof(PAIRWISE)))){
          printf("   Memory error loading PW buffer in load_energies\n");
          return USERERR;
    }


    if (fread(conf, sizeof(CONF_HEAD), ematrix->n, fp2) != n_conf) {
        printf("   Error in loading pairwise interaction headers at position %d\n", i);
        return USERERR;
    }
    for (i=0; i<ematrix->n; i++) {
        fread(pw, sizeof(PAIRWISE), ematrix->n, fp2);
		if ( verbose == 1 ) {
			printf("reading row %05d\r", i); fflush(stdout);
		}
        
        if (conf[i].on == 't') {
           memcpy(&ematrix->conf[i], &conf[i], sizeof(CONF_HEAD));
           memcpy(ematrix->pw[i], pw, ematrix->n*sizeof(PAIRWISE));
        }
        if (ematrix->conf[i].on != 't') // initilize even if not a valid delphi run 
           memcpy(&ematrix->conf[i], &conf[i], sizeof(CONF_HEAD));
    }
    if ( verbose == 1) {
		printf("\n"); fflush(stdout);
    }
    free(conf);
    free(pw);
    fclose(fp2);
    return n_conf;
}

int extract_matrix(EMATRIX *ematrix, char *dir, int verbose)
{
	int i, j; 
	FILE *fp;
	char sbuff[256];
	char fname[256];
	struct stat buf;

	if (stat(dir, &buf)) { /* no such file or dir */
      if (mkdir(dir, 0755)) {
         printf("   FATAL: Failed creating directory %s, no write permission.\n", dir);
         return USERERR;
      }
   }

   for (i=0; i<ematrix->n; i++) {
      sprintf(sbuff, "%s/%s.opp", dir, ematrix->conf[i].uniqID);
      if (!(fp=fopen(sbuff, "w"))) {
         printf("   Open file %s error\n", sbuff);
         return USERERR;
      }
      for (j=0; j<ematrix->n; j++) {
         /*if (fabs(ematrix->pw[i][j].ele) > PW_CUTOFF ||
             fabs(ematrix->pw[i][j].vdw) > PW_CUTOFF ||
             fabs(ematrix->pw[i][j].crt) > PW_CUTOFF ||
             fabs(ematrix->pw[i][j].ori) > PW_CUTOFF)*/
         fprintf(fp, "%05d %s %8.3f%8.3f%8.3f%8.3f %s\n", j+1, ematrix->conf[j].uniqID, ematrix->pw[i][j].ele, ematrix->pw[i][j].vdw, ematrix->pw[i][j].crt, ematrix->pw[i][j].ori, ematrix->pw[i][j].mark);
      }
      fclose(fp);
   }
	
	/*
	if ( verbose == 1 ) {
		printf(" Extracting matrix ...\n"); fflush(stdout);
	}

	for (i=0; i<ematrix->n; i++) {
		//printf("%11s\n",ematrix->conf[i].history);
		if (ematrix->conf[i].on != 't') continue;
		if (ematrix->conf[i].uniqID[3] == 'D' && ematrix->conf[i].uniqID[4] == 'M') continue;
		sprintf(sbuff, "%s/%s.opp", dir, ematrix->conf[i].uniqID);
		if (!(fp=fopen(sbuff, "w"))) {
			printf("   Open file %s error\n", sbuff);
			return USERERR;
		}
		for (j=0; j<ematrix->n; j++) {
			fprintf(fp, "%05d %s %8.3f%8.3f%8.3f%8.3f %s\n", j+1, ematrix->conf[j].uniqID, ematrix->pw[i][j].ele, ematrix->pw[i][j].vdw, ematrix->pw[i][j].crt, ematrix->pw[i][j].ori, ematrix->pw[i][j].mark);
		}
		fclose(fp);
	}

	// write head3.lst 
	sprintf(fname, "%s/%s", dir, FN_CONFLIST3);
	if (!(fp = fopen(fname, "w"))) {
		printf("   Can not open file %s to write. Abort ...\n", fname);
		return USERERR;
	}

	fprintf(fp, "iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    history\n");
	for (i=0; i<ematrix->n; i++) {
		if (ematrix->conf[i].on != 't') continue;
		fprintf(fp, "%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n",
				i+1,
				ematrix->conf[i].uniqID,
				'f', 0.00,
				ematrix->conf[i].netcrg,
				ematrix->conf[i].Em,
				ematrix->conf[i].pKa,
				ematrix->conf[i].e,
				ematrix->conf[i].H,
				ematrix->conf[i].E_vdw0,
				ematrix->conf[i].E_vdw1,
				ematrix->conf[i].E_tors,
				ematrix->conf[i].E_epol,
				ematrix->conf[i].E_dsolv,
				ematrix->conf[i].E_extra,
				ematrix->conf[i].history);
	}

	fclose(fp);*/


	return 0;
}


int head3lst_param(EMATRIX ematrix)
{  int kc;
   float rxn0, extra, em0, pka0;
   char sbuff[6];
   int   H, e;

   for (kc=0; kc<ematrix.n; kc++) {
         strncpy(sbuff, ematrix.conf[kc].uniqID, 5); sbuff[5] = '\0';

         if (param_get("EM", sbuff, "", &em0)) {
            printf("   WARNING: No EM entry for %s, set to 0.\n", sbuff);
            fflush(stdout);
            em0 = 0.0;
            param_sav("EM", sbuff, "", &em0, sizeof(float));
         }
         ematrix.conf[kc].Em = em0;

         if (param_get("PKA", sbuff, "", &pka0)) {
            printf("   WARNING: No PKA entry for %s, set to 0.\n",sbuff);
            fflush(stdout);
            pka0 = 0.0;
            param_sav("PKA",sbuff, "", &pka0, sizeof(float));
         }
         ematrix.conf[kc].pKa = pka0;

         if (param_get("PROTON",sbuff, "", &H)) {
            printf("   WARNING: no PROTON for %s, 0 assumed\n",sbuff);
            H = 0;
            param_sav("PROTON",sbuff, "", &H, sizeof(int));
         }
         ematrix.conf[kc].H = H;

         if (param_get("ELECTRON",sbuff, "", &e)) {
            printf("   WARNING: no ELECTRON for %s, 0 assumed\n",sbuff);
            e = 0;
            param_sav("ELECTRON",sbuff, "", &e, sizeof(int));
         }
         ematrix.conf[kc].e = e;

	 
	 if (!strcmp(env.pbe_solver, "delphi") && !strcmp(env.rxn_method, "surface")) {
	    if (param_get("RXN",sbuff, "", &rxn0)) {
	       printf("   WARNING: No RXN entry for %s, set to 0.\n",sbuff);
	       fflush(stdout);
	       rxn0 = 0.0;
	       param_sav("RXN",sbuff, "", &rxn0, sizeof(float));
	    }
	    ematrix.conf[kc].E_rxn0 = rxn0; /*this is the case with delphi*/
	 }
	 else ematrix.conf[kc].E_rxn0 = 0.0; /*this is the case with apbs*/

         if (param_get("EXTRA",sbuff, "", &extra)) {
            /* printf("   WARNING: No EXTRA entry for %s, set to 0.\n",sbuff); */
            fflush(stdout);
            extra = 0.0;
            param_sav("EXTRA",sbuff, "", &extra, sizeof(float));
         }
         ematrix.conf[kc].E_extra = extra;
		 
		 /*if we calculate rxn with self energies, then we don't need rxn0*/
		 ematrix.conf[kc].E_dsolv = ematrix.conf[kc].E_rxn - ematrix.conf[kc].E_rxn0;
         

         /* label unrealistic values */
         if (ematrix.conf[kc].E_vdw0 > 999.0) ematrix.conf[kc].E_vdw0 = 999.0;
         if (ematrix.conf[kc].E_vdw1 > 999.0) ematrix.conf[kc].E_vdw1 = 999.0;
   }

   return 0;
}


int write_energies(EMATRIX *ematrix, char *dir, int verbose)
{
	FILE *fp, *fp2;
    int i, j;
    char fname[MAXCHAR_LINE];
    

    
    // correction on pairwise interaction 
    for (i=0; i<ematrix->n; i++) {
        char conf1_type, conf2_type;
        if (ematrix->conf[i].uniqID[3] == '+' || ematrix->conf[i].uniqID[3] == '-') conf1_type = 'C'; /* one of "C" or "D" */
        else conf1_type = 'D';

        for (j=i; j<ematrix->n; j++) {
            if (ematrix->conf[j].uniqID[3] == '+' || ematrix->conf[j].uniqID[3] == '-') conf2_type = 'C'; /* one of "C" or "D" */
            else conf2_type = 'D';

            /* in a residue ? */
            if(strncmp(ematrix->conf[i].uniqID,ematrix->conf[j].uniqID, 3) ||
            strncmp(ematrix->conf[i].uniqID+5,ematrix->conf[j].uniqID+5, 6)) { /* not in the same residue */
                if (conf1_type == 'C' && conf2_type == 'C') { /* crg-crg corrected average */
                    /* corrected min */
                    if ((ematrix->pw[i][j].crt * ematrix->pw[i][j].ori < 0 || ematrix->pw[j][i].crt * ematrix->pw[j][i].ori < 0)
                      ||((strchr(ematrix->pw[i][j].mark, '?') && strchr(ematrix->pw[j][i].mark, '?')))) { /* abnormal case, sign flipped or both ? marked */
                       if (fabs(ematrix->pw[i][j].ori)>fabs(ematrix->pw[j][i].ori))
                          ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[j][i].ori/1.5;
                       else
                          ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[i][j].ori/1.5;
                    }
                    else if (strchr(ematrix->pw[i][j].mark, '?') && !strchr(ematrix->pw[j][i].mark, '?')) {
                       ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[j][i].crt;
                    }
                    else if (!strchr(ematrix->pw[i][j].mark, '?') && strchr(ematrix->pw[j][i].mark, '?')) {
                       ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[i][j].crt;
                    }
                    else {
                       if (fabs(ematrix->pw[i][j].crt)>fabs(ematrix->pw[j][i].crt))
                          ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[j][i].crt;
                       else
                          ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[i][j].crt;
                    }
                }

                else if (conf1_type == 'D' && conf2_type == 'D') { /* dip-dip corrected average*/
                    ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = (ematrix->pw[i][j].ori+ematrix->pw[j][i].ori)/2.0;
                }
                else { /* crg-dip */
                    /* when sign flipped or both "?" marked, use uncorrected value divided by 1.5 */
                    if ((ematrix->pw[i][j].crt * ematrix->pw[i][j].ori < 0 || ematrix->pw[j][i].crt * ematrix->pw[j][i].ori < 0)
                      ||((strchr(ematrix->pw[i][j].mark, '?') && strchr(ematrix->pw[j][i].mark, '?')))) {
                       ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = (ematrix->pw[i][j].ori+ematrix->pw[j][i].ori)/2.0/1.5;
                    }
                    else if (strchr(ematrix->pw[i][j].mark, '?') && !strchr(ematrix->pw[j][i].mark, '?')) {
                       ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[j][i].crt;
                    }
                    else if (!strchr(ematrix->pw[i][j].mark, '?') && strchr(ematrix->pw[j][i].mark, '?')) {
                       ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = ematrix->pw[i][j].crt;
                    }
                    else {
                       ematrix->pw[i][j].ele = ematrix->pw[j][i].ele = (ematrix->pw[i][j].crt+ematrix->pw[j][i].crt)/2.0;
                    }
                }
            }
            else { /* within a residue */
                ematrix->pw[i][j].ele = ematrix->pw[i][j].vdw = ematrix->pw[j][i].ele = ematrix->pw[j][i].vdw = 0.0;
            }
        }
    }
    
    
    
    // write out the matrices 
    fp2 = tmpfile();
    // The firest line is a two-field record, number of comformers and version number separated by white space 
    fprintf(fp2, "%d %s\n", ematrix->n, VERSION);
    fwrite(ematrix->conf, sizeof(CONF_HEAD), ematrix->n, fp2);
    for (i=0; i<ematrix->n; i++) {
	    if (verbose == 1) {
			printf("writing row %05d\r", i);
		}
        fwrite(ematrix->pw[i], sizeof(PAIRWISE), ematrix->n, fp2);
    }
    fputc(EOF, fp2); rewind(fp2);
	
    sprintf(fname, "%s/%s", dir, ENERGY_TABLE);
    if (!(fp = fopen(fname, "w"))) {
        printf("   Can not open file %s to write. Abort ...\n", fname);
        return USERERR;
    }
	
    if (def(fp2, fp, 9) != Z_OK) {
        printf("Compress file %s error\n", fname);
        fclose(fp);
        fclose(fp2);
        return USERERR;
    }
    fclose(fp2); fclose(fp);
    
    sprintf(fname, "%s/%s", dir, PW_DIR);
    if (extract_matrix(ematrix, fname,0)) return USERERR;

    // write head3.lst
    
    sprintf(fname, "%s/%s", dir, FN_CONFLIST3);
    if (!(fp = fopen(fname, "w"))) {
        printf("   Can not open file %s to write. Abort ...\n", fname);
        return USERERR;
    }

   fprintf(fp, "iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    history\n");
   for (i=0; i<ematrix->n; i++) {
         fprintf(fp, "%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10s%c\n",
                                                                    i+1,
                                                                    ematrix->conf[i].uniqID,
                                                                    'f', 0.00,
                                                                    ematrix->conf[i].netcrg,
                                                                    ematrix->conf[i].Em,
                                                                    ematrix->conf[i].pKa,
                                                                    ematrix->conf[i].e,
                                                                    ematrix->conf[i].H,
                                                                    ematrix->conf[i].E_vdw0,
                                                                    ematrix->conf[i].E_vdw1,
                                                                    ematrix->conf[i].E_tors,
                                                                    ematrix->conf[i].E_epol,
                                                                    ematrix->conf[i].E_dsolv,
                                                                    ematrix->conf[i].E_extra,
                                                                    ematrix->conf[i].history,
                                                                    ematrix->conf[i].on);
   }

   fclose(fp);
   return 0;
}

int free_ematrix(EMATRIX *ematrix)
{  int i;

   for (i=0; i<ematrix->n; i++) {
      free(ematrix->pw[i]);
   }
   if (ematrix->n > 0) {
      free(ematrix->pw);
      free(ematrix->conf);
   }
   ematrix->n = 0;
   
   return 0;
}

int refresh_prot(PROT prot)
{  int kr, kc;
   char fname[256];
   FILE *fp;
   int success = 0;
   char sbuff[MAXCHAR_LINE];
   char *sptr;

   /* update vdw0, vdw1, tors, epol, rxn */
   for (kr=0; kr<prot.n_res; kr++) {
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
         sprintf(fname, "%s.opp", prot.res[kr].conf[kc].uniqID);
         prot.res[kr].conf[kc].E_vdw0 = 0.0;
         prot.res[kr].conf[kc].E_vdw1 = 0.0;
         prot.res[kr].conf[kc].E_epol = 0.0;
         prot.res[kr].conf[kc].E_tors = 0.0;
         prot.res[kr].conf[kc].E_rxn  = 0.0;
         if (!(fp=fopen(fname, "r"))) { /* not done yet */
            success = -1;
         }
         else {
            while (fgets(sbuff, sizeof(sbuff), fp)) {
               sptr = strtok(sbuff, " \n");
               if (sptr == NULL) continue;
               if (!strcmp(sptr, "VDW_SELF")) {
                  sptr = strtok(NULL, " \n");
                  if (sptr) {
                     prot.res[kr].conf[kc].E_vdw0 = atof(sptr);
                  }
                  else prot.res[kr].conf[kc].E_vdw0 = 0.0;
               }
               else if (!strcmp(sptr, "VDW_BKBN")) {
                  sptr = strtok(NULL, " \n");
                  if (sptr) prot.res[kr].conf[kc].E_vdw1 = atof(sptr);
                  else prot.res[kr].conf[kc].E_vdw1 = 0.0;
               }
               else if (!strcmp(sptr, "ELE_BKBN")) {
                  sptr = strtok(NULL, " \n");
                  if (sptr) prot.res[kr].conf[kc].E_epol = atof(sptr);
                  else prot.res[kr].conf[kc].E_epol = 0.0;
               }
               else if (!strcmp(sptr, "TORSION")) {
                  sptr = strtok(NULL, " \n");
                  if (sptr) prot.res[kr].conf[kc].E_tors = atof(sptr);
                  else prot.res[kr].conf[kc].E_tors = 0.0;
               }
               else if (!strcmp(sptr, "RXNSINGLE")) {
                   int del_runs, i_run;
                   float rxn[100]; /*rxn at focusing runs*/
		   
		   /*delphi*/
                   del_runs = 0;
		   if (!strcmp(env.pbe_solver, "delphi") && !strcmp(env.rxn_method, "surface")) {
		      while ((sptr=strtok(NULL, " \n"))) {
			  del_runs++;
			  rxn[del_runs-1] = atof(sptr);
		      }
		      prot.res[kr].conf[kc].E_rxn = 0.0;
		      if (del_runs < 3) i_run=0;
		      else i_run = del_runs-3;
		      for (; i_run<del_runs; i_run++) {
			  if (prot.res[kr].conf[kc].E_rxn > rxn[i_run])
			      prot.res[kr].conf[kc].E_rxn = rxn[i_run]; /* report the most negative in the last 3 runs*/
		      }
		   }
		   /*apbs or delphi self energies*/
		   else {
		      sptr = strtok(NULL, " \n");
		      if (sptr) prot.res[kr].conf[kc].E_rxn = atof(sptr);
		      else prot.res[kr].conf[kc].E_rxn = 0.0;
		   }
               }
            }
            fclose(fp);
         }
      }
   }
   return success;
}
