#include <glib.h>
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

//set to skip monte carlo, only do MFE (needs pre-existing fort.38)
int MFE_ONLY = 0;

#define mev2Kcal 0.0235  // to keep consistent with mfe.py, use this constant instead of PH2KCAL/58, just for mfe part
#define TINY 1.0E-10
#define NMAX 5000
#define GET_PSUM for (j=0; j<ndim; j++) {\
                        for (sum=0.0, i=0; i<mpts; i++) sum += p[i][j];\
                         psum[j] = sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define NUNIQS 1000000




/* normal pw is garanteed to be smaller than 2000. When it is is bigger than 5000, it is 
 * represents a value refenced to a clashed pair of conformers.
 */

/* self energy in post_run_prot, prot_free, prot_fixed and copied to post_run_conflist
 * before sampling
 * post_run_pairwise in post_run_pairwise
 * counter in post_run_conflist
 */
/* This Monte Carlo Program adds entropy correction to each conformer type different by
 * proton and electron, so that the ionized and neutral will get equal chance regardless of
 * the population of conformers.
 * It has been found that the drift of calculated pKa with rotamer number is caused by this
 * entropy effect
 */
typedef struct {
    float vdw0;
    float vdw1;
    float tors;
    float ebkb;
    float dsol;
    float offset;
    float epws;
    float vpws;
    float pHpK0;
    float EhEm0;
    float TS;
    float residues;
    float total;
    float *mfePair;
    float *crg;
    float *vdw;
    float *ele;
} MFE;  // all the mfe energy terms

typedef struct  {
    int n;      /* Number of free confs on this residue */
    int on;
    int *conf;  /* idices of free confs */
    MFE mfe;   // let residue include mfe
} RESIDUE;

struct STAT {
    float a;
    float b;
    float chi2;
};

typedef struct  {
    int n;
    int *res;
} BIGLIST;

typedef struct {
   //int header;
   float x;
   float y;
   float z;
   float rad;
   float crg;
   //int trailer;
} UNF;


struct OccInRes {
    char OccInRes_chainID;
    int  OccInRes_resSeq;
    int  OccInRes_iConf;
    char *key_resSeq;
    float most_occupancy;
    //struct OccInRes *next;
};

// create a Linked List
struct node {
    int key_resSeq;
    float the_occupancy;
    char name[100];
    char node_chainID;
    int  node_resSeq;
    int  node_iConf;
    struct node *next;
};


struct hash {
    struct node *head;
    int count;
};

/* public variables */
PROT     post_run_prot;
RES      post_run_conflist;
float    **post_run_pairwise, **post_run_pairwise_vdw, **post_run_pairwise_ele;
float    post_run_E_base, post_run_E_state;
float    post_run_ph, post_run_eh;
int post_run_Nx;       /* number of titration points */
float *post_run_xp, *post_run_yp;   /* titration points */
RESIDUE  *post_run_free_res, *post_run_fixed_res, *post_run_all_res, *post_run_mfe_res; // post_run_mfe_res is all the mfe residue
BIGLIST  *post_run_biglist;   /* same length as free residues */
int      post_run_n_free, post_run_n_fixed, post_run_n_all, post_run_n_mfe;
FILE     *post_run_fp;
float    **post_run_occ_table;   
float    **most_occ_table; 
char **post_run_shead;
struct hash *hashTable = NULL;
int eleCount = 0;





float  postrun_fround3(float x);
int    postrun_print_mfe(int i_res, float mfeP, FILE *pK_fp,  FILE *res_fp);
int    postrun_get_mfe(int i_res, int t_point);     
int    postrun_load_pairwise();
int    postrun_load_pairwise_fround3();
int    postrun_load_pairwise_vdw();
int    postrun_load_occupancy();
void   postrun_group_confs();
int    postrun_fitit();
int    postrun_load_conflist();
int    potential_map();
int    write_amber_siz(FILE *siz); 
struct STAT postrun_fit(float a, float b);
float  postrun_score(float v[]);
void   postrun_dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk);
int    most_occ();
struct OccInRes searchInHash(int key_resSeq);


int postrun()
{   
    int i;
    /* Load conformer list from FN_CONFLIST3 */
    printf("   Load conformer list from file \"%s\" ...\n", FN_CONFLIST3);
    fflush(stdout);
    post_run_prot = new_prot();
    postrun_load_conflist();
    if (post_run_conflist.n_conf == 0) {
        printf("   FATAL: error in reading conformer list %s", FN_CONFLIST3);
        return USERERR;
    }

    
    /* load post_run_pairwise */
    printf("   Load post_run_pairwise interactions ...\n"); fflush(stdout);
    if (postrun_load_pairwise()) {
        printf("   FATAL: post_run_pairwise interaction not loaded\n");
        return USERERR;
    }
    printf("   Done\n\n");
    fflush(stdout);

    

    post_run_occ_table = (float **) malloc(post_run_conflist.n_conf*sizeof(float *));
    for (i=0; i<post_run_conflist.n_conf; i++) post_run_occ_table[i] = (float *) malloc(env.titr_steps * sizeof(float));

    most_occ_table = (float **) malloc(post_run_conflist.n_conf*sizeof(float *));
    for (i=0; i<post_run_conflist.n_conf; i++) most_occ_table[i] = (float *) malloc(env.titr_steps * sizeof(float));

    post_run_all_res = NULL;
    post_run_free_res = NULL;
    post_run_fixed_res = NULL;
    post_run_biglist  = NULL;
    

    for (i=0; i<env.titr_steps; i++) {
        postrun_group_confs();
    }

    printf("   Loading occupancy from fort.38 ...\n");
    if (postrun_load_occupancy()) {
        printf("   Fatal error detected in postrun_load_occupancy().\n");
        return USERERR;
    }
    printf("   Done\n\n"); fflush(stdout);
    
    
    
    /* curve fitting */
    printf("   Fit titration curves to get pKa/Em ...\n");
    fflush(stdout);
    if (postrun_fitit()) {
        printf("   Fatal error detected in fitting program.\n");
        return USERERR;
    }
    printf("   Done\n\n"); fflush(stdout);

   
    /* free post_run_pairwise */
    for (i=0; i<post_run_conflist.n_conf; i++)
        free(post_run_pairwise[i]);
    free(post_run_pairwise);
    
    //printf("   Creating connectivity table...\n"); fflush(stdout);
    //get_connect12(post_run_prot);
    //printf("   Done\n\n"); fflush(stdout);


    if (env.display_potential_map){
        if (potential_map()) {
            printf("   Fatal error detected in potential_map().\n");
            return USERERR;
        }
    }
    printf("   Output files:\n");
    printf("      %-16s: pKa or Em from titration curve fitting.\n", CURVE_FITTING);
    printf("      %-16s: Summary of residue charges.\n", "sum_crg.out");
    if (env.display_potential_map){
        printf("      %-16s: PDB for the most occupied conformer for each residue.\n", "step5_out.pdb");
        printf("      %-16s: DelPhi atom radii.\n", "fort.11");
        printf("      %-16s: DelPhi atom charges.\n", "fort.12");
        printf("      %-16s: DelPhi parameters.\n", "fort.10");
        printf("      %-16s: DelPhi potential map.\n", "phimap01.cube");
        printf("      %-16s: DelPhi dielectric map.\n", "epsmap01.eps");
    }
    
    return 0;
}

int postrun_load_conflist()
{   FILE *post_run_fp;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    CONF conf_temp;
    char notfound;
    int iconf, ires;
    int kr;
    int counter;
    
    post_run_conflist.n_conf = 0;
    post_run_conflist.conf   = NULL;
       
    if (!(post_run_fp=fopen(FN_CONFLIST3, "r"))) {
        printf("   FATAL: Can't open file %s\n", FN_CONFLIST3);
        return USERERR;
    }
    fgets(sbuff, sizeof(sbuff), post_run_fp); /* skip the first line */
    counter = 0;
    while(fgets(sbuff, sizeof(sbuff), post_run_fp)) {
        /* load this line to a conf template */
        if (strlen(sbuff) < 20) continue;
        sscanf(sbuff, "%d %s %c %f %f %f %f %d %d %f %f %f %f %f %f %s", &conf_temp.iConf,
        conf_temp.uniqID,
        &conf_temp.on,
        &conf_temp.occ,
        &conf_temp.netcrg,
        &conf_temp.Em,
        &conf_temp.pKa,
        &conf_temp.e,
        &conf_temp.H,
        &conf_temp.E_vdw0,
        &conf_temp.E_vdw1,
        &conf_temp.E_tors,
        &conf_temp.E_epol,
        &conf_temp.E_dsolv,
        &conf_temp.E_extra,
        conf_temp.history);

        conf_temp.E_TS = 0.0; /* initialize entropy effect at the time of loading post_run_conflist */

        /* rescale */
        conf_temp.E_vdw0 *= env.scale_vdw0;
        conf_temp.E_vdw1 *= env.scale_vdw1;
        conf_temp.E_epol *= env.scale_ele;
        conf_temp.E_tors *= env.scale_tor;
        conf_temp.E_dsolv*= env.scale_dsolv;

        strncpy(conf_temp.resName, conf_temp.uniqID, 3); conf_temp.resName[3] = '\0';
        strncpy(conf_temp.confName, conf_temp.uniqID, 5); conf_temp.confName[5] = '\0';
        conf_temp.chainID = conf_temp.uniqID[5];
        strncpy(stemp, conf_temp.uniqID+6, 4); stemp[4] = '\0';
        conf_temp.resSeq = atoi(stemp);
        conf_temp.iCode = conf_temp.uniqID[10];
        conf_temp.n_atom = 0;
        if (conf_temp.on == 't' || conf_temp.on == 'T') conf_temp.on = 't';
        else conf_temp.on = 'f';
        conf_temp.iConf = counter;
        /* creating post_run_conflist */
        iconf = ins_conf(&post_run_conflist, post_run_conflist.n_conf, 0);
        cpy_conf(&post_run_conflist.conf[iconf], &conf_temp);
        counter++;
        
        notfound = 1;
        for (kr=0; kr<post_run_prot.n_res; kr++) {
            if (conf_temp.chainID != post_run_prot.res[kr].chainID ||
                strncmp(conf_temp.resName, post_run_prot.res[kr].resName, 3) ||
            conf_temp.resSeq != post_run_prot.res[kr].resSeq ||
            conf_temp.iCode != post_run_prot.res[kr].iCode)
            continue;
            /* residue found */
            iconf = ins_conf(&post_run_prot.res[kr], post_run_prot.res[kr].n_conf, 0);
            cpy_conf(&post_run_prot.res[kr].conf[iconf], &conf_temp);
            
            notfound = 0;
            break;
        }
        if (notfound) { /* belongs to new residue */
            ires  = ins_res(&post_run_prot, post_run_prot.n_res);
            iconf = ins_conf(&post_run_prot.res[kr], post_run_prot.res[kr].n_conf, 0);

            cpy_conf(&post_run_prot.res[kr].conf[iconf], &conf_temp);

            strncpy(post_run_prot.res[ires].resName, conf_temp.resName,3);
            post_run_prot.res[ires].chainID = conf_temp.chainID;
            post_run_prot.res[ires].resSeq = conf_temp.resSeq;
            post_run_prot.res[ires].iCode = conf_temp.iCode;
        }
    }
    
    return 0;
}

void postrun_group_confs()
{   int  kr, kc;
    int free_conf;
    int ir, ic;
    
    
    /* delete existing big list when regrouping */
    if (post_run_biglist != NULL) {
        for (ir=0; ir<post_run_n_free; ir++) {
            if (post_run_biglist[ir].res != NULL)
                free(post_run_biglist[ir].res);
        }
        free(post_run_biglist);
    }


    /* create post_run_all_res */
    if (post_run_all_res != NULL) {
        for (ir=0; ir<post_run_n_all; ir++) {
            if (post_run_all_res[ir].conf != NULL)
                free(post_run_all_res[ir].conf);
        }
        free(post_run_all_res);
    }
    post_run_all_res = (RESIDUE *) malloc(post_run_prot.n_res * sizeof(RESIDUE));
    post_run_n_all = post_run_prot.n_res;
    for (ir=0; ir<post_run_n_all; ir++) {
        post_run_all_res[ir].n    = post_run_prot.res[ir].n_conf;
        post_run_all_res[ir].conf = (int *) malloc(post_run_prot.res[ir].n_conf * sizeof(int));
        post_run_all_res[ir].mfe.crg = (float *) calloc(env.titr_steps, sizeof(float));
    }
    for (ir=0; ir<post_run_n_all; ir++) {
        for (ic=0; ic<post_run_all_res[ir].n; ic++) {
            post_run_all_res[ir].conf[ic] = post_run_prot.res[ir].conf[ic].iConf;
        }
    }
    
    /* create free residues */
    if (post_run_free_res != NULL) {
        for (ir=0; ir<post_run_n_free; ir++) {
            free(post_run_free_res[ir].conf);
        }
        free(post_run_free_res);
        post_run_free_res = NULL;
    }
    
    if (post_run_fixed_res != NULL) {
        for (ir=0; ir<post_run_n_fixed; ir++) {
            free(post_run_fixed_res[ir].conf);
        }
        free(post_run_fixed_res);
        post_run_fixed_res = NULL;
    }
    
    post_run_n_free = 0; post_run_n_fixed = 0;
    for (kr=0; kr<post_run_n_all; kr++) {
        free_conf = 0;
        for (kc=0; kc<post_run_all_res[kr].n; kc++) {
            if (post_run_conflist.conf[post_run_all_res[kr].conf[kc]].on == 'f') free_conf ++;
            /*
            printf("%s, %c, %f\n", post_run_conflist.conf[post_run_all_res[kr].conf[kc]].uniqID,
            post_run_conflist.conf[post_run_all_res[kr].conf[kc]].on,
            post_run_conflist.conf[post_run_all_res[kr].conf[kc]].occ);
            */
        }
        if (free_conf) {       /* free_conf stores number of free conformers in this residue */
            post_run_n_free++;
            post_run_free_res = (RESIDUE*) realloc(post_run_free_res, post_run_n_free * sizeof(RESIDUE));
            post_run_free_res[post_run_n_free-1].n = free_conf;
            post_run_free_res[post_run_n_free-1].conf = (int *) malloc(free_conf*sizeof(int));
            
            ic = 0;
            for (kc=0; kc<post_run_all_res[kr].n; kc++) {
                if (post_run_conflist.conf[post_run_all_res[kr].conf[kc]].on == 'f') {
                    post_run_free_res[post_run_n_free-1].conf[ic] = post_run_conflist.conf[post_run_all_res[kr].conf[kc]].iConf;
                    ic++;
                }
                else { /* 't' must be with occ = 0.0, ignored */
                    if (fabs(post_run_conflist.conf[post_run_all_res[kr].conf[kc]].occ) > 0.000000001)  {
                        printf("   WARNING: %s has a fixed occ %.3f in a free residue\n",
                        post_run_conflist.conf[post_run_all_res[kr].conf[kc]].uniqID, post_run_conflist.conf[post_run_all_res[kr].conf[kc]].occ);
                        fflush(stdout);
                    }
                }
            }
        }
        else {
            post_run_n_fixed++;
            post_run_fixed_res = (RESIDUE*) realloc(post_run_fixed_res, post_run_n_fixed * sizeof(RESIDUE));
            post_run_fixed_res[post_run_n_fixed-1].n = post_run_all_res[kr].n;
            post_run_fixed_res[post_run_n_fixed-1].conf = (int *) malloc(post_run_all_res[kr].n*sizeof(int));
            
            ic = 0;
            for (kc=0; kc<post_run_all_res[kr].n; kc++) {
                post_run_fixed_res[post_run_n_fixed-1].conf[ic] = post_run_conflist.conf[post_run_all_res[kr].conf[kc]].iConf;
                ic++;
            }
        }
    }
    
    return;
}

int postrun_load_pairwise()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, ".", 0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    if (!(post_run_pairwise = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    

    for (kc=0; kc<ematrix.n; kc++) {
        if (!(post_run_pairwise[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
          
           post_run_pairwise[i][j] = post_run_pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele
                                            +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           
        }
    }

    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

int postrun_load_pairwise_vdw()
{
    int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, ".", 0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    if (!(post_run_pairwise_vdw = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    for (kc=0; kc<ematrix.n; kc++) {
        if (!(post_run_pairwise_vdw[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           post_run_pairwise_vdw[i][j] = post_run_pairwise_vdw[j][i] = (ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw/2.0 ;
        }
    }
    
    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

int postrun_fitit()
{   
    int i, j;
    int N_crg;
    int N_res;
    int Counter;
    char line[20];
    struct STAT stat;   /* a structure of statistcs */
    char **head, **mhead;  // add mhead for mfe header
    float *netcrg;
    float **ypp, **ysp;
    float a, b, mid;
    float *crg, *protons, *electrons;    /* total net charge at pHs */
    float n;
    char sbuff[21];
    int n_protons, n_electrons, n_crg;
    int n_protons_grnd, n_electrons_grnd, n_crg_grnd;
    int k; // just for cycle
    int tifound; // whether do postrun_fitit
    FILE *blist_fp;
    float mfePoint;
 
    /*<<< Get titration points >>>*/
    /*--- Determine the titration points, x values ---*/
    post_run_Nx = env.titr_steps;
    
    /*--- Assign titration points to an array ---*/
    post_run_xp = (float *) malloc(post_run_Nx * sizeof(float));      /* store x points */
    post_run_yp = (float *) malloc(post_run_Nx * sizeof(float));      /* store y points */
    ypp = (float **) malloc(post_run_conflist.n_conf * sizeof(float *));
    for (i=0; i<post_run_conflist.n_conf; i++) ypp[i] = (float *) malloc(post_run_Nx * sizeof(float));
    ysp = (float **) malloc(post_run_conflist.n_conf * sizeof(float *));
    for (i=0; i<post_run_conflist.n_conf; i++) ysp[i] = (float *) malloc(post_run_Nx * sizeof(float));

    /*--- Convert x points to float ---*/
    for (i=0; i<post_run_Nx; i++) {
        if (env.titr_type == 'p') post_run_xp[i] = env.titr_ph0 + i*env.titr_phd;
        else  post_run_xp[i] = env.titr_eh0 + i*env.titr_ehd;
    }
    
    head = (char **) malloc(post_run_conflist.n_conf * sizeof(char *));
    for (i=0; i<post_run_conflist.n_conf; i++) head[i] = (char *) malloc(20 * sizeof(char));
    post_run_shead = (char **) malloc(post_run_conflist.n_conf * sizeof(char *));
    for (i=0; i<post_run_conflist.n_conf; i++) post_run_shead[i] = (char *) malloc(20 * sizeof(char));
    mhead = (char **) malloc(post_run_conflist.n_conf * sizeof(char *));
    for (i=0; i<post_run_conflist.n_conf; i++) mhead[i] = (char *) malloc(20 * sizeof(char));
    netcrg = (float *) malloc(post_run_conflist.n_conf*sizeof(float));

    
    /*<<< Group into residues >>>*/
    /* keep only charged conformers */
    Counter = 0;
    for (i=0; i <post_run_conflist.n_conf; i++) {
        if (strchr(post_run_conflist.conf[i].uniqID, '+') || strchr(post_run_conflist.conf[i].uniqID, '-')) {
            strncpy(head[Counter], post_run_conflist.conf[i].uniqID, 4); head[Counter][4] = '\0';
            strncat(head[Counter], post_run_conflist.conf[i].uniqID+5, 6);
            head[Counter][10] = '\0';
            for (j=0; j<post_run_Nx; j++) {
                ypp[Counter][j] = post_run_occ_table[i][j];
            }
            Counter++;
        }
    }
    N_crg = Counter;

    
    /* group residues */
    Counter = 0;
    strncpy(line, head[0], 10);
    for (j=0; j<post_run_Nx; j++) ysp[0][j] = ypp[0][j];
    for (i=1; i<N_crg; i++) {
        if (strncmp(head[i], line, 10)) {      /* not equal, a new residue */
            strncpy(post_run_shead[Counter], line, 10); post_run_shead[Counter][10] = '\0';
            Counter++;
            strncpy(line, head[i], 10);
            for (j=0; j<post_run_Nx; j++) ysp[Counter][j] = ypp[i][j];
        }
        else {                                 /* same residue */
            for (j=0; j<post_run_Nx; j++) ysp[Counter][j] += ypp[i][j];
        }
    }
    strncpy(post_run_shead[Counter], line, 10); post_run_shead[Counter][10] = '\0';
    N_res = Counter + 1;
    post_run_n_mfe = N_res;
    for (i=0; i<post_run_n_all; i++) {
        for (j=0; j<post_run_Nx; j++) post_run_all_res[i].mfe.crg[j] = 0.0;
    }
    // N_res are the number we do mfe
    for (i=0; i<post_run_n_mfe; i++) {
        strncpy(mhead[i], post_run_shead[i], 3); mhead[i][3] = '\0';
        strncat(mhead[i], post_run_shead[i]+4, 6); mhead[i][9] = '\0';
    }
    if (env.mfe_flag) {
        if (env.titr_type == 'p') printf("       Doing mfe at pH %.2f for all the residues\n", env.mfe_point);
        else printf("       Doing mfe at Eh %.3f for all the residues\n", env.mfe_point);
        if ((env.mfe_point>post_run_xp[0] && env.mfe_point>post_run_xp[post_run_Nx-1]) || (env.mfe_point<post_run_xp[0] && env.mfe_point<post_run_xp[post_run_Nx-1]))
            printf("         MFE: mfe point not in the titration range, do mfe at pKa or Em\n");
    }
    else printf("   MFE: didn't specify mfe point, do mfe at pKa or Em \n");
    
    /* free post_run_pairwise */
    for (i=0; i<post_run_conflist.n_conf; i++)
        free(post_run_pairwise[i]);
    free(post_run_pairwise);

    // load the post_run_pairwise interaction again, round the ele and vdw to keep consistent with mfe.py  
    if (postrun_load_pairwise_fround3()) {
        printf("   FATAL: mfe post_run_pairwise interaction not loaded\n");
        return USERERR;
    }

    // load all the conformers of mfe residues
    post_run_mfe_res = (RESIDUE *) calloc(post_run_n_mfe, sizeof(RESIDUE));
    for (k=0; k<post_run_n_mfe; k++) {
        post_run_mfe_res[k].n = 0;
        for (i=0; i<post_run_conflist.n_conf; i++) {
            strncpy(sbuff, post_run_conflist.conf[i].uniqID, 3); sbuff[3] = '\0';
            strncat(sbuff, post_run_conflist.conf[i].uniqID+5, 6); sbuff[9] = '\0';
            if (!strcmp(mhead[k], sbuff)) post_run_mfe_res[k].n++;
        }
    }
    for (i=0; i<post_run_n_mfe; i++) {
        post_run_mfe_res[i].conf = (int *) calloc(post_run_mfe_res[i].n, sizeof(int));
        post_run_mfe_res[i].mfe.mfePair = (float *) calloc(post_run_n_all, sizeof(float));
        post_run_mfe_res[i].mfe.vdw = (float *) calloc(post_run_n_all, sizeof(float));
        post_run_mfe_res[i].mfe.ele = (float *) calloc(post_run_n_all, sizeof(float));
        post_run_mfe_res[i].mfe.crg = (float *) calloc(env.titr_steps, sizeof(float));
    }

    for (i=0; i<post_run_n_mfe; i++) {
        post_run_mfe_res[i].n = 0;
        for (j=0; j<post_run_conflist.n_conf; j++) {
            strncpy(sbuff, post_run_conflist.conf[j].uniqID, 3); sbuff[3] = '\0';
            strncat(sbuff, post_run_conflist.conf[j].uniqID+5, 6); sbuff[9] = '\0';
            if (!strcmp(mhead[i], sbuff)) {
                post_run_mfe_res[i].n++;
                post_run_mfe_res[i].conf[post_run_mfe_res[i].n-1] = post_run_conflist.conf[j].iConf;
            }
        }
    }
    for (i=0; i<post_run_n_mfe; i++) { 
        for (j=0; j<post_run_mfe_res[i].n; j++) {
            post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_self0 = post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_vdw0
                                                      + post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_vdw1 
                                                      + post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_epol 
                                                      + post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_tors 
                                                      + post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_dsolv
                                                      + post_run_conflist.conf[post_run_mfe_res[i].conf[j]].E_extra;
        }
    }

    /* Write net charge ,protons, electrons to sum_crg.out, to be implemented */
    crg       = (float *) malloc(post_run_Nx * sizeof(float));
    protons   = (float *) malloc(post_run_Nx * sizeof(float));
    electrons = (float *) malloc(post_run_Nx * sizeof(float));

    memset(crg, 0, post_run_Nx * sizeof(float));
    memset(protons, 0, post_run_Nx * sizeof(float));
    memset(electrons, 0, post_run_Nx * sizeof(float));
    /*
    if ((post_run_fp = fopen("sum_crg.out", "w")) == NULL) {
        printf("   Can not open \"sum_crg.out\" to write, abort ...\n");
        return USERERR;
    }
    if (env.titr_type == 'p') {   // pH titration 
        fprintf(post_run_fp, "  pH      ");
    }
    else {      // Eh titration assumed
        fprintf(post_run_fp, "  Eh      ");
    }
    for (i=0; i<post_run_Nx; i++) fprintf(post_run_fp, " %5d", (int) post_run_xp[i]);
    fprintf(post_run_fp, "\n");

    for(i=0; i<N_res; i++) {
        //printf("post_run_shead: %s\n", post_run_shead[i]);
        fprintf(post_run_fp, "%s", post_run_shead[i]);
        strncpy(sbuff, post_run_shead[i], 4); sbuff[4] = '1'; sbuff[5] = '\0';
        if (param_get( "PROTON", sbuff, "", &n_protons)) n_protons = 0;
        if (param_get( "ELECTRON", sbuff, "", &n_electrons)) n_electrons = 0;
        // n_crg = n_protons-n_electrons; 
        n_crg = n_protons-n_electrons; 
        
        // number of protons on ground conformer type 
        strncpy(sbuff, post_run_shead[i], 3); sbuff[3] = '0'; sbuff[4] = '1'; sbuff[5] = '\0';
        if (param_get( "PROTON", sbuff, "", &n_protons_grnd)) n_protons_grnd = 0;
        if (param_get( "ELECTRON", sbuff, "", &n_electrons_grnd)) n_electrons_grnd = 0;
        n_crg_grnd = n_protons_grnd - n_electrons_grnd;
        
        for(j=0; j<post_run_Nx; j++) {
            fprintf(post_run_fp, " %5.2f", n_crg*ysp[i][j]);
            post_run_mfe_res[i].mfe.crg[j] = n_crg*ysp[i][j] + n_crg_grnd*(1.0-ysp[i][j]);
            for (k=0; k<post_run_n_all; k++) {
                strncpy(sbuff, post_run_conflist.conf[post_run_all_res[k].conf[0]].uniqID, 3); sbuff[3] = '\0';
                strncat(sbuff, post_run_conflist.conf[post_run_all_res[k].conf[0]].uniqID+5, 6); sbuff[9] = '\0';
                // get all the charges of all the residues
                if (!strcmp(sbuff, mhead[i])) post_run_all_res[k].mfe.crg[j] = n_crg*ysp[i][j] + n_crg_grnd*(1.0-ysp[i][j]);
            }    
            crg[j] += n_crg*ysp[i][j] + n_crg_grnd*(1.0-ysp[i][j]);
            protons[j] += n_protons*ysp[i][j] + n_protons_grnd*(1.0-ysp[i][j]);
            electrons[j] += n_electrons*ysp[i][j] + n_electrons_grnd*(1.0-ysp[i][j]);
        }
        fprintf(post_run_fp, "\n");
    }

    fprintf(post_run_fp, "----------\n");

    fprintf(post_run_fp, "Net_Charge");
    for(j=0; j<post_run_Nx; j++) {
        fprintf(post_run_fp, " %5.2f", crg[j]);
    }
    fprintf(post_run_fp, "\n");

    fprintf(post_run_fp, "Protons   ");
    for(j=0; j<post_run_Nx; j++) {
        fprintf(post_run_fp, " %5.2f", protons[j]);
    }
    fprintf(post_run_fp, "\n");

    fprintf(post_run_fp, "Electrons ");
    for(j=0; j<post_run_Nx; j++) {
        fprintf(post_run_fp, " %5.2f", electrons[j]);
    }
    fprintf(post_run_fp, "\n");

    fclose(post_run_fp);*/

    /*<<< Loop over y values >>>*/
    if (!(post_run_fp = fopen(CURVE_FITTING, "w"))) {
        printf("   FATAL: can not write to file \"%s\".", CURVE_FITTING);
        return USERERR;
    }
    if (env.titr_type == 'p') {   /* pH titration */
        fprintf(post_run_fp, "%10s", "pH");
    }
    else {      /* Eh titration assumed */
        fprintf(post_run_fp, "%10s", "Eh");
    }
    if (env.mfe_pka){
        fprintf(post_run_fp, "    pKa/Em  n(slope) 1000*chi2      vdw0    vdw1    tors    ebkb    dsol   offset  pHpK0   EhEm0    -TS   residues   total");
    } else {
        fprintf(post_run_fp, "%10s%10s%10s", 
                "pKa/Em", 
                "n(slope)", 
                "1000*chi2");
        fprintf(post_run_fp, "%10s%10s%10s%10s%10s%10s%10s%10s%10s", 
                "offset", 
                "tors", 
                "v_s", 
                "v_b",
                "v_r",
                "dsol", 
                "e_b", 
                "e_r", 
                "total");
    }
    fprintf(post_run_fp, "\n");

    if (!(blist_fp = fopen("respair.lst", "w"))) {
        printf("can't write file respair.lst\n");
    }
    fprintf(blist_fp, " residue    partner         vdw     ele  post_run_pairwise  charge\n");



    for (i=0; i <N_res; i++) {
        /*--- Convert y points to float ---*/

        /* a reasonable guess makes optimization easier */
        a = 0.0;
        mid = 0.6;
        for (j=0; j<post_run_Nx; j++) {
            post_run_yp[j] = ysp[i][j];
            if (fabs(post_run_yp[j] - 0.5) < mid) {
                mid = fabs(post_run_yp[j] - 0.5);
                b = post_run_xp[j];
            }
        }
        tifound=1;
        if (mid >= 0.485) {
            tifound=0;
            if (fabs(post_run_yp[0]-post_run_yp[post_run_Nx-1])>0.5) {  /* jumps form <0.015 to >0.985 */
               fprintf(post_run_fp, "%s     %-25s", post_run_shead[i],"titration curve too sharp");
               mfePoint=post_run_xp[post_run_Nx/2];
            }  
            else {    // assume either all the post_run_yp >0.985 or all of them <0.015                    
                if (post_run_yp[0] > 0.985) {
                    if (strchr(post_run_shead[i], '+')) {
                        fprintf(post_run_fp, "%s     >%-24.1f", post_run_shead[i], post_run_xp[post_run_Nx-1]);
                        mfePoint=post_run_xp[post_run_Nx-1];
                    }
                    else {
                        fprintf(post_run_fp, "%s     <%-24.1f", post_run_shead[i], post_run_xp[0]);
                        mfePoint=post_run_xp[0];
                    }
                }
                else {
                    if (strchr(post_run_shead[i], '+')) {
                        fprintf(post_run_fp, "%s     <%-24.1f", post_run_shead[i], post_run_xp[0]);
                        mfePoint=post_run_xp[0];
                    }
                    else {
                        fprintf(post_run_fp, "%s     >%-24.1f", post_run_shead[i], post_run_xp[post_run_Nx-1]);
                        mfePoint=post_run_xp[post_run_Nx-1];
                    }
                }
            }            
        }
        if (tifound == 0) {
            postrun_print_mfe(i, mfePoint, post_run_fp, blist_fp);
            continue;
        }
         //   fprintf(post_run_fp, "
        /* printf("a=%.3f; b=%.3f\n",a,b); */
        stat = postrun_fit(a, b);

        if (env.titr_type == 'p') n = fabs(stat.a * 8.617342E-2 * env.monte_temp / 58.0);
        else if (env.titr_type == 'e') n = fabs(stat.a * 8.617342E-2 * env.monte_temp);
        else n = stat.a;


        //if (stat.b < post_run_xp[0] || stat.b > post_run_xp[post_run_Nx-1]) fprintf(post_run_fp, "%s        pKa or Em out of range   \n", post_run_shead[i]);
        if (stat.b < post_run_xp[0]) {
            fprintf(post_run_fp, "%s     <%-24.1f", post_run_shead[i], post_run_xp[0]);
            mfePoint = post_run_xp[0];
        }
    else if (stat.b > post_run_xp[post_run_Nx-1]) {
            fprintf(post_run_fp, "%s     >%-24.1f", post_run_shead[i], post_run_xp[post_run_Nx-1]);
            mfePoint = post_run_xp[post_run_Nx-1];
        }
        else {
            fprintf(post_run_fp, "%10s %9.3f %9.3f %9.3f", post_run_shead[i], stat.b, n, 1000*stat.chi2);
            mfePoint = stat.b;
        }
        postrun_print_mfe(i, mfePoint, post_run_fp, blist_fp); 
    }     

    fclose(post_run_fp);
    free(crg);
    free(protons);
    free(electrons);
    free(post_run_xp);
    free(post_run_yp);
    for (i=0; i<post_run_conflist.n_conf; i++) free(ypp[i]);
    free(ypp);
    for (i=0; i<post_run_conflist.n_conf; i++) free(ysp[i]);
    free(ysp);
    for (i=0; i<post_run_conflist.n_conf; i++) free(head[i]);
    free(head);
    for (i=0; i<post_run_conflist.n_conf; i++) free(post_run_shead[i]);
    free(post_run_shead);
    
    return 0;
}

struct node * createNode(int key_resSeq, char *name, float the_occupancy, char node_chainID, int node_resSeq, int node_iConf) {
    struct node *newnode;
    newnode = (struct node *) malloc(sizeof(struct node));
    newnode->key_resSeq = key_resSeq;
    newnode->the_occupancy = the_occupancy;
    newnode->node_chainID = node_chainID;
    newnode->node_resSeq = node_resSeq;
    newnode->node_iConf = node_iConf;
    strcpy(newnode->name, name);
    newnode->next = NULL;
    return newnode;
}

void insertToHash(int key_resSeq, char *name, float the_occupancy, char node_chainID, int node_resSeq) {
    int hashIndex = key_resSeq % eleCount;
    int int_confId;
    char *confId = (char*) malloc(3);
    strncpy(confId, name+11, 15);
    int_confId = (int) strtol(confId, (char **)NULL, 10);
    
    struct node *newnode = createNode(key_resSeq, name, the_occupancy, node_chainID, node_resSeq, int_confId);

    /* head of list for the bucket with index "hashIndex" */
    if (!hashTable[hashIndex].head) {
        hashTable[hashIndex].head = newnode;
        hashTable[hashIndex].count = 1;
        return;
    }
    /* adding new node to the list */
    newnode->next = (hashTable[hashIndex].head);
    /*
     * update the head of the list and no of
     * nodes in the current bucket
     */
    hashTable[hashIndex].head = newnode;
    hashTable[hashIndex].count++;
    return;
}

struct OccInRes searchInHash(int key_resSeq) {
    int hashIndex = key_resSeq % eleCount, flag = 0, i = 0;
    struct node *myNode;
    myNode = hashTable[hashIndex].head;
    float max=0.0;
    char *most_occ_confname;
    char chainID;
    int resSeq;
    int iConf;

    if (!myNode) {
        printf("Search element unavailable in hash table\n");
        struct OccInRes result;
        result.key_resSeq = "most_occ_confname";
        result.most_occupancy = 999.00;
        result.OccInRes_chainID = 'A';
        result.OccInRes_resSeq = 0;
        result.OccInRes_iConf = 0;
        return result;
    }
    while (myNode != NULL) {
        if (myNode->key_resSeq == key_resSeq) {
            //printf("VoterID  : %d\n", myNode->key_resSeq);
            //printf("Name     : %s\n", myNode->name);
            //printf("Age      : %5.2f\n", myNode->the_occupancy);
            flag = 1;
            if (myNode->the_occupancy > max){
                max = myNode->the_occupancy;
                most_occ_confname = myNode->name;
                chainID = myNode->node_chainID;
                resSeq = myNode->node_resSeq;
                iConf = myNode->node_iConf;
                //break;
            }
        }
        myNode = myNode->next;
    }

    if (!flag)
        printf("Search element unavailable in hash table\n");
    
    struct OccInRes result;
    result.key_resSeq = most_occ_confname;
    result.most_occupancy = max;
    result.OccInRes_chainID = chainID;
    result.OccInRes_resSeq = resSeq;
    result.OccInRes_iConf = iConf;
    //printf("%d %s max: %5.2f\n",key_resSeq, most_occ_confname, max);    
    return result;
}

int potential_map(){

    FILE *post_run_fp, *fp2, *fp3, *fort11, *fort12, *step5_out, *fort10;
    char sbuff[MAXCHAR_LINE];
    int i, j, k,c, iConf, n_retry, ic, a, most_occ_conf_list;
    int conf_counter;
    int titration_point;
    char *tok;
    float protein_charge;
    float residue_charge;
    int occ_value;
    protein_charge = 0.0;
    residue_charge = 0.0;
    n_retry = 0; /* reset delphi failure counter for this conformer */
    a = 0;
    c = 0;
    struct OccInRes result;

    eleCount = post_run_conflist.n_conf*sizeof(int *);
    hashTable = (struct hash *) calloc(eleCount, sizeof(struct hash));

    if (!(post_run_fp=fopen(STEP2_OUT, "r"))) {
       printf(ANSI_COLOR_RED "   FATAL: potential_map(): \"No step 2 output \"%s\".\n" ANSI_COLOR_RESET, STEP2_OUT);
       return USERERR;
    }
    post_run_prot = load_pdb(post_run_fp);
    fclose(post_run_fp);
    printf("   Sreaching for the most occupied conformer for each residue in fort.38 at %cH %d\n",env.titr_type,env.column_number);
    
    
    // Writing radii file for delphi
    fort11 = fopen("fort.11", "w");
    fprintf(fort11,"!\n");
    fprintf(fort11,"!   This file was created from MCCE3.5 topology files\n");
    fprintf(fort11,"!    Radius of all hydrogens with\n");
    fprintf(fort11,"! Column 1-6: Atom name; column 7-12: Residuename;\n");
    fprintf(fort11,"! column 13-18: atom size.\n");
    fprintf(fort11,"atom__res_radius_\n");
    fprintf(fort11,"1H          1.250\n");
    fprintf(fort11,"2H          1.250\n");
    fprintf(fort11,"3H          1.250\n");
    fprintf(fort11,"H1          1.250\n");
    fprintf(fort11,"H2          1.250\n");
    fprintf(fort11,"H3          1.250\n");
    fprintf(fort11,"OXT         1.480\n");
    fprintf(fort11,"\n");


    fort12 = fopen("fort.12", "w");
    fprintf(fort12,"!\n");
    fprintf(fort12,"!   This file was created from MCCE3.5 topology files\n");
    fprintf(fort12,"!    Radius of all hydrogens with\n");
    fprintf(fort12,"! Column 1-6: atom name; column 7-9: residue name;\n");
    fprintf(fort12,"! column 10-12: residue number; column 14: chain name;\n");
    fprintf(fort12,"! column 15-23: charge magnitude. For example:\n");
    fprintf(fort12,"! CA    ARG130 B -0.2637 In is a valid format\n");
    fprintf(fort12,"atom__resnumbc_charge_\n");
    fprintf(fort12,"1H              0.0906 !!! Hydrogen of CA backbone\n");
    fprintf(fort12,"2H              0.0906\n");
    fprintf(fort12,"3H              0.0906\n");
    fprintf(fort12,"\n");
    fprintf(fort12,"H1              0.4240 !!! Hydrogen for N-terminal residue\n");
    fprintf(fort12,"H2              0.4240 !!!\n");
    fprintf(fort12,"H3              0.4240 !!!\n");
    fprintf(fort12,"\n");

    step5_out = fopen("step5_out.pdb", "w");
    

    // group conf. into res.
    for (i=0; i<post_run_conflist.n_conf; i++){
        insertToHash(post_run_conflist.conf[i].resSeq, post_run_conflist.conf[i].uniqID, post_run_occ_table[i][env.column_number], post_run_conflist.conf[i].chainID, post_run_conflist.conf[i].resSeq);
    }

    if (env.only_backbone){
        for (i=0; i<post_run_prot.n_res; i++){
            iConf =1;
            result = searchInHash(post_run_prot.res[i].resSeq);
            if (strchr(post_run_prot.res[i].conf[0].confName, 'BK')){
                for (k=0; k<post_run_prot.res[i].conf[0].n_atom; k++){
                    if (strlen(post_run_prot.res[i].conf[0].atom[k].name) != 0){
                        protein_charge += post_run_prot.res[i].conf[0].atom[k].crg;
                        fprintf(fort11,"%-05s %s      %5.2f\n", post_run_prot.res[i].conf[0].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].conf[0].atom[k].rad);
                        fprintf(fort12,"%-05s %s %d      %5.2f\n", post_run_prot.res[i].conf[0].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].resSeq, post_run_prot.res[i].conf[0].atom[k].crg);
                        fprintf(step5_out, "ATOM %5d %4s %4s %4d    %8.3f%8.3f%8.3f%7.3f%7.3f           %1s\n",c, post_run_prot.res[i].conf[0].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].resSeq, post_run_prot.res[i].conf[0].atom[k].xyz.x, post_run_prot.res[i].conf[0].atom[k].xyz.y, post_run_prot.res[i].conf[0].atom[k].xyz.z, 1.00, 0.00, post_run_prot.res[i].conf[0].atom[k].name);
                    }
                }
            }
            for (j=0; j<post_run_prot.res[i].n_conf; j++){            
                        if (result.OccInRes_chainID == post_run_prot.res[i].chainID && result.OccInRes_resSeq == post_run_prot.res[i].resSeq && result.OccInRes_iConf == iConf){
                            for (k=0; k<post_run_prot.res[i].conf[j].n_atom; k++){
                                if (strlen(post_run_prot.res[i].conf[j].atom[k].name) != 0){
                                    fprintf(fort11,"%-05s %s      %5.2f\n", post_run_prot.res[i].conf[j].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].conf[j].atom[k].rad);
                                    if (c<99999) c++;
                                    fprintf(step5_out, "ATOM %5d %4s %4s %4d    %8.3f%8.3f%8.3f%7.3f%7.3f           %1s\n",c, post_run_prot.res[i].conf[j].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].resSeq, post_run_prot.res[i].conf[j].atom[k].xyz.x, post_run_prot.res[i].conf[j].atom[k].xyz.y, post_run_prot.res[i].conf[j].atom[k].xyz.z, 1.00, 0.00, post_run_prot.res[i].conf[j].atom[k].name); 
                                }
                            }
                        }
                        iConf++;
            }
            fprintf(fort11,"\n");
            fprintf(fort12,"\n");
        }
    }
    else {
        for (i=0; i<post_run_prot.n_res; i++){
            iConf =1;
            result = searchInHash(post_run_prot.res[i].resSeq);
            if (strchr(post_run_prot.res[i].conf[0].confName, 'BK')){
                //printf("%c %c %d %d %d %d %d ====\n", result.OccInRes_chainID, post_run_prot.res[i].chainID, result.OccInRes_resSeq, post_run_prot.res[i].resSeq, result.OccInRes_iConf, iConf, post_run_prot.res[i].conf[0].n_atom);
                for (k=0; k<post_run_prot.res[i].conf[0].n_atom; k++){
                    if (strlen(post_run_prot.res[i].conf[0].atom[k].name) != 0){
                        fprintf(fort11,"%-05s %s      %5.2f\n", post_run_prot.res[i].conf[0].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].conf[0].atom[k].rad);
                        fprintf(fort12,"%-05s %s %d      %5.2f\n", post_run_prot.res[i].conf[0].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].resSeq, post_run_prot.res[i].conf[0].atom[k].crg);
                        fprintf(step5_out, "ATOM %5d %4s %4s %4d    %8.3f%8.3f%8.3f%7.3f%7.3f           %1s\n",c, 
                            post_run_prot.res[i].conf[0].atom[k].name, 
                            post_run_prot.res[i].resName, 
                            post_run_prot.res[i].resSeq, 
                            post_run_prot.res[i].conf[0].atom[k].xyz.x, 
                            post_run_prot.res[i].conf[0].atom[k].xyz.y, 
                            post_run_prot.res[i].conf[0].atom[k].xyz.z, 
                            1.00, 
                            0.00, 
                            post_run_prot.res[i].conf[0].atom[k].name); 
                        //printf("%c %c %d %d %d %d %s\n", result.OccInRes_chainID, post_run_prot.res[i].chainID, result.OccInRes_resSeq, post_run_prot.res[i].resSeq, result.OccInRes_iConf, iConf, post_run_prot.res[i].conf[0].atom[k].name);
                    }
                }
            }
            for (j=0; j<post_run_prot.res[i].n_conf; j++){    
            //printf("%c %c %d %d %d %d \n", result.OccInRes_chainID, post_run_prot.res[i].chainID, result.OccInRes_resSeq, post_run_prot.res[i].resSeq, result.OccInRes_iConf, iConf);        
                        if (result.OccInRes_chainID == post_run_prot.res[i].chainID && result.OccInRes_resSeq == post_run_prot.res[i].resSeq){
                            if (atoi(post_run_prot.res[i].conf[j].confID) == result.OccInRes_iConf){
                            //printf("testing %d %s %s %s %c %c %d %d %d %d: %f \n",j,post_run_prot.res[i].conf[j].confID, post_run_prot.res[i].conf[j].confName,  post_run_prot.res[i].resName, result.OccInRes_chainID, post_run_prot.res[i].chainID, result.OccInRes_resSeq, post_run_prot.res[i].resSeq, result.OccInRes_iConf, iConf, result.most_occupancy);
                                for (k=0; k<post_run_prot.res[i].conf[j].n_atom; k++){
                                    if (strlen(post_run_prot.res[i].conf[j].atom[k].name) != 0){
                                        fprintf(fort11,"%-05s %s      %5.2f\n", post_run_prot.res[i].conf[j].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].conf[j].atom[k].rad);
                                        protein_charge += post_run_prot.res[i].conf[j].atom[k].crg;
                                        fprintf(fort12,"%-05s %s %d      %5.2f\n", post_run_prot.res[i].conf[j].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].resSeq, post_run_prot.res[i].conf[j].atom[k].crg);
                                        if (c<99999) c++;
                                        fprintf(step5_out, "ATOM %5d %4s %4s %4d    %8.3f%8.3f%8.3f%7.3f%7.3f           %1s\n",c, post_run_prot.res[i].conf[j].atom[k].name, post_run_prot.res[i].resName, post_run_prot.res[i].resSeq, post_run_prot.res[i].conf[j].atom[k].xyz.x, post_run_prot.res[i].conf[j].atom[k].xyz.y, post_run_prot.res[i].conf[j].atom[k].xyz.z, 1.00, 0.00, post_run_prot.res[i].conf[j].atom[k].name); 
                                        //printf("%c %c %d %d %d %d %d\n", result.OccInRes_chainID, post_run_prot.res[i].chainID, result.OccInRes_resSeq, post_run_prot.res[i].resSeq, result.OccInRes_iConf, iConf, post_run_prot.res[i].conf[j].n_atom);        
                                    }
                                }
                            }
                        }
                        iConf++;
            }
            fprintf(fort11,"\n");
            fprintf(fort12,"\n");
        }
    }

    printf("   Protein total charge = %5.2f\n", protein_charge);
    fclose(fort11);
    fclose(fort12);
    fclose(step5_out);

    fort10 = fopen("fort.10", "w");
    //fprintf(post_run_fp, "gsize=%d\n", env.grids_delphi);
    fprintf(fort10, "scale=%.6f\n", (env.grids_per_ang));
    fprintf(fort10, "perfil=70.0\n");
    fprintf(fort10, "in(pdb,file=\"step5_out.pdb\")\n");
    fprintf(fort10, "in(siz,file=\"fort.11\")\n");
    fprintf(fort10, "in(crg,file=\"fort.12\")\n");
    fprintf(fort10, "indi=%.1f\n", env.epsilon_prot);
    fprintf(fort10, "exdi=%.1f\n", env.epsilon_solv);
    fprintf(fort10, "ionrad=%.1f\n", env.ionrad);
    fprintf(fort10, "prbrad=%.1f\n", env.radius_probe);
    fprintf(fort10, "salt=%.2f\n", env.salt);
    fprintf(fort10, "bndcon=2\n");
    fprintf(fort10, "out(phi,file=phimap01.cube, format=cube)\n");
    fprintf(fort10, "out(eps,file=epsmap01.cube)\n");
    fprintf(fort10, "energy(s,c,g)\n");
    fclose(fort10);
    
    sprintf(sbuff, "%s>delphi%02d.log 2>/dev/null", env.delphi_potential_exe, 1);
    printf("   Running DelPhi to get potential map\n");
    system(sbuff);

    printf("   Done\n");

    return 0;
}

float postrun_fround3(float x)
{
    // round off the float point number
    char sbuff[128];
    float y;

    sprintf(sbuff, "%.3f", x);
    y=atof(sbuff);
    return y;
}

int postrun_load_pairwise_fround3()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, ".", 0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }
    
    if (!(post_run_pairwise = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    } 
    post_run_pairwise_vdw = (float **) malloc(ematrix.n * sizeof(float *));
    post_run_pairwise_ele = (float **) malloc(ematrix.n * sizeof(float *));
    
    for (kc=0; kc<ematrix.n; kc++) {
        if (!(post_run_pairwise[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
        if (!(post_run_pairwise_vdw[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
        if (!(post_run_pairwise_ele[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
          // post_run_pairwise[i][j] = post_run_pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele
          //                                 +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
          // if (ematrix.pw[i][j].vdw > 998.0) post_run_pairwise[i][j] = 999.0;
            post_run_pairwise[i][j] = postrun_fround3(ematrix.pw[i][j].ele) * env.scale_ele + postrun_fround3(ematrix.pw[i][j].vdw) * env.scale_vdw;
            post_run_pairwise_vdw[i][j] = postrun_fround3(ematrix.pw[i][j].vdw) * env.scale_vdw;
            post_run_pairwise_ele[i][j] = postrun_fround3(ematrix.pw[i][j].ele) * env.scale_ele;
        }
    }

    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

int postrun_print_mfe(int i_res, float mfeP, FILE *pK_fp, FILE *res_fp)
{
    int i_low, i_high;
    int j, k;
    float cut_off, ratio;
    char head1[20], head2[20];
    RESIDUE old_mfe;

    if (env.mfe_flag) {
        if (! ((env.mfe_point>post_run_xp[0] && env.mfe_point>post_run_xp[post_run_Nx-1]) || (env.mfe_point<post_run_xp[0] && env.mfe_point<post_run_xp[post_run_Nx-1]))) mfeP = env.mfe_point;
    }
    
    cut_off = env.mfe_cutoff;
    // if (cut_off < 0.0001) cut_off = 0.5; 

    for (k=0; k<=post_run_Nx-1; k++) {
        if (post_run_xp[k] >= mfeP) {
            i_high=k; break;
        }
    }
    for (k=post_run_Nx-1; k>=0; k--) {
        if (post_run_xp[k] <= mfeP) {
            i_low=k; break;
        }
    }

    if (abs(i_high-i_low) == 0) { // do mfe at only one point
        postrun_get_mfe(i_res, i_low); 
        if (env.titr_type == 'p')
            if (env.mfe_pka){
                fprintf(pK_fp, "  %8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f\n", 
                   post_run_mfe_res[i_res].mfe.vdw0/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.vdw1/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.tors/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.ebkb/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.dsol/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.offset/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.pHpK0/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.EhEm0/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.TS/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.residues/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.total/PH2KCAL);
            } else {
                fprintf(pK_fp, "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                    post_run_mfe_res[i_res].mfe.offset/PH2KCAL,
                    post_run_mfe_res[i_res].mfe.tors/PH2KCAL,
                    post_run_mfe_res[i_res].mfe.vdw0/PH2KCAL, 
                    post_run_mfe_res[i_res].mfe.vdw1/PH2KCAL,
                    post_run_mfe_res[i_res].mfe.vpws/PH2KCAL,
                    post_run_mfe_res[i_res].mfe.dsol/PH2KCAL, 
                    post_run_mfe_res[i_res].mfe.ebkb/PH2KCAL,
                    post_run_mfe_res[i_res].mfe.epws/PH2KCAL,
                    post_run_mfe_res[i_res].mfe.total/PH2KCAL);
            }        
        else 
            if (env.mfe_pka){
                fprintf(pK_fp, "  %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%10.1f\n", 
                   post_run_mfe_res[i_res].mfe.vdw0/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.vdw1/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.tors/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.ebkb/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.dsol/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.offset/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.pHpK0/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.EhEm0/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.TS/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.residues/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.total/mev2Kcal);
            }else {
                fprintf(pK_fp, "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                       post_run_mfe_res[i_res].mfe.offset/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.tors/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.vdw0/PH2KCAL, 
                       post_run_mfe_res[i_res].mfe.vdw1/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.vpws/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.dsol/PH2KCAL, 
                       post_run_mfe_res[i_res].mfe.ebkb/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.epws/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.total/PH2KCAL);
            }
        for (j=0; j<post_run_n_all; j++) {
            head1[0]='\0'; head2[0]='\0';
            if (fabs(post_run_mfe_res[i_res].mfe.mfePair[j]) < cut_off) continue;
            strcpy(head1, post_run_shead[i_res]);
            strncpy(head2, post_run_conflist.conf[post_run_all_res[j].conf[0]].uniqID, 3); head2[3] = '\0';
            strncat(head2, post_run_conflist.conf[post_run_all_res[j].conf[0]].uniqID+5, 6); head2[9] = '\0';
            if (env.titr_type == 'p') 
                fprintf(res_fp, "%s  %s   %8.2f%8.2f%8.2f%8.2f\n", head1, head2, 
                post_run_mfe_res[i_res].mfe.vdw[j]/PH2KCAL, post_run_mfe_res[i_res].mfe.ele[j]/PH2KCAL, post_run_mfe_res[i_res].mfe.mfePair[j]/PH2KCAL, post_run_all_res[j].mfe.crg[i_low]);
            else
                fprintf(res_fp, "%s  %s   %8.2f%8.2f%8.2f%8.2f\n", head1, head2, 
                post_run_mfe_res[i_res].mfe.vdw[j]/mev2Kcal, post_run_mfe_res[i_res].mfe.ele[j]/mev2Kcal, post_run_mfe_res[i_res].mfe.mfePair[j]/mev2Kcal, post_run_all_res[j].mfe.crg[i_low]);
        }
    }
    else {    /* get mfe from two titration points */
        ratio = (mfeP - post_run_xp[i_low])/(post_run_xp[i_high]-post_run_xp[i_low]);
        postrun_get_mfe(i_res, i_low);
        old_mfe.n = post_run_mfe_res[i_res].n;
        old_mfe.conf = (int *) calloc(old_mfe.n, sizeof(int));
        old_mfe.mfe.mfePair = (float *) calloc(post_run_n_all, sizeof(float));
        old_mfe.mfe.vdw = (float *) calloc(post_run_n_all, sizeof(float));
        old_mfe.mfe.ele = (float *) calloc(post_run_n_all, sizeof(float));
        old_mfe.mfe.vdw0 = post_run_mfe_res[i_res].mfe.vdw0;
        old_mfe.mfe.vdw1 = post_run_mfe_res[i_res].mfe.vdw1;
        old_mfe.mfe.ebkb = post_run_mfe_res[i_res].mfe.ebkb;
        old_mfe.mfe.tors = post_run_mfe_res[i_res].mfe.tors;
        old_mfe.mfe.dsol = post_run_mfe_res[i_res].mfe.dsol;
        old_mfe.mfe.offset = post_run_mfe_res[i_res].mfe.offset;
        old_mfe.mfe.pHpK0 = post_run_mfe_res[i_res].mfe.pHpK0;
        old_mfe.mfe.EhEm0 = post_run_mfe_res[i_res].mfe.EhEm0;
        old_mfe.mfe.TS = post_run_mfe_res[i_res].mfe.TS;
        old_mfe.mfe.residues = post_run_mfe_res[i_res].mfe.residues;
        old_mfe.mfe.total = post_run_mfe_res[i_res].mfe.total;

        for (k=0; k<post_run_n_all; k++) {
            old_mfe.mfe.mfePair[k] = post_run_mfe_res[i_res].mfe.mfePair[k];                       
            old_mfe.mfe.vdw[k] = post_run_mfe_res[i_res].mfe.vdw[k];                       
            old_mfe.mfe.ele[k] = post_run_mfe_res[i_res].mfe.ele[k];                       
        }
        old_mfe.mfe.vpws = post_run_mfe_res[i_res].mfe.vpws;
        old_mfe.mfe.epws = post_run_mfe_res[i_res].mfe.epws;

        postrun_get_mfe(i_res, i_high);
        post_run_mfe_res[i_res].mfe.vdw0 = old_mfe.mfe.vdw0*(1-ratio) + post_run_mfe_res[i_res].mfe.vdw0*ratio;
        post_run_mfe_res[i_res].mfe.vdw1 = old_mfe.mfe.vdw1*(1-ratio) + post_run_mfe_res[i_res].mfe.vdw1*ratio;
        post_run_mfe_res[i_res].mfe.ebkb = old_mfe.mfe.ebkb*(1-ratio) + post_run_mfe_res[i_res].mfe.ebkb*ratio;
        post_run_mfe_res[i_res].mfe.tors = old_mfe.mfe.tors*(1-ratio) + post_run_mfe_res[i_res].mfe.tors*ratio;
        post_run_mfe_res[i_res].mfe.dsol = old_mfe.mfe.dsol*(1-ratio) + post_run_mfe_res[i_res].mfe.dsol*ratio;
        post_run_mfe_res[i_res].mfe.offset = old_mfe.mfe.offset*(1-ratio) + post_run_mfe_res[i_res].mfe.offset*ratio;
        post_run_mfe_res[i_res].mfe.pHpK0 = old_mfe.mfe.pHpK0*(1-ratio) + post_run_mfe_res[i_res].mfe.pHpK0*ratio;
        post_run_mfe_res[i_res].mfe.EhEm0 = old_mfe.mfe.EhEm0*(1-ratio) + post_run_mfe_res[i_res].mfe.EhEm0*ratio;
        post_run_mfe_res[i_res].mfe.TS = old_mfe.mfe.TS*(1-ratio) + post_run_mfe_res[i_res].mfe.TS*ratio;
        post_run_mfe_res[i_res].mfe.residues = old_mfe.mfe.residues*(1-ratio) + post_run_mfe_res[i_res].mfe.residues*ratio;
        post_run_mfe_res[i_res].mfe.total = old_mfe.mfe.total*(1-ratio) + post_run_mfe_res[i_res].mfe.total*ratio;

        for (k=0; k<post_run_n_all; k++) {
            post_run_mfe_res[i_res].mfe.mfePair[k] = old_mfe.mfe.mfePair[k]*(1-ratio) + post_run_mfe_res[i_res].mfe.mfePair[k]*ratio;
            post_run_mfe_res[i_res].mfe.vdw[k] = old_mfe.mfe.vdw[k]*(1-ratio) + post_run_mfe_res[i_res].mfe.vdw[k]*ratio;
            post_run_mfe_res[i_res].mfe.ele[k] = old_mfe.mfe.ele[k]*(1-ratio) + post_run_mfe_res[i_res].mfe.ele[k]*ratio;
        }
        post_run_mfe_res[i_res].mfe.vpws = old_mfe.mfe.vpws*(1-ratio) + post_run_mfe_res[i_res].mfe.vpws*ratio;
        post_run_mfe_res[i_res].mfe.epws = old_mfe.mfe.epws*(1-ratio) + post_run_mfe_res[i_res].mfe.epws*ratio;
        
        if (env.titr_type == 'p') 
            if (env.mfe_pka){
                fprintf(pK_fp, "  %8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f\n", 
                   post_run_mfe_res[i_res].mfe.vdw0/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.vdw1/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.tors/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.ebkb/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.dsol/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.offset/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.pHpK0/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.EhEm0/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.TS/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.residues/PH2KCAL, 
                   post_run_mfe_res[i_res].mfe.total/PH2KCAL);
            }else {
                fprintf(pK_fp, "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                       post_run_mfe_res[i_res].mfe.offset/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.tors/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.vdw0/PH2KCAL, 
                       post_run_mfe_res[i_res].mfe.vdw1/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.vpws/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.dsol/PH2KCAL, 
                       post_run_mfe_res[i_res].mfe.ebkb/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.epws/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.total/PH2KCAL);
            }   
        else 
            if (env.mfe_pka){
                fprintf(pK_fp, "  %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%10.1f\n", 
                   post_run_mfe_res[i_res].mfe.vdw0/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.vdw1/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.tors/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.ebkb/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.dsol/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.offset/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.pHpK0/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.EhEm0/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.TS/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.residues/mev2Kcal, 
                   post_run_mfe_res[i_res].mfe.total/mev2Kcal);
            }else {
                fprintf(pK_fp, "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                       post_run_mfe_res[i_res].mfe.offset/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.tors/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.vdw0/PH2KCAL, 
                       post_run_mfe_res[i_res].mfe.vdw1/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.vpws/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.dsol/PH2KCAL, 
                       post_run_mfe_res[i_res].mfe.ebkb/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.epws/PH2KCAL,
                       post_run_mfe_res[i_res].mfe.total/PH2KCAL);
            }

        for (j=0; j<post_run_n_all; j++) {
            head1[0]='\0'; head2[0]='\0';
            if (fabs(post_run_mfe_res[i_res].mfe.mfePair[j]) < cut_off) continue;
            strcpy(head1, post_run_shead[i_res]);
            strncpy(head2, post_run_conflist.conf[post_run_all_res[j].conf[0]].uniqID, 3); head2[3] = '\0';
            strncat(head2, post_run_conflist.conf[post_run_all_res[j].conf[0]].uniqID+5, 6); head2[9] = '\0';
            if (env.titr_type == 'p') 
                fprintf(res_fp, "%s  %s   %8.2f%8.2f%8.2f%8.2f\n", head1, head2, 
                       post_run_mfe_res[i_res].mfe.vdw[j]/PH2KCAL, post_run_mfe_res[i_res].mfe.ele[j]/PH2KCAL, post_run_mfe_res[i_res].mfe.mfePair[j]/PH2KCAL, 
                       (post_run_all_res[j].mfe.crg[i_low]*(1-ratio) + post_run_all_res[j].mfe.crg[i_high]*ratio));
            else
                fprintf(res_fp, "%s  %s   %8.1f%8.1f%8.1f%8.2f\n", head1, head2, 
                       post_run_mfe_res[i_res].mfe.vdw[j]/mev2Kcal, post_run_mfe_res[i_res].mfe.ele[j]/mev2Kcal, post_run_mfe_res[i_res].mfe.mfePair[j]/mev2Kcal, 
                       (post_run_all_res[j].mfe.crg[i_low]*(1-ratio) + post_run_all_res[j].mfe.crg[i_high]*ratio));
        }
     }

     return 0;
}
          
int postrun_get_mfe(int i_res, int t_point)
{
    // get the mfe of the ith post_run_mfe_res at titration point t_point
    int j, k, q;
    float post_run_ph, post_run_eh;
    double *nocc, *rocc, socc[2]={0.0, 0.0};
    double Eref[2];
    MFE E_ground, E_ionize;
    float *mfe, *epws, *vpws;

    memset(&E_ground, 0, sizeof(MFE));
    memset(&E_ionize, 0, sizeof(MFE));
    memset(&socc, 0, 2*sizeof(float));
    E_ground.mfePair = (float *) calloc(post_run_n_all, sizeof(float));
    E_ground.vdw = (float *) calloc(post_run_n_all, sizeof(float));
    E_ground.ele = (float *) calloc(post_run_n_all, sizeof(float));
    E_ionize.mfePair = (float *) calloc(post_run_n_all, sizeof(float));
    E_ionize.vdw = (float *) calloc(post_run_n_all, sizeof(float));
    E_ionize.ele = (float *) calloc(post_run_n_all, sizeof(float));
    mfe = (float *) calloc(post_run_mfe_res[i_res].n, sizeof(float));
    epws = (float *) calloc(post_run_mfe_res[i_res].n, sizeof(float));
    vpws = (float *) calloc(post_run_mfe_res[i_res].n, sizeof(float));
    nocc = (double *) calloc(post_run_mfe_res[i_res].n, sizeof(double));
    rocc = (double *) calloc(post_run_mfe_res[i_res].n, sizeof(double));

    post_run_ph = env.titr_ph0;
    post_run_eh = env.titr_eh0;
    if (env.titr_type == 'p') post_run_ph = env.titr_ph0 + t_point*env.titr_phd;
    else post_run_eh = env.titr_eh0 + t_point*env.titr_ehd;

    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (!(strncmp(post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID+3, "DM", 2))) {  // dummy conformer
            post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self = post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self0;
            continue;
        }

        post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_ph =  post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].H * (post_run_ph-post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].pKa) * PH2KCAL;
        post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_eh =  post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].e * (post_run_eh-post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].Em) * mev2Kcal;

        for (k=0; k<post_run_conflist.n_conf; k++) {
                mfe[j] += post_run_pairwise[post_run_mfe_res[i_res].conf[j]][post_run_conflist.conf[k].iConf] * postrun_fround3(post_run_occ_table[k][t_point]);
                epws[j] += post_run_pairwise_ele[post_run_mfe_res[i_res].conf[j]][post_run_conflist.conf[k].iConf] * postrun_fround3(post_run_occ_table[k][t_point]);
                vpws[j] += post_run_pairwise_vdw[post_run_mfe_res[i_res].conf[j]][post_run_conflist.conf[k].iConf] * postrun_fround3(post_run_occ_table[k][t_point]);
        }

        post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self = post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self0
                                                     + post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_ph
                                                     + post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_eh
                                                     + mfe[j];
     }

    /* initial Eref */
    // '0' and 'DM' go to neutral form
    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == '0' || post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            Eref[0] = post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self;
            break;
        }
    }
    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] != '0' && post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] != 'D') {
            Eref[1] = post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self;
            break;
        }
    }
    
      /* get Eref for ground and ionized state */
      /* if there is only one kind of conformer, this may be a problem */
    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == '0' || post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            if (Eref[0] > post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self) Eref[0] = post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self;
        }
        else {
            if (Eref[1] > post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self) Eref[1] = post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self;
        }
    }

    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == '0' || post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            rocc[j] = exp(-(post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self-Eref[0]) * KCAL2KT);
        }
        else {
            rocc[j] = exp(-(post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_self-Eref[1]) * KCAL2KT);
        }
    }

    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == '0' || post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == 'D') socc[0] += rocc[j];
        else socc[1] += rocc[j];
    }

    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == '0' || post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            nocc[j] = rocc[j]/socc[0];
        }
        else {
            nocc[j] = rocc[j]/socc[1];
        }
    }


    for (j=0; j<post_run_mfe_res[i_res].n; j++) {
        if (strchr(post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID, '+') || strchr(post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].uniqID, '-')) {
            E_ionize.vdw0 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_vdw0; 
            E_ionize.vdw1 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_vdw1; 
            E_ionize.ebkb += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_epol; 
            E_ionize.tors += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_tors; 
            E_ionize.dsol += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_dsolv; 
            E_ionize.offset += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_extra; 
            E_ionize.pHpK0 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_ph; 
            E_ionize.EhEm0 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_eh; 
            if (nocc[j] > 0.000001) E_ionize.TS += nocc[j] * log(nocc[j])/KCAL2KT;
            E_ionize.residues += nocc[j] * mfe[j];
            E_ionize.epws += nocc[j] * epws[j];
            E_ionize.vpws += nocc[j] * vpws[j];
            for (k=0; k<post_run_n_all; k++) {
                for (q=0; q<post_run_all_res[k].n; q++) {
                    E_ionize.mfePair[k] += nocc[j] * (post_run_pairwise[post_run_mfe_res[i_res].conf[j]][post_run_all_res[k].conf[q]] * postrun_fround3(post_run_occ_table[post_run_all_res[k].conf[q]][t_point])) ;
                    E_ionize.vdw[k] += nocc[j] * (post_run_pairwise_vdw[post_run_mfe_res[i_res].conf[j]][post_run_all_res[k].conf[q]] * postrun_fround3(post_run_occ_table[post_run_all_res[k].conf[q]][t_point])) ;
                    E_ionize.ele[k] += nocc[j] * (post_run_pairwise_ele[post_run_mfe_res[i_res].conf[j]][post_run_all_res[k].conf[q]] * postrun_fround3(post_run_occ_table[post_run_all_res[k].conf[q]][t_point])) ;
                }
                //E_ionize.vpws += E_ionize.vdw[k];
                //E_ionize.epws += E_ionize.ele[k];
            }
            if (env.mfe_pka){
                E_ionize.total = E_ionize.vdw0 + E_ionize.vdw1 + E_ionize.ebkb + E_ionize.tors + E_ionize.dsol + E_ionize.offset + E_ionize.pHpK0 + E_ionize.EhEm0 + E_ionize.TS + E_ionize.residues;
            }
            else {
                E_ionize.total = E_ionize.offset + E_ionize.tors + E_ionize.vdw0 + E_ionize.vdw1 + E_ionize.vpws + E_ionize.dsol + E_ionize.ebkb + E_ionize.epws;
            }
        }
        else {
            E_ground.vdw0 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_vdw0;                      
            E_ground.vdw1 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_vdw1;
            E_ground.ebkb += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_epol;
            E_ground.tors += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_tors;
            E_ground.dsol += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_dsolv;
            E_ground.offset += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_extra;
            E_ground.pHpK0 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_ph;
            E_ground.EhEm0 += nocc[j] * post_run_conflist.conf[post_run_mfe_res[i_res].conf[j]].E_eh;
            if (nocc[j] > 0.000001) E_ground.TS += nocc[j] * log(nocc[j])/KCAL2KT;
            E_ground.residues += nocc[j] * mfe[j];
            E_ground.epws += nocc[j] * epws[j];
            E_ground.vpws += nocc[j] * vpws[j];
            for (k=0; k<post_run_n_all; k++) {
                for (q=0; q<post_run_all_res[k].n; q++) {
                   E_ground.mfePair[k] += nocc[j] * post_run_pairwise[post_run_mfe_res[i_res].conf[j]][post_run_all_res[k].conf[q]] * postrun_fround3(post_run_occ_table[post_run_all_res[k].conf[q]][t_point]);
                   E_ground.vdw[k] += nocc[j] * (post_run_pairwise_vdw[post_run_mfe_res[i_res].conf[j]][post_run_all_res[k].conf[q]] * postrun_fround3(post_run_occ_table[post_run_all_res[k].conf[q]][t_point])) ;
                   E_ground.ele[k] += nocc[j] * (post_run_pairwise_ele[post_run_mfe_res[i_res].conf[j]][post_run_all_res[k].conf[q]] * postrun_fround3(post_run_occ_table[post_run_all_res[k].conf[q]][t_point])) ;
                }
                //E_ground.vpws += E_ground.vdw[k];
                //E_ground.epws += E_ground.ele[k];
            }
            if (env.mfe_pka){
                E_ground.total = E_ground.vdw0 + E_ground.vdw1 + E_ground.ebkb + E_ground.tors + E_ground.dsol + E_ground.offset + E_ground.pHpK0 + E_ground.EhEm0 + E_ground.TS + E_ground.residues;
            }else {
                E_ground.total = E_ground.offset + E_ground.tors + E_ground.vdw0 + E_ground.vdw1 + E_ground.vpws + E_ground.dsol + E_ground.ebkb + E_ground.epws;
            }
        } 
    }

    post_run_mfe_res[i_res].mfe.vdw0 = E_ionize.vdw0 - E_ground.vdw0;
    post_run_mfe_res[i_res].mfe.vdw1 = E_ionize.vdw1 - E_ground.vdw1;
    post_run_mfe_res[i_res].mfe.ebkb = E_ionize.ebkb - E_ground.ebkb;
    post_run_mfe_res[i_res].mfe.tors = E_ionize.tors - E_ground.tors;
    post_run_mfe_res[i_res].mfe.dsol = E_ionize.dsol - E_ground.dsol;
    post_run_mfe_res[i_res].mfe.offset = E_ionize.offset - E_ground.offset;
    post_run_mfe_res[i_res].mfe.pHpK0 = E_ionize.pHpK0 - E_ground.pHpK0;
    post_run_mfe_res[i_res].mfe.EhEm0 = E_ionize.EhEm0 - E_ground.EhEm0;
    post_run_mfe_res[i_res].mfe.TS = E_ionize.TS - E_ground.TS;
    post_run_mfe_res[i_res].mfe.residues = E_ionize.residues - E_ground.residues;
    post_run_mfe_res[i_res].mfe.total = E_ionize.total - E_ground.total;
    for (j=0; j<post_run_n_all; j++) {
        post_run_mfe_res[i_res].mfe.mfePair[j] = E_ionize.mfePair[j] - E_ground.mfePair[j];
        post_run_mfe_res[i_res].mfe.vdw[j] = E_ionize.vdw[j] - E_ground.vdw[j];
        post_run_mfe_res[i_res].mfe.ele[j] = E_ionize.ele[j] - E_ground.ele[j];
    }
    post_run_mfe_res[i_res].mfe.vpws = E_ionize.vpws - E_ground.vpws;
    post_run_mfe_res[i_res].mfe.epws = E_ionize.epws - E_ground.epws;

    free(E_ground.mfePair); free(E_ground.vdw); free(E_ground.ele);
    free(E_ionize.mfePair); free(E_ionize.vdw); free(E_ionize.ele);
    free(mfe); free(epws); free(vpws); free(nocc); free(rocc);

    return 0;
}

struct STAT postrun_fit(float a, float b)
/* initialize the simplex and set up optimization */
{  float **p;
    float y[3];
    int neval;
    struct STAT result;
    
    p = (float **) malloc(3 * sizeof(float *));
    p[0] = (float *) malloc(2 * sizeof(float));
    p[1] = (float *) malloc(2 * sizeof(float));
    p[2] = (float *) malloc(2 * sizeof(float));
    
    p[0][0] = a;      p[0][1] = b;      y[0] = postrun_score(p[0]);
    p[1][0] = a+0.01; p[1][1] = b;      y[1] = postrun_score(p[1]);
    p[2][0] = a;      p[2][1] = b+1.0;  y[2] = postrun_score(p[2]);
    
    postrun_dhill(p, y, 2, 0.0001, &postrun_score, &neval);
    
    /* DEBUG
    printf("(%8.3f%8.3f)=%8.3f exit at %d\n", p[0][0], p[0][1], y[0], neval);
    */
    
    /* postrun_fit this eq. y = exp(a(x-b))/(1+exp(a(x-b))) */
    
    result.a = p[0][0];
    result.b = p[0][1];
    result.chi2 = y[0];
    
    free(p[0]);
    free(p[1]);
    free(p[2]);
    free(p);
    
    return result;
}


float postrun_score(float v[])
{  float S = 0.0;
    float yt;
    float T;
    int i;
    
    for (i=0; i<post_run_Nx; i++) {
        T = exp(v[0]*(post_run_xp[i]-v[1]));
        yt = T/(1.0+T);
        S += (yt-post_run_yp[i])*(yt-post_run_yp[i]);
    }
    
    return S;
}

                 
void postrun_dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk)
{
    float postrun_dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac);
    int i, ihi, ilo, inhi,j, mpts = ndim+1;
    float rtol,sum,swap,ysave,ytry,*psum;
    
    
    psum = (float *)malloc(ndim * sizeof(float));
    *nfunk = 0;
    GET_PSUM
    
    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1,0):(inhi = 0,1);
        for(i=0; i<mpts;i++) {
            if(y[i] <= y[ilo]) ilo=i;
            if(y[i] > y[ihi]) {
                inhi = ihi;
                ihi  = i;
            }
            else if (y[i] > y[inhi] && i != ihi) inhi = i;
        }
        
        /* DEBUG
        for (i=0; i<mpts; i++) {
            printf("%8.3f at (%8.3f, %8.3f)\n", y[i], p[i][0], p[i][1]);
        }
        printf("\n");
        */
        
        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        
        if(rtol < ftol || *nfunk >= NMAX) {
            SWAP(y[0],y[ilo])
            for (i=0; i<ndim; i++) SWAP(p[0][i],p[ilo][i])
                break;
        }
        
        *nfunk += 2;
        
        ytry = postrun_dhtry(p,y,psum,ndim,funk,ihi,-1.0);
        
        if (ytry <= y[ilo]) ytry=postrun_dhtry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry=postrun_dhtry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0; i<mpts; i++) {
                    if (i !=ilo) {
                        for (j=0; j<ndim; j++)
                            p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);
                        y[i] = (*funk)(psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        }
        else --(*nfunk);
    }
    free(psum);
}

float postrun_dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac)
{ 
    int j;
    float fac1, fac2, ytry, *ptry;
    
    ptry = (float *)malloc(ndim*sizeof(float));
    fac1 = (1.0-fac)/ndim;
    fac2 = fac1 - fac;
    for (j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry = (*funk)(ptry);   /* evaluate the function at the trial point */
    if (ytry < y[ihi]) {    /* if it's better than the highest, then replace the highest */
        y[ihi] = ytry;
        for (j=0; j<ndim; j++) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    free(ptry);
    return ytry;
}


int postrun_load_occupancy() 
{
   FILE *post_run_fp;
   int conf_counter, titr_step_counter, j, titr_step_counter_2;
   char *tok, *tok_2;
   char sbuff[MAXCHAR_LINE];
   
   post_run_fp = fopen(OCC_TABLE, "r");
   
   /* skip 1 lines */
   for (j=0; j<1; j++) fgets(sbuff, sizeof(sbuff), post_run_fp);

   conf_counter = 0;
   while (fgets(sbuff, sizeof(sbuff), post_run_fp)) {
      //printf("%s\n", sbuff);
      //if (strstr(sbuff, "DM") != NULL){
          //printf("%s\n", sbuff);
          //tok_2 = strtok (sbuff," \0");
          //printf("tok:  %d\n", tok);
          //titr_step_counter_2 = 0;
          //while (tok_2 != NULL && titr_step_counter_2 < env.titr_steps)
          //{
             //tok_2 = strtok (NULL, " \0");
             //most_occ_table[conf_counter][titr_step_counter_2] = atof(tok_2);
             //++titr_step_counter_2;
          //}
      //}

      tok = strtok (sbuff," \0");
      titr_step_counter = 0;
      titr_step_counter_2 = 0;
      while (tok != NULL && titr_step_counter < env.titr_steps)
      {
         tok = strtok (NULL, " \0");
         //printf("%d %f\n",conf_counter,atof(tok));
         //fflush(stdout);
         //if (strstr(sbuff, "DM") == NULL){
            //most_occ_table[conf_counter][titr_step_counter_2] = atof(tok);
            //++titr_step_counter_2;
         //}
         post_run_occ_table[conf_counter][titr_step_counter] = atof(tok);
         //printf("tok: %5.2f\n", atof(tok));
         ++titr_step_counter;
      }
      ++conf_counter;
   }

   conf_counter = 0;
   

   
   
   fclose(post_run_fp);
   return 0;
}

