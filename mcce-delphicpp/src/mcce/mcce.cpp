#include <stdio.h>
#include "mcce.h"
#include <time.h>

void welcome();
#define LEN 150

int main(int argc, char *argv[])
{
   /* Welcome */
   welcome();


   /* Do step 0, initialization */
   db_open();
   printf("Step 0. Initialize enviroment\n"); fflush(stdout);
   if (init()) {
      db_close();
      printf("Help message: double check file \"run.prm\" in current directory.\n");
      return USERERR;
   }
   else printf("Step 0 Done.\n\n");


   /* Do step 1, premcce */
   if (env.do_premcce) {
      printf("Step 1. Test and format structral file\n"); fflush(stdout);
      if (premcce()) {db_close(); return USERERR;}
      else printf("Step 1 Done.\n\n");
   }
   else printf("Not doing \"Step 1. Test and format structral file\"\n\n");


   /* Do step 2. rotamers */
   if (env.do_rotamers) {
      printf("Step 2. Make multi side chain conformers\n"); fflush(stdout);
      if (rotamers()) {
		db_close(); return USERERR;
      }
      else printf("Step 2 Done.\n\n");
   }
   else printf("Not doing \"Step 2. Make multi side chain conformers\"\n\n");

   /* Do step 3. energies */
   if (env.do_energies) {
      printf("Step 3. Compute energy lookup table\n"); fflush(stdout);
      if (energies()) {db_close(); return USERERR;}
      else printf("Step 3 Done.\n\n");
   }
   else printf("Not doing \"Step 3. Compute energy lookup table\"\n\n");

   /* Do step 4. Monte Carlo */
   if (env.do_monte) {
      if (!env.monte_adv_opt) {
      printf("Step 4. Standard Monte Carlo Sampling\n"); fflush(stdout);
      if (monte()) {db_close(); return USERERR;}
           else printf("Step 4 Done.\n\n");
      }
       else {
         printf("Step 4. Yifan Monte Carlo Sampling\n"); fflush(stdout);
         if (monte2()) {db_close(); return USERERR;}
         else printf("Step 4 Done.\n\n");
      }
     /* if (env.monte_adv_opt) {               // if monte_adv_opt is ture, run Yifan Monte Carlo         
         printf("Step 4. Yifan Monte Carlo Sampling\n"); fflush(stdout);
         if (monte2()) {db_close(); return USERERR;}
         else printf("Step 4 Done.\n\n");
      }      else if (env.monte_ms) {               // If monte_ms is true, run MicroState Monte Carlo Sampling
         printf("Step 4. MicroState Monte Carlo Sampling\n"); fflush(stdout);
         if (monte3()) {db_close(); return USERERR;}
         else printf("Step 4 Done.\n\n");
      }      else{                                  // If all flags are false, run standard Monte Carlo        
         printf("Step 4. Standard Monte Carlo Sampling\n"); fflush(stdout);
         if (monte()) {db_close(); return USERERR;}
         else printf("Step 4 Done.\n\n");
      } */
   }
   else printf("Not doing \"Step 4. Monte Carlo Sampling\"\n\n");
   
   if (env.do_postmcce){
      printf("Step 5. Post Monte Carlo\n"); fflush(stdout);
      if (postrun()) {db_close(); return USERERR;}
      else printf("Step 5 Done.\n\n");
   }
   else printf("Not doing \"Step 5. Post Monte Carlo\"\n\n");

   /* Do step 6. hydrogen bond analysis---Cai */
   if (env.do_analysis){
      printf("Step 6. Hydrogen Bond Network Analysis\n"); fflush(stdout);
      if (analysis()) {db_close(); return USERERR;}
      else printf("Step 6 Done.\n\n");
   }
   else printf("Not doing \"Step 6. Hydrogen Bond Network Analysis\"\n\n");


   db_close();
   return 0;
}

void welcome()
{ 
   printf(" _________________________MCCE 3.5____________________________\n");
   printf("|	   						      |\n");
   printf("|    MCCE (Multi-Conformation Continuum Electrostatics)       |\n");
   printf("|	is a program developed at Marilyn Gunner's lab.	      |\n");
   printf("|	MCCE is a biophysics simulation program combining     |\n");
   printf("|	continuum electrostatics and molecular mechanics.     |\n");
   printf("|	In this program, the protein side chain motions are   |\n");
   printf("|	simulated explicitly while the dielectric effect of   |\n");
   printf("|	solvent and bulk protein material is modeled by	      |\n");
   printf("|	continuum electrostatics.			      |\n");
   printf("|	MCCE can calculate residue pka, cofactor Em and       |\n");
   printf("|	protein PI in protein-solvent systems, and more:      |\n");
   printf("|							      |\n");
   printf("|	 - Protein structural responses to changes in charge  |\n");
   printf("|	 - Changes in charge state of ionizable residues due  |\n");
   printf("|	   to structural changes in the protein		      |\n");
   printf("|	 - The structural and ionization changes caused by    |\n");
   printf("|	   changes in solution pH or Eh			      |\n");
   printf("|	 - Find the location and stoichiometry of proton      |\n");
   printf("|	   transfers coupled to electron transfer	      |\n");
   printf("|	 - Make side chain rotomer packing predictions as     |\n");
   printf("|	   a function of pH				      |\n");
   printf("|							      |\n");
   printf("|	For questions and help, visit                         |\n");
   printf("|		https://sites.google.com/site/mccewiki/home   |\n");
   printf("|		  October 2015, by MCCE Development Team      |\n");
   printf("|							      |\n");
   printf("|_________________________        ____________________________|\n");
   printf("                          MCCE 3.5                          \n\n");
   printf("Latest Updates:                                              \n");
   printf("July 2016, Yifan's monte carlo no longer needs step3_out.pdb directory (fixed)\n");
   printf("July 2016, Removed PASCAL COMTE GENETIC ALGORITHM from this version\n");
   printf("Sept 2016, Added charged ligands to sum_crg.out.\n");
   printf("Oct  2016, Changed Yifan MC pK.out output format to match mfe format.\n");
   printf("Oct  2016, Added step 5 that creates pK.out and sum_crg.out.\n");
   printf("Jan  2017, Added DelPhi potential map.\n\n");
   printf("Mar  2017, Added DelPhi dielectric map.\n\n");
   printf("July 2017, Added monte3.c to output MicroState at step 4.\n\n");
   printf("July 2017, Added step 6 to do Hydrogen Bond Network analysis.\n\n");
   fflush(stdout);
	
	// Added by Jessica on Nov. 2015
   char buf[LEN];
   time_t curtime;
   struct tm *loc_time;

   //Getting current time of system
   curtime = time (NULL);

   // Converting current time to local time
   loc_time = localtime (&curtime);

   // Displaying date and time in standard format
   printf("%s", asctime (loc_time));
   return;
}
