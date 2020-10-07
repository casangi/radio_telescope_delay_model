#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*  #include "a_out.h"  */

      extern struct  { 
       double  delay_vac, rate_vac, dry_atm1[2], wet_atm1[2],
            dry_atm2[2], wet_atm2[2], elev1[2], elev2[2], u_v_w[3],
            az1[2], az2[2];
        }  alma_out_;
  

int almaout_ (int iref, int iremot, int isourc, double xsec,
             int jtag[5], double delvac, double dryatm1, double dryatm2,
             double wetatm1, double wetatm2, char statn1[8], 
             char statn2[8], char sourc20[20] )
{
/*       int iref, iremot, isourc, jtag[5];
         double xsec, delvac; */

         printf(" C: iref,iremot,isourc: %d, %d, %d\n", iref, iremot, isourc); 
         printf(" C: Jtag, xsec: %d %d %d %d %d %f\n", jtag[0],jtag[1],jtag[2],jtag[3],jtag[4], xsec); 
         printf(" C: delvac: %25.20f\n", delvac); 
         printf(" C: dryatm1, dryatm2: %25.20f %25.20f\n", dryatm1, dryatm2); 
         printf(" C: wetatm1, wetatm2: %25.20f %25.20f\n", wetatm1, wetatm2); 
         printf(" C: statn1: %8s\n", statn1);
         printf(" C: statn2: %8s\n", statn2);
         printf(" C: source: %20s\n", sourc20);
/*       printf("delay_vac: %20.10f\n", delay_vac); */


         return 0;
}
