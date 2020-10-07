/*
 *  Include file a_out.h
 *
 *   Include file for all ALMA outputs from the delay computations.
 *
*/
#define alma_out 

      extern struct  { 
       double  delay_vac, rate_vac, dry_atm1[2], wet_atm1[2], 
            dry_atm2[2], wet_atm2[2], elev1[2], elev2[2], u_v_w[3], 
            az1[2], az2[2];
        }  alma_out_;
