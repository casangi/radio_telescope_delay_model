!
!  Include file a_out.i
!
!   Include file for all ALMA outputs from the delay computations.
!
      !DEC$ ATTRIBUTES ALIAS:'alma_out' :: alma_out
      Real*8 delay_vac, rate_vac, dry_atm1(2), wet_atm1(2),             &
     &       dry_atm2(2), wet_atm2(2), elev1(2), elev2(2), u_v_w(3),    &
     &       az1(2), az2(2), tph1(3), tph2(3)
!
      COMMON /alma_out/ delay_vac, rate_vac, dry_atm1, wet_atm1,        &
     &       dry_atm2, wet_atm2, elev1, elev2, u_v_w,                   &
     &       az1, az2, tph1, tph2
