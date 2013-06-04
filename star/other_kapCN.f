      subroutine kapCN_get_Type1( &
            id, k, handle, zbar, X, Zbase, log10_rho, log10_T,  &
            species, chem_id, net_iso, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         
         use kap_lib, only: kap_get_Type1

         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ion charge
         real(dp), intent(in) :: X ! the hydrogen mass fraction
         real(dp), intent(in) :: Zbase ! the metallicity
         real(dp), intent(in) :: log10_rho ! the density
         real(dp), intent(in) :: log10_T ! the temperature
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         double precision, intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         
         ! OUTPUT
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.

         call kap_get_Type1( handle, zbar, X, Zbase, log10_rho, log10_T, &
                             lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                             kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         
      end subroutine kapCN_get_Type1



      subroutine kapCN_get_Type2( &
            id, k, handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
            log10_rho, log10_T, species, chem_id, net_iso, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

         use chem_def
         use kap_lib, only: kap_get_Type2

         implicit none
         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ion charge
         real(dp), intent(in) :: X, Z, Zbase, XC, XN, XO, XNe ! composition    
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         double precision, intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons

         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         ! OUTPUT
         real(dp), intent(out) :: frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
                  
         logical, parameter :: debug=.false.
         real(dp) :: log10_R

         !for converting sp <-> dp
         real(sp) :: Z_sp, X_sp, fC, fN, logRho_sp, logT_sp
         real(sp) :: kap_sp, dlnkap_dlnRho_sp, dlnkap_dlnT_sp

         log10_R = log10_Rho - 3*log10_T + 18d0

         if(log10_T > 4.05d0 .or. log10_R > 1d0)then !use MESA Type2

            call kap_get_Type2( &
            handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
            log10_rho, log10_T,&
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

         else  !use C,N-enhanced low-T opacities:

            Z_sp=real(Zbase,sp); X_sp=real(X,sp); 
            logRho_sp=real(log10_Rho,sp); logT_sp=real(log10_T,sp)

            !Lederer & Aringer tables use Lodders 2003 
            fC=real(XC/(Zbase*L03_element_zfrac(e_C)),sp)
            fN=real(XN/(Zbase*L03_element_zfrac(e_N)),sp)

            call kapCN_get(Z_sp,X_sp,fC,fN,logRho_sp,logT_sp,kap_sp, &
                           dlnkap_dlnRho_sp, dlnkap_dlnT_sp, ierr)

            if(debug.or.kap_sp==1.0)then
               write(*,*) 'logRho=',log10_Rho
               write(*,*) 'logR=',log10_Rho - 3*log10_T +18d0
               write(*,*) 'logT=',log10_T
               write(*,*) 'XC=', XC
               write(*,*) 'XN=', XN
               write(*,*) 'Zbase=', Zbase
               write(*,*) 'Zsp=', Z_sp
               write(*,*) 'Xsp=', X_sp
               write(*,*) 'fC=', fC
               write(*,*) 'fN=', fN
               write(*,*) 'kap_sp=', kap_sp
               write(*,*)
            endif

            kap = real(kap_sp,dp)
            dln_kap_dlnRho = real(dlnkap_dlnRho_sp,dp)
            dln_kap_dlnT = real(dlnkap_dlnT_sp,dp)
            
         endif

         frac_Type2 = 1

      end subroutine kapCN_get_Type2
