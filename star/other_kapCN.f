     subroutine kapCN_get_Type1( &
            id, k, handle, zbar, X, Zbase, log10_rho, log10_T,  &
            species, chem_id, net_iso, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         
         use kap_lib, only: kap_get_Type1
         
         implicit none
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
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, use_Zbase_for_Type1, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

         !use chem_def
         use kap_lib, only: kap_get_Type2

         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ion charge
         real(dp), intent(in) :: X, Z, Zbase, XC, XN, XO, XNe ! composition    
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         logical, intent(in) :: use_Zbase_for_Type1

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
         logical :: try_kapCN
         real(dp) :: log10_R

         !for converting sp <-> dp
         real(sp) :: Z_sp, X_sp, fC, fN, logRho_sp, logT_sp
         real(sp) :: kap_sp, dlnkap_dlnRho_sp, dlnkap_dlnT_sp

         ierr=99
         log10_R = log10_Rho - 3*log10_T + 18d0
         
         try_kapCN = X_C_init > 0d0 .and. X_N_init > 0d0 .and. log10_T < 4d0

         if(try_kapCN)then

            Z_sp=real(Zbase,sp); X_sp=real(X,sp); 
            logRho_sp=real(log10_Rho,sp); logT_sp=real(log10_T,sp)

            !Lederer & Aringer tables use Lodders (2003) abundance scale
            fC=real(XC/X_C_init,sp)
            fN=real(XN/X_N_init,sp)

            call kapCN_get(Z_sp,X_sp,fC,fN,logRho_sp,logT_sp,kap_sp, &
                           dlnkap_dlnRho_sp, dlnkap_dlnT_sp, ierr)

            if(debug.or.kap_sp==1.0)then
               write(*,*) 'logRho=',log10_Rho
               write(*,*) 'logR=',log10_R
               write(*,*) 'logT=',log10_T
               write(*,*) 'XC=', XC
               write(*,*) 'XN=', XN
               write(*,*) 'Zbase=', Zbase
               write(*,*) 'Zsp=', Z_sp
               write(*,*) 'Xsp=', X_sp
               write(*,*) 'fC=', fC
               write(*,*) 'fN=', fN
               write(*,*) 'kap_sp=', kap_sp
               write(*,*) 'X_C_init=', X_C_init
               write(*,*) 'X_N_init=', X_N_init
               write(*,*) 'ierr=', ierr
               write(*,*)
            endif

            if(ierr==0)then
               kap = real(kap_sp,dp)
               dln_kap_dlnRho = real(dlnkap_dlnRho_sp,dp)
               dln_kap_dlnT = real(dlnkap_dlnT_sp,dp)
            endif

         endif

         if(.not.try_kapCN .or. ierr/=0)then
            call kap_get_Type2( &
            handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
            log10_rho, log10_T,&
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, use_Zbase_for_type1, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         endif

         frac_Type2 = 1d0

      end subroutine kapCN_get_Type2
