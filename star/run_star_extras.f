!     ***********************************************************************
!     
!     Copyright (C) 2010  Bill Paxton
!     
!     this file is part of mesa.
!     
!     mesa is free software; you can redistribute it and/or modify
!     it under the terms of the gnu general library public license as published
!     by the free software foundation; either version 2 of the license, or
!     (at your option) any later version.
!     
!     mesa is distributed in the hope that it will be useful, 
!     but without any warranty; without even the implied warranty of
!     merchantability or fitness for a particular purpose.  see the
!     gnu library general public license for more details.
!     
!     you should have received a copy of the gnu library general public license
!     along with this software; if not, write to the free software
!     foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!     
!     ***********************************************************************
      
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      use kapCN
      
      implicit none

      double precision :: X_C_init, X_N_init
      
! these routines are called by the standard run_star check_model
      contains

      !include 'other_kapCN.f90'
      
      subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

! this is the place to set any procedure pointers you want to change
! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
      s% other_kap_get_Type1 => kapCN_get_Type1
!s% other_kap_get_Type2 => kapCN_get_Type2
      s% other_kap_get_Type2 => null_other_kap_get_Type2
      
! Uncomment these lines if you wish to use the functions in this file,
! otherwise we use a null_ version which does nothing.
      s% extras_startup => extras_startup

! Once you have set the function pointers you want,
! then uncomment this (or set it in your star_job inlist)
! to disable the printed warning message,
      s% job% warn_run_star_extras =.false.       

      ! the following values are used in the kapCN routines to set the ratios fC and fN
      X_C_init = s% x_ctrl(1)
      X_N_init = s% x_ctrl(2)
     
      end subroutine extras_controls
      
! None of the following functions are called unless you set their
! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
      use chem_def, only: ic12, in14
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_startup = 0
      if (.not. restart) then
         call alloc_extra_info(s)
      else                      ! it is a restart
         call unpack_extra_info(s)
      end if

      call kapCN_init(ierr)
      
      end function extras_startup

      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info


      
      subroutine null_other_kap_get_Type2( &
            id, k, handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
            log10_rho, log10_T, species, chem_id, net_iso, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, use_Zbase_for_Type1, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
 
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
                  
         frac_Type2 = 0; kap = 0; dln_kap_dlnRho = 0; dln_kap_dlnT = 0
         
         write(*,*) 'no implementation for other_kap_get_Type2'
         ierr = -1

      end subroutine null_other_kap_get_Type2

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


      
      end module run_star_extras
      
