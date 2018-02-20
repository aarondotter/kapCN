! Aaron Dotter -- May 2013
! kapCN is a module that reads in and interpolates the low-T
! opacities by M.T. Lederer & B. Aringer, 2009, A&A, 494, 403
! see doc/ReadMe for details of the data files; it uses MESA
! modules for interpolation and a few other things

      module kapCN

      use const_def, only: sp, dp, mesa_data_dir
      use crlibm_lib, only: exp10_cr
      use utils_lib, only: alloc_iounit, free_iounit
      use num_lib, only: binary_search_sg
      
      implicit none
     
      private

      public :: kapCN_init, kapCN_shutdown, kapCN_interp, kapCN_get

      integer, parameter :: num_kapCN_Xs = 3
      integer, parameter :: num_kapCN_Zs = 14
      integer, parameter :: num_kapCN_fCs = 7
      integer, parameter :: num_kapCN_fNs = 3
      integer, parameter :: num_falpha = 1 !not useful yet
      integer, parameter :: num_tbl = num_kapCN_Xs*num_kapCN_fCs*num_kapCN_fNs ! 63
      integer, parameter :: num_logT = 18
      integer, parameter :: num_logR = 17
      integer, parameter :: tbl_size = num_logR*num_logT
     
      real(sp), target :: kapCN_Z(num_kapCN_Zs)
      real(sp), target :: kapCN_fN(num_kapCN_fNs,num_kapCN_Zs)
      real(sp), target :: kapCN_fC(num_kapCN_fCs,num_kapCN_Zs)
      real(sp), target :: kapCN_X(num_kapCN_Xs)  = [ 0.5, 0.7, 0.8 ]
      real(sp), target :: kapCN_logT(num_logT), kapCN_logR(num_logR)
      real(sp) :: kapCN_min_logR, kapCN_max_logR, kapCN_min_logT, kapCN_max_logT

      !Z mass fractions of C and N, unused
      !real(sp) :: zC=0.1644, zN=0.0532 

      logical, parameter :: debug = .false., use_cache = .true.
      logical :: kapCN_is_initialized = .false.

      character(len=32), parameter :: kapCN_param_file = 'kR_Z_fCN.dat'
      character(len=256) :: kapCN_data_dir

      !for 2-D interpolation
      integer :: ibcx=0, ibcy=0, ilinx=1, iliny=1
      integer :: ict(6) = [ 1, 1, 1, 0, 0, 0 ]
      real(sp), parameter :: bc(num_logT) = 0.0

      !stores a table for one {Z,X,fC,fN}
      type KapCN_Table
         real(sp) :: X
         real(sp) :: Z, Zbase
         real(sp) :: fC
         real(sp) :: fN
         real(sp) :: falpha
         real(sp), pointer :: kap(:)
      end type KapCN_Table
      
      !stores a set of tables for one Z
      type KapCN_Set
         character(len=12) :: filename
         logical :: not_loaded_yet
         real(sp) :: Zbase
         real(sp) :: fC(num_kapCN_fCs)
         real(sp) :: fN(num_kapCN_fNs)
         type(KapCN_Table) :: t(num_tbl)
      end type KapCN_Set
      type(KapCN_Set), target :: kCN(num_kapCN_Zs)
      
      contains

      subroutine kapCN_init(ierr)
      integer, intent(out) :: ierr
      integer :: iZ
      type(kapCN_set), pointer :: k

      if(len_trim(mesa_data_dir)==0)then
         write(0,*) ' Need to set data_dir with call to const_init!'
         ierr=-1
         return
      else
         kapCN_data_dir = mesa_data_dir
      endif
      
      do iZ=1,num_kapCN_Zs
         k => kCN(iZ)
         k% not_loaded_yet = .true.
         k% Zbase = kapCN_Z(iZ)
         k% fC = kapCN_fC(:,iZ)
         k% fN = kapCN_fN(:,iZ)
      enddo

      call read_kapCN_tables(ierr)

      if(ierr==0) kapCN_is_initialized = .true.
      
      end subroutine kapCN_init


      subroutine read_kapCN_tables(ierr)
      integer, intent(out) :: ierr
      character(len=256) :: filename
      character(len=4) :: prefix, string_Z(num_kapCN_Zs), suffix
      character(len=6) :: tmp_Z
      integer :: i, io
      
      io=alloc_iounit(ierr)
      if(ierr/=0) write(*,*) 'problem allocating io unit number'

      prefix = 'kR_Z'
      suffix = '.dat'

      string_Z = ''
      
      filename = trim(kapCN_data_dir) // '/kap_data/' // trim(kapCN_param_file)

      if(debug) write(*,*) '   parameter file = ', trim(filename)

      open(io,file=trim(filename))
      read(io,*) !skip first line
      do i=num_kapCN_Zs,1,-1 !read backwards into array so Z is increasing
         read(io,'(f7.5,11f9.1)') kapCN_Z(i), kapCN_fC(:,i), kapCN_fN(:,i)
      enddo
      close(io)
      call free_iounit(io)

      do i=1,num_kapCN_Zs
         write(tmp_Z,'(1p,1e6.1e1)') kapCN_Z(i)
         string_Z(i) = tmp_Z(1:1) // tmp_Z(4:6)
      enddo

      kCN(:)% filename = prefix // string_Z // suffix
      
      do i=1,num_kapCN_Zs
         call read_one_table(i,ierr)
      enddo

      kapCN_min_logR = minval(kapCN_logR)
      kapCN_max_logR = maxval(kapCN_logR)
      kapCN_min_logT = minval(kapCN_logT)
      kapCN_max_logT = maxval(kapCN_logT)

      end subroutine read_kapCN_tables


      subroutine read_one_table(iZ,ierr)
      use interp_2d_lib_sg, only: interp_mkbicub_sg
      integer, intent(in) :: iZ
      integer, intent(out) :: ierr
      character(len=256) :: ascii_file, cache_file
      character(len=10) :: my_Z
      integer :: i, io, j
      real(sp) :: c_div_o, y
      real(sp) :: table(num_logR,num_logT)
      real(sp), pointer :: logR(:), logT(:)
      logical :: have_cache
      type(kapCN_set), pointer :: k

      ierr=0

      k => kCN(iZ)
      logR => kapCN_logR
      logT => kapCN_logT

      ! data files
      ascii_file = trim(kapCN_data_dir) // '/kap_data/' // trim(k% filename)
      cache_file = trim(kapCN_data_dir) // '/kap_data/cache/' // trim(k% filename(1:9)) // 'bin'

      if(debug) write(*,*) ' ascii file = ', trim(ascii_file)
      if(debug) write(*,*) ' cache file = ', trim(cache_file)

      !check for cached bin file; if it exists, read it and exit
      if(use_cache)then
         inquire(file=trim(cache_file),exist=have_cache)
         if(have_cache)then
            !read cache
            open(io,file=trim(cache_file),form='unformatted',iostat=ierr)
            if(ierr/=0) write(*,*) &
                ' read_one_table: problem opening cache file for read'
            read(io) logR, logT
            read(io) k% t(:)% X
            read(io) k% t(:)% Z
            read(io) k% t(:)% fC
            read(io) k% t(:)% fN
            read(io) k% t(:)% falpha
            do i=1,num_tbl
               allocate(k% t(i)% kap(4*tbl_size))
               read(io) k% t(i)% kap
            enddo
            !table is now fully loaded
            k% not_loaded_yet = .false.
            return
         endif
      endif

      !no cache file, so read ascii file
      io=alloc_iounit(ierr)
      open(io,file=trim(ascii_file),iostat=ierr)
      if(ierr/=0) write(*,*) &
           ' read_one_table: problem opening ascii file for read'
      do i=1,4
         read(io,*) !skip header
      enddo
      read(io,'(6x,a10)') my_Z
      if(debug) write(*,*) ' has Z = ', my_Z, ' Zbase = ', k% Zbase

      do i=1,14
         read(io,*)
      enddo
      
      do i=1,num_tbl
         read(io,*) j, k% t(i)% x, y, k% t(i)% z, c_div_o, &
                    k% t(i)% fc, k% t(i)% fn, k% t(i)% falpha
      enddo

      read(io,*) !last line of #'s

      do i=1,num_tbl
         do j=1,14
            read(io,*) !skip individual headers
         enddo
         read(io,'(5x,17f7.3)') logR
         do j=1,num_logT
            read(io,'(f5.3,17f7.3)') logT(j), table(:,j)
         enddo
         allocate(k% t(i)% kap(4*tbl_size))
         k% t(i)% kap(1:4*tbl_size:4) = reshape(table,[tbl_size])
      enddo

      close(io)

      !create interpolants for tables
      do i=1,num_tbl
         call interp_mkbicub_sg( logR, num_logR, &
                                 logT, num_logT, &
                                 k% t(i)% kap, num_logR, &
                                 ibcx, bc, ibcx, bc, &
                                 ibcy, bc, ibcy, bc, &
                                 ilinx, iliny, ierr)
      end do

      !table is now fully loaded
      k% not_loaded_yet = .false.
      
      !write cache file
      open(io,file=trim(cache_file),form='unformatted',iostat=ierr)
      if(ierr/=0) write(*,*) &
             ' read_one_table: problem opening cache file for write'
      write(io) logR, logT
      write(io) k% t(:)% x
      write(io) k% t(:)% z
      write(io) k% t(:)% fC
      write(io) k% t(:)% fN
      write(io) k% t(:)% falpha
      do i=1,num_tbl
         write(io) k% t(i)% kap
      enddo
      close(io)
      call free_iounit(io)

      end subroutine read_one_table


      !this one is specifically for use with MESA
      subroutine kapCN_get(Z,X,fC,fN,logRho,logT,kappa, &
                           dlnkap_dlnRho,dlnkap_dlnT,ierr)
      real(sp), intent(in) :: Z, X, fC, fN, logRho, logT
      real(sp), intent(out) :: kappa
      real(sp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
      integer, intent(out) :: ierr
      real(sp) :: logR, result(3)

      ierr = 0
      logR = logRho - 3.0*logT + 18.0

      call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      if(ierr==0)then
         kappa = real(exp10_cr(real(result(1),kind=dp)),kind=sp)
         dlnkap_dlnRho = result(2)
         dlnkap_dlnT = result(3) - 3*result(2)
      else
         kappa = -1d0
         dlnkap_dlnRho = 0d0
         dlnkap_dlnT = 0d0
      endif
      
      end subroutine kapCN_get


      subroutine kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      !result=(logKappa, dlogKappa/dlogR, dlogKappa/dlogT)
      use interp_1d_def, only: pm_work_size
      use interp_1d_lib_sg, only: interpolate_vector_sg, interp_pm_sg
      real(sp), intent(in) :: Z, X, fC, fN, logR, logT
      real(sp), intent(out) :: result(3)
      integer, intent(out) :: ierr
      real(sp), pointer :: Z_ary(:), work1(:)
      integer, parameter :: nZ = 4
      integer :: i,iZ
      real(sp) :: my_Z, res(3,nZ), x_new(1), v_new(1)
      character(len=32) :: junk
      
      result=0.0; iZ=0; ierr=0

      if(.not.kapCN_is_initialized)then
         write(*,*) ' kapCN is not initialized; call kapCN_init()'
         ierr=-99
         return
      endif

      if(outside_R_and_T_bounds(logR,logT))then
         !write(*,*) 'kapCN_interp: logR, logT outside of table bounds'
         ierr=-1
         return
      endif

      Z_ary => kapCN_Z
      my_Z = max(min(Z,Z_ary(num_kapCN_Zs)),Z_ary(1))
      iZ=binary_search_sg(num_kapCN_Zs,Z_ary, iZ, Z)
      iZ = max(nz/2,min(iZ,num_kapCN_Zs-nz/2))
      
      !check to see if exact match in Z, then just need one call
      if(Z==Z_ary(iZ))then
         call kapCN_interp_fixedZ(iZ,X,fC,fN,logR,logT,result)
         return
      endif

      !else do the full Z interpolation
      do i=1,nZ
         call kapCN_interp_fixedZ(iZ+i-2,X,fC,fN,logR,logT,res(:,i))
      enddo

      nullify(work1)
      allocate(work1(nZ*pm_work_size))
      x_new(1)=log10(my_Z)

      do i=1,3
         call interpolate_vector_sg(nZ,log10(Z_ary(iZ-1:iZ+2)), 1, &
                                    x_new, res(i,:), v_new, &
                                    interp_pm_sg, pm_work_size, &
                                    work1, junk, ierr )
         result(i)=v_new(1)
      enddo
      
      deallocate(work1)

      end subroutine kapCN_interp


      logical function outside_R_and_T_bounds(logR,logT)
      real(sp), intent(in) :: logR, logT     
      outside_R_and_T_bounds = &
         logR < kapCN_min_logR .or. logR > kapCN_max_logR .or. &
         logT < kapCN_min_logT .or. logT > kapCN_max_logT
      end function outside_R_and_T_bounds


      subroutine kapCN_interp_fixedZ(iZ,X,fC,fN,logR,logT,result)
      !using simple quadratic interpolation in each of X, fC, fN
      integer, intent(in) :: iZ
      real(sp), intent(in) :: X, fC, fN, logR, logT
      real(sp), intent(out) :: result(3)
      integer, parameter :: my_num_kapCN_fCs = 3, f1=num_kapCN_fCs*num_kapCN_Xs, f2=num_kapCN_fCs
      real(sp) :: my_X, my_fC, my_fN, wX(num_kapCN_Xs), wfN(num_kapCN_fNs), wfC(my_num_kapCN_fCs)
      real(sp) :: res(3)
      integer :: iX, ifC, ifN, i,j,tbl,n,nlo,nhi,ierr
      real(sp), pointer :: X_ary(:), fC_ary(:), fN_ary(:)

      result = 0.0; iX=0; ifC=0; ifN=0
      X_ary => kapCN_X; fC_ary => kapCN_fC(:,iZ); fN_ary => kapCN_fN(:,iZ)

      !limit input quantities to lie within tabulated range
      my_X = max(min(X, X_ary(num_kapCN_Xs)),X_ary(1))
      my_fC = max(min(fC,fC_ary(num_kapCN_fCs)),fC_ary(1))
      my_fN = max(min(fN,fN_ary(num_kapCN_fNs)),fN_ary(1))

      !locate inputs in arrays
      iX = binary_search_sg(num_kapCN_Xs, X_ary, iX, my_X)
      ifC= binary_search_sg(num_kapCN_fCs, fC_ary, ifC, my_fC)
      ifN= binary_search_sg(num_kapCN_fNs, fN_ary, ifN, my_fN)

      !make sure 1 < ifC < num_kapCN_fCs
      ifC = max(2,min(ifC,num_kapCN_fCs-1))

      !interpolation coefficients in X
      call interp(X_ary,wX,my_X,num_kapCN_Xs)

      !interpolation coefficients in fN
      call interp(fN_ary(:),wfN,my_fN,num_kapCN_fNs)

      !interpolation coefficients in fC
      call interp(fC_ary(ifC-1:ifC+1),wfC,my_fC,my_num_kapCN_fCs)

      do i=1,num_kapCN_fNs
         do j=1,num_kapCN_Xs
            do n=1,my_num_kapCN_fCs
               res=0.0
               tbl = f1*(i-1) + f2*(j-1) + (ifC+n-2)
               nlo = 3*(n-1)+1
               nhi = nlo+2
               call kapCN_interp_RT(iZ,tbl,logR,logT,res,ierr)
               result = result + wfN(i)*wX(j)*wfC(n)*res
            enddo
         enddo
      enddo
      
      if(debug) write(*,'(1p9e12.4)') my_X, my_fC, my_fN, result
      
      contains

      subroutine interp(a,b,x,n)
! {a} are the tabulated values for use in interpolation
! {b} are coefficients of the interpolating polynomial
!  x  is the abscissa to be interpolated
!  n  is the number of points to be used, interpolating polynomial
!     has order n-1 
      integer, intent(in) :: n
      real(sp), intent(in) :: a(n), x
      real(sp), intent(out) :: b(n)
      integer :: i,j
      do i=1,n
         b(i)=1.0
         do j=1,n
            if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
      end subroutine interp

      end subroutine kapCN_interp_fixedZ


      subroutine kapCN_interp_RT(iZ,tbl,logR,logT,result,ierr)
      use interp_2d_lib_sg, only: interp_evbicub_sg
      integer, intent(in) :: iZ, tbl
      real(sp), intent(in) :: logR, logT
      real(sp), intent(out) :: result(3)
      integer, intent(out) :: ierr
      real(sp) :: res(6)
      real(sp), pointer :: logR_ary(:), logT_ary(:)
      type(kapCN_set), pointer :: k

      k => kCN(iZ)
      logR_ary => kapCN_logR
      logT_ary => kapCN_logT

      call interp_evbicub_sg( logR, logT, &
                              logR_ary, num_logR, &
                              logT_ary, num_logT, &
                              ilinx, iliny, &
                              k% t(tbl)% kap, num_logR, &
                              ict, res, ierr)

      result = res(1:3)

      end subroutine kapCN_interp_RT


      subroutine kapCN_shutdown
      kapCN_is_initialized = .false.
      end subroutine kapCN_shutdown

      end module kapCN
