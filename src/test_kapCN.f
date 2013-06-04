      program test_kapCN

      use const_def, only: sp
      use const_lib, only: const_init
      use kapCN

      implicit none

      integer :: ierr
      real(sp) :: logR, logT, result(3), X, fC, fN, Z
      character(len=256) :: data_dir

      data_dir = '' !means use mesa_data_dir
      call const_init(data_dir,ierr)
      call kapCN_init(ierr)
      if(ierr/=0) write(*,*) ' ierr = ', ierr

      !this one requires full Z interpolation
      X=0.6
      Z=0.009999
      fC=1.11
      fN=1.22
      logR=-3.44343
      logT=3.377
      call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      call write_results

      !this one should only need one fixed Z call
      Z=0.01
      call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      call write_results


      X=0.75
      Z=0.01
      fC=1.9
      fN=1.4
      call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      call write_results

      X=0.75
      Z=0.031
      fC=1.9
      fN=1.4
      call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      call write_results

      logR=1.05
      call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
      call write_results

      call kapCN_shutdown

      contains 

      subroutine write_results
 1    format(a25,f11.6)
      write(*,*)
      write(*,1) ' X = ', X
      write(*,1) ' Z = ', Z
      write(*,1) 'fC = ', fC
      write(*,1) 'fN = ', fN
      write(*,1) 'logR=', logR
      write(*,1) 'logT=', logT
      write(*,1) 'log(kap)=', result(1)
      write(*,1) 'dlog(kap)/dlogT=', result(2)
      write(*,1) 'dlog(kap)/dlogR=', result(3)
      write(*,*) 'ierr = ', ierr
      write(*,*)
      end subroutine write_results

      end program test_kapCN
