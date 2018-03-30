#include "macro.h"
subroutine construct_matrix(iint)
  use openarray
  use variables
  use config
  use pvariable
  !use petscvec
  use petscsys
  implicit none
!#include "petsc.h"
!petsc variable(global value)---------------
!  PetscInt  :: p_size,p_rank,p_N,p_numperrow,p_i,p_j,p_nn
!  PetscScalar :: p_aprA,p_aprb
!  Mat       :: p_A
!  Vec       :: p_b,p_x
!  KSP       :: P_ksp
!  PC        :: pc
!  PetscInt,dimension(:),allocatable :: ix
!-------------------------------------------

  integer :: i,j,k,m,ierr
  integer :: iint
  type(array) :: aa1,aa2,aa3,aa4,bb1,bb2,bb3,bb4,&
    ga1,ga2,gb1,gb2,gc1,gc2,gen
  type(array) :: r,r1,r2
  type(array) :: aaf,bbf,ccf,fff
  type(array) :: ddx,ddy,dq
  type(array) :: fsm_tmp
  integer :: ind(3)
!  integer, allocatable :: ix(:)
  real(kind=8) :: tmp0,tmp1,tmp2,tmp3,tmp4
  real(kind=8),pointer :: paa1(:,:,:),paa2(:,:,:),paa3(:,:,:),&
                    paa4(:,:,:),pbb1(:,:,:),pbb2(:,:,:),&
                    pbb3(:,:,:),pbb4(:,:,:),pga1(:,:,:),&
                    pga2(:,:,:),pgb1(:,:,:),pgb2(:,:,:),&
                    pgc1(:,:,:),pgc2(:,:,:),pgen(:,:,:),&
                    paaf(:,:,:),pbbf(:,:,:),pdx(:,:,:), &
                    pddx(:,:,:),pddy(:,:,:),pdy(:,:,:), &
                    pdq(:,:,:) ,pfsm(:,:,:),pr(:,:,:),pq(:,:,:)
  integer,pointer:: pnumnelt(:,:,:)
  aa1=mat_zeros;gen=mat_zeros;aa2=mat_zeros
  aa3=mat_zeros;aa4=mat_zeros;bb1=mat_zeros
  bb2=mat_zeros;bb3=mat_zeros;bb4=mat_zeros
  ga1=mat_zeros;ga2=mat_zeros;gb1=mat_zeros
  gc1=mat_zeros;gc2=mat_zeros;aaf=mat_zeros
  gb2=mat_zeros;gc1=mat_zeros;gen=mat_zeros
  bbf=mat_zeros;ccf=mat_zeros;fff=mat_zeros
  apr=0; ia=0; ja=0 
  
  ddx=AXF(dx); ddy=AXF(dy)
  fsm_tmp=fsm; dq=dtf
  call get_local_buffer(paa1,aa1);call get_local_buffer(paa2,aa2)
  call get_local_buffer(paa3,aa3);call get_local_buffer(paa4,aa4)
  call get_local_buffer(pbb1,bb1);call get_local_buffer(pbb2,bb2)
  call get_local_buffer(pbb3,bb3);call get_local_buffer(pbb4,bb4)
  call get_local_buffer(pga1,ga1);call get_local_buffer(pga2,ga2)
  call get_local_buffer(pgb1,gb1);call get_local_buffer(pgb2,gb2)
  call get_local_buffer(pgc1,gc1);call get_local_buffer(pgc2,gc2)
  call get_local_buffer(pgen,gen);call get_local_buffer(paaf,aaf)
  call get_local_buffer(pbbf,bbf);call get_local_buffer(pddx,ddx)
  call get_local_buffer(pddy,ddy);call get_local_buffer(pdq,dq)
  call get_local_buffer(pfsm,fsm);call get_local_buffer(pdx,dx)
  call get_local_buffer(pdy,dy)  ;call get_local_buffer(pnumnelt,numnelt)
  
  call update_ghost(ddx)
  call update_ghost(ddy)
  call set(sub(fsm,1,':'),0.d0)
  call set(sub(fsm,im,':'),0.d0)
  call set(sub(fsm,':',1),0.d0)
  call set(sub(fsm,':',jm),0.d0)
  call set(sub(dq,1,':'),sub(dtf,2,':'))
  call set(sub(dq,im,':'),sub(dtf,imm1,':'))
  call set(sub(dq,':',1),sub(dtf,':',2))
  call set(sub(dq,':',jm),sub(dtf,':',jmm1))
  
  call update_ghost(dq)
  call update_ghost(fsm)
  call coef2(aaf,bbf,ccf,fff)
  call disp(aaf,'aaf=')
  call disp(bbf,'bbf=')
  call update_ghost(aaf)
  call update_ghost(bbf)
  aaf=bbf
  call disp(aaf,'aaf=')
  call update_ghost(numnelt)

!petsc init--------------------------------------------------------
  PETSC_COMM_WORLD=MPI_COMM_WORLD
  p_N=total_nums
  p_nn=cnn
  p_numperrow=15
  if(iint==2)then
  print*,"allocate petsc memory"
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,p_rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,p_size,ierr)
  call MatCreate(PETSC_COMM_WORLD,p_A,ierr)
  call MatSetSizes(p_A,PETSC_DECIDE,PETSC_DECIDE,p_N,p_N,ierr)  
  !call MatSetSizes(p_A,p_nn,p_nn,p_N,p_N,ierr)  
  call MatSetFromOptions(p_A,ierr)
  if (p_size>1) then
  call MatMPIAIJSetPreallocation(p_A,p_numperrow,&
      PETSC_NULL_INTEGER,p_numperrow,PETSC_NULL_INTEGER,ierr)
  else
  call MatSeqAIJSetPreallocation(p_A,p_numperrow,&
      PETSC_NULL_INTEGER,ierr)
  endif
!  endif
  call VecCreate(PETSC_COMM_WORLD,p_b,ierr)
  !call VecSetType(p_b,VECMPI,ierr)
  !call VecSetSizes(p_b,p_nn,p_N,ierr)
  call VecSetSizes(p_b,PETSC_DECIDE,p_N,ierr)
  call VecSetFromOptions(p_b,ierr)
  call VecDuplicate(p_b,p_x,ierr)
  call KSPCreate(PETSC_COMM_WORLD,p_ksp,ierr)
  call KSPSetOperators(p_ksp,p_A,p_A,ierr)
  call KSPSetTolerances(p_ksp,1.e-20,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,5000,ierr)
!  endif
!  print*,'cnn=',cnn,'nums_s=',nums_s,'total_nums=',total_nums,'rank=',get_rank()
  allocate(ix(0:total_nums-1))
  allocate(qq(1:total_nums))
  do p_i=0,total_nums-1
    ix(p_i)=p_i
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,p_N,p_x_all,ierr)
  call ISCreateGeneral(PETSC_COMM_WORLD,p_N,ix,PETSC_COPY_VALUES,p_from,ierr)
  call ISCreateGeneral(PETSC_COMM_WORLD,p_N,ix,PETSC_COPY_VALUES,p_to,ierr)
  endif
!  if(get_rank()==0)then
!    print*,ix
!  endif
!------------------------------------------------------------------
  if(p_rank==0) then
    print*,iint
  endif

  ind = shape(pfsm)
  m=0
  do k=2,kb
    do j=2,ind(2)-1
      do i=2,ind(1)-1
        if(pfsm(i,j,2)>0.9d0)then
          paa1(i+1,j,k+1)=0.25d0*paaf(i+1,j,k+1)/pdx(i,j,2)
          paa2(i-1,j,k+1)=-0.25d0*paaf(i-1,j,k+1)/pdx(i,j,2)
          pbb1(i,j+1,k+1)=0.25d0*pbbf(i,j+1,k+1)/pdy(i,j,2)
          pbb2(i,j-1,k+1)=-0.25d0*pbbf(i,j-1,k+1)/pdy(i,j,2)
          pga1(i+1,j,k)=dz1(k-1)/pdx(i,j,2)/pddx(i,j,2)*pdq(i+1,j,2)
          pga2(i-1,j,k)=dz1(k-1)/pdx(i,j,2)/pddx(i-1,j,2)*pdq(i-1,j,2)
          pgb1(i,j+1,k)=dz1(k-1)/pdy(i,j,2)/pddy(i,j,2)*pdq(i,j+1,2)
          pgb2(i,j-1,k)=dz1(k-1)/pdy(i,j,2)/pddy(i,j-1,2)*pdq(i,j-1,2)
          pgc1(i,j,k+1)=1.e0/(dz1(k-1)*pdq(i,j,2))*ramp
          if(k>2)then
            paa3(i+1,j,k-1)=-0.25d0*paaf(i+1,j,k-1)/pdx(i,j,2)
            paa4(i-1,j,k-1)=0.25d0*paaf(i-1,j,k-1)/pdx(i,j,2)
            pbb3(i,j+1,k-1)=-0.25d0*pbbf(i,j+1,k-1)/pdy(i,j,2)
            pbb4(i,j-1,k-1)=0.25d0*pbbf(i,j-1,k-1)/pdy(i,j,2)
            pgc2(i,j,k-1)=1.d0/(dzz1(k-2)*pdq(i,j,2))*ramp
            pgen(i,j,k) = - pdq(i,j,2)*dz1(k-1)*((1.d0/pddx(i,j,2)+1.d0/pddx(i-1,j,2))/ &
                  pdx(i,j,2)+(1.d0/pddy(i,j,2)+1.d0/pddy(i,j-1,2))/pdy(i,j,2))- &
                  (1.d0/dzz1(k-2)+1.d0/dzz1(k-1))/(pdq(i,j,2))*ramp
          else
            pgen(i,j,k) = - pdq(i,j,2)*dz1(k-1)*((1.d0/pddx(i,j,2)+1.d0/pddx(i-1,j,2))/ &
                  pdx(i,j,2)+(1.d0/pddy(i,j,2)+1.d0/pddy(i,j-1,2))/pdy(i,j,2))-&
                  (1.d0/dzz1(k-1)+1.d0/dzz1(k-1))/(pdq(i,j,2))*ramp
          endif
          if(k==kb)then
            pga1(i+1,j,k)=pga1(i+1,j,k)+paa1(i+1,j,k+1)
            pga2(i-1,j,k)=pga2(i-1,j,k)+paa2(i-1,j,k+1)
            pgb1(i,j+1,k)=pgb1(i,j+1,k)+pbb1(i,j+1,k+1)
            pgb2(i,j-1,k)=pgb2(i,j-1,k)+pbb2(i,j-1,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pgc1(i,j,k+1)
          endif
          if(k==3)then
            paa2(i+1,j,k-1)=0.d0
            pbb2(i,j+1,k-1)=0.d0 !need to debug
            pbb4(i,j-1,k-1)=0.d0
            pgc2(i,j,k-1)=0.d0
          endif

          if(i==3)then
            if(k>2)then
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+paa4(i-1,j,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+paa2(i-1,j,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pga2(i-1,j,k)
          else
            if(pfsm(i-1,j,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pga2(i-1,j,k)
            endif
          endif
     
          if(i==im)then
            if(k>2)then 
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+paa3(i+1,j,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+paa1(i+1,j,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pga1(i+1,j,k)
          else
            if(pfsm(i+1,j,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pga1(i+1,j,k)
            endif
          endif

          if(j==3)then
            if(k>2)then
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+pbb4(i,j-1,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+pbb2(i,j-1,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pgb2(i,j-1,k)
          else
            if(pfsm(i,j-1,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pgb2(i,j-1,k)
            endif
          endif

          if(j==jm)then
            if(k>2)then
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+pbb3(i,j+1,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+pbb1(i,j+1,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pgb1(i,j+1,k)
          else
            if(pfsm(i,j+1,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pgb1(i,j+1,k)
            endif
          endif
        if(k<kb)then 
          if(pfsm(i+1,j,2)>0.9d0)then
            p_aprA=paa1(i+1,j,k+1)
            p_j=pnumnelt(i+1,j,k+1)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(pfsm(i,j+1,2)>0.9d0)then
            p_aprA=pbb1(i,j+1,k+1)
            p_j=pnumnelt(i,j+1,k+1)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(pfsm(i,j,2)>0.9d0)then
            p_aprA=pgc1(i,j,k+1)
            p_j=pnumnelt(i,j,k+1)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(pfsm(i,j-1,2)>0.9d0)then
            p_aprA=pbb2(i,j-1,k+1)
            p_j=pnumnelt(i,j-1,k+1)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(pfsm(i-1,j,2)>0.9d0)then
            p_aprA=paa2(i-1,j,k+1)
            p_j=pnumnelt(i-1,j,k+1)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif
        endif

          if(pfsm(i+1,j,2)>0.9d0)then
            p_aprA=pga1(i+1,j,k)
            p_j=pnumnelt(i+1,j,k)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(pfsm(i,j+1,2)>0.9d0)then
            p_aprA=pgb1(i,j+1,k)
            p_j=pnumnelt(i,j+1,k)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif
 
          p_aprA=pgen(i,j,k)
          p_j=pnumnelt(i,j,k)-1
          p_i=pnumnelt(i,j,k)-1
          call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)

          if(pfsm(i,j-1,2)>0.9d0)then
            p_aprA=pgb2(i,j-1,k)
            p_j=pnumnelt(i,j-1,k)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(pfsm(i-1,j,2)>0.9d0)then
            p_aprA=pga2(i-1,j,k)
            p_j=pnumnelt(i-1,j,k)-1
            p_i=pnumnelt(i,j,k)-1
            call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
          endif

          if(k>=3)then
            if(pfsm(i+1,j,2)>0.9d0)then
              p_aprA=paa3(i+1,j,k-1)
              p_j=pnumnelt(i+1,j,k-1)-1
              p_i=pnumnelt(i,j,k)-1
              call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
            endif

            if(pfsm(i,j+1,2)>0.9d0)then
              p_aprA=pbb3(i,j+1,k-1)
              p_j=pnumnelt(i,j+1,k-1)-1
              p_i=pnumnelt(i,j,k)-1
              call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
            endif

            if(pfsm(i,j,2)>0.9d0)then
              p_aprA=pgc2(i,j,k-1)
              p_j=pnumnelt(i,j,k-1)-1
              p_i=pnumnelt(i,j,k)-1
              call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
            endif

            if(pfsm(i,j-1,2)>0.9d0)then
              p_aprA=pbb4(i,j-1,k-1)
              p_j=pnumnelt(i,j-1,k-1)-1
              p_i=pnumnelt(i,j,k)-1
              call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
           endif

           if(pfsm(i-1,j,2)>0.9d0)then
             p_aprA=paa4(i-1,j,k-1)
             p_j=pnumnelt(i-1,j,k-1)-1
             p_i=pnumnelt(i,j,k)-1
             call MatSetValues(p_A,1,p_i,1,p_j,p_aprA,INSERT_VALUES,ierr)
           endif
         endif
       endif
     enddo
   enddo
 enddo
 call MatAssemblyBegin(p_A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(p_A,MAT_FINAL_ASSEMBLY,ierr)
 !call MatView(p_A,PETSC_VIEWER_STDOUT_WORLD,ierr)
 call update_ghost(aa1)
 call update_ghost(aa2)
 call update_ghost(aa3)
 call update_ghost(aa4)
 call update_ghost(bb1)
 call update_ghost(bb2)
 call update_ghost(bb3)
 call update_ghost(bb4)
 call update_ghost(ga1)
 call update_ghost(ga1)
 call update_ghost(gb1)
 call update_ghost(gb2)
 call update_ghost(gc1)
 call update_ghost(gc2)
 call update_ghost(gen)
! call disp(aa1,'aa1=')
! call disp(aa2,'aa2=')
! call disp(aa3,'aa3=')
! call disp(aa4,'aa4=')
! call disp(bb1,'bb1=')
! call disp(bb2,'bb2=')
!! call disp(bb3,'bb3=')
!! call disp(bb4,'bb4=')
! call disp(ga1,'ga1=')
! call disp(ga2,'ga2=')
! call disp(gb1,'gb1=')
! call disp(gb2,'gb2=')
!! call disp(gc1,'gc1=')
! call disp(gc2,'gc2=')
! call disp(gen,'gen=')


  r1=uf*AXB(dq)*dz
  r2=vf*AYB(dq)*dz
  call coef2(aaf,bbf,ccf,fff)
  wq=-fff+wwf+ccf
  r=1.d0/dti2*(DXF(r1)+DYF(r2)-DZF(wq)*dz)
  !call disp(r,'r=')
  !stop
  call get_local_buffer(pr,r)    
  call get_local_buffer(pq,q) 
  do k=2,kb
    do j=2,ind(2)-1
      do i=2,ind(1)-1
        if(pfsm(i,j,2)>0.9d0)then
          m=pnumnelt(i,j,k)-1
          p_aprb=pr(i,j,k)
          call VecSetValues(p_b,1,m,p_aprb,INSERT_VALUES,ierr)
        endif
      enddo
    enddo
  enddo  
  call VecAssemblyBegin(p_b,ierr)
  call VecAssemblyEnd(p_b,ierr)
!  call VecView(p_b,PETSC_VIEWER_STDOUT_WORLD,ierr)

  !ksp solve-------------------:-----------------------------------
  call KSPSetOperators(p_ksp,p_A,p_A,ierr)
  call KSPGetPC(p_ksp,pc,ierr)
  call KSPSetFromOptions(p_ksp,ierr)
  call KSPSolve(p_ksp,p_b,p_x,ierr)
  !call KSPView(p_ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call VecGetValues(p_x,cnn,ix,qq,ierr)
  !call VecView(p_x,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
  call VecScatterCreate(p_x,p_from,p_x_all,p_to,p_ctx,ierr)
  call VecScatterBegin(p_ctx,p_x,p_x_all,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(p_ctx,p_x,p_x_all,INSERT_VALUES,SCATTER_FORWARD,ierr)
 ! call VecView(p_x_all,PETSC_VIEWER_STDOUT_SELF,ierr)
  call VecGetValues(p_x_all,p_N,ix,qq,ierr)
!  if(p_rank==0) then
!  do i=1,378
!  print*,'qq=',qq(i)
!  enddo
!  endif
!  stop
  !---------------------------------------------------------------

  !ksp destroy----------------------------------------------------
  if(iint==100)then
  print*,"deallocate petsc memory"
  deallocate(ix)
  deallocate(qq)
  call VecScatterDestroy(p_ctx,ierr)
  call VecDestroy(p_b,ierr)
  call VecDestroy(p_x,ierr)
  call MatDestroy(p_A,ierr)
  call KSPDestroy(p_ksp,ierr)
  endif
  !---------------------------------------------------------------
  m=0
  do k=2,kb
    do j=2,ind(2)-1
      do i=2,ind(1)-1
        if(pfsm(i,j,2)>0.9d0)then
          m=pnumnelt(i,j,k)
          pq(i,j,k)=qq(m)
        endif
      enddo
    enddo
  enddo
!  call disp(q,'q=')
  stop
  call set(sub(q,1,':',':'),sub(q,2,':',':'))
  call set(sub(q,im,':',':'),sub(q,imm1,':',':'))
  call set(sub(q,':',1,':'),sub(q,':',2,':'))
  call set(sub(q,':',jm,':'),sub(q,':',jmm1,':'))

  call set(sub(q,':',':',kb),sub(q,':',':',kbm1))
  call set(sub(q,':',':',1),0.d0)
!  call disp(q,'q=')

  fsm=fsm_tmp
end subroutine
