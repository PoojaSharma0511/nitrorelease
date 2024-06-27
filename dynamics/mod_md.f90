Module mod_md
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! PES
integer,parameter::ntrain=196,nclass=2
real*8 alpha(ntrain),x_train(ntrain,nclass),xmin(nclass)
real*8 mass(nclass),omg(nclass)
real*8 gamma_B

!! dynamic variables
real*8 x(nclass),v(nclass),acc(nclass)
integer nsteps,ntraj
real*8,allocatable :: pop(:)
real*8 dt,tot_time,curr_time
real*8 path1,path2,path3
real*8 temperature
integer ifriction,iwrite,if_react,if_ts

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer,allocatable :: seed(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar

  open(10,file="md.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) dt
  read(10,*) tot_time
  read(10,*) ntraj
  read(10,*) temperature
  read(10,*) ifriction
  read(10,*) iwrite
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 
  nsteps=int(tot_time/dt)
  allocate(pop(nsteps))

  !-----------------------------------------------------------------  
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    Ntraj=Ntraj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*Ntraj
        seed=seed+1
!        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i
  real*8 t1,t2

  call cpu_time(t1)

  call setup_parameters
  call init_avg
  !call draw_pes

  do i=1,ntraj
    call init_cond
    call evolve
  enddo
  pop=pop/real(ntraj)
  open(100,file="pop.out")
  do i=1,nsteps
    write(100,*) i*dt*1.d15,pop(i)
  enddo
  close(100)
  open(101,file="cnts.out")
  write(101,*)path1,path2,path3
  close(101)
  call cpu_time(t2)
  open(10,file="output")
  write(10,*) "Total time  ",t2-t1
  close(10)

end subroutine main
!---------------------------------------------------------- 

subroutine setup_parameters
  implicit none
  integer i

mass(1)=6.85d0*amu2kg
mass(2)=6.46d0*amu2kg
!print*, temperature, beta
!stop


  open(30,file="alpha_NA_196_24jan.inp")
  open(31,file="kernel_NA_196_24jan.inp")
  do i=1,ntrain
    read(30,*)alpha(i)
    read(31,*)x_train(i,:)
  enddo
  close(30)
  close(31)

  xmin(1)=2.3d-10
  xmin(2)=1.3d-10

  omg(1)=1973.d0*2*pi*clight
  omg(2)=1034.d0*2*pi*clight

!  gamma_B=1000.d0*2*pi*clight
   gamma_B=omg(1)

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine init_avg
  implicit none

  pop=0.d0
  path1=0.d0
  path2=0.d0
  path3=0.d0

end subroutine init_avg
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i
  real*8 ak,sig_x,sig_p,rnd

  do i=1,nclass
    !ak=2/(hbar*omg(i))*dtanh(beta*hbar*omg(i)/2.d0) !! Wigner
    ak=1.d0/(kb*temperature)    !! Classical
    sig_x=1.d0/dsqrt(ak*mass(i)*omg(i)**2)
    sig_p=dsqrt(mass(i)/ak)
    call gaussian_random_number(rnd)
    x(i)=rnd*sig_x + xmin(i)
    call gaussian_random_number(rnd)
    v(i)=(1.d0/mass(i)*(rnd*sig_p))
  enddo
  
    call check_init_cond


end subroutine init_cond
!-----------------------------------------------------------------

subroutine check_init_cond
implicit none

  if(x(1)>2.4d-10 .or. x(2)<1.4d-10) then
    call init_cond
  endif
  
end subroutine check_init_cond
!_________________________________________________________________
subroutine evolve
  implicit none
  integer i
  real*8 pot,acc(nclass),energy,acc_old(nclass)
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass),vnew(nclass),xnew(nclass)

  call compute_acc(x,pot,acc)
  energy = pot+0.5*sum(mass*v*v)
  curr_time=0.d0

  if(ifriction==0) then
    do i=1,nsteps
      write(10,*) curr_time*1.d15,energy/wave_to_J
      write(11,*) curr_time*1.d15,x*1.d10,v
      x=x+v*dt+0.5*acc*dt**2
      v=v+0.5*acc*dt
      call compute_acc(x,pot,acc)
      v=v+0.5*acc*dt
      energy = pot+0.5*sum(mass*v*v)
      curr_time=i*dt
    enddo
  endif


  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
    do i=1,nsteps
      if(iwrite==1) then
        write(10,*) curr_time*1.d15,energy/wave_to_J
        write(11,*) curr_time*1.d15,x*1.d10,v
      endif
      call compute_pops(i)
      call compute_paths
      call stochastic_force(delta_r,delta_v,dt)
      xnew=x+c1*dt*v+c2*dt*dt*acc+delta_r
      acc_old=acc
      call compute_acc(xnew,pot,acc)
      vnew=c0*v+(c1-c2)*dt*acc_old+c2*dt*acc+delta_v
      call reflect(vnew,v,xnew,x)
      energy = pot+0.5*sum(mass*v*v)
      curr_time=i*dt
    enddo
  endif

end subroutine evolve
!-----------------------------------------------------------------  

subroutine compute_paths
  implicit none

  if(x(1)>2.1d-10 .and. x(2)<1.8d-10) then
    if_react=1
    if_ts=0
  endif

!  if(x(1)<1.8d-10 .and. x(2)>2.5d-10) then
!   if(if_ts==1) path1=path1+1
!    if(if_ts==2) path2=path2+1
 !   if(if_ts==0.and.if_react==1) path3=path3+1
 !   if_react=0
 !   if_ts=0
 ! endif

  if(if_ts==0) then

   if(x(1)<2.d-10 .and. x(2)<2.d-10 .and.  if_react==1) then
!      path1=path1+1
      if_ts=1
!write(6,*) "path1",curr_time*1.d15,x*1.d10
    endif

    if(x(1)>2.2d-10 .and. x(2)>2.d-10  .and. if_react==1) then
!      path2=path2+1
      if_ts=2
!write(6,*) "path2",curr_time*1.d15,x*1.d10
    endif

  endif

  if(x(1)<1.8d-10 .and. x(2)>2.5d-10) then
   if(if_ts==1) path1=path1+1
    if(if_ts==2) path2=path2+1
    if(if_ts==0.and.if_react==1) path3=path3+1
    if_react=0
    if_ts=0
  endif



end subroutine compute_paths
!-----------------------------------------------------------------  

subroutine reflect(vnew,v,xnew,x)
  implicit none
  real*8,intent(in)::vnew(nclass),xnew(nclass)
  real*8,intent(inout)::v(nclass),x(nclass)
  integer ireflect,doub_ref

  ireflect =0
 doub_ref=0
if(xnew(1)<1.25d-10 .or. xnew(1)>2.4d-10) then
if(xnew(2)<1.4d-10 .or. xnew(2)>3.0d-10) then
v=-v
ireflect=1
doub_ref=1
endif
endif

if(doub_ref==0) then
  if(xnew(1)<1.25d-10 .or. xnew(1)>2.4d-10) then
!write(6,*) "Before 1",curr_time*1.d15,x*1.d10,v
    v(1)=-v(1);v(2)=vnew(2)
    x(2)=xnew(2)
    ireflect=1
!write(6,*) "After 1",curr_time*1.d15,x*1.d10,v
  endif
  if(xnew(2)<1.4d-10 .or. xnew(2)>3.0d-10) then
!write(6,*) "Before 2",curr_time*1.d15,xnew*1.d10,vnew
    v(2)=-v(2);v(1)=vnew(1)
    x(1)=xnew(1)
    ireflect=1
!write(6,*) "After 2",curr_time*1.d15,x*1.d10,v
  endif
endif

  if(ireflect==0) then
    v=vnew
    x=xnew
  endif

end subroutine reflect
!-----------------------------------------------------------------  

subroutine compute_pops(i)
  implicit none
  integer,intent(in):: i
  real*8 m,c,y

  m=1.d0/1.2d0
  c=2.d-10
  y=m*(x(1)-1.2d-10)+c
  if(x(2)<y) pop(i)=pop(i)+1.d0


end subroutine compute_pops
!-----------------------------------------------------------------  

subroutine compute_acc(x,pot,acc)
  implicit none
  integer i,j
  real*8,intent(in) :: x(nclass)
  real*8,intent(out):: pot,acc(nclass)
  real*8 xx(2),kk(ntrain),kk_deriv(ntrain,nclass)!,kk_dd(ntrain,nclass)
  real*8 p2,cons,y_avg
  real*8 distsq,VV,delv_dels(nclass)

  p2=0.114     !196
  cons=29.935135941210703
  y_avg=0.d0!-551.3186522959184

  xx=x*1.d10
  do j=1,ntrain
     distsq=sum((xx-x_train(j,:))**2)/p2**2
     kk(j)=exp(-0.5*distsq)
     kk_deriv(j,1)=-kk(j)*(xx(1)-x_train(j,1))/p2**2*1.d10
     kk_deriv(j,2)=-kk(j)*(xx(2)-x_train(j,2))/p2**2*1.d10

     !kk_dd(j,1)=-kk(j)*(1-((xx(1)-x_train(j,1))/p2)**2)/p2**2*1.d20
     !kk_dd(j,2)=-kk(j)*(1-((xx(2)-x_train(j,2))/p2)**2)/p2**2*1.d20
   enddo
   VV=sum(alpha*kk)/cons+y_avg
   pot=VV* au2J

   delV_dels(1)=sum(alpha*kk_deriv(:,1))/cons*au2J
   delV_dels(2)=sum(alpha*kk_deriv(:,2))/cons*au2J

  acc=-delV_dels/mass

end subroutine compute_acc
!-----------------------------------------------------------------  

subroutine draw_pes
  implicit none
  integer i,j,n
  real*8 pot,acc(nclass)

  open(511, file='pes.out')
  do i = 1,100
     x(1)=1.d-10+2.d-10*i/100.d0
     do j = 1,100
         x(2)=1.d-10+2.5d-10*j/100.d0
         call compute_acc(x,pot,acc)
           write(511,'(4f19.7)') x*1.d10,pot/wave_to_J
     enddo
  write(511,'(4f19.7)')
  enddo
  close(511)
  !write(511,'(4f19.7)')

end subroutine draw_pes
!-------------------------------------------------------------------------

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!----------------------------------------------------------

subroutine stochastic_force(delr,delv,dt)
  !! stoachastic forces for langevin equation
  !! Not used for the Holstein model results 
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass

    sig_r=dt*dsqrt(kb*temperature/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

End Module mod_md
