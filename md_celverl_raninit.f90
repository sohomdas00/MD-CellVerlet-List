!---------------------------3D LJ Mixture at RHO=.3---------------------------------
!-----------------------------implementing NVE MD-----------------------------------
!----------------using a combination of cell and verlet lists-----------------------
!-------------------------particles are randomly placed-----------------------------
    module parameters
      integer*4,parameter :: n=5000 !no. of particles
      integer*4 :: cntmax,idum
      integer*4 :: cnt,npart
      real*8 :: tmax,rho,en,rcut,en2,en1,en_var,en_sdev,e_drift
      real*8 :: l,sigma,ecut,fcorr,dt,t,etot
      real*8 :: fx(n),fy(n),fz(n),x(n),y(n),z(n),vx(n),vy(n),vz(n)
      integer*4 :: ncel
      integer*4, allocatable :: hoc(:,:,:),l1(:)
      integer*4,dimension(n) :: nlist    ! neighbours of each particle
      integer*4 :: list(n,5000)        ! assume max neigh per part = 500
      real*8 :: xv(n),yv(n),zv(n)
      real*8 :: rn,rc,rv
    end module parameters

    program molecular_dynamics
    use parameters
    implicit none
    integer*4 :: i,j1
    real*8 :: dx1,dy1,dz1,dist1,zl,yl,xl
    character(len=100) :: name1
    npart = n
    idum = 45678

    rho = 0.3
    l = REAL(int((dble(n)/rho)**(1d0/3d0))) !length of cubical box
    sigma = 1.d0   !diameter of particles
    rcut = 2.5d0*sigma
    t = 0.0d0
    ecut = 4.d0*((1.d0/(rcut**12d0))-(1.d0/(rcut**6d0)))
    fcorr = 48.d0*((0.5d0/(rcut**7d0))-(1.d0/(rcut**13d0)))    
    rc = 2.5*sigma
    rv = 2.7*sigma
    xl = l
    yl = l
    zl = l   
    
    ncel = int(xl/rc)
    rn = xl/ncel
    !print*, l,rv,rc,rn,ncel
    !rho = REAL(npart)/(xl**2)   ! ideal gas number densityssss
    allocate(l1(n),hoc(ncel,ncel,ncel))

    dt=0.005

    call init
    tmax = 50.0d0
    !cntmax = int(tmax/dt)  ! energy variance is calculated in all the runs for equal time
    !print *, cntmax
    cntmax = 3000
    en2 = 0.d0
    en1 = 0.d0
    en_var = 0.d0
    en_sdev = 0.d0


    call new_nvlist
    call force
    cnt = 0       !number of time steps progressed
    do while (cnt.le.cntmax)
     
       if (mod(cnt,500).eq.0) then
        name1='snap_'
        call addnumtostring(name1,cnt)
        open (51,file=name1,status='unknown')
        do i=1,n
          write(51,'(6(f13.10,2X))') x(i),y(i),z(i),vx(i),vy(i),vz(i)
        enddo
       endif
       call integrate
       call sample      
       print*,cnt
       cnt=cnt+1
       t=t+dt
    enddo
      !dt = dt + 1d-03
      !print *, dt
      !if (dt .le. 1.0001d-02) goto 15
    end program molecular_dynamics
!----------------------------------------------------------------------------------------------------- 
!*****************************************************************************************************
!--------------------------------------------------------------------
      SUBROUTINE new_nlist   ! creates a new cell list
      use parameters
      implicit none
      integer*4 :: i,icel,jcel,kcel

      do icel = 1,ncel
        do jcel = 1,ncel
          do kcel = 1,ncel
          hoc(icel,jcel,kcel) = 0    ! 0 implies empty cell
          enddo
        enddo
      enddo

      do i=1,n
        icel = int(x(i)/rn)+1
        if (icel.gt.ncel)icel=1
        jcel = int(y(i)/rn)+1
        if (jcel.gt.ncel)jcel=1
        kcel = int(z(i)/rn)+1
        if (kcel.gt.ncel)kcel=1
        l1(i) = hoc(icel,jcel,kcel)  ! link list the head of chain
        hoc(icel,jcel,kcel) = i      ! update the head of cell
        !write(48,*)l1(i)
      enddo
      return
      end
!-------------------------------------------------------------------------------------------------------
     SUBROUTINE new_nvlist     ! makes Verlet list using cell list
      use parameters  
      implicit none
      integer*4 :: i,j,k,ll,m,icel,jcel,kcel,ii,jj,kk
      real*8 :: xr,yr,zr,dr2,dr
      real*8 :: xl,yl,zl
      xl = l
      yl = l
      zl = l

      call new_nlist            
      do i=1,n
        nlist(i) = 0
        xv(i) = x(i)    ! store prev posn
        yv(i) = y(i)
        zv(i) = z(i)
      enddo
      do i=1,n
        icel = int(x(i)/rn)+1
        if (icel.gt.ncel)icel=1
        jcel = int(y(i)/rn)+1
        if (jcel.gt.ncel)jcel=1
        kcel = int(z(i)/rn)+1
        if (kcel.gt.ncel)kcel=1

        do k = -1,1       ! loop over neighbours of each cell
          do ll = -1,1
            do m = -1,1
             ii = icel+k
             jj = jcel+ll
             kk = kcel+m
             if (ii.gt.ncel)ii = 1
             if (jj.gt.ncel)jj = 1
             if (kk.gt.ncel)kk = 1
             if (ii.lt.1)ii = ncel
             if (jj.lt.1)jj = ncel
             if (kk.lt.1)kk = ncel
             j = hoc(ii,jj,kk)
             !print*, ii,jj,kk,j
             do while (j .ne. 0) 
              if (j.gt.i) then
                xr = x(i) - x(j)
                yr = y(i) - y(j)
                zr = z(i) - z(j)
                if (xr.gt.xl/2.0)xr = xr - xl
                if (yr.gt.yl/2.0)yr = yr - yl
                if (zr.gt.zl/2.0)zr = zr - zl
                if (xr.lt.-xl/2.0)xr = xr + xl
                if (yr.lt.-yl/2.0)yr = yr + yl
                if (zr.lt.-zl/2.0)zr = zr + zl
                dr2 = xr**2 + yr**2 + zr**2
                dr = sqrt(dr2)
                if (dr.lt.rv) then
                  nlist(i) = nlist(i) + 1
                  !nlist(j) = nlist(j) + 1
                  list(i,nlist(i)) = j
                  !list(j,nlist(j)) = i
                endif
              endif
              j = l1(j)   ! linked list element points to next particle
             enddo       ! looping until no other particles
          enddo         
        enddo           ! neighbor loop ends
      enddo             ! particle loop ends
      enddo
      return
      end
!------------------------------------------------------------------------------------------------------
      SUBROUTINE force
      use parameters
      implicit none
      integer*4 :: i,j,jj
      real*8 :: xr,yr,zr,r2,r1,rcutsq,r2i,r6i,ff
      real*8 :: xl,yl,zl
      rcutsq = rcut**2
      xl = l
      yl = l
      zl = l
      
        fx = 0.0
        fy = 0.0
        fz = 0.0
      en = 0.d0
      do i=1,n
        do jj=1,nlist(i)
          j = list(i,jj)
          xr = x(i)-x(j)
          yr = y(i)-y(j)
          zr = z(i)-z(j)
          if (xr.gt.xl/2.0)xr=xr-xl
          if (yr.gt.yl/2.0)yr=yr-yl
          if (zr.gt.zl/2.0)zr=zr-zl
          if (xr.lt.-xl/2.0)xr=xr+xl
          if (yr.lt.-yl/2.0)yr=yr+yl
          if (zr.lt.-zl/2.0)zr=zr+zl
          r2 = xr**2 + yr**2 + zr**2
          if (r2.lt.rcut*rcut) then
            r1 = sqrt(r2)
            r2i = 1.0/r2
            r6i = r2i**3
            ff = 48.0*r2i*r6i*(r6i-0.50)
            fx(i)=fx(i)+(ff*xr+fcorr*xr/r1)
	          !write(*,'(9(es10.3,X))')xr,yr,zr,ff,fcorr,xr,r2,r1,xr/r1
            fx(j)=fx(j)-(ff*xr+fcorr*xr/r1)
            fy(i)=fy(i)+(ff*yr+fcorr*yr/r1)
            fy(j)=fy(j)-(ff*yr+fcorr*yr/r1)
            fz(i)=fz(i)+(ff*zr+fcorr*zr/r1)
            fz(j)=fz(j)-(ff*zr+fcorr*zr/r1)
            en=en+(4.0*r6i*(r6i-1.0))-ecut-fcorr*(r1-rcut)
          endif
        enddo
      enddo
      return
      end

!-----------------------------------------------------------------------------------------------------
  subroutine integrate
  use parameters
  implicit none
  integer*4 :: i,j
  real*8 :: sumvx,sumvx2,sumvy,sumvy2,sumvz,sumvz2,vxn,vyn,vzn,temp1,ketot,entot,sumv
  real*8 :: fxi(n),fyi(n),fzi(n),dx,dy,dz,dist
  character(len=100) :: name1

  name1 = 'energies_time_5e-03_raninit'
  open(15,file=name1,status='unknown')
  open(24,file='vel-centreofmass',status='unknown')

  sumvx = 0.d0
  sumvy = 0.d0
  sumvz = 0.d0
  sumvx2 = 0.d0
  sumvy2 = 0.d0
  sumvz2 = 0.d0

  do i=1,n
    fxi(i)=fx(i)
    fyi(i)=fy(i)  ! stores forces at the ith time step
    fzi(i)=fz(i)
 enddo
  do i=1,n
    x(i) = x(i) + vx(i)*dt + 0.5d0*fx(i)*(dt**2d0)
    y(i) = y(i) + vy(i)*dt + 0.5d0*fy(i)*(dt**2d0)
    z(i) = z(i) + vz(i)*dt + 0.5d0*fz(i)*(dt**2d0)
  end do

  do i=1,n
!_______periodic boundary conditions___________________
    if (x(i) .gt. l) x(i) = x(i) - l
    if (x(i) .lt. 0.d0)x(i) = x(i) + l
    if (y(i) .gt. l) y(i) = y(i) - l
    if (y(i) .lt. 0.d0) y(i) = y(i) + l
    if (z(i).gt. l) z(i) = z(i) - l
    if (z(i).lt. 0.d0) z(i) = z(i) + l
!______________________________________________________
  enddo
  do j=1,n
    dx = x(j) - xv(j)
    dy = y(j) - yv(j)
    dz = z(j) - zv(j)
    dx = dx - l*nint(dx/l)
    dy = dy - l*nint(dy/l)
    dz = dz - l*nint(dz/l)
    dist = sqrt(dx**2 + dy**2 + dz**2)
    if (dist.gt.rv-rc) then
      print*, 'b'
      call new_nvlist
   endif
  enddo 
  call force    !stores forces at (i+1)th time step
  do i=1,n
    vxn = vx(i) + 0.5d0*(fx(i)+fxi(i))*dt
    vyn = vy(i) + 0.5d0*(fy(i)+fyi(i))*dt
    vzn = vz(i) + 0.5d0*(fz(i)+fzi(i))*dt
    sumvx = sumvx + vxn
    sumvy = sumvy + vyn
    sumvz = sumvz + vzn
    sumvx2 = sumvx2 + vxn**2d0
    sumvy2 = sumvy2 + vyn**2d0
    sumvz2 = sumvz2 + vzn**2d0
    vx(i) = vxn
    vy(i) = vyn
    vz(i) = vzn
  enddo
  temp1 = (sumvx2 + sumvy2 + sumvz2)/(3.d0*dble(n))   ! instantaneous temperature
  ketot = 0.5d0*(sumvx2+sumvy2+sumvz2)/dble(n)        ! total ke per particle	
  etot = (en + 0.5d0*(sumvx2+sumvy2+sumvz2))/dble(n)  ! total energy per particle
  entot = en/dble(n)                                  ! total pe per particle
  write(15,'(i5,X,5(f13.10,2X))') cnt,ketot,entot,etot,temp1
  write(24, '(3(f13.10,2X))') sumvx,sumvy,sumvz
  return
  end
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
  subroutine sample
  use parameters
  real*8 :: e0
  open(16,file='dt_energy-variance_raninit',status='unknown')
	
  if (cnt.eq.1)e0 = etot
  e_drift = e_drift + abs((e0-etot)/e0)
  en2 = en2 + etot**2d0
  en1 = en1 + etot
	!print *, en2,en1
  if (cnt.eq.cntmax) then
    en2 = en2/dble(cntmax)
    en1 = en1/dble(cntmax)
    e_drift = e_drift/dble(cntmax)
    print *, en2,en1
    en_var = en2 - en1*en1
    en_sdev = dsqrt(en_var)
    write(16,'(4(f13.10X))') dt,en_var,en_sdev,e_drift
  endif
  return
  end
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
      SUBROUTINE init
      use parameters
      implicit none
      real*8 :: dist,ran3
      real*8 :: fs,sumv2,sumvx,sumvy,sumvz,sumvx2,sumvy2,sumvz2,temp,fsx,fsy,fsz
      integer*4, dimension(:,:,:), allocatable :: cell
      integer*4 :: i, j, k, i1, i2, p, q, r, ll, p1, q1
      real*8 :: x0,y0,z0,x1,y1,z1,xr,yr,zr,dr
      !npart = dble(n)      ! number of particles
      !rho = 0.3d0          ! number density of particles   
      temp = 0.5          ! reduced temperature   
      !ll = int(l)
      !allocate (cell(ll,ll,ll))
      i1=1  !dummy index stores the number of particles in box
      dist = sigma**2d0 !square of the minimum permissible distance b/w two particles
      	  
      open(13,file='rho_0.3_raninit',status='unknown')    

!------ randomly generating 1st particle -----------------------------------------------------
      x0 = l*ran3(idum)
      y0 = l*ran3(idum)
      z0 = l*ran3(idum)
      x(i1) = x0
      y(i1) = y0
      z(i1) = z0
15   x1 = l*ran3(idum)
     y1 = l*ran3(idum)
     z1 = l*ran3(idum)
     do i=1,i1
       xr = x1 - x(i)
       yr = y1 - y(i)
       zr = z1 - z(i)
       xr = xr - l*nint(xr/l)
       yr = yr - l*nint(yr/l)
       zr = zr - l*nint(zr/l)
       dr = sqrt(xr**2 + yr**2 +zr**2)
       if (dr.lt.sigma) goto 15
      enddo
    i1 = i1 + 1
    x(i1) = x1
    y(i1) = y1
    z(i1) = z1
	!print *,p,q,r,cell(p,q,r)
     if (i1.lt.n) goto 15

!-------write particle co-ordinates in file ------------------------------------------------
     do i=1,n
      write(13,*) x(i),y(i),z(i)
    enddo
	
        sumvx=0.d0
        sumvy=0.d0
        sumvz=0.d0
        sumvx2=0.d0
        sumvy2=0.d0
        sumvz2=0.d0
        do i=1,n
         vx(i)=ran3(idum)-0.5d0
         vy(i)=ran3(idum)-0.5d0
         vz(i)=ran3(idum)-0.5d0
         sumvx=sumvx+vx(i)
         sumvy=sumvy+vy(i)
         sumvz=sumvz+vz(i)
         sumvx2=sumvx2+vx(i)**2d0
         sumvy2=sumvy2+vy(i)**2d0
         sumvz2=sumvz2+vz(i)**2d0
       enddo
       sumvx=sumvx/dble(n)
       sumvy=sumvy/dble(n)
       sumvz=sumvz/dble(n)
       sumvx2=sumvx2/dble(n)
       sumvy2=sumvy2/dble(n)
       sumvz2=sumvz2/dble(n)
       sumv2=sumvx2+sumvy2+sumvz2
       fsx=dsqrt(temp/sumvx2)
       fsy=dsqrt(temp/sumvy2)
       fsz=dsqrt(temp/sumvz2)
       do i=1,n
        vx(i)=(vx(i)-sumvx)*fsx
        vy(i)=(vy(i)-sumvy)*fsy
        vz(i)=(vz(i)-sumvz)*fsz
	!print*, vx(i),vy(i),vz(i)
       enddo
       return
       end
!*********************************************************************************************
!---------------------------------------------------------------------------------------------      	  
!---------------------------------------------------------------------------------------------
!-----Random Number Generator taken from Num. Rec. -------------------------------------------    
       FUNCTION ran3(idum)
        INTEGER*4 idum
        INTEGER*4 MBIG,MSEED,MZ
  !      REAL MBIG,MSEED,MZ
        REAL*8 ran3,FAC
        PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
  !     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
        INTEGER*4 i,iff,ii,inext,inextp,k
        INTEGER*4 mj,mk,ma(55)
  !     REAL mj,mk,ma(55)
        SAVE iff,inext,inextp,ma
        DATA iff /0/
        if(idum.lt.0.or.iff.eq.0)then
          iff=1
          mj=MSEED-iabs(idum)
          mj=mod(mj,MBIG)
          ma(55)=mj
          mk=1
          do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
  11      continue
          do 13 k=1,4
            do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
  12        continue
  13      continue
          inext=0
          inextp=31
          idum=1
        endif
        inext=inext+1
        if(inext.eq.56)inext=1
        inextp=inextp+1
        if(inextp.eq.56)inextp=1
        mj=ma(inext)-ma(inextp)
        if(mj.lt.MZ)mj=mj+MBIG
        ma(inext)=mj
        ran3=mj*FAC
        return
        END
!__________________________________________________________________________________________________
! Subr: Generating files
       SUBROUTINE addnumtostring(string,number)
       implicit none
       integer*4 i,strlen,number,nodig,num,snum
       character*(*) string
       snum=number
       do i=len(string),1,-1
       if (string(i:i) .ne. ' ') goto 10
       enddo
   10  strlen=i
       nodig=int(log10(1.0*snum+0.1))+1
       do i=nodig,1,-1
       num=snum/10**(i-1)
       string(strlen+1:strlen+1)=char(48+num)
       strlen=strlen+1
       snum=snum-num*10**(i-1)
       enddo
       return
       end
       integer*4 function strlen(string)
       character*(*) string
       integer*4 i
       do i=len(string),1,-1
       if (string(i:i).ne. ' ') goto 100
       enddo ! i
 100   strlen=i
       end
!------------------------------------------------------------------------------------------------------  
 
