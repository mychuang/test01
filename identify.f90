      program identify
      use netcdf
      implicit none
      integer :: i,j,k,l,m
      !
      integer,parameter :: NX=300,NY=300,NZ=60
      integer, parameter :: dx = 1000,dy = 1000, dz = 250
      real, parameter :: rlatd = 22.25, rlond = 119.75
      real, parameter :: lon_l = 118.29151, lat_l = 20.900108
      real, parameter :: lon_g = 9.72325634E-03, lat_g = 8.99927784E-03
!~background
      real :: u0(NX,NY,NZ),v0(NX,NY,NZ),qv0(NX,NY,NZ),t0(NX,NY,NZ)
!~ascat
      integer,parameter :: sea_info=358
      real :: sea_site(3,sea_info)
      integer :: sea_site_id(sea_info)
      character*13 :: infile_sea(sea_info)
      !
      INTEGER :: ID,YY,MM,DD,HH,NN,itime
      REAL :: pre,temp,td,RH,WIND,WINDIR
      real :: alt(sea_info),lon(sea_info),lat(sea_info)
      real :: spd(sea_info),dir(sea_info)
      !
      real,parameter :: badpt = -999.98999
      real :: ua(sea_info),va(sea_info)
      real,parameter :: u_hfactor=-7.7368259E-02
      real,parameter :: v_hfactor=0.612216
      real,parameter :: u_mfactor=-0.6058777
      real,parameter :: v_mfactor=-0.1411584
      real,parameter :: u_lfactor=6.8603404E-02
      real,parameter :: v_lfactor=0.1272687
      !
      real :: power_law
!~sitecheck
      real :: xxsta(sea_info),yysta(sea_info)
      integer :: ii,jj
      real :: uu(NX,NY),vv(NX,NY)
!~verify_bg
      real :: rmse_ub,rmse_vb
!~This is the name of the data file
      character (len = *), &
      parameter :: FILE_NAME = "./ncfdataf_00000000.nc"
      integer :: ncid
      integer :: lenth, dimid
      !
      character (len = *), parameter :: time = "time"
      character (len = *), parameter :: x0 = "x0"
      character (len = *), parameter :: y0 = "y0"
      character (len = *), parameter :: z0 = "z0"
      character (len = *), parameter :: lat0 = "lat0"
      character (len = *), parameter :: lon0 = "lon0"
      integer :: x0_varid, y0_varid, z0_varid
      integer :: lat0_varid, lon0_varid, time_varid
      integer :: x(NX),y(NY),z(NZ)
      integer :: time0(1)
      real :: lat1(NX,NY),lon1(NX,NY)
      !
      character (len = *), parameter :: ufct = "ufct"
      character (len = *), parameter :: vfct = "vfct"
      integer :: ufct1_varid, vfct1_varid
      real :: u1(NX,NY,NZ,1),v1(NX,NY,NZ,1)
!~verify analysis
      real :: rmse_ua,rmse_va,mvd
      real :: fx,f1,f2,sum_dir
      real :: check_dir(102)
!--------------------------------------------------------------------------
!~read background
      open(10,file='gradsout_calc',status='OLD',&
              form='unformatted',access='direct',recl=4*NX*NY)
        do m = 1,NZ
          read(10,rec=m) ((u0(l,k,m),k=1,NX),l=1,NY)
          read(10,rec=m+NZ) ((v0(l,k,m),k=1,NX),l=1,NY)
          read(10,rec=m+NZ*2) ((qv0(l,k,m),k=1,NX),l=1,NY)
          read(10,rec=m+NZ*3) ((t0(l,k,m),k=1,NX),l=1,NY)
        enddo
      close(10)
!~read ascat
      open(17,file='seasite_pro',form='formatted')
        do i=1,sea_info
          read(17,*) sea_site_id(i),sea_site(1,i),sea_site(2,i),&
                     sea_site(3,i)
        enddo
      close(17)

      open(12,file='./ascat/list',status='old',form='formatted')
        do i=1,sea_info
          read(12,"(a13)") infile_sea(i)
        enddo
      close(12)
!~process ascat
     
      do i=1,sea_info
        open(23,file='./ascat/'//infile_sea(i),&
                status='old',form='formatted')
          do while(.true.)
            read(23,110,end=3100)ID,YY,MM,DD,HH,NN,PRE,TEMP,TD,RH,&
                                 WIND,WINDIR
110         format(I6.6,I4.4,I2.2,I2.2,I2.2,I2.2,&
                   F6.1,F5.1,F5.1,F4.0,F4.1,F5.1)
            do k = 1,sea_info
              if(id.eq.sea_site_id(k))then
                alt(i) = sea_site(1,k)
                lon(i) = sea_site(2,k)
                lat(i) = sea_site(3,k)
                spd(i) = WIND
                dir(i) = WINDIR
              endif
            enddo
          enddo
3100      continue
        close(23)
      enddo

!~convert u;v & power law
      ua = 0.0 ; va = 0.0

      do i=1,sea_info
        if (dir(i).ne.999.9.and.dir(i).lt.361.) then
          ua(i) = spd(i)*SIN(dir(i)*0.01745329)
          va(i) = spd(i)*COS(dir(i)*0.01745329)
          if(spd(i).ge.10.)then
            ua(i) = ua(i) + u_hfactor
            va(i) = va(i) + v_hfactor
          elseif(spd(i).lt.10..and.spd(i).ge.5.)then
            ua(i) = ua(i) + u_mfactor
            va(i) = va(i) + v_mfactor
          elseif(spd(i).lt.5..and.spd(i).ge.0.)then
            ua(i) = ua(i) + u_lfactor
            va(i) = va(i) + v_lfactor
          endif
        else
          ua(i) = badpt
          va(i) = badpt
        endif
      enddo

      power_law = (10.0/(0.5*dz))**0.106
      ua = ua / power_law 
      va = va / power_law

!~site check
      print*,"site check"

      do k = 1,sea_info
        xxsta(k) = 1.+(lon(k) - lon_l)/lon_g
        yysta(k) = 1.+(lat(k) - lat_l)/lat_g
      enddo

      m = 0
      uu = badpt ; vv = badpt 
      open(58,file='stn_sea.dat',form='unformatted',&
              access='direct',recl=4*2)
        do k = 1,sea_info
          ii = nint(xxsta(k)) ; jj = nint(yysta(k))
          if(ii.ge.0.and.ii.le.NX-1.and.jj.ge.0.and.jj.le.NY-1)then
            m = m + 1
            uu(ii,jj) = ua(k) ; vv(ii,jj) = va(k)
            write(58,rec=m) lat(k),lon(k)
          endif
        enddo
      close(58)

      open(59,file='seasite.bin',form='unformatted'&
             ,access='direct',recl=4*NX*NY)
        write(59,rec=1)((uu(ii,jj),ii=1,NX),jj=1,NY)
        write(59,rec=2)((vv(ii,jj),ii=1,NX),jj=1,NY)
      close(59)

!~BG verify
      m = 0
      rmse_ub = 0.0 ; rmse_vb = 0.0
      do i = 1, NX
      do j = 1, NY
        if(uu(i,j).gt.badpt.and.vv(i,j).gt.badpt)then
          m = m + 1
          rmse_ub = rmse_ub + (uu(i,j) - u0(i,j,1))**2.
          rmse_vb = rmse_vb + (vv(i,j) - v0(i,j,1))**2.
        endif
      enddo
      enddo
      rmse_ub = sqrt(rmse_ub/float(m))
      rmse_vb = sqrt(rmse_vb/float(m))
      print*,"BG RMSE:",rmse_ub,rmse_vb

      m = 0
      fx = 0.0 ; sum_dir = 0.0
      do i = 1, nx
      do j = 1, ny
         if(uu(i,j).gt.badpt.and.vv(i,j).gt.badpt)then
            m = m + 1
           f1 = uu(i,j)*u0(i,j,1) + vv(i,j)*v0(i,j,1)
           f2 = sqrt(uu(i,j)**2.+vv(i,j)**2.)*sqrt(u0(i,j,1)**2.+v0(i,j,1)**2.)
           fx = acos(f1 / f2)*180./acos(-1.0)
           sum_dir = sum_dir + fx
         endif
      enddo
      enddo
      !print*,sum_dir/float(m)
     
      m = 0
      fx = 0.0 ; check_dir = 0.0
      do i = 1, nx
      do j = 1, ny
         if(uu(i,j).gt.badpt.and.vv(i,j).gt.badpt)then
           m = m + 1
           f1 = uu(i,j)*u0(i,j,1) + vv(i,j)*v0(i,j,1)
           f2 = sqrt(uu(i,j)**2.+vv(i,j)**2.)*sqrt(u0(i,j,1)**2.+v0(i,j,1)**2.)
           fx = acos(f1 / f2)*180./acos(-1.0)
           check_dir(m) = fx
         endif
      enddo
      enddo
     
      open(69,file='check_dir0.bin',form='unformatted'&
             ,access='direct',recl=4*102)
        write(69,rec=1)(check_dir(ii),ii=1,102)
      close(69)
!~read files-----------------------------------------------
      print*,"read analysis file "
      call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid, time, time_varid) )
      call check( nf90_get_var(ncid, time_varid, time0) )
      call check( nf90_inq_varid(ncid, x0, x0_varid) )
      call check( nf90_get_var(ncid, x0_varid, x) )
      call check( nf90_inq_varid(ncid, y0, y0_varid) )
      call check( nf90_get_var(ncid, y0_varid, y) )
      call check( nf90_inq_varid(ncid, z0, z0_varid) )
      call check( nf90_get_var(ncid, z0_varid, z) )
      call check( nf90_inq_varid(ncid, lat0, lat0_varid) )
      call check( nf90_get_var(ncid, lat0_varid, lat1) )
      call check( nf90_inq_varid(ncid, lon0, lon0_varid) )
      call check( nf90_get_var(ncid, lon0_varid, lon1) )
      !
      call check( nf90_inq_varid(ncid, ufct, ufct1_varid) )
      call check( nf90_get_var(ncid, ufct1_varid, u1) )
      call check( nf90_inq_varid(ncid, vfct, vfct1_varid) )
      call check( nf90_get_var(ncid, vfct1_varid, v1) )
      !
!~verify
      m = 0
      rmse_ua = 0.0 ; rmse_va = 0.0
      mvd = 0.0
      do i = 1, nx
      do j = 1, ny
        if(uu(i,j).gt.badpt.and.vv(i,j).gt.badpt)then
          m = m + 1
          rmse_ua = rmse_ua + (uu(i,j) - u1(i,j,1,1))**2.
          rmse_va = rmse_va + (vv(i,j) - v1(i,j,1,1))**2.
          mvd = mvd + sqrt((uu(i,j) - u1(i,j,1,1))**2. +&
                      (vv(i,j) - v1(i,j,1,1))**2.)
        endif
      enddo
      enddo
!
      rmse_ua = sqrt(rmse_ua/float(m))
      rmse_va = sqrt(rmse_va/float(m))
      mvd = mvd / float(m)
      print*,"analysis RMSE:",rmse_ua,rmse_va
      print*,"analysis mvd:",mvd


!~verify for wind dir 
     
      m = 0
      fx = 0.0 ; sum_dir = 0.0
      do i = 1, nx
      do j = 1, ny
         if(uu(i,j).gt.badpt.and.vv(i,j).gt.badpt)then
            m = m + 1
           f1 = uu(i,j)*u1(i,j,1,1) + vv(i,j)*v1(i,j,1,1)
           f2 = sqrt(uu(i,j)**2.+vv(i,j)**2.)*sqrt(u1(i,j,1,1)**2.+v1(i,j,1,1)**2.)
           fx = acos(f1 / f2)*180./acos(-1.0)
           sum_dir = sum_dir + fx
         endif
      enddo
      enddo
      print*,m 
     
      m = 0
      fx = 0.0 ; check_dir = 0.0
      do i = 1, nx
      do j = 1, ny
         if(uu(i,j).gt.badpt.and.vv(i,j).gt.badpt)then
           m = m + 1
           f1 = uu(i,j)*u1(i,j,1,1) + vv(i,j)*v1(i,j,1,1)
           f2 = sqrt(uu(i,j)**2.+vv(i,j)**2.)*sqrt(u1(i,j,1,1)**2.+v1(i,j,1,1)**2.)
           fx = acos(f1 / f2)*180./acos(-1.0)
           check_dir(m) = fx
         endif
      enddo
      enddo
      
      open(59,file='check_dir.bin',form='unformatted'&
             ,access='direct',recl=4*102)
        write(59,rec=1)(check_dir(ii),ii=1,102)
      close(59)


      stop
      end
!===========================================================
      subroutine check(status)
      use netcdf , only : nf90_noerr, nf90_strerror
      implicit none
      integer, intent(in) :: status
      !
      if(status /= nf90_noerr) then
         print *,trim(nf90_strerror(status))
         stop "Stopped"
      endif
      end subroutine  check
