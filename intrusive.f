c intrusive light calculation program
c
c gfortran intrusive.f intrants2D.f extrants2d.f -o intrusive
c 
      program intrusive
      integer width                                                       ! Matrix dimension in Length/width and height
      parameter (width=1024)

      real Irad(width,width)                                              ! intrusive radiance
      real h_w                                                            ! height of the center of a window relative to the ground
      real pvalno(181),pval(181)                                          ! Light output pattern
      real srei(width,width)                                              ! Modis reflectance
      real d_o(width,width)                                               ! avg distance between obstacles
      real h_o(width,width)                                               ! avg height of obstacles relative to the ground
      real h_l(width,width)                                               ! avg lamp height relative to the ground
      real intlu(width,width)                                             ! lamp flux
      real val2d(width,width)
      real z,pi,dz
      real inteo,integ,pvalto
      real Iradmax,gain,offset,xcell0,ycell0,pixsiz
      integer stype,i,j,k,iw,nwav,nw,nzon,nz,lenbase,valmax
      integer nbx,nby
      character*3 zon(120),wav(400)
      character*72 basenm,ohfile,odfile,lhfile,rfile,lfile,lopf,intrufi
      character*72 pafile
      character*12 nom
      pi=3.14159
      dz=pi/180.
      print*,'Name of the experiment?'
      read*,basenm
      print*,'Height to the center of a window from the ground (m)?'
      read*,h_w
c
c  determine the Length of basenm
c 
      lenbase=index(basenm,' ')-1  
c   reading wav.lst
c
      nwav=0
      open(unit=1,file='wav.lst',status='old')
          do k=1,400
             read(1,*,end=10) wav(k)
             nwav=nwav+1
          enddo
 10   close(unit=1)
c   reading zon.lst
c
      nzon=0
      open(unit=1,file='zon.lst',status='old')
          do k=1,120
             read(1,*,end=20) zon(k)
             nzon=nzon+1
          enddo
 20   close(unit=1)
c reading obstacle heights
        ohfile=basenm(1:lenbase)//'_obsth.pgm'
        call intrants2d(ohfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             h_o(i,j)=val2d(i,j)                                          ! Filling of the array
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c reading obstacle distances
        odfile=basenm(1:lenbase)//'_obstd.pgm'
        call intrants2d(odfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             d_o(i,j)=val2d(i,j)                                          ! Filling of the array
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c readind lamp heights
        lhfile=basenm(1:lenbase)//'_altlp.pgm'
        call intrants2d(lhfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             h_l(i,j)=val2d(i,j)                                          ! Filling of the array 
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c
c 
      do nw=1,nwav
        do i=1,nbx
          do j=1,nby
            Irad(i,j)=0.
          enddo
        enddo
        Iradmax=0.
c reading modis reflectances
        rfile='modis_'//wav(nw)//'.pgm'
        call intrants2d(rfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            srei(i,j)=val2d(i,j)                                          ! Filling of the array 
          enddo                                                           ! end of the loop over all cells along y.
        enddo      

        do nz=1,nzon
c reading LOPs
          lopf='fctem_wl_'//wav(nw)//'_zon_'//zon(nz)//'.dat'
          pvalto=0.
          open(UNIT=7, FILE=lopf,status='OLD')                            ! opening file pa#.dat, angular photometry.
            do k=1,181                                                    ! beginning of the loop for the 181 data points
              read(7,*) pval(k)                                           ! reading of the data in the array pval.
              pvalto=pvalto+pval(k)*2.*pi*                                ! Sum of the values of the  photometric function 
     a        sin(real(k-1)*dz)*dz                                        ! (pvaleur x 2pi x sin z x z) (ou z egale 
c                                                                         ! (i-1) x 1 degrees).
            enddo                                                         ! end of the loop over the 181 donnees du fichier pa#.dat.
          close(7)                                                        ! closing file pa#.dat, angular photometry.
          do k=1,181
            if (pvalto.ne.0.) pvalno(k)=pval(k)/pvalto                    ! Normalisation of the photometric function.
          enddo   

c reading du fluxes
          lfile=basenm(1:lenbase)//'_'//wav(nw)//'_lumlp_'//zon(nz)//
     +'.pgm'
          call intrants2d(lfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
          do i=1,nbx                                                      ! beginning of the loop over all cells along x.
            do j=1,nby                                                    ! beginning of the loop over all cells along y.
              intlu(i,j)=val2d(i,j)   
              if (intlu(i,j).ne.0.) print*,intlu(i,j)                                    ! Filling of the array 
            enddo                                                         ! end of the loop over all cells along y.
          enddo      

          do i=1,nbx
            do j=1,nby
c calculate the basic angles of the geometry
              z_o=pi/2.-atan((h_o(i,j)-h_l(i,j))/d_o(i,j)/2.)
              z_g=atan(d_o(i,j)/2./h_l(i,j))
              z_w=pi/2.+atan((h_l(i,j)-h_w)/d_o(i,j)/2.)

c integrate the LOP from obstacle base to obstacle top (inteo)
c and LOP from nadir to obstacle base (integ)
              inteo=0.
              integ=0.
              do k=1,181
                z=real(k-1)*pi/180.
                if ((z.ge.z_o).and.(z.lt.z_g)) then
                   inteo=inteo+pvalno(k)*sin(z)*dz
                endif
                if ((z.ge.z_g).and.(z.lt.pi)) then
                   integ=integ+pvalno(k)*sin(z)*dz
                endif
                if (z.eq.z_w) iw=k
              enddo
              Irad(i,j)=Irad(i,j)+intlu(i,j)*(pvalno(iw)+2.*srei(i,j)*
     +        integ+srei(i,j)/4.*inteo)
              if (Irad(i,j).gt.Iradmax) then
                print*,'entre',Iradmax,Irad(i,j),i,j
                 Iradmax=Irad(i,j)
              endif
            enddo
          enddo
        enddo                                                             ! end of loop over zones
c writing the instrusive light map for each wavelength
c
        intrufi=basenm(1:lenbase)//'_'//wav(nw)//'_intrus.pgm'
        nom='Intrusive '
        valmax=65535       
        gain=Iradmax/real(valmax)
        offset=0.
        call extrants2d (intrufi,Irad,nom,xcell0,ycell0,pixsiz,
     +  gain,offset,nbx,nby,valmax)
      enddo                                                               ! end of loop over wavelength
      stop
      end
