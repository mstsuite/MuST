c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zeropt(ezpt,tpzpt,omegws,ztotss)
c     ================================================================
c
      implicit   none
c
      integer    idebye(49)
      integer    iz
c
      real*8     ezpt
      real*8     tpzpt
      real*8     omegws
      real*8     ztotss
      real*8     bolts
      real*8     grune(49)
      real*8     expvol(49)
      real*8     ezero
      real*8     tpvzer
c
      data       bolts/6.33870d-06/
      data       idebye/0,0,344,1440,0,0,0,0,0,
     >                  0,158,400,428,0,0,0,0,0,91,
     >                  230,360,420,380,630,410,467,445,450,343,
     >                  327,320,0,0,0,0,0,56,147,280,
     >                  291,275,450,350,600,480,274,225,209,108/
      data       grune/0.0d0,0.0d0,1.18d0,1.18d0,0.00d0,0.00d0,0.00d0,
     >                 0.0d0,0.0d0,0.00d0,1.31d0,1.48d0,2.19d0,0.00d0,
     >                 0.0d0,0.0d0,0.00d0,0.00d0,1.37d0,1.16d0,1.17d0,
     >                 1.18d0,1.05d0,1.30d0,2.07d0,1.66d0,
     >                 1.93d0,1.88d0,2.00d0,2.01d0,2.00d0,
     >                 0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,1.67d0,1.00d0,
     >                 0.89d0,0.83d0,1.58d0,1.60d0,2.60d0,3.20d0,2.23d0,
     >                 2.28d0,2.36d0,2.23d0,2.37d0/
      data       expvol/0.0d0,143.7d0,0.0d0,54.54d0,0.0d0,0.0d0,0.0d0,
     >                  0.0d0,0.0d0,0.0d0,254.5d0,151.4d0,109.9d0,0.0d0,
     >                  0.0d0,0.0d0,0.0d0,  0.0d0,481.3d0,291.1d0,
     >                  168.7d0,120.3d0,93.48d0,80.63d0,82.84d0,78.95d0,
     >                  74.72d0,73.42d0,78.92d0,99.35d0,132.4d0,
     >                  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,598.9d0,373.6d0,
     >                  194.7d0,139.9d0,119.2d0,102.7d0,97.25d0,92.54d0,
     >                  93.70d0,99.67d0,111.9d0,142.9d0,179.2d0/
c
c     ================================================================
      iz=ztotss+.1
c     write(6,'('' zeropt:: iz,deb,grun,expv,omeg:'',2i5,3d12.4)')
c    >                      iz,idebye(iz),grune(iz),expvol(iz),omegws
      if(iz.lt.1 .or. iz.gt.49) then
         ezero=0.0d0
         tpvzer=0.0d0
      else
         if(expvol(iz).ne.0.0d0) then
            ezero=1.125d0*bolts*idebye(iz)*
     >            (expvol(iz)/omegws)**grune(iz)
         endif
         tpvzer =3.0d0*grune(iz)*ezero
      endif
      ezpt=ezero
      tpzpt=tpvzer
c     ================================================================
      return
      end
