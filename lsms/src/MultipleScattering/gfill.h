      integer kmymaxp
      integer iplmax
      
      parameter (iplmax=6, kmymaxp=2*(iplmax+1)*(iplmax+1))

      complex*16 sxcoeff(kmymaxp,kmymaxp),sxbcoeff(kmymaxp,kmymaxp)
      complex*16 sycoeff(kmymaxp,kmymaxp),sybcoeff(kmymaxp,kmymaxp)
      complex*16 szcoeff(kmymaxp,kmymaxp),szbcoeff(kmymaxp,kmymaxp)
      common/sigmat/sxcoeff,sxbcoeff,sycoeff,sybcoeff,szcoeff,szbcoeff

      complex*16 lxcoeff(kmymaxp,kmymaxp),lxbcoeff(kmymaxp,kmymaxp)
      complex*16 lycoeff(kmymaxp,kmymaxp),lybcoeff(kmymaxp,kmymaxp)
      complex*16 lzcoeff(kmymaxp,kmymaxp),lzbcoeff(kmymaxp,kmymaxp)
      common/lmat/lxcoeff,lxbcoeff,lycoeff,lybcoeff,lzcoeff,lzbcoeff

      integer lamp,lammp
      parameter ( lamp=1,lammp=(lamp+1)*(lamp+1) )

      
      complex*16 gacoeff((iplmax+1)*(iplmax+1),(iplmax+1)*(iplmax+1),
     >                   lammp)
      complex*16 rgacoeff(kmymaxp,kmymaxp,lammp)

      common/ggaunt/gacoeff
      common/rggaunt/rgacoeff
