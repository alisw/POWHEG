c -*- Fortran -*-

c The user must set nlegborn to the appropriate value for his process.
      integer nlegborn,nlegreal
      
      parameter (nlegborn=4)
      parameter (nlegreal=nlegborn+1)

c     ndiminteg is the dimensionality of the full real integral
c     ndiminteg=(nlegreal-2)*3-4+2-1
c     if there are undecayed resonances, we need extra variables to pilot
c     the resonance's masses

      integer ndiminteg
      parameter (ndiminteg=(nlegreal-2)*3-4+2-1
     .    + 0 )  ! 0=no resonance, 1=1 resonance

      integer maxprocborn,maxprocreal
      parameter (maxprocborn=283,maxprocreal=246)

      integer maxalr
c     qqbAgg: 3 * 6 + 2 * 6 (QED)
c     qbqAgg: 3 * 6 + 2 * 6 (QED)
c     qqbAq'q'b: 1 * 6 * 5 + 3 * 6 + 4 * 6 (QED)
c     qbqAq'q'b: 1 * 6 * 5 + 3 * 6 + 4 * 6 (QED)
c     qq'Aqq': 2 * 6 * 5 + 2 * 6 + 4 * 6 (QED)
c     qbq'bAqbq'b: 2 * 6 * 5 + 2 * 6 + 4 * 6 (QED)
c     qq'bAqq'b: 2 * 6 * 5 + 4 * 6 (QED)
c     qbq'Aq'qb: 2 * 6 * 5 + 4 * 6 (QED)
c     qgAqg: 3 * 6 + 2 * 6 (QED)
c     gqAqg: 3 * 6 + 2 * 6 (QED)
c     qbgAqbg: 3 * 6 + 2 * 6 (QED)
c     qbgAqbg: 3 * 6 + 2 * 6 (QED)
c     sum = 468 + 216 (QED) = 684
      parameter (maxalr=1250) ! dunno...
