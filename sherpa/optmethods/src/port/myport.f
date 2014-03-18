      subroutine mydn2fb(m,n,x,bounds,rfctol,afctol,partol,
     +     maxfev,nfev,fval,ierr,iprint,calcf,lv,liv,
     +     scalefactor,uiparm,urparm,v,iv)
      implicit none
      integer lv, liv
      integer m, n, maxfev, nfev, ierr, iprint
      double precision x(n), bounds(2,n)
      double precision rfctol, afctol, partol, fval
      double precision urparm(n), v(lv), scalefactor(n)
      integer uiparm(n), iv( liv )
      external calcf, ufparm, dpptri

      integer ii

      do ii = 1, n
         if ( x( ii ) .ne. 0.0d0 ) then
            scalefactor( ii ) = 1.0d0 / x( ii )
         else
            scalefactor( ii ) = 1.0d6
         endif
      end do

C  ***  ALG = 1 MEANS REGRESSION CONSTANTS.
      call divset( 1, iv, liv, lv, v )
C     MXFCAL=17
      iv( 17 ) = maxfev
C     PRUNIT=21
      iv( 21 ) = iprint
C     AFCTOL=31
      v( 31 ) = afctol
C     RFCTOL=32
      v( 32 ) = rfctol
C     XCTOL=33
      v( 33 ) = partol

      call dn2fb(m,n,x,bounds,calcf,iv,liv,lv,v,uiparm,urparm,ufparm)

      ierr = iv(1)
c     NFCALL=6, NGCALL=30
      nfev = iv( 6 ) + iv( 30 )
      fval = v( 10 )

      return
      end

      subroutine mydmnfb(n,x,bounds,rfctol,afctol,partol,
     +     maxfev,nfev,fval,ierr,iprint,calcf,lv,liv,
     +     scalefactor,uiparm,urparm,v,iv,errdpptri,covarerr)
      implicit none
      integer lv, liv, errdpptri
      integer n, maxfev, nfev, ierr, iprint
      double precision x(n), bounds(2,n)
      double precision rfctol, afctol, partol, fval
      double precision urparm(n), v(lv), scalefactor(n)
      double precision covarerr(n)
      integer uiparm(n), iv( liv )
      external calcf, ufparm

      integer ii

      do ii = 1, n
         if ( x( ii ) .ne. 0.0d0 ) then
            scalefactor( ii ) = 1.0d0 / x( ii )
         else
            scalefactor( ii ) = 1.0d6
         endif
      end do
      
C  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS.
      call divset( 2, iv, liv, lv, v )
C     MXFCAL=17
      iv( 17 ) = maxfev
C     PRUNIT=21
      iv( 21 ) = iprint
C     AFCTOL=31
      v( 31 ) = afctol
C     RFCTOL=32
      v( 32 ) = rfctol
C     XCTOL=33
      v( 33 ) = partol

      call dmnfb(n,scalefactor,x,bounds,
     +     calcf,iv,liv,lv,v,uiparm,urparm,ufparm)


      ierr = iv(1)
c     NFCALL=6, NGCALL=30
      nfev = iv( 6 ) + iv( 30 )
      fval = v( 10 )

      call dpptri( 'U', n, v( iv(42) ), errdpptri )
      if ( errdpptri .eq. 0 ) then
         do ii = 1, n
            covarerr( ii ) = dsqrt( v( iv(42) + ii*(ii+1)/2 - 1) )
         end do
      end if

      return
      end

      subroutine ufparm( n, x, uiparm, urparm )
      implicit none
      integer n, uiparm(*)
      double precision x(*), urparm(*)
      print *, '!!!!!!!!!!!!!!!!!!ufparm!!!!!!!!!!!!!!!!!'
      print *, '!!!!!Warning: Not supposed to be here!!!!!'
      print *, '!!!!!!!!!!!!!!!!!!ufparm!!!!!!!!!!!!!!!!!'
      return
      end
