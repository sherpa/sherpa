#ifndef syminv_hh
#define syminv_hh

#include <cmath>

namespace appliedstats {

  template< typename Real >
  void chola( std::vector<Real>& a, int n, std::vector<Real>& u, int& nullty,
	      int& ifault, Real& rmax, std::vector<Real>& r ) {
    //
    //      subroutine chola(a, n, u, nullty, ifault, rmax, r)
    //c
    //c     algorithm as6, applied statistics, vol.17, 1968, with
    //c     modifications by a.j.miller
    //c
    //c     arguments:-
    //c     a()     = input, a +ve definite matrix stored in lower-triangular
    //c               form.
    //c     n       = input, the order of a
    //c     u()     = output, a lower triangular matrix such that u*u' = a.
    //c               a & u may occupy the same locations.
    //c     nullty  = output, the rank deficiency of a.
    //c     ifault  = output, error indicator
    //c                     = 1 if n < 1
    //c                     = 2 if a is not +ve semi-definite
    //c                     = 0 otherwise
    //c     rmax    = output, an estimate of the relative accuracy of the
    //c               diagonal elements of u.
    //c     r()     = output, array containing bounds on the relative accuracy
    //c               of each diagonal element of u.
    //c
    //c     latest revision - 18 october 1985
    //c
    //c************************************************************************
    //c
    //      integer*4 ifault,n,nullty
    //      real*8 a(*),u(*),r(n),rmax
    //c
    //c     eta should be set equal to the smallest +ve value such that
    //c     1.0 + eta is calculated as being greater than 1.0 in the accuracy
    //c     being used.
    //c
    //      integer*4 i,icol,irow,j,k,l,m
    //      real*8 eta,five,rsq,w,zero
    //      data eta/1.d-16/, zero/0.d0/, five/5.d0/
    const Real eta = std::numeric_limits< Real >::epsilon();
    const Real zero = 0.0;
    const Real five = 5.0;

    int i__ = 0;
    Real rsq = 0.;
    Real w = 0.;

    //c
    //      ifault=1
    //      if(n.le.0) go to 100
    ifault = 1;
    if ( n <= 0 )
      goto L100;

    //      ifault=2
    //      nullty=0
    //      rmax=eta
    //      r(1)=eta
    //      j=1
    //      k=0
    ifault = 2;
    nullty = 0;
    rmax = eta;
    r[0] = eta;
    int j, k;
    j = 1;
    k = 0;


    //c
    //c     factorize column by column, icol = column no.
    //c
    //      do 80 icol=1,n
    //        l=0
    for ( int icol = 1; icol <= n; ++icol) {
      int l = 0;

      //c
      //c     irow = row number within column icol
      //c
      //        do 40 irow=1,icol
      //          k=k+1
      //          w=a(k)
      //          if(irow.eq.icol) rsq=(w*eta)**2
      //          m=j
      for ( int irow = 1; irow <= icol; ++irow) {
	++k;
	w = a[k - 1];
	if (irow == icol)
	  rsq = pow( w * eta, 2.0 );
	int m = j;


	//          do 10 i=1,irow
	//            l=l+1
	//            if(i.eq.irow) go to 20
	//            w=w-u(l)*u(m)
	//            if(irow.eq.icol) rsq=rsq+(u(l)**2*r(i))**2
	//            m=m+1
	//   10     continue
	for (i__ = 1; i__ <= irow; ++i__) {
	  ++l;
	  if (i__ == irow)
	    goto L20;
	  w -= u[l - 1] * u[m - 1];
	  if (irow == icol)
	    rsq += pow( u[l - 1]*u[l - 1] * r[i__ - 1], 2.0 );
	  ++m;
	  /* L10: */
	}

	
	//   20     if(irow.eq.icol) go to 50
	//          if(u(l).eq.zero) go to 30
      L20:
	if (irow == icol)
	  goto L50;
	if (u[l - 1] == zero)
	  goto L30;

	
	//          u(k)=w/u(l)
	//          go to 40
	//   30     u(k)=zero
	//          if(abs(w).gt.abs(rmax*a(k))) go to 100
	//   40   continue
	u[k - 1] = w / u[l - 1];
	goto L40;
      L30:
	u[k - 1] = zero;
	if ( fabs( w ) > fabs( rmax * a[k - 1] ) )
	  goto L100;
      L40:
	;
      }                            // for ( int irow = 0; irow < icol; ++irow )


    //c
    //c     end of row, estimate relative accuracy of diagonal element.
    //c
    //   50   rsq=sqrt(rsq)
    //        if(abs(w).le.five*rsq) go to 60
    //        if(w.lt.zero) go to 100
    //        u(k)=sqrt(w)
    //        r(i)=rsq/w
    //        if(r(i).gt.rmax) rmax=r(i)
    //        go to 70
    L50:
      rsq = sqrt(rsq);
      if (fabs(w) <= five * rsq)
	goto L60;
      if (w < zero)
	goto L100;
      u[k - 1] = sqrt(w);
      r[i__ - 1] = rsq / w;
      if (r[i__ - 1] > rmax)
	rmax = r[i__ - 1];
      goto L70;

      //   60   u(k)=zero
      //        nullty=nullty+1
    L60:
      u[k - 1] = zero;
      ++nullty;
    L70:
      j += icol;
      /* L80: */
    }
    ifault = 0;

  L100:
    return;
  } /* chola */


  template< typename Real >
  void syminv(std::vector<Real>& a, int n, std::vector<Real>& c, 
	      std::vector<Real>& w, int& nullty, int& ifault, Real& rmax) {
    //
    //      subroutine syminv(a,n,c,w,nullty,ifault,rmax)
    //c
    //c     algorithm as7, applied statistics, vol.17, 1968.
    //c
    //c     arguments:-
    //c     a()     = input, the symmetric matrix to be inverted, stored in
    //c               lower triangular form
    //c     n       = input, order of the matrix
    //c     c()     = output, the inverse of a (a generalized inverse if c is
    //c               singular), also stored in lower triangular.
    //c               c and a may occupy the same locations.
    //c     w()     = workspace, dimension at least n.
    //c     nullty  = output, the rank deficiency of a.
    //c     ifault  = output, error indicator
    //c                     = 1 if n < 1
    //c                     = 2 if a is not +ve semi-definite
    //c                     = 0 otherwise
    //c     rmax    = output, approximate bound on the accuracy of the diagonal
    //c               elements of c.  e.g. if rmax = 1.e-04 then the diagonal
    //c               elements of c will be accurate to about 4 dec. digits.
    //c
    //c     latest revision - 18 october 1985
    //c
    //c************************************************************************
    //c
    //      integer*4 ifault,n,nullty
    //      real*8 a(*),c(*),rmax,w(n)
    //
    //      integer*4 i,icol,irow,j,jcol,k,l,mdiag,ndiag,nn,nrow
    //      real*8 one,x,zero
    //      data zero/0.d0/, one/1.d0/
    //c
    const Real zero = 0.;
    const Real one = 1.;

    //      integer*4 i,icol,irow,j,jcol,k,l,mdiag,ndiag,nn,nrow
    //      real*8 one,x,zero
    //      data zero/0.d0/, one/1.d0/
    //c
    int j, k, l;
    Real x;
    int nn, icol, jcol, irow, nrow, mdiag, ndiag;

    nrow = n;
    ifault = 1;
    if (nrow <= 0) {
      goto L100;
    }
    ifault = 0;

    //c
    //c     cholesky factorization of a, result in c
    //c
    //      call chola(a,nrow,c,nullty,ifault,rmax,w)
    //      if(ifault.ne.0) go to 100
    chola(a, nrow, c, nullty, ifault, rmax, w);
    if ( ifault != 0 )
      goto L100;

    //c
    //c     invert c & form the product (cinv)'*cinv, where cinv is the inverse
    //c     of c, row by row starting with the last row.
    //c     irow = the row number, ndiag = location of last element in the row.
    //c
    //      nn=nrow*(nrow+1)/2
    //      irow=nrow
    //      ndiag=nn
    nn = nrow * (nrow + 1) / 2;
    irow = nrow;
    ndiag = nn;
    
    //   10 if(c(ndiag).eq.zero) go to 60
    //      l=ndiag
    //      do 20 i=irow,nrow
    //        w(i)=c(l)
    //        l=l+i
  L10:
    if (c[ndiag - 1] == zero) {
      goto L60;
    }
    l = ndiag;
    for ( int i__ = irow; i__ <= nrow; ++i__) {
      w[i__ - 1] = c[l - 1];
      l += i__;
      /* L20: */
    }

    //   20 continue
    //      icol=nrow
    //      jcol=nn
    //      mdiag=nn
    icol = nrow;
    jcol = nn;
    mdiag = nn;

    //   30 l=jcol
    //      x=zero
    //      if(icol.eq.irow) x=one/w(irow)
    //      k=nrow
  L30:
    l = jcol;
    x = zero;
    if ( icol == irow )
      x = one / w[irow - 1];
    k = nrow;

    //   40 if(k.eq.irow) go to 50
    //      x=x-w(k)*c(l)
    //      k=k-1
    //      l=l-1
    //      if(l.gt.mdiag) l=l-k+1
    //      go to 40
  L40:
    if ( k == irow )
      goto L50;
    x -= w[k - 1] * c[l - 1];
    --k;
    --l;
    if (l > mdiag) {
      l = l - k + 1;
    }
    goto L40;

    //   50 c(l)=x/w(irow)
    //      if(icol.eq.irow) go to 80
    //      mdiag=mdiag-icol
    //      icol=icol-1
    //      jcol=jcol-1
    //      go to 30
    
  L50:
    c[l - 1] = x / w[irow - 1];
    if (icol == irow) {
      goto L80;
    }
    mdiag -= icol;
    --icol;
    --jcol;
    goto L30;


    //c
    //c     special case, zero diagonal element.
    //c
    //   60 l=ndiag
    //      do 70 j=irow,nrow
    //        c(l)=zero
    //        l=l+j
    //   70 continue
  L60:
    l = ndiag;
    for (j = irow; j <= nrow; ++j) {
      c[l - 1] = zero;
      l += j;
      /* L70: */
    }
    //c
    //c      end of row.
    //c

    //   80 ndiag=ndiag-irow
    //      irow=irow-1
    //      if(irow.ne.0) go to 10
    //  100 return
    //      end
  L80:
    ndiag -= irow;
    --irow;
    if (irow != 0) {
      goto L10;
    }
  L100:
    return;

  }                                                                   // syminv

}                                                     // namespace appliedstats
#endif
