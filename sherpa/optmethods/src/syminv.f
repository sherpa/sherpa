      subroutine syminv(a,n,c,w,nullty,ifault,rmax)
c
c     algorithm as7, applied statistics, vol.17, 1968.
c
c     arguments:-
c     a()     = input, the symmetric matrix to be inverted, stored in
c               lower triangular form
c     n       = input, order of the matrix
c     c()     = output, the inverse of a (a generalized inverse if c is
c               singular), also stored in lower triangular.
c               c and a may occupy the same locations.
c     w()     = workspace, dimension at least n.
c     nullty  = output, the rank deficiency of a.
c     ifault  = output, error indicator
c                     = 1 if n < 1
c                     = 2 if a is not +ve semi-definite
c                     = 0 otherwise
c     rmax    = output, approximate bound on the accuracy of the diagonal
c               elements of c.  e.g. if rmax = 1.e-04 then the diagonal
c               elements of c will be accurate to about 4 dec. digits.
c
c     latest revision - 18 october 1985
c
c*************************************************************************
c
      integer*4 ifault,n,nullty
      real*8 a(*),c(*),rmax,w(n)

      integer*4 i,icol,irow,j,jcol,k,l,mdiag,ndiag,nn,nrow
      real*8 one,x,zero
      data zero/0.d0/, one/1.d0/
c
      nrow=n
      ifault=1
      if(nrow.le.0) go to 100
      ifault=0
c
c     cholesky factorization of a, result in c
c
      call chola(a,nrow,c,nullty,ifault,rmax,w)
      if(ifault.ne.0) go to 100
c
c     invert c & form the product (cinv)'*cinv, where cinv is the inverse
c     of c, row by row starting with the last row.
c     irow = the row number, ndiag = location of last element in the row.
c
      nn=nrow*(nrow+1)/2
      irow=nrow
      ndiag=nn
   10 if(c(ndiag).eq.zero) go to 60
      l=ndiag
      do 20 i=irow,nrow
        w(i)=c(l)
        l=l+i
   20 continue
      icol=nrow
      jcol=nn
      mdiag=nn
   30 l=jcol
      x=zero
      if(icol.eq.irow) x=one/w(irow)
      k=nrow
   40 if(k.eq.irow) go to 50
      x=x-w(k)*c(l)
      k=k-1
      l=l-1
      if(l.gt.mdiag) l=l-k+1
      go to 40
   50 c(l)=x/w(irow)
      if(icol.eq.irow) go to 80
      mdiag=mdiag-icol
      icol=icol-1
      jcol=jcol-1
      go to 30
c
c     special case, zero diagonal element.
c
   60 l=ndiag
      do 70 j=irow,nrow
        c(l)=zero
        l=l+j
   70 continue
c
c      end of row.
c
   80 ndiag=ndiag-irow
      irow=irow-1
      if(irow.ne.0) go to 10
  100 return
      end





      subroutine chola(a, n, u, nullty, ifault, rmax, r)
c
c     algorithm as6, applied statistics, vol.17, 1968, with
c     modifications by a.j.miller
c
c     arguments:-
c     a()     = input, a +ve definite matrix stored in lower-triangular
c               form.
c     n       = input, the order of a
c     u()     = output, a lower triangular matrix such that u*u' = a.
c               a & u may occupy the same locations.
c     nullty  = output, the rank deficiency of a.
c     ifault  = output, error indicator
c                     = 1 if n < 1
c                     = 2 if a is not +ve semi-definite
c                     = 0 otherwise
c     rmax    = output, an estimate of the relative accuracy of the
c               diagonal elements of u.
c     r()     = output, array containing bounds on the relative accuracy
c               of each diagonal element of u.
c
c     latest revision - 18 october 1985
c
c*************************************************************************
c
      integer*4 ifault,n,nullty
      real*8 a(*),u(*),r(n),rmax
c
c     eta should be set equal to the smallest +ve value such that
c     1.0 + eta is calculated as being greater than 1.0 in the accuracy
c     being used.
c
      integer*4 i,icol,irow,j,k,l,m
      real*8 eta,five,rsq,w,zero
      data eta/1.d-16/, zero/0.d0/, five/5.d0/

      i=0
      rsq=0.0d0
      w=0.0d0
c
      ifault=1
      if(n.le.0) go to 100
      ifault=2
      nullty=0
      rmax=eta
      r(1)=eta
      j=1
      k=0
c
c     factorize column by column, icol = column no.
c
      do 80 icol=1,n
        l=0
c
c     irow = row number within column icol
c
        do 40 irow=1,icol
          k=k+1
          w=a(k)
          if(irow.eq.icol) rsq=(w*eta)**2
          m=j
          do 10 i=1,irow
            l=l+1
            if(i.eq.irow) go to 20
            w=w-u(l)*u(m)
            if(irow.eq.icol) rsq=rsq+(u(l)**2*r(i))**2
            m=m+1
   10     continue
   20     if(irow.eq.icol) go to 50
          if(u(l).eq.zero) go to 30
          u(k)=w/u(l)
          go to 40
   30     u(k)=zero
          if(abs(w).gt.abs(rmax*a(k))) go to 100
   40   continue
c
c     end of row, estimate relative accuracy of diagonal element.
c
   50   rsq=sqrt(rsq)
        if(abs(w).le.five*rsq) go to 60
        if(w.lt.zero) go to 100
        u(k)=sqrt(w)
        r(i)=rsq/w
        if(r(i).gt.rmax) rmax=r(i)
        go to 70
   60   u(k)=zero
        nullty=nullty+1
   70   j=j+icol
   80 continue
      ifault=0
c
  100 return
      end

