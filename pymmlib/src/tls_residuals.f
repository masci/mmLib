C Copyright 2009-2010 by PyMMLib Development Group (see AUTHORS file)
C This code is part of the PyMMLib distribution and governed by
C its license.  Please see the LICENSE file that should have been
C included as part of this package.
C
C DESCRIPTION: tls_residuals - Fortran77 subroutines and functions
C used for validation PDB structures.
C AUTHORS:
C  - P. Christoph Champ (PCC),
C  - Frank Zucker,
C  - Ethan A. Merritt (EAM)
C REQUIREMENTS: f2py (part of the Numpy package)
C
************************************************************************
* UPDATES
* PCC Sep 2009  - added sum_square_diff() + diff_trace_uv()
* PCC Aug 2009  - initial version
*
************************************************************************

C     ******************************************************************
      SUBROUTINE calc_ccuij(corwij,U,V)
C     ******************************************************************
C --- Calculates the correlation coefficient of Uij and Vij
c      implicit none
      integer MAXATOM
      parameter (MAXATOM=10000)
      real inv3x3
      real U(3,3),  UINV(3,3)
      real V(3,3),  VINV(3,3)
      real W(3,3),  WINV(3,3)
      real          ZINV(3,3)
      real VV(3,3), VVINV(3,3)
      real det, det_U, det_V, det_W, det_VV
      real corwij
      real uijcor(MAXATOM)
Cf2py intent(in) u, v
Cf2py intent(out) corwij
Cf2py depend(u, v) corwij

c Take inverses of the two matrices
      det_U = inv3x3(Uinv, U)
      det_V = inv3x3(Vinv, V)

      if (det_U .le. 0 .or. det_V .le. 0 .or. 
     *     U(1,1) .le. 0 .or. V(1,1) .le. 0 .or.
     *     U(1,1)*U(2,2) .le. U(1,2)*U(2,1) .or.
     *     V(1,1)*V(2,2) .le. V(1,2)*V(2,1)) then
          corwij = 0
          return
      endif

c And their sum (for the covariance)
      call add3x3( Winv, Uinv, Vinv )
      det_W = 1. / inv3x3( W, Winv )

      if (det_W .le. 0) then
          corwij = 0
          return
      endif

c Now the correlation coefficient
      corwij = (8.0 * det_W) / sqrt( det_U * det_V )
      corwij = sqrt(corwij)
c     if (corwij .gt. 1.00001) write (0,*) ' BOGUS CCuij = ',corwij

      return
      end

C     ******************************************************************
      SUBROUTINE sum_square_diff(sqd,U,V)
C     ******************************************************************
C --- Calculates the sum of the differences of anisotropic ADP tensors
C --- U and V squared.
C --- sqd = abs(numpy.sum(numpy.subtract(U,V)**2))
      real U(3,3)
      real V(3,3)
      real W(3,3)
      real sqd, WSUM
Cf2py intent(in) u, v
Cf2py intent(out) sqd
Cf2py depend(u, v) sqd

      call subtract3x3(W, U, V)
      call sum3x3(WSUM, W**2)

      sqd = abs(WSUM)
      return
      end

C     ******************************************************************
      SUBROUTINE diff_trace_uv(dtuv,U,V)
C     ******************************************************************
C --- Calculates the trace difference of anisotropic ADP tensors U and V.
C --- dtuv = abs((numpy.trace(U) - numpy.trace(V))/ 3.0)
      real U(3,3)
      real V(3,3)
      real W(3,3)
      real dtuv, Utrace, Vtrace
Cf2py intent(in) u, v
Cf2py intent(out) dtuv
Cf2py depend(u, v) dtuv

      call trace3x3(Utrace, U)
      call trace3x3(Vtrace, V)

C      write(*,*) 'U-V = ', (Utrace - Vtrace)

      dtuv = abs((Utrace - Vtrace) / 3.0)
      return
      end

C     ******************************************************************
      SUBROUTINE rosenfeld(r,Upos,Vpos,U,V)
C     ******************************************************************
C --- Calculates the "Rosenfeld" of anisotropic ADP tensors U and V.
C --- r = abs(Un - Vn)
C
C     a/b = x,y,z of U/b
C     d   = distance between a, b
C
C     n = numpy.array([(a[0] - b[0])/d, (a[1] - b[1])/d, (a[2] - b[2])/d])
C     Un = numpy.dot(numpy.dot(n, U), numpy.transpose(n))
C     Vn = numpy.dot(numpy.dot(n, V), numpy.transpose(n))
C     return abs(Un - Vn)
C----------------------
      real U(3,3), Upos(3), Un, U13(1,3)
      real V(3,3), Vpos(3), Vn, V13(1,3)
      real d, n(1,3), nT(3,1), r
Cf2py intent(in) upos, vpos, u, v
Cf2py intent(out) r
Cf2py depend(u, v) r

      d = sqrt( (Upos(1)-Vpos(1))**2 +
     &          (Upos(2)-Vpos(2))**2 +
     &          (Upos(3)-Vpos(3))**2 )

      do j=1,3
          n(1,j) = (Upos(j)-Vpos(j))/d
          nT(j,1) = n(1,j)
      enddo

      call matmlt1x3(U13, n, U)
      call matmlt1x3(V13, n, V)
      Un = dot(U13, nT)
      Vn = dot(V13, nT)

      r = abs(Un - Vn)

      return
      end

C     ******************************************************************
      SUBROUTINE project_U(UMAT,VEC,UPROJ)
C     ******************************************************************
C --- Calculate projection of UMAT parallel to vector VEC
C --- Assume orthogonal coordinates

      REAL UMAT(6),VEC(3),UPROJ
      REAL V1(3)

      VEC2 = VEC(1)**2 + VEC(2)**2 + VEC(3)**2

      V1(1) = UMAT(1)*VEC(1) + UMAT(4)*VEC(2) + UMAT(5)*VEC(3)
      V1(2) = UMAT(4)*VEC(1) + UMAT(2)*VEC(2) + UMAT(6)*VEC(3)
      V1(3) = UMAT(5)*VEC(1) + UMAT(6)*VEC(2) + UMAT(3)*VEC(3)

      UPROJ = (VEC(1)*V1(1) + VEC(2)*V1(2) + VEC(3)*V1(3))/VEC2

      END


********************************************************************************
*                 MATRIX MANIPULATION SUBROUTINES AND FUNCTIONS                *
********************************************************************************

C     ******************************************************************
      SUBROUTINE array(A,N,M)
C     ******************************************************************
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
      INTEGER N,M,I,J
      REAL*8 A(N,M)
Cf2py intent(in,out,copy) a
Cf2py integer intent(hide),depend(a) :: n=shape(a,0), m=shape(a,1)
      DO J=1,M
         A(1,J) = A(1,J)
      ENDDO
      DO I=1,N
         A(I,1) = A(I,1)
      ENDDO
      END

C     ******************************************************************
      SUBROUTINE array3x1(A,N,M)
C     ******************************************************************
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
      INTEGER N,M,I,J
      REAL*8 A(N,M)
C***Cf2py intent(in,out,copy) a
C***Cf2py integer intent(hide),depend(a) :: n=shape(a,0), m=shape(a,1)
      DO J=1,M
         A(1,J) = A(1,J)
      ENDDO
      DO I=1,N
         A(I,1) = A(I,1)
      ENDDO
      END

C     ******************************************************************
      SUBROUTINE mul3x3(C,A,B)
C     multiply a 3x3 matrix/array
C     ******************************************************************
C --- A SUBROUTINE FOR MULTIPLYING TWO 3*3 MATRICES TOGETHER.
      real A(3,3), B(3,3), C(3,3)
      integer i,j
C***Cf2py intent(in) a, b
C***Cf2py intent(out) c
C***Cf2py depend(a) c
      do i = 1,3
      do j = 1,3
          C(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j)
      enddo
      enddo
      return
      end

C     ******************************************************************
      SUBROUTINE matmlt3x3(C,A,B)
C     multiply a 3x3 matrix/array
C     ******************************************************************
C --- A SUBROUTINE FOR MULTIPLYING TWO 3*3 MATRICES TOGETHER.
      real A(3,3), B(3,3), C(3,3)
      integer i,j
C***Cf2py intent(in) a, b
C***Cf2py intent(out) c
C***Cf2py depend(a) c
      do i = 1,3
      do j = 1,3
          C(i,j) = A(i,1)*B(j,1) + A(i,2)*B(j,2) + A(i,3)*B(j,3)
      enddo
      enddo
      return
      end

C     ******************************************************************
      SUBROUTINE matmlt1x3(C,A,B)
C     ******************************************************************
C --- A SUBROUTINE FOR MULTIPLYING (1x3)*(3x3) MATRICES TOGETHER.
      real A(1,3), B(3,3), C(1,3)
      integer i,j

      i = 1
      do j = 1,3
          C(i,j) = A(i,1)*B(j,1) + A(i,2)*B(j,2) + A(i,3)*B(j,3)
      enddo
      return
      end

C     ******************************************************************
      SUBROUTINE trn3x3(B,A)
C     ******************************************************************
C --- transpose a 3x3 matrix/array
      real    A(3,3), B(3,3)
      integer i, j
C***Cf2py intent(in) a
C***Cf2py intent(out) b
C***Cf2py depend(a) b
      do i = 1, 3
      do j = 1, 3
          B(i,j) = A(j,i)
      enddo
      enddo
      return
      end

C     ******************************************************************
      subroutine add3x3(C,A,B)
C     ******************************************************************
C --- add a 3x3 matrix/array
      real A(3,3), B(3,3), C(3,3)
      integer i,j
C***Cf2py intent(in) a, b
C***Cf2py intent(out) c
C***Cf2py depend(a) c
      do i = 1,3
      do j = 1,3
          C(i,j) = A(i,j) + B(i,j)
      enddo
      enddo
      return
      end

C     ******************************************************************
      subroutine sum3x3(B,A)
C     ******************************************************************
C --- sum a 3x3 matrix/array
      real A(3,3), B
      integer i,j
C***Cf2py intent(in) a
C***Cf2py intent(out) b
C***Cf2py depend(a) b
      B = 0.0
      do i = 1,3
      do j = 1,3
          B = B + A(i,j)
      enddo
      enddo
      return
      end

C     ******************************************************************
      subroutine subtract3x3(C,A,B)
C     ******************************************************************
C --- subtract a 3x3 matrix/array
      real A(3,3), B(3,3), C(3,3)
      integer i,j
C***Cf2py intent(in) a, b
C***Cf2py intent(out) c
C***Cf2py depend(a) c
      do i = 1,3
      do j = 1,3
          C(i,j) = A(i,j) - B(i,j)
      enddo
      enddo
      return
      end

C     ******************************************************************
      subroutine trace3x3(T,A)
C     ******************************************************************
C --- trace of a 3x3 matrix/array
      real A(3,3), T
      integer i
C***Cf2py intent(in) a
C***Cf2py intent(out) t
C***Cf2py depend(a) t
      T = 0.0
      do i = 1,3
          T = T + A(i,i)
      enddo
      return
      end

C     ******************************************************************
      SUBROUTINE CROSS(B,C,A)
C     ******************************************************************
C --- cross product of two 3x3 matrices/arrays
      REAL A(3),B(3),C(3)
      A(1)=B(2)*C(3)-C(2)*B(3)
      A(2)=B(3)*C(1)-C(3)*B(1)
      A(3)=B(1)*C(2)-C(1)*B(2)
      RETURN
      END

C     ******************************************************************
      FUNCTION DOT(A,B)
C     ******************************************************************
C --- dot product of two 3x3 matrices/arrays
      REAL A(3),B(3)
      DOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      RETURN
      END

C     ******************************************************************
      FUNCTION inv3x3( A, B )
C     ******************************************************************
C --- inverse of a 3x3 matrix/array
      REAL A(3,3),B(3,3),C(3,3),D
      CALL CROSS(B(1,2),B(1,3),C(1,1))
      CALL CROSS(B(1,3),B(1,1),C(1,2))
      CALL CROSS(B(1,1),B(1,2),C(1,3))
      D=DOT(B(1,1), C(1,1))
      IF (ABS(D) .LE. 1.E-30) THEN
          inv3x3 = 0.0
          RETURN
      ENDIF
      DO I=1,3
      DO J=1,3
          A(I,J) = C(J,I)/D
      ENDDO
      ENDDO
      inv3x3 = D
      RETURN
      END

C     ******************************************************************
      subroutine assert (logic, message)
C     ******************************************************************
      logical logic
      character*(*) message
      if (logic) return
      write(0,*) '>>> ',message,' <<<'
      stop
      end


************************************************************************
*              Matrix inversion via LU decomposition                   *
*       adapted from Numerical Recipes in Fortran (1986)               *
************************************************************************
*
CCC     input  NxN matrix A is replaced by its LU decomposition
CC      output index(N) records row permutation due to pivoting
C       output D is parity of row permutations
c
C     ******************************************************************
        subroutine ludcmp( A, n, index, D )
C     ******************************************************************
        implicit  NONE
        integer   n,index(n)
        real      A(n,n)
        real      D
 
        integer   i,imax,j,k
        real      aamax,dum,sum
        integer    NMAX
        parameter (NMAX=10)
        real      vv(NMAX)
 
        d = 1.
        do i=1,n
            aamax = 0.
            do j=1,n
                if (abs(A(i,j)).gt.aamax) aamax = abs(A(i,j))
            enddo
            call assert(aamax.ne.0.,'Singular matrix')
            vv(i) = 1. / aamax
        enddo
        do j=1,n
            if (j.gt.1) then
                do i=1,j-1
                    sum = A(i,j)
                    if (i.gt.1) then
                        do k=1,i-1
                            sum = sum - A(i,k)*A(k,j)
                        enddo
                        A(i,j) = sum
                    endif
                enddo
            endif
            aamax = 0.
            do i=j,n
                sum = A(i,j)
                if (j.gt.1) then
                    do k=1,j-1
                        sum = sum - A(i,k)*A(k,j)
                    enddo
                    A(i,j) = sum
                endif
                dum = vv(i) * abs(sum)
                if (dum.ge.aamax) then
                    imax = i
                    aamax = dum
                endif
            enddo
            if (j.ne.imax) then
                do k=1,n
                    dum = A(imax,k)
                    A(imax,k) = A(j,k)
                    A(j,k) = dum
                enddo
                d = -d
                vv(imax) = vv(j)
            endif
            index(j) = imax
            call assert(A(j,j).ne.0.,'Singular matrix')
            if (j.ne.n) then
                dum = 1. / A(j,j)
                do i=j+1,n
                    A(i,j) = A(i,j) * dum
                enddo
            endif
        enddo
        return
        end

C     ******************************************************************
        subroutine lubksb( A, N, index, B )
C     ******************************************************************
CCC     corresponding back-substitution routine
CC
C
        implicit NONE
        integer  n, index(n)
        real     A(n,n), B(n)
 
        integer  i,ii,j,ll
        real     sum
 
        ii = 0
        do i=1,n
            ll = index(i)
            sum = B(ll)
            B(ll) = B(i)
            if (ii.ne.0) then
                do j=ii,i-1
                    sum = sum - A(i,j)*B(j)
                enddo
            else if (sum.ne.0.) then
                ii = i
            endif
            B(i) = sum
        enddo
 
        do i=n,1,-1
            sum = B(i)
            if (i.lt.n) then
                do j=i+1,n
                    sum = sum - A(i,j)*B(j)
                enddo
            endif
            B(i) = sum / A(i,i)
        enddo
        return
        end


************************************************************************
* Find eigenvalues and eigenvectors of nxn symmetric matrix using      *
* cyclic Jacobi method. NP is (Fortran) physical storage dimension.    *
* On return A is destroyed, D contains eigenvalues,and each column of  *
* V is a normalized eigenvector                                        *
* Adapted from Numerical Recipes in Fortran (1986)                     *
* We only need it for 3x3 symmetric matrices, so it's overkill.        *
************************************************************************
CCC
CC
      subroutine eigenvalues(eigens,U)
C --- sqd = abs(numpy.sum(numpy.subtract(U,V)**2))
      real U(3,3)
      real EIGENS(3), EVECS(4,4)
Cf2py intent(in) u
Cf2py intent(out) eigens
Cf2py depend(u) eigens

      call jacobi(U,3,4,EIGENS,EVECS)

c     Units of this matrix are A**2; we want values in A
      do i = 1,3
          if (EIGENS(i).gt.0.) then
              EIGENS(i) = sqrt(EIGENS(i))
          else
C             write (noise,*) 'Non-positive definite ellipsoid!'
              EIGENS(i) = 0.
          endif
      enddo

      return
      end
C
        subroutine jacobi( A, n, np, D, V )
C                  jacobi(UU, 3, 4, EIGENS, EVECS )
        PARAMETER (NMAX=4)
c       Machine dependent! (converge when off-diagonal sum is less than this)
        PARAMETER (TINY = 1.e-37)
        real A(np,np), D(n), V(np,np)
        real B(NMAX), Z(NMAX)
C *** Cf2py intent(in) a, n, np
C *** Cf2py intent(out) eigens, evecs
C *** Cf2py depend(a) eigens, evecs
c
        call assert(n.le.NMAX,'Matrix too big for eigenvector routine')

c
c       Initialize V to identity matrix, B and D to diagonal of A
        do ip = 1,n
        do iq = 1,n
            V(ip,iq) = 0.0
        enddo
        V(ip,ip) = 1.0
        B(ip) = A(ip,ip)
        D(ip) = B(ip)
        Z(ip) = 0.0
        enddo
        nrot = 0

c
c       50 sweeps is never expected to happen; Press et al (1986) claim
c       6 to 10 are typical for moderate matrix sizes
c       Empirical trials of rastep show 6 sweeps, 8-10 rotations
        do i = 1,50
            sm = 0.0
            do ip = 1,n-1
            do iq = ip+1,n
                sm = sm+abs(A(ip,iq))
            enddo
            enddo
            if (sm .lt. TINY) return
c           After 4 sweeps skip rotation if the off-diagonal is small
            if (i .lt. 4) then
                thresh = 0.2*sm/(n*n)
            else
                thresh = 0.0
            endif
            do ip = 1,n-1
            do iq = ip+1,n
                g = 100. * abs(A(ip,iq))
                if ((i.gt.4) .and.
     &              (abs(D(ip)) + g .eq. abs(D(ip))) .and.
     &              (abs(D(iq)) + g .eq. abs(D(iq)))) then
                    A(ip,iq) = 0.0
                else if (abs(A(ip,iq)).gt.thresh) then
                    h = D(iq) - D(ip)
                    if (abs(h) + g .eq. abs(h)) then
                        t = A(ip,iq) / h
                    else
                        theta = 0.5 * h / A(ip,iq)
                        t = 1. / (abs(theta) + sqrt(1.+theta*theta))
                        if (theta.lt.0.) t = -t
                    endif
                    c = 1./sqrt(1.+t*t)
                    s = t * c
                    tau = s/(1.+c)
                    h = t * A(ip,iq)
                    Z(ip) = Z(ip) - h
                    Z(iq) = Z(iq) + h
                    D(ip) = D(ip) - h
                    D(iq) = D(iq) + h
                    A(ip,iq) = 0.0
                    do j = 1,ip-1
                        g = A(j,ip)
                        h = A(j,iq)
                        A(j,ip) = g - s*(h+g*tau)
                        A(j,iq) = h + s*(g-h*tau)
                    enddo
                    do j = ip+1,iq-1
                        g = A(ip,j)
                        h = A(j, iq)
                        A(ip,j) = g - s*(h+g*tau)
                        A(j,iq) = h + s*(g-h*tau)
                    enddo
                    do j = iq+1,n
                        g = A(ip,j)
                        h = A(iq,j)
                        A(ip,j) = g - s*(h+g*tau)
                        A(iq,j) = h + s*(g-h*tau)
                    enddo
                    do j = 1,n
                        g = V(j,ip)
                        h = V(j,iq)
                        V(j,ip) = g - s*(h+g*tau)
                        V(j,iq) = h + s*(g-h*tau)
                    enddo
                    nrot = nrot + 1
                endif
            enddo
            enddo
c
c           Update D with sum and reinitialize Z
            do ip = 1,n
                B(ip) = B(ip) + Z(ip)
                D(ip) = B(ip)
                Z(ip) = 0.0
            enddo
        enddo
        call assert(.false.,'Failed to find eigenvectors')
        return
        end
c

