      subroutine psort7(v,a,ii,jj)
      implicit none
!
! puts into a the permutation vector which sorts v into
! increasing order. the array v is not modified.
! only elements from ii to jj are considered.
! arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
!
! this is a modification of cacm algorithm #347 by r. c. singleton,
! which is a modified hoare quicksort.
!
      dimension a(jj),v(jj),iu(20),il(20)
      integer t,tt
      integer a,m,i,j,ii,jj,k,l,ij,il,iu
      double precision v,vt,vtt
      m=1
      i=ii
      j=jj
   10 if (i.ge.j) go to 80
   20 k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
   30 l=j
      if (v(a(j)).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
      go to 50
   40 a(l)=a(k)
      a(k)=tt
   50 l=l-1
      if (v(a(l)).gt.vt) go to 50
      tt=a(l)
      vtt=v(tt)
   60 k=k+1
      if (v(a(k)).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
   70 il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
   80 m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
   90 if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
  100 i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(t)
      if (v(a(i)).le.vt) go to 100
      k=i
  110 a(k+1)=a(k)
      k=k-1
      if (vt.lt.v(a(k))) go to 110
      a(k+1)=t
      go to 100
      end subroutine psort7
