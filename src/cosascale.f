      subroutine cosascale(no,ni,x,lx,ip,dp,nv,ds,tg)
       implicit none

c
c       Declare variables
c
        integer ni,no,nv
        integer lx(ni),stan
        double precision x(no,ni),ds(ni),tg(ni,2)

        integer ip(6)

        double precision dp(10)
        double precision xmiss, flo, fhi

        logical robust
        integer m(no),ktrg,ks,nc,nm
        integer i,j,i1,i3,j1,j2
        double precision psq,r,sum,var,xk
c
c       Start code
c
        ktrg=ip(5)
        stan=ip(3)
        robust=(stan.eq.1)

        xmiss=dp(1)
        flo=dp(3)
        fhi=dp(4)

        ks=0
        nv=0
        do j=1,ni
c          print *, j
          if (lx(j).eq.0) goto 10101
          sum=0.0
          nv=nv+1
          do i=1,no
            m(i)=i
            sum=sum+x(i,j)
          end do
cccccc    Standardize data column-wise
          if (stan.ne.0) then
cccccc      if (.not.robust) then
            sum=sum/float(no)
            var=0.0
            do i=1,no
              x(i,j)=x(i,j)-sum
              var=var+x(i,j)**2
            end do
            var=sqrt(var/(float(no)-1))
            if (.not.robust) then
              do i=1,no
                x(i,j)=x(i,j)/var
              end do
            end if
          end if
          call psort7(x(1,j),m,1,no)
          if (x(m(1),j).ge.x(m(no),j)) then
            ks=1
c         Needs fixing to work in R!!
c            write(*,
c     +        '('' All values of attribute'',i6,'' are the same.'')') j
            goto 10101
          end if
          nm=no
          do while (nm.gt.1)
            if (x(m(nm),j).lt.xmiss) exit
            nm=nm-1
          end do
          if (nm.le.1) then
            ks=1
c         Needs fixing to work for R!!
c            write(*,
c     +        '('' Only one nonmissing value on attribute'',i6)') j
c            goto 10101
          end if
          if (x(m(1),j).ge.x(m(nm),j)) then
            ks=1
c         Needs fixing to Work in R!!
c           write(*,
c     +        '('' All nonmissing values of attribute'',
c     +        i6,'' are the same.'')') j
            goto 10101
          end if
          if (lx(j).le.3) then
            i1=int(0.25*nm)
            i3=int(0.75*nm)
            r=x(m(i3),j)-x(m(i1),j)
            do while (r.le.0.0)
              i3=min(i3+1,nm)
              i1=max(i1-1,1)
              r=x(m(i3),j)-x(m(i1),j)
            end do
            if (ktrg.ne.0) then
              j1=max(int(flo*nm+0.5),1)
              j2=min(int(fhi*nm+0.5),nm)
              tg(j,1)=x(m(j2),j)
              tg(j,2)=x(m(j1),j)
            end if
            ds(j)=1.35/r
            if (.not.robust) ds(j)=1
            if (stan.eq.0) ds(j)=1
cc          if(j. eq. 1) write(6,601) ds(j),r,nm
cc  601       format('ds(j),r ',2F8.3,i5)
            goto 10101
          end if
          xk=x(m(1),j)
          nc=1
          psq=0.0
          do i=2,nm
            if (x(m(i),j).le.xk) then
              nc=nc+1
            else
              psq=psq+(float(nc)/nm)**2
              nc=1
              xk=x(m(i),j)
            end if
          end do
          psq=psq+(float(nc)/nm)**2
          ds(j)=1.0/(1.0-psq)
10101     continue
        end do
        if (robust) go to 9
cccccc  Standardize data row-wise
        do i=1,no
          sum=0.0
          do j=1,ni
            sum=sum+x(i,j)
          end do
          sum=sum/float(ni)
          var=0.0
          do j=1,ni
            x(i,j)=x(i,j)-sum
            var=var+x(i,j)**2
          end do
          var=sqrt(var/(float(ni)))
          do j=1,ni
            x(i,j)=x(i,j)/var
          end do
        end do
cccccc  Standardize the data again columnwise
        do j=1,ni
          sum=0.0
          do i=1,no
            sum=sum+x(i,j)
          end do
          sum=sum/float(no)
          var=0.0
          do i=1,no
            x(i,j)=x(i,j)-sum
            var=var+x(i,j)**2
          end do
          var=sqrt(var/(float(no)-1))
          do i=1,no
            x(i,j)=x(i,j)/var
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccc
c    9   print *,ds
9        continue


c
c       Commented out because it should have been checked in R already!
c
c        if (ks.ne.0) stop
c        if (nv.eq.0) then
cc       Needs Fixing to work in R!
cc          write(*,'('' All attributes deleted (lx(j)=0)'')')
c          stop
c        end if


        return
      end subroutine cosascale
