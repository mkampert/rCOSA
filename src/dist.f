      subroutine dist(no,ni,x,lx,xmiss,ds,tg,cst,wi,
     +                  ipwr,pwr,eps,d)

      implicit none

c
c       Declare variables
c

      integer no,ni,ipwr
      double precision x(no,ni),ds(ni),wi(no,ni),d(no*(no-1)/2)
      double precision tg(ni,2)
      integer lx(ni), kp, k
      integer nv, j, l
      double precision s,s2,sw,sw2
      double precision djkl, d0, d1, d2, d12, d3, u
      double precision eps, t, pwr, xmiss
      double precision big /9.9e35/

c
c        Variables for Soft Targeting
c

      double precision ds0, ds12
      double precision extg(2,ni), sqtg(2,ni), cst(ni) ! "c"orrection on "s"pread due to "t"argeting
      integer nps(ni)

c
c       Avoid warnings regarding initialization
c
      
      l = 0
      djkl = big


c
c       Start code
c

      kp=0
      do j=1,ni
         nps(j)=0
         extg(1,j)=0
         extg(2,j)=0
         sqtg(1,j)=0
         sqtg(2,j)=0
      end do

      do k=1,no-1
        do l=k+1,no
          kp=kp+1
          s=0.0
          s2=s
          sw=s2
          sw2=sw
          nv=0
          do j=1,ni
            if ((lx(j).ne.0).and.
     +          (x(k,j).lt.xmiss).and.
     +          (x(l,j).lt.xmiss)) then
              if (lx(j).eq.1) then
                djkl=abs(x(k,j)-x(l,j))
c                if(j.eq.1) then
c                  call printint(kp)
c                  call traceprint(x(k,j), x(l,j),
c     +                      ds(j), x(k,j)-x(l,j))
c                end if
                if (djkl.ne.0) then
                  djkl= exp(pwr*log(ds(j)*djkl))
                end if
c                if(j.eq.1) then
c                  call printdouble(djkl)
c                end if
              elseif (lx(j).eq.2) then
                djkl=max(abs(x(k,j)-tg(j,1)),abs(x(l,j)-tg(j,1)))
                if (djkl.ne.0) then
                  djkl= exp(pwr*log(ds(j)*djkl))
                end if
              elseif (lx(j).eq.3) then
                d1=max(abs(x(k,j)-tg(j,1)),abs(x(l,j)-tg(j,1)))
                d2=max(abs(x(k,j)-tg(j,2)),abs(x(l,j)-tg(j,2)))
                !djkl=min(d1,d2,abs(x(k,j)-x(l,j)))
                djkl=min(d1,d2)
                if (djkl.ne.0) then
                  djkl= exp(pwr*log(ds(j)*djkl))
                end if
              elseif (lx(j).eq.333) then

                d0=cst(j)*abs(x(k,j)-x(l,j))
                d1=max(abs(x(k,j)-tg(j,1)),abs(x(l,j)-tg(j,1)))
                d2=max(abs(x(k,j)-tg(j,2)),abs(x(l,j)-tg(j,2)))
                d3=max(abs(x(k,j)),abs(x(l,j))) ! - 0.0, 3rd target = 0.0
                d12=min(d1,d2,d3)

                extg(1,j)=extg(1,j)+d0 ! for the expectation
                extg(2,j)=extg(2,j)+d12 ! for the expectation
                sqtg(1,j)=sqtg(1,j)+d0*d0 ! for the spread
                sqtg(2,j)=sqtg(2,j)+d12*d12 ! for the spread

c                if(j.eq.1457) print *, sqtg(1,j), sqtg(2,j), kp
                djkl=min(d12,d0)
                nps(j)=nps(j) + 1

                if (djkl.ne.0) then
                  djkl= exp(pwr*log(ds(j)*djkl))
                end if
            elseif (lx(j).eq.4) then
                if (x(k,j).eq.x(l,j)) then
                  djkl=0.0
                else
                  djkl=ds(j)
                end if
              elseif (lx(j).eq.5) then
                if (x(k,j).ne.tg(j,1).or.x(l,j).ne.tg(j,1)) then
                  djkl=ds(j)
                else
                  djkl=0.0
                end if
              elseif (lx(j).eq.6) then
                if ((x(k,j).ne.tg(j,1).or.x(l,j).ne.tg(j,1)).and.
     +              (x(k,j).ne.tg(j,2).or.x(l,j).ne.tg(j,2))) then
                  djkl=ds(j)
                else
                  djkl=0.0
                end if
              end if
              u=exp(djkl/eps)
              if (ipwr.eq.1) then
                t=max(wi(k,j),wi(l,j))
                s=s+t/u
                sw=sw+t
              else
                s=s+wi(k,j)/u
                sw=sw+wi(k,j)
                s2=s2+wi(l,j)/u
                sw2=sw2+wi(l,j)
              end if
              nv=nv+1
            end if
          end do
          if(sw .gt. 0.0)goto 10521
          d(kp)=big
          goto 10511
10521     if(ipwr .ne. 2)goto 10531
          s=min(s/sw,s2/sw2)
          d(kp)=-eps*log(s)
          goto 10541
10531     continue
          d(kp)=-ni*eps*sw*log(s/sw)/nv
10541     continue
10511     continue
        end do
      end do


c     Update targeting spread corrections for the next iteration

      ds12=1.0
      ds0=1.0
      do j=1,ni
        if ((lx(j).eq.333).and.
     +          (x(k,j).lt.xmiss).and.
     +          (x(l,j).lt.xmiss)) then
          ds12=sqtg(2,j) -
     +      (extg(2,j)/(nps(j)))*(extg(2,j)/(nps(j)))

          ds0=sqtg(1,j) -
     +      (extg(1,j)/(nps(j)))*(extg(1,j)/(nps(j)))
          cst(j)=sqrt(ds12/ds0)
c          if(j.eq.1457) then
c            print *, cst(j), j, sqtg(1,j),
c     +        (extg(1,j)/(nps(j)))*(extg(1,j)/(nps(j))),
c     +        sqtg(2,j),
c     +        (extg(2,j)/(nps(j)))*(extg(2,j)/(nps(j)))
c          end if
        end if
      end do

      return
      end subroutine dist
