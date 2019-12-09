      subroutine COSAcrit(no,ni,x,lx,d,kn,dp,ds,wi,tg,crit)
      implicit none

c NOTE: only for the L1 distances!!

c
c       Declare variables
c

        integer no,ni,kn
        integer lx(ni),mn(kn)
        double precision x(no,ni),wi(no,ni),d(no*(no-1)/2)
        double precision ds(ni),sc(kn),tg(ni,2), dp(10)

        double precision crit,d1,d2, sigma, xmiss
        integer i,j,k, kq

        double precision big /9.9e35/
        integer kmx, lp,ld
        double precision s


c
c       Start the code

        xmiss = dp(1)
        sigma = dp(2)
        crit = 0.0
        kq = 0

        !write(6, '(2g14.4)') eps, xmiss ! sigma appears to be 0?!?



        do k=1,no

          do j=1,kn
            sc(j)=big
            mn(j)=0
          end do

          kmx=1
          s=big

          if (k.gt.1) then
            lp=k-1
            ld=1
            do i=1,k-1
              if (d(lp).lt.s) then
                sc(kmx)=d(lp)
                mn(kmx)=i
                s=0.0
                do j=1,kn
                  if (sc(j).gt.s) then
                    s=sc(j)
                    kmx=j
                  end if
                end do
              end if
              ld=ld+1
              lp=lp+no-ld
            end do
          end if

          ! kq = 0  c this should not be here?!
          if (k.lt.no) then
            do i=k+1,no
              kq=kq+1
              if (d(kq).lt.s) then
                sc(kmx)=d(kq)
                mn(kmx)=i
                s=0.0
                do j=1,kn
                  if (sc(j).gt.s) then
                    s=sc(j)
                    kmx=j
                  end if
                end do
              end if
            end do
          end if

          !write(666, '("crit = ", 1g10.5)' ) crit
          !write(666, '("new k = ", 1i3)' ) k
          do i=1,kn
            !write(666, '("index = ", 1i3)' ) mn(i)
            do j=0,ni
              !write(666, '("crit = ", 1g10.5)' ) crit
              !write(666, '("index = ", 1i3)' ) mn(i)
              if (
     +             (lx(j).ne.0).and.
     +             (x(k,j).lt.xmiss).and.
     +             (x(mn(i),j).lt.xmiss).and.
c     line 96 needed because of underflow problems... still needs to be solved with double precision??
     +             (wi(k,j).gt.1.0d-45)
c     line 96 needed because of underflow problems... still needs to be solved with double precision??
     +            ) then
                if (lx(j).eq.1) then
                  crit= crit + wi(k,j)*ds(j)*abs(x(k,j)-x(mn(i),j))
                  crit = crit + sigma*wi(k,j)*log(wi(k,j))
                elseif (lx(j).eq.2) then
                  d1=max(abs(x(k,j)-tg(j,1)),abs(x(mn(i),j)-tg(j,1)))
                  crit = crit + wi(k,j)*ds(j)*d1
                  crit = crit + sigma*wi(k,j)*log(wi(k,j))
                elseif (lx(j).eq.3) then
                  d1=max(abs(x(k,j)-tg(j,1)),abs(x(mn(i),j)-tg(j,1)))
                  d2=max(abs(x(k,j)-tg(j,2)),abs(x(mn(i),j)-tg(j,2)))
                  crit = crit + wi(k,j)*ds(j)*min(d1,d2)
                  crit = crit + sigma*wi(k,j)*log(wi(k,j))
                elseif (lx(j).eq.4) then
                  if (x(k,j).eq.x(mn(i),j)) then
                    crit = crit + sigma*wi(k,j)*log(wi(k,j))
                  else
                    crit = crit + wi(k,j)*ds(j)
                    crit = crit + sigma*wi(k,j)*log(wi(k,j))
                  end if
                elseif (lx(j).eq.5) then
                  if (x(k,j).ne.tg(j,1).or.x(mn(i),j).ne.tg(j,1)) then
                    crit = crit + wi(k,j)*ds(j)
                    crit = crit + sigma*wi(k,j)*log(wi(k,j))
                  else
                    crit = crit + sigma*wi(k,j)*log(wi(k,j))
                  end if
                elseif (lx(j).eq.6) then
                   if (x(k,j).ne.tg(j,1).or.x(mn(i),j).ne.tg(j,1)) then
                   if (x(k,j).ne.tg(j,2).or.x(mn(i),j).ne.tg(j,2)) then
c                  if (((x(k,j).ne.tg(j,1)).or.(x(mn(i),j).ne.tg(j,1))).and.
c                  +              (((x(k,j).ne.tg(j,2)).or.(x(mn(i),j).ne.tg(j,2))))) then
                       crit = crit + wi(k,j)*ds(j)
                       crit = crit + sigma*wi(k,j)*log(wi(k,j))
                   else
                       crit = crit + sigma*wi(k,j)*log(wi(k,j))
                   end if
                   end if
                end if


              end if
            end do
            crit = crit + sigma*log(float(ni))
            !write(6, '("crit = ", 1g10.5)' ) crit
          end do ! do i, kn

        end do ! do k, no
        !crit = crit/kn + sigma*no*log(float(no))
        crit = crit/kn

      return
      end subroutine COSAcrit

