      subroutine wt_scl(no,ni,x,lx,xmiss,d,tg,cst,
     +                      eps,kn,ds,wi,pwr,critrb,wtchg)
      implicit none


c
c       Declare variables
c

        integer no,ni,kn
        double precision x(no,ni),wi(no,ni),d(no*(no-1)/2),ds(ni)
        double precision tg(ni,2),sc(kn),w0(ni)
        integer lx(ni),mn(kn),m(kn)
        integer kq, k, j, kmx, lp, ld, i, nnm, nni, iq
        double precision qt, eps, s, sw, pwr
        double precision s1, s2, s3
        double precision cst(ni), xmiss, wtchg
        double precision big /9.9e35/
        double precision critrb, dm(ni)
c
c       Start code
c
        qt=0.5
        wtchg=0.0
        critrb=0.0
        kq=0
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

          !print *, k,mn
          !print *, k,sc

          sw=0.0
          s1 = 0.0
          s2 = 0.0
          s3 = 0.0
          do j=1,ni
            w0(j)=wi(k,j)
          end do

          do j=1,ni
            if ((x(k,j).lt.xmiss).and.(lx(j).ne.0)) then
              if (lx(j).eq.1) then
                nnm=0
                do i=1,kn
                  if (x(mn(i),j).lt.xmiss) then
                    nnm=nnm+1
                    sc(nnm)=abs(x(mn(i),j)-x(k,j))
                    if(sc(nnm).ne.0) then
                      sc(nnm)=exp(pwr*log(ds(j)*sc(nnm)))
                    end if
                    m(nnm)=nnm
                  end if
                end do
                if (nnm.eq.0) then
                  wi(k,j)=0.0
                  dm(j) = 0.0
                  goto 10681
                else
                  iq=int(qt*nnm+1)
                  call psort7(sc,m,1,nnm)
                  wi(k,j)=exp(-sc(m(iq))/eps)! underflow occurs when ds*sc/eps < 1.0.e-45 ?
                  dm(j)=sc(m(iq))
                  !print *, sc(m(iq))
                end if
              elseif(lx(j).eq.2) then
                nnm=0
                do i=1,kn
                  if (x(mn(i),j).lt.xmiss) then
                    nnm=nnm+1
                    sc(nnm)=max(abs(x(k,j)-tg(j,1)),
     +                          abs(x(mn(i),j)-tg(j,1)))
                    if(sc(nnm).ne.0) then
                      sc(nnm)=exp(pwr*log(ds(j)*sc(nnm)))
                    end if
                    m(nnm)=nnm
                  end if
                end do
                if (nnm.eq.0) then
                  wi(k,j)=0.0
                  dm(j) = 0.0
                  goto 10681
                else
                  iq=int(qt*nnm+1)
                  call psort7(sc,m,1,nnm)
                  wi(k,j)=exp(-sc(m(iq))/eps)! underflow occurs when ds*sc/eps < 1.0.e-45 ?
                  dm(j)=sc(m(iq))
                  !print *, sc(m(iq))
                end if
              elseif(lx(j).eq.3) then
                nnm=0
                do i=1,kn
                  if (x(mn(i),j).lt.xmiss) then
                    nnm=nnm+1
                    s1=max(abs(x(k,j)-tg(j,1)),abs(x(mn(i),j)-tg(j,1)))
                    s2=max(abs(x(k,j)-tg(j,2)),abs(x(mn(i),j)-tg(j,2)))
                    sc(nnm)=min(s1,s2)
                    if(sc(nnm).ne.0) then
                      sc(nnm)=exp(pwr*log(ds(j)*sc(nnm)))
                    end if
                    m(nnm)=nnm
                  end if
                end do
                if (nnm.eq.0) then
                  wi(k,j)=0.0
                  dm(j) = 0.0
                  goto 10681
                else
                  iq=int(qt*nnm+1)
                  call psort7(sc,m,1,nnm)
                  wi(k,j)=exp(-sc(m(iq))/eps)! underflow occurs when ds*sc/eps < 1.0.e-45 ?
                  dm(j)=sc(m(iq))
                  !print *, sc(m(iq))
                end if
              elseif(lx(j).eq.333) then
                nnm=0
                do i=1,kn
                  if (x(mn(i),j).lt.xmiss) then
                    nnm=nnm+1
                    s1=max(abs(x(k,j)-tg(j,1)),abs(x(mn(i),j)-tg(j,1)))
                    s2=max(abs(x(k,j)-tg(j,2)),abs(x(mn(i),j)-tg(j,2)))
                    s3=max(abs(x(k,j)),abs(x(mn(i),j))) != - 0.0, target = 0.0
                    sc(nnm)=min(s1,s2,s3,cst(j)*abs(x(mn(i),j)-x(k,j)))
                    if(sc(nnm).ne.0) then
                      sc(nnm)=exp(pwr*log(ds(j)*sc(nnm)))
                    end if
                    m(nnm)=nnm
                  end if
                end do
                if (nnm.eq.0) then
                  wi(k,j)=0.0
                  dm(j) = 0.0
                  goto 10681
                else
                  iq=int(qt*nnm+1)
                  call psort7(sc,m,1,nnm)
                  wi(k,j)=exp(-sc(m(iq))/eps)! underflow occurs when ds*sc/eps < 1.0.e-45 ?
                  dm(j)=sc(m(iq))
                  !print *, sc(m(iq))
                end if
              else ! lx(j).gt.4
                nnm=0
                nni=nnm
                do i=1,kn
                  if (x(mn(i),j).lt.xmiss) then
                    nnm=nnm+1
                    if (x(mn(i),j).ne.x(k,j)) nni=nni+1
                  end if
                end do
                if (nnm.eq.0) then
                  wi(k,j)=0.0
                  dm(j) = 0.0
                  goto 10681
                else
                  wi(k,j)= exp(-(ds(j)*nni)/(nnm*eps))
                  ! underflow occures when ds*sc/eps < 1.0.e-45
                  dm(j) = ds(j)*nni/nnm
                end if
              end if
              sw=sw+wi(k,j)
            else
              wi(k,j)= 0.0
              dm(j) = 0.0
            end if
10681       continue

          end do


          do j=1,ni
            wi(k,j)=wi(k,j)/sw
            wtchg=wtchg+abs(wi(k,j)-w0(j))

            if (wi(k,j).gt.0.0) then
               critrb = critrb + wi(k,j)*dm(j) +
     +                 eps*wi(k,j)*log(wi(k,j))
            end if

          end do
          critrb = critrb + eps*log(float(ni))

        end do ! for every observation k...
        !wtchg=wtchg/(2*no*(1.0-1.0/ni))
        wtchg=wtchg/(2*no*(1.0-1.0/ni))
        return
      end subroutine wt_scl
