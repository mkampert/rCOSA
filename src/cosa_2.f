c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))
c     adapted by J.J. Meulman 05-15-2009
c      subroutine cosa(no,ni,x,lx,tg,kn,ip,dp,d,wi) ! when using r_cosa_2.f
      subroutine cosa(no,ni,kn,x,ip,dp,lx,tg,d,wi)
       
        implicit none

c
c       Declare and define variables
c
        integer no, ni, kn, nv
c       ipams <- c(niter, noit, stand, wtcom, ktarg, ltarg, itot)
        integer ip(7)
        integer nit, noit, stan, ipwr, ktarg!, ltarg

c       dpams <- c(xmisst, lambda, qntls, relax, conv, pwr)
        double precision dp(10)
        double precision xmiss, sigma, flo, fhi, rate
        double precision dwt, pwr, wtchg, wtot

c        x, tg, sp, ms, lx, d, dl1, wi,
        double precision x(no, ni), tg(ni, 2)
        integer lx(ni)
        double precision sc(kn+ni+1)
        double precision w0,wtchgold,eps
        double precision d(no*(no-1)/2)
        double precision wi(no, ni), cst(ni)
        double precision crit, critrb

        integer j,i,itot,it,oit,ii,icount


c       ipams <- c(niter, noit, stand, wtcom, ktarg, ltarg)
      nit=ip(1)
      noit=ip(2)
      stan=ip(3)
      ipwr=ip(4)
      ktarg=ip(5)
c     ltarg=ip(6) ! not needed?

c     dpams <- c(xmisst, lambda, qntls, relax, conv, pwr)
      xmiss=dp(1)
      sigma=dp(2)
      flo=dp(3)
      fhi=dp(4)
c        iprt=1 ! not used at all?
      rate=dp(5)
      dwt=dp(6)
      pwr=dp(7)

      wtchg = 0.0


c
c       Start code
c

      call runningcosa()


      call cosaheadlog()

      !scale the data
      call cosascale(no,ni,x,lx,ip,dp,nv,sc,tg)

      !initialize weights
      w0=1.0/ni
      do j=1,ni
        cst(j)=1.0
        do i=1,no
          wi(i,j)=w0
        end do
      end do

      ! call dissimilarities
c      call dist(no,ni,x,lx,xmiss,sc,tg,wi,ipwr,pwr,sigma,d)
      call dist(no,ni,x,lx,xmiss,sc,tg,cst,wi,ipwr,pwr,sigma,d)
      ! print *, d

      if (nit.le.0.or.nv.eq.1) return
      wtchgold=dwt
      wtot=0.0
      itot=0

      ! COSA iterations:
      do oit=1, noit
        do it=1, nit
          itot=itot+1

          ! calculate weights

          call wt_scl(no,ni,x,lx,xmiss,d,tg,cst,
     +            sigma,kn,sc,wi,
     +            pwr,critrb, wtchg)

          !call checkpoint()

c         Calculate criterion
          call COSAcrit(no,ni,x,lx,d,kn,dp,sc,wi,tg,crit)


          !call checkpoint()

          eps=sigma*(rate*oit+1.0)

          ! calculate dissimilarities
          call dist(no,ni,x,lx,xmiss,sc,tg,cst,wi,ipwr,pwr,eps,d)

          !call checkpoint()

          ii=11
          icount=MOD(itot,ii)
          wtchgold=wtchgold+wtchg

          call oldcosatracelog(wtchg, it, oit, itot, eps, critrb, crit)

          if(icount .eq. 0) wtchgold=wtchg
          if (wtchg.le.dwt) goto 162

          call rchkusr() ! allow user interrupt

        end do

  162   continue

        wtot=wtchgold/FLOAT(icount + 1)
        if ((wtot.le.dwt.and.icount.gt.9)) goto 163

      end do
  163   continue

      dp(8) = crit
      dp(9) = critrb
      dp(10) = eps
      ip(7) = itot

      return
      end subroutine cosa
