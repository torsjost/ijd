      Program hgw4

* Two dimensional hexagonal plane!! 2 bands!!
* Find Sigma GW and new G !!
* Write to file 117 !!
* Work in tau-space !!


*  beta    = 1/kT = 157681.2/T [1/Ry] and T in Kelvin.
*  beta    = 1/kT =  11594.2/T [1/eV] and T in Kelvin.
*
* T = 0.025 eV gives beta=1/T=40. OK because T = 0.025 corresponds
* to room temperature 300 K (beta=11594.2/300 is also around 40).
* Integrate over the full BZ.

* Run the DOS calculation first with many k-points!!

      implicit real*8 (a-h,o-z)
      parameter (maxp=1500,kdim=2)
      parameter (ngrp=6)
      parameter (nq0=1700)
      parameter (nwm= 550)
      parameter (nwT= 550)
* Old !
*      parameter (ntaum= 2*52+1)
      parameter (ntaum= 2*80+1)   !OK for nsimp=2 and 4
*      parameter (ntaum= 106)

*      parameter (nsimp= 2)
      parameter (nsimp= 4)

      double precision eA,eB,t1,t2,phi
      double precision rtemp1,rtemp2,rtemp3 

      double precision nos(maxp),dos(maxp)
      double precision tnos(maxp),tdos(maxp),x1(2),f1(2)
      double precision emin,emax,ebot,etop,e(2,8),ee(4),bb(4)
      double precision freq(maxp),rdos(maxp),cdos(maxp),
     .                 v2(maxp),v2d(maxp),rgamma(maxp)

* For G(tau).
* Odd Matsubara !!
      double precision rge(0:nwT-1),cge(0:nwT-1),
     .                 rgo(0:nwT-1),cgo(0:nwT-1)
      double precision rge11(0:nwT-1),cge11(0:nwT-1),
     .                 rgo11(0:nwT-1),cgo11(0:nwT-1)
      double precision rge12(0:nwT-1),cge12(0:nwT-1),
     .                 rgo12(0:nwT-1),cgo12(0:nwT-1)
      double precision rge21(0:nwT-1),cge21(0:nwT-1),
     .                 rgo21(0:nwT-1),cgo21(0:nwT-1)
      double precision rge22(0:nwT-1),cge22(0:nwT-1),
     .                 rgo22(0:nwT-1),cgo22(0:nwT-1)
      double precision rgv(-nwT:nwT-1),cgv(-nwT:nwT-1)

* For W !!
      double precision rwe(0:nwm),cwe(0:nwm),
     .                 wvi(0:nwm), wcoswt(0:nwm,ntaum)

      complex*16 wktau(2,2,nq0,ntaum)

      double precision tau(ntaum), rg0(ntaum), cg0(ntaum)
      double precision rgtau(ntaum), cgtau(ntaum)
      double precision rgtau11(ntaum), cgtau11(ntaum)
      double precision rgtau12(ntaum), cgtau12(ntaum)
      double precision rgtau21(ntaum), cgtau21(ntaum)
      double precision rgtau22(ntaum), cgtau22(ntaum)

      double precision rginf(ntaum), cginf(ntaum),
     .                 rgt(ntaum),cgt(ntaum),
     .                 rgt11(ntaum),cgt11(ntaum),
     .                 rgt12(ntaum),cgt12(ntaum),
     .                 rgt21(ntaum),cgt21(ntaum),
     .                 rgt22(ntaum),cgt22(ntaum)

* For tau to iw space !!
      double precision wcos(ntaum), wsin(ntaum)

* Odd Matsubara !!
      double precision coswt(0:nwm-1,ntaum)
      double precision sinwt(0:nwm-1,ntaum)
      double precision vi(0:nwm-1)

      double precision rge0(2), cge0(2)
      double precision rgo0(2), cgo0(2)


      double precision pi,twopi,numk
      double precision uxm,uym,uzm
      double precision kx,ky,kx0,ky0,dkx,dky
      double precision ux,uy,uz,dux,duy,duz, kmaxGK,kmaxKKP

      double precision rsp(maxp),csp(maxp),sumr(maxp),sumc(maxp)
      double precision h1(2,2,2), o1(2,2,2), z1(2,2,2), w1(2,11), 
     .                 e1(2)

      double precision ffk(nq0,2) 
      double precision ffkpq(nq0*nq0,2)  
      double precision phiC(maxp), tx(2), ty(2), ekn(nq0,2), 
     .                                           ekpqn(nq0*nq0,2) 

      double precision wkP(-nwm:nwm)
      double precision wkS(-nwm:nwm-1)
      double precision wkT(-nwT:nwT-1)
      double precision kBZ(2,nq0)

      complex*16 H0, Hx, Hy, Hz
      complex*16 expphi, ctemp, ctemp1, ctemp2

* Even Matsubara use (-nwm:nwm) !!
      complex*16 P11(nq0,-nwm:nwm) 
      complex*16 P12(nq0,-nwm:nwm) 
      complex*16 P21(nq0,-nwm:nwm) 
      complex*16 P22(nq0,-nwm:nwm) 
      complex*16 G0(2,2,nq0,-nwm:nwm)
      complex*16 Wc(2,2,nq0,-nwm:nwm)

* Odd Matsubara use (-nwm:nwm-1) !!
      complex*16 ScTau(2,2,nq0,ntaum)

      complex*16 Sigma(2,2,nq0,-nwm:nwm-1)

      complex*16 green11(nq0,-nwm:nwm-1) 
      complex*16 green12(nq0,-nwm:nwm-1) 
      complex*16 green21(nq0,-nwm:nwm-1) 
      complex*16 green22(nq0,-nwm:nwm-1) 

      complex*16 ukn(nq0,2,2), xi, csum, csum1, csum2, csum3, csumx, csumy
      complex*16 csum5, csum6, csum7
      complex*16 ukpqn(nq0*nq0,2,2), b(2,2,8)

      integer nkabc(3), p(4,6), iw1(2)
      integer itest 

* For inversion !!
      dimension rw1(2,2),cw1(2,2)
      dimension rw2(2,2),cw2(2,2)
      dimension rw3(2,2),cw3(2,2)
      dimension rw4(2,2),cw4(2,2)
      dimension rw7(2)

      dimension reps(2,2),ceps(2,2)
      dimension rpi(2,2),cpi(2,2)
      dimension rw(2,2),cw(2,2)

      dimension ww1(2,2),ww2(2,2)
      dimension wwork(2),ipvt(2)


      data p/1,2,4,5,4,5,7,8,2,4,5,7,1,3,4,5,4,5,6,8,3,4,5,6/

      pi    = 4.d0*atan(1.d0)
      twopi = 8.d0*atan(1.d0)
      xi    = dcmplx(0.d0,1.d0)

* Open files.
      open(unit=1,file='gw4.d',status='OLD')
      open(unit=2,file='W11.d',status='UNKNOWN')
      open(unit=3,file='W12.d',status='UNKNOWN')
      open(unit=4,file='W21.d',status='UNKNOWN')
      open(unit=5,file='W22.d',status='UNKNOWN')
      open(unit=11,file='P11.d',status='UNKNOWN')
      open(unit=12,file='P12.d',status='UNKNOWN')
      open(unit=13,file='P21.d',status='UNKNOWN')
      open(unit=14,file='P22.d',status='UNKNOWN')

      open(unit=15,file='S11.d',status='UNKNOWN')
      open(unit=16,file='S12.d',status='UNKNOWN')
      open(unit=17,file='S21.d',status='UNKNOWN')
      open(unit=18,file='S22.d',status='UNKNOWN')

      open(unit=151,file='Stau11.d',status='UNKNOWN')
      open(unit=161,file='Stau12.d',status='UNKNOWN')
      open(unit=171,file='Stau21.d',status='UNKNOWN')
      open(unit=181,file='Stau22.d',status='UNKNOWN')

      open(unit=19,file='G11.d',status='UNKNOWN')
      open(unit=20,file='G12.d',status='UNKNOWN')
      open(unit=21,file='G21.d',status='UNKNOWN')
      open(unit=22,file='G22.d',status='UNKNOWN')

      open(unit=32,file='band_h1.d',status='UNKNOWN')
      open(unit=33,file='dosnos_h1.d',status='UNKNOWN')

      open(unit=49,file='Gtau11.d',status='UNKNOWN')
      open(unit=50,file='Gtau12.d',status='UNKNOWN')
      open(unit=51,file='Gtau21.d',status='UNKNOWN')
      open(unit=52,file='Gtau22.d',status='UNKNOWN')

      open(unit=53,file='Wtau11.d',status='UNKNOWN')
      open(unit=54,file='Wtau12.d',status='UNKNOWN')
      open(unit=55,file='Wtau21.d',status='UNKNOWN')
      open(unit=56,file='Wtau22.d',status='UNKNOWN')

      open(unit=57,file='nume.d',status='UNKNOWN')

      read(1,*)nkabc
      read(1,*)npts,emin,emax
      read(1,*)t1,t2
      read(1,*)phi
      read(1,*)dx,dy
      read(1,*)nphi
      read(1,*)nw
      read(1,*)beta
      read(1,*)utilde, uprim, c4
      read(1,*)ndel, stepdelmu, delmu
      close(1,status='keep')

* Control parameters
      write(6,*)'hgw4.f !!'
      write(6,*)'kdim:',kdim
      write(6,*)'nkabc(i)',nkabc(1),nkabc(2),nkabc(3)
      write(6,'(a20,f10.3)')'t1:',t1
      write(6,'(a20,f10.3)')'t2:',t2
      write(6,'(a20,f10.3)')'phi (units of pi):',phi
      write(6,'(a20,f10.7)')'dx:',dx
      write(6,'(a20,f10.7)')'dy:',dy

      write(6,'(a20,f10.3)')'U tilde:',utilde
      write(6,'(a20,f10.3)')'C:',c4

* Define !!
      uhubb1 = utilde/c4

      write(6,'(a20,f10.3)')'Hubbard U:',uhubb1
      write(6,'(a50,f10.3)')'U prim (A to B):',uprim
      write(6,'(a20,i4)')'nphi:',nphi
      write(6,'(a20,i4)')'nw (Matsubara):',nw
      write(6,'(a30,5f12.6)')'Start delmu (first phi):',delmu
      write(6,'(a30,5f12.6)')'Step in delmu:',stepdelmu
      write(6,'(a30,i4)')'ndel:',ndel




      if(nw.ne.nwm) then
      write(6,*)'nw not equal nwm !!'
      stop
      endif
      if(nwT.ne.nw) then
      write(6,*)'nwT not equal nw !!'
      stop
      endif

      phi = phi*pi
      write(6,'(a20,f10.3)')'phi (in radians):',phi

*      volwgt = ( (8.d0*pi*pi)/(3.d0*dsqrt(3.d0)) )*   !Volyme BZ divided by #k-points
*     .         ( 1.d0/((nkabc(1)+1)*(nkabc(2)+1)) )

*      volwgt = ( 2.d0/(3.d0*dsqrt(3.d0)) )*   !One over volyme times #k-points
*     .         ( 1.d0/((nkabc(1)+1)*(nkabc(2)+1)) )

      volwgt = ( 1.d0 )*   ! #k-points
     .         ( 1.d0/((nkabc(1)+1)*(nkabc(2)+1)) )

      numk   = 1.d0/( (nkabc(1)+1)*(nkabc(2)+1) ) 
 
      numk   = (nkabc(1)+1)*(nkabc(2)+1)

      write(6,*)'numk (number k):',numk
      write(6,'(a20,f10.3)')'volwgt:',volwgt

      det    = (8.d0*pi*pi)/(3.d0*dsqrt(3.d0)) 

      unitcell = 3.d0*dsqrt(3.d0)/2.d0
      write(6,'(a20,f10.4)')'Unit cell area:',unitcell


      uxm=(2.d0/3.d0)*2.d0*pi
      uym=uxm
      dux=dble(uxm/nkabc(1))
      duy=dble(uym/nkabc(2))



* Choose q-point for plotting !!
*      iqp = 1   !q=0
      iqp =  2 


* Zero arrays !!
      do i = 1,nq0
      do j = 1,2
      ekn(i,j)   = 0.d0
      ffk(i,j)   = 0.d0
      enddo
      enddo

      do i = 1,nq0
      do j = 1,2
      do k = 1,2
      ukn(i,j,k)   = dcmplx(0.d0,0.d0)
      enddo
      enddo

      enddo
      do i = 1,nq0*nq0
      do j = 1,2
      ekpqn(i,j) = 0.d0
      ffkpq(i,j) = 0.d0
      enddo
      enddo

      do i = 1,nq0*nq0
      do j = 1,2
      do k = 1,2
      ukpqn(i,j,k) = dcmplx(0.d0,0.d0)
      enddo
      enddo
      enddo

      if(nq0.lt.numk) stop 'Error nq0 !!'

* Choose yy = Delta/t2 !!
      yy = 6.00d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 5.25d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 3.67d0
*      yy = 3.00d0     !Ceresoli PRB74, 024408
*      yy = 1.00d0
*      yy = 0.00d0

* Read in beta instead !!
*      beta = 1.d0/0.05d0  !Inverse temperature; see Ceresoli PRB74, 024408
*      beta = 11594.2/300.d0   !   [1/eV] and T=300 in Kelvin.
*      beta =  20.d0     !Just try

      dw   = 2.d0*pi/beta
 
      eA = -t2*yy
      eB = -eA


      write(6,'(a20,f10.3)')'eA:',eA
      write(6,'(a20,f10.3)')'eB:',eB
      write(6,'(a20,f10.3)')'Delta/t2:',dabs(eA)/t2
      write(6,'(a20,f10.3)')'beta [1/eV]:',beta
      write(6,'(a20,f10.3)')'Temperature [K]:',1.d0/beta
      write(6,'(a60,f10.3)')
     .'Integrals of 4 atomic orbitals over the unit cell:', c4

* tau-mesh !!

* Set here also energy-parameter !!
      niv    = nwT
      write(6,*)'niv:',niv
      write(6,*)'nwT:',nwT

      ntau  = ntaum

      write(6,*)'ntau:',ntau
      write(6,*)'nsimp:',nsimp

      call gentau4 (beta, ntau/nsimp, nsimp,
     o              tau )
      write(6,*)'tau mesh !'
      do itau = 1,ntau
      write(6,'(i4,5f12.6)')itau,tau(itau),beta
      enddo
*      stop


* Test numerically G(iw) to G(tau) !!
      eu  = -0.5d0    !Relative rmu
      ffu = 1.d0/(1.d0+dexp(beta*eu))

      do i  = -niv,niv-1
      wkT(i) = (dfloat(i)+0.5d0)*dw
      enddo

* Exact G(tau).
      do      it = 1, ntau
      cg0(it)= 0.d0
      enddo

      if (eu .gt. 0.d0) then   !Unoccupied state
      do      it = 1, ntau
      rg0(it)= (ffu - 1.d0) * dexp(-eu*tau(it))
      enddo
      endif

      if (eu .le. 0.d0) then   !Occupied state
      do      it = 1, ntau
      rg0(it)= -ffu * dexp(eu*(beta-tau(it)))
      enddo
      endif


      write(6,*)'G0 on imaginary axis !!'
      do      iv = -niv, niv-1
      call g0iv    (wkT(iv), eu,
     o              rgv(iv), cgv(iv) )
      write(6,'(5f12.6)'),wkT(iv),rgv(iv), cgv(iv),
     .1.d0/(xi*wkT(iv)-eu)
      enddo

* Construct even and odd G !!
      do      iv = 0, niv-1
      rge(iv)    = rgv(iv) + rgv(-iv-1)
      cge(iv)    = cgv(iv) + cgv(-iv-1)
      rgo(iv)    = rgv(iv) - rgv(-iv-1)
      cgo(iv)    = cgv(iv) - cgv(-iv-1)
      enddo

* Weight for G(tau) = (1/beta) S[n] exp(ivn*tau) G(ivn).
* This gtau.f remove and add the asymptotic part so we don't need
* gtinf.f !!
      call gtau4   (beta, tau,
     d              ntau, niv,
     o              vi, coswt, sinwt )


      write (*,*)'G0(tau) (from gtau4.f)!!'
      do      it = 1, ntau
      rgtau(it)  = dot_product (coswt(0:niv-1,it), rge)
     .           + dot_product (sinwt(0:niv-1,it), cgo)
      cgtau(it)  = dot_product (coswt(0:niv-1,it), cge)
     .           - dot_product (sinwt(0:niv-1,it), rgo)
      write (6,'(5f12.6)') tau(it), rg0(it), rgtau(it),
     .                              cg0(it), cgtau(it)
      enddo
*      stop

* Infinite correction.
* cos(vn*tau)/beta and sin(vn*tau)/beta
*      do      it = 1, ntau
*      do      iv = 0, niv-1
*      coswt(iv,it) = dcos(vi(iv)*tau(it)) / beta
*      sinwt(iv,it) = dsin(vi(iv)*tau(it)) / beta
*      enddo
*      enddo
*
*      rge0(1)    = rge(niv-1)
*      rge0(2)    = rge(niv-2)
*      cge0(1)    = cge(niv-1)
*      cge0(2)    = cge(niv-2)
*      rgo0(1)    = rgo(niv-1)
*      rgo0(2)    = rgo(niv-2)
*      cgo0(1)    = cgo(niv-1)
*      cgo0(2)    = cgo(niv-2)
*
*      call gtinf   (beta, tau, 1.d-8,
*     i              rge0, cge0, rgo0, cgo0,
*     i              vi, coswt, sinwt,
*     d              ntau, niv,
*     o              rginf, cginf )
*
*      write (*,*)'G0(tau) (with correction)!!'
*      do      it = 1, ntau
*      rgtau(it)  = dot_product (coswt(0:niv-1,it), rge)
*     .           + dot_product (sinwt(0:niv-1,it), cgo)
*      cgtau(it)  = dot_product (coswt(0:niv-1,it), cge)
*     .           - dot_product (sinwt(0:niv-1,it), rgo)
*      write (6,'(5f12.6)') tau(it), rg0(it), rgtau(it)+rginf(it),
*     .                              cg0(it), cgtau(it)+cginf(it)
*      enddo
*
*      stop
* FT G0(tau) -> G0(iv)
      do      it = 1, ntau
      rgt(it) = rg0(it)
      cgt(it) = cg0(it)
      enddo


      write (6,*)'FT G0(tau) -> G0(iv) !!'

      do      iv = -niv, niv-1
      vn         = (2*iv+1)*(pi/beta)

* Fit to exponential function + quadratic instead of polynomial fit
* -> exact for G0
*      call filonx  (vn, tau, rgt, ntau/2,
*     o              rcos, rsin)
*
*      call filonx  (vn, tau, cgt, ntau/2,
*     o              ccos, csin)
*
*      rsum       = rcos - csin
*      csum       = rsin + ccos

* 2nd order Simpson
*      if(iv.eq.-niv) write(6,*)'niv before filong:',niv
*      call filong  (vn, tau, ntau/2,
*     o              wcos, wsin)

* 4th order Simpson
      if(iv.eq.-niv) write(6,*)'niv before filong4:',niv
      call filong4 (vn, tau, ntau/4,
     o              wcos, wsin)

      rsum       = dot_product (rgt, wcos) - dot_product (cgt, wsin)
      csum       = dot_product (rgt, wsin) + dot_product (cgt, wcos)

      write (6,'(i4,5f12.6)') iv,vn,rgv(iv),rsum,cgv(iv),dreal(csum)

      enddo ! iv
      write (6,*)'End FT G0(tau) -> G0(iv) !!'



* New test !!
      write (6,*)'FT f(tau) = tau -> f(iv) !!'

      do      it = 1, ntau
      rgt(it) = tau(it)  !f(tau) = tau
      enddo

      do      iv = -niv, niv-1
      vn         = (2*iv+1)*(pi/beta)


* 4th order Simpson
      call filong4 (vn, tau, ntau/4,
     o              wcos, wsin)

      rsum       = dot_product (rgt, wcos) 
      csum       = dot_product (rgt, wsin)

      rgv(iv) =  beta*dsin(vn*beta)/vn + (dcos(vn*beta)-1.d0)/(vn**2)
      cgv(iv) = -beta*dcos(vn*beta)/vn + (dsin(vn*beta)     )/(vn**2)


      write (6,'(i4,5f12.6)') iv,vn,rgv(iv),rsum,
     .                              cgv(iv),dreal(csum)

      enddo ! iv
      write (6,*)'End FT f(tau) = tau -> f(iv) !!'

*      stop


* DOS.
      dosfac =  ( 1.d0/(nkabc(1)*nkabc(2)*nkabc(3)) )*(1.d0/12.d0)

      ux=-dux

      do ix = 1, nkabc(1) +1

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) +1

      uy=uy+duy
      uz=-duz

      do iz = 1, nkabc(3) +1
      uz=uz+duz


* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)

      call diagno(kdim,h1,o1,w1,iw1,z1,e1)


      do ll = 1, kdim
      e(ll,1) = e1(ll)
      e(ll,4) = e1(ll)
      enddo
      do i1 = 1, kdim
      do i2 = 1, kdim
      b(i1,i2,1) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      b(i1,i2,4) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      enddo
      enddo

* Energy at ux+dux, uy, uz.
* duy = 0.
      dkx =  dux*dsqrt(3.d0)/2.d0
      dky = -dux/2.d0

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1  = dcos( dsqrt(3.d0)*(kx+dkx) - phi)
      rtemp2  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 + 3.d0*(ky+dky)/2.d0 - phi)
      rtemp3  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - 3.d0*(ky+dky)/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1  = dcos( dsqrt(3.d0)*(kx+dkx) + phi)
      rtemp2  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 + 3.d0*(ky+dky)/2.d0 + phi)
      rtemp3  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - 3.d0*(ky+dky)/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky+dky)
      rtemp2     = dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky+dky)
      rtemp2     = dsin(-dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


      call diagno(kdim,h1,o1,w1,iw1,z1,e1)


      do ll = 1, kdim
      e(ll,2) = e1(ll)
      e(ll,7) = e1(ll)
      enddo
      do i1 = 1, kdim
      do i2 = 1, kdim
      b(i1,i2,2) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      b(i1,i2,7) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      enddo
      enddo

* Energy at ux, uy+duy, uz.
      dkx =  0.d0
      dky =  duy

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1  = dcos( dsqrt(3.d0)*(kx+dkx) - phi)
      rtemp2  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 + 3.d0*(ky+dky)/2.d0 - phi)
      rtemp3  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - 3.d0*(ky+dky)/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1  = dcos( dsqrt(3.d0)*(kx+dkx) + phi)
      rtemp2  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 + 3.d0*(ky+dky)/2.d0 + phi)
      rtemp3  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - 3.d0*(ky+dky)/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky+dky)
      rtemp2     = dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky+dky)
      rtemp2     = dsin(-dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      do ll = 1, kdim
      e(ll,3) = e1(ll)
      e(ll,6) = e1(ll)
      enddo
      do i1 = 1, kdim
      do i2 = 1, kdim
      b(i1,i2,3) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      b(i1,i2,6) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      enddo
      enddo


* Energy at ux+dux, uy+duy, uz.
      dkx =  dux*dsqrt(3.d0)/2.d0
      dky = -dux/2.d0 + duy

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1  = dcos( dsqrt(3.d0)*(kx+dkx) - phi)
      rtemp2  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 + 3.d0*(ky+dky)/2.d0 - phi)
      rtemp3  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - 3.d0*(ky+dky)/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1  = dcos( dsqrt(3.d0)*(kx+dkx) + phi)
      rtemp2  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 + 3.d0*(ky+dky)/2.d0 + phi)
      rtemp3  =
     .     dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - 3.d0*(ky+dky)/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky+dky)
      rtemp2     = dcos(-dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky+dky)
      rtemp2     = dsin(-dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*(kx+dkx)/2.d0 - (ky+dky)/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      do ll = 1, kdim
      e(ll,5) = e1(ll)
      e(ll,8) = e1(ll)
      enddo
      do i1 = 1, kdim
      do i2 = 1, kdim
      b(i1,i2,5) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      b(i1,i2,8) = dcmplx(z1(i1,i2,1),z1(i1,i2,2))
      enddo
      enddo

* Construct total DOS and NOS.
      do ll = 1, kdim
      do it = 1, 6
      do ie = 1, 4
      ee(ie)=e(ll,p(ie,it))
      enddo
      wgt = dosfac 
      ebot=dmin1(ee(1),ee(2),ee(3),ee(4))
      etop=dmax1(ee(1),ee(2),ee(3),ee(4))
      if(ebot .lt. emax)then
      call slinz1(wgt,ee,emin,emax,tnos,npts)
      if(etop .gt. emin)
     .call sliny1(wgt,ee,emin,emax,tdos,npts)
      endif
      enddo
      enddo

      enddo
      enddo
      enddo

      de=(emax-emin)/(npts-1)
      do i =1,npts
      freq(i) = emin+de*(i-1)
      enddo



* Fermi energy
      nel = 1.d0  !Number electrons (spinless)

      do i = 1,npts
      if(tnos(i).gt.nel)goto 2222
      enddo
2222  continue
      x1(1) = tnos(i-1)
      x1(2) = tnos(i)
      f1(1) = freq(i-1)
      f1(2) = freq(i)
      ef    = alagr2(dble(nel),x1,f1)
      write(6,*)
      write(6,*)'Number electrons:',nel
      write(6,*)'Haldane model!'
      write(6,'(a20,f10.3)')'Fermi energy (eV):',ef
      write(2,'(a20,f10.3)')'eA (eV):',eA
      write(2,'(a20,f10.3)')'eB (eV):',eB
      write(3,'(a20,f10.3)')'eA (eV):',eA
      write(3,'(a20,f10.3)')'eB (eV):',eB
      write(4,'(a20,f10.3)')'eA (eV):',eA
      write(4,'(a20,f10.3)')'eB (eV):',eB
      write(5,'(a20,f10.3)')'eA (eV):',eA
      write(5,'(a20,f10.3)')'eB (eV):',eB
      write(33,'(a20,f10.3)')'Fermi energy (eV):',ef
      write(6,'(a20,f10.3)')'Fermi energy (Ry):',ef/13.6d0
      write(33,'(3f12.3)')
     .     (freq(i),tdos(i),tnos(i),i=1,npts)
      close(33,status='keep')

* End DOS.

* Bands.
      do ix = 1, nkabc(1) + 1

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) + 1

      uy=uy+duy

* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy

      dkx =  dux*dsqrt(3.d0)/2.d0
      dky = -dux/2.d0 + duy

      rsum = rsum + dabs(dux*duy)*dsqrt(3.d0)/2.d0

      enddo
      enddo

      fact = 8.d0*pi*pi/( 3.d0*dsqrt(3.d0) )   !Area 1BZ
      write(6,'(a30,f10.6)')'Area 1BZ (should be one):',rsum/fact

* G-K.
      uxm=(2.d0/3.d0)*( dsqrt(1.d0/3.d0) )*2*pi
      dkx=dble(uxm/nkabc(1))
      dky=0.d0

      kx=-dkx
      ky=-dky

      do ix = 1, nkabc(1) + 1

      kx=kx+dkx
      ky=ky+dky

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


789   format(14f15.6)

      call diagno(kdim,h1,o1,w1,iw1,z1,e1)


      if(ix.eq.1) write(6,*)' G-K:'
*      if(ix.eq.1) write(6,*)' kx ky k Ek :'
*      write(6,789)kx/twopi,ky/twopi,dsqrt(kx**2 + ky**2)/twopi,
*     .            e1(1),e1(2)

      write(32,789)dsqrt(kx**2 + ky**2)/twopi,e1(1),e1(2)

* End k-loops
      enddo

      write(32,*)
      kmaxGK = dsqrt(kx**2 + ky**2)/twopi

* K-K'.
      uxm =  (1.d0/3.d0)*( dsqrt(1.d0/3.d0) )*2*pi
      dkx = -dble(uxm/nkabc(1))
      dky = -dsqrt(3.d0)*dkx

      kx  = (2.d0/3.d0)*( dsqrt(1.d0/3.d0) )*2*pi  - dkx
      ky  =                                        - dky

      do ix = 1, nkabc(1) + 1


      kx=kx+dkx
      ky=ky+dky


* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)



      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      if(ix.eq.1) write(6,*)' K-Kprim:'
      write(32,789)kmaxGK + dble(ix-1)*dsqrt(dkx**2 + dky**2)/twopi,
     .             e1(1),e1(2)

* End k-loops
      enddo

      write(32,*)
      kmaxKKP = kmaxGK + dble(nkabc(1))*dsqrt(dkx**2 + dky**2)/twopi

* K'-G.
      uxm =  (1.d0/3.d0)*( dsqrt(1.d0/3.d0) )*2*pi
      dkx = -dble(uxm/nkabc(1))
      dky =  dsqrt(3.d0)*dkx

      kx  = (1.d0/3.d0)*( dsqrt(1.d0/3.d0) )*2*pi  - dkx
      ky  = (1.d0/3.d0)*2*pi                       - dky

      do ix = 1, nkabc(1) + 1


      kx=kx+dkx
      ky=ky+dky

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      if(ix.eq.1) write(6,*)' Kprim-G:'
      write(32,789)kmaxKKP + dble(ix-1)*dsqrt(dkx**2 + dky**2)/twopi,
     .             e1(1),e1(2)

* End k-loops
      enddo

* End bands.



* Matsubara mesh !!
* Even for P0 and W !!
* Odd for Sigma !!
      dw = 2.d0*pi/beta
      write(6,'(a45,f10.3)')'Matsubara mesh (even) in eV units!'
* Even Matsubara !!
      do i = -nw,nw
      wkP(i)=dfloat(i)*dw
      write(6,'(i4,2f15.5)')i,wkP(i)
      enddo
      write(6,'(a45,f10.3)')'Matsubara mesh (odd) in eV units!'
* Odd Matsubara !!
      do i = -nw,nw-1
      wkS(i)=(dfloat(i)+0.5d0)*dw
      write(6,'(i4,2f15.5)')i,wkS(i)
      enddo



* phi-mesh !!
      dphi = pi/dble(nphi-1)
      do i = 1,nphi
      phiC(i) = dphi*dble(i-1)
      enddo


      ip0 =  21 
      write(6,*)'ip0:',ip0
      write(6,*)'iqp:',iqp

      if(ip0.gt.nphi) stop 'Error nphi !!'

*      do ip = 1, nphi
      do ip = ip0, ip0 
*      do ip = ip0, ip0 + 1 


      phi = phiC(ip)
      write (117,'(i4,5f12.6)') ip, phi
 
      call cputid(0)


* Find rmu0 for each phi !!
* Gap is always at K-point!
      kx = twopi*(2.d0/3.d0)*(1.d0/dsqrt(3.d0))
      ky = 0.d0

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =   h1(1,2,1)
      h1(2,1,2)  =  -h1(1,2,2)

      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

* Middle in gap !!
* -0.5 scale < 0.5
      scale =  0.48d0 
      gap0 = (e1(2) - e1(1))
*      rmu0 = (e1(2) + e1(1))/2.d0 - 0.5d0*gap0     !rmu top band 1
*      rmu0 = (e1(2) + e1(1))/2.d0 + 0.5d0*gap0     !rmu bottom band 2
*      rmu0 = (e1(2) + e1(1))/2.d0 + scale*gap0     !rmu between 1 and 2  
      rmu0 = (e1(2) + e1(1))/2.d0                  !rmu in the middle

*      rmu0 = rmu0 - 1.d0

      write(6,'(a20,5f10.6)')'phi (units pi)',phi/pi
      write(6,'(a20,5f12.6)')'VB top:',e1(1)
      write(6,'(a20,5f12.6)')'CB bottom:',e1(2)
      write(6,'(a20,5f12.6)')'rmu0:',rmu0
      write(6,'(a20,5f12.6)')'Gap:',gap0
      write(6,*)
      write(6,*)

      write(79,'(5f10.6)'),phi/pi, rmu0, e1(1),e1(2),gap0


* Fix q !!
* q-vectors !!
      uqx=-dux

      do iqx = 1, nkabc(1) +1

      uqx=uqx+dux
      uqy=-duy

      do iqy = 1, nkabc(2) +1

      uqy=uqy+duy


      qx  =  uqx*dsqrt(3.d0)/2.d0
      qy  = -uqx/2.d0 + uqy


* Combined index.
      icc = (iqx-1)*(nkabc(2)+1) + iqy 
      kBZ(1,icc) = qx
      kBZ(2,icc) = qy
      write(6,'(a20,i4,2f10.6)')'qx,qy:',icc,qx/(2.d0*pi),qy/(2.d0*pi)



      ikk = 0 
* Diagonalize for k = (kx,ky).
      ux=-dux

      do ix = 1, nkabc(1) + 1 

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) + 1 

      uy=uy+duy

      ikk = ikk + 1  !Number of k = (nkabc(1)+1)*(nkabc(2)+1)

      if(icc.eq.1)then
* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy


* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA


      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =   h1(1,2,1)
      h1(2,1,2)  =  -h1(1,2,2)



* For k = (kx,ky) !!
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)


      do ib      = 1,kdim  !Band
      ekn(ikk,ib)    = e1(ib) - rmu0 
      ffk(ikk,ib)    = 1.d0/(1.d0+dexp( beta*(e1(ib)-rmu0) ))
      do i1      = 1,kdim
      ukn(ikk,i1,ib) = dcmplx(z1(i1,ib,1),z1(i1,ib,2))
      enddo
      enddo

      write(6,'(a15,i4,5f10.6)')'Eigenvalue:',ikk,ekn(ikk,1),ffk(ikk,1),
     .                                            ekn(ikk,2),ffk(ikk,2)
      endif   !icc=1 (only for one q-point)



* Diagonalize for k + q.

* For kx+qx,ky+qy !!
      kx  =  ux*dsqrt(3.d0)/2.d0 + qx 
      ky  = -ux/2.d0 + uy        + qy

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =   h1(1,2,1)
      h1(2,1,2)  =  -h1(1,2,2)


* Diagonalize for k + q.
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)


* Index k+q !!
      ikpq = (icc-1)*(nkabc(1)+1)*(nkabc(2)+1) + ikk

      do ib      = 1,kdim  !Band
      ekpqn(ikpq,ib)  = e1(ib) - rmu0                        !Used also for Sigma
      ffkpq(ikpq,ib)  = 1.d0/(1.d0+dexp( beta*(e1(ib)-rmu0) ))
      do i1        = 1,kdim
      ukpqn(ikpq,i1,ib) = dcmplx(z1(i1,ib,1),z1(i1,ib,2))   !Used also for Sigma
      enddo
      enddo

*      if(icc.eq.1)then
      write(66,'(a10,i6,5f10.6)')'ffkpq:',ikpq,ekpqn(ikpq,1),
     .                                         ffkpq(ikpq,1),
     .                                         ekpqn(ikpq,2),
     .                                         ffkpq(ikpq,2)
*      endif


      enddo
      enddo     !End k-loops

*      write(66,*)'*'

      nkk = ikk
      if(nkk.ne.(nkabc(1)+1)*(nkabc(2)+1)) stop


* Loop energy !!
* q is fixed given by icc index!!
      do iw = -nw, nw   !Even 


      do it  = 1,2  !tau and tau' 
      do itp = 1,2 


* For P0!
      csum1 = dcmplx(0.d0,0.d0)

      do iq = 1, nkk  !Sum k

      ikpq = (icc-1)*(nkabc(1)+1)*(nkabc(2)+1) + iq


*      do ib  = 1,kdim
*      do ibp = 1,kdim
*      csum1 = csum1+(-1.d0)*dconjg(ukpqn(ikpq,it,ib))*ukn(iq,it,ibp)*  !No 2 spinless electrons
*     .                    dconjg(ukn(iq,itp,ibp))*ukpqn(ikpq,itp,ib)*
*     .(ffkpq(ikpq,ib)-ffk(iq,ibp))
*     ./(xi*wkP(iw)-ekpqn(ikpq,ib)+ekn(iq,ibp))
*      enddo
*      enddo

* Insulator case !!
* The code below gives almost the same as above but no divergence at w=0 !!
* ib=1 and ibp=2 (ffk=0 ffkpq=1) or ib=2 and ibp=1 (ffk=1 ffkpq=0) !!


      csum1 = csum1 + (-1.d0)*dconjg(ukpqn(ikpq,it,1))*ukn(iq,it,2)*
     .                        dconjg(ukn(iq,itp,2))*ukpqn(ikpq,itp,1)*
     .                ( 1.d0)/(xi*wkP(iw)-ekpqn(ikpq,1)+ekn(iq,2))
     .      +
     .                ( 1.d0)*dconjg(ukpqn(ikpq,it,2))*ukn(iq,it,1)*
     .                        dconjg(ukn(iq,itp,1))*ukpqn(ikpq,itp,2)*
     .                ( 1.d0)/(xi*wkP(iw)-ekpqn(ikpq,2)+ekn(iq,1))

      enddo  !Sum k

      if(it.eq.1.and.itp.eq.1) P11(icc,iw) = csum1*volwgt
      if(it.eq.1.and.itp.eq.2) P12(icc,iw) = csum1*volwgt
      if(it.eq.2.and.itp.eq.1) P21(icc,iw) = csum1*volwgt
      if(it.eq.2.and.itp.eq.2) P22(icc,iw) = csum1*volwgt

      enddo  !it and itp loops
      enddo



* Invert for fix q and iw !!
* U !!
*      rw1(1,1) = uhubb1 
*      rw1(1,2) = uprim 
*      rw1(2,1) = uprim 
*      rw1(2,2) = uhubb1 
*      cw1(1,1) = 0.d0 
*      cw1(1,2) = 0.d0 
*      cw1(2,1) = 0.d0 
*      cw1(2,2) = 0.d0 

      rw1(1,1) = uhubb1 
      rw1(1,2) = uprim*dreal( cdexp(-xi*qy) ) 
      rw1(2,1) = uprim*dreal( cdexp( xi*qy) ) 
      rw1(2,2) = uhubb1 
      cw1(1,1) = 0.d0 
      cw1(1,2) = uprim*dimag( cdexp(-xi*qy) ) 
      cw1(2,1) = uprim*dimag( cdexp( xi*qy) ) 
      cw1(2,2) = 0.d0 


      rw3(1,1) = dreal(P11(icc,iw))
      cw3(1,1) = dimag(P11(icc,iw))

      rw3(1,2) = dreal(P12(icc,iw))
      cw3(1,2) = dimag(P12(icc,iw))

      rw3(2,1) = dreal(P21(icc,iw))
      cw3(2,1) = dimag(P21(icc,iw))

      rw3(2,2) = dreal(P22(icc,iw))
      cw3(2,2) = dimag(P22(icc,iw))

* Put ImP = 0.
*      cw3(1,1) = dcmplx(0.d0,0.d0)
*      cw3(1,2) = dcmplx(0.d0,0.d0)
*      cw3(2,1) = dcmplx(0.d0,0.d0)
*      cw3(2,2) = dcmplx(0.d0,0.d0)


* v*P0  (P0 total)
      call mmulc   (rw1,cw1,2,rw3,cw3,
     i              2,2,2,2,2,
     o              reps,ceps)
* eps = (1 - v*P0)
      call cv       (-1.d0,reps,2*2,
     o               reps )
      call cv       (-1.d0,ceps,2*2,
     o               ceps )
      do      i = 1,2
      reps(i,i) = 1.d0 + reps(i,i)
      enddo
* eps^(-1) = [1 - v*P0]^(-1)
      call minvc   (reps,ceps,
     d              2,2,
     w              wwork,ipvt,ww1,ww2,
     o              rpi,cpi)


* Full W = eps^{-1}*v (stored in rw, cw)
      call mmulc   (rpi,cpi,kdim,rw1,cw1,
     i              kdim,kdim,kdim,kdim,kdim,
     o              rw,cw)

* Wc=rw-rw1.
      Wc(1,1,icc,iw) = dcmplx( rw(1,1), cw(1,1)) - 
     .                 dcmplx(rw1(1,1),cw1(1,1))
      Wc(1,2,icc,iw) = dcmplx( rw(1,2), cw(1,2)) - 
     .                 dcmplx(rw1(1,2),cw1(1,2))
      Wc(2,1,icc,iw) = dcmplx( rw(2,1), cw(2,1)) - 
     .                 dcmplx(rw1(2,1),cw1(2,1))
      Wc(2,2,icc,iw) = dcmplx( rw(2,2), cw(2,2)) - 
     .                 dcmplx(rw1(2,2),cw1(2,2))

* Make Wc real !!
* Only use real part when going to time space so don't need this !!
      Wc(1,1,icc,iw) = dreal(Wc(1,1,icc,iw))
      Wc(1,2,icc,iw) = dreal(Wc(1,2,icc,iw)) 
      Wc(2,1,icc,iw) = dreal(Wc(2,1,icc,iw))
      Wc(2,2,icc,iw) = dreal(Wc(2,2,icc,iw))


* Choose q-point here for W to print !!
* icc=iqp=1 i.e q = 0 !!
      if(icc.eq.iqp.and.ip.eq.ip0)then
      if(iw.eq.-nw) then
      write(2,*)'W element 11 (hgw4.f) !!'
      write(3,*)'W element 12 (hgw4.f) !!'
      write(4,*)'W element 21 (hgw4.f) !!'
      write(5,*)'W element 22 (hgw4.f) !!'
      write(11,*)'P0 element 11 (hgw4.f) !!'
      write(12,*)'P0 element 12 (hgw4.f) !!'
      write(13,*)'P0 element 21 (hgw4.f) !!'
      write(14,*)'P0 element 22 (hgw4.f) !!'
      write(2,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(3,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(4,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(5,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(11,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(12,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(13,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(14,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      endif
      write(2,'(3f12.6)'),wkP(iw), Wc(1,1,icc,iw)
      write(3,'(3f12.6)'),wkP(iw), Wc(1,2,icc,iw)
      write(4,'(3f12.6)'),wkP(iw), Wc(2,1,icc,iw)
      write(5,'(3f12.6)'),wkP(iw), Wc(2,2,icc,iw)
      write(11,'(4f12.6)'),wkP(iw),P11(icc,iw)
      write(12,'(3f12.6)'),wkP(iw),P12(icc,iw)
      write(13,'(3f12.6)'),wkP(iw),P21(icc,iw)
      write(14,'(3f12.6)'),wkP(iw),P22(icc,iw)
      endif


      enddo     !End iw-loop


      enddo     
      enddo     !End  q-loops

      write(6,*)'nkk for W(tau):',nkk

      write(6,*)'Done P0 and full W in iw-space!'
*      stop

* Make Wc exactly even !!
*      do iw = 1, nw   !Even
*      Wc(1,1,icc,-iw) = Wc(1,1,icc,iw)
*      Wc(1,2,icc,-iw) = Wc(1,2,icc,iw)
*      Wc(2,1,icc,-iw) = Wc(2,1,icc,iw)
*      Wc(2,2,icc,-iw) = Wc(2,2,icc,iw)
*      enddo

* Check symmetry !!
* Maybe ImW = 0 ?? Gets smaller if number k is increased !!
      tol    = 0.0000000001d0
      do  iv = 1, nw
      rtemp1 = dreal(Wc(1,1,iqp,-iv)) - dreal(Wc(1,1,iqp,iv))  !Even
      rtemp2 = dimag(Wc(1,1,iqp,-iv)) + dimag(Wc(1,1,iqp,iv))  !Odd 
      if(dabs(rtemp1).gt.tol) stop 'Re(W11) error in symmetry!!'
      if(dabs(rtemp2).gt.tol) stop 'Im(W11) error in symmetry!!'
      rtemp1 = dreal(P11(iqp,-iv)) - dreal(P11(iqp,iv))  !Even
      rtemp2 = dimag(P11(iqp,-iv)) + dimag(P11(iqp,iv))  !Odd 
      if(dabs(rtemp1).gt.tol) stop 'Re(P11) error in symmetry!!'
      if(dabs(rtemp2).gt.tol) stop 'Im(P11) error in symmetry!!'
      enddo
      do  iv = 1, nw
      rtemp1 = dreal(Wc(2,2,iqp,-iv)) - dreal(Wc(2,2,iqp,iv))  !Even
      rtemp2 = dimag(Wc(2,2,iqp,-iv)) + dimag(Wc(2,2,iqp,iv))  !Odd
      if(dabs(rtemp1).gt.tol) stop 'Re(W22) error in symmetry!!'
      if(dabs(rtemp2).gt.tol) stop 'Im(W22) error in symmetry!!'
      rtemp1 = dreal(P22(iqp,-iv)) - dreal(P22(iqp,iv))  !Even
      rtemp2 = dimag(P22(iqp,-iv)) + dimag(P22(iqp,iv))  !Odd 
      if(dabs(rtemp1).gt.tol) stop 'Re(P22) error in symmetry!!'
      if(dabs(rtemp2).gt.tol) stop 'Im(P22) error in symmetry!!'
      enddo

* W12 and W21 no symmetry !!
* Almost symmetry !!
* Re(W12) and Re(W21) almost even (with tol1=0.00001d0) !!
* Im(W12) and Im(W21) almost even (with tol1=0.01d0) !!
*      tol1   = 0.00001d0
      tol1   = 0.0001d0
      do  iv = 1, nw
      rtemp1 = dreal(Wc(1,2,iqp,-iv)) - dreal(Wc(1,2,iqp,iv))  !Even
      rtemp2 = dimag(Wc(1,2,iqp,-iv)) + dimag(Wc(1,2,iqp,iv))  !Odd
      if(dabs(rtemp1).gt.tol1) stop 'Re(W12) error in symmetry!!'
      if(dabs(rtemp2).gt.tol1) stop 'Im(W12) error in symmetry!!'
      enddo

      do  iv = 1, nw
      rtemp1 = dreal(Wc(2,1,iqp,-iv)) - dreal(Wc(2,1,iqp,iv))  !Even
      rtemp2 = dimag(Wc(2,1,iqp,-iv)) + dimag(Wc(2,1,iqp,iv))  !Odd
      if(dabs(rtemp1).gt.tol1) stop 'Re(W21) error in symmetry!!'
      if(dabs(rtemp2).gt.tol1) stop 'Im(W21) error in symmetry!!'
      enddo



* Do W(tau) for all k!!
* Weight for W(tau) = (1/beta) S[n] exp(iwn*tau) W(iwn).
* iwn is even !!
* Routine wtau.f uses that W(iwn) is real and even i.e W(iwn) = W(-iwn) !!
* This means that W in tau-space is a real function.
      call wtau    (beta, tau,
     d              ntau, nw,
     o              wvi, wcoswt )

*      write (*,*)'Even iw mesh (wtau.f)!!'
*      do      iv = 0, nw
*      write(6,'(i4,5f12.6)'),iv,wkP(iv),wvi(iv),(wkP(iv)-wvi(iv))
*      enddo

* Find Wk(tau) !!
* Assume real !!
      do ik      = 1,nkk  !Loop all k

      do i1      = 1,4

      if(i1.eq.1) then

* 11.
      do      iv = 0, nw
      rwe(iv)    = dreal(Wc(1,1,ik,iv))   !ImW odd but small so skip
      enddo

      if(iqp.eq.ik.and.ip.eq.ip0) write (53,*)'W11(tau) !!'
      do      it = 1, ntau

      rsum1      = 0.d0
      do      iv = 0, nw  !Also iw=0
      rsum1      = rsum1 + wcoswt(iv,it)*rwe(iv) 
      enddo
      
      wktau(1,1,ik,it) = dcmplx(rsum1,0.d0) 

      if(iqp.eq.ik.and.ip.eq.ip0) 
     .write (53,'(5f12.5)') tau(it), rsum1 
      enddo

      endif

      if(i1.eq.2) then

* 12.
      do      iv = 0, nw
      rwe(iv)    = dreal(Wc(1,2,ik,iv))
      enddo

      if(iqp.eq.ik.and.ip.eq.ip0) write (54,*)'W12(tau) !!'
      do      it = 1, ntau

      rsum1      = 0.d0
      do      iv = 0, nw  !Also iw=0
      rsum1      = rsum1 + wcoswt(iv,it)*rwe(iv)
      enddo

      wktau(1,2,ik,it) = dcmplx(rsum1,0.d0) 

      if(iqp.eq.ik.and.ip.eq.ip0) 
     .write (54,'(5f12.5)') tau(it), rsum1

      enddo
      endif

      if(i1.eq.3) then

* 21.
      do      iv = 0, nw
      rwe(iv)    = dreal(Wc(2,1,ik,iv))
      enddo

      if(iqp.eq.ik.and.ip.eq.ip0) write (55,*)'W21(tau) !!'
      do      it = 1, ntau

      rsum1      = 0.d0
      do      iv = 0, nw  !Also iw=0
      rsum1      = rsum1 + wcoswt(iv,it)*rwe(iv)
      enddo

      wktau(2,1,ik,it) = dcmplx(rsum1,0.d0) 

      if(iqp.eq.ik.and.ip.eq.ip0) 
     .write (55,'(5f12.5)') tau(it), rsum1 

      enddo
      endif

      if(i1.eq.4) then

* 22.
      do      iv = 0, nw
      rwe(iv)    = dreal(Wc(2,2,ik,iv))   !ImW odd but small so skip
      enddo

      if(iqp.eq.ik.and.ip.eq.ip0) write (56,*)'W22(tau) !!'
      do      it = 1, ntau

      rsum1      = 0.d0
      do      iv = 0, nw  !Also iw=0
      rsum1      = rsum1 + wcoswt(iv,it)*rwe(iv)
      enddo

      wktau(2,2,ik,it) = dcmplx(rsum1,0.d0) 

      if(iqp.eq.ik.and.ip.eq.ip0) 
     .write (56,'(5f12.5)') tau(it), rsum1

      enddo
      endif

      enddo

      enddo  !k loop

      if(icc.gt.nq0) stop

      call cputid(0)

      write(6,*)'Done P0 and full W in tau-space!'
*      stop



      write(6,*)'Start correlation Sigma!'

* Correlation !!
* Fix q !
      uqx=-dux

      do iqx = 1, nkabc(1) +1

      uqx=uqx+dux
      uqy=-duy

      do iqy = 1, nkabc(2) +1

      uqy=uqy+duy

* Combined index.
      iccq = (iqx-1)*(nkabc(2)+1) + iqy   !For q

      qx   =  uqx*dsqrt(3.d0)/2.d0
      qy   = -uqx/2.d0 + uqy
      write(6,'(a20,i4,2f10.6)')'qx,qy:',iccq,qx/(2.d0*pi),qy/(2.d0*pi)


      do ib   = 1,kdim
      do ibp  = 1,kdim
      do itt  = 1,ntau  !time



      csum1 = dcmplx(0.d0,0.d0)
* Sum k, band im and basis tau tau'!
      do ik =  1, nkk   !Sum k

* Index k+q
      ikpq = (iccq-1)*(nkabc(1)+1)*(nkabc(2)+1) + ik

*      do im  = 1,kdim
      do it  = 1,kdim
      do itp = 1,kdim

*      csum1 = csum1 + dconjg(ukpqn(ikpq,itp,im))*ukn(iccq,itp,ibp)*
*     .                dconjg(ukn(iccq,it,ib))*ukpqn(ikpq,it,im)*
*     .        wktau(it,itp,ik,itt)*
*     .       (ffkpq(ikpq,im)-1.d0)*dexp(-ekpqn(ikpq,im)*tau(itt) )

* im = 1 and im=2 seperately !!
      csum1 = csum1 + dconjg(ukpqn(ikpq,itp,1))*ukn(iccq,itp,ibp)*
     .                dconjg(ukn(iccq,it,ib))*ukpqn(ikpq,it,1)*
     .        wktau(it,itp,ik,itt)*
     .       (-ffkpq(ikpq,1))*dexp( ekpqn(ikpq,1)*(beta-tau(itt)) )   !Occupied
     .
     .              + dconjg(ukpqn(ikpq,itp,2))*ukn(iccq,itp,ibp)*
     .                dconjg(ukn(iccq,it,ib))*ukpqn(ikpq,it,2)*
     .        wktau(it,itp,ik,itt)*
     .       (ffkpq(ikpq,2)-1.d0)*dexp(-ekpqn(ikpq,2)*tau(itt) )      !Unoccupied

* ffkpq 0 or one exactly !!
*      csum1 = csum1 + dconjg(ukpqn(ikpq,itp,2))*ukn(iccq,itp,ibp)*
*     .                dconjg(ukn(iccq,it,ib))*ukpqn(ikpq,it,2)*
*     .        wktau(it,itp,ik,itt)*
*     .       (-1.d0)*dexp(-ekpqn(ikpq,2)*tau(itt) )


*      enddo   !m sum
      enddo
      enddo   

      enddo   !Sum k


      ScTau(ib,ibp,iccq,itt)    = -csum1*volwgt

      enddo    !time loop
      enddo
      enddo    !Band ib and ibp

*      if(iccq.eq.iqp) goto 69     !If only iqp points to consider

      enddo
      enddo     !End  q-loops  (iccq index)

      write(6,*)'Number q for new G:',iccq
      nkk = iccq
      if(nkk.ne.(nkabc(1)+1)*(nkabc(2)+1)) stop

      if(iccq.gt.nq0) stop
69    continue

      if(ip.eq.ip0)then
      write(151,*)'Sigma(tau) element 11 (hgw4.f) !!'
      write(161,*)'Sigma(tau) element 12 (hgw4.f) !!'
      write(171,*)'Sigma(tau) element 21 (hgw4.f) !!'
      write(181,*)'Sigma(tau) element 22 (hgw4.f) !!'
      write(151,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      write(161,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      write(171,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      write(181,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      do itt = 1,ntau
      write(151,'(5f12.5)')tau(itt),ScTau(1,1,iqp,itt)
      write(161,'(5f12.5)')tau(itt),ScTau(1,2,iqp,itt)
      write(171,'(5f12.5)')tau(itt),ScTau(2,1,iqp,itt)
      write(181,'(5f12.5)')tau(itt),ScTau(2,2,iqp,itt)
      enddo
      endif

      call cputid(0)

      write(6,*)'Done Sigma!'
*      stop

* Find Sc(iv) from Sc(tau)!!
* Use only Re(Sc) ???

      do iq  = 1, nkk     !If all q-points (same as for P0)

      do itt = 1,ntau

      rgt11(itt) = dreal(ScTau(1,1,iq,itt))
*      cgt11(itt) = dimag(ScTau(1,1,iq,itt))
      cgt11(itt) = dcmplx(0.d0,0.d0) 

      rgt12(itt) = dreal(ScTau(1,2,iq,itt))
*      cgt12(itt) = dimag(ScTau(1,2,iq,itt))
      cgt12(itt) = dcmplx(0.d0,0.d0) 

      rgt21(itt) = dreal(ScTau(2,1,iq,itt))
*      cgt21(itt) = dimag(ScTau(2,1,iq,itt))
      cgt21(itt) = dcmplx(0.d0,0.d0) 

      rgt22(itt) = dreal(ScTau(2,2,iq,itt))
*      cgt22(itt) = dimag(ScTau(2,2,iq,itt))
      cgt22(itt) = dcmplx(0.d0,0.d0) 

      enddo


      do      iv = -niv, niv-1

* Weights for each energy vn! 
      vn         = (2*iv+1)*(pi/beta)


* 2nd order Simpson !!
*      call filong  (vn, tau, ntau/2,
*     o              wcos, wsin)
* 4th order Simpson !!
      call filong4 (vn, tau, ntau/4,
     o              wcos, wsin)

* 11.
      rsum       = dot_product (rgt11,wcos) - dot_product(cgt11,wsin)
      csum       = dot_product (rgt11,wsin) + dot_product(cgt11,wcos)

      Sigma(1,1,iq,iv)  = dcmplx(rsum,dreal(csum))

* 12.
      rsum       = dot_product (rgt12,wcos) - dot_product(cgt12,wsin)
      csum       = dot_product (rgt12,wsin) + dot_product(cgt12,wcos)

      Sigma(1,2,iq,iv)  = dcmplx(rsum,dreal(csum))

* 21.
      rsum       = dot_product (rgt21,wcos) - dot_product(cgt21,wsin)
      csum       = dot_product (rgt21,wsin) + dot_product(cgt21,wcos)

      Sigma(2,1,iq,iv)  = dcmplx(rsum,dreal(csum))

* 22.
      rsum       = dot_product (rgt22,wcos) - dot_product(cgt22,wsin)
      csum       = dot_product (rgt22,wsin) + dot_product(cgt22,wcos)

      Sigma(2,2,iq,iv)  = dcmplx(rsum,dreal(csum))

      enddo ! iv

      enddo  !iq loop

      if(ip.eq.ip0)then
      write(15,*)'Sigma element 11 (hgw4.f) !!'
      write(16,*)'Sigma element 12 (hgw4.f) !!'
      write(17,*)'Sigma element 21 (hgw4.f) !!'
      write(18,*)'Sigma element 22 (hgw4.f) !!'
      write(15,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      write(16,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      write(17,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi
      write(18,'(a20,2f12.6)')'q:',kBZ(1,iqp)/twopi,kBZ(2,iqp)/twopi

      do iwl = -nw,nw-1
      write(15,'(5f12.6)')wkS(iwl),Sigma(1,1,iqp,iwl)
      write(16,'(5f12.6)')wkS(iwl),Sigma(1,2,iqp,iwl)
      write(17,'(5f12.6)')wkS(iwl),Sigma(2,1,iqp,iwl)
      write(18,'(5f12.6)')wkS(iwl),Sigma(2,2,iqp,iwl)
      enddo
      endif

*      stop
      
* New G.
      write(6,*)'Start new G!'

* Loop over delmu !!
* First guess for delm is the read from file !!
      do idel  = 1, ndel    !If only iqp q-points


*      do iq  = 1, iqp    !If only iqp q-points
      do iq  = 1, nkk     !If all q-points (same as for P0)
      do iw  = -nw, nw-1


      rw3(1,1) = dreal(xi*wkS(iw)-ekn(iq,1)-Sigma(1,1,iq,iw)+delmu)    !ekn from above (same q-points)
      cw3(1,1) = dimag(xi*wkS(iw)-ekn(iq,1)-Sigma(1,1,iq,iw)+delmu)
      rw3(1,2) = dreal(                    -Sigma(1,2,iq,iw))
      cw3(1,2) = dimag(                    -Sigma(1,2,iq,iw))
      rw3(2,1) = dreal(                    -Sigma(2,1,iq,iw))
      cw3(2,1) = dimag(                    -Sigma(2,1,iq,iw))
      rw3(2,2) = dreal(xi*wkS(iw)-ekn(iq,2)-Sigma(2,2,iq,iw)+delmu)
      cw3(2,2) = dimag(xi*wkS(iw)-ekn(iq,2)-Sigma(2,2,iq,iw)+delmu)

* Invert (G0^{-1} - Sigma + delmu).
      call minvc   (rw3,cw3,
     d              2,2,
     w              wwork,ipvt,ww1,ww2,
     o              rpi,cpi)


      green11(iq,iw) = dcmplx(rpi(1,1),cpi(1,1))
      green12(iq,iw) = dcmplx(rpi(1,2),cpi(1,2))
      green21(iq,iw) = dcmplx(rpi(2,1),cpi(2,1))
      green22(iq,iw) = dcmplx(rpi(2,2),cpi(2,2))


      enddo  !iw loop
      enddo  !iq loop

      write(6,*)'Done new G(q,iw)!'

* Choose q-point here for G to print !!
* icc=1 i.e q = 0 !!
      icc =iqp

      if(ip.eq.ip0)then
      write(19,*)'New G element 11 (hgw4.f) !!'
      write(20,*)'New G element 12 (hgw4.f) !!'
      write(21,*)'New G element 21 (hgw4.f) !!'
      write(22,*)'New G element 22 (hgw4.f) !!'
      write(19,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(20,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(21,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi
      write(22,'(a20,2f12.6)')'q:',kBZ(1,icc)/twopi,kBZ(2,icc)/twopi

      do iw = -nw,nw-1
      write(19,'(5f12.6)'),wkS(iw), green11(icc,iw),
     .1.d0/(xi*wkS(iw)-ekn(icc,1))
      write(20,'(5f12.6)'),wkS(iw), green12(icc,iw)
      write(21,'(5f12.6)'),wkS(iw), green21(icc,iw)
      write(22,'(5f12.6)'),wkS(iw), green22(icc,iw),
     .1.d0/(xi*wkS(iw)-ekn(icc,2))
      enddo
      endif
      call cputid(0)
*      stop

* New G(q,0) all q !!

      rtemp1 = 0.d0
      rtemp2 = 0.d0
      rtemp3 = 0.d0
      do iq  = 1, nkk
*      do iq  = 1, iqp    !If only iqp q-points

      eu  = ekn(iq,1)    !Band 1
*      eu  = ekn(iq,2)    !Band 2
      ffu = 1.d0/(1.d0+dexp(beta*eu))

* Construct even and odd G !!
      do      iv = 0, niv-1

      rge11(iv)    = dreal( green11(iq,iv) + green11(iq,-iv-1) )
      cge11(iv)    = dimag( green11(iq,iv) + green11(iq,-iv-1) )
      rgo11(iv)    = dreal( green11(iq,iv) - green11(iq,-iv-1) )
      cgo11(iv)    = dimag( green11(iq,iv) - green11(iq,-iv-1) )

      rge12(iv)    = dreal( green12(iq,iv) + green12(iq,-iv-1) )
      cge12(iv)    = dimag( green12(iq,iv) + green12(iq,-iv-1) )
      rgo12(iv)    = dreal( green12(iq,iv) - green12(iq,-iv-1) )
      cgo12(iv)    = dimag( green12(iq,iv) - green12(iq,-iv-1) )

      rge21(iv)    = dreal( green21(iq,iv) + green21(iq,-iv-1) )
      cge21(iv)    = dimag( green21(iq,iv) + green21(iq,-iv-1) )
      rgo21(iv)    = dreal( green21(iq,iv) - green21(iq,-iv-1) )
      cgo21(iv)    = dimag( green21(iq,iv) - green21(iq,-iv-1) )

      rge22(iv)    = dreal( green22(iq,iv) + green22(iq,-iv-1) )
      cge22(iv)    = dimag( green22(iq,iv) + green22(iq,-iv-1) )
      rgo22(iv)    = dreal( green22(iq,iv) - green22(iq,-iv-1) )
      cgo22(iv)    = dimag( green22(iq,iv) - green22(iq,-iv-1) )

      enddo


      write(6,'(a20,2f12.6)')'qx,qy:',kBZ(1,iq)/twopi,kBZ(2,iq)/twopi
      write(6,'(a20,2f12.6)')'ek band 1:',ekn(iq,1)
*      write(6,'(a20,2f12.6)')'ek band 2:',ekn(iq,2)


* Weight for G(tau) = (1/beta) S[n] exp(ivn*tau) G(ivn).
* This gtau.f remove and add the asymptotic part !!
* Thus don't use gtinf.f !!
      call gtau4   (beta, tau,
     d              ntau, niv,
     o              vi, coswt, sinwt )

      write (6,*)'G(tau) element 11 !!'
*      write (6,*)'G(tau) element 12 !!'
*      write (6,*)'G(tau) element 21 !!'
*      write (6,*)'G(tau) element 22 !!'

      do      it = 1, ntau

* Exact G0(tau).
      cg0(it)= 0.d0
      if (eu .gt. 0.d0) then
      rg0(it)= (ffu - 1.d0) * dexp(-eu*tau(it))
      endif
      if (eu .le. 0.d0) then
      rg0(it)= -ffu * dexp(eu*(beta-tau(it)))
      endif

      rgtau11(it)  = dot_product (coswt(0:nw-1,it), rge11)
     .             + dot_product (sinwt(0:nw-1,it), cgo11)
      cgtau11(it)  = dot_product (coswt(0:nw-1,it), cge11)
     .             - dot_product (sinwt(0:nw-1,it), rgo11)

      rgtau12(it)  = dot_product (coswt(0:nw-1,it), rge12)
     .             + dot_product (sinwt(0:nw-1,it), cgo12)
      cgtau12(it)  = dot_product (coswt(0:nw-1,it), cge12)
     .             - dot_product (sinwt(0:nw-1,it), rgo12)

      rgtau21(it)  = dot_product (coswt(0:nw-1,it), rge21)
     .             + dot_product (sinwt(0:nw-1,it), cgo21)
      cgtau21(it)  = dot_product (coswt(0:nw-1,it), cge21)
     .             - dot_product (sinwt(0:nw-1,it), rgo21)

      rgtau22(it)  = dot_product (coswt(0:nw-1,it), rge22)
     .             + dot_product (sinwt(0:nw-1,it), cgo22)
      cgtau22(it)  = dot_product (coswt(0:nw-1,it), cge22)
     .             - dot_product (sinwt(0:nw-1,it), rgo22)

      
      G0(1,1,iq,it) = dcmplx(rgtau11(it),cgtau11(it))
      G0(1,2,iq,it) = dcmplx(rgtau12(it),cgtau12(it))
      G0(2,1,iq,it) = dcmplx(rgtau21(it),cgtau21(it))
      G0(2,2,iq,it) = dcmplx(rgtau22(it),cgtau22(it))


**      write (6,'(5f12.6)') tau(it), rgtau11(it), rg0(it),
**     .                              cgtau11(it), cg0(it)


*      write (6,'(5f12.6)') tau(it), rgtau12(it), rg0(it),
*     .                              cgtau12(it), cg0(it)
*      write (6,'(5f12.6)') tau(it), rgtau21(it), rg0(it),
*     .                              cgtau21(it), cg0(it)
*      write (6,'(5f12.6)') tau(it), rgtau22(it), rg0(it),
*     .                              cgtau22(it), cg0(it)


      enddo    !tau loop

* Infinite correction.
* cos(vn*tau)/beta and sin(vn*tau)/beta
*      do      it = 1, ntau
*      do      iv = 0, niv-1
*      coswt(iv,it) = dcos(vi(iv)*tau(it)) / beta
*      sinwt(iv,it) = dsin(vi(iv)*tau(it)) / beta
*      enddo
*      enddo
*
*      rge0(1)    = rge11(niv-1)
*      rge0(2)    = rge11(niv-2)
*      cge0(1)    = cge11(niv-1)
*      cge0(2)    = cge11(niv-2)
*      rgo0(1)    = rgo11(niv-1)
*      rgo0(2)    = rgo11(niv-2)
*      cgo0(1)    = cgo11(niv-1)
*      cgo0(2)    = cgo11(niv-2)
*
*      call gtinf   (beta, tau, 1.d-8,
*     i              rge0, cge0, rgo0, cgo0,
*     i              vi, coswt, sinwt,
*     d              ntau, niv,
*     o              rginf, cginf )
*
*      write(6,*)'Infinite correction for new G(tau) element 11 !!'
*      do      it = 1, ntau
*
*      rgtau11(it)  = dot_product (coswt(0:nw-1,it), rge11)
*     .             + dot_product (sinwt(0:nw-1,it), cgo11)
*      cgtau11(it)  = dot_product (coswt(0:nw-1,it), cge11)
*     .             - dot_product (sinwt(0:nw-1,it), rgo11)
*
*      write (6,'(5f12.6)') tau(it), rgtau11(it)+rginf(it), 
*     .                              cgtau11(it)+cginf(it)
*
*      enddo
***


      write (117,'(i4,6f12.6)') iq, kBZ(1,iq),kBZ(2,iq),
     .                         -rgtau11(ntau),  
     .                         -rgtau12(ntau),
     .                         -rgtau21(ntau),
     .                         -rgtau22(ntau)

*      rtemp1 = rtemp1 + (-rgtau11(ntau))
*      rtemp2 = rtemp2 + (-rgtau22(1))

      rtemp1 = rtemp1 + (-rgtau11(ntau))
      rtemp2 = rtemp2 + (-rgtau22(1))
      rtemp3 = rtemp3 - (rgtau11(ntau)+rgtau22(ntau))

      enddo  !q-loop

      write(6,'(a40,5f12.2)')'Initial delmu:',delmu

      if(rtemp3/dble(numk).gt.1.d0) delmu = delmu - stepdelmu
      if(rtemp3/dble(numk).lt.1.d0) delmu = delmu + stepdelmu 

      write(6,*)'idel:',idel
      write(6,'(a40,5f12.2)')'New delmu:',delmu
      write(6,'(a40,5f12.4)')'Sum of G11(tau=beta) and G22(tau=beta)!',
     .                        rtemp3/dble(numk)

      told = 0.005d0
      if(rtemp3/dble(numk).lt.(1.d0+told).and.
     .   rtemp3/dble(numk).gt.(1.d0-told)) goto 99

      enddo  !delmu-loop
99    continue


      write(6,*)
      write(6,*)

      write(6,*)'Number iterations for delmu:',idel
      write(6,'(a40,5f12.2)')'phi:',phiC(ip)/pi
      write(6,'(a40,i5)')'phi index:',ip
      write(6,'(a40,5f12.2)')'Sum of -G11(tau=0)!',
     .                        rtemp1/dble(numk)
      write(6,*)
      write(6,'(a40,5f12.2)')'Sum of -G22(tau=0)!',
     .                        rtemp2/dble(numk)
      write(6,*)
      write(6,'(a40,5f12.2)')'Sum of G11(tau=beta) and G22(tau=beta)!',
     .                        rtemp3/dble(numk)

      write(57,'(2i4,5f12.2)')ip,idel,phiC(ip)/pi,
     .                                      rtemp1/dble(numk),
     .                                      rtemp2/dble(numk),
     .                                      rtemp3/dble(numk)

      call cputid(0)
      write(6,*)'Done new G(q,0) for all bands!'

      enddo     !End phi-loop


      stop
      end

      double precision function alagr2 (x,xi,fi)

c 92.03.02
c 92.04.10 from alagr3
c two-point interpolation
c given a function fi at two points xi, the routine interpolates
c the function at x
c f(x) = [ (x-x2)/(x1-x2) ] f1
c      + [ (x-x1)/(x2-x1) ] f2

c x  = the point at which the function is to be interpolated
c xi(2) = points where the function is given
c fi(2) = the function at xi

      implicit real*8 (a-h,o-z)
      dimension xi(2),fi(2)

      xx1        = x-xi(1)
      xx2        = x-xi(2)
      x12        = xi(1)-xi(2)
      alagr2     = (xx2*fi(1) - xx1*fi(2))/x12

      return
      end                       
