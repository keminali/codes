c Unsteady BEM programme to calculate C_T, C_Q, C_P of a rotor of HAWT or HAMT
c TUDK axisymmetric unsteady BEM method [1,7]
C Unsteadiness is determined  by changes in the blade's pitch angle and/or tip speed ratio (TSR)
c Blade's pitch angle time variation is determined by the subroutine find_theta
c TSR time variation is determined by the subroutine find_tsr
c The user should modify those subroutines to fit purposes.
c 
c Steady BEM Options
c
c High axial correction factor, CT is of
c ig = 0 Spera (1994) Eq (6.38) [1], ig =1 Glauert Eq (6.37) [1]
c
C Tip losses:
c iprant=1 - Prandtl of tip and hub Eq. (7,8,13) [2], (8,9) [3]
c iprant=0 - Goldstein approx Eq (A11,A12) [4]
c
c Rotational effect on cl: Snel's approximation Eq. (4.2) [5]
c
c Post stall cl and cd: Viterna equations Eq. (4-10) [6]
c
C Unsteady BEM (TUDK wake model):
c wake model: Eqs.(9.18) to (9.22) [1]
c Unsteady lift&drag, see dyn_surf.f 
c
c Eldad Avital, Profiles data - Kaiming Ai, 2015
c
c Updated for surface waves [8], but with TUDK model
c
c Refs:
C [1] Hansen (2005) Aerodynamics of Wind Turbines, 2nd Ed.
c [2] Moriarty & Hansen (2005) Aerodyn Theory model, NREL/TP-500-36881
c [3] Liu & Janajreh (2012), Development and applictaion of an improved bem 
c model on horizonal axis wind turbines, Int J. Ener Environ Eng 3:30
c [4] Batten et al. (2008) The prediction of the hydrodynamic performnace of 
c marine current turbines, 33, 1085-1096
c [5] Lindenburg (2003) Investigation into rotor aerodynamics ECN-C--03-025
c [6] Tangler & Kocurek (2004), Wind turbine post-stall airfoil performance
c characteristics guidelines for BEM methods, NREL/CP-500-36900
c [7] Snel & Schepers (1995) Join investigation of dynamic inflow effects and
c implementation of an engineering method, ECN-C--94-107
c [8] Galloway et al(2014), Renewable energy, 297-307

       include 'incl.for'
       character*30 ifile, ifilepr(ntb)
       DIMENSION c(nxb),ar(nxb),x(nxb),thet(nxb), thet_st(nxb), psi(nmb)
       DIMENSION a(nxb,nmb),a1(nxb,nmb),fbp(nxb,nmb),alfa(nxb,nmb)
     1 ,cl(nxb,nmb),cd(nxb,nmb),c_ts(nxb,nmb),c_ps(nxb,nmb)
     2, ur(nxb,nmb),wr(nxb,nmb),tpr(nxb,nmb)
       DIMENSION time(ntb),tp(ntb),thet_tim(ntb)
     1  ,c_t(ntb),c_q(ntb),c_p(ntb),t_print(ntb)
       DIMENSION wn_qs(nxb,nmb,2),wt_qs(nxb,nmb,2),wn_int(nxb,nmb,2)
     1 ,wt_int(nxb,nmb,2),wn(ntb,nmb,2),wt(ntb,nmb,2)
       INTEGER ipro(nxb,2)

       COMMON/wave_data/h,z0,ak,omega_nor,wave_amp,omega_w

       open(unit=30,file='data_surf.in')
       read(30,*) thet0, nb, nr, rb, v0, tp0
       read(30,*) ifile, itlf,iprant,ig
       read(30,*) xs, xhub
       read(30,*) h,z0,omega_nor,wave_amp,iroot
       read(30,*) sf, timeb, nt
       read(30,*) imarch, iwake, inl
       read(30,*) nprint
       do i=1,nprint
         read(30,*) t_print(i), ifilepr(i)
       end do
       close(unit=30)

       sf = max(ONE,sf)
       imarch = max(-1,min(1,imarch))
       iroot = max(1,min(2,iroot))

       pai = 4.*atan(1.)
       h = abs(h)
       z0 = -abs(z0)
       wave_amp = abs(wave_amp)
       omega_nor = abs(omega_nor)

       write(*,*)'Blade tilted by',thet0,'(deg)'
       thet0 = thet0*pai/180.
       write(*,*)'Number of blades=',nb
       write(*,*)'Number of BEM divisions along blade=',nr
       write(*,*)'Blade length R=',rb,'(m)'
       write(*,*)'Tip speed ratio=',tp0,'without surface wave'
       write(*,*)'Blade data file is=',ifile
       write(*,*)'Using tip loss factor (Y=1,N=0)',itlf
       write(*,*)'0 - Goldstein, 1 - Prandtl=',iprant
       write(*,*)'C_T correction for high a:'
       write(*,*) '0 - Spera (1994), 1 - Glauert =',ig
       write(*,*)'BEM starts from r/R=',xs
       write(*,*)'r/R of hub =',xhub
       write(*,*) 'water and hub depths =',h,abs(z0),'(m)'
       write(*,*) 'surface wave frequency normalized by turbine freq=',
     1 omega_nor
       write(*,*) 'surface wave amplitude=',wave_amp,'(m)'
       write(*,*) 'Time step safety factor=',sf
       write(*,*) 'Time to stop=',timeb
       write(*,*) 'Number of time steps=',nt
       if(imarch.eq.-1) then 
         write(*,*) 'No unsteady aerodynamics'
       else if (imarch.eq.0) then
         write(*,*)'Unsteady aerodynamics lags behind iterative solver'
       else
         write(*,*)'Unsteady aerodynamics is fully in iterative solver'
       end if

       if(iwake.eq.0) then
         write(*,*) 'No TUDK wake model'
         imarch = -1
       else
         write(*,*) 'With TUDK wake model'
       end if

       if(inl.eq.0) then
         write(*,*) 'Linear surface wave'
       else
         write(*,*) 'Stokes surface wave'
       end if
       

       if(nt.gt.ntb) then
         write(*,*) 'Increase ntb in incl.for to nt+1'
         write(*,*) 'Stopping'
         stop
       end if

c Deciding on the trigger a>ac for C_T empirical correction
       ig = max(0,min(1,ig))
       ac = 0.2*(1-ig) + dble(ig)*0.34
c       ac = 100
       write(*,*) 'Border between momentum and empirical C_T, ac=',ac

       call read_blade(ifile)

       dx = (1.-xs)/dble(nr)
       do i=1,nr
         x(i) = xs + (i-0.5)*dx
       end do
       x(nr+1) = 2. - x(nr)
  
       call station(c,x,ar,thet_st,cd_max,ipro,nr)
 
       do ir=1,nr
         thet(ir) = thet_st(ir) + thet0
       end do
c     Init the unsteady aero model
c v0 incoming water speed    
       call init_theo(imarch)
       omega_t=tp0*v0/rb
       call surf_disp(omega_t,v0,iroot)
       print*,'Surface wave number=',ak,'(1/m)'
       print*,'Turbine and surface wave time periods=',
     1 2*pai/omega_t,2*pai/omega_w

       dpsi = 2*pai/dble(nb)

       dt = 0.13*pai*rb/(v0*sf)*min(ONE,1./omega_nor)
       write(*,*) 'dt=',dt

       do it=1,nt
         time(it) = (it-1)*dt
c         time(it) = 0.
         do ib=1,nb
           psi(ib) = (ib-1)*dpsi + omega_w*time(it)
           inc = psi(ib)/(2*pai)
           psi(ib) = psi(ib) - 2*pai*inc
         end do 

         do ib=1,nb
           do ir=1,nr
             r = x(ir)*rb
             call wave_sp(ur(ir,ib),wr(ir,ib),v0,r
     1       ,psi(ib),time(it),inl)
c             write(101,*) it,ib,ir,ur(ir,ib)
c             read*,dum
c             ur(ir,ib) = 0.1
c             wr(ir,ib) = 0.
             vv = v0+ur(ir,ib)
             vv = max(EPS,vv)
             tpr(ir,ib) = tp0*v0*vv/(vv*vv+1.d-20)
             wr(ir,ib) = wr(ir,ib)*cos(psi(ib))/v0

             a(ir,ib) = 0.
             a1(ir,ib) = 0.
             fbp(ir,ib) = 1.

             sigma = c(ir) * nb / (2*pai*x(ir))
c Linear iteration to find a and a1 (a')
             do ilin=1,3
               c_ts0=1.d20
               c_ps0=1.d20
               do iter=1,300
                 call func_fai(fai,a(ir,ib),a1(ir,ib)
     1           ,tpr(ir,ib),wr(ir,ib),x(ir))
                 call func_fbp(fbp(ir,ib),nb,fai,x(ir),xhub,itlf,iprant)

c Finding time scale of unsteady aerodynamics
                 call find_stime(stime,v0,c(ir),x(ir),a(ir,ib),a1(ir,ib)
     1           ,rb,tpr,it)

                 call cl_cd(cl(ir,ib),cd(ir,ib),alfa(ir,ib),fai
     1           ,thet(ir),ar(ir),c(ir),x(ir),tpr(ir,ib)
     2           ,cd_max,stime,dt,ipro,ir,ib,it)

                 call coeff(cn,ct,fai,cl(ir,ib),cd(ir,ib))
                 call func_a1(a1(ir,ib),fai,fbp(ir,ib),ct,sigma)
                 call func_a(a(ir,ib),fai,fbp(ir,ib),cn,sigma,ac,ig)
 
                 call ct_cp_bem(c_ts(ir,ib),c_ps(ir,ib),a(ir,ib)
     1           ,a1(ir,ib),c(ir),tpr(ir,ib),cn,ct,fai,nb)

                 tar = max(abs(c_ts(ir,ib)),1.d-12)
                 tar1 = max(abs(c_ps(ir,ib)),1.d-12)

                 if(abs(c_ts(ir,ib)-c_ts0).le.EPS*tar1
     1           .and.abs(c_ps(ir,ib)-c_ps0).le.EPS*tar) goto 20

                 c_ts0=c_ts(ir,ib)
                 c_ps0=c_ps(ir,ib)
               end do
               write(*,*) 'Linear iteration did not converge in tp='
     1         ,tpr
               if(ilin.eq.1) then
                 write(*,*) 'Using previous value as initial guess'
                 ir1 = max(1,ir-1)

                 a(ir,ib) = a(ir1,ib)
                 a1(ir,ib) = a1(ir1,ib)
               else if(ilin.eq.2) then
                 write(*,*) 'Use two previous values for initial guess'
                 ir1 = max(1,ir-1)
                 ir2 = max(1,ir-2)

                 a(ir,ib) = 2*a(ir1,ib)-a(ir2,ib)
                 a1(ir,ib) = 2*a1(ir1,ib)-a1(ir2,ib)
               else
                 write(*,*) 'Trying to using high tip speed approx'
                 if(tpe.ge.4.and.abs(a1(ir,ib)).lt.0.1) goto 10
                 write(*,*) 'Tip speed ratio is only ',tpr
                 write(*,*) 'Stopping'
                 stop
               end if
             end do
c Bi-section method, correct only for high tip speed ratio
c The effect of a1 (a') on the flow angle fai is neglected and thus
c only a has to be found.
10           a0 = -2
             a1s = 0.
             call func_fai(fai,a0,a1s,tpr(ir,ib),wr(ir,ib),x(ir))
             call func_fbp(fbp(ir,ib),nb,fai,x(ir),xhub,itlf,iprant)
c Finding time scale of unsteady aerodynamics
             call find_stime(stime,v0,c(ir),x(ir),a0,a1s
     1         ,rb,tpr(ir,ib),it)

             call cl_cd(cl0,cd0,alfa0,fai,thet(ir),ar(ir)
     1       ,c(ir),x(ir),tpr(ir,ib),cd_max,stime,dt,ipro,ir,ib,it)

             call coeff(cn,ct,fai,cl0,cd0)
             call func_af(y0,a0,fai,fbp(ir,ib),cn,sigma,ac,ig)

             a2 = 1
             call func_fai(fai,a2,a1s,tpr(ir,ib),wr(ir,ib),x(ir))
             call func_fbp(fbp(ir,ib),nb,fai,x(ir),xhub,itlf,iprant)
c Finding time scale of unsteady aerodynamics
             call find_stime(stime,v0,c(ir),x(ir),a2,a1s
     1         ,rb,tpr(ir,ib),it)

             call cl_cd(cl2,cd2,alfa2,fai,thet(ir),ar(ir)
     1       ,c(ir),x(ir),tpr(ir,ib),cd_max,stime,dt,ipro,ir,ib,it)
             call coeff(cn,ct,fai,cl2,cd2)
             call func_af(y2,a2,fai,fbp(ir,ib),cn,sigma,ac,ig)

             if(y0*y2.gt.ZERO) then
               write(*,*) 'Bi-section method failed at tp=',tpr(ir,ib)
               write(*,*) 'Try to change values of a in the code'
               write(*,*) 'Stopping'
               stop
             end if

             do iter=1,200
               a(ir,ib) = 0.5*(a0+a2)
               call func_fai(fai,a(ir,ib),a1s,tpr(ir,ib)
     1         ,wr(ir,ib),x(ir))
               call func_fbp(fbp(ir,ib),nb,fai,x(ir),xhub,itlf,iprant)
c Finding time scale of unsteady aerodynamics
               call find_stime(stime,v0,c(ir),x(ir),a(ir,ib),a1s
     1         ,rb,tpr(ir,ib),it)

               call cl_cd(cl(ir,ib),cd(ir,ib),alfa(ir,ib),fai
     1         ,thet(ir),ar(ir),c(ir),x(ir),tpr(ir,ib)
     2         ,cd_max,stime,dt,ipro,ir,ib,it)

               call coeff(cn,ct,fai,cl(ir,ib),cd(ir,ib))
               call func_af(y1,a(ir,ib),fai,fbp(ir,ib),cn,sigma,ac,ig)
               call func_a1(a1(ir,ib),fai,fbp(ir,ib),ct,sigma)

               if (abs(a2-a0).le.EPS*abs(a(ir,ib))) goto 20
               if(y0*y1.le.ZERO) then
                 a2 = a(ir,ib)
                 y2 = y1
               else
                 a0 = a(ir,ib)
                 y0 = y1
               end if
             end do
             write(*,*) 'No result in high TSR approx, tp=',tp(it)
             stop
20           call ct_cp_bem(c_ts(ir,ib),c_ps(ir,ib),a(ir,ib),a1(ir,ib)
     1         ,c(ir),tpr(ir,ib),cn,ct,fai,nb)
           end do
         end do

         if(iwake.eq.0) goto 30
c Wake model
         irat = max(0,min(1,it-1))

         do ib=1,nb
           do ir=1,nr
             wn_qs(ir,ib,2) = a(ir,ib)
             wt_qs(ir,ib,2) = a1(ir,ib)
             wn_qs(ir,ib,1) = wn_qs(ir,ib,1)*irat 
     1       + (1-irat)*wn_qs(ir,ib,2)
             wt_qs(ir,ib,1) = wt_qs(ir,ib,1)*irat 
     1       + (1-irat)*wt_qs(ir,ib,2)

             wn_int(ir,ib,1) = wn_int(ir,ib,1)*irat 
     1       + (1-irat)*wn_qs(ir,ib,2)
             wt_int(ir,ib,1) = wt_int(ir,ib,1)*irat 
     1       + (1-irat)*wt_qs(ir,ib,2)

             wn(ir,ib,1) = wn(ir,ib,1)*irat 
     1       + (1-irat)*wn_qs(ir,ib,2)
             wt(ir,ib,1) = wt(ir,ib,1)*irat 
     1       + (1-irat)*wt_qs(ir,ib,2)

             am = min(a(ir,ib),0.5d0)
             tau1 = 1.1*rb/((1.-1.3*am)*v0)
             tau2 = (0.39 - 0.26*x(ir)*x(ir))*tau1

             akw = 0.6
             hn = wn_qs(ir,ib,2)+akw*(wn_qs(ir,ib,2)-wn_qs(ir,ib,1))/dt
             ht = wt_qs(ir,ib,2)+akw*(wt_qs(ir,ib,2)-wt_qs(ir,ib,1))/dt

             wn_int(ir,ib,2) = hn+(wn_int(ir,ib,1)-hn)*exp(-dt/tau1)
             wt_int(ir,ib,2) = ht+(wt_int(ir,ib,1)-ht)*exp(-dt/tau1)

             wn(ir,ib,2)=wn_int(ir,ib,2)+(wn(ir,ib,1)-wn_int(ir,ib,2))
     1              *exp(-dt/tau2)
             wt(ir,ib,2)=wt_int(ir,ib,2)+(wt(ir,ib,1)-wt_int(ir,ib,2))
     1              *exp(-dt/tau2)
            end do
         end do

         do ib=1,nb
           do ir=1,nr
             a(ir,ib) = wn(ir,ib,2)
             a1(ir,ib) = wt(ir,ib,2)
 
             call func_fai(fai,a(ir,ib),a1(ir,ib),tpr(ir,ib)
     1      ,wr(ir,ib),x(ir))
             call func_fbp(fbp(ir,ib),nb,fai,x(ir),xhub,itlf,iprant)

c Finding time scale of unsteady aerodynamics
             call find_stime(stime,v0,c(ir),x(ir),a(ir,ib),a1(ir,ib)
     1         ,rb,tpr(ir,ib),it)

             call cl_cd(cl(ir,ib),cd(ir,ib),alfa(ir,ib),fai,thet(ir)
     1    ,ar(ir),c(ir),x(ir),tpr(ir,ib),cd_max,stime,dt,ipro,ir,ib,it)
           
c           write(42,*) ir,a(ir),a1(ir)

             call coeff(cn,ct,fai,cl(ir,ib),cd(ir,ib))
 
             call ct_cp_bem(c_ts(ir,ib),c_ps(ir,ib),a(ir,ib),a1(ir,ib)
     1     ,c(ir),tpr(ir,ib),cn,ct,fai,nb)
c Updating         
             wn(ir,ib,1) = wn(ir,ib,2)        
             wt(ir,ib,1) = wt(ir,ib,2)
         
             wn_int(ir,ib,1) = wn_int(ir,ib,2)        
             wt_int(ir,ib,1) = wt_int(ir,ib,2)
         
             wn_qs(ir,ib,1) = wn_qs(ir,ib,2)        
             wt_qs(ir,ib,1) = wt_qs(ir,ib,2)
c Updating unsteady aero
             call updat_th(alfa(ir,ib),stime,dt,ir,ib,it)
           end do
         end do
c         stop

c Calculating c_t,c_q and c_p for this tip speed ratio, tp(it)
30       x0 = xs
c         read*,dum
         c_t(it) = 0.
         c_q(it) = 0.
         c_p(it) = 0.
         do ib=1,nb
           do ir=1,nr-1
             rat = 0.5*(tpr(ir,ib)+tpr(ir+1,ib))
             rat = tp0*rat/(rat*rat+1.d-20)
             rat2 = rat*rat
             rat3 = rat2*rat

             ca = (c_ts(ir+1,ib)-c_ts(ir,ib))/(x(ir+1)-x(ir))
             cb = c_ts(ir,ib) - ca*x(ir)
             c_t(it) = (ca*(x(ir+1)*x(ir+1)-x(ir)*x(ir))/2.d0
     1       +cb*(x(ir+1)-x(ir)))*rat2 + c_t(it)

             ca = (c_ps(ir+1,ib)-c_ps(ir,ib))/(x(ir+1)-x(ir))
             cb = c_ps(ir,ib) - ca*x(ir)
             c_p(it) = (ca*(x(ir+1)**3-x(ir)**3)/3.d0
     1       +cb*(x(ir+1)**2-x(ir)**2)/2.d0)*rat3 + c_p(it)
           end do

           rat = tpr(nr,ib)
           rat = tp0*rat/(rat*rat+1.d-20)
           rat2 = rat*rat
           rat3 = rat2*rat

           ca = (c_ts(nr,ib)*(1-itlf)-c_ts(nr,ib))/(1-x(nr))
           cb = c_ts(nr,ib) - ca*x(nr)
           c_t(it) = (ca*(1-x(nr)*x(nr))/2.d0
     1     +cb*(1-x(nr)))*rat2 + c_t(it)

           ca = (c_ps(nr,ib)*(1-itlf)-c_ps(nr,ib))/(1-x(nr))
           cb = c_ps(nr,ib) - ca*x(nr)
           c_p(it) = (ca*(1-x(nr)**3)/3.d0
     1   +cb*(1-x(nr)**2)/2.d0)*rat3 + c_p(it)

           c_q(it) = c_p(it)/tp0

c Correcting to c_ts=dC_T/dx and c_ps=dC_P/dx
           do ir=1,nr
             c_ts(ir,ib)=c_ts(ir,ib)
             c_ps(ir,ib)=x(ir)*c_ps(ir,ib)
           end do
         end do
         do ii=1,nprint
           if(abs(time(it)-t_print(ii)).lt.dt)
     1      call print_int(time(it),tp0,thet(1),x,alfa,cl,cd
     2      ,a,a1,c_ts,c_ps,nr,it,ifilepr(ii))
         end do
         if(time(it).ge.timeb) goto 50
       end do
       write(*,*) 'Reached only time',time(nt) 
       it = it - 1      

50     nt = it
       open(unit=30,file='coeff_time.plt')
       write(30,*)'# time,C_T, C_Q, C_P'
       do it=1,nt
         write(30,1000) time(it),c_t(it),c_q(it),c_p(it)
       end do

1000   FORMAT(4(1x,1pe14.6))

       stop
       END

C Reads the blade data
       SUBROUTINE read_blade(ifile)
       include 'incl.for'
       character*30 ifile, ifile_pr, idum
         
       DIMENSION c_sta(nxb), x_sta(nxb), thet_sta(nxb)
       COMMON/sta/c_sta, x_sta, thet_sta, nx_sta

       pai = 4.*atan(1.)
       open(unit=40,file=ifile)
       read(40,*) idum
       read(40,*) nx_sta
       do i=1,nx_sta
         read(40,*) x_sta(i),c_sta(i),thet_sta(i)
     1    ,ifile_pr,itype
         thet_sta(i) = thet_sta(i) * pai/180
         if(itype.eq.1) then
           call read_prof(ifile_pr,i)
         else
           call read_prof_cl_cd(ifile_pr,i)
         end if
       end do

       close(unit=40)

       return
       END

C Reads the profile data; alfa, cl, cd (type 1) 
       SUBROUTINE read_prof(ifile,k)
       include 'incl.for'
       character*30 ifile, idum

       DIMENSION cl_p(nxb,npb),cd_p(nxb,npb),alfa_p(nxb,npb)
     1  ,alfa0_p(npb),alfass_p(npb),alfase_p(npb),cd0_p(npb)
       INTEGER nx(npb)
       COMMON/profile/cl_p,cd_p,alfa_p,alfa0_p
     1  ,alfass_p,alfase_p,cd0_p,nx,np 

       pai = 4*atan(1.)

       open(unit=30,file=ifile)
       read(30,*) idum
       read(30,*) nx(k)
       do j=1,nx(k)
         read(30,*) alfa_p(j,k), cl_p(j,k), cd_p(j,k)
         alfa_p(j,k) = alfa_p(j,k)*pai/180.
       end do 
   
       if(alfa_p(nx(k),k).le.alfa_p(1,k)) then
         write(*,*) 'Rearrange alfa in ascending order'
         write(*,*) 'In file ', ifile
         stop
       end if


       close(unit=30)

C Finding alfa of cl=0 and cd at that angle
       sma1 = 1.d20
       do i=1,nx(k)
         if(abs(cl_p(i,k)).le.sma1) then
           is1 = i
           sma1 = abs(cl_p(i,k))
         end if
       end do

       sma2 = 1.d20
       do i=1,nx(k)
         if(abs(cl_p(i,k)).le.sma2.and.i.ne.is1) then
           is2 = i
           sma2 = abs(cl_p(i,k))
         end if
       end do

       ca = (cl_p(is2,k)-cl_p(is1,k))/(alfa_p(is2,k)-alfa_p(is1,k))
       alfa0_p(k) = alfa_p(is1,k) - cl_p(is1,k)/ca  

       cd0_p(k) = cd_p(is1,k)+(alfa0_p(k)-alfa_p(is1,k))
     1 *(cd_p(is2,k)-cd_p(is1,k))/(alfa_p(is2,k)-alfa_p(is1,k)) 

c Finding point of where to start the dynamic stall (lift curve slope smaller than 0.05 1/deg)

       alfass_p(k) = alfa_p(1,k)
       do i=2,is1
         ca = (cl_p(i,k)-cl_p(i-1,k))/(alfa_p(i,k)-alfa_p(i-1,k))
         if(ca.gt.0.05*180/pai) then
           alfass_p(k) = alfa_p(i-1,k)
           goto 10
         end if
       end do

10     alfase_p(k) = alfa_p(nx(k),k)
       do i=is1+1,nx(k)
         ca = (cl_p(i,k)-cl_p(i-1,k))/(alfa_p(i,k)-alfa_p(i-1,k))
         if(ca.lt.0.05*180/pai) then
           alfase_p(k) = alfa_p(i-1,k)
           goto 20
         end if
       end do
20     continue     

       return
       END 

C The following subroutine has not been tested! 
C Use only type 1 data for profile data.
C Reads the profile data, type 2
       SUBROUTINE read_prof_cl_cd(ifile,k)
       include 'incl.for'
       character*30 ifile

       DIMENSION cd(nxb), cl(nxb)
       DIMENSION cl_p(nxb,npb),cd_p(nxb,npb),alfa_p(nxb,npb)
     1  ,alfa0_p(npb),alfass_p(npb),alfase_p(npb),cd0_p(npb)
       INTEGER nx(npb)
       COMMON/profile/cl_p,cd_p,alfa_p,alfa0_p
     1  ,alfass_p,alfase_p,cd0_p,nx,np 

       pai = 4*atan(1.)

       open(unit=30,file=ifile)
       read(30,*) cla, alfa0_p(k), alfa_st
       alfa0_p(k) = alfa0_p(k) *pai/180.
       alfa_st = alfa_st *pai/180.
       cla = cla * 180/pai
       read(30,*) nx(k)
       do j=1,nx(k)
         read(30,*) cl_p(j,k), cd_p(j,k)

c Finding alfa (assuming linear variation of cl with alfa)
         alfa_p(j,k) = cl_p(j,k)/cla+alfa0
       end do

       close(unit=30)

       np = k

       return
       END 

C Calculates the BEM stations chord and profile
       SUBROUTINE station(c,x,ar,thet,cd_max,ipro,nr)
       include 'incl.for'
       DIMENSION c(nxb),ar(nxb),x(nxb),thet(nxb)
       INTEGER ipro(nxb,2)
       
       DIMENSION c_sta(nxb), x_sta(nxb), thet_sta(nxb)
       COMMON/sta/c_sta, x_sta, thet_sta, nx_sta

       s = 0.
       do i=1,nx_sta-1
         s = s + 0.5*(x_sta(i+1)-x_sta(i))*(c_sta(i)+c_sta(i+1))
       end do

       aspect = min(50.d0,1./s)
c Hoerner's estimate for Cd of a flate plate with angle of attack of 90 degs.
       cd_max = 1.11 + 0.018*aspect 

       write(*,*) 'cd_max=',cd_max

       do j=1,nr
         do i=1,nx_sta
           if(x(j).lt.x_sta(i)) goto 10
         end do
         ii = nx_sta
10       ii = min(nx_sta-1,max(i-1,1))
         ipro(j,1) = ii
         ipro(j,2) = ii+1

         ar(j) = (x(j)-x_sta(ii))/(x_sta(ii+1)-x_sta(ii))
         c(j) = c_sta(ii) + ar(j)*(c_sta(ii+1)-c_sta(ii))
         thet(j) = thet_sta(ii) + ar(j)*(thet_sta(ii+1)-thet_sta(ii))
       end do
       return
       END

c Calculate a
       SUBROUTINE func_a(a,fai,fbp,cn,sigma,ac,ig)
       include 'incl.for'
 
       si = sin(fai)
       si2 = si*si
       rat = sigma*cn
       akw = 4*fbp*si2*rat/(rat*rat+1.d-20)
       rat = 4*fbp*si2 + sigma*cn
       an = sigma*cn*rat/(rat*rat+1.d-20)
       
       if(an.le.ac) then
c Momentum method result
         a = an
       else
         if(ig.eq.0) then
c Spera's emperical straight line 
           rat = akw * (1 - 2 * ac) + 2
           det = max(ZERO,rat*rat + 4*(akw*ac*ac-1)) 
           a = 0.5*(2+akw*(1-2*ac)-sqrt(det))
         else
c Glauert emperical line 
           if(si2.lt.10*sigma*cn.and.si2.lt.1d-20) then
             a = 1
             return
           else
c Bi-section solution
             a0 = ac
             call func_af(y0,a0,fai,fbp,cn,sigma,ac,ig)
             a2 = 1
             call func_af(y2,a2,fai,fbp,cn,sigma,ac,ig)
             if(y0*y2.gt.ZERO) then
c               print*,'enlarge in func_a'
c               read*,dum
               a = an
               return
             else
               do k=1,500
                 a = 0.5*(a0+a2)
                 call func_af(y1,a,fai,fbp,cn,sigma,ac,ig)
                 if(abs(a2-a0).lt.EPS*abs(a)) return
                 if(y0*y1.le.ZERO) then
                   a2 = a
                   y2 = y1
                 else
                   a0 = a
                   y0 = y1
                 end if
               end do
             end if
           end if
         end if
       end if
       return
       END

C Match C-T for the bi-section method and high tip speed ratio approx
       SUBROUTINE func_af(y,a,fai,fbp,cn,sigma,ac,ig)
       include 'incl.for'
 
       si = sin(fai)
       si2 = si * si

       if(a.le.ac) then
         y = 4*a*(1-a)*fbp*si2 -(1-a)*(1-a)*sigma*cn
       else
         ys = 4*(ac*ac+(1-2*ac)*a)*fbp*si2 - (1-a)*(1-a)*sigma*cn
         yg = 4*a*(1-0.25*(5-a/ac)*a)*fbp*si2 - (1-a)*(1-a)*sigma*cn
         y = yg*ig + ys*(1-ig)
       end if

       return
       END

C Calculates c_t and c_p 
       SUBROUTINE ct_cp(c_ts,c_ps,a,a1,fbp,tp,ac,ig)
       include 'incl.for'

       c_ps =8*tp*tp*(1-a)*a1*fbp

       if(a.le.ac) then
         c_ts = 4*a*(1-a)*fbp
       else
         if(ig.eq.0) then
           c_ts = 4*(ac*ac+(1-2*ac)*a)*fbp
         else
           c_ts = 4*a*(1-0.25*(5-a/ac)*a)*fbp
         end if
       end if
       
       c_ts = 2 * c_ts
 
       return
       END

       SUBROUTINE ct_cp_bem(c_ts,c_ps,a,a1,c,tp,cn,ct,fai,nb)
       include 'incl.for'
        
       pai = 4.*atan(1.)
       si = sin(fai)
       si2 = si * si

       cs = cos(fai)
       c_ts = 1.d0/pai * (1-a)*(1-a)*c*cn/si2

       c_ps = 1.d0/pai*tp*(1-a)*(1-a)*c*ct/si2

       return
       END


c Calculate a1 (a')
       SUBROUTINE func_a1(a1,fai,fbp,ct,sigma)
       include 'incl.for'

       rat = 2*fbp*sin(2*fai)-sigma*ct
       a1 = sigma*ct * rat/(rat*rat+1.d-20)
c       print*,'a1,fai,fbp,ct,sigma=',a1,fai,fbp,ct,sigma

       return
       END

c Calculate fbp, Prandtl tip correction
       SUBROUTINE func_fbp(fbp,nb,fai,r,rhub,itlf,iprant)
       include 'incl.for'

       pai = 4.*atan(1.)
       si = sin(fai)
       cs = cos(fai)
       fbp = 1.

       if(iprant.eq.1) then
c Prandtl's tip loss factor
         if(abs(si).gt.1.d-8) then
           f = 0.5*nb*(1./r-1)/abs(si)
           f = min(20.d0,max(ZERO,f))
           ft = dexp(-f)
           ft = 2./pai*acos(ft)

           f = 0.5*nb*(1.-rhub/r)/abs(si)
           f = min(20.d0,max(ZERO,f))
           fh = dexp(-f)
           fh = 2./pai*acos(fh)
           fbp = fh * ft
         else
           fbp = 1.d0
         end if
       else
c Goldstein's tip loss factor
         if(abs(si).gt.1.d-8) then
           f = 0.5*nb*cs/(si*r)-0.5
           if(f.gt.0.) then
             fe2 = max(-10.d0,-2*f)
             fer2 = max(-10.d0,-2*r*f)
             fer = max(-10.d0,(r-1)*f)
           else
             fe2 = max(-10.d0,2*f)
             fer2 = max(-10.d0,2*r*f)
             fer = max(-10.d0,-(r-1)*f)
           end if
           rat = dexp(fer)*(1+dexp(fer2))/(1+dexp(fe2))
           rat = max(ZERO,min(ONE,rat))
           fpb = 2./pai*acos(rat)
         else
           fpb = 1.
         end if
       end if

       fbp = fbp*itlf + (1-itlf)
      
       fbp = max(ZERO,min(ONE,fbp))

       return
       END


c Calculate cl and cd
       SUBROUTINE cl_cd(cl,cd,alfa_l,fai,thet,ar,c,x,tp
     1 ,cd_max,stime,dt,ipro,jj,ib,it)
       include 'incl.for'
       DIMENSION clm(2),cdm(2)
       INTEGER ipro(nxb,2)

       DIMENSION cl_p(nxb,npb),cd_p(nxb,npb),alfa_p(nxb,npb)
     1  ,alfa0_p(npb),alfass_p(npb),alfase_p(npb),cd0_p(npb)
       INTEGER nx(npb)
       COMMON/profile/cl_p,cd_p,alfa_p,alfa0_p
     1  ,alfass_p,alfase_p,cd0_p,nx,np 

       COMMON/stall_pro/alfa_e,alfa0,alfa_ss,alfa_se,cd0,cl_s
 
       pai = 4.*atan(1.)      
       alfa_l = fai - thet

       call find_alfai(alfai,alfa_l,stime,dt,jj,ib,it)
c Adding unsteady wake attached flow contribution
       alfa_l = alfa_l - alfai
       alfa_l =max(-0.5*pai,min(0.5*pai,alfa_l))
       tpx = tp*x

       do k=1,2
         if(alfa_l.lt.alfa_p(1,ipro(jj,k))) then
c post stall conditions, Viterna & Corrigan (1981), typo in cb2 corrected
c Snel's 3D rotational correction for cl
           cl_pot = 2*pai*(alfa_p(1,ipro(jj,k))-alfa0_p(ipro(jj,k)))
           cl_s = cl_p(1,ipro(jj,k))
           cl_s = cl_s + 3.1*tpx*tpx*c*c/(x*x)
     1     /(1+tpx*tpx)*(cl_pot-cl_s)

           cd_s = cd_p(1,ipro(jj,k))
           si_s = sin(alfa_p(1,ipro(jj,k)))
           cs_s = cos(alfa_p(1,ipro(jj,k)))

           ca2 = (cl_s-cd_max*si_s*cs_s)*si_s/(cs_s*cs_s)
           cb2 = (cd_s - cd_max*si_s*si_s)/cs_s
           cb1 = cd_max
           ca1 = 0.5 * cb1

           si = sin(alfa_l)
           si2 = si * si
           cs = cos(alfa_l)
           cs2 = cs * cs

           clm(k) = ca1*sin(2*alfa_l) + ca2*cs2/si
           cdm(k) = cb1*si*si + cb2*cs
         else if (alfa_l.gt.alfa_p(nx(ipro(jj,k)),ipro(jj,k))) then
c post stall conditions, Viterna & Corrigan (1981), typo in cb2 corrected
c Snel et al. rotational correction
           cl_pot = 2*pai*(alfa_p(nx(ipro(jj,k)),ipro(jj,k))
     1     -alfa0_p(ipro(jj,k)))
           cl_s = cl_p(nx(ipro(jj,k)),ipro(jj,k))
           cl_s = cl_s + 3.1*tpx*tpx*c*c/(x*x)
     1     /(1+tpx*tpx)*(cl_pot-cl_s)

           cd_s = cd_p(nx(ipro(jj,k)),ipro(jj,k))
           si_s = sin(alfa_p(nx(ipro(jj,k)),ipro(jj,k)))
           cs_s = cos(alfa_p(nx(ipro(jj,k)),ipro(jj,k)))

           ca2 = (cl_s-cd_max*si_s*cs_s)*si_s/(cs_s*cs_s)
           cb2 = (cd_s - cd_max*si_s*si_s)/cs_s
           cb1 = cd_max
           ca1 = 0.5 * cb1

           si = sin(alfa_l)
           si2 = si * si
           cs = cos(alfa_l)
           cs2 = cs * cs

           clm(k) = ca1*sin(2*alfa_l) + ca2*cs2/si
           cdm(k) = cb1*si*si + cb2*cs
         else
           ii = 1
           do i=1,nx(ipro(jj,k))
             if(alfa_l.lt.alfa_p(i,ipro(jj,k))) goto 10
           end do
           i = nx(ipro(jj,k))
10         ii = min(nx(ipro(jj,k))-1,max(i-1,1))
c Snel et al. rotational correction
           cl_pot = 2*pai*(alfa_p(ii,ipro(jj,k))
     1     -alfa0_p(ipro(jj,k)))
           cl_0 = cl_p(ii,ipro(jj,k))
           cl_0 = cl_0 + 3.1*tpx*tpx*c*c/(x*x)
     1     /(1+tpx*tpx)*(cl_pot-cl_0)

           cl_pot = 2*pai*(alfa_p(ii+1,ipro(jj,k))
     1     -alfa0_p(ipro(jj,k)))
           cl_1 = cl_p(ii+1,ipro(jj,k))
           cl_1 = cl_1 + 3.1*tpx*tpx*c*c/(x*x)
     1     /(1+tpx*tpx)*(cl_pot-cl_1)

           clm(k) = cl_0 + (alfa_l-alfa_p(ii,ipro(jj,k)))
     1     *(cl_1-cl_0)/(alfa_p(ii+1,ipro(jj,k))-alfa_p(ii,ipro(jj,k)))

           cdm(k) = cd_p(ii,ipro(jj,k))
     1     + (alfa_l-alfa_p(ii,ipro(jj,k)))
     2     *(cd_p(ii+1,ipro(jj,k))-cd_p(ii,ipro(jj,k)))
     3     /(alfa_p(ii+1,ipro(jj,k))-alfa_p(ii,ipro(jj,k)))
         end if
       end do

       cd = cdm(1) + ar * (cdm(2)-cdm(1))
       cl = clm(1) + ar * (clm(2)-clm(1))
c Going back to the geometric alfa
       alfa_l = alfa_l + alfai
c Preparing the required data for the dynamic stall calculations
       cl_s = cl
       alfa0 = alfa0_p(ipro(jj,1))
     1 +ar*(alfa0_p(ipro(jj,2))-alfa0_p(ipro(jj,1)))
       alfa_ss = alfass_p(ipro(jj,1))
     1 +ar*(alfass_p(ipro(jj,2))-alfass_p(ipro(jj,1)))
       alfa_se = alfase_p(ipro(jj,1))
     1 +ar*(alfase_p(ipro(jj,2))-alfase_p(ipro(jj,1)))
       cd0 = cd0_p(ipro(jj,1))
     1 +ar*(cd0_p(ipro(jj,2))-cd0_p(ipro(jj,1)))
c Adding dynamic stall
       call dyn_cl(cl,cd,alfa_l,stime,dt,jj,ib,it)
c Adding implusive lift and induced drag
       call impluse(cl,cd,alfa_l,alfai,stime,dt,jj,ib,it)
       
 
       return
       END

c Calculate cn and ct
       SUBROUTINE coeff(cn,ct,fai,cl,cd)
       include 'incl.for'

       cs = cos(fai)
       si = sin(fai)

       cn = cl*cs + cd*si
       ct = cl*si - cd*cs

       return
       END

c Calculate fai
       SUBROUTINE func_fai(fai,a,a1,tp,wr,x)
       include 'incl.for'

       pai = 4. * atan(1.)
       rat = (1 + a1) * tp * x + wr
       if(abs(rat).gt.1.d-20) then
         fai = (1-a)/rat
         fai = min(1.d3,max(-1.d3,fai))
         fai = atan(fai)
       else
         if(abs(1-a).gt.1.d-19) then
           if((1-a)*rat.gt.ZERO) then
             fai = 0.5*pai
           else
             fai = -0.5*pai
           end if
           fai = 0.
         end if
       end if

       return
       END

       SUBROUTINE print_int(time,tp,thet,x,alfa,cl,cd
     1 ,a,a1,c_ts,c_ps,nr,it,ifile)
       include 'incl.for'
       character*30 ifile
       DIMENSION a(nxb,nmb),a1(nxb,nmb),x(nxb),alfa(nxb,nmb)
     1 ,cl(nxb,nmb),cd(nxb,nmb),c_ts(nxb,nmb),c_ps(nxb,nmb)

       pai = 4.*atan(1.)
       open(unit=31,file=ifile)
         write(31,*)'#time=',time,'tip speed ratio=',tp
         write(31,*)'#hub pitch angle=',thet
     1   ,'(rad), x,alfa,cl,cd,a,a1,d(c_t)/dx,d(c_p)/dx'
         do i=1,nr
           write(31,1000) real(x(i)),real(alfa(i,1)*180./pai)
     1    ,real(cl(i,1)),real(cd(i,1)),real(a(i,1)),real(a1(i,1))
     2    ,real(c_ts(i,1)),real(c_ps(i,1))
         end do
       close(unit=31)

1000   format(8(1x,1pe12.4))
       return
       END

       SUBROUTINE surf_disp(omega_t,v0,iroot)
       include 'incl.for'
       DIMENSION ak_root(2)
       COMMON/wave_data/h,z0,ak,omega_nor,wave_amp,omega_w

       omega_w = omega_t * omega_nor
       print*,'Wave frequency (rad/s)=',omega_w
       pai = 4.*atan(1.)
       g = 9.81
       ak0 = 0.
       call disper(y0,omega_w,ak0,h,v0)
       do ii=1,2
         dk = sqrt(omega_w)/(g*100)
         do i=1,100000
           ak = ak0+dk
           call disper(y,omega_w,ak,h,v0)
           if(y*y0.lt.ZERO) then
             if(abs(dk).le.abs(ak)*EPS) goto 10
             dk = 0.5*dk
           else
             ak0 = ak
             y0 = y
           end if
         end do
         print*,'Could not find the length of the surface wave'
         print*,'Change omega_nor'
         stop

10       ak0 = ak
         ak_root(ii) = ak
         y0 = y
       end do

       write(*,*) 'Possible wave numbers are',ak_root(1),ak_root(2)
       if(iroot.eq.1) then
         write(*,*) 'Choosing long wave (low wave number)'
         ak = ak_root(1)
       else
         write(*,*) 'Choosing short wave (high wave number)'
         ak = ak_root(2)
       end if
       return
       END

       SUBROUTINE disper(y,omega,ak,h,v0)
       include 'incl.for'

       dum = omega - ak*v0
c       dum = omega
       g = 9.81
       at = min(99.d0,ak*h)
       at = exp(-at)
       at = (1-at)/(1+at)
       y = dum*dum-g*ak*at

       return
       END

       SUBROUTINE wave_sp_lin(ur,vr,v0,r,psi,time)
       include 'incl.for'
       COMMON/wave_data/h,z0,ak,omega_nor,wave_amp,omega_w

       z = z0+r*cos(psi)
       g = 9.81
       dum = omega_w-ak*v0
       rat = dum/(dum*dum+1.d-20)
       amp = ak*g*wave_amp*rat
       akh2 = min(99.d0,2*ak*h)
       akhz2 = min(99.d0,2*ak*(h+z))

       ur = amp*exp(ak*z)*(1+exp(-akhz2))/(1+exp(-akh2))
     1    *cos(omega_w*time)
       vr = -amp*exp(ak*z)*(1-exp(-akhz2))/(1+exp(-akh2))
     1    *sin(omega_w*time)


       return
       END

       SUBROUTINE wave_sp(ur,wr,v0,r,psi,time,inl)
       include 'incl.for'
       COMMON/wave_data/h,z0,ak,omega_nor,wave_amp,omega_w

       z = z0+r*cos(psi)
       dum = omega_w-ak*v0

       akh2 = min(10.d0,2*ak*h)
       akhz2 = max(ZERO,min(10.d0,2*ak*(z+h)))
       akhz4 = max(ZERO,min(10.d0,4*ak*(z+h)))
       akz = max(-10.d0,min(ZERO,ak*z))
       ak2zh = max(-10.d0,min(ZERO,ak*(2*z-h)))

       theta = -omega_w*time
       st = sin(theta)
       st2 = sin(2*theta)
       ct = cos(theta)
       ct2 = cos(2*theta)

       aa = exp(akz)*(1+exp(-akhz2))/(1-exp(-akh2))
       bb = 3*ak*wave_amp/8.*exp(ak2zh)
     1    *(1+exp(-akhz4))/(1-exp(-akh2))**3
       ur = wave_amp*dum*aa*(ct+2*inl*bb*ct2)

       aav = exp(akz)*(1-exp(-akhz2))/(1-exp(-akh2))
       bbv = 6*ak*wave_amp/8.*exp(ak*(2*z-h))
     1    *(1-exp(-akhz4))/(1-exp(-akh2))**3
       wr = wave_amp*dum*aav*(st+inl*bb*st2)
     1    + wave_amp*dum*aa*inl*bbv*st2

c       ur = amp
c       write(102,*) ur,wr

       return
       END


