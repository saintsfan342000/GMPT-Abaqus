c     user amplitude subroutine
      Subroutine uamp(
C          passed in for information and state variables
     *     ampName, time, ampValueOld,  dt,  nProps, props, nSvars, svars, lFlagsInfo,
     *     nSensor, sensorValues, sensorNames, jSensorLookUpTable, 
C          to be defined
     *     ampValueNew, 
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral,
     *     AmpIncDoubleIntegral)
      
      include 'aba_param.inc'

C     svars - additional state variables, similar to (V)UEL
      dimension sensorValues(nSensor), svars(nSvars), props(nProps)
      character*80 sensorNames(nSensor)
      character*80 ampName



C     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
C     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           nFlagsInfo        = 2)
C     optional flags to be defined
      parameter (iComputeDeriv       = 1,
     *           iComputeSecDeriv    = 2,
     *           iComputeInteg       = 3,
     *           iComputeDoubleInteg = 4,
     *           iStopAnalysis       = 5,
     *           iConcludeStep       = 6,
     *           nFlagsDefine        = 6)
      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)

C     User code to compute  ampValue = F(sensors)

         if (lFlagsInfo(iInitialization).eq.1) then 
            ampValueNew      = zero
C            write (*,*) 'initialization'
         else
            if (lFlagsInfo(iRegularInc).eq.1) then 
c     get sensor values first; note that
c               write (*,*) meanradius
               iR_Torque  = iGetSensorID('PRESS', jSensorLookUpTable)
               valueR_Torque  = sensorValues(iR_Torque)               
              ampValueNew = valueR_Torque
C               write (*,*) valueR_Torque,ampValueNew               
            else
               ampValueNew = ampValueOld
C            write (*,*) 'bad increment'
            end if
         end if

      return
      end
C
C
C
C
C
C
C
C
C
C
C
C
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA) 
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C
C The subroutine gives the material orientation direction as output
C First Material direction components in the global basis directions
      UVAR(1) = DIRECT(1,1)
      UVAR(2) = DIRECT(2,1)
      UVAR(3) = DIRECT(3,1)
C Second Material direction components in the global basis directions
      UVAR(4) = DIRECT(1,2)
      UVAR(5) = DIRECT(2,2)
      UVAR(6) = DIRECT(3,2)
C Third Material direction components in the global basis directions
      UVAR(7) = DIRECT(1,3)
      UVAR(8) = DIRECT(2,3)
      UVAR(9) = DIRECT(3,3)
      RETURN
      END
C
C
C
C
C
C
C
C
C
C
C
C
C
        SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1  RPL,DDSDDT,DRPLDE,DRPLDT,
     1  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     1  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     1  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        CHARACTER*80 CMNAME
C
        DIMENSION DROT(3,3),DSTRAN(NTENS),STRESS(NTENS),STRAN(NTENS),
     1  DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),COORDS(3),
     1  STATEV(NSTATV),PROPS(NPROPS),PREDEF(1),DPRED(1),TIME(2),
     1  DFGRD0(3,3),DFGRD1(3,3)
C 
       dimension esp(6),cmatl(6,6),cfelas(6,6),sn(6),a_sn(6),dep(6)
c cfelas is matrix for full 3D calculations, neglecting two out-of-plane shears
       common/ki_select/i_hard,i_yield
       common/ksh_mtl/ams1,ams2,ams3,y_0,abac1,abac2,abac12,delk
       common/elas/ ekp,esh
       common/khill/r11,r22,r33,r12,r23,r31,a1,a2,a3,a4,a5,a6
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
       common /yld2_2 / oc11,oc22,oc12,oc21,oc33
       common /yld2_3 / tc11,tc22,tc12,tc21,tc33
	 common /rate/ eratlim,erexp
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c      *** beginning of user input in "*.inp" file ***
c      -------------------------------------
c      USER INPUT IN "*.inp" FILE : PROPS(*)
c      -------------------------------------
c      PROPS(1) : Young's modulus
c      PROPS(2) : Poisson's ratio
c
c      PROPS(3) : Hardening Option 
c                 |1| = Power law hardening
c                 |2| = Voce Hardening
c                 |3| = Table
c	if PROPS(3) < 0 include rate sensitivity
c
c      PROPS(4) : Yield Function Option
c                 1 = Hill's(1948) 
c                 2 not used
c                 3 = Yld2004 Old Version1
c                 4 = Yld2004 Old Version2
c                 5 = Yld2004 Anisotropy-Final version (Please use Option 5)
c
c      PROPS(5) : delk(kinemtaic hardening portion) 
c                  delk =0 : pure kinematic hardening with Chaboche Model
c                  delk =1 & Chaboche coefficients >0  isotropic+kin hardening
c                  delk=1 & Chaboche coefficients =0: isotropic hardening 
c
c	PROPS(6) : # of interpolation points through thickness for properties
c		   PROPS(6) = nint
c
c	PROPS(7) : field variable (temp) used for interpolation associated
c			with this material set of constants
c      
c      PROPS(8) : Initial Yield Stress
c                 PROPS(8) = Sigmabar at Epbar = 0 in Stress-Strain Curve
c
c      PROPS(9-11) : For power-law hardening, Sigmabar = a(b + Epbar)**c
c                   For voce hardening, Sigmabar = a - bEXP(-c*Epbar)
c                   PROPS(9) = a
c                   PROPS(10) = b
c                   PROPS(11) = c
c
c      PROPS(12-17) : Anisotropic coefficients for Hill's(1948) yield function
c                    PROPS(12) = r11
c                        thru
c                    PROPS(17) = r31
c
c      PROPS(12-30) : Anisotropic coefficients for Yld2004 function
c                     PROPS(12) = alpa1
c                     PROPS(13) = alpa2
c                     PROPS(14) = alpa3
c                     PROPS(15) = alpa4
c                     PROPS(16) = alpa5
c                     PROPS(17) = alpa6
c                     PROPS(18) = alpa7
c                     PROPS(19) = alpa8
c                     PROPS(20) = alpa9
c                     PROPS(21) = alpa10
c                     PROPS(22) = alpa11
c                     PROPS(23) = alpa12
c                     PROPS(24-29) = alpa13-18
c                     PROPS(30) = m
c
c	PREOPS(31-32) Chaboche hardening constants C (abac1) and gamma (abac2)
c		      PROPs(31) = abac1
c		      PROPS(32) = abac2
c	PROPS(33-34) use for rate sensitivity,threshold limit & exponent
c		      PROPS(33) = eratlim
c		      PROPS(34) = erexp
c
c	Repeat next 28 PROPS for next set of material constants. 
c	Do this nint times
c	Interpolate using temp and PROPS(7)
c
c      -----------------------------------
c      USER INPUT IN "*.inp" FILE : NSTATV
c      -----------------------------------
c      NSTATV : Number of State Variables
c               NSTATV = 8
c
c      *** end of user input in "*.inp" file ***
c
c
c      *** Further Definition of ABAQUS variables ***
c
c      *P.S) Below part is not user input 
c            (It is just explanation for ABAQUS variables)
c      --------------------------------------- 
c      DEFINITION OF STATE VARIABLE : STATEV(*)
c      ----------------------------------------
c      STATEV(1) : equivalent plastic strain  
c      STATEV(2) : equivalent stress
c      STATEV(3) : back-stress component 11  
c      STATEV(4) : back-stress component 22  
c      STATEV(5) : back-stress component 33
c      STATEV(6) : back-stress component 12
c      STATEV(7) : back-stress component 23  
c      STATEV(8) : back-stress component 31
c
c      -------------------------------------------- 
c      DEFINITION OF INCREMENTAL STRAIN : DSTRAN(*)
c      --------------------------------------------
c      DSTRAN(1-6) : INCREMENTAL STRAIN FROM reference state ((N)-step) 
c                    TO current state ((N+1)-step and (I)-iteration)
c                    DSTRAN(1) = Epsilon(xx)               
c                    DSTRAN(2) = Epsilon(yy)               
c                    DSTRAN(3) = Epsilon(zz)               
c                    DSTRAN(4) = 2*Epsilon(xy)               
c                    DSTRAN(5) = 2*Epsilon(yz)               
c                    DSTRAN(6) = 2*Epsilon(zx)               
c      *P.S) DSTRAN(4 thru 6) = 2*(Tensor Strain)              
c
c      -------------------------------- 
c      DEFINITION OF STRESS : STRESS(*)
c      --------------------------------
c      STRESS(1-6) : BEFORE Umat - Value of reference state (N-step)
c                    AFTER Umat  - Value of current state 
c                                  ((N+1)-step and (I)-iteration)
c                    STRESS(1) = Sigma(xx)           
c                    STRESS(2) = Sigma(yy)           
c                    STRESS(3) = Sigma(zz)           
c                    STRESS(4) = Sigma(xy)           
c                    STRESS(5) = Sigma(yz)           
c                    STRESS(6) = Sigma(zx)           
c
c      ------------------------------------------------------- 
c      DEFINITION OF JACOBIAN FOR ABAQUS IMPLICIT: DDSDDE(*,*)
c      -------------------------------------------------------
c      DDSDDE(6,6) : BEFORE Umat - Value of previous iteration
c                    AFTER Umat - Value of current state
c------------------------------------------------------------------------------
c      assign of ABAQUS variables
c      --------------------------
c	# of properties associated with each field variable
       data nlprp /28/
       data zero,one,two,three,smtol /0.d0,1.d0,2.d0,3.d0,1.d-8/
       iprint = 0
c initialize incremental plastic strain
	do i = 1,6
	dep(i) = zero
	enddo
       young = PROPS(1) 
       xnu = PROPS(2) 
       ekp = young/(3.0*(1.0-2.0*xnu))
       esh = young/(2.0*(1.0+xnu))
c
       i_hard = NINT(PROPS(3))
       i_yield = NINT(PROPS(4))
c
       nint1 = PROPS(6)
c	index for first property value
       n1 = 7
c interpolate to find properties
       nist = 1
       fac = zero
       do 20 i = 1,nint1
       ind = n1 + (i-1)*nlprp
       if(TEMP .gt. PROPS(ind)) then
       nist = nist + 1
c TEMP is greater than highest supplied field variable
	if(nist .gt. nint1) then
        nist = nint1
        indm1 = ind - nlprp
        go to 25
        endif
       elseif(i .eq. 1) then
c TEMP is less than or equal to lowest supplied field variable
       indm1 = ind
       go to 25
       else
       indm1 = ind - nlprp
       fac = (PROPS(ind)-TEMP)/(PROPS(ind)-PROPS(indm1))
       if(fac .lt. zero) fac = zero
       go to 25
       endif
20	continue
25	continue
c
c	compute properties using interpolation
c
      y_0 = PROPS(ind+1) - fac*(PROPS(ind+1)-PROPS(indm1+1))
      ams1 = PROPS(ind+2) - fac*(PROPS(ind+2)-PROPS(indm1+2))
      ams2 = PROPS(ind+3) - fac*(PROPS(ind+3)-PROPS(indm1+3))
	ams3 = PROPS(ind+4) - fac*(PROPS(ind+4)-PROPS(indm1+4))
c
	if(i_yield .eq. 1) then
      	r11 = PROPS(ind+5) - fac*(PROPS(ind+5)-PROPS(indm1+5))
		r22 = PROPS(ind+6) - fac*(PROPS(ind+6)-PROPS(indm1+6))
		r33 = PROPS(ind+7) - fac*(PROPS(ind+7)-PROPS(indm1+7))
		r12 = PROPS(ind+8) - fac*(PROPS(ind+8)-PROPS(indm1+8))
		r23 = PROPS(ind+9) - fac*(PROPS(ind+9)-PROPS(indm1+9))
		r31 = PROPS(ind+10) - fac*(PROPS(ind+10)-PROPS(indm1+10))
       elseif(i_yield.ge.3) then
		al1 = PROPS(ind+5) - fac*(PROPS(ind+5)-PROPS(indm1+5))
       	al2 = PROPS(ind+6) - fac*(PROPS(ind+6)-PROPS(indm1+6))
       	al3 = PROPS(ind+7) - fac*(PROPS(ind+7)-PROPS(indm1+7))
       	al4 = PROPS(ind+8) - fac*(PROPS(ind+8)-PROPS(indm1+8))
       	al5 = PROPS(ind+9) - fac*(PROPS(ind+9)-PROPS(indm1+9))
       	al6 = PROPS(ind+10) - fac*(PROPS(ind+10)-PROPS(indm1+10))
       	al7 = PROPS(ind+11) - fac*(PROPS(ind+11)-PROPS(indm1+11))
       	al8 = PROPS(ind+12) - fac*(PROPS(ind+12)-PROPS(indm1+12))
       	al9 = PROPS(ind+13) - fac*(PROPS(ind+13)-PROPS(indm1+13))
       	al10 = PROPS(ind+14) - fac*(PROPS(ind+14)-PROPS(indm1+14))
       	al11 = PROPS(ind+15) - fac*(PROPS(ind+15)-PROPS(indm1+15))
       	al12 = PROPS(ind+16) - fac*(PROPS(ind+16)-PROPS(indm1+16))
       	al13 = PROPS(ind+17) - fac*(PROPS(ind+17)-PROPS(indm1+17))
       	al14 = PROPS(ind+18) - fac*(PROPS(ind+18)-PROPS(indm1+18))
       	al15 = PROPS(ind+19) - fac*(PROPS(ind+19)-PROPS(indm1+19))
       	al16 = PROPS(ind+20) - fac*(PROPS(ind+20)-PROPS(indm1+20))
       	al17 = PROPS(ind+21) - fac*(PROPS(ind+21)-PROPS(indm1+21))
       	al18 = PROPS(ind+22) - fac*(PROPS(ind+22)-PROPS(indm1+22))
       	pem = PROPS(ind+23) - fac*(PROPS(ind+23)-PROPS(indm1+23))
	endif
c
      abac1 = PROPS(ind+24) - fac*(PROPS(ind+24)-PROPS(indm1+24))
      abac2 = PROPS(ind+25) - fac*(PROPS(ind+25)-PROPS(indm1+25))
c
      if(abs(abac2) .lt. smtol) then
		abac12 = zero
       else
       abac12 = abac1/abac2
       endif
       eratlim=PROPS(ind+26) - fac*(PROPS(ind+26)-PROPS(indm1+26))
       erexp = PROPS(ind+27) - fac*(PROPS(ind+27)-PROPS(indm1+27))
       delk = PROPS(5)   
       if(i_yield .eq. 1) then
C NEW DEFINITIONS FOR HILL'S CONSTANTS-MATCHES ABAQUS AND HOSFORD & CADDELL
	r112 = one/(r11*r11)
	r222 = one/(r22*r22)
	r332 = one/(r33*r33)
         a1 = (r222 + r332 - r112)/two
         a2 = (r332 + r112 - r222)/two
         a3 = (r112 + r222 - r332)/two
         a4 = three/(two*r12*r12)
         a5 = three/(two*r23*r23)
         a6 = three/(two*r31*r31)
       endif
       esp(1) = DSTRAN(1)
       esp(2) = DSTRAN(2)
       esp(3) = DSTRAN(3)
       esp(4) = DSTRAN(4)
       esp(5) = DSTRAN(5)
       esp(6) = DSTRAN(6)
c
       xp_ep = STATEV(1)
	ri = STATEV(2)
       a_sn(1) = STATEV(3)   !! change
       a_sn(2) = STATEV(4)   !! change
       a_sn(3) = STATEV(5)   !! change
       a_sn(4) = STATEV(6)   !! change
       a_sn(5) = STATEV(7)   !! change
       a_sn(6) = STATEV(8)   !! change
c
       sn(1) = STRESS(1) - a_sn(1)
       sn(2) = STRESS(2) - a_sn(2)
       sn(3) = STRESS(3) - a_sn(3)
       sn(4) = STRESS(4) - a_sn(4)
       sn(5) = STRESS(5) - a_sn(5)
       sn(6) = STRESS(6) - a_sn(6)
c      ***************************************************************
c      -------------------------------------------------
c      define of elastic modulus tensor for plane stress
c      -------------------------------------------------
       do 10 ii=1,6
	 do 10 jj=1,6
10        cmatl(ii,jj)=zero
c
       p1 = young/(1.0 - xnu*xnu)
       p2 = young/(1.0 + xnu)
       p3 = p2/2.0
       pf0 = young/(1.0 + xnu)/(1.0 - 2.0*xnu)
       pf1 = pf0*(1.0 - xnu)
       pf2 = pf0*xnu
       pf3 = p3
c
       cmatl(1,1) = pf1
       cmatl(1,2) = pf2
	cmatl(1,3) = pf2
       cmatl(2,1) = pf2
       cmatl(2,2) = pf1
	cmatl(2,3) = pf2
       cmatl(3,1) = pf2
       cmatl(3,2) = pf2
	cmatl(3,3) = pf1
       cmatl(4,4) = pf3
       cmatl(5,5) = pf3
       cmatl(6,6) = pf3
c
	do ii = 1,6
		do jj = 1,6
	       cfelas(ii,jj) = cmatl(ii,jj)
		enddo
	enddo
c      -------------------------------------------------------------
c      **compute stresses based on Incremental Deformation Theory** 
c        (core routine)
c      -------------------------------------------------------------
       call  kup_jaman(cmatl,cfelas,esp,dep,sn,a_sn,xp_ep,ri,dtime,
     1	pnewdt,kstep,kinc,noel,kspt,layer,iprint,props,nprops)
c      ------------------------------------------
c      return of new variables to ABAQUS implicit
c      ------------------------------------------
c	compute strain energy density variables
	do i = 1,6
	delas = esp(i) - dep(i)
	sav = (stress(i)+sn(i))/2.
	sse = sse + sav*delas
C plastic dissipation is inaccurate when stress at the beginning of an 
c increment is elastic
	spd = spd + sav*dep(i)
	enddo
       do i = 1,6
         STRESS(i) = sn(i)
       enddo
       STATEV(1) = xp_ep
       STATEV(2) = ri
       STATEV(3) = a_sn(1)   !! change
       STATEV(4) = a_sn(2)   !! change
       STATEV(5) = a_sn(3)   !! change 
       STATEV(6) = a_sn(4)   !! change
       STATEV(7) = a_sn(5)   !! change
       STATEV(8) = a_sn(6)   !! change 
c tangent matrix
	do ii = 1,6
		do jj = 1,6
	       ddsdde(ii,jj) = cmatl(ii,jj)
		enddo
	enddo
c
       return
       end
c---------------------------------------------------------------
c---------------------------------------------------------------
        subroutine kup_jaman(d,cfe,de,dep,sn,a_sn,xp_ep,ri,dtime,
     1		pnewdt,kstep,kinc,noel,kspt,layer,iprint,props,nprops)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        dimension dd(6,6),d(6,6),sn(6),si(6),de(6),xni(6),xnin(6),
     1       cfe(6,6),a_sn(6),a_si(6),dsi(6),atemp(6),sit(6),asit(6),
     1		x2n(6,6),conv1(3),conv2(3),res(7),siv(7),xnsv(7),
     1		a_sif(6),dsi2(6),dep(6),props(nprops)
        common/ki_select/i_hard,i_yield
       common /intp/ aif,aifi
        common/ksh_mtl/ams1,ams2,ams3,y_0,abac1,abac2,abac12,delk
       data qtol,xntol,delim /2.5d-4,2.5d-4,1.d-12/
       data zero,one,two,three /0.0d0,1.d0,2.d0,3.d0/
C       data iqlim,iqlim2 /5,10/
       data iqlim,iqlim2 /3,3/
C alpha intergation factor set to 0.5 for average values over increment
       data aif /0.5d0/
c       ---------------------------------------
c       Stress integration with Chaboche Model
c       ---------------------------------------
	mspt = noel
        xp_rn= zero
        pi = xp_ep
       aifi = one - aif
        do i=1,6
          si(i) = sn(i)
          a_si(i) = a_sn(i) 
       conv1(i) = zero
       conv2(i) = zero
        enddo
c trial stress using full strain increment and elastic stiffnes matrix
        do  i=1,6
          do  j=1,6
          si(i) = si(i)+d(i,j)*de(j)
          dd(i,j) = d(i,j)
          enddo
       sit(i) = si(i)
       asit(i) = a_si(i)
       enddo
c check for no incremnt in strain- use elastic props and return elastic matrix
       desum = zero
              do i=1,6
		fac = one
c shear strains are gamma's (factor of two times epsilom) but gam_ij = gamji
c so there are two contributions
		if(i .gt. 3) fac = one/two
              desum = desum + de(i)*de(i)*fac
              enddo
c effective strain increment for first trial
       desum = sqrt(two*desum/three)
	erate = desum/dtime
       if (desum .lt. delim) go to 333
c no strain rate sensitivity
        call kyshr(pi,zero,dtime,xka,hk,props,nprops)
c ignore strain rate sensitivity on first call
c	if(i_hard .lt. zero) call kyshr(pi,erate,dtime,xkar,hkr,props,nprops)
c set flag to calculate 2nd derivatives
		xka_old = xka
		eprate_old = zero
		hkr_old = hk
       i2der = 0
        if(i_yield.eq.1)then
          call khill_pot(si,xni,ri,x2n,i2der,iprint)
c        elseif(i_yield.eq.2) then 
c          call kbarlat_yld96(si,xni,ri,x2n,i2der)
        elseif(i_yield.ge.3) then
          call barlat_yld2000(si,xni,ri,x2n,i2der,mspt,iprint)
        endif
       qi = ri - xka
        qt = qi
        if(qt.lt.zero) then !elastic
          goto 333
	endif
c ignore rate sensitivity on first call
c	if(i_hard .lt. 0) then
c        qi = ri - xkar
cc strain rate overestimates yield strength
cc use rate insensitive yield strength
c	if(qi .lt. zero)  qi = ri - xka
c        qt = qi
c        endif
        xka_iso = delk*(xka - y_0) + y_0 
c flow stress at beginning of increment
       xka_isot = xka_iso
        hk_iso = delk*hk
        hk_kin = dabs((one-delk))*hk
C following lines are inserted for comparison with ABAQUS solution
C        hk_al = hk_kin
C        hk_be = y_0*dabs((1.0-delk))
        hk_al = abac1
        hk_be = abac2
        ninc = qt/(y_0)
        if(ninc.lt. 1) ninc = 1
       if(ninc .gt. 10) then
C use for problem with incipient yield
C       if(ninc .gt. 50) then
       PNEWDT = 0.56d0
       ninc = 1
       go to 333
       endif
c ninc has been set equal to one permanently
       ninc = 1
        rninc = ninc
        dq = qt/rninc
        do 25 k=1,ninc
c follwoing line incrments applies strees
       cc = k*dq
       	dqi = zero
C initialize subiteration counter
       iq = 0
60	continue
C r2 is redefined in this approach
          r2 = hk_iso
          cc = cc + dqi
c
          do 35 ij=1,6
          do 35 ik=1,6
C original
C            r2 = r2 + xni(ij)*d(ij,ik)*xni(ik)
C using tangent stiffnes matrix for d sig eq 7.1
            r2 = r2 + xni(ij)*dd(ij,ik)*xni(ik)
35        continue
          do ii=1,6
c       dum1 = si(ii)
c       dum2 = a_si(ii)
c       duma = hk_al*si(ii)/xka_iso
c       dumb = hk_be*a_si(ii)
C derivatives taken at latest estimate of yield location
            atemp(ii) = hk_al*si(ii)/xka_iso - hk_be*a_si(ii)
            r2 = r2 + xni(ii)*atemp(ii)
          enddo
c
          dl = cc/r2
C xp_rn is not incremnted by dl but are dl because cc increases to 
C full dq not incremntal
C          xp_rn = xp_rn + dl
          xp_rn =  dl
c plastic strain rate 
	eprate = dl/dtime
C pi is current estimate of total plastic strain
          pi = pi + dl
C          pi = dl
          do ii=1,6
C deviatoric stresses updated from values at beginning of increment which are
C based on trial stresses
       		si(ii) = sit(ii)
            do ij=1,6
       		dgam = dl*dd(ii,ij)*xni(ij)
              si(ii) = si(ii) - dl*dd(ii,ij)*xni(ij)
       		if(abs(de(ij)) .gt. 1.d-8) then
       		ptanij = dgam/de(ij)
       		else
       		ptanij = zero
       		endif
            enddo
C compute weighhted average stress - alpha in increment
       		siav = aif*si(ii) + aifi*sn(ii)
C weighted average flow stress at beginning and end of increment
       		xkav = aif*xka_iso + aifi*xka_isot
       		sxkav = siav/xkav
       		sxki = sn(ii)/xka_isot
C            dsi(ii) = dl*(hk_al*si(ii)/xka_iso - hk_be*a_si(ii))
c following dsi is from integration of evolution equation assuming
c that si/xka_iso is constant over the interval
c       	if(abac2 .gt. zero) then
c       	si2c = a_sn(ii) - abac12*sxkav
c       	dsi(ii) = si2c*exp(-dl*abac2) - si2c
c       	else
C incorporates the weighted average alpha over the increment
       		dlf = dl/(one + dl*abac2*aif)
            dsi(ii) = dlf*(abac1*sxkav - ABAC2*a_sn(ii))
c       	endif
            a_si(ii) = a_sn(ii) + dsi(ii)
C sig-alpha must include d alpha term in update
              si(ii) = si(ii) - dsi(ii)
C not sure believe sig-alpha is correctly calculated initially, only alpha
C needs to be updated
          enddo
c
c establish limits on pi, if out-of-bounds go to gensec
       	if(pi.lt.zero .or. pi.gt.three) then
       		if(iq .eq. 0) then
c very bad first guess, return with smaller time increment
       		ifail = 1
       		pnewdt = 0.56d0
       		go to 333
       		else
c bad current guess, take best estimate and make call to gensec
       		go to 120
       		endif
       endif
c no strain rate sensitivity
        call kyshr(pi,zero,dtime,xka,hk,props,nprops)
c strain rate sensitivity
	if(i_hard .lt. zero) call kyshr(pi,eprate,dtime,xkar,hkr,props,nprops)
          if(i_yield.eq.1)then
            call khill_pot(si,xnin,ri,x2n,i2der,iprint)
c          elseif(i_yield.eq.2) then 
c            call kbarlat_yld96(si,xnin,ri,x2n,i2der)
          elseif(i_yield.ge.3) then 
          call barlat_yld2000(si,xnin,ri,x2n,i2der,mspt,iprint)
          endif
          qi = ri - xka
          xka_iso = delk*(xka - y_0) + y_0 
          hk_iso = delk*hk
          hk_kin = dabs((one-delk))*hk
c account for rate sensitivity
	if(i_hard .lt. 0) then
        qi = ri - xkar
c strain rate overestimates yield strength
c use rate insensitive yield strength & strain rate
		if((cc+qi) .lt. zero)  then
c save current estimate of stress state as starting point in gensec
	       do ii = 1,6
	       siv(ii) = si(ii)
	       xnsv(ii) = xnin(ii)
	       enddo
		epratsv = eprate
		dlsv = dl
		qisv = qi
          xka_iso = delk*(xkar - y_0) + y_0 
          hk_iso = delk*hkr
          hk_kin = dabs((one-delk))*hkr
		go to 120
		else
		xka_old = xkar
		eprate_old = eprate
		hkr_old = hkr
		endif
          xka_iso = delk*(xkar - y_0) + y_0 
          hk_iso = delk*hkr
          hk_kin = dabs((one-delk))*hkr
	endif
c check for convergence between old and new normals
c
c magnitude of old normal
       xnmag = zero
       xndif = zero
       do ik=1,6
       xnmag = xnmag + xni(ik)*xni(ik)
c differences between old and new normals
       xndif = xndif + (xnin(ik)-xni(ik))*(xnin(ik)-xni(ik))
       enddo
       xncon = sqrt(xndif/xnmag)
       qicon = abs(qi/y_0)
       if(qicon .lt. qtol .and. xncon .lt. xntol) then
       go to 25
       else
c save best estimate of stress state as starting point in gensec
       if(iq .eq. 0) then
       do ii = 1,6
       siv(ii) = si(ii)
       xnsv(ii) = xnin(ii)
       enddo
       dlsv = dl
c plastic strain rate 
	epratsv = eprate
       conmin = qicon + xncon
       qisv = qi
       else
       conc = qicon + xncon
       	if(conc .lt. conmin) then
       do ii = 1,6
       siv(ii) = si(ii)
       xnsv(ii) = xnin(ii)
       enddo
       dlsv = dl
c plastic strain rate 
	epratsv = eprate
       qisv = qi
       conmin = conc
       	endif
       endif
C Incremental approach to convergence, accepting first solution but 
C modifying it       
c       cc = qi
C Complete approach with redefined normals based on latest solution
C following line not needed because xp_rn is  not solved for incremntally
C          xp_rn = xp_rn - dl
C pi is current estimate of total plastic strain and must be rezeroed
          pi = pi - dl
c redefine normals as average of old and new values
       do ik = 1,6
      xni(ik) = (xni(ik) + xnin(ik))/two
       enddo
C add differences of flow stress and yield function to qt
       dqi = qi
       iq = iq + 1
c redefine latest convergence values
       		do im = 1,2
       		im3 = 4 - im
       		conv1(im3) = conv1(im3-1)
       		conv2(im3) = conv2(im3-1)
       		enddo
       		conv1(1) = xncon
       		conv2(1) = qicon
       		if(iq .le. iqlim) then
c use gensec for rate dependent problems 
		if(i_hard .lt. 0 .and. iq .ge. 2) go to 120
       		go to 60
              elseif(iq.gt.iqlim .and. iq.le.iqlim2) then
c check rate of convergence
       		if(conv1(1) .gt. xntol) then
       		call conest(conv1,xntol,xnit,mspt)
       		else 
       		xnit = two
       	        endif
       		if(conv2(1) .gt. qtol) then
       		call conest(conv2,qtol,qit,mspt)
       		else
       		qit = two
       		endif
       		if(xnit.gt.three .or. qit.gt.three) go to 120
c resume bisection method if rate of convergence is ok
       		go to 60
       		else
c call general secant method
120		res(7) = qisv
          do ii=1,6
C deviatoric stresses updated from values at beginning of increment which are
C based on trial stresses
C use best estimate from bisection method
       		delst = zero
            do ij=1,6
              delst = delst - dlsv*dd(ii,ij)*xnsv(ij)
       	    enddo
C compute weighted average stress - alpha in increment
       		siav = aif*si(ii) + aifi*sn(ii)
C weighted average flow stress at beginning and end of increment
       		xkav = aif*xka_iso + aifi*xka_isot
       		sxkav = siav/xkav
       		sxki = sn(ii)/xka_isot
C uses integration of evolution equation assuming si/xka=constant
C       	if(abac2 .gt. zero) then
C       	si2c = a_sn(ii) - abac12*sxkav
C       	delal = si2c*exp(-dl*abac2) - si2c
C       	else
C incorporates the weighted average alpha over the increment
C chanhe on 2003-05-01 dlsv in denom
       		dlf = dlsv/(one + dlsv*abac2*aif)
C       		dlf = dlsv/(one + dl*abac2*aif)
C initial guess for delta alpha
            delal = dlf*(abac1*sxkav - ABAC2*a_sn(ii))
C       	endif
C            delal = dlsv*(hk_al*siv(ii)/xka_iso - hk_be*asit(ii))
C sig-alpha must include d alpha term in update
              res(ii) = siv(ii) -(sit(ii)+delst-delal)
       	    enddo
       		ifail = 0
       		call gensec(dd,siv,sit,sn,asit,a_sn,res,dlsv,xp_ep,
     1		xka_isot,epratsv,dtime,ifail,mspt,iprint,props,nprops)
             if(ifail .gt. 0) then
       PNEWDT = 0.56d0
       go to 333
		elseif(ifail .lt. 0) then
       go to 334
       endif
c update quantities after gensec
       		dl = dlsv
       		pi = xp_ep + dl
          xp_rn =  dl
		eprate = epratsv
c update qi with best estimate
          call kyshr(pi,eprate,dtime,xka,hk,props,nprops)
          xka_iso = delk*(xka - y_0) + y_0 
          hk_iso = delk*hk
          hk_kin = dabs((1.0-delk))*hk
       	    do ii = 1,6
       		si(ii) = siv(ii)
       	enddo
       	    do 315 ii = 1,6
C compute weighted average stress - alpha in increment
       		siav = aif*si(ii) + aifi*sn(ii)
C weighted average flow stress at beginning and end of increment
       		xkav = aif*xka_iso + aifi*xka_isot
       		sxkav = siav/xkav
       		sxki = sn(ii)/xka_isot
C            dsi(ii) = dl*(hk_al*si(ii)/xka_iso - hk_be*a_si(ii))
c following dsi is from integration of evolution equation assuming
c that si/xka_iso is constant over the interval
C       	if(abac2 .gt. zero) then
C       	si2c = a_sn(ii) - abac12*sxkav
C       	dsi(ii) = si2c*exp(-dl*abac2) - si2c
C       	else
C incorporates the weighted average alpha over the increment
       		dlf = dl/(one + dl*abac2*aif)
            dsi(ii) = dlf*(ABAC1*sxkav - ABAC2*a_sn(ii))
C       	endif
            a_si(ii) = a_sn(ii) + dsi(ii)
       		denom = one + dl*hk_be
       	a_sif(ii) = a_sn(ii) + dsi(ii)/denom 
315	continue
c            enddo
         if(i_yield.eq.1)then
            call khill_pot(si,xnin,ri,x2n,i2der,iprint)
c          elseif(i_yield.eq.2) then 
c            call kbarlat_yld96(si,xnin,ri,x2n,i2der)
          elseif(i_yield.ge.3) then 
          call barlat_yld2000(si,xnin,ri,x2n,i2der,mspt,iprint)
          endif
c
          qi = ri - xka
             if(ifail .gt. 0) go to 333
              go to 25
              	endif
       endif
c
25      continue
c
c Calculate quantity to update tangent modulus
       	scal2 = zero
       scal2b = zero
          do ii=1,6
C derivatives taken at latest estimate of yield location
            alphar = hk_al*si(ii)/xka_iso - hk_be*a_si(ii)
            scal2 = scal2 + xnin(ii)*alphar
       	    scal2b = scal2b + xnin(ii)*(si(ii)-sn(ii))
          enddo
       		scal2b = scal2b/dl
c set flag to calculate 2nd derivatives
       i2der = 1
          if(i_yield.eq.1)then
            call khill_pot(si,xnin,ri,x2n,i2der,iprint)
c          elseif(i_yield.eq.2) then 
c            call kbarlat_yld96(si,xnin,ri,x2n,i2der)
          elseif(i_yield.ge.3) then 
          call barlat_yld2000(si,xnin,ri,x2n,i2der,mspt,iprint)
          endif
        call kal_tangent(dd,cfe,xnin,x2n,si,a_si,ri,dl,scal2,
     1		scal2b,pi,eprate,dtime,iprint,mspt,props,nprops)    
        do i=1,6
          do j=1,6
              d(i,j) = dd(i,j)
          enddo
        enddo
        xp_ep = xp_ep + xp_rn
c
c calculate plastic strain increments
	do i = 1,6
	dep(i) = xp_rn*xnin(i)
	enddo	
333     continue
c
         do i=1,6
C sn now becones stress only was stres-backstress at beginning of increment
          sn(i) = si(i) + a_si(i)
          a_sn(i) = a_si(i)     
         enddo
c
         return
334	continue
c elastic solution is indicated after initial plastic estimate
         do i=1,6
          sn(i) = sit(i) + asit(i)
          a_sn(i) = asit(i)     
         enddo
         return
         end
c---------------------------------------------------------------
       subroutine conest(y,tol,x,mspt)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       dimension y(3)
       data zero,one,two,four,hun /0.d0,1.d0,2.d0,4.d0,100.d0/
c determine # of additional iterations using quadratic convergence
	 xlin = hun
       xquad = hun
       den = y(1)-y(2)
       if(abs(den) .gt. 1.d-10)  xlin = (tol-y(2))/den
       if(xlin .lt. 1.d0) xlin = hun
       a = y(3)
       c = (-two*y(2)+y(1)+y(3))/two
       b = y(2)-y(3)-c
       if(abs(c) .lt. 1.d-10) then
       go to 100
       else
       rad = b*b - four*(a-tol)*c
              if(rad .lt. zero) then
c convergence unlikely
              xquad = hun
       		go to 100
              endif
       xquad = (-b+sqrt(rad))/two/c
c convergence unlikely
       if(xquad .lt. two) xquad = hun
       endif
c compute number of additional iterations estimated
100       xquad = xquad - two
        xlin = xlin - one
c take minimum of linear and quadratic convergence estimates
       x = xlin
       if(xquad .lt. x) x = xquad
       return
       end       
c---------------------------------------------------------------
       subroutine gensec(DD,SI,SIT,SN,ASIT,ASN,RES,CHI,PI,XKA_ISOT,
     1	eprate,dtime,IFAIL,mspt,iprint,props,nprops)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C THIS SUBROUTINE USES A GENERALIZED SECANT METHOD TO SOLVE FOR THE STRESS 
C STATE OF A MATERIAL UNDERGOING 3D DEFORMATION.
C THERE ARE seven UNKNOWNS, A,B,C,d,e,f AND CHI. Seven INITIAL ESTIMATES ARE READ IN. 
C THE FIRST six EQUATIONS COME FROM THE NORMALITY FLOW RULE. THE LAST
C EQUATION ENFORCES EQUALITY BETWEEN THE YIELD STRESS AND FLOW STRESS.
        common/ki_select/i_hard,i_yield
        common/ksh_mtl/ams1,ams2,ams3,y_0,abac1,abac2,abac12,delk
       common /intp/ aif,aifi
	DIMENSION DD(6,6),SI(6),SIT(6),ASIT(6),RES(7),YES(8,7),
     1		HRES(8,7),AM(28),BM(28),CM(28),JADIAG(7),YNEW(7),
     1		xni(6),x2n(6,6),sn(6),asn(6),props(nprops)
	LOGICAL AFAC,BACK
C       DATA ERLIM,RFAC,zero,one,two/2.5d-3,0.01d0,0.d0,1.d0,2.d0/
C for rate dependence
C       DATA ERLIM,RFAC,zero,one,two/5.0d-5,0.01d0,0.d0,1.d0,2.d0/
c for non-rate problems
       DATA ERLIM,RFAC,zero,one,two/1.d-7,0.01d0,0.d0,1.d0,2.d0/
C       DATA ERLIM,RFAC,zero,one,two/1.d-10,0.01d0,0.d0,1.d0,2.d0/
	DATA CHIFAC,CHITOL,BIGNUM /1.D-1,5.D-1,1.D+3/
       data ILIM /20/
c       data ILIM /25/
	DATA JADIAG/1,3,6,10,15,21,28/
C INITIALIZE
	AFAC = .TRUE.
	BACK = .TRUE.
c reduce erlim for rate dependent problems
	if(i_hard .lt. 0) ERLIM = 5.0d-5
	I = 0
	ISNFLG = 0
       ry0 = rfac*y_0
C BASE VECTOR
	DO 90 K = 1,6
90	YES(1,K) = SI(K)
       YES(1,7) = CHI/CHIFAC
C       YES(1,7) = DLOG(CHI)
	DO 95 K = 1,7
95	HRES(1,K) = RES(K)
C COMPUTE S MATRIX
100	I = I + 1
C ON FIRST CALL, CALCULATE AN INITIAL GUESS FOR CHI.  THE NORMALIZATION TERM
C WILL REMAIN CONSTANT DURING THE INITIAL SEARCH ROUTINE.  THE ESTIMATE FOR
C CI IS COMPUTED BY DIVIDING THE APPROXIMATE WORK RATE BY THE INNER PRODUCT
C OF THE STRESS TENSOR TIMES THE DERIVATIVE OF THE PLASTIC POTENTIAL WRT STRESS.
	IF(I .EQ. 1) GO TO 115
C COMPUTE RESIDUAL VECTOR
C limit estimate of plastic strain to positive numbers, if negative, return 
c with reduction in time step
       if(epn .lt. zero ) then
c call hardening law with zero strain and strain rate
        call kyshr(zero,zero,dtime,xka,hk,props,nprops)
	elseif (eprate .lt. zero) then
c call hardening law with zero strain rate
        call kyshr(EPN,zero,dtime,xka,hk,props,nprops)
	else
        call kyshr(EPN,eprate,dtime,xka,hk,props,nprops)
	endif
        xka_iso = delk*(xka - y_0) + y_0 
       i2der = 0
        if(i_yield.eq.1)then
          call khill_pot(SI,XNI,ri,x2n,i2der,iprint)
c        elseif(i_yield.eq.2) then 
c          call kbarlat_yld96(SI,xni,ri,x2n,i2der)
        elseif(i_yield.ge.3) then
          call barlat_yld2000(si,xni,ri,x2n,i2der,mspt,iprint)
        endif
C DIFFERENCE BETWEEN FLOW STRESS AND YIELD STRESS
       RES(7) = RI - XKA
	IF(ISNFLG .GT. 0) THEN
C PENALTY NUMBER FOR CONSTRAINT ON NEGATIVE PLATIC STRAIN INCREMENTS
C	RES(7) = RES(7) - (eprate/chitol)*BIGNUM
C USE QUADRATIC PENALTY FACTOR
	RES(7) = RES(7) + BIGNUM*(abs(eprate)/chitol)**2.0
	ISNFLG = 0
	ENDIF
C DIFFERENCE BETWEEN CURRENT STRESS VALUES AND CALCULATED ONES 
          do ii=1,6
C deviatoric stresses updated from values at beginning of increment which are
C based on trial stresses
       		delst = zero
            do ij=1,6
              delst = delst - CHI*dd(ii,ij)*xni(ij)
       	    enddo
C CALCULATE CHANGE IN BACK STRESS
C            delal = CHI*(ABAC1*SI(ii)/xka_iso - ABAC2*asn(ii))
C compute weighted average stress - alpha in increment
       		siav = aif*si(ii) + aifi*sn(ii)
C weighted average flow stress at beginning and end of increment
       		xkav = aif*xka_iso + aifi*xka_isot
       		sxkav = siav/xkav
       		sxki = sn(ii)/xka_isot
c following dsi is from integration of evolution equation assuming
c that si/xka_iso is constant over the interval
C       	if(abac2 .gt. zero) then
C       	si2c = asn(ii) - abac12*sxkav
C       	delal = si2c*exp(-chi*abac2) - si2c
C       	else
C incorporates the weighted average alpha over the increment
       		chif = chi/(one + chi*abac2*aif)
            delal = CHIF*(ABAC1*sxkav - ABAC2*asn(ii))
C       	endif
C sig-alpha must include d alpha term in update
              RES(ii) = SI(ii) -(sit(ii)+delst-delal)
       	    enddo
C DETERMINE CONVERGENCE
115	ERSUM = ZERO
	DO 110 J = 1,7
C SET UPPER LIMIT ON RESIUDAL DUE TO MACHINE TOLERANCE
	IF(RES(J) .GT. 1.D19) GO TO 270
110	ERSUM = ERSUM + RES(J)*RES(J)
	ERSUM = DSQRT(ERSUM)/y_0
	IF(ERSUM .LE. ERLIM) THEN
	RETURN
	ELSE
		IF(I .EQ. 1) THEN
		ERSUMIN = ERSUM
		AMIN = YES(1,1)
		BMIN = YES(1,2)
		CMIN = YES(1,3)
		DMIN = YES(1,4)
		EMIN = YES(1,5)
		FMIN = YES(1,6)
		CHIMIN = YES(1,7)*CHIFAC
C		CHIMIN = EXP(YES(1,7))
		epratsv = chimin/dtime
		ELSE
		IF(ERSUM .LT. ERSUMIN) THEN
		ERSUMIN = ERSUM
		AMIN = SI(1)
		BMIN = SI(2)
		CMIN = SI(3)
		DMIN = SI(4)
		EMIN = SI(5)
		FMIN = SI(6)
		CHIMIN = CHI
		epratsv = chimin/dtime
		ENDIF
		ENDIF
		IF(I .LE. ILIM) GO TO 260
270		IFAIL = 1
C		GO TO 270
		SI(1) = AMIN
		SI(2) = BMIN
		SI(3) = CMIN
		SI(4) = DMIN
		SI(5) = EMIN
		SI(6) = FMIN
C SOLUTION INDICATES NO YIELDING
		IF(RES(7) .LT. ZERO .OR. EPRATE .LT. ZERO) IFAIL = -1
		CHI = CHIMIN
C		eprate = epratsv
		RETURN
260		CONTINUE
	ENDIF
	IF(I-8) 10,20,30
10       IF(I .EQ. 1) GO TO 125
	DO 120 J = 1,7
120	HRES(I,J) = RES(J)
C RFAC IS USED TO MAKE SUBSEQUENT ESTIMATES, CURRENTLY RFAC IS LIMITED
C BY  .001<RFAC<.1 BUT THESE LIMITS MAY BE ADJUSTED
125	DO 130 J = 1,7
130	YES(I+1,J) = YES(1,J)
	YES(I+1,I) = YES(1,I)*(ONE+RFAC)
c if a stress component is zero, then the new vector is not changed by the 
c above method
       if(i.le.6 .and. abs(yes(i+1,i)).lt. ry0) yes(i+1,i)=ry0
       CHI = YES(I+1,7)*CHIFAC
C       CHI = EXP(YES(I+1,7))
       EPN = PI + CHI
	eprate = chi/dtime
       DO K = 1,6
       SI(K) = YES(I+1,K)
       ENDDO
	GO TO 100
30	IF(I .GT. 9) THEN
C REARRANGE COLUMNS
	DO 140 IR = 8,2,-1
	DO 140 J = 1,7
140	HRES(IR,J) = HRES(IR-1,J)
	ENDIF
	DO 150 J = 1,7
150	HRES(1,J) = RES(J)
	GO TO 159
20	DO 155 J = 1,7
155	HRES(I,J) = RES(J)
C CALCULATE NEW YES VECTOR
159	JCNT = 0
	DO 160 IR = 1,7
	DO 170 J = 1,IR
	JCNT = JCNT + 1
	AM(JCNT) = HRES(9-IR,J) - HRES(1,J)
	IF(IR .EQ. J) THEN
	CM(JCNT) = ONE
	ELSE
	CM(JCNT) = HRES(9-J,IR) - HRES(1,IR)
	ENDIF
170	CONTINUE
160	BM(IR) = -HRES(1,IR)
C SOLVE EQUATIONS FOR FACTORS
	CALL UACTCL(AM,CM,BM,JADIAG,7,AFAC,BACK,mspt,iprint)
C COMPUTE NEW SOLUTION VECTOR
	DO 210 IR = 1,7
	YNEW(IR) = YES(1,IR)
	DO 210 J = 1,7
210	YNEW(IR) = YNEW(IR) + BM(J)*(YES(9-J,IR)-YES(1,IR))
C REARRANGE Y VECTORS
	IF( I .GT. 8) THEN
	DO 220 IR = 1,7
	DO 220 J = 1,7
220	YES(9-J,IR) = YES(8-J,IR)
	ENDIF
C SAVE NEW ESTIMATE
	DO 230 IR = 1,7
230	YES(1,IR) = YNEW(IR)
C REDEFINE NEW SOLUTION ESTIMATE
       CHI = YES(1,7)*CHIFAC
C       CHI = EXP(YES(1,7))
c for slightly small negative strain increments, reset CHI, activate penalty
c constraint 
	if(chi .lt. zero)then
	ISNFLG = 1
	endif
       EPN = PI + CHI
	eprate = chi/dtime
       DO K = 1,6
       SI(K) = YES(1,K)
       ENDDO
	GO TO 100
	END
C***********************************************************************
        SUBROUTINE UACTCL(A,C,B,JDIAG,NEQ,AFAC,BACK,mspt,iprint)
C***********************************************************************
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C THIS SUBROUTINE SOLVES THE SYSTEM OF NON-SYMMETRIC EQUATIONS
C TAKEN FROM F. STASA'S "APPLIED FINITE ELEMENT ANALYSIS FOR
C ENGINEERS", HOLT,RINEHART,&WINSTON, 1985.
        LOGICAL AFAC,BACK
        DIMENSION A(1),C(1),B(1),JDIAG(1)
       DATA ZERO /0.D0/
C FACTOR A TO UT*D*U, REDUCE B TO Y
        JR = 0
        DO 300 J = 1,NEQ
        JD = JDIAG(J)
        JH = JD - JR
        IF(JH .LE. 1) GO TO 300
        IS = J + 1 - JH
        IE = J - 1
        IF(.NOT. AFAC) GO TO 250
        K = JR + 1
        ID = 0
C       REDUCE ALL EQUATIONS EXCEPT DIAGONAL
        DO 200 I = IS,IE
        IR = ID
        ID = JDIAG(I)
        IH = JMIN0(ID-IR-1,I-IS)
        IF(IH .EQ. 0) GO TO 150
        A(K) = A(K) - DOT(A(K-IH),C(ID-IH),IH)
        C(K) = C(K) - DOT(C(K-IH),A(ID-IH),IH)
150     IF(A(ID) .NE. ZERO) C(K) = C(K)/A(ID)
200     K = K + 1
C       REDUCE DIAGONAL TERM
        A(JD) = A(JD) - DOT(A(JR+1),C(JR+1),JH-1)
C       FORWARD REDUCE THE R.H.S.
250     IF(BACK) B(J) = B(J) - DOT(C(JR+1),B(IS),JH-1)
300     JR = JD
        IF(.NOT.BACK) RETURN
C       BACK SUBSTITUTION
        J = NEQ
        JD = JDIAG(J)
500     IF(A(JD) .NE. ZERO) B(J) = B(J)/A(JD)
        D = B(J)
        J = J - 1
        IF(J .LE. 0) RETURN
        JR = JDIAG(J)
        IF((JD-JR) .LE. 1) GO TO 700
        IS = J - JD + JR + 2
        K = JR - IS + 1
        DO 600 I = IS,J
600     B(I) = B(I) - A(I+K)*D
700     JD = JR
        GO TO 500
        END
C***********************************************************************
        FUNCTION DOT(A,B,N)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(1),B(1)
C THIS FUNCTION COMPUTES THE DOT PRODUCT
        DOT = 0.D0
        DO 100 I = 1,N
        DOT = DOT + A(I)*B(I)
100     CONTINUE
        RETURN
        END
c---------------------------------------------------------------------------
      subroutine kal_tangent(c,cfe,xn,x2n,s,a_si,q,dl,scal2,
     1		scal2b,pepsbr,eprate,dtime,iprint,mspt,props,nprops)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        common/ksh_mtl/ams1,ams2,ams3,y_0,abac1,abac2,abac12,delk
       common/elas/ ekp,esh
      dimension c(6,6),ts(6,6),tsi(6,6),
     1		cfe(6,6),temp(6,6),xn(6),temp2(6,6),
     1		z(6,6),x2n(6,6),s(6),sd(6),xnmat(6,6),xmu(6),
     1 		a_si(6),xhvec(7),tempv(6),xnvec(6),xmvec(6),
     1		dumn(6),dumm(6),props(nprops)
       data zero,one,two,three /0.d0,1.d0,2.0d0,3.0d0/
c     -------------------------------
c     Jacobian based on current state
c     -------------------------------
      call kyshr(pepsbr,eprate,dtime,ys,hr,props,nprops)
c adjust moduli to include alpha terms
       scal2 = scal2 + hr
C !!NOTE!! the definition of the inverse matrix in the ABAQUS theory manual for 
C H_hat and w_hat is wrong; the one term should be dha/dHb, not with a and b 
C subscripts reversed
C define addtional terms for N matrix- xnmat
C psi2 term comes from dh/dsig
       psi2 = abac1/ys
C psi term comes from dh/dH where H is vecotr of internal variables
C (eps_pl,alp1,alp2,alp3,alp12,alp23,alp31)
       psi = -psi2 - abac2
       xhfacd = one-dl*psi
       phi = psi2/xhfacd
C minus sign in factor comes from the fact that d^2g/dsds=-d^2g/dsdH
       facn = one - dl*phi
       do ii = 1,6
       do jj = 1,6
       xnmat(ii,jj) = x2n(ii,jj)*facn       
       enddo
       enddo
C define addtional terms for H_hat vector- inverse matrix*h_vector
C xhvec is H_hat in abaqus
       xmucon = dl*psi2*hr/ys
       do ii = 1,6
       xmu(ii) = s(ii)*xmucon + abac1*s(ii)/ys - abac2*a_si(ii)
       xhvec(ii+1) = xmu(ii)/xhfacd
       enddo
       xhvec(1) = one
C define additonal terms for n vector
c comnpute d^2g/dsdH*H
c first column of the matrix d^2g/dsdH is zero since g independent of eps,
c which is the first hardening variable
       do ii = 1,6
       tempv(ii) = zero
       do jj = 1,6
       tempv(ii) = tempv(ii) - x2n(ii,jj)*xhvec(jj+1)
       enddo
       xnvec(ii) = xn(ii) + dl*tempv(ii)
       enddo
C define additonal terms for m vector
       xmfac = one - dl*phi
       do ii = 1,6
       xmvec(ii) = xn(ii)*xmfac
       enddo
C define additonal terms for scalar d
C following minus sign comes from df = d(ys(sig)-fp(eps)) -d(fp) since
C xhvec(1) is for eps only
       dadd = -hr*xhvec(1)
       do ii = 1,6
C following minus sign comes from df/ds = d(ys(sig)-fp(eps))/dalp= -d(ys)/dsig
       dadd = dadd - xn(ii)*xhvec(ii+1)
       enddo
C dot product of elastic matrix and xnvec
C then xmvec with the result
       xd = zero
       do ii = 1,6
       tempv(ii) = zero
       do jj = 1,6
       tempv(ii) = tempv(ii) + cfe(ii,jj)*xnvec(jj)
       enddo
       xd = xd + tempv(ii)*xmvec(ii)
       enddo
C following minus sign comes from ABAQUS definition for "d"
       xd = xd - dadd
c define Z matrix in ABAQUS
       do ii = 1,6
       tempv(ii) = zero
C dot product of xmvec and Del
       do jj = 1,6
       tempv(ii) = tempv(ii) + cfe(jj,ii)*xmvec(jj)
       enddo
       enddo
       do 20 i = 1,6
       do 20 m = 1,6
       ck = zero
       if(i .eq. m) ck = one
       z(i,m) = ck - xnvec(i)*tempv(m)/xd
20	continue
C combine terms to define I + lambda*Del:Z:N matrix from ABAQUS
C original method
CC first do inner product of Z:N=ZN
       call matmul6(z,xnmat,temp)
CC next do inner product of Del:ZN
       call matmul6(cfe,temp,temp2)
       do 70 i = 1,6
       do 70 j = 1,6
       ck = zero
       if(i .eq. j) ck = one
	delam =dl
70     ts(i,j) = ck + delam*temp2(i,j)
C invert this matrix
       call inver6(ts,tsi)
C combine terms to define tsi:Del:Z matrix from ABAQUS
C do inner product of Del:Z
       call matmul6(cfe,z,temp2)
C do inner product of tsi:DelZ
       call matmul6(tsi,temp2,c)
      return
      end
c------------------------------------------------------------------
	SUBROUTINE INVER6(A,AI)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C THIS SUBROUTINE INVERTS A 6 BY 6 MATRIX
	DIMENSION A(6,6),BS(3,3),ES(3,3),BI(3,3),EI(3,3),CS(3,3),DS(3,3)
     1		,XI(3,3),UI(3,3),TMP(3,3),DUM(3,3),B(6,6),AI(6,6)
       data zero,one,two,three,four /0.d0,1.d0,2.0d0,3.0d0,4.d0/
C ALTERNATIVE INVERSION SCHEME
	DO 120 I = 1,6
	 DO 110 J = 1,6
	  B(I,J) = A(I,J)
110	 CONTINUE
120	CONTINUE
	DO 170 K = 1,6
	 B(K,K) = ONE/B(K,K)
	 DO 130 J = 1,6
	 IF(J.EQ.K) GO TO 130
	 B(K,J) = -B(K,K)*B(K,J)
130	 CONTINUE
	DO 150 I = 1,6
	 DO 140 J = 1,6
	  IF((I.EQ.K) .OR. (J.EQ.K)) GO TO 140
	  B(I,J) = B(I,J) + B(K,J)*B(I,K)
140	 CONTINUE
150	CONTINUE
	DO 160 I = 1,6
	 IF(I .EQ. K) GO TO 160
	 B(I,K) = B(K,K)*B(I,K)
160	CONTINUE
170	CONTINUE
	DO 180 I = 1,6
	 DO 190 J = 1,6
	  AI(I,J) = B(I,J)
190	 CONTINUE
180	CONTINUE
	IF(AI(1,1) .NE. ZERO) RETURN
	DO 10 I = 1,3
	DO 10 J = 1,3
10	BS(I,J) = A(I,J)
	CALL INVER3(BS,BI)
	DO 20 I = 1,3
	DO 20 J = 1,3
20	ES(I,J) = A(I+3,J+3)
	CALL INVER3(ES,EI)
	DO 30 I = 1,3
	DO 30 J = 1,3
30	DS(I,J) = A(I+3,J)
	DO 40 I = 1,3
	DO 40 J = 1,3
40	CS(I,J) = A(I,J+3)
	CALL MM3(EI,DS,DUM)
	CALL MM3(CS,DUM,TMP)
	DO 100 I = 1,3
	DO 100 J = 1,3
100	BS(I,J) = BS(I,J) - TMP(I,J)
	CALL MM3(BI,CS,DUM)
	CALL MM3(DS,DUM,TMP)
	DO 200 I = 1,3
	DO 200 J = 1,3
200	ES(I,J) = ES(I,J) - TMP(I,J)
	CALL INVER3(BS,XI)
	CALL INVER3(ES,UI)
	DO 210 I = 1,3
	DO 210 J = 1,3
210	AI(I,J) = XI(I,J)
	DO 220 I = 1,3
	DO 220 J = 1,3
220	AI(I+3,J+3) = UI(I,J)
	CALL MM3(CS,UI,DUM)
	CALL MM3(BI,DUM,TMP)
	DO 230 I = 1,3
	DO 230 J = 1,3
230	AI(I,J+3) = - TMP(I,J)
	CALL MM3(DS,XI,DUM)
	CALL MM3(EI,DUM,TMP)
	DO 240 I = 1,3
	DO 240 J = 1,3
240	AI(I+3,J) = - TMP(I,J)
	RETURN
	END
C****************************************************************************
	SUBROUTINE INVER3(A,AI)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C THIS SUBROUTINE INVERTS A 3 BY 3 MATRIX
	DIMENSION A(3,3),AI(3,3),U(2,2),EI(2,2)
       data zero,one,two,three,four /0.d0,1.d0,2.0d0,3.0d0,4.d0/
	DET = A(2,2)*A(3,3) - A(2,3)*A(3,2)
	EI(1,1) = A(3,3)/DET
	EI(2,2) = A(2,2)/DET
	EI(1,2) = -A(2,3)/DET
	EI(2,1) = -A(3,2)/DET
	X = A(1,1) - A(1,2)*(EI(1,1)*A(2,1)+EI(1,2)*A(3,1))
     1		- A(1,3)*(EI(2,1)*A(2,1)+EI(2,2)*A(3,1))
	AI(1,1) = ONE/X
	DO 100 I = 1,2
	DO 100 J = 1,2
100	U(I,J) = A(I+1,J+1) - A(I+1,1)*A(1,J+1)/A(1,1)
	DET = U(1,1)*U(2,2) - U(1,2)*U(2,1)
	AI(2,2) = U(2,2)/DET
	AI(3,3) = U(1,1)/DET
	AI(2,3) = -U(1,2)/DET
	AI(3,2) = -U(2,1)/DET
	AI(1,2) = -(A(1,2)*AI(2,2)+A(1,3)*AI(3,2))/A(1,1)
	AI(1,3) = -(A(1,2)*AI(2,3)+A(1,3)*AI(3,3))/A(1,1)
	AI(2,1) = -AI(1,1)*(EI(1,1)*A(2,1)+EI(1,2)*A(3,1))
	AI(3,1) = -AI(1,1)*(EI(2,1)*A(2,1)+EI(2,2)*A(3,1))
	RETURN
	END
C****************************************************************************
	SUBROUTINE MM3(A,B,C)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES
	DIMENSION A(3,3),B(3,3),C(3,3)
	DO 100 I = 1,3
	DO 100 J = 1,3
100	C(I,J) = A(I,1)*B(1,J)+A(I,2)*B(2,J)+A(I,3)*B(3,J)
	RETURN
	END
c------------------------------------------------------------------
         subroutine khill_pot(si,xn,q,x2n,i2der,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         dimension si(6),xn(6),x2n(6,6)
       common/khill/r11,r22,r33,r12,r23,r31,a1,a2,a3,a4,a5,a6
	data two /2.d0/
c        ----------------------------
c        Hill's (1948) Yield Function
c        ----------------------------
	sd12 = si(1) - si(2)
	sd23 = si(2) - si(3)
	sd31 = si(3) - si(1)
         q = (a3*sd12*sd12 + a1*sd23*sd23 + a2*sd31*sd31 
     1       + two*a4*si(4)*si(4) + two*a5*si(5)*si(5)) 
     2		+ two*a6*si(6)*si(6)
         q = sqrt(q)

        small = 10.d-7
        if((abs(q)).lt.small)q= small
C NORMAL DEFINITIONS USING SIX COMPONENTS AND INCOMPRESSIBILITY
         xn(1) = (a3*sd12 - a2*sd31)/q
         xn(2) = (a1*sd23 - a3*sd12)/q
         xn(3) = -(xn(1) + xn(2))
         xn(4) = two*a4*si(4)/q
         xn(5) = two*a5*si(5)/q
         xn(6) = two*a6*si(6)/q
c if flag is sero return
       if(i2der .eq. 0) return
c calculate 2nd derivatives in the 6x6 array with shear being 4th-6th components
       x2n(1,1) = (a3+a2 - xn(1)*xn(1))/q
       x2n(1,2) = (-a3 - xn(1)*xn(2))/q
       x2n(1,3) = -(x2n(1,1)+x2n(1,2))
       x2n(1,4) = -xn(1)*xn(4)/q
       x2n(1,5) = -xn(1)*xn(5)/q
       x2n(1,6) = -xn(1)*xn(6)/q
       x2n(2,2) = (a1+a3 - xn(2)*xn(2))/q
       x2n(2,3) = -(x2n(2,2)+x2n(1,2))
       x2n(2,4) = -xn(2)*xn(4)/q
       x2n(2,5) = -xn(2)*xn(5)/q
       x2n(2,6) = -xn(2)*xn(6)/q
       x2n(3,3) = (a1+a2 - xn(3)*xn(3))/q
       x2n(3,4) = -xn(3)*xn(4)/q
       x2n(3,5) = -xn(3)*xn(5)/q
       x2n(3,6) = -xn(3)*xn(6)/q
       x2n(4,4) = (two*a4 - xn(4)*xn(4))/q
       x2n(4,5) = -xn(4)*xn(5)/q
       x2n(4,6) = -xn(4)*xn(6)/q
       x2n(5,5) = (two*a5 - xn(5)*xn(5))/q
       x2n(5,6) = -xn(5)*xn(6)/q
       x2n(6,6) = (two*a6 - xn(6)*xn(6))/q
	do i = 2,6
		do j = 1,i-1
		       x2n(i,j) = x2n(j,i)
		enddo
	enddo
         return
         end
c---------------------------------------------------------------
        subroutine barlat_yld2000(si,xn,q,x2n,i2der,mspt,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
        dimension si(6),xn(6),x2n(6,6)
! #############################################################################
! STRESS POTENTIAL - YLD2000-3D 
! AUGUST 28, 2002 
! #############################################################################
! ATTENTION, IN THE 3D VERSION, THE FOLLOWING NOTATION APPLIES
! L1 = C prime
! L2 = C double prime
! L3 = C triple prime- not used at this time
c calculate the L matrices, L=C*T, where T is the deviator matrix, and C 
C contains the alpha parameters
       call rdani3d(mspt,iprint)
c calculate the yield strength, normals and 2nd derivatives, if necessary
       call skew3d(si,xn,q,x2n,i2der,mspt,iprint)
       return
       end
c------------------------------------------------------------------
      SUBROUTINE RDANI3D(mspt,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common/ki_select/i_hard,i_yield
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
      common/yld2_2/ RL(2,6,6),C(2,6,6),T(2,6,6)
C      INTEGER N,I,J,P,NVER
C      CHARACTER*50 COM
      DATA EPSILON /1D-10/
       data zero,one,two,three /0.d0,1.d0,2.0d0,3.0d0/
       onth = one/three
	DO N = 1 , 2
		DO I = 1 , 6
			DO J = 1 , 6
				RL(N,I,J) = zero
				C(N,I,J) = zero
			    T(N,I,J) = (I/J)*(J/I)
			END DO
		END DO
	END DO
	DO N = 1 , 2
		DO I = 1 , 3
			DO J = 1 , 3
				T(N,I,J) = T(N,I,J) - onth
			END DO
		END DO
	END DO
	if(i_yield .eq. 3) then
C use the third version of the C matrix
C the fourth version of the C matrix is commented
		C(1,1,1) = al1
C		C(1,1,2) = one - al14
C		C(1,2,1) = one - al13
		C(1,2,2) = al2
C		C(1,3,1) = -al15
C		C(1,3,2) = -al16
c s-xy or s-12
		C(1,4,4) = al7
c s-zx or s-31
		C(1,5,5) = al9
c s-yz or s-23
		C(1,6,6) = al11
C		C(2,1,1) = al3
C		C(2,1,2) = one - al18
C		C(2,2,1) = one - al17
C		C(2,2,2) = al4
C		C(2,3,1) = -al5
C		C(2,3,2) = -al6
		C(2,1,1) = two*al5
		C(2,1,2) = al6
		C(2,2,1) = al3
		C(2,2,2) = two*al4
c s-xy or s-12
		C(2,4,4) = al8
c s-zx or s-31
		C(2,5,5) = al10
c s-yz or s-23
		C(2,6,6) = al12
	elseif (i_yield .eq. 4) then
		C(1,1,1) = al1
		C(1,2,2) = al2
		C(1,3,3) = al3
c s-xy or s-12
		C(1,4,4) = al7
c s-zx or s-31
		C(1,5,5) = al9
c s-yz or s-23
		C(1,6,6) = al11
		C(2,1,2) = -al5
		C(2,1,3) = -al6
		C(2,2,1) = -al4
		C(2,2,3) = -al6
		C(2,3,1) = -al4
		C(2,3,2) = -al5
c s-xy or s-12
		C(2,4,4) = al8
c s-zx or s-31
		C(2,5,5) = al10
c s-yz or s-23
		C(2,6,6) = al12
	elseif (i_yield .eq. 5) then
		C(1,1,2) = -al1
		C(1,1,3) = -al2
		C(1,2,1) = -al3
		C(1,2,3) = -al4
		C(1,3,1) = -al5
		C(1,3,2) = -al6
c s-xy or s-12
		C(1,4,4) = al9
c s-zx or s-31
		C(1,5,5) = al8
c s-yz or s-23
		C(1,6,6) = al7
		C(2,1,2) = -al10
		C(2,1,3) = -al11
		C(2,2,1) = -al12
		C(2,2,3) = -al13
		C(2,3,1) = -al14
		C(2,3,2) = -al15
c s-xy or s-12
		C(2,4,4) = al18
c s-zx or s-31
		C(2,5,5) = al17
c s-yz or s-23
		C(2,6,6) = al16
	endif
	DO N = 1 , 2
		DO I = 1 , 6
			DO J = 1 , 6
			DO K = 1 , 6
			RL(N,I,J) = RL(N,I,J) + C(N,I,k)*T(N,k,J)
			END DO
			END DO
		END DO
	END DO
      RETURN
      END
! #############################################################################
! L1 = L prime
! L2 = L double prime
! L3 = L triple prime
      SUBROUTINE SKEW3D(SI,XN,SIG0,X2N,I2DER,mspt,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
      common/yld2_2/ RL(2,6,6),C(2,6,6),T(2,6,6)
      dimension D(6,6),x2n(6,6),xn(6),SI(6),S(6)
      dimension X1(6),X2(6),RKS1(3),RKS2(3),H1(3),H2(3)
       data zero,one,two,three /0.d0,1.d0,2.0d0,3.0d0/
C SAVE STRESSES
	DO I =1,6
	S(I) = SI(I)
	ENDDO
! CALCULATE FIRST LINEAR TRANSFORMATION
      DO I = 1 , 6
         X1(I) = ZERO
         DO J = 1 , 6
            X1(I) = X1(I) + RL(1,I,J)*S(J)
         END DO
      END DO
! CALCULATE SECOND LINEAR TRANSFORMATION
      DO I = 1 , 6
         X2(I) = ZERO
         DO J = 1 , 6
            X2(I) = X2(I) + RL(2,I,J)*S(J)
         END DO
      END DO
! CALCULATE PRINCIPAL TRANSFORMED STRESSES
      CALL CAR3D(X1,RKS1,H1,T1,EXPR1,Q1,mspt,iprint)
      CALL CAR3D(X2,RKS2,H2,T2,EXPR2,Q2,mspt,iprint)
! CALCULATE YIELD FUNCTION
      CALL YLD3D(X1,X2,RKS1,RKS2,H1,H2,S,F,SIG0,iprint)
! CALCULATE FIRST DERIVATIVES AND 2ND IF REQUIRED
      CALL RAT3D(SIG0,X1,X2,RKS1,RKS2,H1,H2,
     &           EXPR1,EXPR2,Q1,Q2,XN,X2N,I2DER,mspt,iprint)
      RETURN
      END

! #############################################################################
      SUBROUTINE CAR3D(X,RKS,H,TH,EXPR,Q,mspt,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
      common/yld2_2/ RL(2,6,6),C(2,6,6),T(2,6,6)
      DIMENSION RKS(3),H(3),X(6)
       data zero,one,two,three,four /0.d0,1.d0,2.0d0,3.0d0,4.d0/
      DATA EPSILON/1D-10/
      PI = TWO*DASIN(ONE)
       THHF = THREE/TWO
! H1, H2, H3 CALCULATION
      H(1) = (X(1)+X(2)+X(3))/THREE
      H(2) = (X(4)*X(4)+X(5)*X(5)+X(6)*X(6)
     1		-X(2)*X(3)-X(3)*X(1)-X(1)*X(2))/THREE
      H(3) = X(1)*X(2)*X(3) + TWO*X(4)*X(5)*X(6)
      H(3) = (H(3) - X(3)*X(4)*X(4)-X(1)*X(6)*X(6)-X(2)*X(5)*X(5))/TWO
! CALCULATION OF PRINCIPAL VALUES
      P = -THREE*(H(1)**TWO+H(2)) 
      Q = -(TWO*H(1)**THREE+THREE*H(1)*H(2)+TWO*H(3))
       QT = Q/TWO
C PT>0
       PT = -P/THREE
      EXPR = QT*QT - PT*PT*PT
      IF(DABS(EXPR).LT.EPSILON) THEN
         TH = DACOS(DSIGN(ONE,-Q/2))
      ELSE
         ARG = -QT/(PT**THHF)
       ARG1 =  DABS(ARG) - ONE
       IF(DABS(ARG1) .LT. EPSILON) 
     1	ARG = (ONE - EPSILON)*DSIGN(ONE,ARG)
       TH = DACOS(ARG)
C         TH = DACOS((-QT)/((PT)**(THHF)))
      END IF
       SQH = TWO*DSQRT(PT)
      RKS(1) = H(1)+SQH*DCOS(TH/THREE)
      RKS(2) = H(1)+SQH*DCOS((TH+FOUR*PI)/THREE)
      RKS(3) = H(1)+SQH*DCOS((TH+TWO*PI)/THREE)
	DO I = 1 , 3
         IF (DABS(RKS(I)).LT.1D-16) THEN
            RKS(I) = ZERO
         END IF
      END DO
      RETURN
      END
! #############################################################################
      SUBROUTINE YLD3D(X1,X2,RKS1,RKS2,H1,H2,S,F,SIG0,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common/ki_select/i_hard,i_yield
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
      DIMENSION X1(6),X2(6),RKS1(3),RKS2(3),H1(3),H2(3),S(6)
       data zero,one,two,three,four /0.d0,1.d0,2.0d0,3.0d0,4.d0/
! RKS(3) EIGENVALUES OF Xij
! K YLD SECOND MEMBER
! INPUT ISOTROPIC YIELD FUNCTION (ALSO CHANGE IN RAT)                         
! YIELD FUNCTION               
C select largest rks value for normalization
	yn1 = max(abs(rks1(1)),abs(rks1(2)),abs(rks1(3)))
	yn2 = max(abs(rks2(1)),abs(rks2(2)),abs(rks2(3)))
	ynorm = max(yn1,yn2)
      DO I = 1 , 3
         RKS1(I) = RKS1(I)/ynorm
         RKS2(I) = RKS2(I)/ynorm
      END DO
	if(i_yield .eq. 3) then
	F = DABS(RKS1(1)-RKS1(2))**PEM - DABS(RKS1(1))**PEM                                    
     &  + DABS(RKS1(2)-RKS1(3))**PEM - DABS(RKS1(2))**PEM                                     
     &  + DABS(RKS1(3)-RKS1(1))**PEM - DABS(RKS1(3))**PEM 
     &  + DABS(RKS2(1))**PEM  + DABS(RKS2(2))**PEM
     &  + DABS(RKS2(3))**PEM                                     
C CT IS 1/SIG_0
       SIG0 = ynorm*(F/two)**(one/PEM)
	CT = ONE/SIG0
C	CT = (TWO/F)**(ONE/PEM)
	elseif(i_yield .eq. 4) then
	F = DABS(RKS1(1)-RKS1(2))**PEM 
     &  + DABS(RKS1(2)-RKS1(3))**PEM 
     &  + DABS(RKS1(3)-RKS1(1))**PEM 
     &  + DABS(RKS2(1))**PEM  + DABS(RKS2(2))**PEM
     &  + DABS(RKS2(3))**PEM                                     
C CT IS 1/SIG_0
	CX = (TWO+(TWO/THREE)**PEM+TWO*(ONE/THREE)**PEM)/TWO
C CT IS 1/SIG_0
       SIG0 = ynorm*(F/two/CX)**(one/PEM)
	CT = ONE/SIG0
C	CT = (TWO*CX/F)**(ONE/PEM)
	elseif(i_yield .eq. 5) then
	F = ZERO
	DO I = 1,3
	DO J = 1,3
c change on 04-29-09 to see if compatible with linux
c	F = F + DABS(RKS1(I) - RKS2(J))**PEM
	ARGF = DABS(RKS1(I) - RKS2(J))
        IF(ARGF .GT. 1.E-12) THEN
	F = F + ARGF**PEM
        ENDIF
	ENDDO
	ENDDO
C CT IS 1/SIG_0
       SIG0 = ynorm*(F/four)**(one/PEM)
	CT = ONE/SIG0
	endif
! SCALE VARIABLES TO SIZE OF YIELD FUNCTIONS
! PRINCIPAL TRANSFORM STRESSES
      DO I = 1 , 3
         RKS1(I) = ynorm*RKS1(I)*CT
         RKS2(I) = ynorm*RKS2(I)*CT
      END DO
! INVARIANTS
       CT2 = CT*CT
       CT3 = CT2*CT
	H1(1) = H1(1)*CT
      H1(2) = H1(2)*CT2
      H1(3) = H1(3)*CT3
	H2(1) = H2(1)*CT
      H2(2) = H2(2)*CT2
      H2(3) = H2(3)*CT3
! TRANSFORMED AND REAL STRESSES
      DO I = 1 , 6
         X1(I) = X1(I)*CT
         X2(I) = X2(I)*CT
         S(I)  = S(I)*CT
      END DO
      RETURN
      END
! #############################################################################
! FIRST & SECOND DERIVATIVE
! #############################################################################
      SUBROUTINE RAT3D(SIG0,X1,X2,RKS1,RKS2,H1,H2,
     &                 EXPR1,EXPR2,Q1,Q2,E,X2N,I2DER,mspt,iprint)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       common/yld2_1 /al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,
     1	al12,al13,al14,al15,al16,al17,al18,pem
      common/ki_select/i_hard,i_yield
      common/yld2_2/ RL(2,6,6),C(2,6,6),T(2,6,6)
      DIMENSION RL1(6,6),RL2(6,6),E(6),X2N(6,6)
      DIMENSION X1(6),X2(6),RKS1(3),RKS2(3),H1(3),H2(3)
      DIMENSION DFL1(3),DHX1(3,6),DLH1(3,3),DFH1(3)
      DIMENSION DFL2(3),DHX2(3,6),DLH2(3,3),DFH2(3)
       DIMENSION DF2L1(3,3),DF2L2(3,3),D2LH1(3,3,3),D2LH2(3,3,3)
       DIMENSION DF2L12(3,3),DF2L21(3,3)
       DIMENSION D2HX1(3,6,6),D2HX2(3,6,6)
       DIMENSION D1A(6,6),D2A(6,6),D3A(6,6)
       DIMENSION D1B(6,6),D2B(6,6),D3B(6,6)
       DIMENSION DD2A(6,6),DD2B(6,6)
       DIMENSION T13A(3),T3A(3,3),T36A(3,6),T63A(6,3),T63A2(6,3),
     1	T36A2(3,6),T6A(6,6),T6A2(6,6),T13A1(3),T13A2(3),T13A3(3),
     2	T36A12(3,6),T36AA(3,6),T3A12(3,3)
       DIMENSION T13B(3),T3B(3,3),T36B(3,6),T63B(6,3),T63B2(6,3),
     1	T36B2(3,6),T6B(6,6),T6B2(6,6),T13B1(3),T13B2(3),T13B3(3),
     2	T36B12(3,6),T36BB(3,6),T3B12(3,3)
      DIMENSION P(9),Y1(6),Y2(6)
      DATA EPSILON/1D-10/
       data zero,one,two,three,four /0.d0,1.d0,2.0d0,3.0d0,4.d0/
! INIT
       PEM1 = PEM - ONE
       PEM2 = PEM - TWO
	if(i_yield .eq. 4) 	CX = (TWO+(TWO/THREE)**PEM+
     1		TWO*(ONE/THREE)**PEM)/TWO

      DO I = 1 , 6
         Y1(I) = ZERO
         Y2(I) = ZERO
      END DO
! DERIVATIVES FOR FIRST LINEAR TRANSFORMATION 
! F = DABS(KS1(1)-KS1(2))**PEM - DABS(KS1(1))**PEM                                    
!   + DABS(KS1(2)-KS1(3))**PEM - DABS(KS1(2))**PEM                                     
!   + DABS(KS1(3)-KS1(1))**PEM - DABS(KS1(3))**PEM 
! DERIVATIVES DF/DL1i  
	IF(I_YIELD.EQ.3 .OR. I_YIELD.EQ.4) THEN
      P(1) = RKS1(2)-RKS1(3)                            
	P(2) = RKS1(3)-RKS1(1)                            
	P(3) = RKS1(1)-RKS1(2)
	P(4) = RKS1(1)
	P(5) = RKS1(2)
	P(6) = RKS1(3)
	ELSEIF(I_YIELD.EQ.5) THEN
	DO I = 1,3
	DO J = 1,3
	IJ = (I-1)*3 + J
	P(IJ) = RKS1(I) - RKS2(J)
	ENDDO
	ENDDO
	ENDIF
	if(i_yield .eq. 3) then
      DFL1(1) = (-DSIGN(ONE,P(2))*DABS(P(2))**PEM1
     &          + DSIGN(ONE,P(3))*DABS(P(3))**PEM1
     &          - DSIGN(ONE,P(4))*DABS(P(4))**PEM1)*PEM
      DFL1(2) = (-DSIGN(ONE,P(3))*DABS(P(3))**PEM1
     &          + DSIGN(ONE,P(1))*DABS(P(1))**PEM1
     &          - DSIGN(ONE,P(5))*DABS(P(5))**PEM1)*PEM
      DFL1(3) = (-DSIGN(ONE,P(1))*DABS(P(1))**PEM1
     &          + DSIGN(ONE,P(2))*DABS(P(2))**PEM1
     &          - DSIGN(ONE,P(6))*DABS(P(6))**PEM1)*PEM
	elseif(i_yield .eq. 4) then
      DFL1(1) = (-DSIGN(ONE,P(2))*DABS(P(2))**PEM1
     &          + DSIGN(ONE,P(3))*DABS(P(3))**PEM1)*PEM
      DFL1(2) = (-DSIGN(ONE,P(3))*DABS(P(3))**PEM1
     &          + DSIGN(ONE,P(1))*DABS(P(1))**PEM1)*PEM
      DFL1(3) = (-DSIGN(ONE,P(1))*DABS(P(1))**PEM1
     &          + DSIGN(ONE,P(2))*DABS(P(2))**PEM1)*PEM
	elseif(i_yield .eq. 5) then
      DFL1(1) = (DSIGN(ONE,P(1))*DABS(P(1))**PEM1
     &          + DSIGN(ONE,P(2))*DABS(P(2))**PEM1
     &          + DSIGN(ONE,P(3))*DABS(P(3))**PEM1)*PEM
      DFL1(2) = (DSIGN(ONE,P(4))*DABS(P(4))**PEM1
     &          + DSIGN(ONE,P(5))*DABS(P(5))**PEM1
     &          + DSIGN(ONE,P(6))*DABS(P(6))**PEM1)*PEM
      DFL1(3) = (DSIGN(ONE,P(7))*DABS(P(7))**PEM1
     &          + DSIGN(ONE,P(8))*DABS(P(8))**PEM1
     &          + DSIGN(ONE,P(9))*DABS(P(9))**PEM1)*PEM
	endif
! DERIVATIVE DH1i/dX1kl CALCULATION
      DHX1(1,1) = ONE/THREE
      DHX1(1,2) = ONE/THREE
      DHX1(1,3) = ONE/THREE
      DHX1(1,4) = ZERO
      DHX1(1,5) = ZERO
      DHX1(1,6) = ZERO
      DHX1(2,1) = -(X1(2)+X1(3))/THREE
      DHX1(2,2) = -(X1(3)+X1(1))/THREE
      DHX1(2,3) = -(X1(1)+X1(2))/THREE
      DHX1(2,4) = TWO*X1(4)/THREE
      DHX1(2,5) = TWO*X1(5)/THREE
      DHX1(2,6) = TWO*X1(6)/THREE
C FOR FULL 3D
      DHX1(3,1) = (X1(2)*X1(3)-X1(6)**TWO)/TWO
      DHX1(3,2) = (X1(3)*X1(1)-X1(5)**TWO)/TWO
      DHX1(3,3) = (X1(1)*X1(2)-X1(4)**TWO)/TWO
C PLANE STRESS
C      DHX1(3,1) = X1(2)*X1(3)/TWO
C      DHX1(3,2) = X1(3)*X1(1)/TWO
C      DHX1(3,3) = (X1(1)*X1(2)-X1(4)**TWO)/TWO
C PLANE STRESS X_z*X_xy
C      DHX1(3,4) = -X1(3)*X1(4)
C** CAUTION !!! MAKE SURE X1(4)=SXY,X1(5)=SZX,X1(6)=SYZ
      DHX1(3,4) = X1(5)*X1(6)-X1(3)*X1(4)
      DHX1(3,5) = X1(6)*X1(4)-X1(2)*X1(5)
      DHX1(3,6) = X1(4)*X1(5)-X1(1)*X1(6)
      IF(DABS(EXPR1).GT.EPSILON) THEN
       ! DERIVATIVES dL1k/dH1j
         DLH1(1,1)=(RKS1(1)**TWO)/(RKS1(1)**TWO-TWO*RKS1(1)*H1(1)-H1(2))     
         DLH1(2,1)=(RKS1(2)**TWO)/(RKS1(2)**TWO-TWO*RKS1(2)*H1(1)-H1(2)) 
         DLH1(3,1)=(RKS1(3)**TWO)/(RKS1(3)**TWO-TWO*RKS1(3)*H1(1)-H1(2)) 
         DLH1(1,2) =      RKS1(1)/(RKS1(1)**TWO-TWO*RKS1(1)*H1(1)-H1(2))
         DLH1(2,2) =      RKS1(2)/(RKS1(2)**TWO-TWO*RKS1(2)*H1(1)-H1(2))
         DLH1(3,2) =      RKS1(3)/(RKS1(3)**TWO-TWO*RKS1(3)*H1(1)-H1(2))
         DLH1(1,3) =    TWO/(RKS1(1)**TWO-TWO*RKS1(1)*H1(1)-H1(2))/THREE
         DLH1(2,3) =    TWO/(RKS1(2)**TWO-TWO*RKS1(2)*H1(1)-H1(2))/THREE
         DLH1(3,3) =    TWO/(RKS1(3)**TWO-TWO*RKS1(3)*H1(1)-H1(2))/THREE       
       ! DERIVATIVES dF/dH1j
         DO I = 1 , 3
            DFH1(I) = ZERO
            DO J = 1 , 3
               DFH1(I) = DFH1(I) + DFL1(J)*DLH1(J,I)
            END DO
         END DO
      ELSE
         IF(-Q1.GT.ZERO) THEN		! STRESS CASE THETA1 = 0 
         ! DERIVATIVES dL11/dH1i
       DLH1(1,1)=(RKS1(1)**TWO)/(RKS1(1)**TWO-TWO*RKS1(1)*H1(1)-H1(2))
       DLH1(1,2) =   RKS1(1)/(RKS1(1)**TWO-TWO*RKS1(1)*H1(1)-H1(2))
       DLH1(1,3) =   TWO/(RKS1(1)**TWO-TWO*RKS1(1)*H1(1)-H1(2))/THREE
         ! DERIVATIVES dF/dH1i
            DFH1(1) = (DFL1(1)-DFL1(2))*DLH1(1,1)+3*DFL1(2)
            DFH1(2) = (DFL1(1)-DFL1(2))*DLH1(1,2)
            DFH1(3) = (DFL1(1)-DFL1(2))*DLH1(1,3)
         ELSE					! STRESS CASE THETA1 = 180
         ! DERIVATIVES DL13/dH1i
         DLH1(3,1)=(RKS1(3)**TWO)/(RKS1(3)**TWO-TWO*RKS1(3)*H1(1)-H1(2))
         DLH1(3,2)=    RKS1(3)/(RKS1(3)**TWO-TWO*RKS1(3)*H1(1)-H1(2))
         DLH1(3,3)= TWO/(RKS1(3)**TWO-TWO*RKS1(3)*H1(1)-H1(2))/THREE
         ! DERIVATIVES dF/dH1i
            DFH1(1) = (DFL1(3)-DFL1(2))*DLH1(3,1)+THREE*DFL1(2)
            DFH1(2) = (DFL1(3)-DFL1(2))*DLH1(3,2)
            DFH1(3) = (DFL1(3)-DFL1(2))*DLH1(3,3)
         END IF
      END IF
      ! DERIVATIVES dF1/dX1
      DO I = 1 , 6
         Y1(I)=DFH1(1)*DHX1(1,I)+DFH1(2)*DHX1(2,I)+DFH1(3)*DHX1(3,I)
C         J = I + 3
C         Y1(J)=(DFH1(1)*DHX1(1,J)+DFH1(2)*DHX1(2,J)+DFH1(3)*DHX1(3,J))/2
      END DO
C       Y1(4)=(DFH1(1)*DHX1(1,4)+DFH1(2)*DHX1(2,4)+DFH1(3)*DHX1(3,4))
	IF(I_YIELD.EQ.3 .OR. I_YIELD.EQ.4) THEN
! DERIVATIVES FOR SECOND LINEAR TRANSFORMATION 
!  F =  DABS(KS2(1))**PEM  + DABS(KS2(2))**PEM  + DABS(KS2(3))**PEM                                     
	P(7) = RKS2(1)
	P(8) = RKS2(2)
	P(9) = RKS2(3)
	F2 =  DABS(P(7))**PEM + DABS(P(8))**PEM + DABS(P(9))**PEM
      DFL2(1) = PEM*DSIGN(ONE,P(7))*DABS(P(7))**PEM1
      DFL2(2) = PEM*DSIGN(ONE,P(8))*DABS(P(8))**PEM1
      DFL2(3) = PEM*DSIGN(ONE,P(9))*DABS(P(9))**PEM1
	elseif(i_yield .eq. 5) then
      DFL2(1) = -(DSIGN(ONE,P(1))*DABS(P(1))**PEM1
     &          + DSIGN(ONE,P(4))*DABS(P(4))**PEM1
     &          + DSIGN(ONE,P(7))*DABS(P(7))**PEM1)*PEM
      DFL2(2) = -(DSIGN(ONE,P(2))*DABS(P(2))**PEM1
     &          + DSIGN(ONE,P(5))*DABS(P(5))**PEM1
     &          + DSIGN(ONE,P(8))*DABS(P(8))**PEM1)*PEM
      DFL2(3) = -(DSIGN(ONE,P(3))*DABS(P(3))**PEM1
     &          + DSIGN(ONE,P(6))*DABS(P(6))**PEM1
     &          + DSIGN(ONE,P(9))*DABS(P(9))**PEM1)*PEM
	endif
! DERIVATIVE DH1i/dX1kl CALCULATION
      DHX2(1,1) = ONE/THREE
      DHX2(1,2) = ONE/THREE
      DHX2(1,3) = ONE/THREE
      DHX2(1,4) = ZERO
      DHX2(1,5) = ZERO
      DHX2(1,6) = ZERO
      DHX2(2,1) = -(X2(2)+X2(3))/THREE
      DHX2(2,2) = -(X2(3)+X2(1))/THREE
      DHX2(2,3) = -(X2(1)+X2(2))/THREE
      DHX2(2,4) = TWO*X2(4)/THREE
      DHX2(2,5) = TWO*X2(5)/THREE
      DHX2(2,6) = TWO*X2(6)/THREE
C FULL 3D
      DHX2(3,1) = (X2(2)*X2(3)-X2(6)**TWO)/2
      DHX2(3,2) = (X2(3)*X2(1)-X2(5)**TWO)/2
      DHX2(3,3) = (X2(1)*X2(2)-X2(4)**TWO)/2
C PLANE STRESS
C      DHX2(3,1) = X2(2)*X2(3)/TWO
C      DHX2(3,2) = X2(3)*X2(1)/TWO
C      DHX2(3,3) = (X2(1)*X2(2)-X2(4)**TWO)/TWO
C      DHX2(3,4) = -X2(3)*X2(4)
C** CAUTION !!! MAKE SURE X1(4)=SXY,X1(5)=SXZ,X1(6)=SYZ
      DHX2(3,4) = (X2(5)*X2(6)-X2(3)*X2(4))
      DHX2(3,5) = (X2(6)*X2(4)-X2(2)*X2(5))
      DHX2(3,6) = (X2(4)*X2(5)-X2(1)*X2(6))
      IF(DABS(EXPR2).GT.EPSILON) THEN
       ! DERIVATIVES dL1k/dH1j
       DLH2(1,1) = (RKS2(1)**TWO)/(RKS2(1)**TWO-TWO*RKS2(1)*H2(1)-H2(2)) 
       DLH2(2,1) = (RKS2(2)**TWO)/(RKS2(2)**TWO-TWO*RKS2(2)*H2(1)-H2(2)) 
       DLH2(3,1) = (RKS2(3)**TWO)/(RKS2(3)**TWO-TWO*RKS2(3)*H2(1)-H2(2)) 
       DLH2(1,2) =      RKS2(1)/(RKS2(1)**TWO-TWO*RKS2(1)*H2(1)-H2(2))
       DLH2(2,2) =      RKS2(2)/(RKS2(2)**TWO-TWO*RKS2(2)*H2(1)-H2(2))
       DLH2(3,2) =      RKS2(3)/(RKS2(3)**TWO-TWO*RKS2(3)*H2(1)-H2(2))
       DLH2(1,3) =      TWO/(RKS2(1)**TWO-TWO*RKS2(1)*H2(1)-H2(2))/THREE
       DLH2(2,3) =      TWO/(RKS2(2)**TWO-TWO*RKS2(2)*H2(1)-H2(2))/THREE
       DLH2(3,3) =      TWO/(RKS2(3)**TWO-TWO*RKS2(3)*H2(1)-H2(2))/THREE
       ! DERIVATIVES dF/dH1j
         DO I = 1 , 3
            DFH2(I) = ZERO
            DO J = 1 , 3
               DFH2(I) = DFH2(I) + DFL2(J)*DLH2(J,I)
            END DO
         END DO
      ELSE
         IF(-Q2.GT.ZERO) THEN       ! STRESS CASE THETA1 = 0 
C         ! DERIVATIVES dL11/dH1i
         DLH2(1,1)=(RKS2(1)**TWO)/(RKS2(1)**TWO-TWO*RKS2(1)*H2(1)-H2(2))
         DLH2(1,2) =    RKS2(1)/(RKS2(1)**TWO-TWO*RKS2(1)*H2(1)-H2(2))
         DLH2(1,3) =    TWO/(RKS2(1)**TWO-TWO*RKS2(1)*H2(1)-H2(2))/THREE       
C         ! DERIVATIVES dF/dH1i
            DFH2(1) = (DFL2(1)-DFL2(2))*DLH2(1,1)+THREE*DFL2(2)
            DFH2(2) = (DFL2(1)-DFL2(2))*DLH2(1,2)
            DFH2(3) = (DFL2(1)-DFL2(2))*DLH2(1,3)
         ELSE                                ! STRESS CASE THETA1 = 180
C         ! DERIVATIVES DL13/dH1i
         DLH2(3,1)=(RKS2(3)**TWO)/(RKS2(3)**TWO-TWO*RKS2(3)*H2(1)-H2(2))
            DLH2(3,2) = RKS2(3)/(RKS2(3)**TWO-TWO*RKS2(3)*H2(1)-H2(2))
            DLH2(3,3) = TWO/(RKS2(3)**TWO-TWO*RKS2(3)*H2(1)-H2(2))/THREE
C         ! DERIVATIVES dF/dH1i
            DFH2(1) = (DFL2(3)-DFL2(2))*DLH2(3,1)+THREE*DFL2(2)
            DFH2(2) = (DFL2(3)-DFL2(2))*DLH2(3,2)
            DFH2(3) = (DFL2(3)-DFL2(2))*DLH2(3,3)
         END IF
      END IF
C      ! DERIVATIVES dF1/dX1
      DO I = 1 , 6
         Y2(I)=DFH2(1)*DHX2(1,I)+DFH2(2)*DHX2(2,I)+DFH2(3)*DHX2(3,I)
C         J = I + 3
C         Y2(J)=(DFH2(1)*DHX2(1,J)+DFH2(2)*DHX2(2,J)+DFH2(3)*DHX2(3,J))/2
      END DO
C         Y2(4)=(DFH2(1)*DHX2(1,4)+DFH2(2)*DHX2(2,4)+DFH2(3)*DHX2(3,4))
! CALCULATE STRAIN RATES
	if(i_yield .eq. 3) then
       TWOM = TWO*PEM
	elseif(i_yield .eq. 4) then
       TWOM = TWO*CX*PEM
	elseif(i_yield .eq. 5) then
       TWOM = FOUR*PEM
	endif
      DO J = 1 , 6
         E(J) = ZERO
         DO I = 1 , 6
            E(J) = E(J) + Y1(I)*RL(1,I,J) + Y2(I)*RL(2,I,J)
         END DO
C CONVERT D(PHI)/DSIGIJ TO D(SIGO)/DSIGIJ
	   E(J) = E(J)/TWOM
      END DO
c       endif
       IF(I2DER .EQ. 0) RETURN
C CALCULATE SECOND DERIVATIVES
C 
       DO I = 1,3
       DO J = 1,3
       DF2L1(I,J) = ZERO
       DF2L2(I,J) = ZERO
       DF2L12(I,J) = ZERO
       DF2L21(I,J) = ZERO
              DO K = 1,3
              D2LH1(I,J,K) = ZERO
              D2LH2(I,J,K) = ZERO
       	      ENDDO
       ENDDO
       ENDDO
       DO I = 1,3
       DO J = 1,6
       DO K = 1,6
       D2HX1(I,J,K) = ZERO
       D2HX2(I,J,K) = ZERO
       ENDDO
       ENDDO
       ENDDO
       DO J = 1,6
       DO K = 1,6
       RL1(J,K) = RL(1,J,K)
       RL2(J,K) = RL(2,J,K)
       ENDDO
       ENDDO
C CALCULATE D2(PHI)/D(XPRIN)**2
C
	PEMM1 = PEM*PEM1
	if(i_yield .eq. 3) then
       DF2L1(1,1) =  PEMM1*(DABS(P(2))**PEM2+DABS(P(3))**PEM2-
     1			DABS(P(4))**PEM2)
       DF2L1(1,2) = -PEMM1*DABS(P(3))**PEM2
       DF2L1(1,3) = -PEMM1*DABS(P(2))**PEM2
       DF2L1(2,2) =  PEMM1*(DABS(P(1))**PEM2+DABS(P(3))**PEM2-
     1			DABS(P(5))**PEM2)
       DF2L1(2,3) = -PEMM1*DABS(P(1))**PEM2
       DF2L1(3,3) =  PEMM1*(DABS(P(1))**PEM2+DABS(P(2))**PEM2-
     1			DABS(P(6))**PEM2)
	elseif(i_yield .eq. 4) then
       DF2L1(1,1) =  PEMM1*(DABS(P(2))**PEM2+DABS(P(3))**PEM2)
       DF2L1(1,2) = -PEMM1*DABS(P(3))**PEM2
       DF2L1(1,3) = -PEMM1*DABS(P(2))**PEM2
       DF2L1(2,2) =  PEMM1*(DABS(P(1))**PEM2+DABS(P(3))**PEM2)
       DF2L1(2,3) = -PEMM1*DABS(P(1))**PEM2
       DF2L1(3,3) =  PEMM1*(DABS(P(1))**PEM2+DABS(P(2))**PEM2)
	elseif(i_yield .eq. 5) then
c change on 4-29-09 to accomodate raising zero to a power
c       DF2L1(1,1) =  PEMM1*(DABS(P(1))**PEM2+DABS(P(2))**PEM2+
c     1			DABS(P(3))**PEM2)
          darg1 = zero
          if(DABS(P(1)) .gt. ztol) darg1 = darg1 + DABS(P(1))**pem2
          if(DABS(P(2)) .gt. ztol) darg1 = darg1 + DABS(P(2))**pem2
          if(DABS(P(3)) .gt. ztol) darg1 = darg1 + DABS(P(3))**pem2
          df2l1(1,1) = pemm1*darg1
c       endif
c       DF2L1(2,2) =  PEMM1*(DABS(P(4))**PEM2+DABS(P(5))**PEM2+
c     1			DABS(P(6))**PEM2)
          darg1 = zero
          if(DABS(P(4)) .gt. ztol) darg1 = darg1 + DABS(P(4))**pem2
          if(DABS(P(5)) .gt. ztol) darg1 = darg1 + DABS(P(5))**pem2
          if(DABS(P(6)) .gt. ztol) darg1 = darg1 + DABS(P(6))**pem2
          df2l1(2,2) = pemm1*darg1
c       DF2L1(3,3) =  PEMM1*(DABS(P(7))**PEM2+DABS(P(8))**PEM2+
c     1			DABS(P(9))**PEM2)
          darg1 = zero
          if(DABS(P(7)) .gt. ztol) darg1 = darg1 + DABS(P(7))**pem2
          if(DABS(P(8)) .gt. ztol) darg1 = darg1 + DABS(P(8))**pem2
          if(DABS(P(9)) .gt. ztol) darg1 = darg1 + DABS(P(9))**pem2
          df2l1(3,3) = pemm1*darg1
	endif
	if(i_yield.eq.3 .OR. I_YIELD.EQ.4) then
       DF2L1(2,1) = DF2L1(1,2)
       DF2L1(3,1) = DF2L1(1,3)
       DF2L1(3,2) = DF2L1(2,3)
       DF2L2(1,1) = PEMM1*DABS(P(7))**PEM2
       DF2L2(2,2) = PEMM1*DABS(P(8))**PEM2
       DF2L2(3,3) = PEMM1*DABS(P(9))**PEM2
	elseif(i_yield .eq. 5) then
c change on 4-29-09 to accomodate raising zero to a power
c       DF2L2(1,1) =  PEMM1*(DABS(P(1))**PEM2+DABS(P(4))**PEM2+
c     1			DABS(P(7))**PEM2)
          darg1 = zero
          if(DABS(P(1)) .gt. ztol) darg1 = darg1 + DABS(P(1))**pem2
          if(DABS(P(4)) .gt. ztol) darg1 = darg1 + DABS(P(4))**pem2
          if(DABS(P(7)) .gt. ztol) darg1 = darg1 + DABS(P(7))**pem2
          df2l2(1,1) = pemm1*darg1
c       DF2L2(2,2) =  PEMM1*(DABS(P(2))**PEM2+DABS(P(5))**PEM2+
c     1			DABS(P(8))**PEM2)
          darg1 = zero
          if(DABS(P(2)) .gt. ztol) darg1 = darg1 + DABS(P(2))**pem2
          if(DABS(P(5)) .gt. ztol) darg1 = darg1 + DABS(P(5))**pem2
          if(DABS(P(8)) .gt. ztol) darg1 = darg1 + DABS(P(8))**pem2
          df2l2(2,2) = pemm1*darg1
c       DF2L2(3,3) =  PEMM1*(DABS(P(3))**PEM2+DABS(P(6))**PEM2+
c     1			DABS(P(9))**PEM2)
          darg1 = zero
          if(DABS(P(3)) .gt. ztol) darg1 = darg1 + DABS(P(3))**pem2
          if(DABS(P(6)) .gt. ztol) darg1 = darg1 + DABS(P(6))**pem2
          if(DABS(P(9)) .gt. ztol) darg1 = darg1 + DABS(P(9))**pem2
          df2l2(3,3) = pemm1*darg1
	DO I = 1,3
	DO J = 1,3
C CROSS DERIVATIVE WRT S1(I) AND THEN S2(J)
	IJ = (I-1)*3 + J
c change on 4-29-09 to accomodate raising zero to a power
       if(iprint .eq. 1) then
             OPEN(98,FILE='/usr2/karabin/frct/dbg.sav',
     1	STATUS='OLD',ACCESS='APPEND')
        write(98,1085) ij,P(ij)
 1085   format('ij,P'i4,e12.5)
       CLOSE(98)
       endif
        if (dabs(P(ij)) .gt. ztol) then
	DF2L12(I,J) = -PEMM1*DABS(P(IJ))**PEM2
        else
	DF2L12(I,J) = zero
        endif
	ENDDO
	ENDDO
	DO I = 1,3
	DO J = 1,3
C CROSS DERIVATIVE WRT S2(I) AND THEN S1(J)
	DF2L21(I,J) = DF2L12(J,I)
	ENDDO
	ENDDO
	endif
C COMPUTE D2(H)/DXab**2
       OTM = -ONE/THREE
       TTP =  TWO/THREE
       D2HX1(2,1,2) = OTM
       D2HX1(2,1,3) = OTM
       D2HX1(2,2,3) = OTM
       D2HX1(2,2,1) = D2HX1(2,1,2)
       D2HX1(2,3,1) = D2HX1(2,1,3)
       D2HX1(2,3,2) = D2HX1(2,2,3)
       D2HX1(2,4,4) = TTP
       D2HX1(2,5,5) = TTP
       D2HX1(2,6,6) = TTP
       D2HX1(3,1,2) = X1(3)/TWO
       D2HX1(3,1,3) = X1(2)/TWO
       D2HX1(3,1,6) = -X1(6)
       D2HX1(3,2,3) = X1(1)/TWO
       D2HX1(3,2,5) = -X1(5)
       D2HX1(3,3,4) = -X1(4)
       D2HX1(3,4,4) = -X1(3)
       D2HX1(3,4,5) = X1(6)
       D2HX1(3,4,6) = X1(5)
       D2HX1(3,5,5) = -X1(2)
       D2HX1(3,5,6) = X1(4)
       D2HX1(3,6,6) = -X1(1)
       D2HX1(3,2,1) = D2HX1(3,1,2)
       D2HX1(3,3,1) = D2HX1(3,1,3)
       D2HX1(3,6,1) = D2HX1(3,1,6)
       D2HX1(3,3,2) = D2HX1(3,2,3)
       D2HX1(3,5,2) = D2HX1(3,2,5)
       D2HX1(3,4,3) = D2HX1(3,3,4)
       D2HX1(3,5,4) = D2HX1(3,4,5)
       D2HX1(3,6,4) = D2HX1(3,4,6)
       D2HX1(3,6,5) = D2HX1(3,5,6)
C FOR X2 VARIABLES
       D2HX2(2,1,2) = D2HX1(2,1,2)
       D2HX2(2,1,3) = D2HX1(2,1,3)
       D2HX2(2,2,3) = D2HX1(2,2,3)
       D2HX2(2,2,1) = D2HX1(2,1,2)
       D2HX2(2,3,1) = D2HX1(2,1,3)
       D2HX2(2,3,2) = D2HX1(2,2,3)
       D2HX2(2,4,4) = TTP
       D2HX2(2,5,5) = TTP
       D2HX2(2,6,6) = TTP
       D2HX2(3,1,2) = X2(3)/TWO
       D2HX2(3,1,3) = X2(2)/TWO
       D2HX2(3,1,6) = -X2(6)
       D2HX2(3,2,3) = X2(1)/TWO
       D2HX2(3,2,5) = -X2(5)
       D2HX2(3,3,4) = -X2(4)
       D2HX2(3,4,4) = -X2(3)
       D2HX2(3,4,5) = X2(6)
       D2HX2(3,4,6) = X2(5)
       D2HX2(3,5,5) = -X2(2)
       D2HX2(3,5,6) = X2(4)
       D2HX2(3,6,6) = -X2(1)
       D2HX2(3,2,1) = D2HX2(3,1,2)
       D2HX2(3,3,1) = D2HX2(3,1,3)
       D2HX2(3,6,1) = D2HX2(3,1,6)
       D2HX2(3,3,2) = D2HX2(3,2,3)
       D2HX2(3,5,2) = D2HX2(3,2,5)
       D2HX2(3,4,3) = D2HX2(3,3,4)
       D2HX2(3,5,4) = D2HX2(3,4,5)
       D2HX2(3,6,4) = D2HX2(3,4,6)
       D2HX2(3,6,5) = D2HX2(3,5,6)
       IBEG = 1
       IEND = 3
      IF(DABS(EXPR1).LT.EPSILON) THEN
         IF(-Q1.GT.ZERO) THEN		! STRESS CASE THETA1 = 0 
              IEND = 1
         ELSE				! STRESS CASE THETA1 = 180
       		IBEG = 3
         ENDIF
       ENDIF
C FIRST FIND DERIVATIVE OF D2(XPRIN')/DH**2 FOR STANDARD CASE
       DO I = IBEG,IEND
       DENOM1 = RKS1(I)*RKS1(I)-TWO*RKS1(I)*H1(1)-H1(2)
       D1S = DENOM1*DENOM1
       R12 = RKS1(I)*RKS1(I)
       R13 = R12*RKS1(I)
       D2LH1(I,1,1) = (FOUR*R13+TWO*H1(1)*R12*R12/DENOM1-
     1			TWO*R12*R13/DENOM1)/D1S
       D2LH1(I,1,2) = (THREE*R12+TWO*H1(1)*R13/DENOM1-
     1			TWO*R12*R12/DENOM1)/D1S
       D2LH1(I,1,3) = FOUR*(RKS1(I)+H1(1)*R12/DENOM1-
     1			R13/DENOM1)/(THREE*D1S)
       D2LH1(I,2,2) = TWO*(RKS1(I)+H1(1)*R12/DENOM1-
     1			R13/DENOM1)/D1S
       D2LH1(I,2,3) = (TWO+FOUR*H1(1)*RKS1(I)/DENOM1-
     1			FOUR*R12/DENOM1)/(THREE*D1S)
       D2LH1(I,3,3) = FOUR*TWO*(H1(1)-RKS1(I))/(THREE*DENOM1*THREE*D1S)
       D2LH1(I,2,1) = D2LH1(I,1,2)
       D2LH1(I,3,1) = D2LH1(I,1,3)
       D2LH1(I,3,2) = D2LH1(I,2,3)
       ENDDO
       IBEG = 1
       IEND = 3
      IF(DABS(EXPR2).LT.EPSILON) THEN
         IF(-Q2.GT.ZERO) THEN		! STRESS CASE THETA1 = 0 
              IEND = 1
         ELSE				! STRESS CASE THETA1 = 180
       		IBEG = 3
         ENDIF
       ENDIF
C FIRST FIND DERIVATIVE OF D2(XPRIN")/DH**2 FOR STANDARD CASE
       DO I = IBEG,IEND
       DENOM2 = RKS2(I)*RKS2(I)-TWO*RKS2(I)*H2(1)-H2(2)
       D2S = DENOM2*DENOM2
       R22 = RKS2(I)*RKS2(I)
       R23 = R22*RKS2(I)
       D2LH2(I,1,1) = (FOUR*R23+TWO*H2(1)*R22*R22/DENOM2-
     1			TWO*R22*R23/DENOM2)/D2S
       D2LH2(I,1,2) = (THREE*R22+TWO*H2(1)*R23/DENOM2-
     1			TWO*R22*R22/DENOM2)/D2S
       D2LH2(I,1,3) = FOUR*(RKS2(I)+H2(1)*R22/DENOM2-
     1			R23/DENOM2)/(THREE*D2S)
       D2LH2(I,2,2) = TWO*(RKS2(I)+H2(1)*R22/DENOM2-
     1			R23/DENOM2)/D2S
       D2LH2(I,2,3) = (TWO+FOUR*H2(1)*RKS2(I)/DENOM2-
     1			FOUR*R22/DENOM2)/(THREE*D2S)
       D2LH2(I,3,3) = FOUR*TWO*(H2(1)-RKS2(I))/(THREE*DENOM2*THREE*D2S)
       D2LH2(I,2,1) = D2LH2(I,1,2)
       D2LH2(I,3,1) = D2LH2(I,1,3)
       D2LH2(I,3,2) = D2LH2(I,2,3)
       ENDDO
      IF(DABS(EXPR1).LT.EPSILON) GO TO 100
C PERFORM MATRIX MULTIPLICATION TO COMPUTE D2(PHI')/DSIG**2 IN 
C NONDIMENSIONAL FORM
C FIRST TERM HAS D2(PHI)/D(XPRIN)**2
C 3'X3 TIME 3'X3 MATRICES-ROW PRODUCT
       CALL MATMUL3(DLH1,DF2L1,T3A)
C 3X3' (T3A) TIMES A 3'X6-INNER PRODUCT
       CALL MATMUL36(T3A,DHX1,T36A)
C 3X6' (T36A) TIMES A 6'X6(RL1)-INNER PRODUCT
       CALL MATMUL36I(T36A,RL1,T36A2)
	IF(I_YIELD .EQ. 5) THEN
C NEED ADDITIONAL TERMS FOR CROSS DERIVATIVES WRT S1 AND S2
C 3'X3 TIME 3'X3 MATRICES-ROW PRODUCT
       CALL MATMUL3(DLH2,DF2L12,T3A)
C 3X3' (T3A) TIMES A 3'X6-INNER PRODUCT
       CALL MATMUL36(T3A,DHX2,T36A)
C 3X6' (T36A) TIMES A 6'X6(RL2)-INNER PRODUCT
       CALL MATMUL36I(T36A,RL2,T36AA)
C ADD TWO MATRICES BEFORE PROCEEDING WITH THE CHAIN RULE
	DO I = 1,3
	DO J = 1,6
	T36A2(I,J) = T36A2(I,J) + T36AA(I,J)
	ENDDO
	ENDDO
	ENDIF	
C 3'X6 TIMES A 3'X3 (D(XPRIN)/DH)-ROW PRODUCT
       CALL MATMUL63(T36A2,DLH1,T63A)
C 6X3'  TIMES A 3'X6 (DH/DX)-INNER PRODUCT
       CALL MATMUL66(T63A,DHX1,T6A)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6A,RL1,D1A)
C SECOND TERM HAS D2(XPRIN)/DH**2
C 1X3' (DPHI/DXPRIN) X 3'X3X3 (D2(XPRIN)/DHDH) MATRICES-PRIMES DENOTE MULTIPLIED
C INDICES
       CALL MATMUL33(DFL1,D2LH1,T3A)
C 3X3' (T3A) TIMES A 3'X6(DH/DX)-INNER PRODUCT
       CALL MATMUL36(T3A,DHX1,T36A)
C 3X6' (T36A)TIMES  (DX/DSIG) A 6'X6-INNER PRODUCT
       CALL MATMUL36I(T36A,RL1,T36A2)
C 3'X6  TIMES A 3'X6 (DH/DX)-ROW PRODUCT
       CALL MATMUL36R(T36A2,DHX1,T6A)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6A,RL1,D2A)
C THIRD TERM HAS D2H/DSIG**2
C 1X3' (DPHI/DXPRIN) TIMES 3'X3 (D(XPRIN)/DH) MATRICES
       CALL MATMUL313(DFL1,DLH1,T13A)
C 1X3' (DPHI/DXPRIN) X 3X6X6 (D2H/DX2) MATRICES
       CALL MATMUL336(T13A,D2HX1,T6A)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6A,RL1,T6A2)
C 6'X6  TIMES A 6'X6 (DX/DSIG)-ROW PRODUCT
       CALL MATMUL6R(T6A2,RL1,D3A)
C COMBINE ALL TERMS FOR SECOND DEVIATIVES
       DO I = 1,6
       DO J = 1,6
       DD2A(I,J) = D1A(I,J) + D2A(I,J) + D3A(I,J)
       ENDDO
       ENDDO
       GO TO 200
100         IF(-Q1.GT.ZERO) THEN		! STRESS CASE THETA1 = 0 
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = (DF2L1(1,1)-DF2L1(1,2))*DLH1(1,I) + DF2L1(1,2)*DT3
       T13A1(I) = DLH1(1,I)
       T13A2(I) = (DF2L1(1,2)-DF2L1(2,2))*DLH1(1,I) + DF2L1(2,2)*DT3
       T13A3(I) = DT3 - DLH1(1,I)
       ENDDO
C COMBINE TERMS TO CALCULATE D2(PHI)/DH2 A 3X3 MATRIX
       FACA = DFL1(1) - DFL1(2)       
       DO I = 1,3
       DO J = 1,3
       T3A(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J) + 
     1		FACA*D2LH1(1,I,J)       
       ENDDO
       ENDDO
		IF(I_YIELD .EQ. 5) THEN
			IF(DABS(EXPR2).LT.EPSILON) THEN
C FOR MULTIPLE ROOTS IN RKS2 ALSO
			IF(-Q2.GT.ZERO) THEN	! STRESS CASE THETA2 = 0 
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = (DF2L12(1,1)-DF2L12(1,2))*DLH2(1,I)+DF2L12(1,2)*DT3
       T13A1(I) = DLH1(1,I)
       T13A2(I) = (DF2L12(2,1)-DF2L12(2,2))*DLH2(1,I)+DF2L12(2,2)*DT3
       T13A3(I) = DT3 - DLH1(1,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3A12(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J)
       ENDDO
       ENDDO
		         ELSE			! STRESS CASE THETA2 = 180
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = (DF2L12(1,3)-DF2L12(1,2))*DLH2(3,I)+DF2L12(1,2)*DT3
       T13A1(I) = DLH1(1,I)
       T13A2(I) = (DF2L12(2,3)-DF2L12(2,2))*DLH2(3,I)+DF2L12(2,2)*DT3
       T13A3(I) = DT3 - DLH1(1,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3A12(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J)
       ENDDO
       ENDDO
			ENDIF
C ONLY MULTIPLE ROOTS IN RKS1
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = DF2L12(1,1)*DLH2(1,I)+DF2L12(1,2)*DLH2(2,I)
     1		+DF2L12(1,3)*DLH2(3,I)
       T13A1(I) = DLH1(1,I)
       T13A2(I) = DF2L12(2,1)*DLH2(1,I)+DF2L12(2,2)*DLH2(2,I)
     1		+DF2L12(2,3)*DLH2(3,I)
       T13A3(I) = DT3 - DLH1(1,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3A12(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J)
       ENDDO
       ENDDO
			ENDIF
		ENDIF
         ELSE					! STRESS CASE THETA1 = 180
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = (DF2L1(3,3)-DF2L1(3,2))*DLH1(3,I) + DF2L1(3,2)*DT3
       T13A1(I) = DLH1(3,I)
       T13A2(I) = (DF2L1(3,2)-DF2L1(2,2))*DLH1(3,I) + DF2L1(2,2)*DT3
       T13A3(I) = DT3 - DLH1(3,I)
       ENDDO
C COMBINE TERMS TO CALCULATE D2(PHI)/DH2 A 3X3 MATRIX
       FACA = DFL1(3) - DFL1(2)       
       DO I = 1,3
       DO J = 1,3
       T3A(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J) + 
     1		FACA*D2LH1(3,I,J)       
       ENDDO
       ENDDO
		IF(I_YIELD .EQ. 5) THEN
			IF(DABS(EXPR2).LT.EPSILON) THEN
C FOR MULTIPLE ROOTS IN RKS2 ALSO
			IF(-Q2.GT.ZERO) THEN	! STRESS CASE THETA2 = 0 
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = (DF2L12(3,1)-DF2L12(3,2))*DLH2(1,I)+DF2L12(3,2)*DT3
       T13A1(I) = DLH1(3,I)
       T13A2(I) = (DF2L12(2,1)-DF2L12(2,2))*DLH2(1,I)+DF2L12(2,2)*DT3
       T13A3(I) = DT3 - DLH1(3,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3A12(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J)
       ENDDO
       ENDDO
		         ELSE			! STRESS CASE THETA2 = 180
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = (DF2L12(3,3)-DF2L12(3,2))*DLH2(3,I)+DF2L12(3,2)*DT3
       T13A1(I) = DLH1(3,I)
       T13A2(I) = (DF2L12(2,3)-DF2L12(2,2))*DLH2(3,I)+DF2L12(2,2)*DT3
       T13A3(I) = DT3 - DLH1(3,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3A12(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J)
       ENDDO
       ENDDO
			ENDIF
C ONLY MULTIPLE ROOTS IN RKS1
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13A(I) = DF2L12(3,1)*DLH2(1,I)+DF2L12(3,2)*DLH2(2,I)
     1		+DF2L12(3,3)*DLH2(3,I)
       T13A1(I) = DLH1(3,I)
       T13A2(I) = DF2L12(2,1)*DLH2(1,I)+DF2L12(2,2)*DLH2(2,I)
     1		+DF2L12(2,3)*DLH2(3,I)
       T13A3(I) = DT3 - DLH1(3,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3A12(I,J) = T13A(I)*T13A1(J) + T13A2(I)*T13A3(J)
       ENDDO
       ENDDO
			ENDIF
		ENDIF
       ENDIF
C MULTIPLY MATRICES (CHAIN RULE) TO GET D2PHI/DH**2
C 3x3' (T3A) TIMES 3'X6 (DH/DX) INNER PRODUCT
       CALL MATMUL36(T3A,DHX1,T36A)
C 3X6' (T36A) TIMES  (DX/DSIG) A 6'X6
       CALL MATMUL36I(T36A,RL1,T36A2)
	IF(I_YIELD .EQ. 5) THEN
C NEED ADDITIONAL TERMS FOR CROSS DERIVATIVES WRT S1 AND S2
C 3x3' (T3A12) TIMES 3'X6 (DH/DX) INNER PRODUCT
       CALL MATMUL36(T3A12,DHX2,T36A)
C 3X6' (T36A) TIMES  (DX/DSIG) A 6'X6
       CALL MATMUL36I(T36A,RL2,T36AA)
C ADD TWO MATRICES BEFORE PROCEEDING WITH THE CHAIN RULE
	DO I = 1,3
	DO J = 1,6
	T36A2(I,J) = T36A2(I,J) + T36AA(I,J)
	ENDDO
	ENDDO
	ENDIF	
C 3'X6 (T36A2) TIMES A 3'X6 (DH/DX) ROW X ROW
       CALL MATMUL36R(T36A2,DHX1,T6A)
C 6X6'  TIMES A 6'X6 (DX/DSIG)
       CALL MATMUL6(T6A,RL1,D1A)
C MULTIPLY MATRICES (CHAIN RULE) TO GET D2H/DX**2
C 1X3' (DPHI/DH) X 3'X6X6 (D2H/DX2) MATRICES
       CALL MATMUL336(DFH1,D2HX1,T6A)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6A,RL1,T6A2)
C 6'X6  TIMES A 6'X6 (DX/DSIG)-ROW PRODUCT
       CALL MATMUL6R(T6A2,RL1,D3A)
C COMBINE ALL TERMS FOR SECOND DEVIATIVES
       DO I = 1,6
       DO J = 1,6
       DD2A(I,J) = D1A(I,J) + D3A(I,J)
       ENDDO
       ENDDO
200      IF(DABS(EXPR2).LT.EPSILON) GO TO 300
C PERFORM MATRIX MULTIPLICATION TO COMPUTE D2(PHI")/DSIG**2 IN 
C NONDIMENSIONAL FORM
C FIRST TERM HAS D2(PHI)/D(XPRIN)**2
C 3'X3 TIME 3'X3 MATRICES-ROW PRODUCT
       CALL MATMUL3(DLH2,DF2L2,T3B)
C 3X3' (T3B) TIMES A 3'X6(DHX2)-INNER PRODUCT
       CALL MATMUL36(T3B,DHX2,T36B)
C 3X6' (T36B) TIMES A 6'X6(RL2)-INNER PRODUCT
       CALL MATMUL36I(T36B,RL2,T36B2)
	IF(I_YIELD .EQ. 5) THEN
C NEED ADDITIONAL TERMS FOR CROSS DERIVATIVES WRT S2 AND S1
C 3'X3 TIME 3'X3 MATRICES-ROW PRODUCT
       CALL MATMUL3(DLH1,DF2L21,T3B)
C 3X3' (T3B) TIMES A 3'X6-INNER PRODUCT
       CALL MATMUL36(T3B,DHX1,T36B)
C 3X6' (T36B) TIMES A 6'X6(RL1)-INNER PRODUCT
       CALL MATMUL36I(T36B,RL1,T36BB)
C ADD TWO MATRICES BEFORE PROCEEDING WITH THE CHAIN RULE
	DO I = 1,3
	DO J = 1,6
	T36B2(I,J) = T36BB(I,J) + T36B2(I,J)
	ENDDO
	ENDDO
	ENDIF	
C 3'X6 TIMES A 3'X3 (D(XPRIN)/DH)-ROW PRODUCT
       CALL MATMUL63(T36B2,DLH2,T63B)
C 6X3'  TIMES A 3'X6 (DH/DX)-INNER PRODUCT
       CALL MATMUL66(T63B,DHX2,T6B)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6B,RL2,D1B)
C SECOND TERM HAS D2(XPRIN)/DH**2
C 1X3' (DPHI/DXPRIN) X 3'X3X3 (D2(XPRIN)/DHDH) MATRICES-PRIMES DENOTE MULTIPLIED
C INDICES
       CALL MATMUL33(DFL2,D2LH2,T3B)
C 3X3' (T3B) TIMES A 3'X6(DH/DX)-INNER PRODUCT
       CALL MATMUL36(T3B,DHX2,T36B)
C 3X6' (T36B)TIMES  (DX/DSIG) A 6'X6-INNER PRODUCT
       CALL MATMUL36I(T36B,RL2,T36B2)
C 3'X6  TIMES A 3'X6 (DH/DX)-ROW PRODUCT
       CALL MATMUL36R(T36B2,DHX2,T6B)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6B,RL2,D2B)
C THIRD TERM HAS D2H/DSIG**2
C 1X3' (DPHI/DXPRIN) TIMES 3'X3 (D(XPRIN)/DH) MATRICES
       CALL MATMUL313(DFL2,DLH2,T13B)
C 1X3' (DPHI/DXPRIN) X 3X6X6 (D2H/DX2) MATRICES
       CALL MATMUL336(T13B,D2HX2,T6B)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6B,RL2,T6B2)
C 6'X6  TIMES A 6'X6 (DX/DSIG)-ROW PRODUCT
       CALL MATMUL6R(T6B2,RL2,D3B)
C COMBINE ALL TERMS FOR SECOND DEVIATIVES
       DO I = 1,6
       DO J = 1,6
       DD2B(I,J) = D1B(I,J) + D2B(I,J) + D3B(I,J)
       ENDDO
       ENDDO
       GO TO 600
300         IF(-Q2.GT.ZERO) THEN		! STRESS CASE THETA2 = 0 
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = (DF2L2(1,1)-DF2L2(1,2))*DLH2(1,I) + DF2L2(1,2)*DT3
       T13B1(I) = DLH2(1,I)
       T13B2(I) = (DF2L2(1,2)-DF2L2(2,2))*DLH2(1,I) + DF2L2(2,2)*DT3
       T13B3(I) = DT3 - DLH2(1,I)
       ENDDO
C COMBINE TERMS TO CALCULATE D2(PHI)/DH2 A 3X3 MATRIX
       FACB = DFL2(1) - DFL2(2)       
       DO I = 1,3
       DO J = 1,3
       T3B(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J) + 
     1		FACB*D2LH2(1,I,J)       
       ENDDO
       ENDDO
		IF(I_YIELD .EQ. 5) THEN
			IF(DABS(EXPR1).LT.EPSILON) THEN
C FOR MULTIPLE ROOTS IN RKS1 ALSO
			IF(-Q2.GT.ZERO) THEN	! STRESS CASE THETA1 = 0 
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = (DF2L21(1,1)-DF2L21(1,2))*DLH1(1,I)+DF2L21(1,2)*DT3
       T13B1(I) = DLH2(1,I)
       T13B2(I) = (DF2L21(2,1)-DF2L21(2,2))*DLH1(1,I)+DF2L21(2,2)*DT3
       T13B3(I) = DT3 - DLH2(1,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3B12(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J)
       ENDDO
       ENDDO
		         ELSE			! STRESS CASE THETA1 = 180
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = (DF2L21(1,3)-DF2L21(1,2))*DLH1(3,I)+DF2L21(1,2)*DT3
       T13B1(I) = DLH2(1,I)
       T13B2(I) = (DF2L21(2,3)-DF2L21(2,2))*DLH1(3,I)+DF2L21(2,2)*DT3
       T13B3(I) = DT3 - DLH2(1,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3B12(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J)
       ENDDO
       ENDDO
			ENDIF
C ONLY MULTIPLE ROOTS IN RKS2
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = DF2L21(1,1)*DLH1(1,I)+DF2L21(1,2)*DLH1(2,I)
     1		+DF2L21(1,3)*DLH1(3,I)
       T13B1(I) = DLH2(1,I)
       T13B2(I) = DF2L21(2,1)*DLH1(1,I)+DF2L21(2,2)*DLH1(2,I)
     1		+DF2L21(2,3)*DLH1(3,I)
       T13B3(I) = DT3 - DLH2(1,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3B12(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J)
       ENDDO
       ENDDO
			ENDIF
		ENDIF

         ELSE					! STRESS CASE THETA2 = 180
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = (DF2L2(3,3)-DF2L2(3,2))*DLH2(3,I) + DF2L2(3,2)*DT3
       T13B1(I) = DLH2(3,I)
       T13B2(I) = (DF2L2(3,2)-DF2L2(2,2))*DLH2(3,I) + DF2L2(3,2)*DT3
       T13B3(I) = DT3 - DLH2(3,I)
       ENDDO
C COMBINE TERMS TO CALCULATE D2(PHI)/DH2 A 3X3 MATRIX
       FACB = DFL2(3) - DFL2(2)       
       DO I = 1,3
       DO J = 1,3
       T3B(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J) + 
     1		FACB*D2LH2(3,I,J)       
       ENDDO
       ENDDO
		IF(I_YIELD .EQ. 5) THEN
			IF(DABS(EXPR1).LT.EPSILON) THEN
C FOR MULTIPLE ROOTS IN RKS1 ALSO
			IF(-Q2.GT.ZERO) THEN	! STRESS CASE THETA1 = 0 
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = (DF2L21(3,1)-DF2L21(3,2))*DLH1(1,I)+DF2L21(3,2)*DT3
       T13B1(I) = DLH2(3,I)
       T13B2(I) = (DF2L21(2,1)-DF2L21(2,2))*DLH1(1,I)+DF2L21(2,2)*DT3
       T13B3(I) = DT3 - DLH2(3,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3B12(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J)
       ENDDO
       ENDDO
		         ELSE			! STRESS CASE THETA1 = 180
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = (DF2L21(3,3)-DF2L21(3,2))*DLH1(3,I)+DF2L21(3,2)*DT3
       T13B1(I) = DLH2(3,I)
       T13B2(I) = (DF2L21(2,3)-DF2L21(2,2))*DLH1(3,I)+DF2L21(2,2)*DT3
       T13B3(I) = DT3 - DLH2(3,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3B12(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J)
       ENDDO
       ENDDO
			ENDIF
C ONLY MULTIPLE ROOTS IN RKS2
       DO I = 1,3
       DT3 = ZERO
       IF(I .EQ. 1) DT3 = THREE
       T13B(I) = DF2L21(3,1)*DLH1(1,I)+DF2L21(3,2)*DLH1(2,I)
     1		+DF2L21(3,3)*DLH1(3,I)
       T13B1(I) = DLH2(3,I)
       T13B2(I) = DF2L21(2,1)*DLH1(1,I)+DF2L21(2,2)*DLH1(2,I)
     1		+DF2L21(2,3)*DLH1(3,I)
       T13B3(I) = DT3 - DLH2(3,I)
       ENDDO
       DO I = 1,3
       DO J = 1,3
       T3B12(I,J) = T13B(I)*T13B1(J) + T13B2(I)*T13B3(J)
       ENDDO
       ENDDO
			ENDIF
		ENDIF
       ENDIF
C MULTIPLY MATRICES (CHAIN RULE) TO GET D2PHI/DH**2
C 3x3' (T3B) TIMES 3'X6 (DH/DX) INNER PRODUCT
       CALL MATMUL36(T3B,DHX2,T36B)
C 3X6' (T36B) TIMES  (DX/DSIG) A 6'X6
       CALL MATMUL36I(T36B,RL2,T36B2)
	IF(I_YIELD .EQ. 5) THEN
C NEED ADDITIONAL TERMS FOR CROSS DERIVATIVES WRT S1 AND S2
C 3x3' (T3B12) TIMES 3'X6 (DH/DX) INNER PRODUCT
       CALL MATMUL36(T3B12,DHX1,T36B)
C 3X6' (T36B) TIMES  (DX/DSIG) A 6'X6
       CALL MATMUL36I(T36B,RL1,T36BB)
C ADD TWO MATRICES BEFORE PROCEEDING WITH THE CHAIN RULE
	DO I = 1,3
	DO J = 1,6
	T36B2(I,J) = T36B2(I,J) + T36BB(I,J)
	ENDDO
	ENDDO
	ENDIF	
C 3'X6 (T36B2) TIMES A 3'X6 (DH/DX) ROW X ROW
       CALL MATMUL36R(T36B2,DHX2,T6B)
C 6X6'  TIMES A 6'X6 (DX/DSIG)
       CALL MATMUL6(T6B,RL2,D1B)
C MULTIPLY MATRICES (CHAIN RULE) TO GET D2H/DX**2
C 1X3' (DPHI/DH) X 3'X6X6 (D2H/DX2) MATRICES
       CALL MATMUL336(DFH2,D2HX2,T6B)
C 6X6'  TIMES A 6'X6 (DX/DSIG)-INNER PRODUCT
       CALL MATMUL6(T6B,RL2,T6B2)
C 6'X6  TIMES A 6'X6 (DX/DSIG)-ROW PRODUCT
       CALL MATMUL6R(T6B2,RL2,D3B)
C COMBINE ALL TERMS FOR SECOND DEVIATIVES
       DO I = 1,6
       DO J = 1,6
       DD2B(I,J) = D1B(I,J) + D3B(I,J)
       ENDDO
       ENDDO
C600       DFAC1 = -PEM1/(FOUR*PEM*PEM*SIG0)
C NORMALS ALREADY HAVE 1/2M IN THEM
600       DFAC1 = -PEM1/SIG0
	if(i_yield .eq. 3) then
       DFAC2 = ONE/(TWO*PEM*SIG0)
	elseif(i_yield .eq. 4) then
       DFAC2 = ONE/(TWO*CX*PEM*SIG0)
	elseif(i_yield .eq. 5) then
       DFAC2 = ONE/(FOUR*PEM*SIG0)
	endif
C SECOND DERIVATIVE MATRIX
       DO I = 1,6
       DO J = 1,6
       X2N(I,J) = DFAC1*E(I)*E(J) + DFAC2*(DD2A(I,J)+DD2B(I,J))
       ENDDO
       ENDDO
	RETURN
      END
c------------------------------------------------------------------
      subroutine matmul3(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3,3),b(3,3),c(3,3)
       data zero /0.d0/
c
       do 10 i = 1,3
       do 10 j = 1,3
       c(i,j) = zero
       do 20 k = 1,3
c note reverse order for indices of a matrix, because it is dXm/dHn
20       c(i,j) = c(i,j) + a(k,j)*b(k,i)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul36(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3,3),b(3,6),c(3,6)
       data zero /0.d0/
c
       do 10 i = 1,3
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(i,k)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul63(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3,6),b(3,3),c(6,3)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,3
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(k,i)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul36i(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3,6),b(6,6),c(3,6)
       data zero /0.d0/
c
       do 10 i = 1,3
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,6
20       c(i,j) = c(i,j) + a(i,k)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul63i(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(6,3),b(3,3),c(6,3)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,3
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(i,k)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul63b(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(6,6),b(6,3),c(6,3)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,3
       c(i,j) = zero
       do 20 k = 1,6
20       c(i,j) = c(i,j) + a(k,i)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul313(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3),b(3,3),c(3)
       data zero /0.d0/
c
       do 10 i = 1,3
       c(i) = zero
       do 20 k = 1,3
20       c(i) = c(i) + a(k)*b(k,i)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul33(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3),b(3,3,3),c(3,3)
       data zero /0.d0/
c
       do 10 i = 1,3
       do 10 j = 1,3
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(k)*b(k,i,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul36r(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3,6),b(3,6),c(6,6)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(k,i)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul66(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(6,3),b(3,6),c(6,6)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(i,k)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul6(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(6,6),b(6,6),c(6,6)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,6
20       c(i,j) = c(i,j) + a(i,k)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul6r(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(6,6),b(6,6),c(6,6)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,6
20       c(i,j) = c(i,j) + a(k,i)*b(k,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine matmul336(a,b,c)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension a(3),b(3,6,6),c(6,6)
       data zero /0.d0/
c
       do 10 i = 1,6
       do 10 j = 1,6
       c(i,j) = zero
       do 20 k = 1,3
20       c(i,j) = c(i,j) + a(k)*b(k,i,j)
10	continue
       return
       end              
c------------------------------------------------------------------
      subroutine kyshr(pi,epr,dtime,ys,hr,props,nprops)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common/ki_select/i_hard,i_yield
      common/ksh_mtl/ams1,ams2,ams3,y_0,abac1,abac2,abac12,delk
	common /rate/ eratlim,erexp
      data zero,one /0.d0,1.d0/
	dimension props(nprops)
c Hardening
       if(iabs(i_hard).eq.1) then
          ys = ams1*(ams2+pi)**ams3
          hr = ams1*ams3*(ams2+pi)**(ams3-1.0)
        elseif(iabs(i_hard).eq.2) then
          ys = ams1 -ams2*exp(-ams3*pi)
          hr = ams2*ams3*exp(-ams3*pi)
        elseif(iabs(i_hard).eq.3) then
          nvalue=(nprops-34)/2
	    call ahard(ys,hr,pi,props(35),nvalue)

        endif
c include rate sensitivity provided plastic rate is greater than threshold
	if(i_hard .lt. 0) then
c determine type of rate sensitivity
		if(eratlim .gt. zero) then
c ABAQUS rate factor method for positive strain rate
			if(epr .gt. zero) then			
			efac = one + (epr/eratlim)**erexp
c additional factor in hardening to account for rate 
c derivative of rate factor wrt de_pl
			efacr = erexp*(efac - one)/(epr*dtime)
			ys = ys*efac
			hr = hr*efac + ys*efacr/efac
			endif
		else
c power law rate method for strain rate greater than threshold
			perlim = abs(eratlim)
			if(epr .gt. perlim) then			
			efac = (epr/perlim)**erexp
c additional factor in hardening to account for rate 
c derivative of rate factor wrt de_pl
			ys = ys*efac
			hr = hr*efac + ys*erexp/(epr*dtime)
			endif
		endif
	endif
      return
      end
c
c---------------------------------------------------------------------------
c---------------------------------------------------------------------------
c
      SUBROUTINE AHARD(SYIELD,hard,eqplas,TABLE,nvalue)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION TABLE(2,nvalue)    
c
c    If more than one entry, search table
c     if(nvalue.GT.1) then
        do 10 k1=1,nvalue-1
          EQPL1=TABLE(2,k1+1)
          if(eqplas.LT.EQPL1) then
            EQPL0=TABLE(2,k1)
 
c            if(EQPL1.LE.EQPL0) then
c                WRITE(6,1) 
c 1              FORMAT(//,30X,'***ERROR -PLASTIC STRAIN MUST BE ENTERED IN ASCENDING ORDER')
c                call XIT
c            endif
C         CURRENT YIELD STRESS AND HARDENING
            DEQPL=EQPL1-EQPL0
            SYIEL0=TABLE(1,k1)
            SYIEL1=TABLE(1,k1+1)
            DSYIEL=SYIEL1-SYIEL0
            hard=DSYIEL/DEQPL
            SYIELD=SYIEL0+(eqplas-EQPL0)*hard
            GOTO 20
          endif
 10     continue
 20     continue
c     endif
      RETURN
      END
