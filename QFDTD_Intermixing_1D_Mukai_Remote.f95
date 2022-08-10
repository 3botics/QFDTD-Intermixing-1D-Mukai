program QFDTD_Intermixing_1D_Mukai_Remote
    implicit none
!*****************************************************************************************
!***************************************** Variables *************************************
    character(len=30):: Title
!   ---------------Constant---------------------------------------------------------------
    double precision :: muo,pi,epso,co,qe,h,h_rad,m0
!   ---------------Domain and Iteration---------------------------------------------------
    integer          :: k,m,n,t,tmax,EN,S_Max_Itr,NFFT,NF,iMov
    integer          :: ie,ib

    double precision :: dt,xp
    double precision :: xx,ddx
    double precision :: x
!   ---------------PML--------------------------------------------------------------------
    integer          :: pml_grid,minv
    double precision :: sigma0,pml_thick
    complex          :: Rpml
!   ---------------Material Properties----------------------------------------------------
    double precision :: mat_PML,mat_QW,mat_Bar
    double precision :: mat_InGaAsP,mat_InP,mat_GaAs,mat_AlGaAs

    double precision :: t_QW,t_Bar
    double precision :: t_InGaAsP,t_InP,t_GaAs,t_AlGaAs

    double precision :: Eg_QW,Eg_Bar
    double precision :: me_QW,me_Bar,mLH_QW,mLH_Bar,mHH_QW,mHH_Bar
    !--------InGaAsP-------
    double precision :: Eg_InGaAsP
    double precision :: me_InGaAsP,mLH_InGaAsP,mHH_InGaAsP
    !--------InGaAsP2-------
    double precision :: Eg_InGaAsP2
    double precision :: me_InGaAsP2,mLH_InGaAsP2,mHH_InGaAsP2
    !--------InP-------
    double precision :: Eg_InP
    double precision :: me_InP,mLH_InP,mHH_InP
    !--------GaAs-------
    double precision :: Eg_GaAs
    double precision :: me_GaAs,mLH_GaAs,mHH_GaAs
    !--------AlGaAs-------
    double precision :: Eg_AlGaAs
    double precision :: me_AlGaAs,mLH_AlGaAs,mHH_AlGaAs
!   ---------------Schrodinger Equation---------------------------------------------------
    double precision :: ra,psi_norm,Fs,phi_norm,phi0_norm,arg
!   ---------------Test Function----------------------------------------------------------
    double precision :: L2,L4a,L4b,L6a,L6b,L8a,L8b,L38a,L38b,L10a,L10b,L310a,L310b
    double precision :: time
    double precision :: nlambda,np,n0
!   ---------------Diffusion Equation-----------------------------------------------------
    double precision :: Comp_QW,Comp_Bar,F0,D1,D2,aDif,kDif,tDif,xDif,max_dEg
    integer          :: nDif

!*****************************************************************************************
!***************************************** Universal Constants ***************************
    parameter (Title='QFDTD_Intermixing_1D_Mukai')
    parameter (pi=4.0*atan(1.0),        muo=4.0*pi*1.0d-7,       epso=8.854187817d-12)
    parameter (co=1.0/dsqrt(muo*epso),  qe=1.602176487d-19,      m0=9.10938215d-31)
    parameter (h_rad=1.054571628d-34,   h=6.62606896d-34)

!***************************************** Computational Domain **************************
!   ---------------Iteration--------------------------------------------------------------
    parameter (tmax=3, iMov=2**10)
    parameter (NF=2**6)
    parameter (NFFT=23,S_Max_Itr=2**NFFT, EN=3)
    parameter (pml_grid = 100)
!   ---------------Domain-----------------------------------------------------------------
    parameter (t_InGaAsP = 5d-9, t_InP = 50d-9)
    parameter (t_QW = t_InGaAsP, t_Bar = t_InP)
    parameter (xx = t_Bar + t_QW + t_Bar, ddx = 0.25d-9)
    parameter (ie = nint(xx/ddx) + 2*pml_grid, ib = ie + 1)
    parameter (pml_thick = float(pml_grid)*ddx)
!   ---------------Initial Wave-----------------------------------------------------------
    parameter (xp = 2.0d-9)
    parameter (nlambda=50*ddx,   np=10*ddx,  n0=3*np)
!   ---------------PML--------------------------------------------------------------------
    parameter (sigma0=0.0001, Rpml=cmplx(cos(pi/4),sin(pi/4)))
!   ---------------Material Parameters----------------------------------------------------
    !--------InGaAsP-------
    parameter (Eg_InGaAsP = 0.1)
    parameter (me_InGaAsP = 0.041*m0, mHH_InGaAsP = 0.5*m0, mLH_InGaAsP = 0.052*m0)
    !--------InGaAsP2-------
    parameter (Eg_InGaAsP2 = 0.1)
    parameter (me_InGaAsP2 = 0.044*m0, mHH_InGaAsP2 = 0.49*m0, mLH_InGaAsP2 = 0.052*m0)
    !--------InP-------
    parameter (Eg_InP= 0.1)
    parameter (me_InP = 0.08*m0, mHH_InP = 0.56*m0, mLH_InP = 0.052*m0)
    !--------GaAs-------
!    parameter (Eg_GaAs = 0.1)
    parameter (me_GaAs = 0.067*m0,mHH_GaAs = 0.42*m0, mLH_GaAs = 0.052*m0)
    !--------AlGaAs-------
!    parameter (Eg_AlGaAs = 0.1)
    parameter (me_AlGaAs = 0.088*m0,mHH_AlGaAs = 0.53*m0, mLH_AlGaAs = 0.052*m0)
!   ---------------Interdiffusion---------------------------------------------------------
    parameter (Comp_QW = 0, Comp_Bar = 0.25)
    parameter (F0 = Comp_Bar)
    parameter (D1 = 6.9e-23, D2 = 6.9e-23, kDif = 1 ) ! Al[x]Ga[1-x]As
!    parameter (D1 = 2.1e-23, D2 = 2.1e-21, kDif = 30) ! In[1-x]Ga[x]As[y]P[1-y]
    parameter (aDif = sqrt(D1/D2))
    parameter (nDif = 1)

!*****************************************************************************************
!***************************************** Vectors ***************************************
    double precision, dimension(ib) :: U,Eg,dEg,dEc,Vs,me,mHH,mLH,mat
    double precision, dimension(ib) :: F_comp,sumD1,sumD2,x_comp,y_comp
    double precision, dimension(ib) :: psir,psim,phir,phim
    double precision, dimension(ib) :: ca_psim,cb_psim,cc_psim,ca_psir,cb_psir,cc_psir
    double precision, dimension(ib) :: sigmax,gamma_re,gamma_im,delta_psir,delta_psim
    double precision, dimension(EN) :: Eig_Energy,Eig_Freq
    double precision, dimension(ib,EN) :: psir_init,phirf,phimf
    double precision, dimension(S_Max_Itr) :: psirt,psimt,Han_Win,freq
    integer, dimension(EN) :: Probe

!*****************************************************************************************
!***************************************** Open Files ************************************
    open(20,file = 'QFDTD_Intermixing_1D_Mukai_param',status='replace')
    open(21,file = 'QFDTD_Intermixing_1D_Mukai_U_init',status='replace')
    open(22,file = 'QFDTD_Intermixing_1D_Mukai_psir_init',status = 'replace')
    open(23,file = 'QFDTD_Intermixing_1D_Mukai_psim_init',status = 'replace')
    open(24,file = 'QFDTD_Intermixing_1D_Mukai_psirt',status = 'replace')
    open(25,file = 'QFDTD_Intermixing_1D_Mukai_psimt',status = 'replace')
    open(26,file = 'QFDTD_Intermixing_1D_Mukai_psirtf',status = 'replace')
    open(27,file = 'QFDTD_Intermixing_1D_Mukai_phirf',status = 'replace')
    open(28,file = 'QFDTD_Intermixing_1D_Mukai_phimf',status = 'replace')
    open(29,file = 'QFDTD_Intermixing_1D_Mukai_phir',status = 'replace')

    open(30,file = 'QFDTD_Intermixing_1D_Mukai_me',status = 'replace')
    open(31,file = 'QFDTD_Intermixing_1D_Mukai_HanWin',status = 'replace')
    open(32,file = 'QFDTD_Intermixing_1D_Mukai_dEc_init',status = 'replace')
    open(33,file = 'QFDTD_Intermixing_1D_Mukai_dEg_init',status = 'replace')
    open(35,file = 'QFDTD_Intermixing_1D_Mukai_dEc',status = 'replace')
    open(36,file = 'QFDTD_Intermixing_1D_Mukai_dEg',status = 'replace')

    open(41,file = 'QFDTD_Intermixing_1D_Mukai_psirt_mov',access = 'stream', form = 'unformatted',status = 'replace')
    open(42,file = 'QFDTD_Intermixing_1D_Mukai_psimt_mov',access = 'stream', form = 'unformatted',status = 'replace')
    open(43,file = 'QFDTD_Intermixing_1D_Mukai_U_mov',access = 'stream', form = 'unformatted',status = 'replace')
    open(44,file = 'QFDTD_Intermixing_1D_Mukai_Atom_Com_mov',access = 'stream', form = 'unformatted',status = 'replace')
    open(45,file = 'QFDTD_Intermixing_1D_Mukai_phir_mov',access = 'stream', form = 'unformatted',status = 'replace')

    open(51,file = 'QFDTD_Intermixing_1D_Mukai_sigma',status = 'replace')
    open(52,file = 'QFDTD_Intermixing_1D_Mukai_gamma_re',status = 'replace')
    open(53,file = 'QFDTD_Intermixing_1D_Mukai_gamma_im',status = 'replace')

    open(61,file = 'QFDTD_Intermixing_1D_Mukai_mat', access = 'stream',FORM = 'unformatted')

    open(81,file = 'QFDTD_Intermixing_1D_Mukai_Atom_Com_Init',status='replace')
    open(82,file = 'QFDTD_Intermixing_1D_Mukai_Atom_Com_Final',status='replace')

    open(101,file = 'QFDTD_Intermixing_1D_Mukai_Eig_Energy',status = 'replace')

!*****************************************************************************************
!***************************************** Initialization ********************************
    print*,'pml_thick',pml_thick
    dt = 0.125d-19
    ra = (0.5*h_rad/m0)*(dt/ddx**2)
    Fs = 1/dt;
!   ---------------Interdiffusion Initialization------------------------------------------
    sumD1 = 0;  sumD2 = 0;

    F_comp = 0;
    F_comp(1:nint((pml_thick+t_Bar)/ddx)) = Comp_Bar;
    F_comp(nint((pml_thick+t_Bar+t_QW)/ddx):ib) = Comp_Bar;

    write(81,*)F_comp; close(81);
!   ---------------PML Initialization-----------------------------------------------------
    sigmax = 0; gamma_re = 0; gamma_im = 0;
    do m = 1,(pml_grid+1)
        minv = ib - (m-1)
        sigmax(m) = sigma0*(float(m - (pml_grid+1)))**2
        sigmax(minv) = sigma0*(float(m - (pml_grid+1)))**2
    enddo
    gamma_re = real((1/(1 + Rpml*sigmax))**2)
    gamma_im = aimag((1/(1 + Rpml*sigmax))**2)

    write(51,*)sigmax;      close(51);
    write(52,*)gamma_re;    close(52);
    write(53,*)gamma_im;    close(53);

!   ---------------Material Initialization------------------------------------------------
    me = m0;                mHH = m0;           mLH = m0;
    me_QW = me_GaAs;        me_Bar = me_AlGaAs;
    mHH_QW = mHH_GaAs;      mHH_Bar = mHH_AlGaAs;
    mLH_QW = mLH_GaAs;      mLH_Bar = mLH_AlGaAs;

    Eg = 0; dEg = 0; dEc = 0;
!   ---------------In[1-x]Ga[x]As[y]P[1-y]------------------------------------------------
!    x_comp = 0.47
!    y_comp = 1 - 0
!    Eg_InGaAsP = 1.35 + 0.672*x_comp - 1.091*y_Comp + 0.758*(x_comp**2) + 0.101*(y_Comp**2) + &
!         0.111*x_comp*y_Comp - 0.58*(x_comp**2)*y_Comp - 0.159*x_comp*(y_Comp**2) + 0.268*(x_comp**2)*(y_Comp**2)
!    Eg_InP = 1.35

!   ---------------Al[x]Ga[1-x]As---------------------------------------------------------
    Eg_AlGaAs = 1.424 + 1.247*0.25
    Eg_GaAs = 1.424

    Eg_QW = Eg_GaAs;     Eg_Bar = Eg_AlGaAs;

    do m = 1,ib
        x = float(m)*ddx
        if(m.le.nint((pml_thick + t_Bar)/ddx).or.m.ge.nint((pml_thick + t_Bar + t_QW)/ddx)) then
            Eg(m) = Eg_Bar;
            me(m) = me_Bar;
            mHH(m) = mHH_Bar;
            mLH(m) = mLH_Bar;
        elseif(m.gt.nint((pml_thick + t_Bar)/ddx).and.m.lt.nint((pml_thick + t_Bar +t_QW)/ddx)) then
            Eg(m) = Eg_QW;
            me(m) = me_QW;
            mHH(m) = mHH_QW;
            mLH(m) = mLH_QW;
        endif
    enddo

!   ---------------Domain Initialization--------------------------------------------------
    mat_PML = 1.0;
    mat_InGaAsP = 0.5;      mat_InP = 0.755;
    mat_QW = mat_InGaAsP;   mat_Bar = mat_InP;
    mat = mat_PML;

    mat (1:pml_grid) = mat_PML;
    mat (pml_grid + 1 :nint((pml_thick + t_Bar)/ddx)) = mat_Bar;
    mat (nint((pml_thick + t_Bar)/ddx) + 1 : nint((pml_thick + t_Bar + t_QW)/ddx)-1) = mat_QW;
    mat (nint((pml_thick + t_Bar + t_QW)/ddx) : nint((pml_thick + 2*t_Bar + t_QW)/ddx)) = mat_Bar;

!   ---------------Writing Parameters-----------------------------------------------------
    write(20,*)ib,ddx,dt,tmax,iMov,pml_grid,S_Max_Itr,NFFT; close(20);
    write(30,*)me;      close(30);
    write(61)mat;       close(61);

!   ---------------Schrodinger Initialization---------------------------------------------
!    me = 0.067*m0;
    ca_psir = 1.0; cb_psir = -h_rad*dt/(2*me*ddx**2); cc_psir =  dt/h_rad;
    ca_psim = 1.0; cb_psim =  h_rad*dt/(2*me*ddx**2); cc_psim = -dt/h_rad;

!   ---------------Psi Initialization-----------------------------------------------------
    L2 = (2*pml_thick + xx)/2.0
!    L4a = pml_thick + t_Bar + t_QW/4.0 ;        L4b = pml_thick + t_Bar + t_QW - t_QW/4.0
!    L6a = pml_thick + t_Bar + t_QW/6.0 ;        L6b = pml_thick + t_Bar + t_QW - t_QW/6.0
!    L8a = pml_thick + t_Bar + t_QW/8.0 ;        L8b = pml_thick + t_Bar + t_QW - t_QW/8.0
!    L38a = pml_thick + t_Bar + 3.0*t_QW/8.0 ;   L38b = pml_thick + t_Bar + t_QW - 3.0*t_QW/8.0
!    L10a = pml_thick + t_Bar + t_QW/10.0 ;      L10b = pml_thick + t_Bar + t_QW - t_QW/10.0
!    L310a = pml_thick + t_Bar + 3.0*t_QW/10.0 ; L310b = pml_thick + t_Bar + t_QW - 3.0*t_QW/10.0

    L4a = pml_thick + (t_Bar + t_QW)/4.0 ;        L4b = pml_thick + 3*(t_Bar + t_QW)/4.0
    L6a = pml_thick + (t_Bar + t_QW)/6.0 ;        L6b = pml_thick + 5*(t_Bar + t_QW)/6.0
    L8a = pml_thick + t_Bar + t_QW/8.0 ;        L8b = pml_thick + t_Bar + t_QW - t_QW/8.0
    L38a = pml_thick + t_Bar + 3.0*t_QW/8.0 ;   L38b = pml_thick + t_Bar + t_QW - 3.0*t_QW/8.0
    L10a = pml_thick + t_Bar + t_QW/10.0 ;      L10b = pml_thick + t_Bar + t_QW - t_QW/10.0
    L310a = pml_thick + t_Bar + 3.0*t_QW/10.0 ; L310b = pml_thick + t_Bar + t_QW - 3.0*t_QW/10.0

!    Probe = (/nint(L2/ddx),nint(L4a/ddx),nint(L6a/ddx),nint(L8a/ddx),nint(L10a/ddx)/)
    Probe = (/nint(L2/ddx),nint(L4a/ddx),nint(L6a/ddx)/)
    print*,'Probe 1',Probe(1),Probe(1)*ddx/1.0d-9,L2
    print*,'Probe 2',Probe(2),Probe(2)*ddx/1.0d-9,L4a
    print*,'Probe 3',Probe(3),Probe(3)*ddx/1.0d-9,L6a
!    print*,'Probe 4',Probe(3),Probe(3)*ddx/1.0d-9,L8a
!    print*,'Probe 5',Probe(3),Probe(3)*ddx/1.0d-9,L10a

    psir_init = 0.0
    psir = 0.0; psim = 0.0;

    do n = 1,ib
        x = float(n)*ddx

        psir_init(n,1) = exp(-((x-L2)/xp)**2)
        psir_init(n,2) = exp(-((x-L4a)/xp)**2) - exp(-((x-L4b)/xp)**2)
        psir_init(n,3) = exp(-((x-L6a)/xp)**2) - exp(-((x-L2)/xp)**2) + exp(-((x-L6b)/xp)**2)
!        psir_init(n,4) = exp(-((x-L8a)/xp)**2) - exp(-((x-L38a)/xp)**2) + exp(-((x-L38b)/xp)**2) - &
!                       exp(-((x-L8b)/xp)**2)
!        psir_init(n,5) = exp(-((x-L10a)/xp)**2)- exp(-((x-L310a)/xp)**2) + exp(-((x-L2)/xp)**2) - &
!                        exp(-((x-L310b)/xp)**2) + exp(-((x-L10b)/xp)**2)
    enddo
    write(22,*)psir_init; close(22);

!   ---------------Band Offset and Potential Energy-------------------------------------------------------
    dEg = Eg - Eg_QW
    dEc = 0.65*dEg ! Al[x]Ga[1-x]As
!    dEc = 0.4*dEg ! In[1-x]Ga[x]As[y]P[1-y]

    max_dEg = maxval(dEg)
    Vs = 0
    U = -qe*Vs + qe*dEc
    write(21,*)U/qe;    close(21);
    write(32,*)dEc;  close(32);
    write(33,*)dEg;  close(33);

!    print*,'stop'; stop;

!   ---------------Hanning Window---------------------------------------------------------
    do m = 1,S_Max_Itr
        Han_Win(m) = 0.5*(1-cos(2*pi*m/S_Max_Itr))
    enddo
    write(31,*)Han_Win; close(31);

!   ---------------Frequency Axis---------------------------------------------------------
    do m = 1,S_Max_Itr
        freq(m) = Fs*float(m-1)/S_Max_Itr - Fs/2.0
    enddo
    write(26,*)freq*h/qe;  close(26);

!*****************************************************************************************
!***************************************** Printing Monitor ******************************
    print*, '------------------------------------------------------'
    print*, 'Title : ',Title
    print*, 'Note  : ...'
    print*, '------------------------------------------------------'
    print*, 'dt : ',dt
    print*, 'ra : ',ra
    print*, '------------------------------------------------------'
    print*, 'Ie : ',ie, ' points'
    print*, 'Ib : ',ib, ' points'
    print*, '------------------------------------------------------'
    print*, 'x  : ',xx,  ' m'
    print*, 'dx : ',ddx, ' m'
    print*, '------------------------------------------------------'
    print*, 'NFFT : ',NFFT
    print*, '------------------------------------------------------'
    print*, 'Bar : ',nint((pml_thick + t_Bar)/ddx)
    print*, 'QW  : ',nint((pml_thick + t_Bar + t_QW)/ddx)
    print*, '------------------------------------------------------'
    print*, 'D1 : ',D1
    print*, 'D2  : ',D2
    print*, '------------------------------------------------------'
!    print*, 'stop'; stop;

!*****************************************************************************************
!***************************************** Time Loop *************************************
    print*,'..start main loop..'

!    do t=1,1
!        tDif = float(t)*60*60*2
!        if(mod(t,1).eq.0) then
!            print*,'Diffuse time : ',tDif,'second'
!        endif

!***************************************** Initial Function *******************************
    do k = 1,EN
        psir = 0.0 ; psim = 0.0;
        psir(1:ib) = psir_init(1:ib,k)

        psi_norm = sqrt(sum(psir(1:ib)**2 + psim(1:ib)**2))
        psir = psir/psi_norm;
        psim = psim/psi_norm;
        psi_norm = sqrt(sum(psir(1:ib)**2 + psim(1:ib)**2))
        print*,'Psi Norm :', psi_norm

!***************************************** Main Loop **************************************
    do n = 1,S_Max_Itr
        if(mod(n,S_Max_Itr/20).eq.0) print*,n
        time=dt*float(n)

        delta_psir(2:ie) = psir(1:ie-1) - 2*psir(2:ie) + psir(3:ib)
        delta_psim(2:ie) = psim(1:ie-1) - 2*psim(2:ie) + psim(3:ib)

!        psir(2:ie) = psir(2:ie) + cb_psir(2:ie)*delta_psim(2:ie) + cc_psir(2:ie)*U(2:ie)*psim(2:ie)
!        psim(2:ie) = psim(2:ie) + cb_psim(2:ie)*delta_psir(2:ie) + cc_psim(2:ie)*U(2:ie)*psir(2:ie)

        psir(2:ie) = psir(2:ie) + cb_psir(2:ie)*(delta_psim(2:ie)*gamma_re(2:ie) + &
                                                 delta_psir(2:ie)*gamma_im(2:ie)) + &
                                  cc_psir(2:ie)*U(2:ie)*psim(2:ie)
        psim(2:ie) = psim(2:ie) + cb_psim(2:ie)*(delta_psir(2:ie)*gamma_re(2:ie) - &
                                                 delta_psim(2:ie)*gamma_im(2:ie)) + &
                                  cc_psim(2:ie)*U(2:ie)*psir(2:ie)

        psirt(n) =  psir(Probe(k))*Han_Win(n)
        psimt(n) =  psim(Probe(k))*Han_Win(n)

        write(24,*)psirt(n)
        write(25,*)psimt(n)

        if(mod(n,iMov).eq.0) write(41)psir
        if(mod(n,iMov).eq.0) write(42)psim
    enddo

!***************************************** Eigen Energy ***********************************
        call FFTC(psirt,psimt,S_Max_Itr,NFFT)
        psirt = dsqrt(psirt**2 + psimt**2)
        write(28,*)psirt;           close(28);

        Eig_Freq(1) = 2*pi*freq(maxloc(psirt,1))
        Eig_Energy(1) = abs(h_rad * Eig_Freq(1))
        arg = Eig_Freq(1) * dt

        print*,'Eigen Energy', k ,':' ,Eig_Energy(k)/qe * 1000.0 , 'meV'
        write(101,*)Eig_Energy;     close(101);

!***************************************** Eigen Function *********************************
!       psir = 0.0 ; psim = 0.0; phir = 0.0 ; phim = 0.0
!       psir(1:ib) = psir_init(1:ib,k)
!
!       psi_norm = sqrt(sum(psir(1:ib)**2 + psim(1:ib)**2))
!       psir = psir/psi_norm;
!       psim = psim/psi_norm;
!       psi_norm = sqrt(sum(psir(1:ib)**2 + psim(1:ib)**2))
!       print*,'Psi Norm :', psi_norm
!
!   do m=1,S_Max_Itr
!       delta_psir(2:ie) = psir(1:ie-1) - 2*psir(2:ie) + psir(3:ib)
!       delta_psim(2:ie) = psim(1:ie-1) - 2*psim(2:ie) + psim(3:ib)
!
!        psir(2:ie) = psir(2:ie) + cb_psir(2:ie)*delta_psim(2:ie) + cc_psir(2:ie)*U(2:ie)*psim(2:ie)
!        psim(2:ie) = psim(2:ie) + cb_psim(2:ie)*delta_psir(2:ie) + cc_psim(2:ie)*U(2:ie)*psir(2:ie)
!
!       psir(2:ie) = psir(2:ie) + cb_psir(2:ie)*(delta_psim(2:ie)*gamma_re(2:ie) + &
!                                                delta_psir(2:ie)*gamma_im(2:ie)) + &
!                                 cc_psir(2:ie)*U(2:ie)*psim(2:ie)
!       psim(2:ie) = psim(2:ie) + cb_psim(2:ie)*(delta_psir(2:ie)*gamma_re(2:ie) - &
!                                                delta_psim(2:ie)*gamma_im(2:ie)) + &
!                                 cc_psim(2:ie)*U(2:ie)*psir(2:ie)
!
!       phir(2:ie) = phir(2:ie) + Han_Win(m)*(psir(2:ie)*cos(arg*m) + psim(2:ie)*sin(arg*m));
!       phim(2:ie) = phim(2:ie) + Han_Win(m)*(psim(2:ie)*cos(arg*m) - psir(2:ie)*sin(arg*m));
!
!       if(mod(m,NF).eq.0) write(45)phir
!   enddo
!       phi_norm = sqrt(sum(phir**2 + phim**2))
!       phir = phir/phi_norm;
!       phim = phim/phi_norm;
!       phi_norm = sqrt(sum(phir**2 + phim**2))
!       print*, 'Phi Norm : ', phi_norm
!
!       phirf(1:ib,k) = phir(1:ib)
!       phimf(1:ib,k) = phim(1:ib)
   enddo
!
!   write(27,*)phirf;
!   write(28,*)phimf;

    do m = 1,EN
        print*,'Eigen Energy', m ,':' ,Eig_Energy(m)/qe*1000.0 , 'meV'
    enddo

    print*, '------------------------------------------------------'
    print*, 'Diffusion : '
!***************************************** Atomic Composition *****************************
    do m = 1,ib
        x = float(m)*ddx;
        xDif = x - xx/2 - pml_thick;

        if (m .le. nint((pml_thick+t_Bar)/ddx) .or. m >= nint((pml_thick+t_Bar+t_QW)/ddx )) then
            do n = 1,nDif
            sumD1(m) = sumD1(m) + (((1 - aDif*kDif)/(1 + aDif*kDif))**(2*(n-1)))* &
                    (1 - erf((abs(xDif) + t_QW/2 + (4*n - 2)*aDif*t_QW/2)/(2*sqrt(D1*tDif))) + &
                    ((1 - aDif*kDif)/(1 + aDif*kDif)) * &
                    (1 - erf((abs(xDif) - t_QW/2 + 4*n*aDif*t_QW/2)/(2*sqrt(D1*tDif)))));
            enddo
            F_comp(m) = F0/2*(1 + erf((abs(xDif) - t_QW/2)/(2*sqrt(D1*tDif))) - &
                    ((1 - aDif*kDif)/(1 + aDif*kDif))*(1 - erf((abs(xDif) - t_QW/2)/(2*sqrt(D1*tDif)))) +&
                    (4*aDif*kDif/((1 + aDif*kDif)**2))*sumD1(m));

        elseif (m .gt. nint((pml_thick+t_Bar)/ddx) .and. m .lt. nint((pml_thick+t_Bar+t_QW)/ddx)) then

            do n = 1,nDif
            sumD2(m) = sumD2(m) + (((1 - aDif*kDif)/(1 + aDif*kDif))**(2*(n-1))) * &
                    (2 + erf((abs(xDif) - (4*n - 3)*t_QW/2)/(2*sqrt(D2*tDif))) - &
                    erf((abs(xDif) + (4*n - 3)*t_QW/2)/(2*sqrt(D2*tDif))) + &
                    ((1 - aDif*kDif)/(1 + aDif*kDif)) * &
                    (2 - erf((abs(xDif )+ (4*n - 1)*t_QW/2)/(2*sqrt(D1*tDif))) + &
                    erf((abs(xDif) - (4*n - 1)*t_QW/2)/(2*sqrt(D2*tDif)))));
            enddo
            F_comp(m) = F0*aDif/(1 + aDif*kDif)*sumD2(m);
        endif
    enddo
    write(44) F_comp
    write(82,*) F_comp

!***************************************** Band Offset and Potential Energy ***************
    dEg = F_comp/maxval(F_comp)*max_dEg;
!    dEc = 0.40*dEg ! In[1-x]Ga[x]As[y]P[1-y]
    dEc = 0.65*dEg ! Al[x]Ga[1-x]As

    U = -qe*Vs + qe*dEc
    write(43) U/qe
    write(35,*) dEc
    write(36,*) dEg

    sumD1 = 0;  sumD2 = 0;
    F_comp = 0;
    F_comp(1:nint((pml_thick+t_Bar)/ddx)) = Comp_Bar;
    F_comp(nint((pml_thick+t_Bar+t_QW)/ddx):ib) = Comp_Bar;
!***************************************** End of Time Loop *******************************
!    enddo

    close(27); close(28);
    close(44); close(43);

!**************************************************** Program Ended ***********************
    print*, '------------------------------------------------------'
    print*, 'Title : ',Title
    print*, 'Note  : ...'
    print*, '----------------FINISH, ALHAMDULILLAAH!!!-------------'

end program QFDTD_Intermixing_1D_Mukai_Remote


!******************************************************************************************
!*********************************************** FFT **************************************
SUBROUTINE FFTC (AR,AI,N,M)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: N,M
    INTEGER :: N1,N2,I,J,K,L,L1,L2
    double precision :: PI,A1,A2,Q,U,V
    double precision, INTENT (INOUT), DIMENSION (N) :: AR,AI

    PI = 4.0*ATAN(1.0)
    N2 = N/2
    N1 = 2**M
    IF(N1.NE.N) STOP 'Indices do not match'
!   --------------------------------------------------------------------------------------
!   Rearrange the data to the bit reversed order
    L = 1
    DO K = 1, N-1
        IF (K.LT.L) THEN
            A1    = AR(L)
            A2    = AI(L)
            AR(L) = AR(K)
            AR(K) = A1
            AI(L) = AI(K)
            AI(K) = A2
        END IF
        J   = N2
        DO WHILE (J.LT.L)
            L = L-J
            J = J/2
        END DO
        L = L+J
    END DO
!   --------------------------------------------------------------------------------------
!   Perform additions at all levels with reordered data
    L2 = 1
    DO L = 1, M
        Q  =  0.0
        L1 =  L2
        L2 =  2*L1
        DO K = 1, L1
            U   =  COS(Q)
            V   = -SIN(Q)
            Q   =  Q + PI/L1
            DO J = K, N, L2
                I     =  J + L1
                A1    =  AR(I)*U-AI(I)*V
                A2    =  AR(I)*V+AI(I)*U
                AR(I) =  AR(J)-A1
                AR(J) =  AR(J)+A1
                AI(I) =  AI(J)-A2
                AI(J) =  AI(J)+A2
            END DO
        END DO
    END DO
!   --------------------------------------------------------------------------------------
    AR = cshift(AR,-N2)
    AI = cshift(AI,-N2)
END SUBROUTINE FFTC

