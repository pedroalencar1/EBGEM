!programa com modelo de Foster and Lane (Erosion by concentred flow in farm fields - 1983)
!e Haan et.al. (Desing hydrology and sedimentology for small catchments  -1994)


program FLSM
IMPLICIT NONE
CHARACTER arquivo*30,arquivo1*30,arquivo2*30
integer No,i
REAL Q(1500),n,S,Kr,Ds,Pa,Tc,Xnc,Gc,Rn,WP,Rh,Ta,Wn,Weq,Le,Delta,Param,Tmax,tm(1500),t(1500),tn
REAL NEL,Tb,Gcf,Xncf,Xnc1,Xnc2,Gc1,Gc2,Rn1,Rn2,WP1,Lim,Ch,Phi,Ang,Dvcr,g,pi,es,Ws
REAL WP2,Rh1,Rh2,Ta1,Ta2,Dg1,Dg2,Tnc1,Tnc2,Xncf1,Xncf2,gcf1,gcf2,tnf1,tnf2,dgf1,dgf2
REAL wi(1500),wf,dp(1500),depth,DeltaW,PerimF,AreaI,AreaF,Rh_S,Max_wi,Max_W,dW
COMMON /grand1/Gcf,Xncf1,Xncf2
COMMON /grand2/Ch,g,Ds,es,Pa,Phi,depth,Ang,pi

!Phi = internal friction angle of the soil (degrees) > the data can be obtained in the NAVFAC 7.02, based in texture; in Ortiz et al 1989 based in texture, Plasticity intex and porosity or via laboratory experiments
!Ch = Cohesion (Pa) > the data can be obtained in the NAVFAC 7.02, based in texture; in Ortiz et al 1989 based in texture, Plasticity intex and porosity or via laboratory experiments
!Lim = a threshold for aplication of sidorchuk routine
!Q = peak discharge (m3/s)
!tm  = discharge's duration (min)
!n = Manning's number
!S = Slope (abs) - (not in degrees or percentage)
!Kr = rill erodibility factor (s/m)
!Tc = Critical shear stress (Pa)
!Xnc = Normalized distance in the WP where T=Tc (abs)
!Gc = Conveyance function at Xnc (abs)
!Rn = Normalized Hydraulic Radius (abs)
!WP = Wet Perimeter (m)
!Rh = Hydraulic Radius (m)
!Ta = Average shear stress (Pa)
!Wn = Normalized width (abs)
!Weq = equilibrium width (m)
!Er = Erosion rate (kg/m.s)
!Ve = Velocity of movement down (m/s)
!Pse = Weigth of soil eroded in the event (kg)
!Le = Sheet thickness of eroded soil(m)
!Param = A parameter that is repeated throughout the calculations
!t = time (s)
!NEL = Depth of nonerodible layer(m)
!Tb = Tension when the erosion reach the NEL (Pa)
!Gcf = conveyance funcion final (abs)
!Wf = final width of the channel (m)
!tne = time to reach  the NEL (s)

        WRITE(*,*)'******************************************'
        WRITE(*,*)'           Foster_Lane program            '
        WRITE(*,*)'                                          '
        WRITE(*,*)'    Calculates the final total of soil    '
        WRITE(*,*)'       eroded in an ephemeral gully       '
        WRITE(*,*)'        for a time-serie rainfall         '
        WRITE(*,*)'                                          '
        WRITE(*,*)'      Universidade Federal do Ceara       '
        WRITE(*,*)'    Departamento de Engenharia Agricola   '
        WRITE(*,*)'      Mestrado em Engenharia Agricola     '
        WRITE(*,*)'                                          '
        WRITE(*,*)'            Pedro Alencar, 2017           '
        WRITE(*,*)'                                          '
        WRITE(*,*)'******************************************'


!1. READS INPUT DATA
write(*,*)'Insert the number of observed events'
Read(*,*)No
write(*,'(a)')'Insert the name of the file containing the runoff data (Discharge e Duration): ' !in the absense of measured data we sugest use the 30-min intensity
Read(*,'(a30)')arquivo1
write(*,'(a)')'Insert the name of the file containing the hillslope and soil data: '
Read(*,'(a30)')arquivo2
WRITE(*,'(a)',advance='no') 'Insert the name of the output file: '
READ(*,'(a30)')arquivo

open(50,file=arquivo1)
    do i=1, No
        read(50,*) Q(i),tm(i)
        write(*,*) Q(i),tm(i)
        t(i)=tm(i)*60.
    end do
close(50)

open(60,file=arquivo2)
    read(60,*)n,S,Tc,kr,Ds,es,NEL,Ch,Phi,Lim
    write(*,*)'The number of Manning of the channel is....',n
    write(*,*)'The declivity of the hillslope is..........',S
    write(*,*)'The critical shear stress is...............',Tc
    write(*,*)'The rill erodibility coefficient is........',Kr
    write(*,*)'The soil Bulk density is...................',Ds
    write(*,*)'The porosity of the soil is................',es
    write(*,*)'The depth of the nonerodible layer is......',NEL
    write(*,*)'The soil cohesion is.......................',Ch
    write(*,*)'The internal friction angle is.............',Phi
    write(*,*)'The threshold for wall erosion is..........',Lim
close(60)

Pa = 9803. !N/m³
g = 9.803 !gravity
i = 1
depth = 0
AreaI = 0
Max_wi = 0
Max_W  = 0
pi = 3.1415926536
Phi = Phi*pi/180.
Ws = 0.
i=1

!2. FOSTER AND LANE EQUATIONS FOR INCISION BEFORE REACHING THE NEL
open(40,file=arquivo)

do while(i.le.No)

    Param = (n*Q(i)/(S**(0.5)))**0.375

    if (depth.lt.NEL) then
        Gc = Pa*S*Param/Tc
        if (Gc.lt.1.79) then
            i=i+1
            else

                Delta = sqrt(3.9429**2 - 4*6.9594*(1/Gc))
                Xnc1 = (3.9429-Delta)/(2*6.9594)
                Xnc2 = (3.9429+Delta)/(2*6.9594)
                Rn1 = -0.8834*Xnc1+0.1395*Xnc1+0.151
                Rn2 = -0.8834*Xnc2+0.1395*Xnc2+0.151

        end if

        if (Rn1.gt.0) then
            WP1 = Param/Rn1**(0.625)
            Rh1 = Rn1*WP1
            Ta1 = Pa*Rh1*S
            Tnc1 = Tc/Ta1
            Gc1 = 1./(Tnc1*(Rn1**0.375))
            Dg1 = abs(Gc-Gc1)

            else
                Gc1 = 0
                Dg1 = abs(Gc-Gc1)
        end if

        if (Rn2.gt.0) then
            WP2 = Param/Rn2**(0.625)
            Rh2 = Rn2*WP2
            Ta2 = Pa*Rh2*S
            Tnc2 = Tc/Ta2
            Gc2 = 1/(Tnc2*(Rn2**(0.375)))
            Dg2 = abs(Gc-Gc2)

            else
                Gc2 = 0
                Dg2 = abs(Gc-Gc2)
        end if

        if (Dg1.gt.Dg2) then
            Xnc = Xnc2
            else
                Xnc = Xnc1
        end if
        Write(*,*)'Xnc ',Xnc

        Wn = -1.4873*Xnc+0.7436
        Rn = -0.8834*Xnc**2+0.1395*Xnc+0.151
        WP = Param/Rn**(0.625)
        Weq = WP*Wn
        Rh = Rn*WP
        Ta = Pa*Rh*S
        Tmax = 1.35*Ta

        if (Tmax.lt.Tc) then
            Write(*,*)'The event didn´t cause erosion.'
            wi(i) = 0.
            dp(i) = 0.
            else
                wi(i) = Weq
                Le = Kr*(Tmax-Tc)*t(i)/Ds
                if (Le.gt.NEL) then
                    dp(i) = NEL
                    else
                        dp(i) = Le
                end if
                depth = depth + dp(i)
                if (Max_wi.lt.Weq) then
                    Max_wi = wi(i)
                    Max_W  = wi(i)
                end if
                AreaI = AreaI+wi(i)*dp(i)


            Write(40,*)'Fase1',i,wi(i),dp(i),Max_wi,AreaI
        end if
    end if

!3. FOSTER AND LANE EQUATIONS FOR INCISION AFTER REACHING THE NEL
    if (depth.ge.NEL) then
        Gcf = Pa*S*Param/Tc
        if (Gcf.lt.1.78) then
            i=i+1
            else
                CALL Newton_NEL
        end if

        tnf1 = 1.35*(1-(1-2.*xncf1)**2.9)
        gcf1 = 1/(tnf1*(xncf1*(1-2*xncf1))**(0.375))
        dgf1 = abs(gcf-gcf1)
        tnf2 = 1.35*(1-(1-2.*xncf2)**2.9)
        gcf2 = 1/(tnf2*(xncf2*(1-2.*xncf2))**(0.375))
        dgf2 = abs(gcf-gcf2)

            if (dgf1.gt.dgf2) then
                xncf = xncf2
                else
                    xncf = xncf1
            end if

        Tb = 1.35*Ta*(1-(1-2*xncf)**2.9)*Pa*S*Param*(xncf*(1-2*xncf))**0.375
        if (Tb.le.Tc)then
            dW = 0
            else
            dW = Kr*(Tb-Tc)/Ds
        end if

        wf = Param*((1-2.*xncf)/(xncf**1.667))**0.375

        if (Max_W.ge.Wf) then
            tn = 0.
            else
                tn = t(i)*dW/(Wf-Max_W)
        end if

        Max_W = (1-exp(-tn))*(Wf-Max_W)+Max_W

        write(40,*)'Fase2',i, Max_W, NEL, dw
    end if
    i = i+1
end do

DeltaW = Max_W-Max_Wi
AreaF = AreaI + DeltaW*depth
PerimF = Max_W + 2.*depth
Rh_S = AreaF/PerimF

!4. WALL TRANSFORMATION

if (AreaF.gt.Lim) then
    Dvcr = (2.*Ch*cos(Phi)/(g*Ds))/(sin(0.5*(Phi+pi/2.)))**2.
    if  (Dvcr.le.NEL) then
        CALL Newton_PHI
        Ws = Max_W +2.*NEL/tan(Ang)
        AreaF = NEL*(Max_W+Ws)/2.
    end if
end if

write(40,*)'The cross section area is: ',AreaF,' m2'

if (depth.lt.NEL) then
    write(40,*)'The depth is: ',depth,'m'
    write(40,*)' '
    write(40,*)'The erosion didn´t reach the non erodible layer.'
    else
        write(40,*)'The depth is: ',NEL,'m'
        write(40,*)' '
        write(40,*)'The erosion reached the non erodible layer.'
    end if

if (Ws.gt.0.) then
    write(40,*)'There is erosion of the walls.'
    write(40,*)'The top width is',Ws,'m and the wall slope is',100.*tan(Ang),'%'
    write(40,*)Ang*180/pi
end if

write(*,*)' '
write(*,*)'Please, see the output file'



end program


Subroutine Newton_NEL
!Subrotina para calcular as raizes da equação de Xncf
implicit none
integer i,j
real x,xa,xb,f,df,erro,gcf,xncf1,xncf2
COMMON /grand1/Gcf,Xncf1,Xncf2

Xa = 0.1
Xb = 0.4
i=0
j=0

erro=1000000.
X=Xa
do while (erro.gt.0.000001)

    f = -321.29*x**6+249.02*x**5+9.059*x**4-58.033*x**3+13.012*x**2+1.7199*x-0.004883*COS(x) -1/Gcf
    df = -321.29*6*x**5+249.02*5*x**4+9.059*4*x**3-58.033*3*x**2+13.012*2*x+1.7199+0.004883*sin(x)

    x = xa - f/df
    erro = abs(x-xa)/xa
    xa=x
    i=i+1
end do

erro=1000000.
x=Xb
do while (erro.gt.0.000001)

    f = -321.29*x**6+249.02*x**5+9.059*x**4-58.033*x**3+13.012*x**2+1.7199*x-0.004883*COS(x) -1/Gcf
    df = -321.29*6*x**5+249.02*5*x**4+9.059*4*x**3-58.033*3*x**2+13.012*2*x+1.7199+0.004883*sin(x)

    x = xb - f/df
    erro = abs(x-xb)/xb
    xb=x
    j=j+1
end do

Xncf1 = xa
Xncf2 = xb

end subroutine

Subroutine Newton_PHI
!Subroutine to calculate wal stable angle
    implicit none

    integer i
    real k1,k2,f1,df1,x,xa,erro,Ch,g,Ds,es,Pa,Phi,depth,Ang,pi
    COMMON /grand2/Ch,g,Ds,es,Pa,Phi,depth,Ang,pi

    i=0
    k1 = Ch/(g*Ds*depth)
    k2 = tan(Phi)*(Ds - es*1000.)/1000.

    erro=1000000.
    x = pi/4.
        do while (erro.gt.0.00000001)

        f1 = k2*(cos(x))**2 - sin(2*x)/2. - k1
        df1 = -2*k2*cos(x) - cos(2*x)

        xa = x - f1/df1
        erro = abs(x-xa)
        x=xa
        i=i+1

    end do
Ang = x
end subroutine
