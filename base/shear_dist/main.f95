! A fortran95 program for G95
! By WQY
program dist_shear
  implicit none
  integer op, cont1, cont2, points, aux
  real*16 lbd_w, lbd_b, Tr_w,Tr_b,Tb_med,Tw_med, Tw_max,Tb_max, T0 !shear stress variables
  real*16 yw, yb, P, Pw, Pb, Rh, S, b, h, fp, Lb, Lw, Sfw, Csf, A !geometry variables
  real*16 db, dw, x(500), y(500), T(500), loc,gy, fx, dfx, erro, erro_max, x_w, x_b, x_new, xt


!Calculates the shear stress in each position;
!It assumes rectangular section

10  write(*,*)'Which distribution you want to use?'
  write(*,*)'Press 1 for Sterling and Knight (2002) or 2 for PoMCE'
  write(*,*)'Press 0 to end the program'
  read(*,*) op

  write(*,*)'Please, type in the flow depth.'
  read(*,*) h
  write(*,*)'Please, type in the flow width.'
  read(*,*) b

  write(*,*) 'Please type in the average shear stress in the wall (Pa)'
  read(*,*) Tw_med
  write(*,*) 'Please type in the maximum shear stress in the wall (Pa)'
  read(*,*) Tw_max
  write(*,*) 'Please type in the average shear stress in the bed (Pa)'
  read(*,*) Tb_med
  write(*,*) 'Please type in the maximum shear stress in the bed (Pa)'
  read(*,*) Tb_max


  points = 100 !number of points for calculation of shear stress in each sector
  erro_max = 1E-5 !for calibration of T in the PoMCE method

  !1. Geometric properties of the section
!  S = 0.00076
!  b = 0.316
!  h = 0.158
!
!  A = b*h
!  P = b + 2.*h
!
!  Rh = A/P
!  T0 = 9810*S*Rh

  Pb = b
  Pw = 2.*h

  Lb = b/2.
  Lw = h

  db = Lb/points
  dw = Lw/points

!  !1.1 Calculating shear stress parameters based on Knight and Sterling (2000)
!  fp = Pb/Pw
!
!  if (fp .lt. 4.374) then
!        Csf = 1.
!        else
!            Csf = 0.6603*(fp**0.28125)
!  end if
!
!  Sfw = -3.23*log10(fp/1.38 + 1) + 4.6052
!  Sfw = 0.01*Csf*exp(Sfw)
!
!  Tw_med = T0 * Sfw * (1+fp)
!  Tb_med = T0 * (1-Sfw) * (1 + 1/fp)
!  Tw_max = T0 * Sfw * 2.0372 * (fp**0.7108)
!  Tb_max = T0 * (1-Sfw) * 2.1697 * (fp**(-0.3287))
!
!!  write(*,*)Csf, Sfw, T0, Tw_med,Tb_med,Tw_max,Tb_max

  Tr_w = Tw_med/Tw_max
  Tr_b = Tb_med/Tb_max

  if (op .eq. 0) then
    goto 20

    !2. Calculating distribution based on sterling equations
    elseif (op .eq. 1) then
        write(*,*) ' '
        write(*,*) 'You selected the method of Sterling!'
        write(*,*) ' '
        write(*,*)'please type the calibrated lambda for the wall. It should be a positive number.'
        read(*,*)lbd_w
        write(*,*)'please type the calibrated lambda for the bed. It should be a positive number.'
        read(*,*)lbd_b


        if (lbd_w*Tw_max .le. 40) then
        !2.1 shear stress distribution in the wall
        loc = dw
        cont1 = 0
        do while (loc .le. Lw)
            T(cont1) = (1/lbd_w) * log(1 + (exp(lbd_w*Tw_max) - 1) * (loc/Lw))
            x(cont1) = 0.
            y(cont1) = Lw - loc
!            write(*,*)cont1,loc, T(cont1)
            cont1 = cont1+1
            loc = loc + dw

        end do

        else
            loc = dw
            cont1 = 0
            do while (loc .le. Lw)
                T(cont1) = Tw_max + (1/lbd_w) * log(loc/Lw)
                x(cont1) = 0.
                y(cont1) = Lw - loc
!                write(*,*)cont1,loc, T(cont1)
                cont1 = cont1+1
                loc = loc + dw
            end do
        end if

        if (lbd_b*Tb_max .le. 40) then
        !2.2 shear stress distribution in the wall
        loc = db
        cont2 = cont1
        do while (loc .le. Lb)
            T(cont2) = (1/lbd_b) * log(1 + (exp(lbd_b*Tb_max) - 1) * (loc/Lb))
            x(cont2) = 0.
            y(cont2) = 0.
!            write(*,*)cont2,loc, T(cont2)
            cont2 = cont2 + 1
            loc = loc + db
        end do

        else
            loc = db
            cont2 = cont1
            do while (loc .le. Lb)
                T(cont2) = Tb_max + (1/lbd_b) * log(loc/Lb)
                x(cont2) = 0.
                y(cont2) = Lw - loc
!                write(*,*)cont1,loc, T(cont1)
                cont2 = cont2 + 1
                loc = loc + dw
            end do
        end if

        !3. Calculating distribution based on sterling equations
        elseif (op .eq. 2) then
            write(*,*) ' '
            write(*,*) 'You selected the PoMCE-Based method!'
            write(*,*) ' '
            write(*,*)'please type the calibrated lambda for the wall. It should be a negative number.'
            read(*,*)lbd_w
            write(*,*)'please type the calibrated lambda for the bed. It should be a negative number.'
            read(*,*)lbd_b

            !3.1 Shear stress distribution in the wall
!            if (Tr_w .lt. 0.97) then
            if (lbd_w .gt. -40) then
                loc = dw
                cont1 = 0
                do while (loc .le. Lw)
                    gy = 1. - exp(-lbd_w) * (lbd_w + 1.)
                    gy = 1. - gy * (loc/Lw)
                    erro = 1000.
                    x_w = 0.8
                    do while (erro .gt. erro_max)
                        fx = exp(-lbd_w * x_w) * (lbd_w * x_w + 1.) - gy
                        dfx = -1. * (lbd_w**2) * x_w * exp(-lbd_w * x_w)

                        x_new = x_w - fx/dfx

                        erro = abs(x_new - x_w)
                        x_w = x_new
                    end do
                    T(cont1) = x_w
                    x(cont1) = 0.
                    y(cont1) = Lw - loc
!                    write(*,*)cont1,loc, T(cont1), 'here'
                    cont1 = cont1+1
                    loc = loc + dw
                end do

                else
                    xt =  dw
                    loc = dw
!                    write(*,*)xt,loc,dw
                    cont1 = 0
                    do while (loc .le. xt)
                        gy = 1. - exp(-lbd_w) * (lbd_w + 1.)
                        gy = 1. - gy * (loc/Lw)
!                        write(*,*)loc,gy
                        erro = 1000.
                        x_w = 0.8
                        aux = 0
                        do while (erro .gt. erro_max)
                            fx = exp(-lbd_w * x_w) * (lbd_w * x_w + 1.) - gy
                            dfx = -1. * (lbd_w**2) * x_w * exp(-lbd_w * x_w)

                            x_new = x_w - fx/dfx
                            erro = abs(x_new - x_w)
                            x_w = x_new
                            aux = aux+1
!                            write(*,*) aux, x_w
                        end do
                        T(cont1) = x_w
                        x(cont1) = 0.
                        y(cont1) = Lw - loc
                        write(*,*)cont1,loc, T(cont1), 'here! wwwwww'
                        cont1 = cont1+1
                        loc = loc + dw
!                        write(*,*)loc,'here here!'
                    end do

                    do while (loc .le. Lw)
                        gy = log(-lbd_w*loc/Lw) - lbd_w
                        erro = 1000.
                        x_w = 0.8
                        do while (erro .gt. erro_max)
                            fx = log(-lbd_w*x_w) - lbd_w*x_w - gy
                            dfx = 1/x_w - lbd_w

                            x_new = x_w - fx/dfx
                            erro = abs(x_new - x_w)
                            x_w = x_new
                        end do
                        T(cont1) = x_w
                        x(cont1) = 0.
                        y(cont1) = Lw - loc
                        write(*,*)cont1,loc, T(cont1)
                        cont1 = cont1+1
                        loc = loc + dw
                    end do

            end if

            !3.1 Shear stress distribution in the bed
!            if (Tr_b .lt. 0.97) then
            if (lbd_b .gt. -40) then
                loc = db
                cont2 = 0
                do while (loc .le. Lb)
                    gy = 1. - exp(-lbd_b) * (lbd_b + 1.)
                    gy = 1. - gy * (loc/Lb)
                    erro = 1000.
                    x_b = 0.8
                    do while (erro .gt. erro_max)
                        fx = exp(-lbd_b * x_b) * (lbd_b * x_b + 1.) - gy
                        dfx = -1. * (lbd_b**2) * x_b * exp(-lbd_b * x_b)

                        x_new = x_b - fx/dfx

                        erro = abs(x_new - x_b)
                        x_b = x_new
                    end do
                    T(cont2) = x_b
                    x(cont2) = 0.
                    y(cont2) = Lb - loc
!                    write(*,*)cont2,loc, T(cont1), 'here'
                    cont2 = cont2 + 1
                    loc = loc + db
                end do

                else
                    xt =  db
                    loc = db
!                    write(*,*)xt,loc,dw
                    cont2 = 0
                    do while (loc .le. xt)
                        gy = 1. - exp(-lbd_b) * (lbd_b + 1.)
                        gy = 1. - gy * (loc/Lb)
!                        write(*,*)loc,gy
                        erro = 1000.
                        x_b = 0.8
                        aux = 0
                        do while (erro .gt. erro_max)
                            fx = exp(-lbd_b * x_b) * (lbd_b * x_b + 1.) - gy
                            dfx = -1. * (lbd_b**2) * x_b * exp(-lbd_b * x_b)

                            x_new = x_b - fx/dfx
                            erro = abs(x_new - x_b)
                            x_b = x_new
                            aux = aux+1
!                            write(*,*) aux, x_w
                        end do
                        T(cont2) = x_b
                        x(cont2) = 0.
                        y(cont2) = Lb - loc
                        write(*,*)cont2,loc, T(cont2), 'here! bbbb'
                        cont2 = cont2+1
                        loc = loc + db
!                        write(*,*)loc,'here here!'
                    end do

                    do while (loc .le. Lb)
                        gy = log(-lbd_b*loc/Lb) - lbd_b
                        erro = 1000.
                        x_b = 0.8
                        do while (erro .gt. erro_max)
                            fx = log(-lbd_b*x_b) - lbd_b*x_b - gy
                            dfx = 1/x_b - lbd_b

                            x_new = x_b - fx/dfx
                            erro = abs(x_new - x_b)
                            x_b = x_new
                        end do
                        T(cont2) = x_b
                        x(cont2) = 0.
                        y(cont2) = Lb - loc
                        write(*,*)cont2,loc, T(cont2)
                        cont2 = cont2+1
                        loc = loc + db

                    end do

            end if



  end if

20 end
