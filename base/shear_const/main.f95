! A fortran95 program for G95
! By WQY
program shear_const
  implicit none
  integer op, cont1, cont2, points, aux
  real*16 lbd_w, lbd_b, Tr_w,Tr_b,Tb_med,Tw_med, Tw_max,Tb_max, T0 !shear stress variables
  real*16 yw, yb, P, Pw, Pb, Rh, S, b, h, fp, Lb, Lw, Sfw, Csf, A !geometry variables
  real*16 db, dw, x(500), y(500), T(500), loc,gy, fx, dfx, erro, erro_max, x_w, x_b, x_new, xt


!Calculates the shear stress in each position;
!It assumes rectangular section

  write(*,*)'Please, type in the flow depth.'
  read(*,*) h
  write(*,*)'Please, type in the flow width.'
  read(*,*) b
  write(*,*)'Please, type in the flow slope.'
  read(*,*) S


  !1. Geometric properties of the section
!  S = 0.00076
!  b = 0.316
!  h = 0.158

  A = b*h
  P = b + 2.*h

  Rh = A/P
  T0 = 9810*S*Rh

  Pb = b
  Pw = 2.*h

  Lb = b/2.
  Lw = h

!  db = Lb/points
!  dw = Lw/points

  !1.1 Calculating shear stress parameters based on Knight and Sterling (2000)
  fp = Pb/Pw

  if (fp .lt. 4.374) then
        Csf = 1.
        else
            Csf = 0.6603*(fp**0.28125)
  end if

  Sfw = -3.23*log10(fp/1.38 + 1) + 4.6052
  Sfw = 0.01*Csf*exp(Sfw)

  Tw_med = T0 * Sfw * (1+fp)
  Tb_med = T0 * (1-Sfw) * (1 + 1/fp)
  Tw_max = T0 * Sfw * 2.0372 * (fp**0.7108)
  Tb_max = T0 * (1-Sfw) * 2.1697 * (fp**(-0.3287))

  write(*,*)'Tw_med = ', Tw_med
  write(*,*)'Tw_max = ', Tw_max
  write(*,*)'Tb_med = ', Tb_med
  write(*,*)'Tb_max = ', Tb_max
  write(*,*)'T0 = ',T0
end
