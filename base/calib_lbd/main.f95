!Routine to calibrate the lagrange multipliers to assess shear stress in open channels' boundaries
!Using the method of Newton-Raphson
!References - Sterling and Knight (2002)
!           - Bonakdari et al. (2014)
!           - Nocedal and Wright (2006)


program calib_ldb
  implicit none
  integer n, cont1,cont2, op
  real*16 lbd1_w,lbd_w, lbd1_b, lbd_b, fs, dfs, fp,dfp, erro, erro_max, erro1, erro2, aux1, aux2
  real*16 Tb_med, Tb_max, Tw_med, Tw_max, x1, x2, x1_new, x2_new, Tmax, Tmed,Tr_w,Tr_b
  character arquivo*30

        WRITE(*,*)'********************************************'
        WRITE(*,*)'*           Newton-Raphson Method          *'
        WRITE(*,*)'*                                          *'
        WRITE(*,*)'*    Calculates the Lagrange multipliers   *'
        WRITE(*,*)'*     for PoME- and PoMCE-based methods    *'
        WRITE(*,*)'*                                          *'
        WRITE(*,*)'*      Technisches Universität Berlin      *'
        WRITE(*,*)'*           (Department of Ecology)        *'
        WRITE(*,*)'*        Federal University of Ceará       *'
        WRITE(*,*)'* (Department of Agricultural Engineering) *'
        WRITE(*,*)'*                                          *'
        WRITE(*,*)'*            Pedro Alencar, 2020           *'
        WRITE(*,*)'*                                          *'
        WRITE(*,*)'********************************************'
        WRITE(*,*)' '
        WRITE(*,*)' '
        WRITE(*,*)' '


10  write(*,*)'Which distribution you want to use?'
  write(*,*)'Press 1 for Sterling and Knight (2002) or 2 for PoMCE'
  write(*,*)'Press 0 to end the program'
  read(*,*) op

  write(*,*) 'Please type in the average shear stress in the wall (Pa)'
  read(*,*) Tw_med
  write(*,*) 'Please type in the maximum shear stress in the wall (Pa)'
  read(*,*) Tw_max
  write(*,*) 'Please type in the average shear stress in the bed (Pa)'
  read(*,*) Tb_med
  write(*,*) 'Please type in the maximum shear stress in the bed (Pa)'
  read(*,*) Tb_max
!b/h = 2

!Tw_med = 0.547814
!Tb_med = 0.628975
!Tw_max = 0.558226
!Tb_max = 0.682565


Tr_b = Tb_med/Tb_max
Tr_w = Tw_med/Tw_max

!write(*,*)Tr_w ,Tr_b

!b/h = 4
!Tw_med = 0.597
!Tb_med = 0.749
!Tw_max = 0.663
!Tb_max = 0.862

erro_max = 1E-6 !Defines precision of approximation

! 1. Calculates sterling's lambdas

!1.1 Calibrating lambda for the wall

if (Tr_w .lt. 0.97) then !check if it is possible a direct solution
    cont1 = 0
    x1 = 20 !initial guess
    erro1 = 1000.0
    do while(erro1 .gt. erro_max)
        fs = x1*exp(x1)
        fs = fs/(exp(x1) - 1)
        fs = fs - x1*Tr_w - 1

        dfs = exp(x1)
        dfs = dfs*(dfs - x1 - 1)
        dfs = dfs/((exp(x1) - 1)**2) - Tr_w

        x1_new = x1 - fs/dfs     !método de newton

        erro1 = abs(x1 - x1_new)
        x1 = x1_new
        cont1 = cont1+1
    end do
    lbd1_w = x1/Tw_max

    else
        write(*,*)' '
        write(*,*)' '
        write(*,*) 'Value of lambda is too big and can cause stack overflow. Modified equation is used. Please check documentation.'
        write(*,*)' '
        write(*,*)' '

        x1 = (1.- Tr_w)**(-1)
        cont1 = 1
        erro1 = -1
        lbd1_w = x1/Tw_max
    end if


!1.2 Calibrating lambda for the bed

if  (Tr_b .lt. 0.97) then !check if it is possible a direct solution
    cont2 = 0
    x1 = 20. !initial guess
    erro2 = 1000.0
    do while(erro2 .gt. erro_max)
        fs = x1*exp(x1)
        fs = fs/(exp(x1) - 1)
        fs = fs - x1*Tr_b - 1

        dfs = exp(x1)
        dfs = dfs*(dfs - x1 - 1)
        dfs = dfs/((exp(x1) - 1)**2) - Tr_b

        x1_new = x1 - fs/dfs     !método de newton

        erro2 = abs(x1 - x1_new)
        x1 = x1_new
        cont2 = cont2+1
    end do
    lbd1_b = x1/Tb_max

    else
        write(*,*)' '
        write(*,*)' '
        write(*,*) 'Value of lambda is too big and can cause stack overflow. Modified equation is used. Please check documentation.'
        write(*,*)' '
        write(*,*)' '

        x1 = (1-Tr_b)**(-1)
        erro2 = -1
        lbd1_b = x1/Tb_max

end if

if (op .eq. 0) then
    goto 20

elseif (op .eq.1) then

    write(*,*) ' '
    write(*,*) 'You selected the method of Sterling!'
    write(*,*) ' '

   lbd_w = lbd1_w
   lbd_b = lbd1_b


    write(*,*)'After ',cont1,' iterations and with estimation error equals to ', erro1
    write(*,*)'the Lagrange Multipliers Lambda 1 for the wall is ', lbd_w
    write(*,*)'After ',cont2,' iterations and with estimation error equals to ', erro2
    write(*,*)'the Lagrange Multipliers Lambda 1 for the bed is ', lbd_b
    write(*,*) ' '
    write(*,*) ' '
    write(*,*) ' '
    write(*,*) ' '

    goto 10

!2. POMCE method

    elseif (op .eq. 2) then !POMCE method is selected
        write(*,*) ' '
        write(*,*) 'You selected the PoMCE method!'
        write(*,*) ' '


        !2.1 Calibrating lambda for the wall
        if (Tr_w .lt. 0.97) then !check if it is possible a direct solution
            cont1 = 0
            x1 = -20 !initial guess
            erro1 = 1000.0
            do while(erro1 .gt. erro_max)
                fp = exp(x1) - x1 - 1.
                fp = 2/x1 - x1/fp - Tr_w

                dfp = (exp(x1) - x1 - 1)
                dfp = (x1*exp(x1) - x1)/dfp**2 - 2/(x1**2) - 1/dfp

                x1_new = x1 - fp/dfp     !método de newton

                erro1 = abs(x1 - x1_new)
                x1 = x1_new
                cont1 = cont1+1
            end do
            erro1 = erro
            lbd1_w = x1

            else !for indirect solution
                x1 = ((Tr_w + 2.)**2) - 8.
                x1 = (Tr_w - 2) - sqrt(x1)
                x1 = x1/(2 - 2*Tr_w)

                erro1 = -1
                cont1 = 1
                lbd1_w = x1
        end if

        !2.2 Calibrating lambda for the bed
        if (Tr_b .lt. 0.97) then
            cont2 = 0
            x1= -20 !initial guess
            erro2 = 1000.0

            do while(erro2 .gt. erro_max)
                fp = exp(x1) - x1 - 1.
                fp = 2/x1 - x1/fp - Tr_b

                dfp = (exp(x1) - x1 - 1)
                dfp = (x1*exp(x1) - x1)/dfp**2 - 2/(x1**2) - 1/dfp

                x1_new = x1 - fp/dfp     !método de newton

                erro2 = abs(x1 - x1_new)
                x1 = x1_new
                cont2 = cont2+1
            end do
            lbd1_b = x1

            else !for indirect solution
                x1 = ((Tr_b + 2.)**2) - 8.
                x1 = (Tr_b - 2) - sqrt(x1)
                x1 = x1/(2 - 2*Tr_b)

                erro2 = -1
                cont2 = 1
                lbd1_b = x1
    end if

   lbd_w = lbd1_w
   lbd_b = lbd1_b

    write(*,*)'After ',cont1,' iterations and with estimation error equals to ', erro1
    write(*,*)'the Lagrange Multipliers Lambda 1 for the wall is ', lbd_w
    write(*,*)'After ',cont2,' iterations and with estimation error equals to ', erro2
    write(*,*)'the Lagrange Multipliers Lambda 1 for the bed is ', lbd_b



    goto 10

else
    write(*,*) 'INVALID OPTION'
    write(*,*) '   '
    write(*,*) '   '

    goto 10

  end if

20 end
