    !  Nullstellensuche.f90
    !
    !  FUNCTIONS:
    !  Nullstellensuche - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Nullstellensuche
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program Nullstellensuche

    implicit none

    integer :: j, k, l
    double precision, parameter :: stepwidth = 0.5d-2
    double precision :: posa = -1.5d0, posb = -1.5d0
    integer, parameter :: steps = 3d0/stepwidth
    integer :: iteration_steps
    logical :: error_flag = .false.
    logical :: nullstelle_enthalten
    double precision :: nullstellen((steps + 1)**2, 2)
    integer :: nullstellen_counter = 1
    double precision :: nullstelle_re, nullstelle_im
    complex :: z
    !open File
    open(1, file="nullstellen.dat")
    open(2, file="nullstellen_pos.dat")
    
    print*, 'Steps: ' , steps , ' at stepwidth ' , stepwidth
    do k= 0, steps
        posb = k*stepwidth - 15d-1

        do j = 0, steps
            posa = j*stepwidth - 15d-1
            call newton(posa, posb, nullstelle_re, nullstelle_im, iteration_steps, error_flag)
            !call newtonComplex(posa, posb, z, iteration_steps, error_flag) !Versuch mit der Komplexen Funtkoin zu arbeiten
            if (error_flag .eqv. .false.) then
                write(1,fmt='(I3.3, 1X, F6.3, 1X, F6.3)') iteration_steps, posa, posb !write into file
            else
                write(1,fmt='((A), 1X, F6.3, 1X, F6.3)') "Konv nicht", posa, posb !write into file
            end if
            nullstelle_enthalten = .false.
            do l = 1, nullstellen_counter
                if (abs(nullstellen(l,1) - nullstelle_re) < 10d-2 .and. abs(nullstellen(l,2) - nullstelle_im) < 10d-2) then
                !if (abs(nullstellen(l,1) - real(z)) < 10d-6 .and. abs(nullstellen(l,2) - aimag(z)) < 10d-6) then
                    nullstelle_enthalten = .true.
                    exit
                end if               
            end do
            if (nullstelle_enthalten .eqv. .false.) then
                nullstellen(nullstellen_counter, 1) = nullstelle_re
                !nullstellen(nullstellen_counter, 1) = real(z)
                nullstellen(nullstellen_counter, 2) = nullstelle_im
                !nullstellen(nullstellen_counter, 2) = aimag(z)
                nullstellen_counter = nullstellen_counter + 1
                !write(2,fmt='(F8.4,1X, F8.4)') real(z), aimag(z)
                write(2,fmt='(F8.4,1X, F8.4)') nullstelle_re, nullstelle_im
            end if
        end do
        write(1,*)
        posa = -1.5d0
    end do
    
    close(1)
    close(2)
    
    end program Nullstellensuche

    subroutine newton(aInput,bInput,a ,b,i,error) ! Real und Imaginärteil des Startwerts und Anzahl der Iterationen
        implicit none
        integer :: i
        double precision a, b, rs, is, ra, ia, rb, ib, aInput, bInput
        logical :: error
        a = aInput
        b = bInput
        i = 0
        error = .false.
        do while (.true.)
            i = i + 1
            call f(a,b,rs,is)
            call fa(a,b, ra, ia)
            call fb(a,b, rb, ib)
            !a = (rs - a* ra) / ra !liefert das gleiche Ergebnis
            a = a - rs / ra
            !b = (is - b*ib) / ib
            b = b - is / ib
            
            call f(a,b, rs,is)
            if (i > 200) then
                error =.true.
                exit
            end if
            if (abs(rs) < 10d-8 .and. abs(is) < 10d-8) exit
        end do
    end subroutine
    subroutine newtonComplex(aInput, bInput, z, i ,error)
        implicit none
        integer :: i
        double precision :: aInput, bInput
        complex :: z, c1
        logical :: error
        i = 0
        error = .false.
        z = cmplx(aInput, bInput)
        c1= cmplx(1,0)
        do while(.true.)
            i = i + 1
            z = z - (z**4-c1)/(4*z**3)
            if (abs(real(z)) < 10d-6 .and. abs(aimag(z)) < 10d-6) exit
            if (i > 200) then
                error = .true.
                exit
            end if
        end do
        
    
    end subroutine
    subroutine f(a, b, re, im) !Funktionswert an der Stelle z=a+i b zurückgegeben als real und imaginärteil aufgeteilt
        implicit none
        double precision :: a,b, re, im
        re = a**4 - 6*a**2*b**2 + b**4 -1
        im = 4*a**3*b - 4*a*b**3
        !re = a**2 - b**2 + 1
        !im = 2 * a * b
        return
    end subroutine
    subroutine fa(a,b,re,im) !Funktionswert an der Ableitung nach a an der Stelle z=a+i b zurückgegeben als real und imaginärteil aufgeteilt
        implicit none
        double precision :: a, b, re, im
        re = 4*a**3 - 12*a*b**2
        im = 12*b*a**2 - 4*b**3
        !re = 2*a
        !im = 2 * b
        return
    end subroutine

    subroutine fb(a,b,re,im)
        implicit none
        double precision :: a,b, re, im
        re = 4*b**3 - 12 *b*a**2
        im = 4*a**3 - 12*a*b**2
        !re = -2 * b
        !im = 2* a
        return
    end subroutine



