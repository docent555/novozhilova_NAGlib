program sys15f
    use, intrinsic :: iso_c_binding
    use fun
    use ifport

    implicit none

    integer(c_int) i, j, hours, minutes, seconds
    real(c_double) start_time, stop_time, calc_time
    complex(c_double_complex) pc

    call init()

    do i = 1, ne
        p(i, 1) = real(exp(ic*(i - 1)/dble(ne)*2*pi))
        p(ne + i, 1) = imag(exp(ic*(i - 1)/dble(ne)*2*pi))
        p(2*ne + i, 1) = real(exp(ic*(i - 1)/dble(ne)*2*pi))
        p(3*ne + i, 1) = imag(exp(ic*(i - 1)/dble(ne)*2*pi))
    end do

    write (*, '(/)')

    start_time = dclock()
    call ode4f()
    stop_time = dclock()

    calc_time = stop_time - start_time

    hours = calc_time/3600
    minutes = (calc_time - hours*3600)/60
    seconds = calc_time - hours*3600 - minutes*60

    write (*, '(/)')
    print *, 'Calcualting took:', hours, 'h :', minutes, 'm :', seconds, 's'

    do i = 2, nt
        do j = 1, 3
            w(j, i - 1) = imag(log(f(2*j - 1, i)*exp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*exp(ic*f(2*j, i - 1)))))/dt
        end do
    end do

    phi(:, 1) = 0; 
    do i = 2, nt
        do j = 1, 3
            phi(j, i) = phi(j, i - 1) + imag(log(f(2*j - 1, i)*exp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*exp(ic*f(2*j, i - 1)))))
        end do
    end do

    breaknum(:) = 0
    phitmp0(1) = f(2*1 - 1, 1)*exp(ic*f(2*1, 1))
    phitmp0(2) = f(2*2 - 1, 1)*exp(ic*f(2*2, 1))
    phitmp0(3) = f(2*3 - 1, 1)*exp(ic*f(2*3, 1))
    !phitmp0(:) = datan2(dimag(f(:, 1)), dreal(f(:, 1)))
    phios(:, 1) = phitmp0(:)
    do i = 2, nt
        do j = 1, 3
            fc = f(2*j - 1, i)*exp(ic*f(2*j, i))
            phitmp1(j) = datan2(dimag(fc), dreal(fc))
            if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
            if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
            phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
            !phios(j, i) = phitmp1(j)
            phitmp0(j) = phitmp1(j)
        end do
    end do

    do i = 1, nt - 1
        do j = 1, 3
            wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
        end do
    end do

    write (*, '(/)')

    pause

    open (1, file='F.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), abs(f(1, i)), abs(f(2, i)), abs(f(3, i))
    end do
    close (1)

    open (1, file='FCOMP.dat')
    do i = 1, nt
        fcarr(1) = f(2*1 - 1, i)*exp(ic*f(2*1, i))
        fcarr(2) = f(2*2 - 1, i)*exp(ic*f(2*2, i))
        fcarr(3) = f(2*3 - 1, i)*exp(ic*f(2*3, i))
        write (1, '(7e17.8)') tax(i), real(fcarr(1)), imag(fcarr(1)), real(fcarr(2)), imag(fcarr(2)), real(fcarr(3)), imag(fcarr(3))
    end do
    close (1)

    open (2, file='E.dat')
    do i = 1, nt
        write (2, '(5e17.8)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
    end do
    close (2)

    open (3, file='W.dat')
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), w(1, i), w(2, i), w(3, i)
    end do
    close (3)

    open (1, file='P.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), phi(1, i), phi(2, i), phi(3, i)
    end do
    close (1)

    open (1, file='POS.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), phios(1, i), phios(2, i), phios(3, i)
    end do
    close (1)

    open (3, file='WOS.dat')
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
    end do
    close (3)
    
    stop
101 print *, 'error of file open.'
    pause
    stop
102 print *, 'error of file reading.'
    pause
    stop
103 print *, 'error of file writing.'
    pause
    stop
end program
