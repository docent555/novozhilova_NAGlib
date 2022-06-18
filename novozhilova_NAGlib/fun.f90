module fun
    use, intrinsic :: iso_c_binding
    use ifcore

    integer(c_int) ne, nt, nz, freq_out
    real(c_double) zex, dz, tend, dtr(2), q(3), i(2), th(2), a(2), dcir(2), r(2), f0(3), dt, pitch, f10, f20, f30
    complex(c_double_complex) fp(2)

    integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :)
    complex(c_double_complex), allocatable, target :: mean(:)
    real(c_double), allocatable, target :: tax(:), zax(:), u(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :), p(:, :)

    complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
    real(c_double), parameter :: pi = 2.0D0*dacos(0.0D0)

    private freq_out, zex, tend, dtr, q, i, th, a, dcir, r, f0, pitch

contains

    subroutine init()
        implicit none

        integer(c_int) ii

        call read_param()

        nt = tend/dt + 1
        nz = zex/dz + 1

        call allocate_arrays()

        f(1, 1) = f10
        f(3, 1) = f20
        f(5, 1) = f30

        do ii = 1, nt
            tax(ii) = (ii - 1)*dt
        end do

        do ii = 1, nz
            zax(ii) = (ii - 1)*dz
        end do

        call calc_u(u, zex, nz, zax)

        do ii = 1, 2
            idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
            idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
        end do

        !open (1, file='test.dat')
        !do ii = 1, ne
        !    write (1,'(4i)') idxre(1,ii), idxim(1,ii), idxre(2,ii), idxim(2,ii)
        !end do
        !close (1)
        !stop

    end subroutine init

    subroutine allocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_alloc

        allocate (f(6, nt), p(4*ne, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), w1(3, nt - 1), &
                  idxre(2, ne), idxim(2, ne), stat=err_alloc)

        if (err_alloc /= 0) then
            print *, "allocation error"
            pause
            stop
        end if
    end subroutine allocate_arrays

    subroutine deallocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_dealloc

        deallocate (f, p, u, tax, zax, mean, eta, etag, w, w1, stat=err_dealloc)

        if (err_dealloc /= 0) then
            print *, "deallocation error"
            pause
            stop
        end if
    end subroutine deallocate_arrays

    subroutine read_param() bind(c, name='read_param')
        use, intrinsic :: iso_c_binding
        import
        implicit none

        namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, &
            dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

        real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2

        open (unit=1, file='input_fortran.dat', status='old', err=101)
        read (unit=1, nml=param, err=102)
        close (unit=1)

        q(1) = q1
        q(2) = q2
        q(3) = q3
        i(1) = i1
        i(2) = i2
        th(1) = th1
        th(2) = th2
        a(1) = a1
        a(2) = a2
        dtr(1) = dtr1
        dtr(2) = dtr2
        dcir(1) = dcir1
        dcir(2) = dcir2
        r(1) = r1
        r(2) = r2

        write (*, nml=param)

        return
101     print *, "error of file open"; pause; stop
102     print *, 'error of reading file "input_fortran.dat"'; pause; stop
    end subroutine read_param

    subroutine ode4f(dydt, y, neq, nt, t0, h)
        import
        implicit none

        interface
            function dydt(t, y) result(s)
                use, intrinsic :: iso_c_binding
                implicit none
                real(c_double) t, y(:), s(size(y))
            end function dydt
        end interface

        integer(c_int) nt, neq, i, j
        real(c_double) h, t0, t
        real(c_double), intent(inout) :: y(:, :)
        real(c_double) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27

        !solve eq. at t=0
        fp(1) = y(1, 1)*exp(ic*y(2, 1))
        fp(2) = y(3, 1)*exp(ic*y(4, 1))
        call ode4p(dpdz, p, ne, nz, 0.0d0, dz)
        !CALL D02BGF(0.0d0,zex,4*ne,Y,TOL,HMAX,M,VAL,FCN,W,IFAIL)
        eta(:, 1) = eff(p(:, nz))
        etag(:, 1) = pitch**2/(pitch**2 + 1)*eta(:, 1)

        do i = 1, nt - 1
            v = y(:, i)
            t = t0 + (i - 1)*h
            s1(:) = dydt(t, v)
            s2(:) = dydt(t + h/2, v + h*s1(:)/2)
            s3(:) = dydt(t + h/2, v + h*s2(:)/2)
            s4(:) = dydt(t + h, v + h*s3(:))
            y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6

            eta(:, i + 1) = eff(p(:, nz))
            etag(:, i + 1) = pitch**2/(pitch**2 + 1)*eta(:, i + 1)

            write (*, '(a,f12.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,\,a)') 'Time = ', t + h, '   |F1| = ', abs(y(1, i + 1)), '   |F2| = ', abs(y(3, i + 1)), '   |F3| = ', abs(y(5, i + 1)), &
                '   Eff1 = ', eta(1, i + 1), '   Eff2 = ', eta(2, i + 1), char(13)

            pressed = peekcharqq()
            if (pressed) then
                key = getcharqq()
                if (ichar(key) .eq. esc) then
                    write (*, '(/,a)') 'Quit?'
                    key = getcharqq()
                    if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                        open (1, file='F.dat')
                        do j = 1, i + 1
                            write (1, '(4e17.8)') (j - 1)*h, abs(y(1, j)), abs(y(2, j)), abs(y(3, j))
                        end do
                        close (2)
                        open (2, file='E.dat')
                        do j = 1, i + 1
                            write (2, '(5e17.8)') (j - 1)*h, eta(1, j), etag(1, j), eta(2, j), etag(2, j)
                        end do
                        close (2)
                        stop
                    end if
                end if
            end if
        end do
    end subroutine ode4f

    !function eff(pex, ne) result(eta)
    !    import
    !
    !    integer(c_int) ne
    !    real(c_double) eta(2)
    !    complex(c_double_complex), intent(in) :: pex(:)
    !
    !    eta(1) = 1 - sum(abs(pex(1:ne))**2)/ne
    !    eta(2) = 1 - sum(abs(pex(ne + 1:2*ne))**2)/ne
    !end function eff

    function eff(pex) result(eta)
        use, intrinsic :: iso_c_binding, only: c_double, c_int
        import, only:ne, idxre, idxim

        implicit none

        integer(c_int) i
        real(c_double) eta(2)
        real(c_double), intent(in) :: pex(:)

        do i = 1, 2
            eta(i) = 1 - sum(sqrt(pex(idxre(i, :))**2 + pex(idxim(i, :))**2))/ne
        end do
    end function eff

    subroutine ode4p(dydt, y, neq, nz, z0, h)!, params)
        import
        implicit none

        interface
            subroutine dydt(z, p, prhs)
                use, intrinsic :: iso_c_binding
                implicit none
                real(c_double), intent(in) :: z, p(:)
                real(c_double), intent(out) :: prhs(:)
            end subroutine dydt
        end interface

        integer(c_int) nz, neq, i
        real(c_double) h, z0, z
        real(c_double), intent(inout) :: y(:, :)
        real(c_double) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))

        do i = 1, nz - 1
            v = y(:, i)
            z = z0 + (i - 1)*h
            call dydt(z, v, s1)
            call dydt(z + h/2, v + h*s1(:)/2, s2)
            call dydt(z + h/2, v + h*s2(:)/2, s3)
            call dydt(z + h, v + h*s3(:), s4)
            y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6
        end do
    end subroutine ode4p

    !function dpdz(z, p) result(s)
    !    import
    !    implicit none
    !
    !    integer(c_int) i, idx(ne)
    !    real(c_double) u, z
    !    complex(c_double_complex), parameter :: ic = (0.0D0, 1.0D0)
    !    complex(c_double_complex) p(:), s(size(p, 1))
    !
    !    do i = 1, 2
    !        u = exp(-3*((z - zex/2)/(zex/2))**2)
    !        idx = (/(i - 1)*ne + 1:ne + (i - 1)*ne/)
    !        s(idx) = ic*(fp(i)*u - (dtr(i) + abs(p(idx))**2 - 1)*p(idx))
    !    end do
    !end function dpdz

    subroutine dpdz(z, p, prhs)
        import :: ne, zex, f, ic, dtr
        implicit none

        real(c_double), intent(in) :: z, p(:)
        real(c_double), intent(out) :: prhs(:)

        integer(c_int) i, reidx(ne), imidx(ne)
        real(c_double) u
        complex(c_double_complex) s(ne), ptmp(ne)

        u = exp(-3*((z - zex/2)/(zex/2))**2)

        do i = 1, 2
            ptmp = dcmplx(p(idxre(i, :)), p(idxim(i, :)))

            s = ic*(fp(i)*u - (dtr(i) + abs(ptmp)**2 - 1)*ptmp)

            prhs(idxre(i, :)) = real(s)
            prhs(idxim(i, :)) = imag(s)
        end do
    end subroutine dpdz

    complex(c_double_complex) function xi(p, num)
        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
        import, only:ne, nz, mean, u, dz, idxre, idxim

        implicit none

        integer(c_int) i, num
        real(c_double) p(:, :)

        do i = 1, nz
            mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i)), 1)/ne
        end do

        mean = u*mean
        !mean = dconjg(u)*mean

        xi = (0.5d0*(mean(1) + mean(2)) + sum(mean(2:nz - 1)))*dz
    end function

    function dfdt(t, f) result(s)
        implicit none

        integer(c_int) :: ii, jj, iter_num = 1, time_num = 1
        real(c_double) t, f(:), s(size(f)), &
            x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
            x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
            f1, f2, f3, phi1, phi2, phi3, a1, a2
        complex(c_double_complex) x1, x2

        fp(1) = f(1)*exp(ic*f(2))
        fp(2) = f(3)*exp(ic*f(4))

        call ode4p(dpdz, p, ne, nz, 0.0d0, dz)

        !open (1, file='test.dat')
        !do jj = 1, nz
        !    do ii = 1, 2*ne
        !        write (1, '(2f17.8)') dreal(p(ii, jj)), dimag(p(ii, jj))
        !    end do
        !end do
        !close (1)
        !stop

        x1 = xi(p(1:2*ne, :), 1)
        x2 = xi(p(2*ne + 1:4*ne, :), 1)

        x1r = real(x1)
        x1i = imag(x1)
        x2r = real(x2)
        x2i = imag(x2)

        f1 = f(1)
        phi1 = f(2)
        f2 = f(3)
        phi2 = f(4)
        f3 = f(5)
        phi3 = f(6)

        q31 = q(3)/q(1)
        i1 = i(1)
        r1 = r(1)
        th1 = th(1)
        dcir1 = dcir(1)
        cos1 = cos(f(2))
        sin1 = sin(f(2))

        q32 = q(3)/q(2)
        i2 = i(2)
        r2 = r(2)
        th2 = th(2)
        dcir2 = dcir(2)
        cos2 = cos(f(4))
        sin2 = sin(f(4))

        q3 = q(3)
        a1 = a(1)
        a2 = a(2)

        s(1) = (-f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*f3*cos(phi3 - phi1 - th1))*q31
        s(2) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*(f3/f1)*sin(phi3 - phi1 - th1))*q31

        s(3) = (-f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*f3*cos(phi3 - phi2 - th2))*q32
        s(4) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*(f3/f2)*sin(phi3 - phi2 - th2))*q32

        s(5) = -f3 + a1*f1*cos(phi1 - phi3) + a2*f2*cos(phi2 - phi3)
        s(6) = a1*f1/f3*sin(phi1 - phi3) + a2*f2/f3*sin(phi2 - phi3)

        if (mod(iter_num, 4) .eq. 0) then
            w1(1, time_num) = s(2)
            w1(2, time_num) = s(4)
            w1(3, time_num) = s(6)
            time_num = time_num + 1
        end if
        iter_num = iter_num + 1
    end function dfdt

    subroutine calc_u(u, zex, nz, zax)
        import
        implicit none

        integer(c_int), intent(in) :: nz
        real(c_double), intent(in) :: zex, zax(nz)
        real(c_double), intent(out) :: u(:)

        integer(c_int) i

        do i = 1, nz
            u(i) = exp(-3*((zax(i) - zex/2)/(zex/2))**2)
        end do

    end subroutine

end module fun
