module fun
    use, intrinsic :: iso_c_binding
    use ifcore

    integer(c_int) ne, nt, nz, freq_out, neqp, lenwrk, l, method, ifail, neqf, lwork, liwork, nrd
 real(c_double) zex, dz, tend, dtr(2), q(3), i(2), th(2), a(2), dcir(2), r(2), f0(3), dt, pitch, f10, f20, f30, ptol, hstart, zstart
    complex(c_double_complex) fp(2)
    logical(c_bool) errass

    integer(c_int) breaknum(3)
    real(c_double) phitmp0(3), phitmp1(3)
    complex(c_double_complex) fc, fcarr(3)

    integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :), iwork(:)
    complex(c_double_complex), allocatable, target :: mean(:)
    real(c_double), allocatable, target :: tax(:), zax(:), u(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :), p(:, :), &
                                           phi(:, :), phios(:, :), wos(:, :), &
                                           pgot(:), ppgot(:), pmax(:), thres(:), workp(:), work(:)

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

        neqp = 4*ne
        lenwrk = 32*neqp
        zstart = 0.0d0
        errass = .false.
        hstart = 0.0d0
        ptol = 1.0d-7
        method = 2
        ifail = 0

        neqf = 6
        nrd = 6
        lwork = 8*neqf + 5*nrd + 21
        liwork = nrd + 21

        call allocate_arrays()

        do l = 1, neqp
            thres(l) = 1.0d-8
        end do

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

        allocate (f(6, nt), p(neqp, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), w1(3, nt - 1), &
          idxre(2, ne), idxim(2, ne), pgot(neqp), ppgot(neqp), pmax(neqp), thres(neqp), workp(lenwrk), work(lwork), iwork(liwork), &
                  wos(3, nt - 1), phi(3, nt), phios(3, nt), stat=err_alloc)

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

    subroutine ode4f()
        import
        implicit none

        integer(c_int) neq, i, j, ipar, iout, idid, itol
        real(c_double) :: h, t, zwant, zgot, rpar, rtol, atol
        real(c_double) y(6)
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27

        call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifail)

        !solve eq. at t=0
        fp(1) = f(1, 1)*exp(ic*f(2, 1))
        fp(2) = f(3, 1)*exp(ic*f(4, 1))

        do i = 1, nz - 1
            zwant = i*dz
            call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifail)
            p(:, i + 1) = pgot
        end do

        eta(:, 1) = eff(p(:, nz))
        etag(:, 1) = pitch**2/(pitch**2 + 1)*eta(:, 1)

        !open (1, file='test.dat')
        !do j = 1, nz
        !    write (1, '(5f17.8)') (j-1)*dz, p(14, j), p(ne+14, j), p(128, j), p(ne+128, j)
        !end do
        !close (1)
        !stop

        rpar = 0.0
        ipar = 0
        t = tax(1)
        y = f(:, 1)
        itol = 0
        rtol = 1.0d-7
        atol = rtol
        iout = 6

        iwork(:) = 0
        work(:) = 0.0d0

        iwork(5) = 6

        call dopri5(6, dfdt, t, y, tend, rtol, atol, itol, solout, iout, work, lwork, iwork, liwork, rpar, ipar, idid)

        eta(:, nt) = eff(p(:, nz))
        etag(:, nt) = pitch**2/(pitch**2 + 1)*eta(:, nt)
        do j = 1, neqf
            f(j, nt) = y(j)
        end do

    end subroutine ode4f

    function eff(pex) result(eta)
        use, intrinsic :: iso_c_binding, only: c_double, c_int
        import, only:ne, idxre, idxim

        implicit none

        integer(c_int) i
        real(c_double) eta(2)
        real(c_double), intent(in) :: pex(:)

        do i = 1, 2
            eta(i) = 1 - sum(abs(dcmplx(pex(idxre(i, :)), pex(idxim(i, :))))**2)/ne
        end do
    end function eff

    subroutine dpdz(z, p, prhs)
        import :: ne, zex, f, ic, dtr
        implicit none

        real(c_double) z, p(*), prhs(*)

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

    subroutine dfdt(neqf, t, f, s, rpar, ipar)
        implicit none

        integer(c_int) :: ii, jj, iter_num = 1, time_num = 1, neqf, ipar
        real(c_double) t, f(neqf), s(neqf), &
            x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
            x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
            f1, f2, f3, phi1, phi2, phi3, a1, a2, zwant, zgot, rpar
        complex(c_double_complex) x1, x2

        call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifail)

        fp(1) = f(1)*exp(ic*f(2))
        fp(2) = f(3)*exp(ic*f(4))

        do jj = 1, nz - 1
            zwant = zax(jj + 1)
            call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifail)
            p(:, jj + 1) = pgot
        end do

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
    end subroutine dfdt

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

    subroutine solout(nr, xold, x, y, n, con, icomp, nd, rpar, ipar, irtrn)
        implicit none

        interface
            function contd5(ii, x, con, icomp, nd)
                implicit double precision(a - h, o - z)
                dimension con(5*nd), icomp(nd)
            end
        end interface

        integer(c_int) nr, n, nd, icomp(nd), ipar, irtrn, j, it
        real(c_double) xold, x, con(5*nd), rpar, y(neqf), xout
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/intern/xout, it

        if (nr .eq. 1) then
            write (*, '(a,f12.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,\,a)') 'Time = ', xout, '   |F1| = ', abs(y(1)), '   |F2| = ', abs(y(3)), '   |F3| = ', abs(y(5)), &
                '   Eff1 = ', eta(1, nr), '   Eff2 = ', eta(2, nr), char(13)
            xout = x + dt
            it = 1
            eta(:, it) = eff(p(:, nz))
            etag(:, it) = pitch**2/(pitch**2 + 1)*eta(:, it)
            do j = 1, neqf
                f(j, it) = y(j)
            end do
        else
10          continue
            if (x .ge. xout) then
                write (*, '(a,f12.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,\,a)') 'Time = ', xout, '   |F1| = ', abs(y(1)), '   |F2| = ', abs(y(3)), '   |F3| = ', abs(y(5)), &
                    '   Eff1 = ', eta(1, nr), '   Eff2 = ', eta(2, nr), char(13)
                xout = xout + dt
                it = it + 1
                eta(:, it) = eff(p(:, nz))
                etag(:, it) = pitch**2/(pitch**2 + 1)*eta(:, it)
                do j = 1, neqf
                    f(j, it) = y(j)
                end do
                goto 10
            end if
        end if

        pressed = peekcharqq()
        if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
                write (*, '(/,a)') 'Quit?'
                key = getcharqq()
                if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                    nt = nr
                    return
                end if
            end if
        end if
        return
    end

end module fun
