!
!   Polarizable Embedding (PE) library
!   Copyright (C) 2013-2018, 2020 The PE library developers. See the CONTRIBUTORS file
!                                 in the top-level directory of this distribution.
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jogvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!

module pelib_blas_interfaces

    use pelib_precision

    implicit none

contains

!------------------------------------------------------------------------------

function nrm2(x)

    real(rp), external :: dnrm2

    real(rp) :: nrm2
    real(rp), dimension(:), intent(in) :: x

    integer(ip) :: n, incx

    incx = 1

    n = size(x)

    nrm2 = dnrm2(n, x(1), incx)

end function nrm2

!------------------------------------------------------------------------------

function dot(x,y)

    real(rp), external :: ddot

    real(rp) :: dot
    real(rp), dimension(:), intent(in) :: x, y

    integer(ip) :: n, incx, incy

    incx = 1
    incy = 1

    n = size(x)

    dot = ddot(n, x(1), incx, y(1), incy)

end function dot

!------------------------------------------------------------------------------

subroutine axpy(x, y, a)

    external :: daxpy

    real(rp), intent(in), optional :: a
    real(rp), dimension(:), intent(in) :: x
    real(rp), dimension(:), intent(inout) :: y

    real(rp) :: o_a
    integer(ip) :: n, incx, incy

    if (present(a)) then
        o_a = a
    else
        o_a = 1.0
    end if

    incx = 1
    incy = 1

    n = size(x)

    call daxpy(n, o_a, x(1), incx, y(1), incy)

end subroutine axpy

!------------------------------------------------------------------------------

subroutine gemm(a, b, c, transa, transb, alpha, beta)

    external :: dgemm

    real(rp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: transa, transb
    real(rp), dimension(:,:), intent(in) :: a, b
    real(rp), dimension(:,:) , intent(inout) :: c

    integer(ip) :: m, n, k, lda, ldb, ldc
    character(len=1) :: o_transa, o_transb
    real(rp) :: o_alpha, o_beta

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0
    end if

    if (present(transa)) then
        o_transa = transa
    else
        o_transa = 'N'
    end if

    if (present(transb)) then
        o_transb = transb
    else
        o_transb = 'N'
    end if

    if (o_transa == 'N') then
        k = size(a, 2)
    else
        k = size(a, 1)
    end if

    m = size(c, 1)
    n = size(c, 2)
    lda = max(1, size(a, 1))
    ldb = max(1, size(b, 1))
    ldc = max(1, size(c, 1))

    call dgemm(o_transa, o_transb, m, n, k, o_alpha, a(1,1), lda, b(1,1), ldb, o_beta, c(1,1), ldc)

end subroutine gemm

!------------------------------------------------------------------------------

subroutine gemv(a, x, y, alpha, beta, trans)

    external :: dgemv

    real(rp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: trans
    real(rp), dimension(:,:), intent(in) :: a
    real(rp), dimension(:), intent(in) :: x
    real(rp), dimension(:), intent(inout) :: y

    integer(ip) :: m, n, incx, incy, lda
    character(len=1) :: o_trans
    real(rp) :: o_alpha, o_beta

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0
    end if

    if (present(trans)) then
        o_trans = trans
    else
        o_trans = 'N'
    end if

    incx = 1
    incy = 1
    lda = max(1, size(a, 1))
    m = size(a, 1)
    n = size(a, 2)

    call dgemv(o_trans, m, n, o_alpha, a(1,1), lda, x(1), incx, o_beta, y(1), incy)

end subroutine gemv

!------------------------------------------------------------------------------

subroutine spmv(ap, x, y, uplo, alpha, beta)

    external :: dspmv

    real(rp), dimension(:), intent(in) :: ap
    real(rp), dimension(:), intent(in) :: x
    real(rp), dimension(:), intent(inout) :: y
    character(len=1), intent(in), optional :: uplo
    real(rp), intent(in), optional :: alpha, beta

    integer(ip) :: n, incx, incy
    real(rp) :: o_alpha, o_beta
    character(len=1) :: o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0
    end if

    incx = 1
    incy = 1
    n = size(x)

    call dspmv(o_uplo, n, o_alpha, ap(1), x(1), incx, o_beta, y(1), incy)

end subroutine spmv

!------------------------------------------------------------------------------

end module pelib_blas_interfaces
