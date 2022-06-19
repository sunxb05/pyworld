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

module pelib_lapack_interfaces

    use pelib_precision

    implicit none

contains

!------------------------------------------------------------------------------

subroutine sptrf(ap, uplo, ipiv, info)

    external :: dsptrf, xerbla

    real(rp), dimension(:), intent(inout) :: ap
    character(len=1), optional :: uplo
    integer(ip), dimension(:), optional, target :: ipiv
    integer(ip), intent(out), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn
    integer(ip), dimension(:), pointer :: o_ipiv

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    if (present(ipiv)) then
        o_ipiv => ipiv
    else
        allocate(o_ipiv(n))
    end if

    call dsptrf(o_uplo, n, ap(1), o_ipiv(1), o_info)

    if (.not. present(ipiv)) then
        deallocate(o_ipiv)
    end if

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('SPTRF', - o_info)
    end if

end subroutine sptrf

!------------------------------------------------------------------------------

subroutine sptri(ap, ipiv, uplo, info)

    external :: dsptri, xerbla

    real(rp), dimension(:), intent(inout) :: ap
    integer(ip), dimension(:), intent(in) :: ipiv
    character(len=1), intent(in), optional :: uplo
    integer(ip), intent(out), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn
    real(rp), dimension(:), allocatable :: work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    allocate(work(n))

    call dsptri(o_uplo, n, ap(1), ipiv(1), work, o_info)

    deallocate(work)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('SPTRI', - o_info)
    end if

end subroutine sptri

!------------------------------------------------------------------------------

subroutine sptrs(ap, b, ipiv, uplo, info)

    external :: dsptrs, xerbla

    real(rp), dimension(:), intent(in) :: ap
    real(rp), dimension(:,:), intent(inout) :: b
    integer(ip), dimension(:), intent(in) :: ipiv
    character(len=1), optional :: uplo
    integer(ip), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn, nrhs, ldb

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)
    nrhs = size(b, 2)
    ldb = max(1, size(b, 1))

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    call dsptrs(o_uplo, n, nrhs, ap(1), ipiv(1), b(1,1), ldb, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('SPTRS', - o_info)
    end if

end subroutine sptrs

!------------------------------------------------------------------------------

subroutine pptrf(ap, uplo, info)

    external :: dpptrf, xerbla

    real(rp), dimension(:), intent(inout) :: ap
    character(len=1), optional :: uplo
    integer(ip), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    call dpptrf(o_uplo, n, ap(1), o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('PPTRF', - o_info)
    end if

end subroutine pptrf

!------------------------------------------------------------------------------

subroutine pptri(ap, uplo, info)

    external :: dpptri, xerbla

    real(rp), dimension(:), intent(inout) :: ap
    character(len=1), optional :: uplo
    integer(ip), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    call dpptri(o_uplo, n, ap(1), o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('PPTRI', - o_info)
    end if

end subroutine pptri

!------------------------------------------------------------------------------

subroutine pptrs(ap, b, uplo, info)

    external :: dpptrs, xerbla

    real(rp), dimension(:), intent(inout) :: ap
    real(rp), dimension(:,:), intent(inout) :: b
    character(len=1), optional :: uplo
    integer(ip), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn, nrhs, ldb

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)
    nrhs = size(b, 2)
    ldb = max(1, size(b, 1))

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    call dpptrs(o_uplo, n, nrhs, ap(1), b(1,1), ldb, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('PPTRS', - o_info)
    end if

end subroutine pptrs

!------------------------------------------------------------------------------

subroutine spcon(ap, ipiv, anorm, rcond, uplo, info)

    external :: dspcon, xerbla

    real(rp), intent(in) :: anorm
    real(rp), intent(out) :: rcond
    real(rp), dimension(:), intent(in) :: ap
    integer(ip), dimension(:), intent(in) :: ipiv
    character(len=1), optional :: uplo
    integer(ip), intent(out), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn
    real(rp), dimension(:), allocatable :: iwork, work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    allocate(iwork(n))
    allocate(work(2*n))

    call dspcon(o_uplo, n, ap(1), ipiv(1), anorm, rcond, work, iwork, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('DSPCON', - o_info)
    end if

end subroutine spcon

!------------------------------------------------------------------------------

subroutine ppcon(ap, anorm, rcond, uplo, info)

    external :: dppcon, xerbla

    real(rp), intent(in) :: anorm
    real(rp), intent(out) :: rcond
    real(rp), dimension(:), intent(in) :: ap
    character(len=1), optional :: uplo
    integer(ip), intent(out), optional :: info

    character(len=1) :: o_uplo
    integer(ip) :: o_info
    integer(ip) :: n, nn
    real(rp), dimension(:), allocatable :: iwork, work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    allocate(iwork(n))
    allocate(work(3*n))

    call dppcon(o_uplo, n, ap(1), anorm, rcond, work, iwork, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('DPPCON', - o_info)
    end if

end subroutine ppcon

!------------------------------------------------------------------------------

function lansp(norm, ap, uplo)

    real(rp), external :: dlansp

    real(rp) :: lansp
    real(rp), dimension(:), intent(in) :: ap
    character(len=1), intent(in) :: norm
    character(len=1), optional :: uplo

    character(len=1) :: o_uplo
    integer(ip) :: n, nn
    real(rp), dimension(:), allocatable :: work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1.0 + sqrt(1.0 + 8.0 * real(nn, rp))) * 0.5)
    end if

    allocate(work(n))

    lansp = dlansp(norm, o_uplo, n, ap(1), work)

    deallocate(work)

end function lansp

!------------------------------------------------------------------------------

end module pelib_lapack_interfaces
