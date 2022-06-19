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

module pelib_precision

    use, intrinsic :: iso_fortran_env, only: real32, real64, real128, &
                                             int32, int64

    implicit none

#if INT32
    integer(int32), parameter :: ip = int32
#elif INT64
    integer(int64), parameter :: ip = int64
#else
    integer(int32), parameter :: ip = int32
#endif

#ifdef REAL32
    integer(ip), parameter :: rp = real32
#elif REAL64
    integer(ip), parameter :: rp = real64
#elif REAL128
    integer(ip), parameter :: rp = real128
#else
    integer(ip), parameter :: rp = real64
#endif

    integer(ip), parameter :: lp = int32

end module pelib_precision
