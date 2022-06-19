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
module pelib_utils

    use pelib_precision

    implicit none

contains

!------------------------------------------------------------------------------

function elem2charge(elem) result(charge)

    use pelib_options, only: luout

    character(len=*), intent(in) :: elem

    integer(ip) :: i
    real(rp) :: charge
    character(len=2), dimension(109) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)

    if (elem == 'X') then
        charge = 0.0
        return
    end if

    do i = 1, 109
        if (elem == trim(elements(i))) then
            charge = real(i, rp)
            exit
        else
            charge = 0.0
        end if
    end do

    if (int(charge) == 0 .and. elem /= 'X') then
        write(luout, *) 'WARNING: charge not found for element: ', elem
    end if

end function elem2charge

!------------------------------------------------------------------------------

function charge2elem(charge) result(elem)

    real(rp), intent(in) :: charge

    character(len=2) :: elem
    character(len=2), dimension(109) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)

    if ((nint(charge) == 0) .or. (nint(charge) > 109)) then
        elem = 'X '
        return
    end if

    elem = elements(nint(charge))

end function charge2elem

!------------------------------------------------------------------------------

function charge2mass(charge) result(mass)

    use pelib_options, only: luout

    real(rp), intent(in) :: charge

    integer(ip) :: i
    real(rp) :: mass
    real(rp), dimension(109) :: masses

    masses = (/   1.007825,   4.002603,   7.016005,   9.012183,&
              &  11.009305,  12.000000,  14.003074,  15.994915,&
              &  18.998403,  19.992439,  22.989770,  23.985045,&
              &  26.981541,  27.976928,  30.973763,  31.972072,&
              &  34.968853,  39.962383,  38.963708,  39.962591,&
              &  44.955914,  47.947947,  50.943963,  51.940510,&
              &  54.938046,  55.934939,  58.933198,  57.935347,&
              &  62.929599,  63.929145,  68.925581,  73.921179,&
              &  74.921596,  79.916521,  78.918336,  83.911506,&
              &  84.911800,  87.905625,  88.905856,  89.904708,&
              &  92.906378,  97.905405,  97.907216, 101.904348,&
              & 102.905503, 105.903475, 106.905095, 113.903361,&
              & 114.903875, 119.902199, 120.903824, 129.906229,&
              & 126.904477, 131.904148, 132.905429, 137.905232,&
              & 138.906347, 139.905433, 140.907647, 141.907719,&
              & 144.912743, 151.919728, 152.921225, 157.924019,&
              & 158.925342, 163.929171, 164.930319, 165.930290,&
              & 168.934212, 173.938859, 174.940770, 179.946546,&
              & 180.947992, 183.950928, 186.955744, 191.961467,&
              & 192.962917, 194.964766, 196.966543, 201.970617,&
              & 204.974401, 207.976627, 208.980374, 208.982404,&
              & 209.987126, 222.017571, 223.019733, 226.025403,&
              & 227.027750, 232.038051, 231.035880, 238.050785,&
              & 237.048168, 244.064199, 243.061373, 247.070347,&
              & 247.070300, 251.079580, 252.082944, 257.095099,&
              & 258.098570, 259.100931, 260.105320, 261.108690,&
              & 262.113760, 263.118220, 262.122930, 269.134100,&
              & 267.138000 /)

    if (int(charge) == 0) then
        mass = 0.0
        return
    end if

    do i = 1, 109
        if (int(charge) == i) then
            mass = masses(i)
            exit
        else
            mass = 0.0
        end if
    end do

    if (int(mass) == 0 .and. int(charge) == 0) then
        write(luout, *) 'WARNING: mass not found for element with charge: ', charge
    end if

end function charge2mass

!------------------------------------------------------------------------------

subroutine chcase(string, uplo)

    character(len=*), intent(inout) :: string
    character(len=*), intent(in), optional :: uplo

    integer(ip) :: i, gap
    character(len=1) :: a, z, o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    end if

    gap = iachar('a') - iachar('A')

    if (o_uplo == 'u' .or. o_uplo == 'U') then
        a = 'a'
        z = 'z'
    else if (o_uplo == 'l' .or. o_uplo == 'L') then
        a = 'A'
        z = 'Z'
        gap = - gap
    else
        error stop 'Unknown case specified'
    end if

    do i = 1, len_trim(string)
        if (lge(string(i:i), a) .and. lle(string(i:i), z)) then
            string(i:i) = achar(iachar(string(i:i)) - gap)
        end if
    end do

end subroutine chcase

!------------------------------------------------------------------------------

subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer(ip), intent(out) :: lunit
    integer(ip) :: i
    logical(lp) :: lexist, lopen

    if (stat == 'old') then
        inquire(file=filename, exist=lexist)

        if (.not. lexist) then
            print *, filename, ' not found.'
            error stop 'file not found.'
        end if
    end if

    do i = 21, 99
        inquire(unit=i, opened=lopen)
        if (lopen) then
            cycle
        else
            lunit = i
            open(unit=lunit, file=filename, status=stat, form=frmt)
            exit
        end if
    end do

    return

end subroutine openfile

!------------------------------------------------------------------------------

end module pelib_utils
