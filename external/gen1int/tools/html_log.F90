!!  gen1int: compute one-electron integrals using rotational London atomic-orbitals
!!  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
!!
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file contains routines of accessing an HTML test log.
!!
!!  2012-02-04, Bin Gao
!!  * first version

! maximum length of a string
#include "max_len_str.h"

!> \brief module of HTML test log routines
!> \author Bin Gao
!> \date 2012-02-05
module html_log

  !  kind type parameter of real numbers
  use xkind
  implicit none

  public :: html_log_begin
  public :: html_log_end
  public :: html_log_head_begin
  public :: html_log_head_end
  public :: html_log_title
  public :: html_log_body_begin
  public :: html_log_body_end
  public :: html_log_heading
  public :: html_log_begin_paragraph
  public :: html_log_end_paragraph
  public :: html_log_line_break
  public :: html_log_text
  public :: html_log_href
  public :: html_log_int_number
  public :: html_log_real_number
  public :: html_log_int_array
  public :: html_log_real_array
  public :: html_log_horizontal_line
  public :: html_log_date_and_time
  public :: html_log_table_begin
  public :: html_log_table_header
  public :: html_log_table_row
  public :: html_log_table_end

  interface html_log_table_row
    module procedure html_log_table_int_row
    module procedure html_log_table_real_row
  end interface html_log_table_row

  contains

  !> \brief begins test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_begin(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<html lang=""en"">"
  end subroutine html_log_begin

  !> \brief ends test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_end(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "</html>"
  end subroutine html_log_end

  !> \brief begins the head of test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_head_begin(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<head>"
  end subroutine html_log_head_begin

  !> \brief ends the head of test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_head_end(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "</head>"
  end subroutine html_log_head_end

  !> \brief writes the title of test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param log_title is the log title
  !> \param io_log is the IO unit of log file
  subroutine html_log_title(log_title, io_log)
    character*(*), intent(in) :: log_title
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<title>"//trim(log_title)//"</title>"
  end subroutine html_log_title

  !> \brief begins the body of test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_body_begin(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<body>"
  end subroutine html_log_body_begin

  !> \brief ends the body of test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_body_end(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "</body>"
  end subroutine html_log_body_end

  !> \brief writes a heading in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param heading is the heading text
  !> \param io_log is the IO unit of log file
  !> \param level is the level of heading text
  subroutine html_log_heading(heading, io_log, level)
    character*(*), intent(in) :: heading
    integer, intent(in) :: io_log
    integer, optional, intent(in) :: level
! digit from 0 to 9
#include "digit.h"
    integer p_level  !level of heading text, from 1 to 6 (private)
    if (present(level)) then
      p_level = max(1,min(6,level))
    else
      p_level = 6
    end if
    write(io_log,"(A)") "<h"//DIGIT(p_level)//">"//trim(heading)// &
                        "</h"//DIGIT(p_level)//">"
  end subroutine html_log_heading

  !> \brief begins a paragraph in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_begin_paragraph(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<p>"
  end subroutine html_log_begin_paragraph

  !> \brief ends a paragraph in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_end_paragraph(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "</p>"
  end subroutine html_log_end_paragraph

  !> \brief creates a line break in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_line_break(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<br>"
  end subroutine html_log_line_break

  !> \brief writes text in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param log_text contains the text
  !> \param io_log is the IO unit of log file
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_text(log_text, io_log, font_size, font_color)
    character*(*), intent(in) :: log_text
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    if (present(font_size)) then
      if (present(font_color)) then
        write(io_log,"(A)") "<font size="""//font_size//""" color="""// &
                            font_color//""">"//log_text//"</font><br>"
      else
        write(io_log,"(A)") "<font size="""//font_size//""">"// &
                            log_text//"</font><br>"
      end if
    else
      if (present(font_color)) then
        write(io_log,"(A)") "<font color="""//font_color//""">"// &
                            log_text//"</font><br>"
      else
        write(io_log,"(A)") log_text//"<br>"
      end if
    end if
  end subroutine html_log_text

  !> \brief specifies a URL in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param url is the URL
  !> \param log_text contains the text
  !> \param io_log is the IO unit of log file
  subroutine html_log_href(log_text, url, url_text, io_log)
    character*(*), intent(in) :: log_text
    character*(*), intent(in) :: url
    character*(*), intent(in) :: url_text
    integer, intent(in) :: io_log
    write(io_log,"(A)") log_text//" <a href="""//url//""">"//url_text//"</a>.<br>"
  end subroutine html_log_href

  !> \brief writes text and an integer in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param log_text contains the text
  !> \param log_int is the integer
  !> \param fmt_int is the format of integer
  !> \param io_log is the IO unit of log file
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_int_number(log_text, log_int, fmt_int, io_log, &
                                 font_size, font_color)
    character*(*), intent(in) :: log_text
    integer, intent(in) :: log_int
    character*(*), intent(in) :: fmt_int
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    if (present(font_size)) then
      if (present(font_color)) then
        write(io_log,"(A,"//fmt_int//",A)")                              &
          "<font size="""//font_size//""" color="""//font_color//""">"// &
          log_text, log_int, "</font><br>"
      else
        write(io_log,"(A,"//fmt_int//",A)") &
          "<font size="""//font_size//""">"//log_text, log_int, "</font><br>"
      end if
    else
      if (present(font_color)) then
        write(io_log,"(A,"//fmt_int//",A)") &
          "<font color="""//font_color//""">"//log_text, log_int, "</font><br>"
      else
        write(io_log,"(A,"//fmt_int//",A)") log_text, log_int, "<br>"
      end if
    end if
  end subroutine html_log_int_number

  !> \brief writes text and a real number in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param log_text contains the text
  !> \param log_real is the real number
  !> \param fmt_real is the format of real number
  !> \param io_log is the IO unit of log file
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_real_number(log_text, log_real, fmt_real, io_log, &
                                  font_size, font_color)
    character*(*), intent(in) :: log_text
    real(REALK), intent(in) :: log_real
    character*(*), intent(in) :: fmt_real
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    if (present(font_size)) then
      if (present(font_color)) then
        write(io_log,"(A,"//fmt_real//",A)")                             &
          "<font size="""//font_size//""" color="""//font_color//""">"// &
          log_text, log_real, "</font><br>"
      else
        write(io_log,"(A,"//fmt_real//",A)") &
          "<font size="""//font_size//""">"//log_text, log_real, "</font><br>"
      end if
    else
      if (present(font_color)) then
        write(io_log,"(A,"//fmt_real//",A)") &
          "<font color="""//font_color//""">"//log_text, log_real, "</font><br>"
      else
        write(io_log,"(A,"//fmt_real//",A)") log_text, log_real, "<br>"
      end if
    end if
  end subroutine html_log_real_number

  !> \brief writes text and integers in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param log_text contains the text
  !> \param log_int contains the integers
  !> \param fmt_int is the format of integers
  !> \param io_log is the IO unit of log file
  !> \param separation is the separation sign for integers
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_int_array(log_text, log_int, fmt_int, io_log, &
                                separation, font_size, font_color)
    character*(*), intent(in) :: log_text
    integer, intent(in) :: log_int(:)
    character*(*), intent(in) :: fmt_int
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: separation
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    integer size_int                     !size of integers
    character(MAX_LEN_STR) fmt_array     !format of the array
    character(MAX_LEN_STR) p_separation  !separation sign for integers (private)
    ! gets the size of integers
    size_int = size(log_int)
    if (size_int==1) then
      if (len_trim(fmt_int)>MAX_LEN_STR)                               &
        call error_stop("html_log_int_array",                          &
                        "increases length of string in max_len_str.h", &
                        len_trim(fmt_int))
      fmt_array = fmt_int
    else
      ! sets the separation sign for integers
      if (present(separation)) then
        if (len_trim(separation)>MAX_LEN_STR)                            &
          call error_stop("html_log_int_array",                          &
                          "increases length of string in max_len_str.h", &
                          len_trim(separation))
        p_separation = separation
      else
        p_separation = ","
      end if
      ! converts \var(size_int)-1 to string
      call int_to_str(size_int-1, MAX_LEN_STR, .true., .false., fmt_array)
      ! checks the length of \var(fmt_array)
      size_int = len_trim(fmt_array)
      if (size_int+2*len_trim(fmt_int)+len_trim(p_separation)+6>MAX_LEN_STR) &
        call error_stop("html_log_int_array",                                &
                        "increases length of string in max_len_str.h",       &
                        size_int+2*len_trim(fmt_int)+len_trim(p_separation)+6)
      fmt_array(size_int+1:) = "("//fmt_int//","""//trim(p_separation)//"""),"//fmt_int
    end if
    if (present(font_size)) then
      if (present(font_color)) then
        write(io_log,"(A,"//trim(fmt_array)//",A)")   &
          "<font size="""//font_size//""" color="""// &
          font_color//""">"//log_text//" (", log_int, ")</font><br>"
      else
        write(io_log,"(A,"//trim(fmt_array)//",A)") &
          "<font size="""//font_size//""">"//log_text//" (", log_int, ")</font><br>"
      end if
    else
      if (present(font_color)) then
        write(io_log,"(A,"//trim(fmt_array)//",A)") &
          "<font color="""//font_color//""">"//log_text//" (", log_int, ")</font><br>"
      else
        write(io_log,"(A,"//trim(fmt_array)//",A)") &
          log_text//" (", log_int, ")<br>"
      end if
    end if
  end subroutine html_log_int_array

  !> \brief writes text and real numbers in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param log_text contains the text
  !> \param log_real contains the real numbers
  !> \param fmt_real is the format of real numbers
  !> \param io_log is the IO unit of log file
  !> \param separation is the separation sign for integers
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_real_array(log_text, log_real, fmt_real, io_log, &
                                 separation, font_size, font_color)
    character*(*), intent(in) :: log_text
    real(REALK), intent(in) :: log_real(:)
    character*(*), intent(in) :: fmt_real
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: separation
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    integer size_real                    !size of real numbers
    character(MAX_LEN_STR) fmt_array     !format of the array
    character(MAX_LEN_STR) p_separation  !separation sign for integers (private)
    ! gets the number of integers
    size_real = size(log_real)
    if (size_real==1) then
      if (len_trim(fmt_real)>MAX_LEN_STR)                              &
        call error_stop("html_log_real_array",                         &
                        "increases length of string in max_len_str.h", &
                        len_trim(fmt_real))
      fmt_array = fmt_real
    else
      ! sets the separation sign for integers
      if (present(separation)) then
        if (len_trim(separation)>MAX_LEN_STR)                            &
          call error_stop("html_log_real_array",                         &
                          "increases length of string in max_len_str.h", &
                          len_trim(separation))
        p_separation = separation
      else
        p_separation = ","
      end if
      ! converts \var(size_real)-1 to string
      call int_to_str(size_real-1, MAX_LEN_STR, .true., .false., fmt_array)
      ! checks the length of \var(fmt_array)
      size_real = len_trim(fmt_array)
      if (size_real+2*len_trim(fmt_real)+len_trim(p_separation)+6>MAX_LEN_STR) &
        call error_stop("html_log_real_array",                                 &
                        "increases length of string in max_len_str.h",         &
                        size_real+2*len_trim(fmt_real)+len_trim(p_separation)+6)
      fmt_array(size_real+1:) = "("//fmt_real//","""//trim(p_separation)//"""),"//fmt_real
    end if
    if (present(font_size)) then
      if (present(font_color)) then
        write(io_log,"(A,"//trim(fmt_array)//",A)")   &
          "<font size="""//font_size//""" color="""// &
          font_color//""">"//log_text//" (", log_real, ")</font><br>"
      else
        write(io_log,"(A,"//trim(fmt_array)//",A)") &
          "<font size="""//font_size//""">"//log_text//" (", log_real, ")</font><br>"
      end if
    else
      if (present(font_color)) then
        write(io_log,"(A,"//trim(fmt_array)//",A)") &
          "<font color="""//font_color//""">"//log_text//" (", log_real, ")</font><br>"
      else
        write(io_log,"(A,"//trim(fmt_array)//",A)") &
          log_text//" (", log_real, ")<br>"
      end if
    end if
  end subroutine html_log_real_array

  !> \brief creates a horizontal line in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_horizontal_line(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<hr>"
  end subroutine html_log_horizontal_line

  !> \brief writes current date and time in test log
  !> \author Bin Gao
  !> \date 2012-02-04
  !> \param io_log is the IO unit of log file
  subroutine html_log_date_and_time(io_log)
    integer, intent(in) :: io_log
    character(8) curr_date   !\var(curr_date) and \var(curr_time) are used for
    character(10) curr_time  !\fn(date_and_time) to get current time and date
    call date_and_time(DATE=curr_date, TIME=curr_time)
    write(io_log,"(A)") "Current date and time: "//                    &
      curr_time(1:2)//":"//curr_time(3:4)//":"//curr_time(5:6)//", "// &
      curr_date(1:4)//"-"//curr_date(5:6)//"-"//curr_date(7:8)//"<br>"
  end subroutine html_log_date_and_time

  !> \brief begins a table in test log
  !> \author Bin Gao
  !> \date 2012-02-07
  !> \param io_log is the IO unit of log file
  subroutine html_log_table_begin(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "<table>"
  end subroutine html_log_table_begin

  !> \brief writes the header cells of a table
  !> \author Bin Gao
  !> \date 2012-02-07
  !> \param header_text contains the text of header cells
  !> \param io_log is the IO unit of log file
  subroutine html_log_table_header(header_text, io_log)
    character*(*), intent(in) :: header_text(:)
    integer, intent(in) :: io_log
    integer num_cells                 !number of header cells
    character(MAX_LEN_STR) fmt_cell  !format of the header cells
    integer ishell                    !incremental recoder over shells
    ! gets the number of header cells
    num_cells = size(header_text)
    ! sets the format of the header cells
    if (num_cells==1) then
      fmt_cell = "(4X,""<th>"",A,""</th>"")"
    else
      fmt_cell(1:4) = "(4X,"
      ! converts \var(num_cells) to string
      call int_to_str(num_cells, MAX_LEN_STR-4, .true., .false., &
                      fmt_cell(5:MAX_LEN_STR))
      ! checks the length of \var(fmt_cell)
      ishell = len_trim(fmt_cell)
      if (ishell+19>MAX_LEN_STR)                                       &
        call error_stop("html_log_table_header",                       &
                        "increases length of string in max_len_str.h", &
                        ishell+19)
      fmt_cell(ishell+1:ishell+19) = "(""<th>"",A,""</th>""))"
    end if
    write(io_log,"(2X,A)") "<tr>"
    write(io_log,trim(fmt_cell)) (trim(header_text(ishell)), ishell=1,num_cells)
    write(io_log,"(2X,A)") "</tr>"
  end subroutine html_log_table_header

  !> \brief writes the standard cells (integers) of a table
  !> \author Bin Gao
  !> \date 2012-02-07
  !> \param cell_int contains the integers of cells
  !> \param num_int contains the number of integers in each cell
  !> \param fmt_int contains format of integers in each cell
  !> \param io_log is the IO unit of log file
  !> \param separation is the separation sign for integers
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_table_int_row(cell_int, num_int, fmt_int, io_log, &
                                    separation, font_size, font_color)
    integer, intent(in) :: cell_int(:)
    integer, intent(in) :: num_int(:)
    character*(*), intent(in) :: fmt_int(:)
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: separation
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    character(MAX_LEN_STR) p_separation  !separation sign for integers (private)
    character(MAX_LEN_STR) p_font_size   !size of text (private)
    character(MAX_LEN_STR) p_font_color  !color of text (private)
    integer num_cells                    !number of standard cells
    integer ishell                       !incremental recoder over shells
    character(MAX_LEN_STR) fmt_cell      !format of the standard cells
    integer len_fmt_int, len_fmt_cell    !lengths of the format string
    integer start_int, end_int           !start and end of the integers in a cell
    ! sets the separation sign for integers
    if (present(separation)) then
      if (len_trim(separation)>MAX_LEN_STR)                            &
        call error_stop("html_log_table_int_row",                      &
                        "increases length of string in max_len_str.h", &
                        len_trim(separation))
      p_separation = separation
    else
      p_separation = ","
    end if
    ! sets the size of text
    if (present(font_size)) then
      if (len_trim(font_size)>MAX_LEN_STR)                             &
        call error_stop("html_log_table_int_row",                      &
                        "increases length of string in max_len_str.h", & 
                        len_trim(font_size))
      p_font_size = font_size
    else
      p_font_size = "4"
    end if
    ! sets the color of text
    if (present(font_color)) then
      if (len_trim(font_color)>MAX_LEN_STR)                            &
        call error_stop("html_log_table_int_row",                      &
                        "increases length of string in max_len_str.h", & 
                        len_trim(font_color))
      p_font_color = font_color
    else
      p_font_color = "black"
    end if
    ! gets the number of standard cells
    num_cells = size(num_int)
    if (num_cells/=size(fmt_int))               &
      call error_stop("html_log_table_int_row", &
                      "incorrect number of standard cells", size(fmt_int))
    ! checks the size of \var(cell_int)
    if (sum(num_int)/=size(cell_int))           &
      call error_stop("html_log_table_int_row", &
                      "incorrect number of integers", size(cell_int))
    write(io_log,"(2X,A)") "<tr>"
    ! loops over the standard cells
    start_int = 1
    do ishell = 1, num_cells
      ! sets the end of the integers in current cell
      end_int = start_int+num_int(ishell)-1
      ! gets the length of format string of integers in current cell
      len_fmt_int = len_trim(fmt_int(ishell))
      if (num_int(ishell)==1) then
        if (len_fmt_int>MAX_LEN_STR)                                     &
          call error_stop("html_log_table_int_row",                      &
                          "increases length of string in max_len_str.h", &
                          len_fmt_int)
        fmt_cell = fmt_int(ishell)
      else
        ! converts \var(num_int(ishell))-1 to string
        call int_to_str(num_int(ishell)-1, MAX_LEN_STR, .true., .false., fmt_cell)
        ! checks the length of \var(fmt_array)
        len_fmt_cell = len_trim(fmt_cell)
        if (len_fmt_cell+2*len_fmt_int+len_trim(p_separation)+6>MAX_LEN_STR) &
          call error_stop("html_log_table_int_row",                          &
                          "increases length of string in max_len_str.h",     &
                          len_fmt_cell+2*len_fmt_int+len_trim(p_separation)+6)
        fmt_cell(len_fmt_cell+1:) = "("//fmt_int(ishell)(1:len_fmt_int)// &
                                    ","""//trim(p_separation)//"""),"//   &
                                    fmt_int(ishell)(1:len_fmt_int)
      end if
      ! writes the current cell
      write(io_log,"(4X,A,"//trim(fmt_cell)//",A)")             &
        "<td><font size="""//trim(p_font_size)//""" color="""// &
        trim(p_font_color)//""">", cell_int(start_int:end_int), "</font></td>"
      ! sets the start of the integers in next cell
      start_int = end_int+1
    end do
    write(io_log,"(2X,A)") "</tr>"
  end subroutine html_log_table_int_row

  !> \brief writes the standard cells (real numbers) of a table
  !> \author Bin Gao
  !> \date 2012-02-07
  !> \param cell_real contains the real numbers of cells
  !> \param num_real contains the number of real numbers in each cell
  !> \param fmt_real contains format of real numbers in each cell
  !> \param io_log is the IO unit of log file
  !> \param separation is the separation sign for real numbers
  !> \param font_size specifies the size of text
  !> \param font_color specifies the color of text
  subroutine html_log_table_real_row(cell_real, num_real, fmt_real, io_log, &
                                     separation, font_size, font_color)
    real(REALK), intent(in) :: cell_real(:)
    integer, intent(in) :: num_real(:)
    character*(*), intent(in) :: fmt_real(:)
    integer, intent(in) :: io_log
    character*(*), optional, intent(in) :: separation
    character*(*), optional, intent(in) :: font_size
    character*(*), optional, intent(in) :: font_color
    character(MAX_LEN_STR) p_separation  !separation sign for real numbers (private)
    character(MAX_LEN_STR) p_font_size   !size of text (private)
    character(MAX_LEN_STR) p_font_color  !color of text (private)
    integer num_cells                    !number of standard cells
    integer ishell                       !incremental recoder over shells
    character(MAX_LEN_STR) fmt_cell      !format of the standard cells
    integer len_fmt_real, len_fmt_cell   !lengths of the format string
    integer start_real, end_real         !start and end of the real numbers in a cell
    ! sets the separation sign for real numbers
    if (present(separation)) then
      if (len_trim(separation)>MAX_LEN_STR)                            &
        call error_stop("html_log_table_real_row",                     &
                        "increases length of string in max_len_str.h", &
                        len_trim(separation))
      p_separation = separation
    else
      p_separation = ","
    end if
    ! sets the size of text
    if (present(font_size)) then
      if (len_trim(font_size)>MAX_LEN_STR)                             &
        call error_stop("html_log_table_real_row",                     &
                        "increases length of string in max_len_str.h", & 
                        len_trim(font_size))
      p_font_size = font_size
    else
      p_font_size = "4"
    end if
    ! sets the color of text
    if (present(font_color)) then
      if (len_trim(font_color)>MAX_LEN_STR)                            &
        call error_stop("html_log_table_real_row",                     &
                        "increases length of string in max_len_str.h", & 
                        len_trim(font_color))
      p_font_color = font_color
    else
      p_font_color = "black"
    end if
    ! gets the number of standard cells
    num_cells = size(num_real)
    if (num_cells/=size(fmt_real))               &
      call error_stop("html_log_table_real_row", &
                      "incorrect number of standard cells", size(fmt_real))
    ! checks the size of \var(cell_real)
    if (sum(num_real)/=size(cell_real))          &
      call error_stop("html_log_table_real_row", &
                      "incorrect number of real numbers", size(cell_real))
    write(io_log,"(2X,A)") "<tr>"
    ! loops over the standard cells
    start_real = 1
    do ishell = 1, num_cells
      ! sets the end of the real numbers in current cell
      end_real = start_real+num_real(ishell)-1
      ! gets the length of format string of real numbers in current cell
      len_fmt_real = len_trim(fmt_real(ishell))
      if (num_real(ishell)==1) then
        if (len_fmt_real>MAX_LEN_STR)                                    &
          call error_stop("html_log_table_real_row",                     &
                          "increases length of string in max_len_str.h", &
                          len_fmt_real)
        fmt_cell = fmt_real(ishell)
      else
        ! converts \var(num_real(ishell))-1 to string
        call int_to_str(num_real(ishell)-1, MAX_LEN_STR, .true., .false., fmt_cell)
        ! checks the length of \var(fmt_array)
        len_fmt_cell = len_trim(fmt_cell)
        if (len_fmt_cell+2*len_fmt_real+len_trim(p_separation)+6>MAX_LEN_STR) &
          call error_stop("html_log_table_real_row",                          &
                          "increases length of string in max_len_str.h",      &
                          len_fmt_cell+2*len_fmt_real+len_trim(p_separation)+6)
        fmt_cell(len_fmt_cell+1:) = "("//fmt_real(ishell)(1:len_fmt_real)// &
                                    ","""//trim(p_separation)//"""),"//     &
                                    fmt_real(ishell)(1:len_fmt_real)
      end if
      ! writes the current cell
      write(io_log,"(4X,A,"//trim(fmt_cell)//",A)")             &
        "<td><font size="""//trim(p_font_size)//""" color="""// &
        trim(p_font_color)//""">", cell_real(start_real:end_real), "</font></td>"
      ! sets the start of the real numbers in next cell
      start_real = end_real+1
    end do
    write(io_log,"(2X,A)") "</tr>"
  end subroutine html_log_table_real_row

  !> \brief ends a table in test log
  !> \author Bin Gao
  !> \date 2012-02-07
  !> \param io_log is the IO unit of log file
  subroutine html_log_table_end(io_log)
    integer, intent(in) :: io_log
    write(io_log,"(A)") "</table>"
  end subroutine html_log_table_end

end module html_log
