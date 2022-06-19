  ! created by GeomPartZero.py on 2011-11-07
  ! number of total geometric derivatives, according to
  ! (1) numbers of differentiated centers
  ! (2) indices of differentiated centers
  ! (3) different indices of bra and ket centers
  ! (4) orders of partial geoemtric derivatives on bra and ket centers
  integer, parameter :: NUM_TOT_GEO = 28
  ! numbers of differentiated centers
  integer, parameter :: NUM_CENTS(NUM_TOT_GEO) = (/                      &
    1,       1,       1,       1,       1,       1,       1,       1,    &
    1,       1,       1,       1,       1,       1,       1,       1,    &
    2,       2,       2,       2,       2,       2,       2,       2,    &
    2,       2,       2,       2/)
  ! indices of differentiated centers
  integer, parameter :: IDX_GEO_CENT(NUM_TOT_GEO*2) = (/                 &
    2, 0,    2, 0,    2, 0,    2, 0,    2, 0,    2, 0,    2, 0,    2, 0, &
    2, 0,    2, 0,    2, 0,    2, 0,    2, 0,    2, 0,    2, 0,    2, 0, &
    1, 2,    1, 2,    1, 2,    1, 2,    1, 2,    1, 2,    1, 2,    1, 2, &
    1, 2,    1, 2,    1, 2,    1, 2/)
  ! orders of differentiated centers
  integer, parameter :: ORDER_GEO_CENT(NUM_TOT_GEO*2) = (/               &
    3, 0,    3, 0,    3, 0,    3, 0,    3, 0,    3, 0,    3, 0,    3, 0, &
    3, 0,    3, 0,    3, 0,    3, 0,    3, 0,    3, 0,    3, 0,    3, 0, &
    3, 4,    3, 4,    3, 4,    3, 4,    3, 4,    3, 4,    3, 4,    3, 4, &
    3, 4,    3, 4,    3, 4,    3, 4/)
  ! indices of bra and ket centers
  integer, parameter :: IDX_BRA_KET(NUM_TOT_GEO*2) = (/                  &
    7, 8,    7, 8,    7, 8,    7, 8,    2, 8,    2, 8,    2, 8,    2, 8, &
    7, 2,    7, 2,    7, 2,    7, 2,    2, 2,    2, 2,    2, 2,    2, 2, &
    2, 7,    2, 7,    2, 7,    2, 7,    1, 2,    1, 2,    1, 2,    1, 2, &
    2, 1,    2, 1,    2, 1,    2, 1/)
  ! orders of partial geoemtric derivatives on bra center
  integer, parameter :: ORDER_GEO_BRA(NUM_TOT_GEO) = (/                  &
    0,       3,       0,       3,       0,       3,       0,       3,    &
    0,       3,       0,       3,       0,       3,       0,       3,    &
    0,       3,       0,       3,       0,       3,       0,       3,    &
    0,       3,       0,       3/)
  ! orders of partial geoemtric derivatives on ket center
  integer, parameter :: ORDER_GEO_KET(NUM_TOT_GEO) = (/                  &
    0,       0,       3,       3,       0,       0,       3,       3,    &
    0,       0,       3,       3,       0,       0,       3,       3,    &
    0,       0,       3,       3,       0,       0,       3,       3,    &
    0,       0,       3,       3/)
  ! referenced results indicating if the total geometric derivatives are zero
  logical, parameter :: REF_ZERO_INTS(NUM_TOT_GEO) = (/                     &
    .true.,  .true.,  .true.,  .true.,  .false., .false., .false., .false., &
    .false., .false., .false., .false., .true.,  .true.,  .true.,  .true.,  &
    .true.,  .true.,  .true.,  .true.,  .false., .false., .false., .false., &
    .false., .false., .false., .false./)
  ! referenced final orders of partial geometric derivatives by adding
  ! those from total geometric derivatives
  integer, parameter :: REF_ORDER_PART(NUM_TOT_GEO*2) = (/               &
    0, 0,    0, 0,    0, 0,    0, 0,    3, 0,    6, 0,    3, 3,    6, 3, &
    0, 3,    3, 3,    0, 6,    3, 6,    0, 0,    0, 0,    0, 0,    0, 0, &
    0, 0,    0, 0,    0, 0,    0, 0,    3, 4,    6, 4,    3, 7,    6, 7, &
    4, 3,    7, 3,    4, 6,    7, 6/)
  ! referenced results indicating if scattering the geometric derivatives later on
  logical, parameter :: REF_SCATTER_DER(NUM_TOT_GEO) = (/                   &
    .false., .false., .false., .false., .false., .true.,  .true.,  .true.,  &
    .false., .false., .true.,  .true.,  .false., .false., .false., .false., &
    .false., .false., .false., .false., .false., .true.,  .true.,  .true.,  &
    .true.,  .true.,  .true.,  .true./)
  ! referenced sequences of bra and ket centers for partial derivative terms
  integer, parameter :: REF_SEQ_PART(NUM_TOT_GEO*2) = (/                 &
    0, 0,    0, 0,    0, 0,    0, 0,    1, 0,    1, 0,    1, 0,    1, 0, &
    2, 0,    2, 0,    2, 0,    2, 0,    0, 0,    0, 0,    0, 0,    0, 0, &
    0, 0,    0, 0,    0, 0,    0, 0,    1, 2,    1, 2,    1, 2,    1, 2, &
    2, 1,    2, 1,    2, 1,    2, 1/)
