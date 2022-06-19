  ! spacing of finite difference
  real(REALK), parameter :: FD_SPACING = 10.0_REALK**(-8)
  ! maximum order to test
  integer, parameter :: FD_MAX_ORDER = 30
  ! order of electronic derivatives
  integer, parameter :: FD_ORDER_ELEC = 0
  ! order of Cartesian multipole moments
  integer, parameter :: FD_ORDER_MOM = 1
  ! default order which is not to test
  integer, parameter :: FD_DEFAULT = 2
