MODULE cosp_kinds
  use machine,  only: kind_phys,kind_dbl_prec,kind_sngl_prec
  implicit none
  integer,parameter :: sp = kind_sngl_prec
  integer,parameter :: dp = kind_dbl_prec
  integer,parameter :: wp = kind_phys
END MODULE cosp_kinds
