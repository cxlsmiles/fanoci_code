module globals
    implicit none

!------------------------------
!
! Input data
!
!------------------------------
    character(len=100) :: prop_type
    character(len=100) :: gam
    character(len=100) :: hs_in
    character(len=100) :: ps_in
    integer :: h_in
    integer :: p_in
    integer :: h_f_min
    integer :: h_f_max
    integer :: p_f_min
    integer :: p_f_max
    integer :: n_occ
    integer :: prnt
    integer :: moint
    double precision :: eig_thresh
    
    integer :: n_hs_in, n_ps_in
    integer, dimension(:), allocatable :: holes_arr
    integer, dimension(:), allocatable :: particles_arr
	
!------------------------------
!
!  CI configurations
!
!------------------------------
    integer :: n_ov
    integer :: n_config
    integer, dimension(:,:), allocatable :: config_2h
    integer, dimension(:,:), allocatable :: config_1h1p
    character(len=20), dimension(:), allocatable :: config_2hnp_str
    character(len=20), dimension(:), allocatable :: config_2h_str
    double precision, dimension(:), allocatable :: coupling_arr_el
    double precision, dimension(:,:), allocatable :: coupling_arr_pol

!------------------------------
!
! Hartree-Fock data
!
!------------------------------
    double precision :: E_0
    integer*8 :: n_int2e
    integer :: n_mo
    double precision, dimension(:), allocatable :: int2e
    double precision, dimension(:), allocatable :: emo
    double precision, dimension(:,:), allocatable :: s_mat

!------------------------------
!
! CI matrix
!
!------------------------------
    double precision, dimension(:,:), allocatable :: CI_init_matrix
    double precision, dimension(:,:), allocatable :: CI_matrix
    double precision, dimension(:,:), allocatable :: CI_eigvec
    
    double precision, dimension(:), allocatable :: CI_init_eigval
    double precision, dimension(:), allocatable :: CI_eigval
    integer :: ci_mat_in_size
    integer :: ci_mat_size

    !...................................
    ! 2h CI matrix, partial decay widths
    !...................................
    double precision, dimension (:,:), allocatable :: mat_2h
    double precision, dimension (:,:), allocatable :: mat_2h_s
    double precision, dimension (:,:), allocatable :: mat_2h_t
    double precision, dimension (:), allocatable :: mat_2h_eigval
    double precision, dimension (:), allocatable :: mat_2h_s_eigval
    double precision, dimension (:), allocatable :: mat_2h_t_eigval
    integer :: ci_mat_2h_s_size
    integer :: ci_mat_2h_t_size

end module globals
