module ci_blocks_el
    use globals
!    use misc
    use mat_element_el
    implicit none
    contains
    ! the following functions needed to construct the CI matrix
    !
    !     main blocks:
    ! function main_2h1p (csf_type, a, nrow, ncol) result(mat)
    !
    !     coupling blocks:
    ! function coupling_2h1p (a, type_1, type_2, nrow, ncol) result(mat)
    !

    ! =================
    !    MAIN BLOCK
    ! =================

    ! -----------------
    !   2h1p
    ! -----------------

    function main_2h1p (csf_type, a, nrow, ncol) result(mat)
        integer :: a
        integer :: ov1, ov2, ov3, ov4, row, col
        integer :: i, j
		integer :: nrow, ncol
		integer :: csf_type
        double precision, dimension(nrow, ncol) :: mat
        

        if (csf_type == 1) then !singlet I
            do i = 1, n_config
                ov1 = config_2h(i, 1)
                ov2 = config_2h(i, 2)

                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)

                    mat(i, j) = mat_el_2h1p_sI (ov1, ov2, a, ov3, ov4, a)

                end do

            end do
            

        else if (csf_type == 2) then ! singlet II
            do i = h_f_min, h_f_max-1
                row = i - h_f_min + 1
                mat(row,row) = mat_el_2h1p_sII (i, a, i, a)

                do j = i + 1, h_f_max
                    col = j - h_f_min + 1
                    mat(row,col) = mat_el_2h1p_sII (i, a, j, a)
                    mat(col,row) = mat(row,col)
                end do
            end do

            mat(n_ov,n_ov) = mat_el_2h1p_sII (h_f_max, a, h_f_max, a)


        else if (csf_type == 3) then ! triplet
            do i = 1, n_config-1
                ov1 = config_2h(i, 1)
                ov2 = config_2h(i, 2)

                mat(i, i) = mat_el_2h1p_t (ov1, ov2, a, ov1, ov2, a)

                do j = i+1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)

                    mat(i, j) = mat_el_2h1p_t (ov1, ov2, a, ov3, ov4, a)
                    mat(j, i) = mat(i, j)

                end do

            end do

            ov1 = config_2h(n_config,1)
            ov2 = config_2h(n_config,2)

            mat(n_config, n_config) = mat_el_2h1p_t (ov1, ov2, a, ov1, ov2, a)

        end if

    end function  main_2h1p

    
    ! ====================
    !    COUPLING BLOCK
    ! ====================
    function coupling_2h1p (a, type_1, type_2, nrow, ncol) result(mat)
        integer :: a
        integer :: ov1, ov2, ov3, ov4
        integer :: i, j, row
		integer :: nrow, ncol
		integer :: type_1, type_2
        double precision, dimension(nrow, ncol) :: mat

        if ((type_1 == 1) .and. (type_2 == 2)) then
            do i = 1, n_config
                ov1 = config_2h(i, 1)
                ov2 = config_2h(i, 2)

                do j = 1, n_ov
                    ov3 = j + h_f_min - 1
                    mat(i, j) = mat_el_2h1p_sI_sII (ov1, ov2, a, ov3, a)

                end do
            end do

        else if ((type_1 == 1) .and. (type_2 == 3)) then
            do i = 1, n_config
                ov1 = config_2h(i, 1)
                ov2 = config_2h(i, 2)

                mat(i, i) = mat_el_2h1p_sI_t (ov1, ov2, a, ov1, ov2, a)

                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)

                    mat(i, j) = mat_el_2h1p_sI_t (ov1, ov2, a, ov3, ov4, a)

                end do
            end do

        else if ((type_1 == 2) .and. (type_2 == 3)) then
            do i = 1, n_ov 
                ov1 = i + h_f_min - 1

                do j = 1, n_config
                    ov2 = config_2h(j, 1)
                    ov3 = config_2h(j, 2)

                    mat(i, j) = mat_el_2h1p_sII_t (ov1, a, ov2, ov3, a)

                end do
            end do
        end if

    end function coupling_2h1p

end module ci_blocks_el
