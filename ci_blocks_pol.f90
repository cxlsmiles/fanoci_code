module ci_blocks_pol
    use globals
    use misc
    use mat_element_pol
    implicit none
    contains
    ! the following functions needed to construct the CI matrix
    !
    !     main blocks:
    ! function main_1h1p (a) result(mat)
    ! function main_2h2p (a, b, csf_type, nrow, ncol) result (mat)
    !
    !     coupling blocks:
    ! function coupling_1h1p_2h2p (a, b, csf_type, nrow, ncol) result (mat)
    ! function coupling_2h2p_2h2p (a, b, csf_type_1, csf_type_2, nrow) result (mat)
    !

    ! =================
    !    MAIN BLOCKS
    ! =================

    ! -----------------
    !   1h1p
    ! -----------------

    function main_1h1p (a) result(mat)
        integer :: i, j, a
        integer :: ov1, ov2
        double precision, dimension(n_ov, n_ov) :: mat

        do i = 1, n_ov
            ov1 = i + h_f_min-1

            do j = 1, n_ov
                ov2 = j + h_f_min-1

                mat(i, j) = mat_el_1h1p (ov1, a, ov2, a)
            end do
        end do

    end function main_1h1p


    ! -----------------
    !   2h2p
    ! -----------------

    function main_2h2p (a, b, csf_type, nrow, ncol) result (mat)
        integer :: i, j, k, l
        integer :: a, b
        integer :: ov1, ov2, ov3, ov4, nrow, ncol
        integer :: csf_type
        double precision, dimension(nrow, ncol) :: mat


        if (csf_type == 1) then
            do i = 1, n_ov
                ov1 = i + h_f_min-1
!                print*,mat_el_test(ov1, a, b) 
!		print*,emo(ov1),emo(a),emo(b)
                do j = 1, n_ov
                    ov2 = j + h_f_min-1
                    mat(i, j) = mat_el_2h2p_11 (ov1, a, b, ov2, a, b)
                end do
            end do

        else if (csf_type == 2) then
            do i = 1, n_config
                ov1 = config_2h (i, 1)
                ov2 = config_2h (i, 2)
!                print*,ov1,ov2,a,b,mat_el_2h2p_a_test(ov1,ov2,a,b)
                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)
                    mat(i, j) = mat_el_2h2p_AA (ov1, ov2, a, b, ov3, ov4, a)
                end do
            end do

        else if (csf_type == 3) then

            do i = 1, n_config
                ov1 = config_2h (i, 1)
                ov2 = config_2h (i, 2)

                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)
                    mat(i, j) = mat_el_2h2p_BB (ov1, ov2, a, b, ov3, ov4, a)
                end do
            end do
        end if

    end function main_2h2p


    ! =================
    !    COUPLING BLOCKS
    ! =================

    ! -----------------
    !   1h1p / 2h2p
    ! -----------------

    function coupling_1h1p_2h2p (a, b, csf_type, nrow, ncol) result (mat)
        integer :: i, j, k
        integer :: a, b
        integer :: ov1, ov2, ov3, ov4, nrow, ncol
        integer :: csf_type
        double precision, dimension(nrow, ncol) :: mat

        if (csf_type == 1) then
            do i = 1, n_ov
                ov1 = i + h_f_min - 1

                do j = 1, n_ov
                    ov2 = j + h_f_min - 1
                    mat(i, j) = mat_el_1h1p_2h2p_1 (ov1, b, ov2, a, b)
                end do
            end do

        else if (csf_type == 2) then
            do i = 1, n_ov
                ov1 = i + h_f_min - 1

                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)
                    mat(i, j) = mat_el_1h1p_2h2p_A (ov1, b, ov3, ov4, a, b)
                end do
            end do

        else if (csf_type == 3) then

            do i = 1, n_ov
                ov1 = i + h_f_min - 1

                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)
                    mat(i, j) = mat_el_1h1p_2h2p_B (ov1, b, ov3, ov4, a, b)
                end do
            end do
        end if

    end function coupling_1h1p_2h2p


    ! -----------------
    !   2h2p / 2h2p
    ! -----------------

    function coupling_2h2p_2h2p (a, b, csf_type_1, csf_type_2, nrow) result (mat)
        integer :: i, j
        integer :: a, b
        integer :: ov1, ov2, ov3, ov4, nrow
        integer :: csf_type_1, csf_type_2
        double precision, dimension (nrow, nrow) :: mat

        if (csf_type_1 == 1) then
            if (csf_type_2 == 2) then
                do i = 1, n_ov
                    ov1 = h_f_min + i - 1

                    do j = 1, n_config
                        ov3 = config_2h(j, 1)
                        ov4 = config_2h(j, 2)

                        mat(i, j) = mat_el_2h2p_1A(ov1, a, b, ov3, ov4, a)
                    end do
                end do

            else if (csf_type_2 == 3) then
                do i = 1, n_ov
                    ov1 = h_f_min + i - 1

                    do j = 1, n_config
                        ov3 = config_2h(j, 1)
                        ov4 = config_2h(j, 2)

                        mat(i, j) = mat_el_2h2p_1B(ov1, a, b, ov3, ov4, a)
                    end do
                end do
            end if

        else if (csf_type_1 == 2) then
            do i = 1, n_config
                ov1 = config_2h(i, 1)
                ov2 = config_2h(i, 2)

                do j = 1, n_config
                    ov3 = config_2h(j, 1)
                    ov4 = config_2h(j, 2)

                    mat(i, j) = mat_el_2h2p_AB(ov1, ov2, a, b, ov3, ov4, a)
                end do
            end do

        end if

    end function coupling_2h2p_2h2p


end module ci_blocks_pol
