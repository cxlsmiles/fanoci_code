module build_ci_mat
    use globals
    use ci_blocks_el
    use ci_blocks_pol
    use mat_element_2h
    implicit none
    contains
    
    
    subroutine build_ci_mat_init ()
        integer :: i, j
        integer :: h1, h2, p1, p2
        
        CI_init_matrix(:,:) = 0.0d0
        if (prop_type == "pol") then
            do i = 1, ci_mat_in_size
                h1 = config_1h1p(i,1)
                p1 = config_1h1p(i,2)

                do j = 1, ci_mat_in_size
                    h2 = config_1h1p(j,1)
                    p2 = config_1h1p(j,2)

                    CI_init_matrix(i, j) = mat_el_1h1p (h1, p1, h2, p2)
                end do
            end do
        end if
    end subroutine build_ci_mat_init
    

    subroutine build_cimatrix(a)
        integer, intent(in) :: a
        integer :: offset_2h1p_t
        integer :: offset_2h2p_a, offset_2h2p_b, offset_2h2p_1

        CI_matrix(:,:) = 0.0d0
    
        if (prop_type == "el") then
        !--------------------------------------------------------
        !   Construction of the CI_matrix for a given virtual
        !   orbital k: only the upper triangular part of the
        !   matrix is created
        !
        !   1. main block singlet I
        !   2. main block singlet II
        !   3. main block triplet III
        !--------------------------------------------------------
        offset_2h1p_t = n_config + n_ov
        if (gam == "singlet") then
           CI_matrix(1:n_config, 1:n_config) = main_2h1p(1, a, n_config, n_config)
           CI_matrix((n_config+1):ci_mat_size, (n_config+1):ci_mat_size) = &
               & main_2h1p(2, a, n_ov, n_ov)
    
        else if (gam == "triplet") then
           CI_matrix = main_2h1p(3, a, n_config, n_config)
    
        else
           CI_matrix(1:n_config, 1:n_config) = main_2h1p(1, a, n_config, n_config)
           CI_matrix((n_config+1):offset_2h1p_t, (n_config+1):offset_2h1p_t) =&
               & main_2h1p(2, a, n_ov, n_ov)
           CI_matrix((offset_2h1p_t+1):ci_mat_size, (offset_2h1p_t+1):ci_mat_size) =&
               & main_2h1p(3, a, n_config, n_config)
        end if
        !--------------------------------------------------------
        !   4. coupling block S I - S II
        !   5. coupling block S I - T III
        !   6. coupling block S II - T III
        !--------------------------------------------------------
        if (gam == "singlet") then
           CI_matrix(1:n_config, (n_config+1):offset_2h1p_t) = &
               & coupling_2h1p(a, 1, 2, n_config, n_ov)
    
        else if (gam == "total") then
           CI_matrix(1:n_config, (n_config+1):offset_2h1p_t) = &
               & coupling_2h1p(a, 1, 2, n_config, n_ov)
                
           CI_matrix(1:n_config, (offset_2h1p_t+1):ci_mat_size) = &
               & coupling_2h1p(a, 1, 3, n_config, n_config)
                    
           CI_matrix((n_config+1):offset_2h1p_t, (offset_2h1p_t+1):ci_mat_size) = &
               & coupling_2h1p(a, 2, 3, n_ov, n_config)
        end if
    
    
        else if (prop_type == "pol") then
        !--------------------------------------------------------
        !   Construction of the CI_matrix for a given virtual
        !   orbital k: only the upper triangular part of the
        !   matrix is created
        !
        !   1. main block 1h1p
        !   2. main block 2h2p | ii ab>
        !   3. main block 2h2p | ij ab, A>
        !   4. main block 2h2p | ij ab, B>
        !--------------------------------------------------------
    
        offset_2h2p_1 = n_ov + 1
        offset_2h2p_a = 2*n_ov + 1
        offset_2h2p_b = 2*n_ov + n_config + 1


        CI_matrix(1:n_ov, 1:n_ov) = main_1h1p(a)

        CI_matrix(offset_2h2p_1 : offset_2h2p_a-1, offset_2h2p_1 : offset_2h2p_a-1) =&
              & main_2h2p (p_in, a, 1, n_ov, n_ov)
 
        CI_matrix(offset_2h2p_a : offset_2h2p_b-1, offset_2h2p_a : offset_2h2p_b-1) =&
              & main_2h2p (p_in, a, 2, n_config, n_config)
            
        CI_matrix(offset_2h2p_b : ci_mat_size, offset_2h2p_b : ci_mat_size) =&
              & main_2h2p (p_in, a, 3, n_config, n_config)
    
        !--------------------------------------------------------
        !   5. coupling block 1h1p 2h2p |ii ab>
        !   6. coupling block 1h1p 2h2p |ij ab, A>
        !   7. coupling block 1h1p 2h2p |ij ab, B>
        !   8. coupling block 2h2p |ii ab> / |ij ab, A>
        !   9. coupling block 2h2p |ii ab> / |ij ab, B>
        !   10. coupling block 2h2p |ij ab, A> / |ij ab, B>
        !--------------------------------------------------------
        CI_matrix(1:n_ov, offset_2h2p_1 : offset_2h2p_a-1) = &
              & coupling_1h1p_2h2p(p_in, a, 1, n_ov, n_ov)

        CI_matrix(1:n_ov, offset_2h2p_a : offset_2h2p_b-1) = &
              & coupling_1h1p_2h2p(p_in, a, 2, n_ov, n_config)

        CI_matrix(1:n_ov, offset_2h2p_b : ci_mat_size) = &
              & coupling_1h1p_2h2p(p_in, a, 3, n_ov, n_config)
        
        CI_matrix(offset_2h2p_1 : offset_2h2p_a-1, offset_2h2p_a : offset_2h2p_b-1) = &
              & coupling_2h2p_2h2p (p_in, a, 1, 2, n_config)
        
        CI_matrix(offset_2h2p_1 : offset_2h2p_a-1, offset_2h2p_b : ci_mat_size) = &
              & coupling_2h2p_2h2p (p_in, a, 1, 3, n_config)
        
        CI_matrix(offset_2h2p_a : offset_2h2p_b-1, offset_2h2p_b : ci_mat_size) = &
              & coupling_2h2p_2h2p (p_in, a, 2, 3, n_config)

        end if
    end subroutine build_cimatrix


    subroutine build_ci2hmatrix (type, mat_2h, m_size)
        integer, intent(in) :: type
        integer, intent(in) :: m_size
        integer :: i, j, h1, h2, h3, h4
        double precision, dimension (m_size, m_size) :: mat_2h

        mat_2h(:,:) = 0.d0
        if (type == 0) then  ! singlet

            do i = 1, n_config
               h1 = config_2h(i, 1)
               h2 = config_2h(i, 2)

               do j = 1, n_config
                  h3 = config_2h(j, 1)
                  h4 = config_2h(j, 2)

                  mat_2h(i, j) = mat_el_2h_s (1, h1, h2, h3, h4)

               end do
            end do


            do i = 1, n_ov
               h1 = h_f_min - 1 + i

               do j = 1, n_ov
                  h2 = h_f_min - 1 + j
                  mat_2h(n_config + i, n_config + j) = mat_el_2h_s (2, h1, h1, h2, h2)
               end do
            end do


            do i = 1, n_config
               h1 = config_2h(i, 1)
               h2 = config_2h(i, 2)

               do j = 1, n_ov
                  h3 = h_f_min - 1 + j

                  mat_2h(i, n_config + j) = mat_el_2h_s (3, h3, h3, h1, h2)
                  mat_2h(n_config + j, i) = mat_el_2h_s (3, h3, h3, h1, h2)
               end do
            end do

        else if (type == 1) then ! triplet

            do i = 1, n_config
               h1 = config_2h(i, 1)
               h2 = config_2h(i, 2)
               do j = 1, n_config
                  h3 = config_2h(j, 1)
                  h4 = config_2h(j, 2)
                  mat_2h(i, j) = mat_el_2h_t (h1, h2, h3, h4)
               end do
            end do

        end if
    end subroutine build_ci2hmatrix


end module build_ci_mat
