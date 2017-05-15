module build_ci_mat
	use globals
	use ci_blocks_el
	use ci_blocks_pol
    implicit none
    contains

    subroutine build_cimatrix(a)
		integer*8 :: a
		integer*8 :: offset_2h1p_t
		integer*8 :: offset_2h2p_a, offset_2h2p_b, offset_2h2p_1
		
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
    ! main blocks

            offset_2h2p_1 = n_ov + 1
            offset_2h2p_a = 2 * n_ov + 1
            offset_2h2p_b = 2 * n_ov + n_config + 1

!            CI_matrix(1:n_ov, 1:n_ov) = main_1h1p(a)
            
			CI_matrix(n_ov+1 : offset_2h2p_a-1, n_ov+1 : offset_2h2p_a-1) =&
						& main_2h2p (p_in, a, 1, n_ov, n_ov)
            
			CI_matrix(offset_2h2p_a : offset_2h2p_b-1, offset_2h2p_a : offset_2h2p_b-1) =&
						& main_2h2p (p_in, a, 2, n_config, n_config)
            
			CI_matrix(offset_2h2p_b : ci_mat_size, offset_2h2p_b : ci_mat_size) =&
						& main_2h2p (p_in, a, 3, n_config, n_config)

    ! coupling blocks
!            CI_matrix(1:n_ov, offset_2h2p_1 : offset_2h2p_a-1) = &
!						& coupling_1h1p_2h2p(p_in, a, 1, n_ov, n_ov)
!						
!            CI_matrix(1:n_ov, offset_2h2p_a : offset_2h2p_b-1) = &
!						& coupling_1h1p_2h2p(p_in, a, 2, n_ov, n_config)
!						
!            CI_matrix(1:n_ov, offset_2h2p_b : ci_mat_size) = &
!						& coupling_1h1p_2h2p(p_in, a, 3, n_ov, n_config)
            
			CI_matrix(n_ov+1 : offset_2h2p_a-1, offset_2h2p_a : offset_2h2p_b-1) = &
						& coupling_2h2p_2h2p (p_in, a, 1, 2, n_config)
            
			CI_matrix(n_ov+1 : offset_2h2p_a-1, offset_2h2p_b : ci_mat_size) = &
						& coupling_2h2p_2h2p (p_in, a, 1, 3, n_config)
            
			CI_matrix(offset_2h2p_a : offset_2h2p_b-1, offset_2h2p_b : ci_mat_size) = &
						& coupling_2h2p_2h2p (p_in, a, 2, 3, n_config)


        end if

    end subroutine build_cimatrix

end module build_ci_mat
