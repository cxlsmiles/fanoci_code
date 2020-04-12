program fanoci_main
     use globals
     use setup_variables
     use build_ci_mat
     use lapack
     use clear_resources
     implicit none
     
     
     
     integer :: i, a, j, count, m, temp_int, ov1, ov2
     double precision :: eps_h_in, en_e
     character(len=100) :: file_gamma
     double precision, dimension(:), allocatable :: eigval_temp
     double precision, dimension(:), allocatable :: vec_2h
     double precision, dimension(:), allocatable :: vec_ci_coupling
     double precision, dimension(:), allocatable :: temp_vec
     double precision, dimension(:,:), allocatable :: mat_prod

     !
     !   Time variables
     !
     double precision :: prog_start, prog_end, t_2eint_in, t_2eint_fi
     double precision :: it_begin, it_end, it_time


     call setup ()
     
     !----------------------------------
     ! 2h CI matrix
     !----------------------------------
     if(gam.eq.'singlet'.or.gam.eq.'total')then
         ci_mat_2h_s_size = n_config + n_ov
         allocate(mat_2h_s(ci_mat_2h_s_size,ci_mat_2h_s_size))
         allocate(mat_2h_s_eigval(ci_mat_2h_s_size)) 
     
         call build_ci2hmatrix(0, mat_2h_s, ci_mat_2h_s_size)
         call RealSymm(mat_2h_s, ci_mat_2h_s_size, mat_2h_s_eigval)
     end if

     if(gam.eq.'triplet'.or.gam.eq.'total')then
         ci_mat_2h_t_size = n_config
         allocate(mat_2h_t(ci_mat_2h_t_size,ci_mat_2h_t_size))
         allocate(mat_2h_t_eigval(n_config))

         call build_ci2hmatrix(1, mat_2h_t, ci_mat_2h_t_size)
         call RealSymm(mat_2h_t, ci_mat_2h_t_size, mat_2h_t_eigval)
     end if

     allocate(mat_2h(ci_mat_size,ci_mat_size))
     allocate(mat_2h_eigval(ci_mat_size))
     mat_2h(:,:) = 0.0d0
     mat_2h_eigval(:) = 0.0d0

     if(gam.eq.'total')then
         mat_2h(1:ci_mat_2h_s_size, 1:ci_mat_2h_s_size) = mat_2h_s
         mat_2h((ci_mat_2h_s_size+1):, (ci_mat_2h_s_size+1):) = mat_2h_t

         mat_2h_eigval(1:ci_mat_2h_s_size) = mat_2h_s_eigval
         mat_2h_eigval((ci_mat_2h_s_size+1):) = mat_2h_t_eigval
     else if(gam.eq.'singlet')then
         mat_2h = mat_2h_s
         mat_2h_eigval = mat_2h_s_eigval
     else if(gam.eq.'triplet')then
         mat_2h = mat_2h_t
         mat_2h_eigval = mat_2h_t_eigval
     end if

     call print_2h_matrix()


    !---------------------------------
    !---------------------------------
    ! End of 2h CI matrix diagonalisation
    !---------------------------------

    if (prop_type == "el") then
        open(unit=718, file="plot_all.dat")
        !open(unit=228, file="ww_coupling_matel.dat")
        !write(228,*) ci_mat_size*(p_f_max - p_f_min + 1)

        h_in = holes_arr(1)
        eps_h_in = emo(h_in)
        
        write(file_gamma,'(A10,A3,A1,I1,A4)') "Couplings.",gam,".",h_in,".txt"
        open(251, file=file_gamma)
        write(251,*) ci_mat_size*(p_f_max - p_f_min + 1)

        coupling_arr_el(:) = 0.0d0
        allocate(mat_prod(ci_mat_size, ci_mat_size))


        do i = 1, ci_mat_size
            open(1000+i)
            write(1000+i,*)ci_mat_size*(p_f_max - p_f_min + 1)
        end do

        allocate(vec_ci_coupling(ci_mat_size))
        allocate(temp_vec(ci_mat_size))

        
        do a = p_f_min, p_f_max
            call cpu_time(it_begin)
            
            call build_cimatrix(a)

            call RealSymm(CI_matrix, ci_mat_size, CI_eigval)

            !------------------------------------------------------------------
            !       Coupling array
            !------------------------------------------------------------------

            call coupling_elements(a)


            !-------------------------------------------------------
            ! Write the WW coupling matrix elements

!            do i = 1, n_config
!               ov1 = config_2h(i,1)
!               ov2 = config_2h(i,2)
!               en_e = mat_el_2h1p_sI(ov1, ov2, a, ov1, ov2, a) - abs(eps_h_in)
!               write(228, *)en_e,(coupling_arr_el(i))**2
!            end do
!           
!            do i = 1, n_ov
!               ov1 = i + h_f_min - 1
!               en_e = mat_el_2h1p_sII(ov1, a, ov1, a) - abs(eps_h_in)
!               write(228, *)en_e,(coupling_arr_el(i + n_config))**2
!            end do
! 
!            do i = n_config+n_ov+1, ci_mat_size
!               temp_int = i - n_config-n_ov
!               ov1 = config_2h(temp_int,1)
!               ov2 = config_2h(temp_int,2)
!               en_e = mat_el_2h1p_t (ov1, ov2, a, ov1, ov2, a) - abs(eps_h_in)
!               write(228, *)en_e,(coupling_arr_el(i))**2
!            end do

            ! Finish WW coupling matrix elements
            !----------------------------------------------------------

            !----------------------------------------------------------
            !  Overlap matrix between the 2h CI states and the 2h1p
            !  CI states
            !----------------------------------------------------------
            mat_prod(:,:) = 0.0d0
            mat_prod = matmul(transpose(mat_2h), CI_matrix) ! rows - 2h states, cols - CI 2h1p states

           
            !----------------------------------------------------------
            !  Coupling matrix elements |<chi_a_i|H|psi_d>|2
            !---------------------------------------------------------- 
            vec_ci_coupling(:) = 0.d0
            do i = 1, ci_mat_size
               vec_ci_coupling(i) = dot_product(coupling_arr_el(:),CI_matrix(:,i))
               write(251, *)CI_eigval(i)-abs(eps_h_in), vec_ci_coupling(i)**2
            end do


            do i = 1, ci_mat_size ! loop over the 2h ci states
               do j = 1, ci_mat_size ! loop over the 2h1p ci states
                   write(1000+i,*)CI_eigval(j)-abs(eps_h_in), (mat_prod(i,j) * vec_ci_coupling(j))**2
               end do
            end do

            
            if ( prnt .eq. 1) then
                call print_ci_eigval(a)
                call print_2h_contributions(mat_prod,a)
            end if
            
            call cpu_time(it_end)
            it_time = it_end - it_begin
            write(*,'(A5,I3,A6,F10.5,A3)')"Step ", a - p_f_min + 1, " took ", it_time, " s."
        end do

        deallocate(mat_prod)    
        deallocate(vec_ci_coupling)
        deallocate(temp_vec)


        close(718)
        close(251)


        do i = 1, ci_mat_size
            close(1000+i)
        end do
        
    else if (prop_type == "pol") then

        write(file_gamma,'(A10,A3,A4)') "Couplings.", gam, ".txt"
        open(251, file=file_gamma)
        write(251,*) ci_mat_size*(p_f_max - p_f_min + 1)*ci_mat_in_size

        open(unit=718, file="plot_all.dat")

        call build_ci_mat_init ()

        call print_mat(ci_mat_in_size, ci_mat_in_size,CI_init_matrix)

        call RealSymm(CI_init_matrix, ci_mat_in_size, CI_init_eigval)
       
        coupling_arr_pol(:,:) = 0.0d0
        h_in = config_1h1p(1, 1)
        p_in = config_1h1p(1, 2)
                    
        count = 1
        do a = p_f_min, p_f_max
        
            call cpu_time(it_begin)
            
            call build_cimatrix(a)

            call RealSymm(CI_matrix, ci_mat_size, CI_eigval)

            call coupling_elements(a)
            
            if ( prnt .eq. 1) then
                call print_ci_eigval(a)
            end if

            !call print_2h_contributions(a)
        
        !------------------------------------------------------------------
        !       Coupling array
        !------------------------------------------------------------------

            do j = 1, ci_mat_size
               !write(251, *)CI_eigval(j)-abs(CI_init_eigval(1)), (dot_product(matmul(coupling_arr_pol(1,:),CI_matrix(:,j)),CI_init_matrix(:,1)))**2
               write(251, *)CI_eigval(j)-abs(CI_init_eigval(1)), (dot_product(coupling_arr_pol(1,:),CI_matrix(:,j)))**2
            end do
            
            call cpu_time(it_end)
            it_time = it_end - it_begin
            write(*,'(A5,I3,A6,F10.5,A3)')"Step ", count, " took ", it_time, " s."
            count = count + 1
       end do    
	
    close(718)
    close(251)
    end if
	call cpu_time(prog_end)

    write(*,'(A15,F8.2,A15)')"Diagonalization time ", (prog_end - t_2eint_fi)/(p_f_max - p_f_min + 1), " s / iteration "
    write(*,'(A15,F8.2,A2)')"Execution time ", prog_end - prog_start, " s"

    write(*,*)"THE END! :)"
	call clear_res()
	
end program fanoci_main
