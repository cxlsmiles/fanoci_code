module setup_variables
	use globals
	use read_hf_data
        use intindex_module
        use reduce_moint
	use configurations
	use misc
	use coupling
	contains
	
	subroutine setup ()
		integer*8 :: temp
        integer :: i
		double precision :: t_in, t_fi

      ! Determine number of orbitals from HForbitals.txt
      ! it can be modified later manually in input.txt
		call read_n_mos()
		call read_input()
		write(*,'(A25)')" Input successfully read."
		write(*,'(A21,I5)')" Molecular orbitals: ",n_mo
	
		allocate(emo(n_mo))
		call read_mo_energies(n_mo, emo)
        
        !-----------------------
        ! Two-electron integrals
        !-----------------------
		if (prop_type == "el")  then
			n_int2e = intindex(n_occ, n_occ, n_mo, n_mo)
		else if (prop_type == "pol") then
			temp = n_mo + 1
			n_int2e = temp * (temp+1) * (temp*temp + temp + 2)/8
		end if

        allocate(int2e(n_int2e))
        int2e(:) = 0.0d0
        
        if (moint .eq. 1) then
            if (prop_type == "el") then
                call reduce_int2e()
                stop                
            else
                write(*,*)"You cannot reduce the number of 2e integrals."
                stop
            end if
        else
	    call cpu_time(t_in)
            call read_2e_integrals(n_int2e, int2e)  
            call cpu_time(t_fi)
        end if
        
        
        write(*,'(A25,I15)')" Number of 2e integrals: ",n_int2e
        write(*,'(A30,f5.2,A2)')" Two-electron integrals read: ",t_fi-t_in," s"
		
        !---------------
        ! Overlap matrix
        !---------------
        allocate(s_mat(n_mo,n_mo))
        call initialize_smat()
        
        
        
        n_ov = h_f_max - h_f_min + 1
        n_config = n_ov * (n_ov - 1)/2
        
        
        allocate(config_2h(n_config,2))
        call generate_2h_configurations()
        
        
        
        
        if (prop_type == "el") then
        	if (gam == "singlet") then
        		ci_mat_size = n_config + n_ov
        	else if (gam == "triplet") then
        		ci_mat_size = n_config
        	else
        		ci_mat_size = n_ov * n_ov
        	end if
        	
        else if (prop_type == "pol") then
                ci_mat_in_size = n_hs_in * n_ps_in
        	ci_mat_size = n_ov * (n_ov + 1)
                allocate(config_1h1p(ci_mat_in_size,2))
                call generate_1h1p_configurations()
        end if
        write(*,'(A30,I10,A3,I10)')" The size of the CI matrix is ", ci_mat_size, " x ", ci_mat_size
        
        
        if (prop_type == "el") then
            allocate(coupling_arr_el(ci_mat_size))
        else if (prop_type == "pol") then
            allocate(coupling_arr_pol(ci_mat_in_size, ci_mat_size))
        end if


        allocate(config_2hnp_str(ci_mat_size))
        allocate(config_2h_str(n_ov + 2*n_config))

        call generate_2h_config_str()
        


        allocate(CI_init_matrix(ci_mat_in_size, ci_mat_in_size))
        allocate(CI_init_eigval(ci_mat_in_size))
        
        allocate(CI_matrix(ci_mat_size,ci_mat_size))
        allocate(CI_eigval(ci_mat_size))
        
        
        if (p_f_max == 0) then
        	p_f_max = n_mo
        end if


	end subroutine setup
	
end module setup_variables
