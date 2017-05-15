module coupling
    use globals
    use mat_element_el
    use mat_element_pol
    contains
    
    subroutine coupling_elements(k)
        integer :: k, i, counter
        
        if (prop_type == "el") then
            coupling_arr_el(:) = 0.0d0
            if (gam == "singlet") then
                counter = 1
                do i = 1, n_config
                   coupling_arr_el(counter) = mat_el_1h_2h1p_sI(h_in, config_2h(i,1), config_2h(i, 2), k)
                   write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A3)')"|",config_2h(i, 1), ",",config_2h(i, 2),",", k," s>"
                   counter = counter + 1
                end do
        
                do i = 1, n_ov
                   coupling_arr_el(counter) = mat_el_1h_2h1p_sII(h_in, i+h_f_min-1, k)
                   write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A3)')"|",i+h_f_min-1, ",",i+h_f_min-1,",", k," s>"
                   counter = counter + 1
                end do
        
        else if (gam == "triplet") then
        
            do i = 1, n_config
                coupling_arr_el(i) = mat_el_1h_2h1p_t(h_in, config_2h(i,1), config_2h(i, 2), k)
                write(config_2hnp_str(i),'(A1,I3,A1,I3,A1,I3,A3)')"|",config_2h(i, 1), ",",config_2h(i, 2),",", k," t>"
            end do
        
        else if (gam == "total") then
            counter = 1
            do i = 1, n_config
               coupling_arr_el(counter) = mat_el_1h_2h1p_sI(h_in, config_2h(i,1), config_2h(i, 2), k)
               write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A3)')"|",config_2h(i, 1), ",",config_2h(i, 2),",", k," s>"
               counter = counter + 1
            end do
        
            do i = 1, n_ov
               coupling_arr_el(counter) = mat_el_1h_2h1p_sII(h_in, i + h_f_min-1, k)
               write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A3)')"|",i+h_f_min-1, ",",i+h_f_min-1,",", k," s>"
          counter = counter + 1
       end do
        
       do i = 1, n_config
          coupling_arr_el(counter) = mat_el_1h_2h1p_t(h_in, config_2h(i,1), config_2h(i, 2), k)
           write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A3)')"|",config_2h(i, 1), ",",config_2h(i, 2),",", k," t>"
           counter = counter + 1
       end do
       end if
    
    else if (prop_type == "pol") then
    
        h_in = config_1h1p(1, 1)
        p_in = config_1h1p(1, 2)
        counter = 1
        do i = 1, n_ov
            coupling_arr_pol(1, counter) = mat_el_1h1p(h_in, p_in, h_f_min + i - 1, k)
            write(config_2hnp_str(counter),'(A1,I3,A1,I3,A2)')"|",i+h_f_min-1,",", k," >"
            counter = counter + 1
        end do
    
        do i = 1, n_ov
            coupling_arr_pol(1, counter) = mat_el_1h1p_2h2p_1(h_in, p_in, h_f_min + i - 1, p_in, k)
            write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A1,I3,A2)')"|",i+h_f_min-1, ",",i+h_f_min-1,",", k,",", p_in," >"
            counter = counter + 1
        end do
    
        do i = 1, n_config
            coupling_arr_pol(1, counter) = mat_el_1h1p_2h2p_A(h_in, p_in, config_2h(i, 1), config_2h(i, 2), p_in, k)
            write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A1,I3,A3)')"|",config_2h(i, 1), ",",config_2h(i, 2),",", k,",", p_in," A>"
            counter = counter + 1
        end do

        do i = 1, n_config
            coupling_arr_pol(1, counter) = mat_el_1h1p_2h2p_B(h_in, p_in, config_2h(i, 1), config_2h(i, 2), p_in, k)
            write(config_2hnp_str(counter),'(A1,I3,A1,I3,A1,I3,A1,I3,A3)')"|",config_2h(i, 1), ",",config_2h(i, 2),",", k,",", p_in," B>"
            counter = counter + 1
        end do
        end if

    end subroutine coupling_elements

end module coupling
