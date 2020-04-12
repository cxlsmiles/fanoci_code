module clear_resources
    use globals
    contains
    
    subroutine clear_res ()
        deallocate(emo)
        deallocate(int2e)
        deallocate(s_mat)
        deallocate(config_2h)

        if (prop_type == "el") then
           deallocate(coupling_arr_el)
        else 
           deallocate(config_1h1p)
           deallocate(coupling_arr_pol)
        end if

        deallocate(config_2h_str)
        deallocate(CI_init_matrix)
        deallocate(CI_matrix)
        deallocate(CI_init_eigval)
        deallocate(CI_eigval)
        if(gam.eq.'singlet'.or.gam.eq.'total')then
            deallocate(mat_2h_s)
            deallocate(mat_2h_s_eigval)
        end if
        if(gam.eq.'triplet'.or.gam.eq.'total')then
            deallocate(mat_2h_t)
            deallocate(mat_2h_t_eigval)
        end if
        deallocate(mat_2h)
    end subroutine clear_res

end module clear_resources
