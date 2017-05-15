module mat_element_2h
    use globals 
    use intindex_module 
    contains

    function mat_el_2h_s (mult, i, j, k, l) result (mat_el)
        integer :: i, j, k, l, mult
        double precision :: mat_el
        double precision :: dik, dil, djk, djl

        dik = s_mat(i, k)
        dil = s_mat(i, l)
        djk = s_mat(j, k)
        djl = s_mat(j, l)

        mat_el = 0.d0
        if (mult == 1) then ! Singlet I - singlet I
           mat_el = (dil * djk + dik * djl) * (- emo(i) - emo(j)) &
                   &+ int2e(intindex(k, i, l, j)) + int2e(intindex(k, j, l, i))

        else if (mult == 2) then ! Singlet II - singlet II
           mat_el = - 2.d0 * dik * emo(i) + int2e(intindex(k, i, k, i))


        else if (mult == 3) then ! Singlet I - singlet II
           mat_el = -sqrt(2.d0) * int2e(intindex(k, i, l, i))
!           print*,k, i, l, i
!           print*,intindex(k, i, l, i), int2e(intindex(k, i, l, i))


        end if

    end function mat_el_2h_s
    

    function mat_el_2h_t (i, j, k, l) result (mat_el)
        integer :: i, j, k, l
        double precision :: mat_el
        double precision :: dik, dil, djk, djl

        dik = s_mat(i, k)
        dil = s_mat(i, l)
        djk = s_mat(j, k)
        djl = s_mat(j, l)

        ! Triplet - triplet 
        mat_el = (dil * djk - dik * djl) * (emo(i) + emo(j)) &
                &+ int2e(intindex(k, i, l, j)) - int2e(intindex(k, j, l, i))

    end function mat_el_2h_t


end module mat_element_2h
