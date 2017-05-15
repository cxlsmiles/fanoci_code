module mat_element_el
    use globals
	use intindex_module
    contains	
	
    ! =================
    !    MAIN BLOCKS
    ! =================	

    ! The 2-electron integrals are in the chemists notation
	
    function mat_el_2h1p_sI (i, j, a, k, l, c) result (mat_el)
	! -----------------
	! computes < ^{ds1}\Psi_{ija} | H - E_0 | ^{ds1}\Psi_{klc} >
	! -----------------
        integer :: i, j, a, k, l, c
        double precision :: mat_el
        double precision :: dac, dik, dil, djk, djl

        dac = s_mat(a, c)
        dik = s_mat(i, k)
        dil = s_mat(i, l)
        djk = s_mat(j, k)
        djl = s_mat(j, l)

        mat_el = (dik * djl * dac + dil * djk * dac) * (emo(c) - emo(k) - emo(l)) &
				+ dac * int2e(intindex(k, j, l, i)) + dac * int2e(intindex(k, i, l, j)) &
				+ 0.5d0 * djl * int2e(intindex(a, i, k, c)) - djl * int2e(intindex(a, c, k, i)) &
				+ 0.5d0 * dik * int2e(intindex(a, j, l, c)) - dik * int2e(intindex(a, c, l, j)) &
				+ 0.5d0 * djk * int2e(intindex(a, i, l, c)) - djk * int2e(intindex(a, c, l, i)) &
				+ 0.5d0 * dil * int2e(intindex(a, j, k, c)) - dil * int2e(intindex(a, c, k, j))
		
    end function mat_el_2h1p_sI


    function mat_el_2h1p_sII (i, a, k, c) result (mat_el)
	! -----------------
    ! computes < ^{ds2}\Psi_{iia} | H - E_0 | ^{ds2}\Psi_{kkc} >
	! -----------------
        integer :: i, a, k, c
        double precision :: dik, dac
        double precision :: mat_el

        dac = s_mat(a, c)
        dik = s_mat(i, k)

        mat_el = dik * dac * (emo(c) - emo(k) - emo(k)) + dac * int2e(intindex(k, i, k, i))&
				+ dik * int2e(intindex(a, i, k, c)) -2.0d0 * dik * int2e(intindex(a, c, k, i))
				
    end function mat_el_2h1p_sII


    function mat_el_2h1p_sI_sII (i, j, a, k, c) result (mat_el)
	! -----------------
    ! computes < ^{ds1}\Psi_{ija} | H - E_0 | ^{ds2}\Psi_{kkc} >
	! -----------------
        integer :: i, j, a, k, c
        double precision :: dik, djk, dac
        double precision :: mat_el

        dac = s_mat(a, c)
        dik = s_mat(i, k)
        djk = s_mat(j, k)

        mat_el = -sqrt(0.5d0) *(&
				-2.0d0 * djk * int2e(intindex(a, c, k, i)) + djk * int2e(intindex(a, i, k, c))&
				-2.0d0 * dik * int2e(intindex(a, c, k, j)) + dik * int2e(intindex(a, j, k, c))&
				+ dac * int2e(intindex(k, i, k, j)) + dac * int2e(intindex(k, j, k, i)))

    end function mat_el_2h1p_sI_sII


    function mat_el_2h1p_sI_t (i, j, a, k, l, c) result (mat_el)
	! -----------------
    ! computes < ^{ds1}\Psi_{ija} | H - E_0 | ^{dt}\Psi_{klc} >
	! -----------------
        integer :: i, j, a, k, l, c
        double precision :: dik, dil, djk, djl
        double precision :: mat_el

        dik = s_mat(i, k)
        dil = s_mat(i, l)
        djk = s_mat(j, k)
        djl = s_mat(j, l)

        mat_el = 0.5d0 * sqrt(3.0d0) * &
				(dik * int2e(intindex(a, j, l, c)) - djl * int2e(intindex(a, i, k, c))&
				- dil * int2e(intindex(a, j, k, c)) + djk * int2e(intindex(a, i, l, c)))
				
    end function mat_el_2h1p_sI_t


    function mat_el_2h1p_sII_t (k, c, i, j, a) result (mat_el)
	! -----------------
    ! computes < ^{ds2}\Psi_{kkc} | H - E_0 | ^{dt}\Psi_{ija} >
	! -----------------
        integer :: i, j, a, k, c
        double precision :: dik, djk, dac
        double precision :: mat_el

        dac = s_mat(a, c)
        dik = s_mat(i, k)
        djk = s_mat(j, k)

        mat_el = -sqrt(1.5d0) * (dik * int2e(intindex(c, k, j, a)) - djk * int2e(intindex(c, k, i, a)))
		
    end function mat_el_2h1p_sII_t


    function mat_el_2h1p_t (i, j, a, k, l, c) result (mat_el)
	! -----------------
    ! computes < ^{dt}\Psi_{ija} | H - E_0 | ^{dt}\Psi_{klc}
	! -----------------
        integer :: i, j, a, k, l, c
        double precision :: dik, dil, djk, djl, dac
        double precision :: mat_el

        dik = s_mat(i, k)
        dil = s_mat(i, l)
        djk = s_mat(j, k)
        djl = s_mat(j, l)
        dac = s_mat(a, c)

        mat_el = (dik * djl * dac - dil * djk * dac) * (emo(c) - emo(k) - emo(l))&
				- dac * int2e(intindex(k, j, l, i)) + dac * int2e(intindex(k, i, l, j))&
				+ 1.5d0 * djl * int2e(intindex(a, i, k, c)) - djl * int2e(intindex(a, c, k, i))&
				+ 1.5d0 * dik * int2e(intindex(a, j, l, c)) - dik * int2e(intindex(a, c, l, j))&
				+ djk * int2e(intindex(a, c, l, i)) - 1.5d0 * djk * int2e(intindex(a, i, l, c))& 
				+ dil * int2e(intindex(a, c, k, j)) - 1.5d0 * dil * int2e(intindex(a, j, k, c))

    end function mat_el_2h1p_t


    function mat_el_1h_2h1p_sI (i, k, l, c) result (mat_el)
	! -----------------
    ! computes < ^{d}\Psi_{i} | H | ^{ds1}\Psi_{klc} >
	! -----------------
        integer :: i, k, l, c
        double precision :: mat_el

        mat_el = sqrt(0.5d0) * (int2e(intindex(k, i, l, c)) + int2e(intindex(k, c, l, i)))

    end function mat_el_1h_2h1p_sI


    function mat_el_1h_2h1p_sII (i, k, c) result (mat_el)
	! -----------------
	! computes < ^{d}\Psi_{i} | H | ^{ds2}\Psi_{kkc} >
	! -----------------
        integer :: i, k, c
        double precision :: mat_el

        mat_el = -int2e(intindex(k, i, k, c))

    end function mat_el_1h_2h1p_sII
	

    function mat_el_1h_2h1p_t (i, k, l, c) result (mat_el)
	! -----------------
    ! computes < ^{d}\Psi_{i} | H | ^{dt}\Psi_{klc} >
	! -----------------
        integer :: i, k, l, c
        double precision :: mat_el

        mat_el = sqrt(1.5d0) * (int2e(intindex(k, i, l, c)) - int2e(intindex(k, c, l, i)))

    end function mat_el_1h_2h1p_t

end module mat_element_el
