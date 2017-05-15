module mat_element_pol
	use globals
	use intindex_module
	contains

	function mat_el_1h1p (i, a, j, b) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{i}^{a} | H | ^{1}\Psi_{j}^{b} >
	! -----------------
		integer :: i, j
		integer :: a, b
		double precision :: mat_el
		double precision :: dij, dab

		dij = s_mat(i, j)
		dab = s_mat(a, b)

		mat_el = dij * dab * (emo(b) - emo(j))&
				+ 2.0d0 * int2e(intindex(a, i, j, b))&
				- int2e(intindex(a, b, j, i))
		
	end function mat_el_1h1p


	! -----------------
	!   2h2p type 1
	! -----------------

	function mat_el_2h2p_11 (i, a, b, k, c, d) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{ii}^{ab} | H | ^{1}\Psi_{kk}^{cb} >
	! -----------------
		integer :: i, k
		integer :: a, b, c, d
		double precision :: mat_el
		double precision ::  dad, dbd, dbc, dac, dik

		dbc = s_mat(b, c)
		dac = s_mat(a, c)
		dbd = s_mat(b, d)
		dad = s_mat(a, d)
		dik = s_mat(i, k)

		mat_el = (dik * dac * dbd + dik * dbc * dad)*( - emo(k) + emo(d) - emo(k) + emo(c)) &
			+ dik * int2e(intindex(c, a, d, b)) + dik * int2e(intindex(c, b, d, a)) &
			+ dbc * dad * int2e(intindex(i, k, i, k)) + dac * dbd * int2e(intindex(i, k, i, k)) &
			-2.0d0 * dac * dik * int2e(intindex(d, b, i, k)) + dac * dik * int2e(intindex(d, k, i, b)) &
			-2.0d0 * dbc * dik * int2e(intindex(d, a, i, k)) + dbc * dik * int2e(intindex(d, k, i, a)) &
			-2.0d0 * dad * dik * int2e(intindex(c, b, i, k)) + dad * dik * int2e(intindex(c, k, i, b)) &
			-2.0d0 * dbd * dik * int2e(intindex(c, a, i, k)) + dbd * dik * int2e(intindex(c, k, i, a))
		
	end function mat_el_2h2p_11


    function mat_el_test(i, a, b) result (mat_el)
        integer :: i, a, b
        double precision :: mat_el

        mat_el = emo(a) + emo(b) - 2*emo(i) + &
                + int2e(intindex(b, i, i, b))&
                - 2.0d0* int2e(intindex(b, b, i, i))&
                + int2e(intindex(a, i, i, a))&
                - 2.0d0* int2e(intindex(a, a, i, i))&
                + int2e(intindex(a, a, b, b))&
                + int2e(intindex(a, b, b, a))&
                + int2e(intindex(i, i, i, i))
                
    end function mat_el_test
    
    function mat_el_2h2p_a_test(i, j, a, b) result (mat_el)
        integer :: i, j,a, b
        double precision :: mat_el

        mat_el = emo(a) + emo(b) - emo(j) - emo(i) + &
                + int2e(intindex(a, a, b, b))&
                - int2e(intindex(a, b, b, a))&
                + int2e(intindex(i, i, j, j))&
                - int2e(intindex(i, j, j, i))&
                + 1.5d0* int2e(intindex(a, i, i, a))&
                - int2e(intindex(a, a, i, i))&
                + 1.5d0* int2e(intindex(b, i, i, b))&
                - int2e(intindex(b, b, i, i))&
                + 1.5d0* int2e(intindex(a, j, j, a))&
                - int2e(intindex(a, a, j, j))&
                + 1.5d0* int2e(intindex(b, j, j, b))&
                - int2e(intindex(b, b, j, j))
    
    end function mat_el_2h2p_a_test

	! -----------------
	!   2h2p type A
	! -----------------

	function mat_el_2h2p_AA (i, j, a, b, k, l, c) result (mat_el)
	! -----------------
	! computes < ^{A}\Psi_{ij}^{ab} | H | ^{A}\Psi_{kl}^{cb} >
	! -----------------
		integer :: i, j, k, l
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dac, dbc, dil, dik, djl, djk

		dab = s_mat(a, b)
		dbc = s_mat(b, c)
		dac = s_mat(a, c)
		djl = s_mat(j, l)
		djk = s_mat(j, k)
		dil = s_mat(i, l)
		dik = s_mat(i, k)

		mat_el = (dik * djl * dac - djk * dil * dac) * (emo(b) + emo(c) - emo(k) - emo(l))&
                + dbc * dab * int2e(intindex(i, l, j, k)) - dbc * dab * int2e(intindex(i, k, j, l))& 
                + dik * djl * int2e(intindex(c, a, b, b)) - dik * djl * int2e(intindex(c, b, b, a))&
                + djk * dil * int2e(intindex(c, b, b, a)) - djk * dil * int2e(intindex(c, a, b, b))&
                - dab * dil * int2e(intindex(c, b, j, k)) + 1.5d0 * dab * dil * int2e(intindex(c, k, j, b))&
                - dab * djk * int2e(intindex(c, b, i, l)) + 1.5d0 * dab * djk * int2e(intindex(c, l, i, b))&
                + dab * djl * int2e(intindex(c, b, i, k)) - 1.5d0 * dab * djl * int2e(intindex(c, k, i, b))& 
                + dab * dik * int2e(intindex(c, b, j, l)) - 1.5d0 * dab * dik * int2e(intindex(c, l, j, b))& 
                + dac * dil * int2e(intindex(b, b, j, k)) - 1.5d0 * dac * dil * int2e(intindex(b, k, j, b))& 
                + dac * djk * int2e(intindex(b, b, i, l)) - 1.5d0 * dac * djk * int2e(intindex(b, l, i, b))&
                - dac * dik * int2e(intindex(b, b, j, l)) + 1.5d0 * dac * dik * int2e(intindex(b, l, j, b))&
                - dac * djl * int2e(intindex(b, b, i, k)) + 1.5d0 * dac * djl * int2e(intindex(b, k, i, b))&
                + dbc * djl * int2e(intindex(b, a, i, k)) - 1.5d0 * dbc * djl * int2e(intindex(b, k, i, a))& 			 
                - dbc * djk * int2e(intindex(b, a, i, l)) + 1.5d0 * dbc * djk * int2e(intindex(b, l, i, a))& 
                - dbc * dil * int2e(intindex(b, a, j, k)) + 1.5d0 * dbc * dil * int2e(intindex(b, k, j, a))& 
                + dbc * dik * int2e(intindex(b, a, j, l)) - 1.5d0 * dbc * dik * int2e(intindex(b, l, j, a))& 
                - dik * int2e(intindex(c, a, j, l)) + 1.5d0 * dik * int2e(intindex(c, l, j, a))& 
                + dil * int2e(intindex(c, a, j, k)) - 1.5d0 * dil * int2e(intindex(c, k, j, a))& 
                + djk * int2e(intindex(c, a, i, l)) - 1.5d0 * djk * int2e(intindex(c, l, i, a))&
                - djl * int2e(intindex(c, a, i, k)) + 1.5d0 * djl * int2e(intindex(c, k, i, a))&  
                - dac * int2e(intindex(i, l, j, k)) + dac * int2e(intindex(i, k, j, l))

	end function mat_el_2h2p_AA


	! -----------------
	!   2h2p type B
	! -----------------

	function mat_el_2h2p_BB (i, j, a, b, k, l, c) result (mat_el)
	! -----------------
	! computes < ^{B}\Psi_{ij}^{ab} | H | ^{B}\Psi_{kl}^{cb} >
	! -----------------
		integer :: i, j, k, l
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dac, dbc, dil, dik, djl, djk

		dab = s_mat(a, b)
		dbc = s_mat(b, c)
		dac = s_mat(a, c)
		djl = s_mat(j, l)
		djk = s_mat(j, k)
		dil = s_mat(i, l)
		dik = s_mat(i, k)

		mat_el = (djk * dil * dac + dik * djl * dac) * (emo(b) + emo(c) - emo(k) - emo(l))&
				+ dac * int2e(intindex(i, l, j, k)) + dac * int2e(intindex(i, k, j, l))&
				+ dbc * dab * int2e(intindex(i, l, j, k)) + dbc * dab * int2e(intindex(i, k, j, l))&
				+ djk * dil * int2e(intindex(c, a, b, b)) + dik * djl * int2e(intindex(c, a, b, b))& 
				+ dik * djl * int2e(intindex(c, b, b, a)) + djk * dil * int2e(intindex(c, b, b, a))&
				+ 0.5d0 * dab * djk * int2e(intindex(c, l, i, b)) - dab * djk * int2e(intindex(c, b, i, l))& 
				+ 0.5d0 * dab * djl * int2e(intindex(c, k, i, b)) - dab * djl * int2e(intindex(c, b, i, k))&
				+ 0.5d0 * dab * dik * int2e(intindex(c, l, j, b)) - dab * dik * int2e(intindex(c, b, j, l))& 
				+ 0.5d0 * dab * dil * int2e(intindex(c, k, j, b)) - dab * dil * int2e(intindex(c, b, j, k))& 
				+ 0.5d0 * dac * dik * int2e(intindex(b, l, j, b)) - dac * dik * int2e(intindex(b, b, j, l))&
				+ 0.5d0 * dac * dil * int2e(intindex(b, k, j, b)) - dac * dil * int2e(intindex(b, b, j, k))&			
				+ 0.5d0 * dac * djk * int2e(intindex(b, l, i, b)) - dac * djk * int2e(intindex(b, b, i, l))& 
				+ 0.5d0 * dac * djl * int2e(intindex(b, k, i, b)) - dac * djl * int2e(intindex(b, b, i, k))& 
				+ 0.5d0 * dbc * djk * int2e(intindex(b, l, i, a)) - dbc * djk * int2e(intindex(b, a, i, l))& 
				+ 0.5d0 * dbc * dil * int2e(intindex(b, k, j, a)) - dbc * dil * int2e(intindex(b, a, j, k))&
				+ 0.5d0 * dbc * djl * int2e(intindex(b, k, i, a)) - dbc * djl * int2e(intindex(b, a, i, k))& 
				+ 0.5d0 * dbc * dik * int2e(intindex(b, l, j, a)) - dbc * dik * int2e(intindex(b, a, j, l))&	
				+ 0.5d0 * dik * int2e(intindex(c, l, j, a)) - dik * int2e(intindex(c, a, j, l))&
				+ 0.5d0 * dil * int2e(intindex(c, k, j, a)) - dil * int2e(intindex(c, a, j, k))&
				+ 0.5d0 * djk * int2e(intindex(c, l, i, a)) - djk * int2e(intindex(c, a, i, l))&
				+ 0.5d0 * djl * int2e(intindex(c, k, i, a)) - djl * int2e(intindex(c, a, i, k))				

	end function mat_el_2h2p_BB


	! =====================
	!    COUPLING BLOCKS
	! =====================

	! ----------------------------
	!   2h2p type 1 / 2h2p type A
	! ----------------------------

	function mat_el_2h2p_1A (i, a, b, k, l, c) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{ii}^{ab} | H | ^{A}\Psi_{kl}^{cb} >
	! -----------------
		integer :: i, k, l
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dbc, dac, dil, dik

		dab = s_mat(a, b)
		dbc = s_mat(b, c)
		dac = s_mat(a, c)
		dil = s_mat(i, l)
		dik = s_mat(i, k)

		mat_el = sqrt(1.5d0) * (dab * dil * int2e(intindex(c, k, i, b)) - dab * dik * int2e(intindex(c, l, i, b))&
                + dac * dik * int2e(intindex(b, l, i, b)) - dac * dil * int2e(intindex(b, k, i, b))& 
                + dbc * dik * int2e(intindex(b, l, i, a)) - dbc * dil * int2e(intindex(b, k, i, a))&
                - dik * int2e(intindex(c, l, i, a)) + dil * int2e(intindex(c, k, i, a)))

	end function mat_el_2h2p_1A


	! ----------------------------
	!   2h2p type 1 / 2h2p type B
	! ----------------------------

	function mat_el_2h2p_1B (i, a, b, k, l, c) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{ii}^{ab} | H | ^{B}\Psi_{kl}^{cb} >
	! -----------------
		integer :: i, k, l
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dbc, dac, dil, dik

		dab = s_mat(a, b)
		dbc = s_mat(b, c)
		dac = s_mat(a, c)
		dil = s_mat(i, l)
		dik = s_mat(i, k)

               ! October 15, 2016
               ! a minus sign added due to the difference btw definitions of |ij ab, B>
		mat_el = -sqrt(0.5d0) * &
				(dbc * dab * int2e(intindex(i, l, i, k)) + dbc * dab * int2e(intindex(i, k, i, l))& 
				+ 2.0d0 * dik * dil * int2e(intindex(c, a, b, b)) + 2.0d0 * dik * dil * int2e(intindex(c, b, b, a))& 
				- 2.0d0 * dab * dik * int2e(intindex(c, b, i, l)) + dab * dik * int2e(intindex(c, l, i, b))&
				- 2.0d0 * dab * dil * int2e(intindex(c, b, i, k)) + dab * dil * int2e(intindex(c, k, i, b))& 
				- 2.0d0 * dac * dik * int2e(intindex(b, b, i, l)) + dac * dik * int2e(intindex(b, l, i, b))&
				- 2.0d0 * dac * dil * int2e(intindex(b, b, i, k)) + dac * dil * int2e(intindex(b, k, i, b))&
				- 2.0d0 * dbc * dik * int2e(intindex(b, a, i, l)) + dbc * dik * int2e(intindex(b, l, i, a))&
				- 2.0d0 * dbc * dil * int2e(intindex(b, a, i, k)) + dbc * dil * int2e(intindex(b, k, i, a))&
				+ dac * int2e(intindex(i, k, i, l)) + dac * int2e(intindex(i, l, i, k))& 
				- 2.0d0 * dil * int2e(intindex(c, a, i, k)) + dil * int2e(intindex(c, k, i, a))& 
				- 2.0d0 * dik * int2e(intindex(c, a, i, l)) + dik * int2e(intindex(c, l, i, a)))
		

	end function mat_el_2h2p_1B


	! ----------------------------
	!   2h2p type A / 2h2p type B
	! ----------------------------

	function mat_el_2h2p_AB (i, j, a, b, k, l, c) result (mat_el)
	! -----------------
	! computes < ^{A}\Psi_{ij}^{ab} | H | ^{B}\Psi_{kl}^{cb} >
	! -----------------
		integer :: i, j, k, l
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dac, dbc, dil, dik, djl, djk

		dab = s_mat(a, b)
		dbc = s_mat(b, c)
		dac = s_mat(a, c)
		djl = s_mat(j, l)
		djk = s_mat(j, k)
		dil = s_mat(i, l)
		dik = s_mat(i, k)

               ! October 15, 2016
               ! a minus sign added due to the difference btw definitions of |ij ab, B>
		mat_el = -sqrt(3.0d0) * 0.5d0 * (&
				+ dab * dik * int2e(intindex(c, l, j, b))& 
				+ dab * dil * int2e(intindex(c, k, j, b))&
				- dab * djk * int2e(intindex(c, l, i, b))&
				- dab * djl * int2e(intindex(c, k, i, b))&
				+ dac * dik * int2e(intindex(b, l, j, b))& 
				+ dac * dil * int2e(intindex(b, k, j, b))& 
				- dac * djk * int2e(intindex(b, l, i, b))& 
				- dac * djl * int2e(intindex(b, k, i, b))& 
				- dbc * dik * int2e(intindex(b, l, j, a))&
				- dbc * dil * int2e(intindex(b, k, j, a))& 
				+ dbc * djk * int2e(intindex(b, l, i, a))& 
				+ dbc * djl * int2e(intindex(b, k, i, a))& 
				- dik * int2e(intindex(c, l, j, a))& 
				- dil * int2e(intindex(c, k, j, a))& 
				+ djk * int2e(intindex(c, l, i, a))&
				+ djl * int2e(intindex(c, k, i, a)))


	end function mat_el_2h2p_AB


	! ----------------------------
	!   1h1p / 2h2p type 1
	! ----------------------------

	function mat_el_1h1p_2h2p_1 (i, a, j, b, c) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{i}^{a} | H | ^{1}\Psi_{jj}^{bc} >
	! -----------------
		integer :: i, j
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dac, dij

		dab = s_mat(a, b)
		dac = s_mat(a, c)
		dij = s_mat(i, j)

		mat_el = 0.5d0 * (-2.0d0 * dij * int2e(intindex(a, c, j, b)) -2.0d0 * dij * int2e(intindex(a, b, j, c)) & 
                + dac * int2e(intindex(j, b, j, i)) + dac * int2e(intindex(j, i, j, b))& 
                + dab * int2e(intindex(j, i, j, c)) + dab * int2e(intindex(j, c, j, i)))
					
	end function mat_el_1h1p_2h2p_1


	! ----------------------------
	!   1h1p / 2h2p type A
	! ----------------------------

	function mat_el_1h1p_2h2p_A (i, a, j, k, b, c) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{i}^{a} | H | ^{A}\Psi_{jk}^{bc} >
	! -----------------
		integer :: i, j, k
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dac, dij, dik

		dab = s_mat(a, b)
		dac = s_mat(a, c)
		dij = s_mat(i, j)
		dik = s_mat(i, k)

		mat_el = sqrt(1.5d0) *&
				( dac * int2e(intindex(k, i, j, b)) - dac * int2e(intindex(k, b, j, i))& 
				+ dik * int2e(intindex(a, b, j, c)) - dik * int2e(intindex(a, c, j, b))&
				+ dij * int2e(intindex(a, c, k, b)) - dij * int2e(intindex(a, b, k, c))&
				+ dab * int2e(intindex(k, c, j, i)) - dab * int2e(intindex(k, i, j, c)))
				
	end function mat_el_1h1p_2h2p_A


	! ----------------------------
	!   1h1p / 2h2p type B
	! ----------------------------

	function mat_el_1h1p_2h2p_B (i, a, j, k, b, c) result (mat_el)
	! -----------------
	! computes < ^{1}\Psi_{i}^{a} | H | ^{B}\Psi_{jk}^{bc} >
	! -----------------
		integer :: i, j, k
		integer :: a, b, c
		double precision :: mat_el
		double precision :: dab, dac, dij, dik

		dab = s_mat(a, b)
		dac = s_mat(a, c)
		dij = s_mat(i, j)
		dik = s_mat(i, k)

               ! October 15, 2016
               ! a minus sign added due to the difference btw definitions of |ij ab, B>
		mat_el = -sqrt(0.5d0) * (dac * int2e(intindex(k, i, j, b)) + dac * int2e(intindex(k, b, j, i))& 
				+ dab * int2e(intindex(k, c, j, i)) + dab * int2e(intindex(k, i, j, c))& 
				- dik * int2e(intindex(a, c, j, b)) - dik * int2e(intindex(a, b, j, c))& 
				- dij * int2e(intindex(a, b, k, c)) - dij * int2e(intindex(a, c, k, b)))
		

	end function mat_el_1h1p_2h2p_B

end module mat_element_pol
