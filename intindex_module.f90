module intindex_module
	implicit none
	contains
	
	function intindex(i1,j1,k1,l1) result(tot_index)

        integer, intent(in) :: i1, j1, k1, l1
        integer :: i, j, k, l
        integer :: it, kt
        integer*8 :: ij, kl
        integer*8 :: tot_index

        i = i1; j = j1; k = k1; l = l1

        ! swaps i and j if i < j
        if(i<j) then
            it = i
            i = j
            j = it
        endif

        ! swaps k and l if k < l
        if(k<l) then
            kt = k
            k = l
            l = kt
        endif

        ij = i*(i+1)/2+j
        kl = k*(k+1)/2+l

        ! swaps ij and kl if ij < kl
        if(ij<kl) then
            it = kl
            kl = ij
            ij = it
        endif

        ! computes the integral index
        tot_index = ij*(ij+1)/2+kl

    end function intindex
	
end module intindex_module
