module reduce_moint
    use globals
    use intindex_module
    use read_hf_data
    implicit none
    
    contains
    
    subroutine reduce_int2e ()
        integer :: temp
        integer :: i, j, k, l
        integer*8 :: tot_index, i1
        integer*8 :: n_int2e_total
        double precision :: t_2eint_in, t_2eint_fi
        double precision, dimension(:), allocatable :: int2e_temp

        temp = n_mo + 1
        n_int2e_total = temp * (temp+1) * (temp*temp + temp + 2)/8

        allocate(int2e_temp(n_int2e_total))

        call cpu_time(t_2eint_in)
        call read_2e_integrals(n_int2e_total, int2e_temp)
        call cpu_time(t_2eint_fi)
        
        write(*,'(A27,F10.2,A2)')"Reading 2e integrals took ",t_2eint_fi-t_2eint_in,"s"   
        write(*,'(A30,I15)')"Total number of 2e integrals ",n_int2e_total
        call cpu_time(t_2eint_in)
        
        do i = 1, n_occ
            do j = i, n_occ
                do k = 1, n_occ
                    do l = k, n_mo
                        tot_index = intindex(i,j,k,l)
                        int2e(tot_index) = int2e_temp(tot_index)
                    end do
                end do
            end do
        end do

        
        do i = 1, n_occ
            do j = i, n_occ
                do k = n_occ+1, n_mo
                    do l = k, n_mo
                        tot_index = intindex(i,j,k,l)
                        int2e(tot_index) = int2e_temp(tot_index)
                    end do
                end do
            end do
        end do

        
        do i = 1, n_occ
            do j = n_occ+1, n_mo
                do k = i, n_occ
                    do l = j, n_mo
                        tot_index = intindex(i,j,k,l)
                        int2e(tot_index) = int2e_temp(tot_index)
                    end do
                end do
            end do
        end do

        
        open(329,file="moint_red.txt")
        do i1 = 1, n_int2e
            if ( int2e(i1) .ne. 0.000000000000000) then
                write(329,'(A3,I15,A1,F20.15)')" ",i1," ",int2e(i1)
            end if
        end do
        close(329)

        
        deallocate(int2e_temp)
        call cpu_time(t_2eint_fi)
        write(*,'(A26,F8.2,A2)')"2e integral sorting took ", t_2eint_fi-t_2eint_in, " s"
        
    end subroutine reduce_int2e

end module reduce_moint
