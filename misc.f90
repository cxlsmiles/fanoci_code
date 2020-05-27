module misc
    use globals
    implicit none
    contains

    subroutine initialize_smat()
        integer*8 :: i

        s_mat(:,:) = 0.0d0

        do i = 1, n_mo
            s_mat(i,i) = 1.0d0
        end do

    end subroutine initialize_smat    


    subroutine print_mat (mat_size_1, mat_size_2, mat)
        integer*8 :: mat_size_1, mat_size_2
        integer*8 :: i, j, mat_size
        double precision, dimension(mat_size_1, mat_size_2) :: mat


        do i = 1, mat_size_1
            do j = 1, mat_size_2
                !write(*,'((f20.15))',advance="no")mat(i,j)
                write(*,'(I3,A2,I3,A2,f20.15)')i," ",j," ",mat(i,j)
            end do
            write(*,*)
        end do

    end subroutine print_mat


    subroutine print_ci_eigval (k)
    !--------------------------------------------------------
    !   The eigenvalues for a given k are written in a
    !   file "Eigval.kXX.txt", where XX is the k-value
    !--------------------------------------------------------
        integer :: k
        integer :: i, j
        character(len=100) :: format_string, filename
                
        if (k < 10) then
           format_string = "(A8,I1,A3,I1,A4)"
        else if (k .ge. 10 .and. k < 100) then
           format_string = "(A8,I2,A3,I1,A4)"
        else
           format_string = "(A8,I3,A3,I1,A4)"
        endif


        write (filename,format_string) "Eigval.k", k,".iv",h_in,".txt"
        open(unit=101, file=filename)
        
        write(718, '(F20.15,A8)', advance='no') emo(k)," "
        do i = 1, ci_mat_size
           write(101,'((F20.15))')CI_eigval(i)
           write(718,'((F20.15))', advance='no')CI_eigval(i)-emo(k)
        end do
        
        close(101)
        write(718,*)


        write(filename,format_string) "Eigvec.k",k,".iv",h_in,".txt"
                      
        open(unit=201, file=filename)
        do i = 1, ci_mat_size
            write(201, '(I3,A2,F20.15)')i,": ",CI_eigval(i)
            do j = 1, ci_mat_size
                if (abs(CI_matrix(j,i)) .ge. eig_thresh) then
                    write(201, '(A,A,F10.6)')config_2hnp_str(j)," ",CI_matrix(j,i) 
                    !, advance='no')config_2h_str(j),": ",CI_matrix(j,i),"  "
                end if
            end do
            write(201,*)
        end do
        close(201)
        
    end subroutine print_ci_eigval

    subroutine print_2h_matrix ()
        integer :: i, j

        open(unit=299,file="states.2h.txt")
        do i = 1, ci_mat_size
            write(299, '(I3,A2,F20.15)')i,": ",mat_2h_eigval(i)
            do j = 1, ci_mat_size
                if (abs(mat_2h(j,i)) .gt. 0.000001) then
                    write(299, '(A,A,F20.15)')config_2h_str(j)," ",mat_2h(j,i)
                end if
            end do
            write(299,*)
        end do
        close(299)

        open(unit=312,file="energies.2h.txt")
        do i = 1, ci_mat_size
            write(312, '(F20.15)')mat_2h_eigval(i)
        end do
        close(312)
    end subroutine print_2h_matrix

    subroutine print_2h_contributions (mat_prod,k)
        double precision, dimension(ci_mat_size,ci_mat_size) :: mat_prod
        integer :: i, j
        integer, intent(in) :: k
        character(len=100) :: format_string, filename

        if (k < 10) then
           format_string = "(A12,I1,A4)"
        else if (k .ge. 10 .and. k < 100) then
           format_string = "(A12,I2,A4)"
        else
           format_string = "(A12,I3,A4)"
        endif



        write (filename,format_string) "coeff.2h.k", k,".txt"

        open(354, file=filename)
        do i = 1, ci_mat_size ! loop over the eigenvectors chi_q^a

            write(354,'(I3,A2,F20.15)')i,": ",CI_eigval(i)
            do j = 1, ci_mat_size
                if (abs(mat_prod(j, i)) .gt. 0.d0) then
                    !write(354,*)j, mat_prod(j, i)
                    write(354,*)j, mat_prod(j, i)**2
                end if
            end do
        end do

    end subroutine print_2h_contributions
end module misc
