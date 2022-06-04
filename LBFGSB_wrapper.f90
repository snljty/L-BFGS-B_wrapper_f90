! A wrapper of L-BFGS-B local minimization routine

# define exit_success 0
# define exit_failure 1

# define stdin_unit 5
# define stdout_unit 6
# define stderr_unit 0

module LBFGSB_module
    implicit none

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine lbfgsb                                                                                    !
    !                                                                                                      !
    ! syntax: subroutine lbfgsb(n, x, func_sub, grad_sub, print_sub, info, &                               !
    !    func_tol, grad_tol, num_corr, max_cyc, bound_stat, x_lower, x_upper, print_level)                 !
    !                                                                                                      !
    ! A wrapper of L-BFGS-B                                                                                !
    ! The original version can be found on http://users.iems.northwestern.edu/~nocedal/software.html#lbfgs !
    !                                                                                                      !
    ! This subroutine uses L-BFGS-B method to minimize the value of a function by adjusting the variables, !
    ! with the gradient provided.                                                                          !
    !                                                                                                      !
    ! variables:                                                                                           !
    !                                                                                                      !
    ! n is an integer variable.                                                                            !
    !   On entry n is the dimension of variables, should be positive.                                      !
    !   On exit n is not changed.                                                                          !
    !                                                                                                      !
    ! x is a double precision array of dimension n.                                                        !
    !   On entry x is the initial guess.                                                                   !
    !     the memory of x should be allocated before.                                                      !
    !   On exit x is the final value of variables, and should be the converged values of variables if      !
    !     info == 0                                                                                        !
    !                                                                                                      !
    ! func_sub is a subroutine to computes function value of the function to be minimized, should have the !
    !   syntax func_sub(n, x, f), where                                                                    !
    !                                                                                                      !
    !     n is an integer.                                                                                 !
    !       On entry n is the dimension of variables.                                                      !
    !       On exit n is not changed.                                                                      !
    !                                                                                                      !
    !                                                                                                      !
    !     x is a double precision array of dimension n.                                                    !
    !       On entry x is the current variables.                                                           !
    !       On exit x is not changed.                                                                      !
    !                                                                                                      !
    !     f is a double pricision.                                                                         !
    !       On entry f is not specified.                                                                   !
    !       On exit f is the function value related to the current variables.                              !
    !                                                                                                      !
    !                                                                                                      !
    ! grad_sub is a subroutine to computes gradients of the function to be minimized, should have the      !
    !   syntax grad_sub(n, x, g), where                                                                    !
    !                                                                                                      !
    !     n is an integer.                                                                                 !
    !       On entry n is the dimension of variables.                                                      !
    !       On exit n is not changed.                                                                      !
    !                                                                                                      !
    !                                                                                                      !
    !     x is a double precision arrray of dimension n.                                                   !
    !       On entry x is the current variables.                                                           !
    !       On exit x is not changed.                                                                      !
    !                                                                                                      !
    !     g is an array double pricision of dimension n.                                                   !
    !       On entry g is not specified, and can be assumed to be already allocated.                       !
    !       On exit g is the gradient of the function related to the current variables.                    !
    !                                                                                                      !
    ! print_sub is a subroutine to print intermediate information, shoule have the syntax                  !
    !   print_sub(n, x, f, g, task, ind_cyc), where                                                        !
    !                                                                                                      !
    !     n is an integer.                                                                                 !
    !       On entry n is the dimension of variables.                                                      !
    !       On exit n is not changed.                                                                      !
    !                                                                                                      !
    !     x is a double precision array of dimension n.                                                    !
    !       On entry x is the current variables.                                                           !
    !       On exit x is not changed.                                                                      !
    !                                                                                                      !
    !     f is a double pricision.                                                                         !
    !       On entry f is the current function value.                                                      !
    !       On exit f is not changed.                                                                      !
    !                                                                                                      !
    !     g is a double pricision array of dimension n.                                                    !
    !       On entry g is the current function gradient.                                                   !
    !       On exit g is not changed.                                                                      !
    !                                                                                                      !
    !     task is a character string of length 60.                                                         !
    !       On entry task is the name of the current task.                                                 !
    !       On exit task is not changed.                                                                   !
    !                                                                                                      !
    !     ind_cyc is an integer.                                                                           !
    !       On entry ind_cyc is the index of the current cycle.                                            !
    !       On exit ind_cyc is not changed.                                                                !
    !                                                                                                      !
    ! info is an integer.                                                                                  !
    !   On entry info is not specified.                                                                    !
    !   On exit info indicates the status of exit status of L-BFGS-B algorithm, where                      !
    !      == 0: exited normally.                                                                          !
    !       > 0: exited with an internal error occured.                                                    !
    !       < 0: did not converge within the specified max amount of loops.                                !
    !                                                                                                      !
    ! func_tol is a double pricision, and it is optional.                                                  !
    !   On entry func_tol is the tolerance related to the change of the function value.                    !
    !     If not specified, the default value of func_tol is 1.d-7.                                        !
    !   On exit func_tol is not changed.                                                                   !
    !                                                                                                      !
    ! grad_tol is a double pricision, and it is optional.                                                  !
    !   On entry grad_tol is the tolerance related to the maximum gradient.                                !
    !     If not specified, the default value of grad_tol is 1.d-5.                                        !
    !   On exit grad_tol is not changed.                                                                   !
    !                                                                                                      !
    ! num_corr is an integer, and it is optional.                                                          !
    !   On entry num_corr is the maximum number of variable metric corrections used to define the limited  !
    !     memory matrix. The space complexity is O(num_corr**2) when m is large. It is usually between 3   !
    !     and 10.                                                                                          !
    !     If not specified, the default value of num_corr is 5.                                            !
    !   On exit num_corr is not changed.                                                                   !
    !                                                                                                      !
    ! max_cyc is an integer, and it is optional.                                                           !
    !    On entry max_cyc is the maximum allowed amount of loops, if the number of cycles is greater than  !
    !      max_cyc, the algorithm will exit without a converged solution. If is is zero or negative,       !
    !      the algorithm will never stop until it is converged.                                            !
    !      If not specified, the default value of max_cyc is 0.                                            !
    !                                                                                                      !
    ! bound_stat is an integer array of dimension n, and it is optional.                                   !
    !   On entry bound_stat(i) == 0: x(i) is unbounded,                                                    !
    !                             1: x(i) has only a lower bound,                                          !
    !                             2: x(i) has both lower and upper bound, and                              !
    !                             3: x(i) has only an upper bound.                                         !
    !     If not specified, the default values of bound_stat are 0.                                        !
    !  On exit bound_stat is not changed.                                                                  !
    !                                                                                                      !
    ! x_lower is a double precision array of dimension n, and it is optional.                              !
    !   On entry if bound_stat(i) == 1 or 2, then x_lower is the lower bound of x(i), otherwise useless.   !
    !     If not specified, the default values of x_lower are -1.d0. If bound_stat == 0 or 3, then x_lower !
    !     does not need to be specified.                                                                   !
    !   On exit x_lower is not changed.                                                                    !
    !                                                                                                      !
    ! x_upper is a double precision array of dimension n, and it is optional.                              !
    !   On entry if bound_stat(i) == 2 or 3, then x_upper is the upper bound of x(i), otherwise useless.   !
    !     If not specified, the default values of x_upper are 1.d0. If bound_stat == 0 or 1, then x_upper  !
    !     does not need to be specified.                                                                   !
    !   On exit x_upper is not changed.                                                                    !
    !                                                                                                      !
    ! print_level is an integer.                                                                           !
    !   On entry print_level controls the verbose level of internal printing of subroutine setulb in       !
    !     original L-BFGS-B code, the larger print_level is, the more verbose it is. For details please    !
    !     see the comments of iprint of the original code. If print_level is positive, a file with name    !
    !     "iterate.dat" will be created including the detail information. If print_level is negative,      !
    !     then it will prints nothing.                                                                     !
    !     If not specified, the default value of print_level is -1.                                        !
    !   On exit print_level is not changed.                                                                !
    !                                                                                                      !
    ! End of comments of subroutine lbfgsb.                                                                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lbfgsb(n, x, func_sub, grad_sub, print_sub, info, &
        func_tol, grad_tol, num_corr, max_cyc, bound_stat, x_lower, x_upper, print_level)
        implicit none
        integer, intent(in) :: n
        double precision, intent(inout), dimension(n) :: x
        interface
            subroutine func_val(n, x, f)
                integer, intent(in) :: n
                double precision, dimension(n), intent(in) :: x
                double precision, intent(out) :: f
            end subroutine func_val
            subroutine grad_sub(n, x, g)
                integer, intent(in) :: n
                double precision, dimension(n), intent(in) :: x
                double precision, dimension(n), intent(out) :: g
            end subroutine grad_sub
        end interface
        ! external :: func_sub, grad_sub
        integer, intent(out) :: info
        interface
            subroutine print_sub(n, x, f, g, task, ind_cyc)
                implicit none
                integer(kind=4), intent(in) :: n
                double precision, dimension(n), intent(in) :: x, g
                double precision, intent(in) :: f
                character(len=60), intent(in) :: task
                integer, intent(in) :: ind_cyc
            end subroutine print_sub
        end interface
        ! external :: print_sub
        double precision, intent(in), optional :: func_tol, grad_tol
        integer, intent(in), optional :: num_corr
        integer, intent(in), optional :: max_cyc
        integer, intent(in), dimension(n), optional :: bound_stat
        double precision, intent(in), dimension(n), optional :: x_lower, x_upper
        integer, intent(in), optional :: print_level

        ! for the meaning of all variables used by subtourine setulb, please read comments of the original LBFGSB code.
        double precision :: f ! function value
        double precision, dimension(:), allocatable :: g ! gradient
        character(len=:), allocatable :: task
        double precision :: factr, pgtol
        integer :: max_cyc_use, ind_cyc
        double precision, dimension(:), allocatable :: l, u
        integer, dimension(:), allocatable :: nbd
        integer :: m
        double precision, dimension(:), allocatable :: wa
        integer, dimension(:), allocatable :: iwa
        integer :: iprint
        character(len=:), allocatable :: csave
        logical, dimension(:), allocatable :: lsave
        integer, dimension(:), allocatable :: isave
        double precision, dimension(:), allocatable :: dsave

        integer :: wa_len, iwa_len
        integer, parameter :: task_len = 60, csave_len = 60, lsave_len = 4, isave_len = 44, dsave_len = 29
        integer, parameter :: m_def = 5
        double precision, parameter :: factr_def = 1.d-7, pgtol_def = 1.d-5
        double precision, parameter :: u_def = 1.d0, l_def = -1.d0
        integer, parameter :: iprint_def = -1

        if (n <= 0) then
            write(stderr_unit, "(a,i0,a)") "Error! The dimension of the variables should be positive, but got ", n, "."
            info = 1
            return
        end if

        allocate(g(n))
        allocate(nbd(n))
        allocate(l(n), u(n))

        allocate(character(len=task_len) :: task)
        allocate(character(len=csave_len) :: csave)
        allocate(lsave(lsave_len))
        allocate(isave(isave_len))
        allocate(dsave(dsave_len))

        if (present(func_tol)) then
            factr = func_tol
        else
            factr = factr_def
        end if
        if (present(grad_tol)) then
            pgtol = grad_tol
        else
            pgtol = pgtol_def
        end if
        if (present(num_corr)) then
            m = num_corr
        else
            m = m_def
        end if
        if (present(max_cyc)) then
            max_cyc_use = max_cyc
        else
            max_cyc_use = 0 ! no max cycle limitation
        end if
        if (present(x_lower)) then
            l = x_lower
        else
            l = l_def ! will not be used if nbd == 0
        end if
        if (present(x_upper)) then
            u = x_upper
        else
            u = u_def ! will not be used if nbd == 0
        end if
        if (present(bound_stat)) then
            nbd = bound_stat
        else
            nbd = 0 ! no boundary limitation
        end if
        if (present(print_level)) then
            iprint = print_level
        else
            iprint = iprint_def
        end if

        wa_len = 2 * m * n + 5 * n + 11 * m ** 2 + 8 * m
        iwa_len = 3 * n
        allocate(wa(wa_len))
        allocate(iwa(iwa_len))

        task(:) = "START" ! if task is "character(len=:), allocatable" not "pointer", then (:) cannot be omitted!
        ind_cyc = 0

        do while (index(task, "START") == 1 .or. index(task, "FG") == 1 .or. index(task, "NEW_X") == 1)
            call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
            call print_sub(n, x, f, g, task, ind_cyc)
            if (index(task, "FG") == 1) then ! Need evaluate function value and gradient
                ind_cyc = ind_cyc + 1
                if (max_cyc_use > 0 .and. ind_cyc > max_cyc_use) then
                    exit
                end if
                call func_sub(n, x, f)
                call grad_sub(n, x, g)
                
            end if
        end do

        if (index(task, "ERROR") == 1 .or. index(task, "ABNORMAL") == 1) then
            info = 1
        else if (max_cyc_use > 0 .and. ind_cyc == max_cyc_use + 1) then
            info = -1
        else
            info = 0 ! index(task, "CONVERGENCE") == 1
        end if

        deallocate(g)
        deallocate(nbd)
        deallocate(l, u)
        deallocate(wa)
        deallocate(iwa)

        deallocate(task)
        deallocate(csave)
        deallocate(lsave)
        deallocate(isave)
        deallocate(dsave)

        return      
    end subroutine lbfgsb

end module lbfgsb_module

