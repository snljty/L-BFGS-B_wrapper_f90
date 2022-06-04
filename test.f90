! A test of the wrapper LBFGSB_wrapper of L-BFGS-B local minimization routine

# define exit_success 0
# define exit_failure 1

# define stdin_unit 5
# define stdout_unit 6
# define stderr_unit 0

subroutine func_val(n, x, f)
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: x
    double precision, intent(out) :: f

    f = (x(1) - 2.0d0) ** 2 + (x(2) - 3.5d0) ** 2 + x(1) * x(2)

    return
end subroutine func_val

subroutine func_grad(n, x, g)
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: x
    double precision, dimension(n), intent(out) :: g

    g(1) = 2 * (x(1) - 2.0d0) + x(2)
    g(2) = 2 * (x(2) - 3.5d0) + x(1)

    return
end subroutine func_grad

subroutine print_message(n, x, f, g, task, ind_cyc)
    implicit none
    integer(kind=4), intent(in) :: n
    double precision, dimension(n), intent(in) :: x, g
    double precision, intent(in) :: f
    character(len=60), intent(in) :: task
    integer, intent(in) :: ind_cyc

    if (index(task, "FG") == 1) then
        write(*, "(a,i0)") "Iteration ", ind_cyc
        write(*, "(a)") "Function value:"
        write(*, "(f10.6)") f
        write(*, "(a)") "Gradient:"
        write(*, "(f10.6)") g
        write(*, "()")
    else if (index(task, "NEW_X") == 1) then
        write(*, "(a)") "New variables: "
        write(*, "(f10.6)") x
        write(*, "()")
    end if

    return
end subroutine print_message

program main
    use lbfgsb_module
    implicit none
    interface
        subroutine func_val(n, x, f)
            integer, intent(in) :: n
            double precision, dimension(n), intent(in) :: x
            double precision, intent(out) :: f
        end subroutine func_val
        subroutine func_grad(n, x, g)
            integer, intent(in) :: n
            double precision, dimension(n), intent(in) :: x
            double precision, dimension(n), intent(out) :: g
        end subroutine func_grad
    end interface
    ! external :: func_val, func_grad
    interface
        subroutine print_message(n, x, f, g, task, ind_cyc)
            implicit none
            integer(kind=4), intent(in) :: n
            double precision, dimension(n), intent(in) :: x, g
            double precision, intent(in) :: f
            character(len=60), intent(in) :: task
            integer, intent(in) :: ind_cyc
        end subroutine print_message
    end interface
    ! external :: print_message

    integer, parameter :: n = 2
    double precision, dimension(:), allocatable :: x
    character(len=22) :: print_format
    integer :: info

    write(print_format, "(a,i0,a)") "(", n, "(f11.7,2x))"

    allocate(x(n))
    x = 0.d0 ! initial value

    call lbfgsb(n, x, func_val, func_grad, print_message, info)
    ! call lbfgsb(n, x, func_val, func_grad, info)

    if (info < 0) then
        write(stderr_unit, "(a)") "Error! The L-BFGS-B algorithm did not converged!"
        deallocate(x)
        stop exit_failure
    else if (info > 0) then
        write(stderr_unit, "(a)") "Error! There is an internal error running the L-BFGS-B algorithm."
        deallocate(x)
        stop exit_failure
    else
        write(*, "(a)") "Converged variables:"
        write(*, print_format) x
    end if

    deallocate(x)

    stop
end program main

