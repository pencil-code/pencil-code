module test
        contains
        subroutine add_k(x)
                use module_k
                integer :: k
                integer :: x
                k = 5
                x = x + k
        endsubroutine add_k
        subroutine sub_a()
                use module_j
                integer :: j
                j = j + 1
                add_k(j)
        endsubroutine sub_a
endmodule test
