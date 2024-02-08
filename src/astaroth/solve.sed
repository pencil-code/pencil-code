/singlepass_solve/, $ { 
              s/singlepass_solve/twopass_solve_intermediate/
              t end
              s/rk3/rk3_intermediate/
              p
              s/rk3_intermediate\(.*\)\(, *value.*\),.*)/rk3_final\1\2))/
              t store
              b cont
              : store
              h
              w temp
              : cont
              $ { a Kernel twopass_solve_final()\{
                  r temp
                  a \}
                }
              d
              : end
           }
