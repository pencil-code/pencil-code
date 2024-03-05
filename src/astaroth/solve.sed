/singlepass_solve/, $ { 
              s/singlepass_solve/twopass_solve_intermediate/
              t end
              s/rk3/rk3_intermediate/
              h
              s/ *[vec]*value[^,]*,/ /
              p
              g
              s/rk3_intermediate\(.*,.*value *( *[a-zA-Z0-9_]* *)\),.*)/rk3_final\1, ac_input_step_num) )/
              t store
              b cont
              : store
              h
              w sedtmp
              : cont
              $ { a Kernel twopass_solve_final(int ac_input_step_num)\{
                  r sedtmp
                  a \}
                }
              d
              : end
           }
