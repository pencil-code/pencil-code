input int AC_frequency_bin
ComputeSteps get_source_function_and_opacity(boundconds)
{
        source_function(AC_frequency_bin)
        opacity(AC_frequency_bin)
}

ComputeSteps Qintrinsic_steps(boundconds)
{
        integrate_Q_up()
        integrate_Q_down()
}
ComputeSteps Qextrinsic_steps(boundconds)
{
        revise_Q_up()
        revise_Q_down()
        sum_up_rays()
}

