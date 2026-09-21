# Local change cleanup

Remove an item when it has been resolved and committed. Do not renumber the remaining items.

15. **Repair the `accretion1d` custom-output interface.** The committed `special_spec` path reuses the spectrum callback to write `sigma`, `nu`, and `mdot` profiles. Replace or tighten that interface so it does not pass uninitialized spectrum arguments, allocate unused `nxgrid`-squared work arrays, or overcount profiles when the domain is decomposed in `y` or `z`; also rename or document `special_spec` to reflect that it schedules general Special-module output rather than spectra.

16. **Fix the `hydro_potential` module.** Diagnose and stabilize its potential-flow evolution before recreating a sample. The discarded experimental sample grew from `umax=2.225E-01` at `t=0` to `3.277E+04` at `t=2.082` and wrote crash snapshots on all four MPI ranks; verify the nonlinear potential equation, viscous term, timestep constraints, MPI decomposition behavior, and supported diagnostics with a focused regression test.
