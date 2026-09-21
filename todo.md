# Local change cleanup

Remove an item when it has been resolved and committed. Do not renumber the remaining items.

13. **Revisit the Training import after updating from upstream.** `src/slices.f90` comments out `use Training, only: get_slices_training` while retaining the call. Leave this untouched until the relevant upstream changes have been fetched and compared.

14. **Decide how to preserve the `hydro_potential` sample.** The untracked `samples/hydro_potential/` contains 11 visible configuration/input files plus ignored build products and simulation output. Separate the reusable sample definition from generated executables, objects, snapshots, crash files, and time-series data before committing anything.

15. **Repair the `accretion1d` custom-output interface.** The committed `special_spec` path reuses the spectrum callback to write `sigma`, `nu`, and `mdot` profiles. Replace or tighten that interface so it does not pass uninitialized spectrum arguments, allocate unused `nxgrid`-squared work arrays, or overcount profiles when the domain is decomposed in `y` or `z`; also rename or document `special_spec` to reflect that it schedules general Special-module output rather than spectra.
