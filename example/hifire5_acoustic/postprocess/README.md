# HIFiRE-5 post-processing

These MATLAB scripts reproduce the reduced time-domain, harmonic, spectral,
POD, and DMD products stored under `../output/manuscript_fine/`.

The raw E-2/E-3/E-4 time histories are too large for the ordinary Git tree.
Set `LASTFIT_HIFIRE_RAW_DIR` to an extracted archive with this layout:

```text
raw-root/
  E-2/Pert/
  E-3/Pert/
  E-4/Config.cfg
  E-4/Pert/
  E-4/RESU/
  RESU/Result_00000179.pvts
```

For example, in MATLAB:

```matlab
setenv('LASTFIT_HIFIRE_RAW_DIR', 'D:\data\Un_Hifire_fine');
analyze_e4_unsteady;
analyze_e4_late_periodic;
compare_e_amplitude_effects;
compare_e_modal_analysis;
```

`export_e4_harmonic_3d_paraview` writes the complete three-dimensional
harmonic field to the raw E-4 directory. The Git release retains the complete
wall-only PVTS field because it is the compact dataset needed for the wall
amplitude and phase maps.
