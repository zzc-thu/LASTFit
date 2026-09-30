# Grid files

`unsteady/` contains the 90 x 81 x 160 grid as one PVTS header and four VTS
pieces. `steady/` contains the one-piece 90 x 41 x 40 grid. These correspond to
`Config.cfg` and `Config_steady.cfg`, respectively. The manuscript-grid dataset
remains pending: the 160 x 151 x 160 result and restart archives were found,
but no matching fine-grid `INIT` directory was present beside them. Do not use
the included 90-point grid with `Config_E4_fine.cfg`.
