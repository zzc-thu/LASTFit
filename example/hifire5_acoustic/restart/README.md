# Restart Data

`unsteady/` contains four flow files and four fitted-shock files matching the
90 x 81 x 160, four-rank compact setup.

To continue the compact case, copy these files to `RESU/` and set
`IF_Continue_Calculate=1`. The 192-rank restart associated with
`Config_E4_fine.cfg` is not stored in Git and must not be replaced with these
four-rank files.
