# Reference Grids

[Case guide](../README.md)

| Dataset | Configuration | Grid | VTS pieces |
|---|---|---:|---:|
| [Unsteady grid](unsteady/Initial_grid.pvts) | `Config.cfg` | 90 x 81 x 160 | 4 |
| [Steady grid](steady/Initial_grid.pvts) | `Config_steady.cfg` | 90 x 41 x 40 | 1 |

These files are reference output, not runtime mesh inputs. The production
160 x 151 x 160 grid is not distributed here. Neither compact grid is
compatible with the production configuration or 192-rank restart.
