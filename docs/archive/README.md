# Historical Configurations

[Examples](../../example/README.md) | [Verification](../verification.md)

`Config_E4_fine.cfg` records a 160 x 151 x 160, 192-rank HIFiRE run from
another source revision. Its `ModelType=2` is not HIFiRE in the current
source, where HIFiRE uses `ModelType=4`. Do not copy this file into a run.

Migration requires checking the geometry parameters, generated wall/grid,
base flow, checkpoint partition, time and forcing phase against that
revision. Changing the model number alone is not verified migration.
The production checkpoint and histories are not distributed in this repo.
