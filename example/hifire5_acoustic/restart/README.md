# Restart files

`unsteady/` contains the four flow files and four fitted-shock files used with
the included 90 x 81 x 160 acoustic setup. A steady restart is not included:
the candidate one-rank file was excluded because its hash matches the
parabolic-leading-edge LNS restart.

The raw fine-grid archive contains all 192 flow/shock restart pairs associated
with the 160 x 151 x 160 base-flow result. Their combined size is too large for
the ordinary Git tree, so they should be stored with the versioned data DOI.
They have not been mixed with the four-rank legacy restart included here.
