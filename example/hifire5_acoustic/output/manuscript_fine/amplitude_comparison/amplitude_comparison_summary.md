# E-2 / E-3 / E-4 amplitude-comparison summary

Save interval: `0.0002`; sampling frequency: `5000`.

The harmonic projection frequency is `f0 = 666.666666667`. E-2 and E-3 do not contain a local Config.cfg, so their epsilons are inferred from the folder names.

| Case | epsilon | frames | late window | dominant FFT frequency | mean long | mean short | long/short | long/epsilon | short/epsilon | long 2f0/f0 | short 2f0/f0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| E-2 | 0.05 | 49 | 25-49 | 714.286 | 0.159176 | 0.121389 | 1.31129 | 3.18352 | 2.42778 | 0.238351 | 0.27765 |
| E-3 | 0.005 | 123 | 62-123 | 691.057 | 0.0156419 | 0.0122491 | 1.27699 | 3.12838 | 2.44981 | 0.024924 | 0.0296467 |
| E-4 | 0.0005 | 45 | 23-45 | 666.667 | 0.00168501 | 0.00130097 | 1.29519 | 3.37003 | 2.60195 | 0.0417035 | 0.0428562 |

## Interpretation

- If the response is linear, `|P_f0|` should scale with `epsilon`, and `|P_f0|/epsilon` should remain nearly constant.
- The smallest-amplitude case `E-4` is used as the linear reference for the nonlinearity-index plot.
- Growth of `|P_2f0|/|P_f0|` indicates waveform distortion and nonlinear generation of higher harmonics.
- Long/short-axis differences quantify how the elliptic cone geometry modulates the same incoming fast acoustic disturbance.

## E-2

Mean long-axis `|P_f0| = 0.159176`, mean short-axis `|P_f0| = 0.121389`, so the long/short ratio is `1.311`.
The normalized responses are long `3.18352` and short `2.42778`.
The second-harmonic ratios are long `0.238351` and short `0.27765`.

## E-3

Mean long-axis `|P_f0| = 0.0156419`, mean short-axis `|P_f0| = 0.0122491`, so the long/short ratio is `1.277`.
The normalized responses are long `3.12838` and short `2.44981`.
The second-harmonic ratios are long `0.024924` and short `0.0296467`.

## E-4

Mean long-axis `|P_f0| = 0.00168501`, mean short-axis `|P_f0| = 0.00130097`, so the long/short ratio is `1.295`.
The normalized responses are long `3.37003` and short `2.60195`.
The second-harmonic ratios are long `0.0417035` and short `0.0428562`.
