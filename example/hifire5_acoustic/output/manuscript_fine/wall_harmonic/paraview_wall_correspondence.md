# ParaView 壁面谐波场与二维傅里叶图的对应关系

本目录中的 `Harmonic_f0_666p67_Wall_Smooth.vts` 是将后半段压力扰动 `P_pert` 在 `f0 = 666.667 Hz` 上的 Fourier 复系数重新贴回三维壁面 `j=1` 后得到的单一壁面结构网格。

二维图 `late_wall_pressure_harmonic_amplitude_map.png` 的横坐标 `X` 与该 VTK 壁面上的空间坐标 `X` 对应，纵坐标 circumferential index `k` 与 VTK 数组 `wall_k_index` 对应。因此二维图中的任意点 `(X,k)`，在 ParaView 中就是壁面上相同 `X` 位置和相同 `wall_k_index` 的点。

二维幅值图中的 `|p'_{f0}|` 对应 ParaView 数组 `P_f0_amplitude_smooth`。若需要和二维图色标的 `10^{-3}` 标注一致，也可以显示 `P_f0_amplitude_x1e3_smooth`。

二维相位图中的 phase(rad) 对应 ParaView 数组 `P_f0_phase_rad_smooth`。该数组仍然限制在 `[-pi, pi]`，适合与二维相位图直接对比，但在三维表面上会保留相位包裹跳变。若主要目的是展示连续传播趋势，建议显示 `P_f0_phase_continuous_display`；若需要完整未包裹相位，可显示 `P_f0_phase_unwrapped_x_rad_smooth`。

四条参考轴线的周向索引为：Long axis +: `k=41`，Long axis -: `k=121`，Short axis +: `k=1`，Short axis -: `k=81`。在 ParaView 中可以用 `Threshold` 或 `Plot Over Line` 结合 `wall_k_index` 提取这些轴线，与二维图中的黑色实线/虚线标注对应。

本次新增的平滑不是直接平滑相位角，而是先对 Fourier 复系数 `P_f0_real + i P_f0_imag` 在壁面 `(i,k)` 上做 5x5 加权平滑，再重新计算幅值和相位，因此可以避免 `-pi/pi` 相位跳变处的错误平均。
