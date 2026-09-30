from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

src = OpenDataFile(r"I:/Un_Hifire_fine/E-4/paraview_harmonic_f0_3d/Harmonic_f0_666p67_Wall_Smooth.vts")
view = GetActiveViewOrCreate("RenderView")
disp = Show(src, view)
disp.Representation = "Surface"
ColorBy(disp, ("POINTS", "P_f0_amplitude_smooth"))
ampLUT = GetColorTransferFunction("P_f0_amplitude_smooth")
ampLUT.RGBPoints = [0.0, 0.10, 0.20, 0.70, 0.0015, 0.55, 0.75, 1.00, 0.0030, 1.00, 0.93, 0.60, 0.0045, 0.70, 0.05, 0.05]
ampLUT.ColorSpace = "RGB"
ampLUT.NanColor = [0.6, 0.6, 0.6]
disp.RescaleTransferFunctionToDataRange(False, True)
ampLUT.RescaleTransferFunction(0.0, 0.0045)
ampBar = GetScalarBar(ampLUT, view)
ampBar.Title = "P_f0 amplitude (smooth)"
ampBar.ComponentTitle = ""
ampBar.TitleFontSize = 16
ampBar.LabelFontSize = 14
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 1
view.Update()
ResetCamera(view)

# To display phase with the same smoothed wall file, run these lines in ParaView Python Shell:
# ColorBy(disp, ("POINTS", "P_f0_phase_continuous_display"))
# phaseLUT = GetColorTransferFunction("P_f0_phase_continuous_display")
# phaseLUT.RGBPoints = [-6.0, 0.12, 0.22, 0.85, 0.0, 0.95, 0.95, 0.95, 6.0, 0.75, 0.05, 0.05]
# phaseLUT.ColorSpace = "Diverging"
# phaseLUT.RescaleTransferFunction(-6.0, 6.0)
# GetScalarBar(phaseLUT, view).Title = "P_f0 phase trend (continuous)"
