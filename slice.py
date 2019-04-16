#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'EnSight Reader' and assign it to a variable, 'transientcase'
transientcase = EnSightReader(CaseFileName='/home/kaiming/Documents/ZJU_Projects/Jet/data/transient.case')
transientcase.PointArrays = ['v', 'density', 'pressure', 'temperature']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [788, 837]

# get color transfer function/color map for 'density'
densityLUT = GetColorTransferFunction('density')

# show data in view
transientcaseDisplay = Show(transientcase, renderView1)
# trace defaults for the display properties.
transientcaseDisplay.ColorArrayName = ['POINTS', 'density']
transientcaseDisplay.LookupTable = densityLUT
transientcaseDisplay.GlyphType = 'Arrow'
transientcaseDisplay.ScalarOpacityUnitDistance = 0.0016420380639339577

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
transientcaseDisplay.SetScalarBarVisibility(renderView1, False)

# get opacity transfer function/opacity map for 'density'
densityPWF = GetOpacityTransferFunction('density')

# reset view to fit data
renderView1.ResetCamera()



#################
## slice
################
# create a new 'Slice'
slice1 = Slice(Input=transientcase)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-0.21849990739250558, 0.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.ColorArrayName = ['POINTS', 'density']
slice1Display.LookupTable = densityLUT
slice1Display.GlyphType = 'Arrow'

# hide data in view
Hide(transientcase, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(transientcase)

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [-0.21849990739250558, 0.0, 0.9599974287600782]
renderView1.CameraFocalPoint = [-0.21849990739250558, 0.0, 0.0]
renderView1.CameraParallelScale = 0.2494169143257604

# *****************
# change legend layout, and its font color, position
# *****************

# get color legend for 'densityLUT' in view 'renderView1'
densityLUTColorBar = GetScalarBar(densityLUT, renderView1)

# Properties modified on densityLUTColorBar
densityLUTColorBar.AutoOrient = 0

## legend orientation
densityLUTColorBar.Orientation = 'Horizontal'

## legend normalized position
densityLUTColorBar.Position = [0.3, 0.2]

# change label color to 'black' 
densityLUTColorBar.LabelColor = [0.0, 0.0, 0.0] 

#  change titile color to 'black' 
densityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
# ***************

# set Background color as 'White'
renderView1.Background =[1,1,1]


# get layout
viewLayout1 = GetLayout()


# save screenshot
SaveScreenshot('/home/kaiming/Documents/ZJU_Projects/Jet/paraview/d_64.png', layout=viewLayout1, magnification=1, quality=100)



#################
# pressure contour
##################

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(slice1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'pressure'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'pressure'
pressureLUT = GetColorTransferFunction('pressure')

# get opacity transfer function/opacity map for 'pressure'
pressurePWF = GetOpacityTransferFunction('pressure')


# set active source
SetActiveSource(transientcase)
# *****************
# change legend layout, and its font color, position
# *****************

# get color legend for 'densityLUT' in view 'renderView1'
pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)

# Properties modified on vLUTColorBar
pressureLUTColorBar.AutoOrient = 0

## legend orientation
pressureLUTColorBar.Orientation = 'Horizontal'

## legend normalized position
pressureLUTColorBar.Position = [0.3, 0.2]

# change label color to 'black' 
pressureLUTColorBar.LabelColor = [0.0, 0.0, 0.0] 

#  change titile color to 'black' 
pressureLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
# ***************


# current camera placement for renderView1
renderView1.CameraPosition = [-0.21849990739250558, 0.0, 0.9599974287600782]
renderView1.CameraFocalPoint = [-0.21849990739250558, 0.0, 0.0]
renderView1.CameraParallelScale = 0.2494169143257604

# save screenshot
SaveScreenshot('/home/kaiming/Documents/ZJU_Projects/Jet/paraview/p_64.png', layout=viewLayout1, magnification=1, quality=100)


######################
# #temperature contour
######################

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(slice1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'temperature'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function (color map) for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')

# get opacity transfer function/opacity map for 'temperature'
temperaturePWF = GetOpacityTransferFunction('temperature')

# ******
# legend layout, and its font color, position
# *****

temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)

# Properties modified on vLUTColorBar

temperatureLUTColorBar.AutoOrient = 0

## legend orientation
temperatureLUTColorBar.Orientation = 'Horizontal'

## legend normalized position
temperatureLUTColorBar.Position = [0.3, 0.2]

# label color to 'black' 
temperatureLUTColorBar.LabelColor = [0.0, 0.0, 0.0] 

#  change 'titile' color to 'black' 
temperatureLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
# *****************

# set active source
SetActiveSource(transientcase)

# current camera placement for renderView1
renderView1.CameraPosition = [-0.21849990739250558, 0.0, 0.9599974287600782]
renderView1.CameraFocalPoint = [-0.21849990739250558, 0.0, 0.0]
renderView1.CameraParallelScale = 0.2494169143257604

# save screenshot
SaveScreenshot('/home/kaiming/Documents/ZJU_Projects/Jet/paraview/t_64.png', layout=viewLayout1, magnification=1, quality=100)


# *****************
#  velocity
# ****************

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(slice1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'v'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'v'
vLUT = GetColorTransferFunction('v')

# get opacity transfer function/opacity map for 'v'
vPWF = GetOpacityTransferFunction('v')

# *****************
# change legend layout, and its font color, position
# *****************

# get color legend for 'densityLUT' in view 'renderView1'
vLUTColorBar = GetScalarBar(vLUT, renderView1)

# Properties modified on vLUTColorBar
vLUTColorBar.AutoOrient = 0

## legend orientation
vLUTColorBar.Orientation = 'Horizontal'

## legend normalized position
vLUTColorBar.Position = [0.3, 0.2]

# change label color to 'black' 
vLUTColorBar.LabelColor = [0.0, 0.0, 0.0] 

#  change titile color to 'black' 
vLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
# ***************

# set active source
SetActiveSource(transientcase)

# current camera placement for renderView1
renderView1.CameraPosition = [-0.21849990739250558, 0.0, 0.9599974287600782]
renderView1.CameraFocalPoint = [-0.21849990739250558, 0.0, 0.0]
renderView1.CameraParallelScale = 0.2494169143257604

# save screenshot
SaveScreenshot('/home/kaiming/Documents/ZJU_Projects/Jet/paraview/v_64.png', layout=viewLayout1, magnification=1, quality=100)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-0.21849990739250558, 0.0, 0.9599974287600782]
renderView1.CameraFocalPoint = [-0.21849990739250558, 0.0, 0.0]
renderView1.CameraParallelScale = 0.2494169143257604

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
