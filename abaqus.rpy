# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022 replay file
# Internal Version: 2021_09_15-18.57.30 176069
# Run by engs3052 on Fri May 29 11:13:09 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=185.216873168945, 
    height=63.033332824707)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
o1 = session.openOdb(
    name='/home/engs3052/Downloads/exp_allison_model/AbaqusJob.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/engs3052/Downloads/exp_allison_model/AbaqusJob.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          2
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
leaf = dgo.LeafFromElementSets(elementSets=("PART-1-1.VENT_SURF", ))
dg = session.DisplayGroup(leaf=leaf, name='DisplayGroup-2')
leaf = dgo.LeafFromElementSets(elementSets=("PART-1-1.VENT_SURF", ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.replace(leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=498.215, 
    farPlane=658.256, width=679.911, height=246.317, cameraPosition=(274.13, 
    220.621, -430.607), cameraUpVector=(-0.818669, 0.540477, 0.194077), 
    cameraTarget=(3.02464, 0.750938, 28.0589))
session.viewports['Viewport: 1'].view.setValues(nearPlane=494.181, 
    farPlane=659.938, width=674.406, height=244.323, cameraPosition=(-214.509, 
    340.028, -384.983), cameraUpVector=(-0.616155, 0.000710711, 0.787624), 
    cameraTarget=(1.45776, 1.13384, 28.2052))
session.viewports['Viewport: 1'].view.setValues(nearPlane=488.008, 
    farPlane=664.57, width=665.982, height=241.271, cameraPosition=(-354.315, 
    454.697, 41.3831), cameraUpVector=(-0.234097, -0.628063, 0.742116), 
    cameraTarget=(1.29344, 1.26862, 28.7063))
session.viewports['Viewport: 1'].view.setValues(nearPlane=507.38, 
    farPlane=644.188, width=692.419, height=250.849, cameraPosition=(-10.7441, 
    289.426, 527.68), cameraUpVector=(0.0066804, -0.983129, 0.182793), 
    cameraTarget=(1.23844, 1.29507, 28.6284))
session.viewports['Viewport: 1'].view.setValues(nearPlane=519.035, 
    farPlane=631.895, width=708.324, height=256.611, cameraPosition=(2.23374, 
    -93.5032, 597.021), cameraUpVector=(-0.00335803, -0.874805, -0.484464), 
    cameraTarget=(1.22497, 1.69238, 28.5564))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=20 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=19 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=18 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=17 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=16 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=20 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=499.082, 
    farPlane=652.365, width=681.096, height=246.746, cameraPosition=(548.086, 
    65.3735, 196.88), cameraUpVector=(-0.423486, -0.741599, 0.520278), 
    cameraTarget=(0.355438, 1.43929, 29.1938))
session.viewports['Viewport: 1'].view.setValues(nearPlane=485.433, 
    farPlane=665.484, width=662.469, height=239.998, cameraPosition=(410.682, 
    -381.269, -101.308), cameraUpVector=(-0.539014, -0.338859, 0.771128), 
    cameraTarget=(0.512325, 1.94931, 29.5343))
session.viewports['Viewport: 1'].view.setValues(nearPlane=485.862, 
    farPlane=664.563, width=663.055, height=240.21, cameraPosition=(-288.703, 
    -463.713, 204.985), cameraUpVector=(-0.242003, 0.780858, 0.57593), 
    cameraTarget=(1.63364, 2.08148, 29.0432))
session.viewports['Viewport: 1'].view.setValues(nearPlane=502.668, 
    farPlane=647.734, width=685.99, height=248.519, cameraPosition=(-375.38, 
    -74.4857, 457.716), cameraUpVector=(0.0663872, 0.973829, -0.217372), 
    cameraTarget=(1.80976, 1.29063, 28.5297))
del session.viewports['Viewport: 1']
#* CanvasError: SystemError: the current viewport may not be deleted.
