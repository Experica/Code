# Preparation

- Restart System, turn off unnecessary processes. Pychotoolbox could set the graphcs card color lookup table which is system wide and may not be unset.

- Start Host(`Alt + H`) in **Command**, if **Environment** is used to present stimuli, check **Command** IP address, and connect to **Command**.

- Initialize Parameters in **ConditionTest**, which should be the first one in the experiment list when started. This special experiment has no parameters been marked as inheritable(Green), so we can set global parameters here once, then all subsequent experiments will inherite proper parameters. 

  Check in Experimennt Panel:
    - Experimenter, Subject Info, DataDir, Hemisphere, Eye,  RecordSession, RecordSite, Display_ID, NotifyExperimenter

  Check in Environment Panel:
    - ScreenToEye, ScreenHeight, ScreenAspect, CLUT, MarkerSize, MarkerCorner, Marker On/Off Colors.

- Start Experiment(`Alt + E`) of **ConditionTest**, which will flip marker periodicly. Power/Reset the photodiode logic board. In 5 seconds, the maximun and minimun of photodiode signal will be learned and used to set appropriate threshold to convert analog signal to digital signal in real time. 

  Measure the `RiseLag` and `FallLag` of the digital signal relative to photodiode analog signal(need to re-learn and re-measure every power/reset). Digital IO is connected to **SpikeGLX** `Channel 0` and Photodiode Logic Board to **SpikeGLX** `Channel 1`. These two Channels are identical, to provide extra robustness if any channel is corrupted.

# ISI Procedure

0. Focus and illuminate evenly on cortex, protect optics from movement and external light. Adjust histogram and save `bloodvessel` map.

0. **ISIEpochOri8**, episodic drifting square gratings for ocular dominence map, direction map, and orientation map.
    - Left Eye
    - Right Eye

0. **ISICycle2Color**, temporal modulation of maximum cone isolating colors for cone-opponent functional domains
    - Achromatic(ColorSpace=DKL, Color=X)
    - L cone isolating(ColorSpace=LMS, Color=X)
    - M cone isolating(ColorSpace=LMS, Color=Y)
    - S cone isolating(ColorSpace=LMS, Color=Z)

0. **ISICycleOri**, temporally modulate orientation of drifting square grating for direction map and orientation map.
    - CycleDirection=1
    - CycleDirection=-1

0. **ISICycleColorPlane**, temporally modulated colors on 3 DKL planes, each with CCW and CW cycle direction.
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=1)
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=-1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=-1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=-1)

# ISI Procedure Automation



# EPhys Procedure for each penetration (Turn Off/Avoid Login TeamViewer/AnyDesk)

Here **Command** is used without **Environment**, it should be set FullScreen (`Alt + F`) and FullViewport (`Alt + V`) to present stimulus, and other visual guides (Grid, Eye Fixation, etc.) turned off (`Alt + G`).

0. Based on cortical maps, insert the probe at a chosen location, test probe signals and mark the penetration site on the cortical maps.

0. In **ConditionTest**, manually refine stimulus, so that most of the electrodes showing maximum response.
    - LeftStick = Position
    - RightStick = Size
    - LT/RT = Ori -/+
    - LB/RB = Visible Off/On
    - A + LT/RT = SpatialFreq -/+, not yet
    - B + LT/RT = TemporalFreq -/+, not yet
    - X + LT/RT = Diameter -/+, not yet
    - Y + LT/RT = Diameter -/+, not yet

0. **Flash2Color**, flip fullscreen colors for CSD analysis. Maximum isolating colors are used to seperate 4A, 4B and 4Ca/b. **Remember** to change and set `Eye`.
    - Achromatic(non-dominent eye, ColorSpace=DKL, Color=X)
    - Achromatic(dominent eye, ColorSpace=DKL, Color=X)
    - DKL L-M axis(dominent eye, ColorSpace=DKL, Color=Y)
    - DKL S-(L+M) axis(dominent eye, ColorSpace=DKL, Color=Z) 

0. **Color**, this is hue onset/offset against a background.
    - HSL hues/wps with equal physical luminence(ColorSpace=HSL, Color=HueYm)
    - DKL hues/wps with maximum cone contrast(ColorSpace=DKL, Color=HueL0)

0. **CycleColorPlane**, this is temporal modulated colors on 3 DKL planes, each with CCW and CW cycle direction.
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=1)
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=-1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=-1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=-1)

0. **HartleySubspace**, this is static sinusoidal gratings with Ori, SpatialFreq and SpatialPhase in hartley space.
    - Achromatic(ColorSpace=DKL, Color=X)
    - L cone isolating(ColorSpace=LMS, Color=Xmcc)
    - M cone isolating(ColorSpace=LMS, Color=Ymcc)
    - S cone isolating(ColorSpace=LMS, Color=Zmcc)

0. **Image**, natural image set, 4 color modulations the same as HartleySubspace.

0. **OriSF**, drifting sinusoidal grating of Ori and SF, 4 color modulations the same as HartleySubspace, with an additional Red/Blue square grating(ColorSpace=HSL, Color=RBYm)

# EPhys Procedure Automation for each penetration (Turn Off/Avoid Login TeamViewer/AnyDesk)

There are experiment sessions **ColorEPhys/ColorEPhysLite** to automate sequence of above experiments. To use it properly:

1. In **ConditionTest**, manually map receptive field, and set parameters
    - `Experimenter` to recieve notification
    - `RecordSession` and `RecordSite`
    - dominent `Eye`
    - `Position`
    - `Diameter`

0. In **Flash2Color**
    - mannually switch to non-doniment eye
    - set `Eye` to non-dominent
    - `ColorSpace` = DKL, `Color` = X
    - click `FullViewportSize` to make stimulus cover full viewport 

0. Switch back to dominent eye, then start **ColorEPhys/ColorEPhysLite** experiment session (`Alt + S`).

# Problems
- joystick need to add color bottons
- find ways to get HSL Same Physical luminance hue automatically for each display



scanbox recorder
