# Preparation

- Restart system, configure displays and end unnecessary processes.

- Start Host (`Alt + H`) in **Command**. If **Environment** is used to present stimuli, check and set **Command** IP address in **Environment** so it can connect to **Command**.

- Initialize parameters in **ConditionTest** which should be the first one in experiment list when **Command** started. This experiment has no parameters been marked as inheritable(Green), so we can set global parameters here once, then all subsequent experiments will inherite proper parameters.

    Experiment Panel:
    - Experimenter, Subject, DataDir, Hemisphere, Eye,  RecordSession, RecordSite, Display, etc.

    Environment Panel:
    - ScreenToEye, ScreenHeight, ScreenAspect, CLUT, MarkerSize, MarkerCorner, etc.

- Start Experiment (`Alt + E`) **ConditionTest** and marker will flip periodicly. Restart threshold learning of photodiode logic board.

    Measure the `RiseLag` and `FallLag` of the photodiode digital signal relative to the photodiode analog signal, and set in the configuration of **Command**. Digital GPIO Synchronization is connected to **Record System** as `Channel 0`, and Digital Photodiode Synchronization as `Channel 1`.

# ISI Procedure

0. Focus and illuminate evenly on cortex. Adjust depth of field(shallow) and image dynamic range. Protect optics from movement and external light.  Save `bloodvessel` image under green and red light.

0. **ISIEpochOri8**, episodic drifting square gratings for ocular dominence map, direction map, and orientation map.
    - Left Eye
    - Right Eye
    - Both Eyes

0. **ISICycle2Color**, temporal modulation of maximum cone isolating colors for cone-opponent functional domains.
    - Achromatic(ColorSpace=DKL, Color=X)
    - L cone isolating(ColorSpace=LMS, Color=X)
    - M cone isolating(ColorSpace=LMS, Color=Y)
    - S cone isolating(ColorSpace=LMS, Color=Z)
    - Flip Red/Green (ColorSpace=HSL, Color=RGYm)
    - Flip Blue/Yellow (ColorSpace=HSL, Color=BYYm)

0. **ISIEpochFlash2Color**, episodic maximum cone isolating colors for cone-opponent functional domains.
    - Achromatic(ColorSpace=DKL, Color=X)
    - L cone isolating(ColorSpace=LMS, Color=X)
    - M cone isolating(ColorSpace=LMS, Color=Y)
    - S cone isolating(ColorSpace=LMS, Color=Z)
    - Red/Green (ColorSpace=HSL, Color=RGYm)
    - Blue/Yellow (ColorSpace=HSL, Color=BYYm)

0. **ISICycleOri**, temporally modulate orientation of drifting square grating for direction map and orientation map.
    - CycleDirection=1
    - CycleDirection=-1

0. **ISICycleColorPlane**, temporally modulate colors on 3 DKL planes, each with CCW and CW cycle direction.
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=1)
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=-1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=-1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=-1)

0. After ISI experiments, save `bloodvessel` image under red and green light again, and finally the scale bar image.

# ISI Procedure Automation

# ISI Online Analysis

# Guide Map for EPhys Penetration


# EPhys Procedure

Here **Command** is used without **Environment**, it should be set FullScreen (`Alt + F`) and FullViewport (`Alt + V`) to present stimuli while unnecessary guides (Grid, Eye Trace, etc.) turned off (`Alt + G`).

0. Based on ISI guide map, slowly (<1mm/min) insert probe perpendicular to cortex at a chosen location, check ephys signals and mark the penetration site on the ISI guide map.

0. In **ConditionTest**, turn on **Input** and manually refine stimulus, so that most of the channels recording the region of interest showing maximum responses.
    - LeftStick = Position
    - RightStick = Size
    - LT/RT = Ori -/+
    - LB/RB = Visible Off/On
    - A + LT/RT = SpatialFreq -/+, not yet
    - B + LT/RT = TemporalFreq -/+, not yet
    - X + LT/RT = Diameter -/+, not yet
    - Y + LT/RT = Diameter -/+, not yet

0. **Flash2Color**, flip fullscreen colors for layer identification. Specific colors are used to seperate 4A, 4B and 4Cα/β. **Remember** to check and set `Eye`.
    - Achromatic(dominent eye, ColorSpace=DKL, Color=X)
    - DKL L-M axis(dominent eye, ColorSpace=DKL, Color=Y)
    - DKL S-(L+M) axis(dominent eye, ColorSpace=DKL, Color=Z)

0. **Color**, uniform hue patch onset/offset against a background.
    - HSL hues/wp with equal luminence(ColorSpace=HSL, Color=HueYm)
    - DKL hues/wp with maximum cone contrast(ColorSpace=DKL, Color=HueL0)

0. **CycleColorPlane**, temporal modulate colors on 3 DKL planes, each with CCW and CW cycle direction.
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=1)
    - IsoLum Hues(ModulateParam=DKLIsoLum, CycleDirection=-1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=1)
    - IsoLM Hues(ModulateParam=DKLIsoLM, CycleDirection=-1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=1)
    - IsoSLM Hues(ModulateParam=DKLIsoSLM, CycleDirection=-1)

0. **HartleySubspace**, static sinusoidal gratings with Ori, SpatialFreq and SpatialPhase in hartley subspace.
    - Achromatic(ColorSpace=DKL, Color=X)
    - L cone isolating(ColorSpace=LMS, Color=Xmcc)
    - M cone isolating(ColorSpace=LMS, Color=Ymcc)
    - S cone isolating(ColorSpace=LMS, Color=Zmcc)

0. **Image**, natural image set with its natural colors or modulated colors the same as in HartleySubspace.

0. **OriSF**, drifting sinusoidal gratings of Ori and SF, color modulations the same as HartleySubspace.

# EPhys Procedure Automation

ExperimentSession **ColorEPhys/ColorEPhysLite** is to automate the sequence of above experiments.

0. In **ConditionTest**, manually map receptive field and set parameters:
    - `Experimenter` to recieve notification
    - `RecordSession` and `RecordSite`
    - presenting `Eye`
    - `Position`
    - `Diameter`

0. Start the experimentsession (`Alt + S`).

# Notes
- joystick need to add color bottons

scanbox recorder
