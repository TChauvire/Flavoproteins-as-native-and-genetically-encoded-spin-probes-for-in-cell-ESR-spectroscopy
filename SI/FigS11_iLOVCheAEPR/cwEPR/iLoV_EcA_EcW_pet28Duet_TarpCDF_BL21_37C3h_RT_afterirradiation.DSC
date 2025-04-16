#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	3260.000000
XWID	119.882812
*
*	Documentational Text:
*
TITL	'iLoV_EcA_EcW_pet28Duet_TarpCDF_BL21_37C3h_RT_afterirradiation'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    04/27/19
TIME    22:44:13
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.332
A1SW    0.012
MWFQ    9.321505e+09
MWPW    0.006325
AVGS    4
RESO    shqe1023
SPTP    0.04201
RCAG    60
RCHM    1
B0MA    0.0004
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    1024
RCTC    0.04096
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3320.00 G
Delay              1.0 s
FieldFlyback       On
FieldWait          Given delay
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3254.000 G
SweepDirection     Up
SweepWidth         120.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.321505 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        60.4
Power              6.325 mW
PowerAtten         15.0 dB

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         4
NbScansDone        4
NbScansToDo        4
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       40

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           9.271e-09
ConvTime           42.01 ms
DModAmp            1.00 G
DModCalibrated     True
DModDetectSCT      First
DModEliDelay       1.0 us
DModExtLockIn      False
DModExtTrigger     False
DModFieldMod       First
DModGain           60 dB
DModIntegrator     True
DModModOutput      Internal
DModTimeConst      1.28 ms
DoubleModFreq      3.00 kHz
DoubleModPhase     0.0
DoubleMode         False
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               60 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             4.000 G
ModFreq            100.00 kHz
ModOutput          Internal
ModPhase           0.0
Offset             0.0 %
PolyCof            0,0.99652,0,0,0.00737177,0,0,-0.00559614,0,0,-2.88221e-05,0,0,1.00404e-05,0,0,3.43695e-08,0,0,-5.0404e-09,0,0,-1.4783e-11,0,0,-1.29132e-12,0
QuadMode           False
ResCenter          62.5 mm
ResLength          40 mm
Resolution         1024
Resonator          1
SctNorm            True
SctRevision        AllegroDMod
SpuExtension       False
SweepTime          43.02 s
TimeConst          40.96 ms
TimeExp            False
TuneCaps           50

*
************************************************************
