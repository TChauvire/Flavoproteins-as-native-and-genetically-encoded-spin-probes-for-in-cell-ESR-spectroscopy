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
XMIN	3262.000000
XWID	119.882812
*
*	Documentational Text:
*
TITL	'241210_Tempo50uM_M=1G.DSC'
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
DATE    12/11/24
TIME    12:11:09
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.3322
A1SW    0.012
MWFQ    9.322904e+09
MWPW    0.006325
AVGS    2
RESO    shqe1023
SPTP    0.04224
RCAG    20
RCHM    1
B0MA    0.0001
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
CenterField        3322.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3262.000 G
SweepDirection     Up
SweepWidth         120.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.322904 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        60.4
Power              6.325 mW
PowerAtten         15.0 dB
QValue             3900

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         2
NbScansDone        2
NbScansToDo        2
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       40

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvTime           42.24 ms
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
Gain               20 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             1.000 G
ModFreq            100.00 kHz
ModOutput          Internal
ModPhase           0.0
Offset             0.0 %
QuadMode           False
Resolution         1024
Resonator          1
SctNorm            True
SctRevision        AllegroDMod
SpuExtension       False
SweepTime          43.25 s
TimeConst          40.96 ms
TimeExp            False
TuneCaps           50

*
************************************************************
