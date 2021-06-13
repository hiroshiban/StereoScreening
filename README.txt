# StereoScreening
A software package for stereo vision screening

FULL UPDATED VERSION OF README IS COMING SOON.




function StereoScreening(subjID,acq,:displayfile,:stimlusfile,:gamma_table,:overwrite_flg,:force_proceed_flag)
(: is optional)

Displays rectangular planes with binocular disparities (+/- arcmings) for testing psychophysical depth
discrimination acuity based on binocular disparity.
This script shoud be used with MATLAB Psychtoolbox version 3 or above.
This is a full update of the old NearFarRectangleBehavior.m so as to be compatible with the new PTB3 toolbox
(especially new PsychImaging functions) and new 3D viewing setup with an nVidia graphics card and 3D vision tools.

Currently, you can run this script only on Windows OS since the psignifit tool (expecially executables) linked
from this script is only compatible with Windows. Before running the test, please configure the psignifit executables
so that they can be called from the system, following one of the procedures below.
1. please add a path to the psignifit executables (psignifit-bootstrap.exe, psignifit-diagnostics.exe, and psignifit-mapestimate.exe)
   to you 'PATH' environmental variable.
2. Or you can uncomment setenv() descriptions to set the environmental variable from this script.
3. Or please copy the psignifit executables to your system path such as C:\Windows\System32.

Created    : "Tue Aug 17 12:25:39 2010 ban"
Last Update: "2021-06-13 15:46:11 ban"


[acknowledgment]
This screening uses Psignifit (mpsignifit) tool develped by Ingo Frund, Valentin Hanel and Felix Wichmann,
in computing the subject discrimination threshold and psychmetric functions. 
https://github.com/wichmann-lab/psignifit/wiki
https://github.com/wichmann-lab/psignifit/archive/master.zip
We would like to express our sincere gratitudes to the authors for sharing the great tool.


[input variables]
sujID         : ID of subject, string, such as 's01'
                the directory named as subj_name (e.g. 'HB') should be located in ~/StereoScreening/Presentation/subjects/
                under which configurations files are included. By changing parameters in the configuration files,
                stimulus type, size, colors, moving speed, presentation timing etc can be manipulated as you like.
                For details, please see the files in ~/StereoScreening/Presentation/subjects/_DEFAULT_.
                If subject directory does not exist in the specific directory described above, the parameters in the
                _DEFAULT_ directory would be automatically copied as subj_name and the default parameters are used for
                stimulus presentation. you can modify the default parameters later once the files are copied and the
                script is terminated.
                !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!
                !!! if 'debug' (case insensitive) is included          !!!
                !!! in subjID string, this program runs as DEBUG mode; !!!
                !!! stimulus images are saved as *.png format at       !!!
                !!! ~/CurvatureShading/Presentation/images             !!!
                !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!
acq           : acquisition number (design file number),
                a integer, such as 1, 2, 3, ...
displayfile   : (optional) display condition file, such as 'shadow_display_fmri.m'
                as an example, please see ~/StereoScreening/subjects/_DEFAULT_/nf_display.m
stimulusfile  : (optional) stimulus condition file, such as 'shadow_stimulus_exp1.m'
                all of the stimuli in this script are generated in real-time based on
                the parameters in the stimulus file. For details, please see this
                function and the other function in ../Generation and ../Common directries.
                as an example, please see ~/StereoScreening/subjects/_DEFAULT_/nf_stimulus.m
gamma_table   : (optional) table(s) of gamma-corrected video input values (Color LookupTable).
                256(8-bits) x 3(RGB) x 1(or 2,3,... when using multiple displays) matrix
                or a *.mat file specified with a relative path format. e.g. '/gamma_table/gamma1.mat'
                The *.mat should include a variable named "gamma_table" consists of a 256x3xN matrix.
                if you use multiple (more than 1) displays and set a 256x3x1 gamma-table, the same
                table will be applied to all displays. if the number of displays and gamma tables
                are different (e.g. you have 3 displays and 256x3x!2! gamma-tables), the last
                gamma_table will be applied to the second and third displays.
                if empty, normalized gamma table (repmat(linspace(0.0,1.0,256),3,1)) will be applied.
overwrite_flg : (optional) whether overwriting pre-existing result file. if 1, the previous result
                file with the same acquisition number will be overwritten by the previous one.
                if 0, the existing file will be backed-up by adding a prefix '_old' at the tail
                of the file. 0 by default.
force_proceed_flag : (optional) whether proceeding stimulus presentatin without waiting for
                the experimenter response (e.g. presesing the ENTER key) or a trigger.
                if 1, the stimulus presentation will be automatically carried on.

!!! NOTE !!!!
displayfile & stimulusfile should be located at
./subjects/(subjID)/
like ./subjects/(subjID)/nearfar_display.m
     ./subjects/(subjID)/nearfar_stimulus.m


[output variables]
no output matlab variable.


[output files]
1. behavioral result
   stored ./subjects/(subjID)/results/(today)
   as ./subjects/(subjID)/results/(today)/(subjID)_StereoScreening_results_run_(run_num).mat


[example]
>> StereoScreening('s01',1,'nf_display.m','nf_stimulus_exp1.m')


[About displayfile]
The contents of the displayfile is as below.
(The file includs 6 lines of headers and following display parameters)

(an example of the displayfile)

************************************************************
This_is_the_display_file_for_NEAR_FAR_experiment.
Please_change_the_parameters_below.
MotionParallaxConstant.m
Programmed_by_Hiroshi_Ban___Aug_24_2015
************************************************************

% display mode, one of "mono", "dual", "dualcross", "dualparallel", "cross", "parallel", "redgreen", "greenred",
% "redblue", "bluered", "shutter", "topbottom", "bottomtop", "interleavedline", "interleavedcolumn", "propixxmono", "propixxstereo"
dparam.ExpMode='shutter';%'cross';

dparam.scrID=1; % screen ID, generally 0 for a single display setup, 1 for dual display setup

% a method to start stimulus presentation
% 0:ENTER/SPACE, 1:Left-mouse button, 2:the first MR trigger pulse (CiNet),
% 3:waiting for a MR trigger pulse (BUIC) -- checking onset of pin #11 of the parallel port,
% or 4:custom key trigger (wait for a key input that you specify as tgt_key).
dparam.start_method=1;

% a pseudo trigger key from the MR scanner when it starts, only valid when dparam.start_method=4;
dparam.custom_trigger=KbName(84);

%% keyboard settings
dparam.Key1=37; % key 1 = near
dparam.Key2=39; % key 2 = far

%% whether giving correct/incorrect feedback
% 0 = no feedback, 1 = giving correct (green fixation and high-tone sound) or incorrect(red fixation and low-tone sound) feedbacks
dparam.givefeedback=1;

% screen settings

%%% whether displaying the stimuli in full-screen mode or as is (the precise resolution), 'true' or 'false' (true)
dparam.fullscr=false;

%%% the resolution of the screen height
dparam.ScrHeight=1200;

%% the resolution of the screen width
dparam.ScrWidth=1920;

% whther skipping the PTB's vertical-sync signal test. if 1, the sync test is skipped
dparam.skip_sync_test=0;


[About stimulusfile]
The contents of the stimulusfile is as below.
(The file includs 6 lines of headers and following stimulus parameters)

(an example of the stimulusfile)

************************************************************
This_is_the_stimulus_parameter_file_for_NEAR_FAR_experiment.
Please_change_the_parameters_below.
MotionParallaxConstant.m
Programmed_by_Hiroshi_Ban___Aug_24_2015
***********************************************************

% "sparam" means "stimulus generation parameters"

%%% image size
sparam.outerRectFieldSize=[8,8]; % target stimulus size in deg, [row,col]
sparam.innerRectFieldSize=[4,4]; % target stimulus size in deg, [row,col]
sparam.gapRectFieldSize=[0,0];   % widths [row(top and bottom),col(right and left)] of the gap between inner and outer rectangles in deg (if 0, no gap). gapRectFieldSize + innerRectFieldSize < outerRectFieldSize

%%% disparity magnitude(s)
sparam.base_disparity=0;         % target base disparity in deg (if non-zero, the target plane is located to near/far side compared to the fixation plane)
sparam.disparity=[8,  4,  2,  1, 0.5, -0.5, -1, -2, -4, -8]; % target disparities in arcmin

% when sparam.reference_disparity is NaN, the task is 1AFC in which participants have to judge whether the target
% is located to near or far, while, when this value is set to a specific disparity (arcmins), the task becomes 2AFC
% in which particpants have to judge which of the planes (the first or the second) is located to near.
sparam.reference_disparity=NaN;

%%% the number of trials for each stimulus condition
sparam.numTrials=30;

%%% RDS parameters
sparam.dotRadius=[0.05,0.05]; % radius of RDS's white/black ovals
sparam.dotDens=2; % deinsity of dot in RDS image (1-100)
sparam.colors=[255,0,128]; % RDS colors [dot1,dot2,background](0-255)
sparam.oversampling_ratio=2; % oversampling_ratio for fine scale RDS images, [val]

%%% RDS noise paramters
% when sparam.noise_mode='none' (no noise condition), all parameters related to noise are ignored.
% when sparam.noise_mode='anti' (anticorrelated (left/right dot colors are flipped) RDSs are presented), only sparam.noise_ratio is used in generating RDSs.
% when sparam.noise_mode='snr'  (noises are added in disparities), please set all the noise-related paramters carefully.
sparam.noise_mode='none';  % one of 'none', 'anti', and 'snr'
sparam.noise_ratio=30;     % percentage (0-100) of noises in the RDS, used when sparam.noise_mode='anti' or sparam.noise_mode='snr'
sparam.noise_mean=0;       % noise mean in disparity (arcmin), used only when sparam.noise_mode='snr'
sparam.noise_sd=5;         % noise SD in disparity (arcmin), used only when sparam.noise_mode='snr'
sparam.noise_method='add'; % one of 'add' and 'replace', used only when sparam.noise_mode='snr'

% stimulus display durations in msec

% initial fixation duration in msec
sparam.initial_fixation_time=1000;

%%% duration in msec for each condition
sparam.condition_duration=2000;

%%% duration in msec for simulus ON period for each trial, integer (1)
sparam.stim_on_duration=1000;

%%% duration in msec for each trial
sparam.BetweenDuration=500;

%%% fixation size & color
sparam.fixsize=24; % radius in pixels
sparam.fixlinesize=[12,2]; % [height,width] of the fixation line in pixel
sparam.fixcolor=[255,255,255];

%%% background color
sparam.bgcolor=[128,128,128];

%%% RGB for background patches
sparam.patch_size=[30,30];   % background patch size, [height,width] in pixels
sparam.patch_num=[20,20];    % the number of background patches along vertical and horizontal axis
sparam.patch_color1=[255,255,255]; % 1x3 matrices
sparam.patch_color2=[0,0,0];       % 1x3 matrices

%%% size parameters
run(fullfile(fileparts(mfilename('fullpath')),'sizeparams'));
%sparam.ipd=6.4;
%sparam.pix_per_cm=57.1429;
%sparam.vdist=65;


[note on the stimulus presentation and task]
The presentation will start by pressing the start button you defined as sparam.start_method.
The stimuli are presented as below,

***************************************************
*** 1 AFC (when sparam.reference_disparity=NaN) ***
***************************************************

stim -- blank -- response -- stim -- blank -- response -- stim -- blank -- reseponse ...

The task is to judge whether the targe rectangular plane is near (closer to you compared to the fixational plane) or far.
- press key 1 or left-mouse-click when the stimulus is to near (key 1/2 are defined in the display file)
- press key 2 or right-mouse-click when the stimulus is to far.

***************************************************
*** 2 AFC (when sparam.reference_disparity is set as to be a real value (reference arcmin)) ***
***************************************************

stim(or ref) - blank - stim(or ref) - blank - response - stim(or ref) -blank -- stim(or ref) -- blank -- response ...

The task is to judge which (the first or the second) of the rectangular planes is near (closer to you compared to the fixational plane).
- press key 1 or left-mouse-click when the first stimulus is to near (key 1/2 are defined in the display file)
- press key 2 or right-mouse-click when the second stimulus is to near.


[reference]
for stmulus generation, see ../Generation & ../Common directories.
