function StereoScreening(subjID,acq,displayfile,stimulusfile,gamma_table,overwrite_flg,force_proceed_flag)

% function StereoScreening(subjID,acq,:displayfile,:stimlusfile,:gamma_table,:overwrite_flg,:force_proceed_flag)
% (: is optional)
%
% Displays rectangular planes with binocular disparities (+/- arcmings) for testing psychophysical depth
% discrimination acuity based on binocular disparity.
% This script shoud be used with MATLAB Psychtoolbox version 3 or above.
% This is a full update of the old NearFarRectangleBehavior.m so as to be compatible with the new PTB3 toolbox
% (especially new PsychImaging functions) and new 3D viewing setup with an nVidia graphics card and 3D vision tools.
%
% Currently, you can run this script only on Windows OS since the psignifit tool (expecially executables) linked
% from this script is only compatible with Windows. Before running the test, please configure the psignifit executables
% so that they can be called from the system, following one of the procedures below.
% 1. please add a path to the psignifit executables (psignifit-bootstrap.exe, psignifit-diagnostics.exe, and psignifit-mapestimate.exe)
%    to you 'PATH' environmental variable.
% 2. Or you can uncomment setenv() descriptions to set the environmental variable from this script.
% 3. Or please copy the psignifit executables to your system path such as C:\Windows\System32.
%
% Created    : "Tue Aug 17 12:25:39 2010 ban"
% Last Update: "2021-06-13 15:46:46 ban"
%
%
% [acknowledgment]
% This screening uses Psignifit (mpsignifit) tool develped by Ingo Frund, Valentin Hanel and Felix Wichmann,
% in computing the subject discrimination threshold and psychmetric functions. 
% https://github.com/wichmann-lab/psignifit/wiki
% https://github.com/wichmann-lab/psignifit/archive/master.zip
% We would like to express our sincere gratitudes to the authors for sharing the great tool.
%
%
% [input variables]
% sujID         : ID of subject, string, such as 's01'
%                 the directory named as subj_name (e.g. 'HB') should be located in ~/StereoScreening/Presentation/subjects/
%                 under which configurations files are included. By changing parameters in the configuration files,
%                 stimulus type, size, colors, moving speed, presentation timing etc can be manipulated as you like.
%                 For details, please see the files in ~/StereoScreening/Presentation/subjects/_DEFAULT_.
%                 If subject directory does not exist in the specific directory described above, the parameters in the
%                 _DEFAULT_ directory would be automatically copied as subj_name and the default parameters are used for
%                 stimulus presentation. you can modify the default parameters later once the files are copied and the
%                 script is terminated.
%                 !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!
%                 !!! if 'debug' (case insensitive) is included          !!!
%                 !!! in subjID string, this program runs as DEBUG mode; !!!
%                 !!! stimulus images are saved as *.png format at       !!!
%                 !!! ~/CurvatureShading/Presentation/images             !!!
%                 !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!
% acq           : acquisition number (design file number),
%                 a integer, such as 1, 2, 3, ...
% displayfile   : (optional) display condition file, such as 'shadow_display_fmri.m'
%                 as an example, please see ~/StereoScreening/subjects/_DEFAULT_/nf_display.m
% stimulusfile  : (optional) stimulus condition file, such as 'shadow_stimulus_exp1.m'
%                 all of the stimuli in this script are generated in real-time based on
%                 the parameters in the stimulus file. For details, please see this
%                 function and the other function in ../Generation and ../Common directries.
%                 as an example, please see ~/StereoScreening/subjects/_DEFAULT_/nf_stimulus.m
% gamma_table   : (optional) table(s) of gamma-corrected video input values (Color LookupTable).
%                 256(8-bits) x 3(RGB) x 1(or 2,3,... when using multiple displays) matrix
%                 or a *.mat file specified with a relative path format. e.g. '/gamma_table/gamma1.mat'
%                 The *.mat should include a variable named "gamma_table" consists of a 256x3xN matrix.
%                 if you use multiple (more than 1) displays and set a 256x3x1 gamma-table, the same
%                 table will be applied to all displays. if the number of displays and gamma tables
%                 are different (e.g. you have 3 displays and 256x3x!2! gamma-tables), the last
%                 gamma_table will be applied to the second and third displays.
%                 if empty, normalized gamma table (repmat(linspace(0.0,1.0,256),3,1)) will be applied.
% overwrite_flg : (optional) whether overwriting pre-existing result file. if 1, the previous result
%                 file with the same acquisition number will be overwritten by the previous one.
%                 if 0, the existing file will be backed-up by adding a prefix '_old' at the tail
%                 of the file. 0 by default.
% force_proceed_flag : (optional) whether proceeding stimulus presentatin without waiting for
%                 the experimenter response (e.g. presesing the ENTER key) or a trigger.
%                 if 1, the stimulus presentation will be automatically carried on.
%
% !!! NOTE !!!!
% displayfile & stimulusfile should be located at
% ./subjects/(subjID)/
% like ./subjects/(subjID)/nearfar_display.m
%      ./subjects/(subjID)/nearfar_stimulus.m
%
%
% [output variables]
% no output matlab variable.
%
%
% [output files]
% 1. behavioral result
%    stored ./subjects/(subjID)/results/(today)
%    as ./subjects/(subjID)/results/(today)/(subjID)_StereoScreening_results_run_(run_num).mat
%
%
% [example]
% >> StereoScreening('s01',1,'nf_display.m','nf_stimulus_exp1.m')
%
%
% [About displayfile]
% The contents of the displayfile is as below.
% (The file includs 6 lines of headers and following display parameters)
%
% (an example of the displayfile)
%
% ************************************************************
% This_is_the_display_file_for_NEAR_FAR_experiment.
% Please_change_the_parameters_below.
% MotionParallaxConstant.m
% Programmed_by_Hiroshi_Ban___Aug_24_2015
% ************************************************************
%
% % display mode, one of "mono", "dual", "dualcross", "dualparallel", "cross", "parallel", "redgreen", "greenred",
% % "redblue", "bluered", "shutter", "topbottom", "bottomtop", "interleavedline", "interleavedcolumn", "propixxmono", "propixxstereo"
% dparam.ExpMode='shutter';%'cross';
%
% dparam.scrID=1; % screen ID, generally 0 for a single display setup, 1 for dual display setup
%
% % a method to start stimulus presentation
% % 0:ENTER/SPACE, 1:Left-mouse button, 2:the first MR trigger pulse (CiNet),
% % 3:waiting for a MR trigger pulse (BUIC) -- checking onset of pin #11 of the parallel port,
% % or 4:custom key trigger (wait for a key input that you specify as tgt_key).
% dparam.start_method=1;
%
% % a pseudo trigger key from the MR scanner when it starts, only valid when dparam.start_method=4;
% dparam.custom_trigger=KbName(84);
%
% %% keyboard settings
% dparam.Key1=37; % key 1 = near
% dparam.Key2=39; % key 2 = far
%
% %% whether giving correct/incorrect feedback
% % 0 = no feedback, 1 = giving correct (green fixation and high-tone sound) or incorrect(red fixation and low-tone sound) feedbacks
% dparam.givefeedback=1;
%
% % screen settings
%
% %%% whether displaying the stimuli in full-screen mode or as is (the precise resolution), 'true' or 'false' (true)
% dparam.fullscr=false;
%
% %%% the resolution of the screen height
% dparam.ScrHeight=1200;
%
% %% the resolution of the screen width
% dparam.ScrWidth=1920;
%
% % whther skipping the PTB's vertical-sync signal test. if 1, the sync test is skipped
% dparam.skip_sync_test=0;
%
%
% [About stimulusfile]
% The contents of the stimulusfile is as below.
% (The file includs 6 lines of headers and following stimulus parameters)
%
% (an example of the stimulusfile)
%
% ************************************************************
% This_is_the_stimulus_parameter_file_for_NEAR_FAR_experiment.
% Please_change_the_parameters_below.
% MotionParallaxConstant.m
% Programmed_by_Hiroshi_Ban___Aug_24_2015
%************************************************************
%
% % "sparam" means "stimulus generation parameters"
%
% %%% image size
% sparam.outerRectFieldSize=[8,8]; % target stimulus size in deg, [row,col]
% sparam.innerRectFieldSize=[4,4]; % target stimulus size in deg, [row,col]
% sparam.gapRectFieldSize=[0,0];   % widths [row(top and bottom),col(right and left)] of the gap between inner and outer rectangles in deg (if 0, no gap). gapRectFieldSize + innerRectFieldSize < outerRectFieldSize
%
% %%% disparity magnitude(s)
% sparam.base_disparity=0;         % target base disparity in deg (if non-zero, the target plane is located to near/far side compared to the fixation plane)
% sparam.disparity=[8,  4,  2,  1, 0.5, -0.5, -1, -2, -4, -8]; % target disparities in arcmin
%
% % when sparam.reference_disparity is NaN, the task is 1AFC in which participants have to judge whether the target
% % is located to near or far, while, when this value is set to a specific disparity (arcmins), the task becomes 2AFC
% % in which particpants have to judge which of the planes (the first or the second) is located to near.
% sparam.reference_disparity=NaN;
%
% %%% the number of trials for each stimulus condition
% sparam.numTrials=30;
%
% %%% RDS parameters
% sparam.dotRadius=[0.05,0.05]; % radius of RDS's white/black ovals
% sparam.dotDens=2; % deinsity of dot in RDS image (1-100)
% sparam.colors=[255,0,128]; % RDS colors [dot1,dot2,background](0-255)
% sparam.oversampling_ratio=2; % oversampling_ratio for fine scale RDS images, [val]
%
% %%% RDS noise paramters
% % when sparam.noise_mode='none' (no noise condition), all parameters related to noise are ignored.
% % when sparam.noise_mode='anti' (anticorrelated (left/right dot colors are flipped) RDSs are presented), only sparam.noise_ratio is used in generating RDSs.
% % when sparam.noise_mode='snr'  (noises are added in disparities), please set all the noise-related paramters carefully.
% sparam.noise_mode='none';  % one of 'none', 'anti', and 'snr'
% sparam.noise_ratio=30;     % percentage (0-100) of noises in the RDS, used when sparam.noise_mode='anti' or sparam.noise_mode='snr'
% sparam.noise_mean=0;       % noise mean in disparity (arcmin), used only when sparam.noise_mode='snr'
% sparam.noise_sd=5;         % noise SD in disparity (arcmin), used only when sparam.noise_mode='snr'
% sparam.noise_method='add'; % one of 'add' and 'replace', used only when sparam.noise_mode='snr'
%
% % stimulus display durations in msec
%
% % initial fixation duration in msec
% sparam.initial_fixation_time=1000;
%
% %%% duration in msec for each condition
% sparam.condition_duration=2000;
%
% %%% duration in msec for simulus ON period for each trial, integer (1)
% sparam.stim_on_duration=1000;
%
% %%% duration in msec for each trial
% sparam.BetweenDuration=500;
%
% %%% fixation size & color
% sparam.fixsize=24; % radius in pixels
% sparam.fixlinesize=[12,2]; % [height,width] of the fixation line in pixel
% sparam.fixcolor=[255,255,255];
%
% %%% background color
% sparam.bgcolor=[128,128,128];
%
% %%% RGB for background patches
% sparam.patch_size=[30,30];   % background patch size, [height,width] in pixels
% sparam.patch_num=[20,20];    % the number of background patches along vertical and horizontal axis
% sparam.patch_color1=[255,255,255]; % 1x3 matrices
% sparam.patch_color2=[0,0,0];       % 1x3 matrices
%
% %%% size parameters
% run(fullfile(fileparts(mfilename('fullpath')),'sizeparams'));
% %sparam.ipd=6.4;
% %sparam.pix_per_cm=57.1429;
% %sparam.vdist=65;
%
%
% [note on the stimulus presentation and task]
% The presentation will start by pressing the start button you defined as sparam.start_method.
% The stimuli are presented as below,
%
% ***************************************************
% *** 1 AFC (when sparam.reference_disparity=NaN) ***
% ***************************************************
%
% stim -- blank -- response -- stim -- blank -- response -- stim -- blank -- reseponse ...
%
% The task is to judge whether the targe rectangular plane is near (closer to you compared to the fixational plane) or far.
% - press key 1 or left-mouse-click when the stimulus is to near (key 1/2 are defined in the display file)
% - press key 2 or right-mouse-click when the stimulus is to far.
%
% ***************************************************
% *** 2 AFC (when sparam.reference_disparity is set as to be a real value (reference arcmin)) ***
% ***************************************************
%
% stim(or ref) - blank - stim(or ref) - blank - response - stim(or ref) -blank -- stim(or ref) -- blank -- response ...
%
% The task is to judge which (the first or the second) of the rectangular planes is near (closer to you compared to the fixational plane).
% - press key 1 or left-mouse-click when the first stimulus is to near (key 1/2 are defined in the display file)
% - press key 2 or right-mouse-click when the second stimulus is to near.
%
%
% [reference]
% for stmulus generation, see ../Generation & ../Common directories.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check the input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear global; clear mex;
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(displayfile), displayfile=[]; end
if nargin<4 || isempty(stimulusfile), stimulusfile=[]; end
if nargin<5 || isempty(gamma_table), gamma_table=[]; end
if nargin<6 || isempty(overwrite_flg), overwrite_flg=0; end
if nargin<7 || isempty(force_proceed_flag), force_proceed_flag=0; end

% check the aqcuisition number
if acq<1, error('Acquistion number must be integer and greater than zero'); end

% check the subject directory
if ~exist(fullfile(pwd,'subjects',subjID),'dir'), error('can not find subj directory. check the input variable.'); end

rootDir=fileparts(mfilename('fullpath'));

% check the display/stimulus files
if ~isempty(displayfile)
  if ~strcmpi(displayfile(end-1:end),'.m'), displayfile=[displayfile,'.m']; end
  if ~exist(fullfile(rootDir,'subjects',subjID,displayfile),'file'), error('displayfile not found. check the input variable.'); end
end

if ~isempty(stimulusfile)
  if ~strcmpi(stimulusfile(end-1:end),'.m'), stimulusfile=[stimulusfile,'.m']; end
  if ~exist(fullfile(rootDir,'subjects',subjID,stimulusfile),'file'), error('stimulusfile not found. check the input variable.'); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Add paths to the subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add paths to the subfunctions
addpath(genpath(fullfile(rootDir,'..','Common')));
addpath(fullfile(rootDir,'..','Generation'));
addpath(fullfile(rootDir,'..','mpsignifit'));
%setenv('PATH',[getenv('PATH'),';',fullfile(rootDir,'..','mpsignifit')]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% For a log file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get date
today=datestr(now,'yymmdd');

% result directry & file
resultDir=fullfile(rootDir,'subjects',num2str(subjID),'results',today);
if ~exist(resultDir,'dir'), mkdir(resultDir); end

% record the output window
logfname=fullfile(resultDir,[num2str(subjID),'_stereo_screening_run_',num2str(acq,'%02d'),'.log']);
diary(logfname);
warning off; %#ok warning('off','MATLAB:dispatcher:InexactCaseMatch');

%%%%% try & catch %%%%%
try


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% check the PTB version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PTB_OK=CheckPTBversion(3); % check wether the PTB version is 3
if ~PTB_OK, error('Wrong version of Psychtoolbox is running. ColorAdaptation requires PTB ver.3'); end

% debug level, black screen during calibration
Screen('Preference', 'VisualDebuglevel', 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup random seed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InitializeRandomSeed();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reset display Gamma-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(gamma_table)
  gamma_table=repmat(linspace(0.0,1.0,256),3,1); %#ok
  GammaResetPTB(1.0);
else
  GammaLoadPTB(gamma_table);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Validate dparam (displayfile) and sparam (stimulusfile) structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize dparam
dparam=struct(); % initialize
if ~isempty(displayfile), run(fullfile(rootDir,'subjects',subjID,displayfile)); end % load specific dparam parameters configured for each of the participants
dparam=ValidateStructureFields(dparam,... % validate fields and set the default values to missing field(s)
         'ExpMode','shutter',...
         'scrID',0,...
         'start_method',1,...
         'custom_trigger',KbName(84),...
         'Key1',37,...
         'Key2',39,...
         'givefeedback',1,...
         'fullscr',false,...
         'ScrHeight',1200,...
         'ScrWidth',1920,...
         'skip_sync_test',0);

% organize sparam
sparam=struct(); % initialize
if ~isempty(stimulusfile), run(fullfile(rootDir,'subjects',subjID,stimulusfile)); end % load specific sparam parameters configured for each of the participants
sparam=ValidateStructureFields(sparam,... % validate fields and set the default values to missing field(s)
         'outerRectFieldSize',[8,8],...
         'innerRectFieldSize',[4,4],...
         'gapRectFieldSize',[0,0],...
         'base_disparity',0,...
         'reference_disparity',NaN,...
         'disparity',[8,  4,  2,  1, 0.5, -0.5, -1, -2, -4, -8],...
         'numTrials',30,...
         'dotRadius',[0.05,0.05],...
         'dotDens',2,...
         'colors',[255,0,128],...
         'oversampling_ratio',2,...
         'noise_mode','none',...
         'noise_ratio',30,...
         'noise_mean',0,...
         'noise_sd',5,...
         'noise_method','add',...
         'initial_fixation_time',1000,...
         'condition_duration',2000,...
         'stim_on_duration',1000,...
         'BetweenDuration',500,...
         'fixsize',24,...
         'fixlinesize',[12,2],...
         'fixcolor',[255,255,255],...
         'bgcolor',[128,128,128],...
         'patch_size',[30,30],...
         'patch_num',[20,20],...
         'patch_color1',[255,255,255],...
         'patch_color2',[0,0,0],...
         'ipd',6.4,...
         'pix_per_cm',57.1429,...
         'vdist',65);

% change unit from msec to sec.
sparam.initial_fixation_time = sparam.initial_fixation_time./1000;
sparam.condition_duration    = sparam.condition_duration./1000;
sparam.BetweenDuration       = sparam.BetweenDuration./1000;
sparam.stim_on_duration      = sparam.stim_on_duration./1000;
sparam.stim_off_duration     = sparam.condition_duration - sparam.stim_on_duration;

% set the number of conditions
sparam.numConds=numel(sparam.disparity);

% set the other parameters
dparam.RunScript = mfilename();
sparam.RunScript = mfilename();

% validating the parameters
if sparam.gapRectFieldSize+sparam.innerRectFieldSize >= sparam.outerRectFieldSize
  error('sparam.gapRectFieldSize + sparam.innerRectFieldSize should be smaller than sparam.outerRectFieldSize. check the contents of sparam.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying the presentation parameters you set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('The Presentation Parameters are as below.\n\n');
fprintf('************************************************\n');
fprintf('Date and time          : %s %s\n',datestr(now,'yyyy-mm-dd(DDD) HH:MM:SS'),getusername());
fprintf('****** Script, Subject, Acquistion Number ******\n');
fprintf('Running Script Name    : %s\n',mfilename());
fprintf('Subject ID             : %s\n',subjID);
fprintf('Acquisition Number     : %d\n',acq);
fprintf('********* Run Type, Display Image Type *********\n');
fprintf('Display Mode           : %s\n',dparam.ExpMode);
fprintf('use Full Screen Mode   : %d\n',logical(dparam.fullscr));
fprintf('Start Method           : %d\n',dparam.start_method);
if dparam.start_method==4
  fprintf('Custom Trigger         : %s\n',dparam.custom_trigger);
end
fprintf('*************** Screen Settings ****************\n');
fprintf('Screen Height          : %d\n',dparam.ScrHeight);
fprintf('Screen Width           : %d\n',dparam.ScrWidth);
fprintf('*********** Stimulation Periods etc. ***********\n');
fprintf('Fixation Time(sec)     : %.2f\n',sparam.initial_fixation_time);
fprintf('Cond Duration(sec)     : %.2f\n',sparam.condition_duration);
fprintf('Between Trial Dur(sec) : %.2f\n',sparam.BetweenDuration);
fprintf('Stim ON Duration(sec)  : %.2f\n',sparam.stim_on_duration);
fprintf('Stim OFF Duration(sec) : %.2f\n',sparam.stim_off_duration);
fprintf('*** The number of cond, block, trials, imgs ****\n');
fprintf('#trials in each cond   : %d\n',sparam.numTrials);
fprintf('#conditions            : %d\n',sparam.numConds);
fprintf('Noise mode             : %s\n',sparam.noise_mode);
fprintf('************ Response key settings *************\n');
fprintf('Reponse Key #1         : %d=%s\n',dparam.Key1,KbName(dparam.Key1));
fprintf('Reponse Key #2         : %d=%s\n',dparam.Key2,KbName(dparam.Key2));
fprintf('************************************************\n\n');
fprintf('Please carefully check before proceeding.\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate trial design & response matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create trial design array (stimulus presentation sequences)
if numel(sparam.numConds)>=4
  design=GenerateRandomDesignSequence(sparam.numConds,sparam.numTrials,2,0,1)';
elseif numel(sparam.numConds)<=2
  design=GenerateRandomDesignSequence(sparam.numConds,sparam.numTrials,0,0,0)';
else
  design=GenerateRandomDesignSequence(sparam.numConds,sparam.numTrials,1,0,1)';
end

% create response matrix
% reponse matrix stores the number of trials in which the observer perceive the planes are to far
% [1 x numel(sprama.disparity)] matrix
respmatrix=zeros(1,sparam.numConds);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialize response & event logger objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize MATLAB objects for event and response logs
event=eventlogger();
resps=responselogger([dparam.Key1,dparam.Key2]);
resps.initialize(event); % initialize responselogger


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Wait for user reponse to start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~force_proceed_flag
  [user_answer,resps]=resps.wait_to_proceed();
  if ~user_answer, diary off; return; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialization of Left & Right screens for binocular presenting/viewing mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dparam.skip_sync_test, Screen('Preference','SkipSyncTests',1); end

% open the Screen(s)
[winPtr,winRect,nScr,dparam.fps,dparam.ifi,initDisplay_OK]=InitializePTBDisplays(dparam.ExpMode,sparam.bgcolor,0,[],dparam.scrID);
if ~initDisplay_OK, error('Display initialization error. Please check your exp_run parameter.'); end
HideCursor();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setting the PTB runnning priority to MAX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the priority of this script to MAX
priorityLevel=MaxPriority(winPtr,'WaitBlanking');
Priority(priorityLevel);

% conserve VRAM memory: Workaround for flawed hardware and drivers
% 32 == kPsychDontShareContextRessources: Do not share ressources between
% different onscreen windows. Usually you want PTB to share all ressources
% like offscreen windows, textures and GLSL shaders among all open onscreen
% windows. If that causes trouble for some weird reason, you can prevent
% automatic sharing with this flag.
%Screen('Preference','ConserveVRAM',32);

% Enable OpenGL mode of Psychtoolbox: This is crucially needed for clut animation
InitializeMatlabOpenGL();
AssertOpenGL();

% set OpenGL blend functions
Screen('BlendFunction', winPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying 'Initializing...'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displaying texts on the center of the screen
DisplayMessage2('Initializing...',sparam.bgcolor,winPtr,nScr,'Arial',36);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cm per pix
sparam.cm_per_pix=1/sparam.pix_per_cm;

% pixles per degree
sparam.pix_per_deg=round( 1/( 180*atan(sparam.cm_per_pix/sparam.vdist)/pi ) );

% adjust parameters for oversampling (this adjustment shoud be done after creating heightfields)
dotDens=sparam.dotDens/sparam.oversampling_ratio;
ipd=sparam.ipd*sparam.oversampling_ratio;
vdist=sparam.vdist*sparam.oversampling_ratio;
pix_per_cm_x=sparam.pix_per_cm*sparam.oversampling_ratio;
pix_per_cm_y=sparam.pix_per_cm;

% generate ovals to be used in RDS
dotSize=round(sparam.dotRadius.*[pix_per_cm_y,pix_per_cm_x]*2); % radius(cm) --> diameter(pix)
basedot=double(MakeFineOval(dotSize,[sparam.colors(1:2) 0],sparam.colors(3),1.2,2,1,0,0));
wdot=basedot(:,:,1);     % get only gray scale image (white)
bdot=basedot(:,:,2);     % get only gray scale image (black)
dotalpha=basedot(:,:,4)./max(max(basedot(:,:,4))); % get alpha channel value 0-1.0;

%% prepare the target rectangular stimuli

% generate the rectangular field
rect_field=cell(numel(sparam.disparity),1);
for ii=1:1:numel(sparam.disparity)
  rect_height=[CalcDistFromDisparity(sparam.ipd,sparam.disparity(ii)+sparam.base_disparity,sparam.vdist),...
               CalcDistFromDisparity(sparam.ipd,sparam.base_disparity,sparam.vdist)]; % unit: cm
  rect_field{ii}=nf_CreateRectField(sparam.outerRectFieldSize,sparam.innerRectFieldSize,sparam.gapRectFieldSize,rect_height,sparam.pix_per_deg,sparam.oversampling_ratio); % unit: cm
end

% calculate position shifts for each value of heightfield
pos=cell(numel(sparam.disparity),2); % 2 = left/right
for ii=1:1:numel(sparam.disparity)
  [pos{ii,1},pos{ii,2}]=RayTrace_ScreenPos_X_MEX(rect_field{ii},ipd,vdist,pix_per_cm_x,0);
end

%% prepare the reference stimulus
if ~isnan(sparam.reference_disparity)
  rect_height_ref=[CalcDistFromDisparity(sparam.ipd,sparam.reference_disparity+sparam.base_disparity,sparam.vdist),...
               CalcDistFromDisparity(sparam.ipd,sparam.base_disparity,sparam.vdist)]; % unit: cm
  rect_field_ref=nf_CreateRectField(sparam.outerRectFieldSize,sparam.innerRectFieldSize,sparam.gapRectFieldSize,rect_height_ref,sparam.pix_per_deg,sparam.oversampling_ratio); % unit: cm
  [pos_ref{1},pos_ref{2}]=RayTrace_ScreenPos_X_MEX(rect_field_ref,ipd,vdist,pix_per_cm_x,0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Debug codes
%%%% saving the stimulus images as *.png format files and enter the debug (keyboard) mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%% DEBUG codes start here. The codes below are just to get stimulus images and check the dparam and sparam parameters
% note: debug stimuli have no jitters in binocular disparities
if strfind(upper(subjID),'DEBUG')

  Screen('CloseAll');
  ShowCursor();

  imgL=cell(sparam.numConds,1);
  imgR=cell(sparam.numConds,1);
  posL=cell(sparam.numConds,1);
  posR=cell(sparam.numConds,1);

  % target images
  for ii=1:1:sparam.numConds
    % generate RDS images
    if strcmpi(sparam.noise_mode,'none') % generating an intact RDS stimulus
      [imgL{ii},imgR{ii}]=nf_RDSfastest(pos{ii,1},pos{ii,2},wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,1,1);
    elseif strcmpi(sparam.noise_mode,'anti') % generating an anticorrelated RDS stimulus
      [imgL{ii},imgR{ii}]=nf_RDSfastest_with_acor_noise_MEX(pos{ii,1},pos{ii,2},wdot,bdot,dotalpha,dotDens,sparam.noise_ratio,sparam.colors(3));
    elseif strcmpi(sparam.noise_mode,'snr') % generating a Signal-to-Noise-Ratio(SNR)-modulated RDS stimulus
      % in SNR RDS stimuli, noise dots are assigned in advance of generating the RDS.
      obj_idx=find(rect_field{ii}~=0);
      obj_idx=shuffle(obj_idx);
      obj_idx=obj_idx(1:round(numel(obj_idx)*sparam.noise_ratio/100));
      noise_field=sparam.noise_mean+randn(numel(obj_idx),1).*sparam.noise_sd;
      noise_height=CalcDistFromDisparity(ipd,noise_field,vdist); % unit: cm
      [noiseL,noiseR]=RayTrace_ScreenPos_X_MEX(noise_height,ipd,vdist,pix_per_cm_x,0);

      posL{ii}=pos{ii,1};
      posR{ii}=pos{ii,2};
      if strcmpi(sparam.noise_method,'add') % if you want to add noises on the baseline depths, please use the lines below.
        posL{ii}(obj_idx)=posL{ii}(obj_idx)+noiseL;
        posR{ii}(obj_idx)=posR{ii}(obj_idx)+noiseR;
      elseif strcmpi(sparam.noise_method,'replace') % if you want to replace the baseline depths with generated noises, please use the lines below.
        posL{ii}(obj_idx)=noiseL;
        posR{ii}(obj_idx)=noiseR;
      else
        error('sparam.noise_method should be one of ''add'' or ''replace''.');
      end
      [imgL{ii},imgR{ii}]=nf_RDSfastest(posL{ii},posR{ii},wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,0,1);
    else
      error('sparam.noise_mode should be one of ''none'', ''anti'' or ''snr''. check the stimulus file.');
    end
  end

  % reference image
  if ~isnan(sparam.reference_disparity)
    % generate RDS images
    if strcmpi(sparam.noise_mode,'none') % generating an intact RDS stimulus
      [imgL_ref,imgR_ref]=nf_RDSfastest(pos_ref{1},pos_ref{2},wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,1,1);
    elseif strcmpi(sparam.noise_mode,'anti') % generating an anticorrelated RDS stimulus
      [imgL_ref,imgR_ref]=nf_RDSfastest_with_acor_noise_MEX(pos_ref{1},pos_ref{2},wdot,bdot,dotalpha,dotDens,sparam.noise_ratio,sparam.colors(3));
    elseif strcmpi(sparam.noise_mode,'snr') % generating a Signal-to-Noise-Ratio(SNR)-modulated RDS stimulus
      % in SNR RDS stimuli, noise dots are assigned in advance of generating the RDS.
      obj_idx=find(rect_field_ref~=0);
      obj_idx=shuffle(obj_idx);
      obj_idx=obj_idx(1:round(numel(obj_idx)*sparam.noise_ratio/100));
      noise_field=sparam.noise_mean+randn(numel(obj_idx),1).*sparam.noise_sd;
      noise_height=CalcDistFromDisparity(ipd,noise_field,vdist); % unit: cm
      [noiseL,noiseR]=RayTrace_ScreenPos_X_MEX(noise_height,ipd,vdist,pix_per_cm_x,0);

      posL_ref=pos_ref{1};
      posR_ref=pos_ref{2};
      if strcmpi(sparam.noise_method,'add') % if you want to add noises on the baseline depths, please use the lines below.
        posL_ref(obj_idx)=posL_ref(obj_idx)+noiseL;
        posR_ref(obj_idx)=posR_ref(obj_idx)+noiseR;
      elseif strcmpi(sparam.noise_method,'replace') % if you want to replace the baseline depths with generated noises, please use the lines below.
        posL_ref(obj_idx)=noiseL;
        posR_ref(obj_idx)=noiseR;
      else
        error('sparam.noise_method should be one of ''add'' or ''replace''.');
      end
      [imgL_ref,imgR_ref]=nf_RDSfastest(posL_ref,posR_ref,wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,0,1);
    else
      error('sparam.noise_mode should be one of ''none'', ''anti'' or ''snr''. check the stimulus file.');
    end
  end

  % save stimuli as *.mat
  save_dir=fullfile(resultDir,'images');
  if ~exist(save_dir,'dir'), mkdir(save_dir); end
  save(fullfile(save_dir,'stereo_screening_stimuli.mat'),'design','dparam','sparam','imgL','imgR','posL','posR','wdot','bdot','dotalpha');
  if ~isnan(sparam.reference_disparity)
    save(fullfile(save_dir,'stereo_screening_stimuli.mat'),'-append','imgL_ref','imgR_ref','posL_ref','posR_ref');
  end

  % plotting/saving figures of the target stimuli
  for ii=1:1:sparam.numConds
    % save generated figures as png
    if ~strcmpi(dparam.ExpMode,'redgreen') && ~strcmpi(dparam.ExpMode,'redblue')
      M = [imgL{ii},sparam.bgcolor(3)*ones(size(imgL{ii},1),20),imgR{ii},sparam.bgcolor(3)*ones(size(imgL{ii},1),20),imgL{ii}];
    else
      M=reshape([imgL{ii},imgR{ii},sparam.bgcolor(3)*ones(size(imgL{ii}))],[size(imgL{ii}),3]); % RGB;
    end

    figure; hold on;
    imshow(M,[0,255]); %#ok
    if ~strcmpi(dparam.ExpMode,'redgreen') && ~strcmpi(dparam.ExpMode,'redblue')
      fname=sprintf('stereo_screening_cond_%03d.png',ii);
    else
      fname=sprintf('stereo_screening_red_green_cond_%03d.png',ii);
    end
    imwrite(M,[save_dir,filesep(),fname,'.png'],'png');
  end % for ii=1:1:sparam.numConds

  % plotting/saving figures of the reference stimulus
  if ~isnan(sparam.reference_disparity)
    % save generated figures as png
    if ~strcmpi(dparam.ExpMode,'redgreen') && ~strcmpi(dparam.ExpMode,'redblue')
      M = [imgL_ref,sparam.bgcolor(3)*ones(size(imgL_ref,1),20),imgR_ref,sparam.bgcolor(3)*ones(size(imgL_ref,1),20),imgL_ref];
    else
      M=reshape([imgL_ref,imgR_ref,sparam.bgcolor(3)*ones(size(imgL_ref))],[size(imgL_ref),3]); % RGB;
    end

    figure; hold on;
    imshow(M,[0,255]); %#ok
    if ~strcmpi(dparam.ExpMode,'redgreen') && ~strcmpi(dparam.ExpMode,'redblue')
      fname='stereo_screening_reference.png';
    else
      fname='stereo_screening_red_green_reference.png';
    end
    imwrite(M,[save_dir,filesep(),fname,'.png'],'png');
  end % if ~isnan(sparam.reference_disparity)

  keyboard;

  % clean up
  Priority(0);
  GammaResetPTB(1.0);
  rmpath(genpath(fullfile(rootDir,'..','Common')));
  rmpath(fullfile(rootDir,'..','Generation'));
  rmpath(fullfile(rootDir,'..','mpsignifit'));
  diary off;

  return;

end % if strfind(upper(subjID),'DEBUG')
% %%%%%% DEBUG code ends here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creating background images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the central aperture size of the background image
edgeY=mod(dparam.ScrHeight,sparam.patch_num(1)); % delete exceeded region
p_height=round((dparam.ScrHeight-edgeY)/sparam.patch_num(1)); % height in pix of patch_height + interval-Y

edgeX=mod(dparam.ScrWidth,sparam.patch_num(2)); % delete exceeded region
p_width=round((dparam.ScrWidth-edgeX)/sparam.patch_num(2)); % width in pix of patch_width + interval-X

aperture_size(1)=2*( p_height*ceil(size(pos{1,1},1)/2/p_height) );
aperture_size(2)=2*( p_width*ceil(size(pos{1,1},2)./sparam.oversampling_ratio/2/p_width) );

bgimg = CreateBackgroundImage([dparam.ScrHeight,dparam.ScrWidth],...
          aperture_size,sparam.patch_size,sparam.bgcolor,sparam.patch_color1,sparam.patch_color2,sparam.fixcolor,sparam.patch_num,0,0,0);
background = Screen('MakeTexture',winPtr,bgimg{1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creating the central fixation, cross images (left/right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create fixation cross images
[fix_L,fix_R]=CreateFixationImg(sparam.fixsize,sparam.fixcolor,sparam.bgcolor,sparam.fixlinesize(2),sparam.fixlinesize(1),0,0);
[dark_fix_L,dark_fix_R]=CreateFixationImg(sparam.fixsize,0.2*sparam.fixcolor,sparam.bgcolor,sparam.fixlinesize(2),sparam.fixlinesize(1),0,0);

[correct_fix_L,correct_fix_R]=CreateFixationImg(sparam.fixsize,[0,255,0],sparam.bgcolor,sparam.fixlinesize(2),sparam.fixlinesize(1),0,0);
[incorrect_fix_L,incorrect_fix_R]=CreateFixationImg(sparam.fixsize,[255,0,0],sparam.bgcolor,sparam.fixlinesize(2),sparam.fixlinesize(1),0,0);

fcross=cell(2,2); % the first 2 = bright/dark and the second 2 = left/right
fcross{1,1}=Screen('MakeTexture',winPtr,fix_L);
fcross{1,2}=Screen('MakeTexture',winPtr,fix_R);
fcross{2,1}=Screen('MakeTexture',winPtr,dark_fix_L);
fcross{2,2}=Screen('MakeTexture',winPtr,dark_fix_R);
fcross{3,1}=Screen('MakeTexture',winPtr,correct_fix_L);
fcross{3,2}=Screen('MakeTexture',winPtr,correct_fix_R);
fcross{4,1}=Screen('MakeTexture',winPtr,incorrect_fix_L);
fcross{4,2}=Screen('MakeTexture',winPtr,incorrect_fix_R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare blue lines for stereo image flip sync with VPixx PROPixx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There seems to be a blueline generation bug on some OpenGL systems.
% SetStereoBlueLineSyncParameters(winPtr, winRect(4)) corrects the
% bug on some systems, but breaks on other systems.
% We'll just disable automatic blueline, and manually draw our own bluelines!

if strcmpi(dparam.ExpMode,'propixxstereo')
  SetStereoBlueLineSyncParameters(winPtr, winRect(4)+10);
  blueRectOn(1,:)=[0, winRect(4)-1, winRect(3)/4, winRect(4)];
  blueRectOn(2,:)=[0, winRect(4)-1, winRect(3)*3/4, winRect(4)];
  blueRectOff(1,:)=[winRect(3)/4, winRect(4)-1, winRect(3), winRect(4)];
  blueRectOff(2,:)=[winRect(3)*3/4, winRect(4)-1, winRect(3), winRect(4)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% image size adjusting to match the current display resolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dparam.fullscr
  ratio_wid=( winRect(3)-winRect(1) )/dparam.ScrWidth;
  ratio_hei=( winRect(4)-winRect(2) )/dparam.ScrHeight;
  stimSize =[size(pos{1,1},2)*ratio_wid,size(pos{1,1},1)*ratio_hei].*[1/sparam.oversampling_ratio,1];
  bgSize=[size(bgimg{1},2)*ratio_wid,size(bgimg{1},1)*ratio_hei];
  fixSize=[2*sparam.fixsize*ratio_wid,2*sparam.fixsize*ratio_hei];
else
  stimSize=[size(pos{1,1},2),size(pos{1,1},1)].*[1/sparam.oversampling_ratio,1];
  bgSize=[dparam.ScrWidth,dparam.ScrHeight];
  fixSize=[2*sparam.fixsize,2*sparam.fixsize];
end

% for some display modes in which one screen is splitted into two binocular displays
if strcmpi(dparam.ExpMode,'cross') || strcmpi(dparam.ExpMode,'parallel') || ...
   strcmpi(dparam.ExpMode,'topbottom') || strcmpi(dparam.ExpMode,'bottomtop')
  stimSize=stimSize./2;
  bgSize=bgSize./2;
  fixSize=fixSize./2;
end

stimRect=[0,0,stimSize]; % used to display target stimuli
bgRect=[0,0,bgSize]; % used to display background images;
fixRect=[0,0,fixSize]; % used to display the central fixation point


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Saving the current parameters temporally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saving the current parameters
% (this is required to analyze part of the data obtained even when the experiment is interrupted unexpectedly)
fprintf('saving the stimulus generation and presentation parameters...');
savefname=fullfile(resultDir,[num2str(subjID),'_stereo_screening_results_run_',num2str(acq,'%02d'),'.mat']);

% backup the old file(s)
if ~overwrite_flg
  BackUpObsoleteFiles(fullfile('subjects',num2str(subjID),'results',today),...
                      [num2str(subjID),'_stereo_screening_results_run_',num2str(acq,'%02d'),'.mat'],'_old');
end

% save the current parameters
eval(sprintf('save %s subjID acq design sparam dparam gamma_table;',savefname));
fprintf('done.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying 'Ready to Start'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displaying texts on the center of the screen
DisplayMessage2('Ready to Start',sparam.bgcolor,winPtr,nScr,'Arial',36);
ttime=GetSecs(); while (GetSecs()-ttime < 0.5), end  % run up the clock.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Flip the display(s) to the background image(s)
%%%% and inform the ready of stimulus presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change the screen and wait for the trigger or pressing the start button
for nn=1:1:nScr
  Screen('SelectStereoDrawBuffer',winPtr,nn-1);
  Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
  Screen('DrawTexture',winPtr,fcross{2,nn},[],CenterRect(fixRect,winRect));

  % blue line for stereo sync
  if strcmpi(dparam.ExpMode,'propixxstereo')
    Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
    Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
  end
end
Screen('DrawingFinished',winPtr);
Screen('Flip', winPtr,[],[],[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Wait for the first trigger pulse from fMRI scanner or start with button pressing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add time stamp (this also works to load add_event method in memory in advance of the actual displays)
fprintf('\nWaiting for the start...\n');
event=event.add_event('Experiment Start',strcat([datestr(now,'yymmdd'),' ',datestr(now,'HH:mm:ss')]),NaN);

% waiting for stimulus presentation
resps.wait_stimulus_presentation(dparam.start_method,dparam.custom_trigger);
PlaySound(1);
fprintf('\nExperiment running...\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Event logs and timer (!start here!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[event,the_experiment_start]=event.set_reference_time(GetSecs());
targetTime=the_experiment_start;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial Fixation Period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wait for the initial fixation
if sparam.initial_fixation_time~=0
  event=event.add_event('Initial Fixation',[]);
  fprintf('\nfixation\n\n');

  for nn=1:1:nScr
    Screen('SelectStereoDrawBuffer',winPtr,nn-1);
    Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
    Screen('DrawTexture',winPtr,fcross{1,nn},[],CenterRect(fixRect,winRect));

    % blue line for stereo sync
    if strcmpi(dparam.ExpMode,'propixxstereo')
      Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
      Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
    end
  end
  Screen('DrawingFinished',winPtr);
  Screen('Flip',winPtr,[],[],[]);

  % wait for the initial fixation
  targetTime=targetTime+sparam.initial_fixation_time;
  while (GetSecs() < targetTime), [resps,event]=resps.check_responses(event); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The Trial Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for currenttrial=1:1:numel(design)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Stimulus generation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tBetweenTrial=GetSecs();

  %% generate the target stimulus

  % get the current stimulus parameters
  stimID=design(currenttrial);
  disparity=sparam.disparity(stimID);

  % generate RDS images
  if strcmpi(sparam.noise_mode,'none') % generating an intact RDS stimulus
    [imgL,imgR]=nf_RDSfastest(pos{stimID,1},pos{stimID,2},wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,1,1);
  elseif strcmpi(sparam.noise_mode,'anti') % generating an anticorrelated RDS stimulus
    [imgL,imgR]=nf_RDSfastest_with_acor_noise_MEX(pos{stimID,1},pos{stimID,2},wdot,bdot,dotalpha,dotDens,sparam.noise_ratio,sparam.colors(3));
  elseif strcmpi(sparam.noise_mode,'snr') % generating a Signal-to-Noise-Ratio(SNR)-modulated RDS stimulus
    % in SNR RDS stimuli, noise dots are assigned in advance of generating the RDS.
    obj_idx=find(rect_field{stimID}~=0);
    obj_idx=shuffle(obj_idx);
    obj_idx=obj_idx(1:round(numel(obj_idx)*sparam.noise_ratio/100));
    noise_field=sparam.noise_mean+randn(numel(obj_idx),1).*sparam.noise_sd;
    noise_height=CalcDistFromDisparity(ipd,noise_field,vdist); % unit: cm
    [noiseL,noiseR]=RayTrace_ScreenPos_X_MEX(noise_height,ipd,vdist,pix_per_cm_x,0);

    posL=pos{stimID,1};
    posR=pos{stimID,2};
    if strcmpi(sparam.noise_method,'add') % if you want to add noises on the baseline depths, please use the lines below.
      posL(obj_idx)=posL(obj_idx)+noiseL;
      posR(obj_idx)=posR(obj_idx)+noiseR;
    elseif strcmpi(sparam.noise_method,'replace') % if you want to replace the baseline depths with generated noises, please use the lines below.
      posL(obj_idx)=noiseL;
      posR(obj_idx)=noiseR;
    else
      error('sparam.noise_method should be one of ''add'' or ''replace''.');
    end
    [imgL,imgR]=nf_RDSfastest(posL,posR,wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,0,1);
  else
    error('sparam.noise_mode should be one of ''none'', ''anti'' or ''snr''. check the stimulus file.');
  end
  stim{1}=Screen('MakeTexture',winPtr,imgL);
  stim{2}=Screen('MakeTexture',winPtr,imgR);

  %% generate the reference stimulus

  if ~isnan(sparam.reference_disparity)
    % generate RDS images
    if strcmpi(sparam.noise_mode,'none') % generating an intact RDS stimulus
      [imgL,imgR]=nf_RDSfastest(pos_ref{1},pos_ref{2},wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,1,1);
    elseif strcmpi(sparam.noise_mode,'anti') % generating an anticorrelated RDS stimulus
      [imgL,imgR]=nf_RDSfastest_with_acor_noise_MEX(pos_ref{1},pos_ref{2},wdot,bdot,dotalpha,dotDens,sparam.noise_ratio,sparam.colors(3));
    elseif strcmpi(sparam.noise_mode,'snr') % generating a Signal-to-Noise-Ratio(SNR)-modulated RDS stimulus
      % in SNR RDS stimuli, noise dots are assigned in advance of generating the RDS.
      obj_idx=find(rect_field_ref~=0);
      obj_idx=shuffle(obj_idx);
      obj_idx=obj_idx(1:round(numel(obj_idx)*sparam.noise_ratio/100));
      noise_field=sparam.noise_mean+randn(numel(obj_idx),1).*sparam.noise_sd;
      noise_height=CalcDistFromDisparity(ipd,noise_field,vdist); % unit: cm
      [noiseL,noiseR]=RayTrace_ScreenPos_X_MEX(noise_height,ipd,vdist,pix_per_cm_x,0);

      posL=pos_ref{1};
      posR=pos_ref{2};
      if strcmpi(sparam.noise_method,'add') % if you want to add noises on the baseline depths, please use the lines below.
        posL(obj_idx)=posL(obj_idx)+noiseL;
        posR(obj_idx)=posR(obj_idx)+noiseR;
      elseif strcmpi(sparam.noise_method,'replace') % if you want to replace the baseline depths with generated noises, please use the lines below.
        posL(obj_idx)=noiseL;
        posR(obj_idx)=noiseR;
      else
        error('sparam.noise_method should be one of ''add'' or ''replace''.');
      end
      [imgL,imgR]=nf_RDSfastest(posL,posR,wdot,bdot,dotalpha,dotDens,sparam.colors(3),0,0,1);
    else
      error('sparam.noise_mode should be one of ''none'', ''anti'' or ''snr''. check the stimulus file.');
    end
    stim_ref{1}=Screen('MakeTexture',winPtr,imgL);
    stim_ref{2}=Screen('MakeTexture',winPtr,imgR);
  end % if ~isnan(sparam.reference_disparity)

  % wait for the BetweenDuration
  tBetweenTrial=tBetweenTrial+sparam.BetweenDuration;
  while GetSecs()<tBetweenTrial, resps.check_responses(event); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Stimulus display & animation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf('TRIAL [%03d]: ID %03d, disparity %.2f arcmin, response ',currenttrial,stimID,disparity);

  tStart=GetSecs();
  event=event.add_event('Start Trial',disparity);

  % determine stimulus prsentation order
  if ~isnan(sparam.reference_disparity) % for 2AFC, 1=target, 2=reference
    stim_order=shuffle([1,2]);
  else % if isnan(sparam.reference_disparity) % for 1AFC, 1=target alone
    stim_order=1;
  end

  % presentation
  for pp=stim_order

    %% stimulus ON
    for nn=1:1:nScr
      Screen('SelectStereoDrawBuffer',winPtr,nn-1);
      Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
      if pp==1 % target
        Screen('DrawTexture',winPtr,stim{nn},[],CenterRect(stimRect,winRect));
      else % reference
        Screen('DrawTexture',winPtr,stim_ref{nn},[],CenterRect(stimRect,winRect));
      end
      Screen('DrawTexture',winPtr,fcross{1,nn},[],CenterRect(fixRect,winRect));

      % blue line for stereo sync
      if strcmpi(dparam.ExpMode,'propixxstereo')
        Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
        Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
      end
    end
    Screen('DrawingFinished',winPtr);
    Screen('Flip',winPtr,[],[],[]);
    event=event.add_event('Stimulus ON',disparity);

    % wait for stim_on_duration
    tStart=tStart+sparam.stim_on_duration;
    while GetSecs()<tStart, resps.check_responses(event); end

    %% stimulus OFF
    for nn=1:1:nScr
      Screen('SelectStereoDrawBuffer',winPtr,nn-1);
      Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
      Screen('DrawTexture',winPtr,fcross{1,nn},[],CenterRect(fixRect,winRect));

      % blue line for stereo sync
      if strcmpi(dparam.ExpMode,'propixxstereo')
        Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
        Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
      end
    end
    Screen('DrawingFinished',winPtr);
    Screen('Flip',winPtr,[],[],[]);
    event=event.add_event('Stimulus OFF',disparity);

    % wait for stim_off_duration
    tStart=tStart+sparam.stim_off_duration;
    while GetSecs()<tStart, resps.check_responses(event); end

  end % for pp=stim_order

  %% observer response

  % display response cue
  for nn=1:1:nScr
    Screen('SelectStereoDrawBuffer',winPtr,nn-1);
    Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
    Screen('DrawTexture',winPtr,fcross{2,nn},[],CenterRect(fixRect,winRect));

    % blue line for stereo sync
    if strcmpi(dparam.ExpMode,'propixxstereo')
      Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
      Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
    end
  end
  Screen('DrawingFinished',winPtr);
  Screen('Flip',winPtr,[],[],[]);

  % get a response
  respFlag=0;
  while ~respFlag
    [x,y,button]=GetMouse(); %#ok
    [dummy1,dummy2,keyCode]=resps.check_responses(event); %#ok
    if ~isnan(sparam.reference_disparity)
      disparity_comp=disparity-sparam.reference_disparity;
      if stim_order(1)==2 % reference is presented first
        if button(1) || keyCode(dparam.Key1)==1 % reference is nearer
          respmatrix(stimID)=respmatrix(stimID)+1;
          event=event.add_event('Response','First');
          respFlag=sign(disparity_comp); % if disparity_comp>0 (reference is nearer): respFlag=1 (correct), while if disparity_comp<0: respFlag=-1 (incorrect)
        elseif button(3) || keyCode(dparam.Key2)==1
          event=event.add_event('Response','Second');
          respFlag=-1*sign(disparity_comp); % if disparity_comp>0 (reference is farther): respFlag=1 (correct), while if disparity_comp<0 (near): respFlag=-1 (incorrect)
        else
          respFlag=0;
        end
      else % if stim_order(1)==1 % target is presented first
        if button(1) || keyCode(dparam.Key1)==1 % reference is farther
          event=event.add_event('Response','First');
          respFlag=-1*sign(disparity_comp); % if disparity_comp<0 (target is nearer): respFlag=1 (correct), while if disparity_comp>0: respFlag=-1 (incorrect)
        elseif button(3) || keyCode(dparam.Key2)==1
          respmatrix(stimID)=respmatrix(stimID)+1;
          event=event.add_event('Response','Second');
          respFlag=sign(disparity_comp); % if disparity_comp>0 (target is farther): respFlag=1 (correct), while if disparity_comp<0: respFlag=-1 (incorrect)
        else
          respFlag=0;
        end
      end
    else % if isnan(sparam.reference_disparity)
      if button(1) || keyCode(dparam.Key1)==1
        event=event.add_event('Response','Near');
        respFlag=-1*sign(disparity); % if disparity<0 (near): respFlag=1 (correct), while if disparity>0 (far): respFlag=-1 (incorrect)
      elseif button(3) || keyCode(dparam.Key2)==1
        respmatrix(stimID)=respmatrix(stimID)+1;
        event=event.add_event('Response','Far');
        respFlag=sign(disparity); % if disparity>0 (far): respFlag=1 (correct), while if disparity<0 (near): respFlag=-1 (incorrect)
      else
        respFlag=0;
      end
    end
  end

  if respFlag==1 % correct response
    fprintf('correct\n');
  else % if respFlag==-1 or 0 % incorrect response
    fprintf('incorrect\n');
  end

  % giving a correct/incorrect feedback
  if dparam.givefeedback
    if respFlag==1 % correct response
      event=event.add_event('Feedback','Correct');
    else % if respFlag==-1 or 0 % incorrect response
      event=event.add_event('Feedback','Incorrect');
    end
    for nn=1:1:nScr
      Screen('SelectStereoDrawBuffer',winPtr,nn-1);
      Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
      if respFlag==1 % correct response
        Screen('DrawTexture',winPtr,fcross{3,nn},[],CenterRect(fixRect,winRect));
      else % if respFlag==-1 or 0 % incorrect response
        Screen('DrawTexture',winPtr,fcross{4,nn},[],CenterRect(fixRect,winRect));
      end

      % blue line for stereo sync
      if strcmpi(dparam.ExpMode,'propixxstereo')
        Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
        Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
      end
    end
    Screen('DrawingFinished',winPtr);
    Screen('Flip',winPtr,[],[],[]);
    PlaySound(respFlag>0);

    ctime=GetSecs();
    while GetSecs()<ctime+0.5, resps.check_responses(event); end
  end

  %% back to the default view and wait for sparam.BetweenDuration (duration between trials)

  for nn=1:1:nScr
    Screen('SelectStereoDrawBuffer',winPtr,nn-1);
    Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect));
    Screen('DrawTexture',winPtr,fcross{1,nn},[],CenterRect(fixRect,winRect));

    % blue line for stereo sync
    if strcmpi(dparam.ExpMode,'propixxstereo')
      Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
      Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
    end
  end
  Screen('DrawingFinished',winPtr);
  Screen('Flip',winPtr,[],[],[]);

  % garbage collections, clean up the current texture & release memory
  for ii=1:1:2, Screen('Close',stim{ii}); end
  if ~isnan(sparam.reference_disparity), for ii=1:1:2, Screen('Close',stim_ref{ii}); end, end
  event=event.add_event('End Trial',disparity);

end % for currenttrial=1:1:numel(design)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment & scanner end here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experimentDuration=GetSecs()-the_experiment_start;
event=event.add_event('End',[]);

fprintf(['\nExperiment Duration was: ', num2str(experimentDuration),' secs\n\n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Write data into file for post analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saving the results
fprintf('saving data...');

% save data
savefname=fullfile(resultDir,[num2str(subjID),'_stereo_screening_results_run_',num2str(acq,'%02d'),'.mat']);
eval(sprintf('save -append %s subjID acq sparam dparam design event gamma_table respmatrix;',savefname));
fprintf('done.\n');

% tell the experimenter that the measurements are completed
try
  for ii=1:1:3, Snd('Play',sin(2*pi*0.2*(0:900)),8000); end
catch
  % do nothing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot the disparity discrimination performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displaying texts on the center of the screen
DisplayMessage2('Computing the performance...',sparam.bgcolor,winPtr,nScr,'Arial',36);

% sort the data
[disparities,idx]=sort(sparam.disparity);
data=respmatrix(idx);
data=[disparities',data',repmat(sparam.numTrials,[numel(disparities),1])];

% plotting
f1=figure('Name','Near-Far Plane: Disparity Discrimination Performance','Position','default',...
          'Color',[0.95,0.95,0.95],'NumberTitle','off','ToolBar','none','MenuBar','none');
plot(disparities,data(:,2)./sparam.numTrials,'o','MarkerFaceColor',[0,1,1],'MarkerEdgeColor',[0,0,1],'MarkerSize',10,'LineWidth',2);
hold on;

% fit the psychometric function (gaussian) using psignifit
priors.m_or_a = 'None';
priors.w_or_b = 'None';
priors.lambda = 'Uniform(0,.1)';
priors.gamma  = 'Uniform(0,.1)';
inference=BootstrapInference(data,priors,'nafc',1,'sigmoid','gauss','cuts',[0.25, 0.5, 0.75]);
diag=Diagnostics(inference.data,inference.params_estimate,'sigmoid',inference.sigmoid,'core',inference.core,'nafc',inference.nafc,'cuts',inference.cuts);
plot(diag.pmf(:,1),diag.pmf(:,2),'-','color',[0,0,0.5],'LineWidth', 3 );

% get confidence intervals
for cut=1:1:length(inference.cuts)
  ci=getCI(inference,cut,0.95 );
  if inference.nafc>1
    guess=1./inference.nafc;
  else
    guess=inference.params_estimate(end);
  end
  h=guess+inference.cuts(cut)*(1-guess-inference.params_estimate(3));
  plot(ci,[h,h],'-','color',[0,0,0.5]);
end

% get slope and threshold
slope=getSlope(inference,2);
threshold=getThres(inference,2);

plot([min(disparities)-0.5,max(disparities)+0.5],[0.5,0.5],'b:','LineWidth',1);
plot([threshold,threshold],[0,1],'k:','LineWidth',2);
if ~isnan(sparam.reference_disparity) % plot reference disparity
  plot([sparam.reference_disparity,sparam.reference_disparity],[0,1],'b:','LineWidth',3);
end

% decorate figure
set(gca,'Color',[0.95,0.95,1.0]);
set(gca,'XLim',[min(disparities)-0.5,max(disparities)+0.5]);
set(gca,'XTick',disparities);
set(gca,'XTickLabel',disparities);
set(gca,'YLim',[0-0.02,1+0.02]);
set(gca,'YTick',0:0.2:1);

if ~isnan(sparam.reference_disparity)
  title(sprintf('%s: Disparity discrimination performance (2AFC), Thr: %.2f, slope: %.2f (arcmin)',num2str(subjID),threshold,slope));
else
  title(sprintf('%s: Disparity discrimination performance, Thr: %.2f, slope: %.2f (arcmin)',num2str(subjID),threshold,slope));
end
xlabel('Disparity [arcmin (nagative is near)]');
ylabel('Far-selection ratio [0.0 -1.0]');

set(gcf,'PaperPositionMode','auto');
saveas(f1,fullfile(resultDir,[num2str(subjID),'_near_far_rectangle_results_run_',num2str(acq,'%02d'),'.fig']),'fig');
print(f1,fullfile(resultDir,[num2str(subjID),'_near_far_rectangle_results_run_',num2str(acq,'%02d'),'.png']),'-dpng','-r0');

% displaying texts on the center of the screen
DisplayMessage2('Completed.',sparam.bgcolor,winPtr,nScr,'Arial',36);
pause(0.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Cleaning up the PTB screen, removing path to the subfunctions, and finalizing the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Screen('CloseAll');

% closing datapixx
if strcmpi(dparam.ExpMode,'propixxmono') || strcmpi(dparam.ExpMode,'propixxstereo')
  if Datapixx('IsViewpixx3D')
    Datapixx('DisableVideoLcd3D60Hz');
    Datapixx('RegWr');
  end
  Datapixx('Close');
end

ShowCursor();
Priority(0);
GammaResetPTB(1.0);
rmpath(genpath(fullfile(rootDir,'..','Common')));
rmpath(fullfile(rootDir,'..','Generation'));
rmpath(fullfile(rootDir,'..','mpsignifit'));
%setenv('PATH',strrep(getenv('PATH'),[';',fullfile(rootDir,'..','mpsignifit')],''));
clear all; clear mex; clear global;
diary off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Catch the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catch %#ok
  % this "catch" section executes in case of an error in the "try" section
  % above.  Importantly, it closes the onscreen window if its open.
  Screen('CloseAll');

  if exist('dparam','var')
    if isstructmember(dparam,'ExpMode')
      if strcmpi(dparam.ExpMode,'propixxmono') || strcmpi(dparam.ExpMode,'propixxstereo')
        if Datapixx('IsViewpixx3D')
          Datapixx('DisableVideoLcd3D60Hz');
          Datapixx('RegWr');
        end
        Datapixx('Close');
      end
    end
  end

  ShowCursor();
  Priority(0);
  GammaResetPTB(1.0);
  tmp=lasterror; %#ok
  if exist('event','var'), event=event.get_event(); end %#ok % just for debugging
  diary off;
  fprintf(['\nError detected and the program was terminated.\n',...
           'To check error(s), please type ''tmp''.\n',...
           'Please save the current variables now if you need.\n',...
           'Then, quit by ''dbquit''\n']);
  keyboard;
  rmpath(genpath(fullfile(rootDir,'..','Common')));
  rmpath(fullfile(rootDir,'..','Generation'));
  rmpath(fullfile(rootDir,'..','mpsignifit'));
  clear all; clear mex; clear global;
  return
end % try..catch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% finish the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return;
% function StereoScreening
