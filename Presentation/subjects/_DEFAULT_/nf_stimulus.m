% ************************************************************
% This_is_the_stimulus_parameter_file_for_NEAR_FAR_experiment.
% Please_change_the_parameters_below.
% NearFarRectanglefMRI.m
% Programmed_by_Hiroshi_Ban___August_01_2010
%************************************************************

% "sparam" means "stimulus generation parameters"

%%% image size
sparam.outerRectFieldSize=[12,12]; % target stimulus size in deg
sparam.innerRectFieldSize=[8,8]; % target stimulus size in deg
%sparam.disparity=[8,  4,  2,  1, 0.5, -0.5, -1, -2, -4, -8]; % target disparities in arcmin
sparam.disparity=[-8, -4,  -2,  -1, -0.5, 0.5, 1, 2, 4, 8]; % target disparities in arcmin

%%% the number of trials for each stimulus condition
sparam.numTrials=20;

%%% RDS parameters
sparam.dotRadius=[0.05,0.05]; % radius of RDS's white/black ovals
sparam.dotDens=2; % deinsity of dot in RDS image (1-100)
sparam.colors=[255,0,128]; % RDS colors [dot1,dot2,background](0-255)
sparam.oversampling_ratio=4; % oversampling_ratio for fine scale RDS images, [val]

%%% RDS noise paramters
% when sparam.noise_mode='none' (no noise condition), all parameters related to noise are ignored.
% when sparam.noise_mode='anti' (anticorrelated (left/right dot colors are flipped) RDSs are presented), only sparam.noise_ratio is used in generating RDSs.
% when sparam.noise_mode='snr'  (noises are added in disparities), please set all the noise-related paramters carefully.
sparam.noise_mode='none';  % one of 'none', 'anti', and 'snr'
sparam.noise_ratio=30;     % percentage (0-100) of noises in the RDS, used when sparam.noise_mode='anti' or sparam.noise_mode='snr'
sparam.noise_mean=0;       % noise mean in disparity (arcmin), used only when sparam.noise_mode='snr'
sparam.noise_sd=10;        % noise SD in disparity (arcmin), used only when sparam.noise_mode='snr'
sparam.noise_method='add'; % one of 'add' and 'replace', used only when sparam.noise_mode='snr'

% stimulus display durations in msec

% initial fixation duration in msec
sparam.initial_fixation_time=1000;

%%% duration in msec for each condition
sparam.condition_duration=600;

%%% duration in msec for simulus ON period for each trial, integer (1)
sparam.stim_on_duration=300;

%%% duration in msec for each trial
sparam.BetweenDuration=200;

%%% fixation size & color
sparam.fixsize=24; % radius
sparam.fixlinesize=[12,2]; % [height,width] of the fixation line in pixel
sparam.fixcolor=[255,255,255];

%%% background color
sparam.bgcolor=[128,128,128];

%%% RGB for background patches
sparam.patch_size=[30,30];   % background patch size, [height,width] in pixels
sparam.patch_num=[20,20];    % the number of background patches along vertical and horizontal axis
sparam.color1=[255,255,255]; % 1x3 matrices
sparam.color2=[0,0,0];       % 1x3 matrices

%%% for creating disparity & shadow
run(fullfile(fileparts(mfilename('fullpath')),'sizeparams'));
