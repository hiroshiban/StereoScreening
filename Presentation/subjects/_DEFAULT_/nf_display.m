% ************************************************************
% This_is_the_display_file_for_NEAR_FAR_experiment.
% Please_change_the_parameters_below.
% MotionParallaxConstant.m
% Programmed_by_Hiroshi_Ban___Aug_24_2015
% ************************************************************

% display mode, one of "mono", "dual", "dualparallel", "dualcross", "cross", "parallel", "redgreen", "greenred",
% "redblue", "bluered", "shutter", "topbottom", "bottomtop", "interleavedline", "interleavedcolumn"
dparam.ExpMode='shutter';%'cross';

dparam.scrID=1; % screen ID, generally 0 for a single display setup, 1 for dual display setup

% a method to start stimulus presentation
% 0:ENTER/SPACE, 1:Left-mouse button, 2:the first MR trigger pulse (CiNet),
% 3:waiting for a MR trigger pulse (BUIC) -- checking onset of pin #11 of the parallel port,
% or 4:custom key trigger (wait for a key input that you specify as tgt_key).
dparam.start_method=1;

% a pseudo trigger key from the MR scanner when it starts, only valid when dparam.start_method=4;
dparam.custom_trigger=KbName(32);

%% keyboard settings
dparam.Key1=37; % key 1 = near
dparam.Key2=39; % key 2 = far

%% whether giving correct/incorrect feedback
% 0 = no feedback, 1 = giving correct (green fixation and high-tone sound) or incorrect(red fixation and low-tone sound) feedbacks
dparam.givefeedback=1;

%% screen settings

%%% whether displaying the stimuli in full-screen mode or as is (the precise resolution), 'true' or 'false' (true)
dparam.fullscr=false;

%%% the resolution of the screen height
dparam.ScrHeight=1440;

%% the resolution of the screen width
dparam.ScrWidth=2560;

% whther skipping the PTB's vertical-sync signal test. if 1, the sync test is skipped
dparam.skip_sync_test=0;
