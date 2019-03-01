function OK=run_exp(subj,acq_ids)

% function run_exp(subj,acq_ids)
% (: is optional)
%
% This is a wrapper to run StereoScreening function.
%
% It displays rectangular planes with binocular disparity (+/- arcmings)
% for testing psychophysical disparity accuracy.
% This script shoud be runned with MATLAB Psychtoolbox version 3 or above.
%
% [input]
% subj    : subject's name, e.g. 'HB'
%           if a non-existed subj name is set, the _DEFAULT_ parameters
%           will be copied and the measurement will start with those parameters
% acq_ids : acquisition number (run), 1,2,3,...
%
% [output]
% OK      : whether this script is completed without any error [true/false]
%
% Created    : "2010-06-25 08:37:58 ban"
% Last Update: "2019-03-01 15:15:19 ban"

% constants, you can change these for your own purpose.
stim_fname='nf_stimulus';
disp_fname='nf_display';

% set the program body
run_script='StereoScreening';


% ********************************************************************************************************
%% check the input variables
% ********************************************************************************************************

if nargin<1, help(mfilename()); OK=false; return; end
if nargin<2 || isempty(acq_ids)
  warning('MATLAB:acq_ids_warning','warning: acq_ids is not specified. using 1...\n');
  acq_ids=1;
end


% ********************************************************************************************************
%% check directory with subject name
% ********************************************************************************************************

% [NOTE]
% if the subj directory is not found, create subj directory, copy all condition
% files from DEFAULT and then run the script using DEFAULT parameters

subj_dir=fullfile(fileparts(mfilename('fullpath')),'subjects',subj);
if ~exist(subj_dir,'dir')

  disp('The subject directory was not found.');
  user_response=0;
  while ~user_response
    user_entry = input('Do you want to proceed using DEFAULT parameters? (y/n) : ', 's');
    if(user_entry == 'y')
      fprintf('Generating subj directory using DEFAULT parameters...');
      user_response=1; %#ok
      break;
    elseif (user_entry == 'n')
      disp('quiting the script...');
      if nargout, OK=false; end
      return;
    else
      disp('Please answer y or n!'); continue;
    end
  end

  %mkdir(subj_dir);
  copyfile(fullfile(fileparts(mfilename('fullpath')),'subjects','_DEFAULT_'),subj_dir);
end


% ********************************************************************************************************
% *** set gamma table. please change the line below to use the actual measuments of the display gamma. ***
% ********************************************************************************************************

% loading gamma_table
load(fullfile('..','gamma_table','ASUS_ROG_Swift_PG278Q','181003','cbs','gammatablePTB.mat'));
%load(fullfile('..','gamma_table','ASUS_VG278HE','181003','cbs','gammatablePTB.mat'));
%load(fullfile('..','gamma_table','MEG_B1','151225','cbs','gammatablePTB.mat'));
%gammatable=repmat(linspace(0.0,1.0,256),3,1)'; %#ok % a simple linear gamma


% ********************************************************************************************************
%% run stimulus presentations
% ********************************************************************************************************

% run Experiment
for ii=acq_ids
  main_exp_name=sprintf('%s(''%s'',%d,''%s'',''%s'',gammatable);',run_script,subj,ii,disp_fname,stim_fname);
  eval(main_exp_name);
end
OK=true;

return;
