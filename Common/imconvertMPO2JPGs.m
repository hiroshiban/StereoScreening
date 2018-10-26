function [successfiles,errorfiles]=imconvertMPO2JPGs(MPO_dir,prefix_mpo,delete_error_flg)

% Converts MPO (multiple picture object) files into separated JPEG files
% function [successfiles,errorfiles]=imconvertMPO2JPGs(MPO_dir,:prefix_mpo,:delete_error_flg)
% (: is optional)
%
% This function reads MPO files in the target directory and converts them
% into the separated JPEG files. The generated JPEG will be named with a
% prefix *_%02d.jpg
%
% [input]
% MPO_dir    : the target directory that contains MPO files.
%              the directory should be set with a relative path format
%              so that the origin is the location where this function
%              is called.
% prefix_mpo : (optional) a prefix to specify MPO files to be processed
%              e.g. prefix_mpo='img_'; empty by default.
% delete_error_flg : (optional) whether delete MPO files which can not
%                    be processed correctly. 0 by default.
%
% [output]
% successfiles : the file name(s) that is(are) successfully processed
% errorfiles   : the file name(s) that is(are) not processed correctly
% the generated JPEG files will be save with a prefix *_%02d.jpg
%
% [dependency]
% wildcardsearch in hbtools
%
%
% Created    : "2015-07-08 14:00:57 ban"
% Last Update: "2015-07-08 16:57:34 ban"

% check the input variable
if nargin<1 || isempty(MPO_dir), help(mfilename()); return; end
if nargin<2 || isempty(prefix_mpo), prefix_mpo=''; end
if nargin<3 || isempty(delete_error_flg), delete_error_flg=0; end

MPO_dir=fullfile(pwd,MPO_dir);
if ~exist(MPO_dir,'dir'), error('can not find MPO_dir. check the input variable.'); end

save_dir=fullfile(pwd,'JPEGs');
if ~exist(save_dir,'dir'), mkdir(save_dir); end

% processing
fprintf('Target : %s\n',MPO_dir);
mpofiles=wildcardsearch(MPO_dir,strrep([prefix_mpo,'*.mpo'],'**','*'));
successcounter=0; successfiles={};
errorcounter=0; errorfiles={};
for ii=1:1:length(mpofiles)
  [~,mpofname,ext]=fileparts(mpofiles{ii});
  fprintf('processing #%05d: %s%s, ',ii,mpofname,ext);
  tgt=relativepath(mpofiles{ii});
  try
    [~,imgsize]=imreadmpo(tgt(1:end-1),1); % this is a file so need to omit irrelevant \ in the tgt
    try
      % check whether the generated jpeg files are valid
      if exist([mpofname,'_01.jpg'],'file') && exist([mpofname,'_02.jpg'],'file')
        imread([mpofname,'_01.jpg']);
        imread([mpofname,'_02.jpg']);

        movefile([mpofname,'_01.jpg'],save_dir);
        movefile([mpofname,'_02.jpg'],save_dir);

        fprintf('[x,y]=[% 4d,% 4d]\n',imgsize(2),imgsize(1));
        clear imgs imgsize;
        successcounter=successcounter+1;
        successfiles{successcounter}=mpofiles{ii}; %#ok
      else
        delete([mpofname,'_01.jpg']);
        delete([mpofname,'_02.jpg']);
        errorcounter=errorcounter+1;
        errorfiles{errorcounter}=mpofiles{ii}; %#ok
        fprintf('!ERROR! can not read the MPO file.\n');
      end
    catch
      delete([mpofname,'_01.jpg']);
      delete([mpofname,'_02.jpg']);
      errorcounter=errorcounter+1;
      errorfiles{errorcounter}=mpofiles{ii}; %#ok
      fprintf('!ERROR! can not read the MPO file.\n');
    end
  catch
    delete([mpofname,'_01.jpg']);
    delete([mpofname,'_02.jpg']);
    errorcounter=errorcounter+1;
    errorfiles{errorcounter}=mpofiles{ii}; %#ok
    fprintf('!ERROR! can not read the MPO file.\n');
  end
end
disp('completed');

if delete_error_flg
  for ii=1:1:length(errorfiles), delete(errorfiles{ii}); end
end

fprintf('\nOverview\n');
fprintf('success/all=%d/%d, error/all=%d/%d\n',...
        length(successfiles),length(mpofiles),length(errorfiles),length(mpofiles));

return
