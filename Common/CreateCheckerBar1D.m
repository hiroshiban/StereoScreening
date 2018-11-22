function [checkerboard,mask]=CreateCheckerBar1D(fieldSize,width,angles,steps,pix_per_deg,ndivsL,ndivsS,phase)

% Generates bar-shaped checkerboard patterns masked by a circular aperture (can be used for pRF stimuli) with an individual ID number on each patch.
% function [checkerboard,mask]=CreateCheckerBar1D(:fieldSize,:width,:angles,:steps,:pix_per_deg,:ndivsL,:ndivsS,:phase)
% (: is optional)
%
% This function generates bar-shaped checkerboards with an individual ID number on each patch, which sweeps the visual field.
% The bar stimulus starts at the right horizontal meridian and sweeps the visual field leftwards.
% Each of two checkers have the compensating values of its counterpart.
% Multiple angles can be accepted.
%
% [input]
% fieldSize   : (optional) field size of the stimulus canvas in degree. 12 by default.
%               the bar stimulus is masked by the circular (or ellipse) aperture whose size is fieldSize.
% width       : (optional) bar width in degree. 3 by default.
% angles      : (optional) bar angle in degree, 0 = right horizontal meridian, counter-clockwise. 0 by default.
% steps       : (optional) steps of the bar in sweeping the visual field defined as fieldSize. 16 by default.
% pix_per_deg : (optional) pixels per degree, [val]. 40 by default.
% ndivsL      : (optional) the number of divisions of the bar along long axis (row). 12 by default.
% ndivsS      : (optional) the number of divisions of the bar along short axis (col). 3 by default.
% phase       : (optional) checker's phase along the short axis. 0 by default.
%
% [output]
% checkerboard : output grayscale bar-shaped checkerboard at each of the visual field locations,
%                cell structure, {numel(angles),numel(steps),2}
%                each pixel shows each checker patch's ID or background(0)
% mask        : (optional) checkerboard regional mask, cell structure, logical
%
% Created    : "2018-11-20 11:23:31 ban"
% Last Update: "2018-11-21 19:00:34 ban"

%% check the input variables.
if nargin<1 || isempty(fieldSize), fieldSize=12; end
if nargin<2 || isempty(width), width=3; end
if nargin<3 || isempty(angles), angles=0; end
if nargin<4 || isempty(steps), steps=16; end
if nargin<5 || isempty(pix_per_deg), pix_per_deg=40; end
if nargin<6 || isempty(ndivsL), ndivsL=12; end
if nargin<7 || isempty(ndivsS), ndivsS=3; end
if nargin<8 || isempty(phase), phase=0; end

%% unit conversion

% from degrees to pixels
fieldSize=ceil(fieldSize*pix_per_deg);
width=ceil(width*pix_per_deg);

% from degree to radius
angles=angles.*pi./180;
phase=phase*pi/180;
if phase>2*pi, phase=mod(phase,2*pi); end

% add small lim in checkerboard image, this is to avoid unwanted juggy edges
imsize_ratio=1.01;

%% generate bar checker board

% base xy distance field
[xx,yy]=meshgrid((0:1:imsize_ratio*fieldSize)-imsize_ratio*fieldSize/2,(0:1:imsize_ratio*fieldSize)-imsize_ratio*fieldSize/2);
%if mod(size(xx,1),2), xx=xx(1:end-1,:); yy=yy(1:end-1,:); end
%if mod(size(xx,2),2), xx=xx(:,1:end-1); yy=yy(:,1:end-1); end

% compute distance
r=sqrt(xx.^2+yy.^2);

checkerboard=cell(numel(angles),steps);
mask=cell(numel(angles),steps);

for aa=1:1:numel(angles)

  % !!!NOTICE!!!
  % We need to create each checkerboard following the procedures below
  %  1. generate radian angle field
  %  2. rotate it based on angles & phase
  %  3. generate checkerboard IDs
  % This consumes much CPU power and time, but it is definitely important.
  %
  % To use imrotate after creating one image may look more sophisticated, but we
  % should not do that. This is because when we use imrotate (or fast_rotate)
  % or Screen('DrawTexture',....,rotangle,...), the displayed image will result
  % in low quality with juggy artefact along checker edges.

  done_flag=0;

  % just flip dimension and copy, if the currect checkerboard is one of
  % 180 deg flipped version of previously generated checkerboards.
  % this is to save calculation time
  if aa>=2
    for tt=1:1:aa-1
      %if angles(aa)==mod(angles(tt)+pi,2*pi)
      if abs(angles(aa)-mod(angles(tt)+pi,2*pi))<0.01 % this is to avoid round off error
        %fprintf('#%d checkerboard is generated by just copying/flipping from #%d checkerboard\n',aa,tt); % debug code
        for pp=1:1:size(checkerboard,2)
          checkerboard{aa,pp}=flipdim(flipdim(checkerboard{tt,pp},2),1);
        end
        if nargout>=2
          for pp=1:1:size(checkerboard,2)
            mask{aa,pp}=flipdim(flipdim(mask{tt,pp},2),1);
          end
        end
        done_flag=1;
        break;
      end
    end
  end

  if ~done_flag

    pos=linspace(fieldSize/2,-fieldSize/2,steps);
    xxxx=xx*cos(angles(aa))-yy*sin(angles(aa));
    yyyy=xx*cos(angles(aa)+pi/2)-yy*sin(angles(aa)+pi/2);

    for pp=1:1:numel(pos)
      % get the target bar region
      inidx=find( (-fieldSize/2<=yyyy & yyyy<=fieldSize/2) & (pos(pp)-width/2<=xxxx & xxxx<=pos(pp)+width/2) );

      % assign indices on the patches along the long axis of the bar
      cidl=zeros(size(yyyy)); % checker id along the long axis
      cidl(inidx)=ceil(yyyy(inidx)/fieldSize*ndivsL);
      cidl(inidx)=cidl(inidx)-min(unique(cidl(inidx)))+1; % -min()+1 is required sicne yyyy contains negative (with 0 at the center of the matrix)

      % assign indices on the patches along the short axis of the bar
      if phase~=0
        widi=linspace(pos(pp)-width/2-(width/ndivsS*phase/(2*pi)),pos(pp)+width/2+(width/ndivsS-(width/ndivsS*phase/(2*pi))),ndivsS+2); widi(1)=[]; % checker borders, considering phase shift
      else
        widi=linspace(pos(pp)-width/2,pos(pp)+width/2,ndivsS+1); widi(1)=[];
      end
      cids=zeros(size(xxxx)); % checker id along the short axis
      for kk=length(widi):-1:1
        cids(pos(pp)-width/2<xxxx & xxxx<=min(widi(kk),pos(pp)+width/2))=kk;
      end

      % generate checker's ID
      checkerboard{aa,pp}=zeros(size(xxxx));
      checkerboard{aa,pp}(inidx)=cidl(inidx)+(cids(inidx)-1)*ndivsL;

      % mask the outer retions and delete outliers
      checkerboard{aa,pp}(r>fieldSize/2 | checkerboard{aa,pp}<0)=0;

      % generate mask
      if nargout>=2, mask{aa,pp}=logical(checkerboard{aa,pp}); end
    end % for pp=1:1:size(pos,1)

  end % if ~done_flag

end % for aa=1:1:numel(angles)

return
