function y = clumpvec(x,itype,mode);

%  y = clumpvec(x,itype,mode);
%  
% A function for taking a nearly maximal efficiency m-sequence
% and clumping it randomly to increase the detection power. 
% expects a vector with numbers 0:Q where Q is the number of
% trial types.
%
%  
% Inputs:
%   x:      input vector  
%   itype:  trial type to clump (e.g. a number from 1 to Q) 
%   mode:   clumping mode
%          supported modes:
%
%          1. simple clumping algorithm, fills in 1 hole for specified event type
%
%          2. modified clumping algorithm, when it finds a "hole" it fills it in 
%             completely. 
%
% Outputs:
%    y: a clumped vector  
%
%  revisions
%  030521  -- make changes to fix bug where no filler was found.
%  040323  -- add  comments
%
%  Send comments and questions to ttliu@ucsd.edu
  

  
if ~exist('mode'); mode = 1;end
  

y = x;
x = x(:)';  
npts = length(x);
indi = find(x == itype);
indilen = length(indi);
dind = diff(indi);
  
looping = 1;
tdiff = 2;  % difference between consecutive indices must be >= 2 
            % for this to be considered a hole

noholes = find(dind > 1);
if(isempty(noholes))
  disp('no holes to fill');
  return
end


%find a hole to fill    
while looping  
% start off with the smallest holes
% and work our way up until we find one. 
  
  iok = find(dind == tdiff); 
  if ~isempty(iok)
    looping = 0;
    % assigns a hole at random
    hole = indi(iok(unidrnd(length(iok)))) + 1;
  end
  tdiff = tdiff + 1;
end  

holesize = tdiff-2;
if mode == 1
  holes = hole;
elseif mode == 2
  holes = hole + [0:(holesize-1)];
  holeregion = holes;
  if(holesize > 1)
    if hole-1 > 0
      holeregion = [hole-1 holes];
   end
   if (hole + holesize) <=  npts
     holeregion = [holeregion hole+holesize];
   end
  end
end

% mindiff is the minimum distance of each event to its
% nearest neighbor
  mindiff = min([dind Inf;Inf dind]);
  
maxmind = max(mindiff);
if maxmind == 1
  % there are no singleton blocks so we want to
  % find the smallest block that is furthest away
  % from the nearest block.
  % here is a semi-brute-force algorithm. 

  blockinds = [];blocksizes = [];
  dindlen = length(dind);
  %blocks must start at indices where dind == 1
  iblocklocs = find(dind == 1);
  nlocs = length(iblocklocs);
  blocksearch = 1;
  ipos = 1;imax = 1;
  while blocksearch
    %figure out how long the block is
    blockstart = indi(iblocklocs(ipos));
%    if imax == 11;keyboard;end
    if ipos == nlocs  % a length 2 block at the very end
      blockend = blockstart + 1;
    elseif dind(iblocklocs(ipos)+1) > 1  % length 2 block
      blockend = blockstart + 1;
    else % search to the end of the block
      nextjump = find(dind(iblocklocs(ipos):dindlen) > 1);
      if isempty(nextjump)
	blockend = indi(indilen);
      else
	blockend = blockstart + min(nextjump)-1;
      end	
    end
    
    if blockend == indi(indilen);
      % must be the last block
      blocksearch = 0;
    end
    blocklen = blockend-blockstart+1;
    ipos = ipos + blocklen-1;
    blocksizes = [blocksizes;blocklen];
    blockinds = [blockinds;blockstart];
    imax  = imax + 1;
%    if(imax > 20) blocksearch = 0;end
  end
  numblocks = length(blocksizes);
  blockends = blockinds + blocksizes-1;
  blockgaps = blockinds(2:numblocks)-blockends(1:(numblocks-1));
  %minimum distance to nearest block
  mingaps = min([Inf blockgaps(:)';blockgaps(:)' Inf]);
  minblocksize = min(blocksizes);
  imin = find(blocksizes == minblocksize);
  [biggestgap,ibiggestgap] =  max(mingaps(imin));
  blockind = imin(ibiggestgap);
  fillerpos = blockinds(blockind)+(0:(blocksizes(blockind)-1));
  do_singleton = 0;
else  
  do_singleton = 1;
end
%keyboard;

for k = 1:length(holes)
  thishole = holes(k);
  % find an event to fill the hole
  % grab an event at random which has
  % the farthest nearest neighbor
  % distance to the next nearest neighbor
  if 1
    if do_singleton
      inds = find(mindiff == maxmind);
    else
            
    end
    
    looping = 1;
    while looping
    
    if do_singleton
      fillind = inds(unidrnd(length(inds)));
      filler = indi(fillind);
    else
      filler = fillerpos(unidrnd(length(fillerpos)));
    end
    
    if mode == 1;
      dofill = (filler ~= thishole);
    else
      dofill = (min(abs(filler-holeregion)) > 0);
    end
     if dofill
        y(thishole) = x(filler);
        y(filler) = x(thishole);
	looping = 0;
	if do_singleton
	  mindiff(fillind) = NaN;
	end
     else 
        looping= looping + 1;
     end
     if looping > 10
       disp('looping > 10, exiting');
       return;
     end
    end
  end
end;



  
  
  
  


