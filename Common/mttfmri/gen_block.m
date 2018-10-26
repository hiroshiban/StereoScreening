function bstim = gen_block(Q,npts,numblocks,itype,doshift)

% Generate a block design
%  
%   bstim = gen_block(Q,npts,numblocks,itype,doshift)
%
% Inputs:
%
%   Q = number of non-Null trial types
%   npts = number of points in the design
%   numblocks = number of blocks in the design 
%   itype:    1 =   simple alternating, A,B,C,NULL,A,B,C
%             2 =   %zero-filled alternating. A,NULL,B,NULL,C,NULL
%             default is 1.
%   doshift:  0 = don't shift
%             1 = shift block design to center it
%             2 = shift only if itype  = 1
%             default is 2.
%
% Output
%   bstim: block design  

  if ~exist('itype');itype = 1;end
  if ~exist('doshift');doshift = 2;end;

  nevents = Q;
  numones = npts/(Q+1);	  
  baseveclen = npts;
  basevec = zeros(baseveclen,1);
	

    sizeofblock = round(numones/numblocks);
    blockspacing = npts/(numblocks);
		
    switch(itype)
     case(1) % simple alternating, A,B,C,NULL,A,B,C
      for k = 1:nevents
	thisbasevec =zeros(baseveclen,1);
	startlocation = 1 + (k-1)*sizeofblock;
	blocklocs = startlocation:blockspacing:npts;
	thisbasevec(blocklocs) =1;
	thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
	basevec(find(thisbasevec)) = k;
      end
		
     case(2) %zero-filled alternating. A,NULL,B,NULL,C,NULL
      zerofill = round((npts-nevents*numones)/(nevents* ...
					       numblocks));
      eventspacing = sizeofblock + zerofill;
      for k = 1:nevents
	thisbasevec =zeros(baseveclen,1);
	startlocation = 1 + (k-1)*eventspacing;
	blocklocs = startlocation:blockspacing:npts;
	thisbasevec(blocklocs) =1;
	thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
	basevec(find(thisbasevec)) = k;
      end
      
     otherwise;
      error('invalid choice of itype');
		  
    end
	
    % as a reasonable guess at the best offset
    % shift the base vector so that the number of
    % zeros at the beginning and end are about the same.
    % note that it would probably be better to do this
    % differently for each lorder (i.e. the reason for the
    % offsetmat matrix), but the choice of optimum offset
    % is not clear for multiple event types -- perhaps the
    % offset that makes the detection powers as equal as
    % possible
    % at each lorder.
		
    basevec0 = basevec;
    
    if (doshift == 1) | (doshift == 2 & itype == 1)
      lastnonzero = max(find(basevec > 0));
      nlastzeros = npts-lastnonzero;
		
      thisoffset = round(nlastzeros/2);
      basevec = shiftvec(basevec,thisoffset);
    end
    
    bstim = basevec(:)';
		
