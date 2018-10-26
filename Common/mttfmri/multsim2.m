% jellyfish for multiple event types

%
%  (1) simple alternating:        1 2 0 1 2 0 1 2, etc.
%                                 1 2 3 0 1 2 3 0, etc.
%                                 
%  (2) zero-filled alternating    1 0/2 2 0/2 1 0/2 2 0/2
%                                 1 0/3 2 0/3 3 0/3 1 0/3 2 0/3 3 0/3
%
%  where 0/2 means half the number of zeros used in the blocks
%  of option (1).
%  
%  (3)  zero-filled clustered alternating:    1 0/2 1 0/2 2 0/2 2 0/2
%                                             1 0/3 1 0/3 2 0/3 2 0/3 3 0/3 3 0/3
%        -- this option is probably not that interesting from a cognitive
%           point of view so we won't worry about it. 
% Version History
%
%   020824  TTL  renamed from bs_mult.m, making changes to 
%                be compatible with formulas in the paper. 
%
%   021212  TTL  adding stuff for colored noise


%defaults;
if ~exist('nevents')
	nevents = 2;
end
if ~exist('zeroarea');
	zeroarea = 0;
end
if ~exist('lorders');
	lorders = [0:3];
end


if ~exist('numperm')
	numperm = 100;
end
if ~exist('nummods')
	nummods = 15;
end
if ~exist('doshift')
  doshift = 2;
end
if ~exist('numpermrand')
  numpermrand = numperm;
end
norders = length(lorders); % number of different orders for subspace
                         % interference. 

ntypes = 1; % types of block-starting pattern
%numones = 32;

%rand('state',100);
newrstate = rand('state');
load rstate2mult;
rand('state',newrstate);
if zeroarea
	h = diff(hemoresp(0:nummods,1.2,3,1));
	%h = h/sum(h);
else
	h = hemoresp(0:(nummods-1),1.2,3,1);
end;
h = h(:);
hnorm2 = h'*h;
h_all = repmat(h,nevents,1);
h_all_norm2 = h_all'*h_all;
hmat = kron(eye(nevents),h);


% Have 240 points per experiment, since this number
% is divisible by 2, 3, 4, 5, and 6, so that
% we can show Q = 1, 2, 3, 4, 5
Q = nevents;
L = Q + 1;
if ~exist('npts')
 npts = 240;           % keep length of experiment the same 
end
numones = npts/(Q+1);	
% notes: in numblocksvec
%        0 = do random designs
%        -1 = do m-sequence based designs
nmshifts = 0;
doclump = 0;

switch(numones)
 case 120 % 1 event for 240 point sequence
  numblocksvec = [1 2  4 8 15 0 -1];
  mparam.base = 2;mparam.power = 8;
  mparam.sht = 0;mrepeat = 1;
  nmshifts =return_mtaps(mparam.base,mparam.power)
  doclump = 1;nperm_clump = 100;
 case 80 % 2 events for 240 point long sequence
  numblocksvec = [  1 2 4 8 16 40 0 -1];
%  numblocksvec = [8 16 40 0 -1]
  mparam.base = 3;mparam.power = 5; % gives length 242 m-seq
  mparam.shift = 0;mrepeat = 1;
  nmshifts =return_mtaps(mparam.base,mparam.power)
  doclump = 1;nperm_clump = 100;
 case 81 % 2 events for 243 point long sequence
    numblocksvec = [ 3 9 27 0 -1];
    mparam.base = 3;mparam.power = 5; % gives length 242 m-seq
    mparam.shift = 0;
    nmshifts =return_mtaps(mparam.base,mparam.power)
 case 60 % 3 events for 240 point long sequence
  numblocksvec =  [1 2 4 10 15 30 0 -1  ];
  mparam.base = 2;mparam.power = 8; % gives length 255 m-seq
  mparam.shift = 0;
  mshifts = 1:round(npts/2);% this was an interesting exercise
                            % with effdet decreasing but no
			    % gain in detpower!
  mshifts = [15]; % optimum found empirically
  nmshifts =length(mshifts);
 case 64 % 3 events for 256 point long sequence
  numblocksvec = [ 1 2 4 8 16 32 0 -1];
  mparam.base = 2;mparam.power = 8; % gives length 255 m-seq
  mparam.shift = 0;
  mshifts = [15]; % optimum found empirically
  nmshifts =length(mshifts);
  
 case 48 % 4 events for 240 point sequence
  
  numblocksvec =  [1 2 4 8 12 24 0  -1];
%  numblocksvec = -1;
  mparam.base = 5;mparam.power = 3; % gives length 124 point m-seq
  mparam.shift = 0;
  nmshifts =return_mtaps(mparam.base,mparam.power)
  mrepeat = 2; % how many times to repeat the m-sequence to 
               % get up to the desired length
  doclump = 1;nperm_clump = 100;
  
	       
 case 40 % 5 events for 240 point sequence
  numblocksvec =  [1 2 4 8 10 20 0 ];
  
 case 24 % 4 events for a 120 point long sequence
    numblocksvec = [1 2 3 4 6 8 12 0 -1];
    mrepeat = 1;
    mparam.base = 5;mparam.power = 3; % gives length 124 point m-seq
    mparam.shift = 0;
    nmshifts =return_mtaps(mparam.base,mparam.power)
  
 otherwise
     error('invalid numones')
end
if doclump == 1
  numperm_mseq = nperm_clump;
  nmshifts = 1;
else
  numperm_mseq = nmshifts;
end
nnblocks = length(numblocksvec);
baseveclen = npts;  % don't add zeros to end. 

%initialize design matrix
X = NaN*ones(npts,nummods*nevents);
if nevents > 1
  combos = nchoosek(1:nevents,2); % number of combinations of contrasts
  ncontrasts = size(combos,1);
else
  ncontrasts = 0;
end
neffdet = nevents + ncontrasts + 1;
%initialize effmat and detmat matrices;

effmat = NaN*ones(ntypes,norders,neffdet,numperm,nnblocks);
detmat = NaN*ones(ntypes,norders,neffdet,numperm,nnblocks);

reffmat = NaN*ones(ntypes,norders,neffdet,numpermrand,1);
rdetmat = NaN*ones(ntypes,norders,neffdet,numpermrand,1);

meffmat = NaN*ones(ntypes,norders,neffdet,numperm_mseq,1);
mdetmat = NaN*ones(ntypes,norders,neffdet,numperm_mseq,1);


% make up the matrix to project out the subspace interference.
PS = NaN*ones(norders,baseveclen,baseveclen);
if docolornoise
  % Burock and Dale autocorrlelation model
    rho = 0.88; lambda = 0.75;
    V = lambda*eye(baseveclen,baseveclen)  +  ...
         (1-lambda)*toeplitz(rho.^(0:baseveclen-1));
    Vi = inv(V); %inverse of the autocovariance matrix
  for iorder = 1:norders
    lorder = lorders(iorder);
    S = legendremat(lorder,baseveclen);
    PS(iorder,:,:) = Vi-Vi*S*inv(S'*Vi*S)*S'*Vi;
  end
else
  rho = NaN; lambda = NaN;
  for iorder = 1:norders
    lorder = lorders(iorder);
    S = legendremat(lorder,baseveclen);
    PS(iorder,:,:) = eye(baseveclen,baseveclen)-S*pinv(S);
  end
end



% make up a big matrix to store the stimulus patterns.
stimmat = NaN*ones(ntypes,numperm,nnblocks,npts);
rstimmat = NaN*ones(ntypes,numpermrand,1,npts);
mstimmat = NaN*ones(ntypes,numperm_mseq,1,npts);

rd = 0;%number of rank deficient matrices
for blocknum = 1:nnblocks;
  t0 = clock;

	blocknum
	%make up the stimulus pattern;
	numblocks = numblocksvec(blocknum);
	
	for itype = 1:ntypes
	basevec = zeros(baseveclen,1);
	
	% note, that in contrast to bigsim.m, here basevec
	% does not correspond to the stimulus pattern, but
	% rather codes for event types.
	
	if(numblocks > 0)
		dorand = 0; domseq = 0;
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
		
%		figure(1);plot([basevec basevec0],'-o');
%		keyboard;
		
		
	elseif (numblocks == 0)
		disp('random patterns');
		dorand = 1;domseq = 0;
		randmat = unidrnd(npts,npts,numpermrand);
		[ignore,prandmat] = sort(randmat);
		basevecmat = rem(prandmat,L);
		
	elseif (numblocks == -1)
	        disp('m-sequence');
		dorand = 0;domseq = 1;
		if(nevents == 2) | (nevents == 4) % use ternary or
                                                  % pentary m-sequences
		 if (nevents == 2) fixEv = 1;
		 elseif (nevents == 4) fixEv = 2;
		 end
		 basevecmat =NaN*ones(npts,numperm_mseq);

		 for ws = 1:nmshifts
		   ms=mseq(mparam.base,mparam.power,mparam.shift,ws)'+ ...
		    fixEv;
		   ms = ms(:);
		   ms = repmat(ms,mrepeat,1);
		   
		   lms = length(ms);
		   %zeropad or truncate sequence
		   if(npts > lms) % zeropad
		     ms = [ms;zeros(npts-lms,1)];
		   else %truncate
		     ms = ms(1:npts);
		   end
		   basevecmat(:,ws) = ms;
		 end
		 
		 if doclump
		   ind = 1;thistype = 1;
		   for k = 1:(nperm_clump-1)
		     basevecmat(:,ind+1) = clumpvec(basevecmat(:,ind), ...
						    thistype);
		     ind = ind + 1;
		     thistype = thistype + 1;
		     if thistype > nevents
		       thistype = 1;
		     end
		   end
		 end
		 

		 
		elseif (nevents == 3) % use sum of binary sequences
  		 baesvecmat =NaN*ones(npts,nmshifts);
		   ws = 1;
 		   ms=mseq(mparam.base,mparam.power,mparam.shift,ws)';
   		   ms = m2bin(ms(:));
		   lms = length(ms);
		   %zeropad or truncate sequence
		   if(npts > lms) % zeropad
		     ms = [ms;zeros(npts-lms,1)];
		   else %truncate
		     ms = ms(1:npts);
		   end
		   for shift = 1:nmshifts
		     basevecmat(:,shift) = ms(:) + 2*circshift(ms(:),mshifts(shift));
		   end
		  
		else
		
		  fprintf('m-sequences not supported for nevents = %d',nevents);
                end
		
	end



	if 1
		
	  %now permute the basevector
	  if dorand == 1
	     numperm_all = numpermrand;
	  elseif domseq == 1
	     numperm_all = numperm_mseq;
	  else
	     numperm_all = numperm;
	  end
	  for pnum = 1:numperm_all
	    if (dorand == 1) | (domseq == 1)
	      basevec = basevecmat(:,pnum);	
	      % keep track of stimulus patterns
	      if(dorand == 1) 
		rstimmat(itype,pnum,1,:) = basevec;
	      else
		mstimmat(itype,pnum,1,:) = basevec;
	      end
	    else
		
	      if(pnum > 1)
		looping = 1;
		while looping % loop until we find two points with different values
		  p1 = unidrnd(npts);
		  p2 = unidrnd(npts);
		  if basevec(p1) ~= basevec(p2)
		    tmp = basevec(p1);
		    basevec(p1)=basevec(p2);
		    basevec(p2)=tmp;
		    looping = 0;
		  end
		end
				
	      end
	      
       	      % keep track of stimulus patterns
    	      stimmat(itype,pnum,blocknum,:) = basevec;
	    end
		

	    for iorder = 1:norders;
	        thisPS = squeeze(PS(iorder,:,:));
		[thiseff,thisdet] = calc_meffdet(basevec,nummods,nevents,thisPS,h);
		if (dorand == 1)
		  reffmat(itype,iorder,:,pnum,1) = thiseff;
		  rdetmat(itype,iorder,:,pnum,1) = thisdet;
		elseif (domseq == 1)
		  meffmat(itype,iorder,:,pnum,1) = thiseff;
		  mdetmat(itype,iorder,:,pnum,1) = thisdet;
		else
		  effmat(itype,iorder,:,pnum,blocknum) = thiseff;
		  detmat(itype,iorder,:,pnum,blocknum) = thisdet;
		end
	    end
			
	  end % pnum
	end	%itype
	end
	etime(clock,t0)
end %blocknum

if 1
filename = sprintf('ne%dl%sno%dz%dds%dc%dv2cn%d',nevents, ...
		   no_blanks(num2str(lorders)),numones,zeroarea,doshift,doclump,docolornoise);
varnames = ['reffmat rdetmat effmat detmat stimmat lorders filename'];
varnames = [varnames,' nummods npts numones numblocksvec'];
varnames = [varnames,' meffmat mdetmat rstimmat mstimmat'];
varnames = [varnames,' rho lambda'];

eval(sprintf('save %s %s',filename,varnames));
end;





