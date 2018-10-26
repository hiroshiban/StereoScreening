% jellyfish for multiple event types
% for starting block patterns, a number of choices
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
%
%
% Version History
%
%   020824  TTL  renamed from bs_mult.m, making changes to 
%                be compatible with formulas in the paper. 


%defaults;
if ~exist('nevents')
	nevents = 2;
end
if ~exist('zeroarea');
	zeroarea = 0;
end
if ~exist('lorder');
	lorders = 0;
end

if ~exist('numones');
	numones = 64;
end
if ~exist('numperm')
	numperm = 100;
end
if ~exist('nummods')
	nummods = 15;
end
norders = length(lorders); % number of different orders for subspace
                         % interference. 

ntypes = 2; % types of block-starting pattern
%numones = 32;
width = 1;
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


switch( numones)   % this is the number of points per event type
case(64)
	numblocksvec = [ 1 2 4 8 16 0 ];
	offsetmat = [0 0 0 0 0 0;
		     0 0 0 0 0 0;
    		     0 0 0 0 0 0;
    		     0 0 0 0 0 0;];
	
otherwise
	error('invalid numones');
end



% Designs all assume that p = 1/(Q+1)
L = nevents + 1;   
npts = L*numones;  % TOTAL NUMBER OF POINTS IN THE DESIGN
nnblocks = length(numblocksvec);
baseveclen = npts;  % don't add zeros to end. 

%initialize design matrix
X = NaN*ones(npts,nummods*nevents);

%**************************
%make up contrast matrices
%**************************
%single event contrast
Amat = eye(nummods*nevents);
% between event contrasts
I = eye(nummods);
combos = nchoosek(1:nevents,2); % number of combinations of contrasts
ncontrasts = size(combos,1);
Cmat = zeros(ncontrasts*nummods,nummods*nevents);
for k = 1:ncontrasts;
	span0 = (1:nummods) + (k-1)*nummods;
	span1 = (1:nummods) + (combos(k,1)-1)*nummods;
	span2 =(1:nummods)  + (combos(k,2)-1)*nummods;
	Cmat(span0,span1) = I;
	Cmat(span0,span2) = -I;
end
neffdet = nevents + ncontrasts + 1;
%initialize effmat and detmat matrices;
effmat = NaN*ones(ntypes,neffdet,numperm,nnblocks,norders);
detmat = NaN*ones(ntypes,neffdet,numperm,nnblocks,norders);


% make up the matrix to project out the subspace interference.
PS = NaN*ones(baseveclen,baseveclen*norders);
for iorder = 1:norders
  lorder = lorders(iorder);
  S = legendremat(lorder,baseveclen);
  span = (1:baseveclen) + (iorder-1)*baseveclen;
  PS(:,span) = eye(baseveclen,baseveclen)-S*pinv(S);
end


% make up a big matrix to store the stimulus patterns.
stimmat = NaN*ones(ntypes,numperm,nnblocks,npts);

for blocknum = 1:nnblocks;
	blocknum
	%make up the stimulus pattern;
	numblocks = numblocksvec(blocknum);
	
	for itype = 1:ntypes
	basevec = zeros(baseveclen,1);
	
	% note, that in contrast to bigsim.m, here basevec
	% does not correspond to the stimulus pattern, but
	% rather codes for event types.
	
	if(numblocks > 0)
		dorand = 0;
		sizeofblock = numones/numblocks;
		blockspacing = npts/(numblocks);
		
		switch(itype)
		 case(1) % simple alternating
		         % A,B,C,Null,A,B,C,Null  
		for k = 1:nevents
		  thisbasevec =zeros(baseveclen,1);
		  startlocation = 1 + (k-1)*sizeofblock;
		  blocklocs = startlocation:blockspacing:npts;
		  thisbasevec(blocklocs) =1;
		  thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
		  basevec(find(thisbasevec)) = k;
		end
		
		 case(2) %zero-filled alternating.
		         % A,Null,B,Null,C,Null
		 eventspacing = sizeofblock + round(sizeofblock/nevents);
		 for k = 1:nevents
		   thisbasevec =zeros(baseveclen,1);
		   startlocation = 1 + (k-1)*eventspacing;
		   blocklocs = startlocation:blockspacing:npts;
		   thisbasevec(blocklocs) =1;
		   thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
		   basevec(find(thisbasevec)) = k;
		 end
		 
                case(3) %zero-filled clustered alternating.
		 
		otherwise;
		  
		  
	        end
	
		thisoffset = offsetmat(lorder+1,blocknum);
		basevec = shiftvec(basevec,thisoffset);
	else
		disp('random patterns');
		dorand = 1;
		basevecmat = unidrnd(npts,npts,numperm);
		basevecmat = rem(basevecmat,L);
	end



	if ~(dorand == 1 & itype > 1) 
		
		%now permute the basevector
		for pnum = 1:numperm
			if dorand
				basevec = basevecmat(:,pnum);	
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
			end
		
			% keep track of stimulus patterns
			stimmat(itype,pnum,blocknum,:) = basevec;
			
			
			%make up design matrix;
			for k = 1:nevents
			  thispattern = zeros(baseveclen,1);
			  thispattern(find(basevec == k)) = 1;
			  span = (1:nummods) + (k-1)*nummods;
			  X(:,span) = toeplitz(thispattern,[thispattern(1) zeros(1,nummods-1)]);
			end
		
			for iorder = 1:norders
			  lorder = lorders(iorder);
			  
			  %FULL COVARIANCE MATRIX
			  CX = inv(X'*squeeze(PS(lorder,:,:))*X);
			
			  %EFFICIENCIES
			  % efficiency of all estimates
			  effmat(itype,1,pnum,blocknum,iorder) = 1/trace(CX);
			  detmat(itype,1,pnum,blocknum,iorder) = h_all'*inv(CX)*h_all/h_all_norm2;
			  %efficiency and detection power of each event estimate
			  % note that for all events and contrasts, we assume that
			  % the expected value of A*h_all = h  -- that is 
			  % all events and contrasts are nominally a gamma density function.
			
			  for k = 1:nevents
			    span = (1:nummods) + (k-1)*nummods;
			    A = Amat(span,:);
			    thisCov = A*CX*A';
			    effmat(k+1,pnum,blocknum) = 1/trace(thisCov);
			    detmat(k+1,pnum,blocknum) = h'*inv(thisCov)*h/hnorm2;
			  end
		
			  %efficiency and detection power of contrasts
			  for k = 1:ncontrasts
			    span = (1:nummods) + (k-1)*nummods;
			    A = Cmat(span,:);
			    thisCov = A*CX*A';
			    effmat(k+nevents+1,pnum,blocknum) = 1/trace(thisCov);
			    detmat(k+nevents+1,pnum,blocknum) = h'*inv(thisCov)*h/hnorm2;
			  end
			
			end %iorder
			
			
		end %pnum
	end % if
	keyboard;

	end	%
end % 



