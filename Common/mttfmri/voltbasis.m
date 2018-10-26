function [B1,B2] = voltbasis(vdim,bdim,diffmode,dt)

%  [B1,B2] = voltbasis(vdim,bdim,diffmode,dt)

% return matrices B1 and B2 which can be used to 
% multiply X and X2 to do a volterra analysis using
% basis functions.
%
% MODS
%  020516 TTL initial version
%  040323 TTL added comments
%  
% send questions to ttliu@ucsd.edu
  
  

if ~exist('dt')
	dt = 1; % default is 1 second resolution
end

t = (0:(vdim-1))*dt; % timeaxis

if bdim <= 3
	alpha = [4 8 16];
else
	error('voltbasis; bdim not supported');
end

if diffmode == 1;
 matdim = 2*bdim;
else
 matdim = bdim;
end
	
 b = NaN*ones(vdim,matdim);
% make up the basis functions
for k = 1:bdim
	[b(:,k),db] = hemoresp(t,1,alpha(k)-1,1);
	if diffmode == 1
	  b(:,k+bdim) = db; % derivative of the basis functions. 
	end
end


B1 = b; % B1 is just the matrix of these basis functions as column vectors
B2len = sum(1:vdim);
B2 = NaN*ones(B2len,matdim^2);
index = 1;
for k = 1:matdim
	for l = 1:matdim
		thismat = b(:,k)*b(:,l)';
		thisvec =[];
	    for m = 1:vdim
			if m > 1
				thisvec = [thisvec;2*diag(thismat,m-1)];
			else
				thisvec = [thisvec;diag(thismat,m-1)];
			end
		end
	%	keyboard;pause;
		B2(:,index) = thisvec;
		index = index + 1;
	end
end



