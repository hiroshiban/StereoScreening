function lmat = legendremat(m,n,recur);
 
% lmat = legendremat(n,m);
%
% makes a (n)x(m+1) matrix where the columns are legendre polynomials
% with unit norm.
%
%   n = number of time points;
%   m = maximum order of legendred polynomial
%       (0th order = DC, 1st order = linear trend, etc.)
%   recur: 1: use recursive routine; 0 use hardwired routines up to m = 6
%
% send comments to ttliu@ucsd.edu

if ~exist('recur') & m <=6
	recur = 0;
end

if (mod(n,2) == 0) %even number of points;
	taxis = ((-n/2+.5):(n/2-0.5))/(n/2-0.5);
else
	taxis = (-(n-1)/2:(n-1)/2)/(n/2);
end
taxis = taxis(:);
lmat = NaN*ones(n,m+1);

if recur
	for k = 1:(m+1);
		lmat(:,k) = legpoly(k-1,taxis);
		lmat(:,k) = lmat(:,k)/norm(lmat(:,k));
	end
else
	for k = 1:(m+1)
		switch (k-1)
			case 0
				pl = ones(size(taxis));
				
			case 1
				pl = taxis;
			case 2
				pl = (3*taxis.^2-1)/2;
			case 3
				pl = (5*taxis.^3-3*taxis)/2;
			case 4
				pl = (35*taxis.^4-30*taxis.^2+3)/8;
			case 5
				pl = (63*taxis.^5-70*taxis.^3+15*taxis)/8;
			case 6
				pl=(231*taxis.^6 - 315*taxis.^4+105*taxis.^2-5);

				
			end
			
		lmat(:,k) = pl / norm(pl);
		
		
	end
end


