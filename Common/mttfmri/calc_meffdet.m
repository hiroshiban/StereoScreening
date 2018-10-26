function [effmat,detmat,ueffmat,udetmat] = ...
          calc_meffdet(basevec,nummods,nevents,PS,h,V,uC,mode,approxflag,B);

% Usage:
% -----
%
% [effmat,detmat,ueffmat,udetmat] = ...
%            calc_meffdet(basevec,nummods,nevents,PS,h,V,uC,mode);
%
% Inputs
%-------
%  basevec:  stimulus pattern
%  nummods:  number of model functions
%  events :  [optional]  number of events, not including control state 
%  PS     :  [optional]  matrix to project out nuisance terms
%                        if PS is an integer, then it is interpreted
%                        as the order of the subspace S
%  h      :  [optional]  hemodynamic response function
%  V      :  [optional] noise correlation matrix, default = identity,
%                       if V is empty, V = I.
%                       if approxflag = 1, then V should be the square root
%                       of the inverse of the covariance matrix
%                       i.e., sqrtm(inv(V))
%  uC     :  [optional] user specified matrix of contrasts
%                       ignored if uC is empty matrix
%  mode   :  [optional ]  0: compute overall power and eff
%                         1: compute both overall and user power and eff
%                         2: compute only user power and eff
%                       defaults are:
%                         0 if uC is non-existent or empty
%                         1 if user has specified non-empty uC
%
%  approxflag :   [optional]   0: compute power and eff in the normal manner
%                         1: use approximation to compute power and eff
%                         default = 0;
%          B  :  [optional] matrix with basis functions as columns
%                           ignored if empty
%                           default = [];
%
%                         
% Outputs
%-------
% effmat: efficiency vector for general contrasts, first element is overall efficiency
%         other elements are the variances for each contrast
%
% detmat: detection power vector for general contrasts, first element is overall detection power
%         other elements are the variances for each contrast
%
% ueffmat: efficiency vector for user contrasts, first element is overall efficiency
%         other elements are the variances for each contrast
%
% udetmat: detection power vector for user contrasts, first element is overall detection power
%         other elements are the variances for each contrast
%
% 
%
% Modification History:
%----------------------
%  020929 TTL added user contrasts input
%  021213 TTL added some comments 
%             added a approxflag mode for looking at approximations
%             to the grammian matrix
%  021215 TTL added basis function expansion mode
% 
% 
% Copyright (c) 2002 UC Regents
%
% Send comments and questions to ttliu@ucsd.edu

  if ~exist('B');
    B = [];
    dobasis = 0;
  else 
    if isempty(B)
      dobasis = 0;
    else
      dobasis = 1;
    end
  end

  npts = length(basevec);
  if ~exist('nevents');
   nevents = max(basevec);
  end
  
  if ~exist('approxflag'); approxflag = 0; end

  % hemodynamic response
  if ~exist('h');
    	h = hemoresp(0:(nummods-1),1.2,3,1);
  else
    if length(h) ~= nummods
       error('length of h must be equal to nummods');
    end 
  end
  h = h(:);
  
  if dobasis
    nbasis = size(B,2)/nevents;
    Bm = B(1:nummods,1:nbasis);
    coeff = inv(Bm'*Bm)*Bm'*h;
    h = coeff;
  end


  hnorm2 = h'*h;
  
  % noise correlation matrix
  if ~exist('V')
    noisecor = 0;
  else 
    if ~isempty(V)
      noisecor = 0;
    else
      Vi = inv(V);
      noisecor = 1;
    end
  end
  
  % nuisance matrix, default is DC 
  makePS = 0;
  if ~exist('PS')
      lorder = 0;
      makePS = 1;
  else
    if(length(PS) == 1)
      lorder = PS;
      makePS = 1;
    end
  end
  
  if makePS
    S = legendremat(lorder,npts);
    if noisecor
      % note this is not a projection matrix!
      % but is analagous to K in the mult paper!
      % it's done this way here for convenience
      % so that matrix square roots are not necessary

      PS = Vi-Vi*S*inv(S'*Vi*S)*S'*Vi;
    else
      PS = eye(npts,npts)-S*pinv(S);
    end
  end
  
  % Handling of user specified contrasts
  n_ucontrasts = 0;
  if ~exist('uC')
    mode = 0;
  else
    if isempty(uC)
      mode = 0;
    else
      if size(uC,2) ~= nevents
	error('number of columns in uC must be equal to nevents');
      end	
      n_ucontrasts = size(uC,1);
      if ~exist('mode') mode = 1; end
    end      
  end
  

  hmat = kron(eye(nevents),h);
  baseveclen = length(basevec);
  
  %initialize design matrix
  X = NaN*ones(npts,nummods*nevents);
  if nevents > 1
   combos = nchoosek(1:nevents,2); % number of combinations of contrasts
   ncontrasts = size(combos,1);
  else
   combos = 0;
   ncontrasts = 0;
  end    
  neffdet = nevents + ncontrasts + 1;
  effmat = NaN*ones(neffdet,1);
  detmat = NaN*ones(neffdet,1);
  ueffmat = NaN*ones(n_ucontrasts + 1,1);
  udetmat = NaN*ones(n_ucontrasts + 1,1);

  %make up design matrix;
  for k = 1:nevents
    thispattern = zeros(baseveclen,1);
    thispattern(find(basevec == k)) = 1;
    span = (1:nummods) + (k-1)*nummods;
    X(:,span) = toeplitz(thispattern,[thispattern(1) zeros(1,nummods-1)]);
  end
		
			  
  %CHECK THE RANK

  if dobasis;% basis function expansion
    X = X*B;
  end

  
  K = X'*PS*X;
  if approxflag % note in approxflag mode, PS should be 
           % the complementary projection onto S and not onto sqrtm(Vi)*S
    K = V'*K*V;
  end

  if(rank(K) < size(X,2));		% rank deficient for efficiency
    effmat(1) = 0;disp('rank deficient for eff');doeff = 0;
  else;
    doeff = 1;
    CX = inv(K);		
  end
  if dobasis
    CX = B*CX*B';
  end

  G = hmat'*K*hmat;
  if(rank(G) < nevents)
    detmat(1) = 0;disp('rank deficient for power');dodet = 0;
  else
    dodet = 1;
    CZ = inv(G);
  end
		

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculation of overall efficiency and contrasts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if mode <= 1
    %EFFICIENCIES
                
    %efficiency and detection power of each event estimate
    for k = 1:nevents
      span = (1:nummods) + (k-1)*nummods;
      D = dcon(k,k,nevents);
      if dodet
	detmat(k+1) = D*CZ*D'* ...
	    hnorm2;
      end
		  
      if doeff
	A = kron(D,eye(nummods));
	effmat(k+1) = trace(A*CX*A');
      end
    end
		
    %efficiency and detection power of contrasts
    for k = 1:ncontrasts
      span = (1:nummods) + (k-1)*nummods;
      
      if dodet
	D = dcon(combos(k,1),combos(k,2),nevents);
	detmat(k+nevents+1) = ...
	    D*CZ*D'*hnorm2;
      end
      
      if doeff
	A = kron(D,eye(nummods));
	effmat(k+nevents+1) = trace(A*CX*A');
      end
    end

    % overall efficiency
    if doeff
      effmat(1) = ...
	  1/mean(effmat(2:neffdet));
    end
    if dodet
      detmat(1) = ...
	  1/mean(detmat(2:neffdet));
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculation of user-specified efficiency and contrasts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if mode >= 1
    for k = 1:n_ucontrasts
      span = (1:nummods) + (k-1)*nummods;
      if dodet
	D = uC(k,:);
	udetmat(k+1) = ...
	    D*CZ*D'*hnorm2;
      end
      if doeff
	A = kron(D,eye(nummods));
	ueffmat(k+1) = trace(A*CX*A');
      end
    end
    if doeff
      ueffmat(1) = ...
	  1/mean(ueffmat(2:n_ucontrasts));
    end
    if dodet
      udetmat(1) = ...
	  1/mean(udetmat(2:n_ucontrasts));
    end
  end


