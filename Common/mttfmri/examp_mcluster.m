% example of m-sequence clustering

% MAKE UP THE m-sequence
Q = 4;npts = 240;
mparam.base = 5;mparam.power = 3; % gives length 124 point m-seq
mparam.shift = 0;
ms = gen_mseq(mparam,npts); % m-sequence
ms0 = ms(:,1); % pick the first m-sequence.
	       
	       
numperm = 20;
mstimmat = NaN*ones(numperm+1,npts);
mstimmat(1,:) = ms0(:)';
ind = 1;thistype = 1;
for k = 1:numperm;
  mstimmat(k+1,:) = clumpvec(mstimmat(k,:),thistype);
  ind = ind + 1;
  thistype = thistype + 1;
  if thistype > Q
    thistype = 1;
  end  
end


% EVALUATE THEIR PERFORMANCE
lorder = 3;nummods = 15;
S = legendremat(lorder,npts); % Nuisance functions
V = eye(npts,npts); % noise correlation matrix;
Vi = inv(V); % inverse of correlation matrix
P = Vi-Vi*S*inv(S'*Vi*S)*S'*Vi; % Projection matrix
h = hemoresp(0:(nummods-1),1.2,3,1); % HRF function
effmat = NaN*ones(numperm+1,1);
detmat = NaN*ones(numperm+1,1);
for k = 1:numperm
  k
  [thiseff,thisdet]= calc_meffdet(mstimmat(k,:),nummods,Q,P,h);
  effmat(k) = thiseff(1);detmat(k)=thisdet(1);
end

entmat = calcentvec(mstimmat,2);

% Theoretical curves
[tdetmat,teffmat,tmaxdet,tmaxeff] = tcurve(Q,nummods,npts);

% Plot out results
figure(1);
subplot(2,1,1);%un-normalized
plot(detmat,effmat,'-o',tdetmat,teffmat,'-x');
subplot(2,1,2);% normalized
plot(detmat/tmaxdet,effmat/tmaxeff,'-o',tdetmat/tmaxdet,teffmat/tmaxeff,'-x');


