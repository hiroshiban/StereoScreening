
% example of block permutation

% MAKE UP THE PERMUTED BLOCK DESIGNS
Q = 4;npts = 240;numblocks = 2;itype = 1;doshift = 2;
bstim = gen_block(Q,npts,numblocks,itype,doshift);
numperm = 20;
bstimmat = NaN*ones(numperm+1,npts);
bstimmat(1,:) = bstim;
for k = 1:numperm;
  bstimmat(k+1,:) = permute_block(bstimmat(k,:));
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
  [thiseff,thisdet]= calc_meffdet(bstimmat(k,:),nummods,Q,P,h);
  effmat(k) = thiseff(1);detmat(k)=thisdet(1);
end

entmat = calcentvec(bstimmat,2);

% Theoretical curves
[tdetmat,teffmat,tmaxdet,tmaxeff] = tcurve(Q,nummods,npts);

% Plot out results
figure(1);
subplot(2,1,1);%un-normalized
plot(detmat,effmat,'-o',tdetmat,teffmat,'-x');
subplot(2,1,2);% normalized
plot(detmat/tmaxdet,effmat/tmaxeff,'-o',tdetmat/tmaxdet,teffmat/tmaxeff,'-x');


