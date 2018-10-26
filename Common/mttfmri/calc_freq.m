

nevents = 10;
minmax = NaN*ones(nm,2);
freqmat = NaN*ones(nm,nevents + 1);
for j = 1:nm
  thismin = Inf; thismax = 0;
  for k = 1:nevents
    nk = sum(basevecmat(:,j) == k);
    if nk > thismax; thismax = nk;end
    if nk < thismin; thismin = nk;end
  end
  minmax(j,:) = [thismin thismax];
  for k = 0:nevents;
    freqmat(j,k+1) = sum(basevecmat(:,j) == k);
  end
end
  