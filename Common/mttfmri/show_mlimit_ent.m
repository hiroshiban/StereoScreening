% show how close m-sequences come to the theoretical limit
% for designs with different numbers of events
%
%   nevent 
%     1          binary  m-seq  base 2
%     2          ternary m-seq  base 3
%     3          base 4 equence (extension field)
%     4          pentary m-seq base 5 
%     5          nothing
%     6          septary m-seq base 7, length 342
%     7          base 8 sequence (extension field)
%                
%     8          base 9 sequence (extension field)
%     9          nothing -- need mseq of base 10
%     10         base 11 mseq
%     12         base 13 mseq
%     15         addition of 4 binary sequences

if ~exist('doprint');doprint = 0;end;
nevents = [1  2  3 4  6 7 8 10 12];
base =    [2  3  4 5  7 8 9 11 13];
power =   [8  5  4 4  3 3 3 3 3];
repeats = [1  1  1 1  1 1 1 1 1];
lens = base.^power-1;
lens = lens.*repeats;
nummods = 15;
nm = length(nevents);
if ~exist('compdata');compdata = 1; end;


if compdata
  disp('compdata = 1; computing');
meff = NaN*ones(nm,1);
centropy = NaN*ones(nm,3);
mdet =  meff;
ieff = meff;
idet = mdet;

for k = 1:nm
  mparam.base = base(k);
  mparam.power = power(k);
  mparam.shift = 0;
  msmat = gen_mseq(mparam,lens(k));
  [thiseff,thisdet] = calc_meffdet(msmat(:,1),nummods);
  meff(k) = thiseff(1);mdet(k) = thisdet(1);
  
  npts = length(msmat(:,1));
  Q = nevents(k);
  ieff(k) = npts/(2*(Q+1)*nummods);
  idet(k) = npts*nummods/(2*(Q+1));
  centropy(k,1) = calcent(msmat(:,1),1);
  centropy(k,2) = calcent(msmat(:,1),2);
  centropy(k,3) = calcent(msmat(:,1),3);
end

else
  load showml_ent

end

iapprox = NaN*ones(nm,1)
iapproxd = NaN*ones(nm,1)
Mapprox = NaN*ones(nm,1);
% compute approximate bound
for k = 1:nm
Q = nevents(k);
npts = lens(k);
numones = npts/(nevents(k)+ 1);
atrace = approxtrace(nummods,numones,npts);
atrace = atrace/((1-numones/npts)*numones);
iapprox(k) = npts*atrace/(2*(Q+1)*nummods^2);
iapproxd(k) = npts*atrace/(2*(Q+1));  
Mapprox(k) = atrace;
end

  ltypes = ['ko:';'k- ';'k--'];
figure(1);
subplot(2,1,1);
hp = plot( nevents,meff./ieff,ltypes(1,:),...
	   [0 13],[1 1],ltypes(2,:),...
	   nevents,iapprox./ieff,ltypes(3,:));
set(gca,'FontSize',14,'FontName','Times');
legend('Estimation Efficiency','Theoretical Bound','Approximate Bound',4);
set(hp,'LineWidth',2);
set(gca,'Xlim',[0 13],'Ylim',[0.97 1.01]);
set(hp,'MarkerSize',8);
xlabel('Number of Trial Types Q');
ylabel('Normalized Efficiency');
title('(a) Estimation Efficiency');
grid;
subplot(2,1,2);

maxentropy = log2(nevents(:)+1);
nentropy = centropy ./ (maxentropy*ones(1,3));

hp2 = plot(nevents,nentropy(:,1),'k-o',nevents,nentropy(:,2),'k:x',...
	   nevents,nentropy(:,3),'k--d');
set(hp2,'MarkerSize',8);
%nevents,Mapprox/nummods^2,'k--');
set(gca,'FontSize',14,'FontName','Times');
set(hp2,'LineWidth',2);
set(gca,'Xlim',[0 13]);
set(gca,'Ylim',[-0.1 1.1]);
legend('1st order','2nd order','3rd order',3);
grid;

xlabel('Number of Trial Types Q');
ylabel('Normalized Conditional Entropy');
title('(b) Conditional Entropy');
save showml_ent nevents meff ieff mdet idet centropy

if doprint
print -dtiff NIMG771_fig5.tif
print -deps NIMG771_fig5.eps
end

