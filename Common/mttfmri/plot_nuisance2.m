% plot out effects of nuisance parameters
% for trial types 2,3,4,5, plot the initial detection power 
%  as a function of the number of nuisance parameters.



load nipaper_01_ent;
titlstrs = str2mat('(a) 1-block detection power',...
		   '(b) 2-block detection power',...
		   '(c) 4-block detection power',...
		   '(d) 8 or 10-block detection power');
fsize = 12;
randomspan = 1;
blockspan = [1:100 110:10:200];

filenames = ['ne2l0123no80z0ds2c0v3cn0ba0np1';
	     'ne3l0123no60z0ds2c0v3cn0ba0np1';
     	     'ne4l0123no48z0ds2c0v3cn0ba0np1';
       	     'ne5l0123no40z0ds2c0v3cn0ba0np1';];


fontsize = 12;


nf = size(filenames,1);
ind= 1;dind = 1;itype = 1;
nr = 2;nc = 2;
figure(1);clf;

% READ IN EACH FILE AT A TIME
if ~exist('loaddata');loaddata = 1;end;
if loaddata;
  detpowers = NaN*ones(nf,4,6);
for nfile = 1:nf;




  filename = deblank(filenames(nfile,:));
  
eval(sprintf('load %s',filename));

nevents = str2num(filename(3))
Q = nevents;
% UPPER BOUNDS
atrace = approxtrace(nummods,numones,npts);
atrace = atrace/((1-numones/npts)*numones*nummods);
tmaxeff = npts/(2*(Q+1)*nummods);

tmaxdet = npts*nummods/(2*(Q+1));
neffmat = effmat/tmaxeff;
ndetmat = detmat/tmaxdet;
nreffmat = reffmat/tmaxeff;
nrdetmat = rdetmat/tmaxdet;
nmeffmat = meffmat/tmaxeff;
nmdetmat = mdetmat/tmaxdet;

% THEORETICAL CURVES
k = nummods;
i0 = 1.0;
alpha =[1/k 0.1:0.05:1.0];
thetavec = [45  55 65  80 90]/180*pi;
teff = alpha.*(1-alpha) ./(1+alpha*(k^2-2*k))*i0*nummods*npts/(2*(Q+1));
alphamat = alpha(:)*ones(1,length(thetavec));
teffmat = teff(:)*ones(1,length(thetavec));
xrangemat = ones(length(alpha),1)*thetavec(:)';

tdet =i0*nummods*((1-alphamat).*(sin(xrangemat).^2)/(k-1)+alphamat.* ...
	  cos(xrangemat).^2);
tdet = tdet*npts/(2*(Q+1));
ntdet = tdet/tmaxdet;
nteffmat = teffmat/tmaxeff;


detpowers(nfile,:,:) = squeeze(ndetmat(1,:,1,1,1:6));
meffs(nfile,:,:) = squeeze(nmeffmat(1,:,1));

end;
end;

nr = 2;nc = 2;clf; hold off;

xax = 0:3;
for ind = 1:4;
  subplot(2,2,ind);
  plot(xax,squeeze(detpowers(1,:,ind)),'-o',...
       xax,squeeze(detpowers(2,:,ind)),'-s',...
       xax,squeeze(detpowers(3,:,ind)),'-d',...
       xax,squeeze(detpowers(4,:,ind)),'-p');
  axis([-0.5 3.5  0 0.5]);
  set(gca,'Fontsize',fsize);
  xlabel('Highest Nuisance Order');
  ylabel('Normalized R_{tot}');
  title(titlstrs(ind,:));
  grid;
  if ind == 2
    legend('2 trial types','3 trial types','4 trial types',...
		    ' 5 trial types',4);
    end
end

   

if ~exist('doprint');doprint=0 ;end  
if doprint  
print -dtiff NIMG771_fig2.tif
print -depsc NIMG771_fig2.eps
end
