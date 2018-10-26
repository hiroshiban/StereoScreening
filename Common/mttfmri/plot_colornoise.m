% Show effect of correlated noise
%
% Files are created by batch_color.m
% 


if ~exist('donimg');donimg = 1;end;

cmspan = [2:100];
fsize = 12;

filenames = ['ne2l0123no80z0ds2c1v3cn1ba0    ';
	     'ne4l0123no48z0ds2c1v3cn1ba0    ';];	     

titlstrs = str2mat('(a) Q = 2, correlated noise',...
		   '(b) Q = 4, correlated noise');



noiseon = [1 1];
fspan = [1 2 ];
filenames = filenames(fspan,:);
noiseon = noiseon(fspan);
titlstrs = titlstrs(fspan,:);

nummods = 15;baseveclen = 240;
% as a quick and dirty kludge
% multiply the noise detection powers 
% by a scaling factor of the form h'*h/(h'*Vi*h)
h = hemoresp(0:(nummods-1),1.2,3,1);
rho = 0.88; lambda = 0.75;
V = lambda*eye(baseveclen,baseveclen)  +  ...
         (1-lambda)*toeplitz(rho.^(0:baseveclen-1));
Vi = inv(V); %inverse of the autocovariance matrix
span = 2+(0:(nummods-1));
Vi = Vi(span,span);
hscale = (h'*h)/(h'*Vi*h);


figure(1);clf;
nf = size(filenames,1);
ind= 1;dind = 1;itype = 1;
tind = 1;
nr = 1;nc = 2;
for nfile = 1:nf
  if nfile  == 5
    figure(2);clf;ind = 1;
  end;

if noiseon(nfile);
  dscale = 1.0;
  tdscale = 1/hscale;
  escale = nummods/trace(Vi);
  escale = 1.0;
else
  dscale = 1.0;
  tdscale = 1.0;
  escale = 1.0;
end

  filename = deblank(filenames(nfile,:));
  
eval(sprintf('load %s',filename));

nevents = str2num(filename(3))
Q = nevents;
% UPPER BOUNDS
atrace = approxtrace(nummods,numones,npts);
atrace = atrace/((1-numones/npts)*numones*nummods);
tmaxeff = npts/(2*(Q+1)*nummods);

tmaxdet = npts*nummods/(2*(Q+1));
neffmat = escale*effmat/tmaxeff;
ndetmat = dscale*detmat/tmaxdet;
nreffmat = escale*reffmat/tmaxeff;
nrdetmat = dscale*rdetmat/tmaxdet;
nmeffmat = escale*meffmat/tmaxeff;
nmdetmat = dscale*mdetmat/tmaxdet;

% THEORETICAL CURVES
k = nummods;
i0 = atrace;
i0 = 1.0;
alpha =[1/k 0.1:0.05:1.0];
thetavec = [60   90]/180*pi;
teff = alpha.*(1-alpha) ./(1+alpha*(k^2-2*k))*i0*nummods*npts/(2*(Q+1));
alphamat = alpha(:)*ones(1,length(thetavec));
teffmat = teff(:)*ones(1,length(thetavec));
xrangemat = ones(length(alpha),1)*thetavec(:)';

tdet =i0*nummods*((1-alphamat).*(sin(xrangemat).^2)/(k-1)+alphamat.* ...
	  cos(xrangemat).^2);
tdet = tdet*npts/(2*(Q+1));
ntdet = tdet/tmaxdet;
nteffmat = teffmat/tmaxeff;

if ~exist('docolor');docolor = 0;end
if docolor
  ltype = ['o';
	 's';
	 '^';
	 '+';
	 '*';
	 '.';
	 'd';
	 'x';];
  ttype = '-';
else
  ltype = ['k+ ';
	 'k+ ';
	 'k+ ';
	 'k+ ';
	 'k+ ';
	 'k+ ';
	 'kd ';
	 'k-o';];
  ttype = 'k-';
  
  
end  



  for iorder = 1;
    
    
        % PICK BEST 10 RANDOM DESIGNS
	[s,s_ind] = sort(squeeze(nreffmat(itype,iorder,dind,:,1)));
	s_ind = flipud(s_ind(:));
	subplot(nr,nc,ind);
	if donimg & ind ~= 4 & ind ~=3
	hp = plot(squeeze(ndetmat(itype,iorder,dind,:,1)),squeeze(neffmat(itype,iorder,dind,:,1)),'m+',...
	     squeeze(ndetmat(itype,iorder,dind,:,6)),squeeze(neffmat(itype, ...
						  iorder,dind,:,6)),'g+',...
		  squeeze(nmdetmat(itype,iorder,dind,1,1)),squeeze(nmeffmat(itype, ...
						  iorder,dind,1,1)),'ro',...
		  squeeze(nmdetmat(itype,iorder,dind,cmspan,1)),squeeze(nmeffmat(itype, ...
						  iorder,dind,cmspan,1)),'r+');

	hold on;
	end
	set(hp(3),'MarkerFaceColor','r','MarkerSize',8);
	if ind == 1;
	  	hleg = legend('1-block','40-block','m-sequence',...
			     'clustered m-seq',1);
         else
	  	hleg = legend('1-block','24-block','m-sequence',...
			     'clustered m-seq',1);
	 end
	set(hleg,'Visible','off');set(allchild(hleg),'visible','on');
	
	if donimg & ind ~= 4 & ind ~=3
	plot(tdscale*ntdet,nteffmat,ttype);
	if noiseon(nfile)
	  axis([0 0.16 0 1.1]);	  
	else
	axis([0 0.6 0 1.1]);	  
	end

	else
	  if noiseon(nfile)
	    axis([0.03 0.07 .95 1.05]);
	  else
	  axis([0.06 0.1 .95 1.05]);
	  end
	end
	% plot soft bounds
%	plot([0 0.6],[atrace atrace],'--');
	
	hold off;

	grid;
	set(gca,'Fontsize',fsize);
	ylabel('Normalized \xi_{tot}')
	xlabel('Normalized R_{tot}')
	title(titlstrs(tind,:));

  end
  
  if ind == 3
    hold on;
    arrow([.037,1.00],[.041,1.035]);
    ht = text(.032, .995,'m-sequence');
    set(ht,'FontSize',fsize);
    arrow([.04,.96],[.055,1.0]);
    ht = text(.035, .955,'clustered m-sequences');
    set(ht,'FontSize',fsize);
    hold off;
  end
  if ind == 4
    hold on;
    arrow([.037,1.00],[.042,1.023]);
    ht = text(.032, .995,'m-sequence');
    set(ht,'FontSize',fsize);
    arrow([.04,.96],[.05,.99]);
    ht = text(.035, .955,'clustered m-sequences');
    set(ht,'FontSize',fsize);
    hold off;
  end

ind = ind + 1;
tind = tind + 1;

end