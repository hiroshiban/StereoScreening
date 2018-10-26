% plot out results from multsim1.m
% 03 May 12th making some modifications to look at the entropy of the 
% stimulus patterns.


maxentorder = 3;
if ~exist('docolor') docolor = 0;end
if ~exist('doentcalc');doentcalc = 1;end;
if ~doentcalc
  disp('loading entropies');
  load nipaper_01_ent;
end
if ~exist('doprint');doprint = 0;end
if docolor;
  colordef none;
else
  colordef white
end

filenames = ['ne2l0123no80z0ds2c0v3cn0ba0';
	     'ne3l0123no60z0ds2c0v3cn0ba0';
     	     'ne4l0123no48z0ds2c0v3cn0ba0';
       	     'ne5l0123no40z0ds2c0v3cn0ba0';];


fontsize = 12;
titlstrs = str2mat('(a) Q = 2','(b) Q = 3','(c) Q = 4','(d) Q = 5');
nf = size(filenames,1);
ind= 1;dind = 1;itype = 1;eind = 1;
nr = 2;nc = 2;
figure(1);clf;
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


if docolor;
ltype = ['m+';
	 'm+';
	 'm+';
	 'm+';
	 'm+';
	 'm+';
	 'yd';
	 'yo';];  
acolor = 'g';facecolor = 'y';
else
ltype = ['k+';
	 'k+';
	 'k+';
	 'k+';
	 'k+';
	 'k+';
	 'kd';
	 'ko';];
acolor = 'k';facecolor = 'k';
end



  for iorder = 1;
    

    
        % PICK BEST 10 RANDOM DESIGNS
	figure(1);
	[s,s_ind] = sort(squeeze(nreffmat(itype,iorder,dind,:,1)));
	s_ind = flipud(s_ind(:));
	subplot(nr,nc,ind);
	plot(squeeze(ndetmat(itype,iorder,dind,:,1)),squeeze(neffmat(itype,iorder,dind,:,1)),ltype(1,:),...
	     squeeze(ndetmat(itype,iorder,dind,:,2)),squeeze(neffmat(itype,iorder,dind,:,2)),ltype(2,:),...
	     squeeze(ndetmat(itype,iorder,dind,:,3)),squeeze(neffmat(itype,iorder,dind,:,3)),ltype(3,:),...	     
	     squeeze(ndetmat(itype,iorder,dind,:,4)),squeeze(neffmat(itype,iorder,dind,:,4)),ltype(4,:),...	     
	     squeeze(ndetmat(itype,iorder,dind,:,5)),squeeze(neffmat(itype, ...
						  iorder,dind,:,5)),ltype(5,:),...
	     squeeze(ndetmat(itype,iorder,dind,:,6)),squeeze(neffmat(itype, ...
						  iorder,dind,:,6)),ltype(6,:),...
	     squeeze(nrdetmat(itype,iorder,dind,s_ind(1:10),1)),squeeze(nreffmat(itype, ...
						  iorder,dind,s_ind(1:10),1)),ltype(7,:));
	axis([0 0.6 0 1.1]);
	hold on;
	if(nevents < 5)
	     hm = plot(squeeze(nmdetmat(itype,iorder,dind,1,1)),squeeze(nmeffmat(itype, ...
						  iorder,dind,1,1)),ltype(8,:)); ...

		  set(hm,'MarkerFaceColor',facecolor);
	end
	if docolor
	  plot(ntdet,nteffmat,'w-'); %theoretical curves
	else
	  plot(ntdet,nteffmat,'k-'); %theoretical curves
	end
	hold off;

	grid;
	set(gca,'FontSize',fontsize);
	ylabel('Normalized \xi_{tot}')
	xlabel('Normalized R_{tot}')
	title(deblank(titlstrs(ind,:)));
	if(ind == 3)
	  arrow([0.2,1],[0.08,.96],'Color',acolor);

%	  arrow([0.25,0.85],[0.14,0.83]);
      
	  arrow([0.3,.7],[0.08,0.8],'Color',acolor);
	  arrow([0.4,0.5],[0.3,0.35],'Color',acolor);
	  h1 = text(0.21,1.0,'m-sequence');set(h1,'FontSize',fontsize);
	  	  h1b = text(0.21,.93,'based');set(h1b,'FontSize',fontsize);
	  h3 = text(0.31,0.70,'randomly');set(h3,'FontSize', ...
							   fontsize);
	  h3b = text(0.31,0.63,'generated');set(h3b,'FontSize',fontsize);
	  h4 = text(0.41,0.5,'permuted');set(h4,'FontSize', ...
						      fontsize);
		  h4b = text(0.41,0.43,'block');set(h4b,'FontSize',fontsize);
	elseif (ind == 1);  
  	  arrow([0.2,1],[0.08,.98],'Color',acolor);

%	  arrow([0.25,0.85],[0.14,0.83]);
      
	  arrow([0.3,.7],[0.08,0.93],'Color',acolor);
	  arrow([0.4,0.5],[0.3,0.35],'Color',acolor);
	  h1 = text(0.21,1.0,'m-sequence');set(h1,'FontSize',fontsize);
	  	  h1b = text(0.21,.93,'based');set(h1b,'FontSize',fontsize);
	  h3 = text(0.31,0.70,'randomly');set(h3,'FontSize', ...
							   fontsize);
	  h3b = text(0.31,0.63,'generated');set(h3b,'FontSize',fontsize);
	  h4 = text(0.41,0.5,'permuted');set(h4,'FontSize', ...
						      fontsize);
		  h4b = text(0.41,0.43,'block');set(h4b,'FontSize', ...
							fontsize);
        elseif (ind == 2)		 
    	  arrow([0.2,1],[0.08,.96],'Color',acolor);

%	  arrow([0.25,0.85],[0.14,0.83]);
      
	  arrow([0.3,.7],[0.08,0.88],'Color',acolor);
	  arrow([0.4,0.5],[0.3,0.35],'Color',acolor);
	  h1 = text(0.21,1.0,'m-sequence');set(h1,'FontSize',fontsize);
	  	  h1b = text(0.21,.93,'based');set(h1b,'FontSize',fontsize);
	  h3 = text(0.31,0.70,'randomly');set(h3,'FontSize', ...
							   fontsize);
	  h3b = text(0.31,0.63,'generated');set(h3b,'FontSize',fontsize);
	  h4 = text(0.41,0.5,'permuted');set(h4,'FontSize', ...
						      fontsize);
		  h4b = text(0.41,0.43,'block');set(h4b,'FontSize', ...
							fontsize);
	elseif (ind == 4)		 
	  arrow([0.3,.7],[0.08,0.75],'Color',acolor);
	  arrow([0.4,0.5],[0.3,0.35],'Color',acolor);
	  h3 = text(0.31,0.70,'randomly');set(h3,'FontSize', ...
							   fontsize);
	  h3b = text(0.31,0.63,'generated');set(h3b,'FontSize',fontsize);
	  h4 = text(0.41,0.5,'permuted');set(h4,'FontSize', ...
						      fontsize);
		  h4b = text(0.41,0.43,'block');set(h4b,'FontSize', ...
							fontsize);
	end
	
	
	%calculate entropies of stimulus patterns
	nperm_mseq = size(mstimmat,2);
        if doentcalc		
	disp('calculating entropies');
	if nperm_mseq > 0
	nblocks = size(stimmat,3)-2;
	else
	nblocks = size(stimmat,3)-1;	  
	end
	nperm = size(stimmat,2);

	nperm_rand = 10;
	evec = NaN*ones(nperm,nblocks,maxentorder);
	revec =  NaN*ones(nperm_rand,maxentorder);
	if nperm_mseq > 0
	mevec =  NaN*ones(nperm_mseq,maxentorder);
	end

	for entorder = 1:maxentorder
	   revec(:,entorder) = calcentvec(squeeze(rstimmat(itype,s_ind(1:nperm_rand), ...
						  1,:)),entorder);
	   if nperm_mseq > 0
	   mevec(:,entorder) = calcentvec(squeeze(mstimmat(itype,:, ...
						  1,:)),entorder);
	   end
	  for blocknum = 1:nblocks
	    evec(:,blocknum,entorder) = calcentvec(squeeze(stimmat(itype,:, ...
						  blocknum,:)),entorder);
	  end
	end
	eval(sprintf('evec_%d = evec;',nfile));
	eval(sprintf('revec_%d = revec;',nfile));
	if nperm_mseq > 0
	eval(sprintf('mevec_%d = mevec;',nfile));
	end
        disp('done with entropies');
	else

	  eval(sprintf('evec = evec_%d;',nfile));
	  eval(sprintf('revec = revec_%d;',nfile));
	  if nperm_mseq > 0
	    eval(sprintf('mevec = mevec_%d;',nfile));
	  end
	end
	
	fignum = 2;
	for entorder = 1:maxentorder
	  figure(fignum);
	  subplot(nr,nc,eind);
	  if nperm_mseq > 0
	  plot(squeeze(neffmat(itype,iorder,dind,:,1:6)),squeeze(2.^evec(:,:,entorder)-1),'-o',...
	     squeeze(nreffmat(itype,iorder,dind,s_ind(1:10),1)),2.^revec(:,entorder)-1,'-x',...
	     squeeze(nmeffmat(itype,iorder,dind,:,1)),2.^mevec(:, ...
						  entorder)-1,'-d');
	  else
	    plot(squeeze(neffmat(itype,iorder,dind,:,1:6)),squeeze(2.^evec(:,:,entorder)-1),'-o',...
	     squeeze(nreffmat(itype,iorder,dind,s_ind(1:10),1)),2.^revec(:,entorder)-1,'-x');

	  end
	       fignum = fignum + 1;
	end

  end

%keyboard;  
  
ind = ind + 1;
eind = eind + 1;
end

if doprint
if docolor
  set(gcf,'InvertHardCopy','off');
  print -dtiff nimg01_ent.tiff
else
  print -deps nimg01_ent.eps
end
end


save nipaper_01_ent evec_1 evec_2 evec_3 evec_4 mevec_1 mevec_2 mevec_3 ...
       revec_1 revec_2 revec_3 revec_4
