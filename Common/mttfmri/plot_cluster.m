% Plot out results from multiple realizations of
% permuted block and clustered m-sequences generated
% by batch_cluster
%
% Uses the function compress_clump to reduce the number
% of datapoints plotted. This is necessary when the number of
% realizations is large (e.g. npaths = 1000). 
% Note that the loading of the files and compressing can take
% a while if the number of npaths is large.
% The number of npaths is currently hardwired into the name of the file.
% In this example version of the code, npaths = 100 (in the paper it was 1000).
%
% To analyze data:
%   clear;nevents = 2;plot_cluster
%   clear;nevents = 4;plot_cluster
%
% To make replotting easier, results needed to plot are stored in
% a clump*mat. To replot a previous result, type
%
%  clear;loaddata = 0;calcdata = 0;nevents = 2;loadprevious = 1;plot_cluster
%  clear;loaddata = 0;calcdata = 0;nevents = 4;loadprevious = 1;plot_cluster
%
%
%  030513 -- correcting what looks like a bug in ploting of 2 block data.
%         -- adding stuff for entropy calcs
%
%  030829 -- adding comments and printing of plots
%         -- adding more stuff to save, so we don't need to reload and
%            recalc just to redo plots. 



if ~exist('doprint');doprint = 0;end;
if ~exist('do100');do100 = 0;end;
if ~exist('nevents');nevents = 4;end;
if ~exist('loaddata');loaddata = 1;end
if ~exist('calcdata');calcdata = 1;end
if ~exist('loadprevious');loadprevious = 0;end;

maxentorder = 2;
iorder = 1;
  mode = 1;    
%always LOAD mixed designs
if nevents == 4  
  load mixed_ne4np240;
else  
  load mixed_ne2np240;
end
if loaddata

  if nevents == 4
  % LOAD CLUSTERED M-SEQUENCES
%  load cl4l0npe30np1000
  load cl4l0npe30np100
  ndet = cldet/tmaxdet;
  neff = cleff/tmaxeff;
  

  
  % LOAD PERMUTED 2-BLOCK SEQUENCES
%  load bl4l0npe100np1000nb2.mat
  load bl4l0npe100np100nb2.mat
  nbdet = bldet/tmaxdet;
  nbeff = bleff/tmaxeff;
  
  
  else
    
  % LOAD CLUSTERED M-SEQUENCES
%  load cl2l0npe30np1000
  load cl2l0npe30np100
  ndet = cldet/tmaxdet;
  neff = cleff/tmaxeff;
  
  
  % LOAD PERMUTED 2-BLOCK SEQUENCES
%  load bl2l0npe100np1000nb2.mat
  load bl2l0npe100np100nb2.mat
  nbdet = bldet/tmaxdet;
  nbeff = bleff/tmaxeff;
  

    
    
  end
  
  if do100
    ndet = ndet(:,:,:,1:100,:);
    neff = neff(:,:,:,1:100,:);
    nbdet = nbdet(:,:,:,1:100,:);
    nbeff = nbeff(:,:,:,1:100,:);
  end

end  


if calcdata
  
    % COMPRESS THE DATA FOR PLOTTING

  [ne,nd,pind] = compress_clump(squeeze(neff(mode,iorder,:,:,1)),...
			    squeeze(ndet(mode,iorder,:,:,1)));
   
  isNumber = find(~isnan(ne(:,1)));
  ne = ne(isNumber,:);
  nd = nd(isNumber,:);
  pind = pind(isNumber,:);
  cl_ne = ne;
  cl_nd = nd;
  %calculate entropies
  if ~exist('m_ent')
   nent = size(ne,1);
   clstim_touse = squeeze(clstim(mode,:,:,:));
   n_cperm = size(clstim,2);
   n_cpaths = size(clstim,3);
   n_cpts = size(clstim,4);
   clstim_touse = reshape(clstim_touse,n_cperm*n_cpaths,n_cpts);
   clstim_touse = clstim_touse(pind(:,1),:); %stimuli for max det power

   m_ent = NaN*ones(nent,maxentorder);
    for entorder = 1:maxentorder
     disp('calculating clumped entropies')
     m_ent(:,entorder) = calcentvec(clstim_touse,entorder);
     disp('done calculating clumped entropies')
     end
  end
  
    
  %COMPRESS
  mode = 1;
    [ne,nd,pind] = compress_clump(squeeze(nbeff(mode,iorder,:,:,1)),...
			    squeeze(nbdet(mode,iorder,:,:,1)));
    isNumber = find(~isnan(ne(:,1)));    
    ne = ne(isNumber,:);
    nd = nd(isNumber,:);
    pind = pind(isNumber,:);
    bl_ne = ne;
    bl_nd = nd;
    %Calculation of entropy
    if ~exist('bl_ent');
       nent = size(ne,1);
       blstim_touse = squeeze(blstim(mode,:,:,:));
       n_blperm = size(blstim,2);
       n_blpaths = size(blstim,3);
       n_blpts = size(blstim,4);
       blstim_touse = reshape(blstim_touse,n_blperm*n_blpaths,n_blpts);
       blstim_touse = blstim_touse(pind(:,1),:);
       % this next thing is a kludge because of an error in some versions
       % of tryclumpgen.m
       ix = find(isnan(blstim_touse));
       blstim_touse(ix) = 0;
       
       bl_ent = NaN*ones(nent,maxentorder);
       for entorder = 1:maxentorder
	 disp('calculating block permuted entropies')
	 bl_ent(:,entorder) = calcentvec(blstim_touse,entorder);
	 disp('done calculating blockpermuted entropies')
       end
    end

meff = squeeze(neff(1,1,1,1,1));
mdet = squeeze(ndet(1,1,1,1,1));
blockeff = squeeze(nbeff(1,1,1,1,1));
blockdet = squeeze(nbdet(1,1,1,1,1));  
mseq_stim = squeeze(clstim(1,1,1,:));
m_ent2 = calcent(clstim(1,1,1,:),2);  
  
end


if ~exist('loadprevious') loadprevious = 0;end;
if loadprevious;
filename = sprintf('clump_ne%d_do%d',nevents,do100);
 eval(sprintf('load %s',filename));
end


%SAVE CALCULATED RESULTS FOR FUTURE USE
if (loadprevious == 0) 
filename = sprintf('clump_ne%d_do%d',nevents,do100);
if do100
  cl_nd100 = cl_nd;
  cl_ne100 = cl_ne;
  bl_ne100 = bl_ne;
  bl_nd100 = bl_nd;
  m_ent100 = m_ent;
  bl_ent100 = bl_ent;
  m_ent2_100 = m_ent2;
  clstim_touse100 = clstim_touse;
  varnames = [' cl_nd100 cl_ne100 bl_ne100 bl_nd100 m_ent100 mdet meff blockeff',...
	      ' blockdet m_ent2_100 bl_ent100 clstim_touse100',...
	       ' mseq_stim'];
	        eval(sprintf('save %s %s',filename,varnames));

else
  
  varnames = [' cl_nd cl_ne bl_ne bl_nd m_ent mdet meff blockeff',...
	      ' blockdet m_ent2 bl_ent clstim_touse blstim_touse mseq_stim' ...
	      ' tmaxdet tmaxeff '];
  eval(sprintf('save %s %s',filename,varnames));
end 
end

%PLOT THE DATA
span = 1:3;figure(1);clf;
subplot('position',[0.07 0.08 0.55 0.85]);
cspan = 3:5:147;
bspan = 3:5:148;
%m-sequence
ht = plot(mdet,meff,'ko',cl_nd(cspan,span),cl_ne(cspan,span),'b+',bl_nd(bspan,span),bl_ne(bspan,span),'r+',...
     blockdet,blockeff,'ms');
set(ht(8),'MarkerFaceColor','m','MarkerSize',9);
set(ht(1),'MarkerFaceColor','k','MarkerSize',9);



% PLOT THEORETICAL CURVES
tcolor = 'k';
nummods = 15;npts = 240;Q = nevents;
k = nummods;
i0 = 1.0;
alpha =[1/k 0.1:0.05:1.0];
if nevents == 4
thetavec = [45 47.5 90]/180*pi;
else
thetavec = [45 90]/180*pi;
end
teff = alpha.*(1-alpha) ./(1+alpha*(k^2-2*k))*i0*nummods*npts/(2*(Q+1));
alphamat = alpha(:)*ones(1,length(thetavec));
teffmat = teff(:)*ones(1,length(thetavec));
xrangemat = ones(length(alpha),1)*thetavec(:)';

tdet =i0*nummods*((1-alphamat).*(sin(xrangemat).^2)/(k-1)+alphamat.* ...
	  cos(xrangemat).^2);
tdet = tdet*npts/(2*(Q+1));
ntdet = tdet/tmaxdet;
nteffmat = teffmat/tmaxeff;hold on;
if nevents == 4;
htheory = plot(ntdet(:,1),nteffmat(:,1),'k-',...
     ntdet(:,2),nteffmat(:,2),'k--',...
     ntdet(:,3),nteffmat(:,3),'k-.');
hold off;
else
htheory = plot(ntdet(:,1),nteffmat(:,1),'k-',...
     ntdet(:,2),nteffmat(:,2),'k-.');
hold off;
end
axis([0 0.5 0 1.02]);


%PLOT MIXED DESIGNS
hold on;
hp = plot(ndetmix(:,1),neffmix(:,1),'go',ndetmix(:,2),neffmix(:,2),'co');
hold off;


fontsize = 12;
set(gca,'FontSize',fontsize);
ylabel('Normalized \xi_{tot}')
xlabel('Normalized R_{tot}')
grid;


% ADD LEGENDS
if nevents == 4
hl = legend([ht([1 2:3:7 8]);hp;htheory],'m-sequence','Clustered m-seq',...
	     'Permuted 2-block','2-block','Mixed Design 1-block',...
	    'Mixed Design 2-block','\theta =45 degrees','\theta = 47.5 degrees',...
	    '\theta = 90 degrees');

else
hl = legend([ht([1 2:3:7 8]);hp;htheory],'m-sequence','Clustered m-seq',...
	     'Permuted 2-block','2-block','Mixed Design 1-block',...
	    'Mixed Design 2-block','\theta =45 degrees',...
	    '\theta = 90 degrees');
end
title(sprintf('(a) Efficiency and Power for %d trial types',nevents));



% find something with twice the detection power of m-seq
if 0
n0 = 2*ndet(1,1,1,1,1);
ndet0 = squeeze(ndet(1,1,:,:,1));
neff0 = squeeze(neff(1,1,:,:,1));
zz = find(abs(ndet0-n0) < .005);
[zmax,zmaxi] = max(neff0(zz));
bestmat = [];
bestmat = [bestmat;neff0(zz(zmaxi)) ndet0(zz(zmaxi))];

ndet0 = squeeze(nbdet(1,1,:,:,1));
neff0 = squeeze(nbeff(1,1,:,:,1));
zz = find(abs(ndet0-n0) < .005);
[zmax,zmaxi] = max(neff0(zz));
bestmat = [bestmat;neff0(zz(zmaxi)) ndet0(zz(zmaxi))];
end

% PLOT ENTROPIES
subplot('position',[0.7 0.08 0.25 0.85]);



span = 1;
ht = plot(2^m_ent2,meff,'k-o',2.^m_ent(cspan,2),cl_ne(cspan,span),'b+',2.^bl_ent(bspan,2),bl_ne(bspan,span),'r+',...
     2^bl_ent(1,2),blockeff,'ms',...
     2.^squeeze(entmix(:,1,2)),neffmix(:,1),'go',...
	  2.^squeeze(entmix(:,2,2)),neffmix(:,2),'co');

set(ht(4),'MarkerFaceColor','m','MarkerSize',9);
set(ht(1),'MarkerFaceColor','k','MarkerSize',9);
grid;
set(gca,'FontSize',fontsize);
xlabel('2\^H_2');
ylabel('Normalized \xi_{tot}');
axis([1 nevents+1 0 1.02])
title('(b) Efficiency vs. 2\^(Entropy)');


% PRINT OUT THE FILES
if doprint
if nevents == 4
  print -dtiff NIMG771_fig4.tif
  print -depsc NIMG771_fig4.eps
  
else
  print -dtiff NIMG771_fig3.tif
  print -depsc NIMG771_fig3.eps
  
end
end

% MAKE UP DIAGRAMS OF STIMULUS PATTERNS
figure(2);clf;fsize = 11;
hstim = subplot('position',[0.08 .5 .85 .45]);

if nevents == 4
map = [0 1 0;1 0 0;0 0 1;0 1 1;1 1 0];  
yscale = 1;
else
map = [0 1 0;1 0 0;0 0 1];
yscale = 0.6;
end  
colormap(map);
target = 0.81;
axis([0 240 0 4.3]);

%mixed design

[m,mind] = min(abs(neffmix(:,1)-target));
this_mixedstim = squeeze(mixedstim(mind,:,1));

[x,y] = stimpatch(this_mixedstim);
patch(x,y*yscale,this_mixedstim);

	      
if nevents == 2;prefix = '(d)';else;prefix = 'c';end;
textstr = ...
[prefix, ' Mixed Design   \xi_{tot} = ',...
 sprintf(' %.2f  R_{tot} = %.2f  H_2 = %.2f bits  ',...
	neffmix(mind,1),ndetmix(mind,1),entmix(mind,1,2)),...
		   ' 2\^H_2 = ',num2str(2^(entmix(mind,1,2)),3)];
if 1

ht = text(0.1,1*yscale + .1,textstr);
set(ht,'FontSize',fsize);	      
else
title(textstr);  
end
%plot(xax,this_mixedstim,':ks');set(gca,'Ylim',[0 4.5]);

% clustered m-sequence
[m,mind] = min(abs(cl_ne(:,1)-target));
this_cstim = squeeze(clstim_touse(mind,:));
yincbase = 1.9;
if nevents == 2;yinc = 2*yincbase;else;yinc = 1.5;end
patch(x,yscale*(y+yinc),this_cstim);
set(gca,'Visible','off');
if 1
ht = text(0.1,yscale*(1+yinc)+.1,...
['(b) Clustered m-sequence   \xi_{tot} = ',...
 sprintf(' %.2f  R_{tot} = %.2f  H_2 = %.2f bits ',...
	cl_ne(mind,1),cl_nd(mind,1),m_ent(mind,2)),...
	 ' 2\^H_2 = ',num2str(2^(m_ent(mind,2)),3)]);	 
	      set(ht,'FontSize',fsize);	      
end	      

% permuted block sequence
if nevents == 2;
[m,mind] = min(abs(bl_ne(:,1)-target));  
this_blstim = squeeze(blstim_touse(mind,:));
yinc = 1*yincbase;
patch(x,yscale*(y+yinc),this_blstim(:)');
ht = text(0.1,yscale*(1+yinc)+.1,...
['(c) Permuted 2-block   \xi_{tot} = ',...
 sprintf(' %.2f  R_{tot} = %.2f  H_2 = %.2f bits ',...
	bl_ne(mind,1),bl_nd(mind,1),bl_ent(mind,2)),...
	 ' 2\^H_2 = ',num2str(2^(bl_ent(mind,2)),3)])
	      set(ht,'FontSize',fsize);	      
end




% m-sequence
if nevents == 2;yinc = 3*yincbase;else;yinc = 3;end
patch(x,yscale*(y+yinc),mseq_stim(:)');
hold on;ypos = 4.15;
if nevents == 2
  ypos = 4.2;xadd = 20;
hl = plot([140 144]+xadd,[ypos ypos],'r-');set(hl,'LineWidth',3);
ht = text(146+xadd,ypos,'A');set(ht,'FontSize',fsize);
hl = plot([155 159]+xadd,[ypos ypos],'b-');set(hl,'LineWidth',3);
ht = text(161+xadd,ypos,'B');set(ht,'FontSize',fsize);
hl = plot([168 172]+xadd,[ypos ypos],'g-');set(hl,'LineWidth',3);
ht = text(174+xadd,ypos,'Null');set(ht,'FontSize',fsize);
else
xpos = 160;  
hl = plot([xpos xpos+4],[ypos ypos],'r-');set(hl,'LineWidth',3);
xpos = xpos + 6;
ht = text(xpos,ypos,'A');set(ht,'FontSize',fsize);
xpos = xpos + 9;
hl = plot([xpos xpos+4],[ypos ypos],'b-');set(hl,'LineWidth',3);
xpos = xpos + 6;
ht = text(xpos,ypos,'B');set(ht,'FontSize',fsize);
xpos = xpos + 9;
hl = plot([xpos xpos+4],[ypos ypos],'c-');set(hl,'LineWidth',3);
xpos = xpos + 6;
ht = text(xpos,ypos,'C');set(ht,'FontSize',fsize);  
xpos = xpos+9;
hl = plot([xpos xpos+4],[ypos ypos],'y-');set(hl,'LineWidth',3);
xpos = xpos + 6;
ht = text(xpos,ypos,'D');set(ht,'FontSize',fsize);  
xpos = xpos + 9;
hl = plot([xpos xpos+4],[ypos ypos],'g-');set(hl,'LineWidth',3);
xpos = xpos + 6;
ht = text(xpos,ypos,'Null');set(ht,'FontSize',fsize);  
  
end
if 1
ht = text(0.1,yscale*(1+yinc)+.1,...
['(a) m-sequence   \xi_{tot} = ',...
 sprintf(' %.2f  R_{tot} = %.2f  H_2 = %.2f bits ',...
	 meff,mdet,m_ent2),...
	 ' 2\^H_2 = ',num2str(2^m_ent2,3)]);

	      set(ht,'FontSize',fsize);	      
end;	      
set(gca,'Visible','off');


	      
	      

% A MEASURE OF RUNNING ENTROPY
if ~exist('do_calcent');do_calcent = 1;end;
if do_calcent
winwidth = 40;
nshifts = npts-winwidth;
mseq_mat = rcircmat(mseq_stim(:),nshifts);
cstim_mat = rcircmat(this_cstim(:),nshifts);

mix_mat = rcircmat(this_mixedstim(:),nshifts);
mseq_mat = mseq_mat(1:winwidth,:);
cstim_mat = cstim_mat(1:winwidth,:);
mix_mat = mix_mat(1:winwidth,:);
eorder = 2;
mseq_ment = calcentvec(mseq_mat',eorder);
cstim_ment = calcentvec(cstim_mat',eorder);
mix_ment = calcentvec(mix_mat',eorder);
eorder = 1;
mseq_ment1 = calcentvec(mseq_mat',eorder);
cstim_ment1 = calcentvec(cstim_mat',eorder);
mix_ment1 = calcentvec(mix_mat',eorder);

if nevents == 2
  blstim_mat = rcircmat(this_blstim(:),nshifts);
  blstim_mat = blstim_mat(1:winwidth,:);
  bl_ment1 = calcentvec(blstim_mat',1);
  bl_ment = calcentvec(blstim_mat',2);
end
end


subplot('position',[0.08 .27 .85 .18]);
xax =  (0:(nshifts-1))+ winwidth/2;
if nevents ==2
plot(xax,2.^mseq_ment1,'g',xax,2.^cstim_ment1,'r',xax,2.^mix_ment1,'b',...
     xax,2.^bl_ment1,'k');  
else  
plot(xax,2.^mseq_ment1,'g',xax,2.^cstim_ment1,'r',xax,2.^mix_ment1,'b');
end
if nevents == 2
axis([0 240 1 3]);
else
 axis([0 240 1 5]);  

end
if nevents == 2;ypos = 3.2;prefix = '(e)';else; ypos = 5.4;prefix = '(d)';end;

ht = text(.1,ypos,[prefix,' 1st order local conditional entropy']);
set(ht,'Fontsize',fsize');
ylabel('2\^H_2');
set(gca,'XTickLabel',[]);
if nevents == 2
  hleg = legend('m-seq','clustered m-seq','mixed','permuted 2-block',3);
  set(hleg,'Visible','off');set(allchild(hleg),'visible','on');
end


grid;
subplot('position',[0.08 .05 .85 .18]);
if nevents == 2
plot(xax,2.^mseq_ment,'g',xax,2.^cstim_ment,'r',xax,2.^mix_ment,'b',...
     xax,2.^bl_ment,'k');  
prefix = '(f)';
else
plot(xax,2.^mseq_ment,'g',xax,2.^cstim_ment,'r',xax,2.^mix_ment,'b');
prefix = '(e)';
end

ht = text(110,0.7,'Seconds');set(ht,'Fontsize',fsize');
ht = text(.1,ypos,[prefix,' 2nd order local conditional entropy']);
set(ht,'Fontsize',fsize');
if nevents == 2
axis([0 240 1 3]);
else
  axis([0 240 1 5]);
end
ylabel('2\^H_1');
grid;
if nevents == 4
hleg = legend('m-seq','clustered m-seq','mixed',2);
set(hleg,'Visible','off');set(allchild(hleg),'visible','on');
end



% PRINT OUT THE FILES

if doprint
if nevents == 4
  print -dtiff NIMG771_fig8.tif
  print -depsc NIMG771_fig8.eps
  
else
  print -dtiff NIMG771_fig7.tif
  print -depsc NIMG771_fig7.eps
  
end
end