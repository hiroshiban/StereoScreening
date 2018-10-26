% plot example of clustered m-sequence algorithm
if ~exist('doprint');doprint = 0;end;
load clusterstate;% load same random state ever time
rand('state',r1); 
nevents = 2;
base = 3;
power = 4;
 mparam.base = base;
 mparam.power = power;
 mparam.shift = 0;
 len = base^power-1;
msmat = gen_mseq(mparam,len);
mstim = msmat(:,1);
map = [0 1 0;1 0 0;0 0 1];
colormap(map);



paths = 11;
mstimmat = NaN*ones(len,paths);
mstimmat(:,1) = mstim(:);
ind = 1;
[x,y] = stimpatch(mstimmat(:,ind));
thistype = 1;

for k = 1:(paths-1);
  mstimmat(:,ind+1) = clumpvec(mstimmat(:,ind), ...
						    thistype);

  thistype = thistype + 1;
  if thistype > nevents
    thistype = 1;
  end

  ind = ind + 1;
end

figure(1);clf;
axis([0 80 2.5 11.5]);
patch(x,y+10,mstimmat(:,1)');
hold on;
patch(x,y+8, mstimmat(:,2)');

inds = find((mstimmat(:,1)-mstimmat(:,2))~=0);
inds = inds-0.5;
arrow([inds(1)-0.5 10],[inds(2) 9],'length',10);
arrow([inds(2)+0.5 10],[inds(1) 9],'length',10);

hold on;
patch(x,y+6, mstimmat(:,3)');
inds = find((mstimmat(:,3)-mstimmat(:,2))~=0);
inds = inds-0.5;
arrow([inds(1)-0.5 8],[inds(2) 7],'length',10);
arrow([inds(2)+0.5 8],[inds(1) 7],'length',10);

hold on;
patch(x,y+3, mstimmat(:,11)');
hc = plot([40 40 40],[5.5 5 4.5],'ko');
set(hc,'MarkerFaceColor','k','MarkerSize',8);
fsize = 12;xpos = 65;
ht = text(xpos,11.25,'m-sequence');
set(ht,'Fontsize',fsize);
ht = text(xpos,9.25,'Iteration 1');
set(ht,'Fontsize',fsize);
ht = text(xpos,7.25,'Iteration 2');
set(ht,'Fontsize',fsize);
ht = text(xpos,4.25,'Iteration 10');
set(ht,'Fontsize',fsize);

set(gca,'Visible','off');

hl = plot([1 5],[11.3 11.3],'r');
set(hl,'LineWidth',5);
ht = text(7,11.3,'A');set(ht,'FontSize',fsize);

hl = plot(4+[7 11],[11.3 11.3],'b');
set(hl,'LineWidth',5);
ht = text(4+12,11.3,'B');set(ht,'FontSize',fsize);

hl = plot(13+[7 11],[11.3 11.3],'g');
set(hl,'LineWidth',5);
ht = text(13+12,11.3,'Null');set(ht,'FontSize',fsize);


arrow([5,2.7],[75,2.7],'length',10);
ht = text(40,2.5,'Time');set(ht,'FontSize',fsize);
if doprint
print -dtiff NIMG771_fig6.tif
print -depsc NIMG771_fig6.eps
end
