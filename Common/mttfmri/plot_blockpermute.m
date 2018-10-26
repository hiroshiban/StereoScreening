% example of block permutation algorithm
map = [0 1 0;1 0 0;0 0 1];
colormap(map);
load blockpermute; % load in the same random seed each time
if ~exist('doprint');doprint = 0;end;
nevents = 2;
numblocks = 2;
numones = 30;
npts = 90;

sizeofblock = round(numones/numblocks);
blockspacing = npts/(numblocks);
baseveclen = npts;
basevec = zeros(baseveclen,1);
for k = 1:nevents
  thisbasevec =zeros(baseveclen,1);
  startlocation = 1 + (k-1)*sizeofblock;
  blocklocs = startlocation:blockspacing:npts;
  thisbasevec(blocklocs) =1;
  thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
  basevec(find(thisbasevec)) = k;
end

paths = 11;
rand('state',r2);
for pnum = 1:paths;
  if(pnum > 1)
  looping = 1;
  while looping % loop until we find two points with different values
    p1 = unidrnd(npts);
    p2 = unidrnd(npts);
    if basevec(p1) ~= basevec(p2)
      tmp = basevec(p1);
      basevec(p1)=basevec(p2);
      basevec(p2)=tmp;
      looping = 0;
    end
  end
  end
  mstimmat(:,pnum) = basevec;
end
[x,y] = stimpatch(mstimmat(:,1));
figure(1);clf;
axis([0 90 2.5 11.5]);
patch(x,y+10,mstimmat(:,1)');
hold on;
patch(x,y+8, mstimmat(:,2)');

inds = find((mstimmat(:,1)-mstimmat(:,2))~=0);
inds = inds-0.5;
arrow([inds(1) 10],[inds(2) 9],'length',10);
arrow([inds(2) 10],[inds(1) 9],'length',10);

hold on;
patch(x,y+6, mstimmat(:,3)');
inds = find((mstimmat(:,3)-mstimmat(:,2))~=0);
inds = inds-0.5;
arrow([inds(1) 8],[inds(2) 7],'length',10);
arrow([inds(2) 8],[inds(1) 7],'length',10);

hold on;
patch(x,y+3, mstimmat(:,11)');
hc = plot([45 45 45],[5.5 5 4.5],'ko');
set(hc,'MarkerFaceColor','k','MarkerSize',8);
fsize = 12;xpos = 75;
ht = text(xpos,11.25,'Block Design');
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

arrow([10,2.7],[80,2.7],'length',10);
ht = text(40,2.5,'Time');set(ht,'FontSize',fsize);
if doprint
print -dtiff NIMG771_fig1.tif
print -depsc NIMG771_fig1.eps
end

