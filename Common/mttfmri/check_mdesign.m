function [eff,det,ueff,udet] = check_mdesign(datname,nummods,nevents,PS,h,V,uC,mode);

% [eff,det,ueff,udet] = ...
% check_mdesign(datname,nummods,nevents,PS,h,V,uC,mode);
%
%  utility function for loading in an ASCII file and computing efficiency
%  and power. 

eval(sprintf('load %s',datname));
disp(sprintf('loaded %s',datname))
eval(sprintf('stimpattern = %s;',datname(1:(length(datname)-4))));
stimpattern = stimpattern(:);

[thiseff,thisdet,this_ueff,this_udet] = ...
    calc_meffdet(stimpattern,nummods,nevents,PS,h,V,uC,mode);
eff = thiseff(1);det = thisdet(1);
ueff = this_ueff(1);udet = this_udet(1);