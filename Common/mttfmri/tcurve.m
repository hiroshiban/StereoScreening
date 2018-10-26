function [tdet,teffmat,tmaxdet,tmaxeff] = tcurve(Q,nummods,npts,theta);

% Make up theoretical trade-off curves
%  
%  [tdet,teffmat,tmaxdet,tmaxeff] = tcurve(Q,nummods,npts,theta);
%  
%  Inputs:
%    Q       = number of non-Null trial types  
%    nummods = number of model functions (e.g. points in HRF).
%    npts    = number of points in the design  
%    theta   = angles to plot over, default = [45 55 65 80 90];
%
%  Outputs:
%    tdet   :   detection power
%    teffmat:   efficiency
%    tmaxdet:   theoretical maximum detection power  
%    tmaxeff:   theoretical maximum efficiency
%      
if ~exist('theta');theta =  [45  55 65  80 90];end;

k = nummods;
i0 = 1.0;
alpha =[1/k 0.1:0.05:1.0];
thetavec = theta/180*pi;
teff = alpha.*(1-alpha) ./(1+alpha*(k^2-2*k))*i0*nummods*npts/(2*(Q+1));
alphamat = alpha(:)*ones(1,length(thetavec));
teffmat = teff(:)*ones(1,length(thetavec));
xrangemat = ones(length(alpha),1)*thetavec(:)';

tdet =i0*nummods*((1-alphamat).*(sin(xrangemat).^2)/(k-1)+alphamat.* ...
	  cos(xrangemat).^2);
tdet = tdet*npts/(2*(Q+1));

tmaxeff = npts/(2*(Q+1)*nummods);
tmaxdet = npts*nummods/(2*(Q+1));
ntdet = tdet/tmaxdet;
nteffmat = teffmat/tmaxeff;

