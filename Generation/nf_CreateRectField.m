function field=nf_CreateRectField(outer_fieldSize,inner_fieldSize,inner_height,pix_per_deg,fine_coefficient)

% function field=nf_CreateRectField(outer_fieldSize,inner_fieldSize,inner_height,pix_per_deg,fine_coefficient)
%
% Creates rectangular height field
% 
% [input]
% outer_fieldSize  : the size of the field in degrees, [row,col] (deg)
% inner_fieldSize  : the size of the field in degrees, [row,col] (deg)
% inner_height     : plane height in pix, [val]
%                    if the value is minus, near surface is generated.
% pix_per_deg : pixels per degree, [pixels]
% fine_coefficient : (optional) if larger, the generated oval become finer,
%                    but comsumes much CPU power. [val]
%                    (default=1, as is, no tuning)
%
% [output]
% field       : rectangular height field image, double format, [row,col]
% 
% !!! NOTICE !!!
% for speeding up image generation process,
% I will not put codes for nargin checks.
% Please be careful.
%
% Created    : "2010-06-14 12:20:56 ban"
% Last Update: "2010-08-01 14:32:33 ban"

% convert to pixels and radians
outer_fieldSize=round(outer_fieldSize.*pix_per_deg).*[1,fine_coefficient];
if mod(outer_fieldSize(1),2), outer_fieldSize(1)=outer_fieldSize(1)-1; end
if mod(outer_fieldSize(2),2), outer_fieldSize(2)=outer_fieldSize(2)-1; end

inner_fieldSize=round(inner_fieldSize.*pix_per_deg).*[1,fine_coefficient];
if mod(inner_fieldSize(1),2), inner_fieldSize(1)=inner_fieldSize(1)-1; end
if mod(inner_fieldSize(2),2), inner_fieldSize(2)=inner_fieldSize(2)-1; end

% create rectangular field
field=zeros(outer_fieldSize);

rowidx=(size(field,1)-inner_fieldSize(1))/2+1:(size(field,1)+inner_fieldSize(1))/2;
colidx=(size(field,2)-inner_fieldSize(2))/2+1:(size(field,2)+inner_fieldSize(2))/2;

field(rowidx,colidx)=inner_height;

return
