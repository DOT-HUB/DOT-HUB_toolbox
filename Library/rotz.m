function out = rotz(in,ang,about)

%in = Nx3 points in 3D
%ang = angle in degrees
%about = centre point of rotation (1x3)
R = [cosd(ang) -sind(ang) 0; sind(ang) cosd(ang) 0; 0 0 1];

if ~exist('about','var')  
    out = (R*in')';
else
    tmpIn = in - repmat(about,size(in,1),1);
    tmpRot = (R*tmpIn')';
    out = tmpRot + repmat(about,size(in,1),1);
end