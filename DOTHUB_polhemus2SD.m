function [SD_3D] = DOTHUB_polhemus2SD(polhemusFilename,SD)

%This script takes the polhemus data associated with all optode points and 
%the existing (2D) SD structure and creates an SD_3D 

% INPUTS: ################################################################
% polhemusFilename =    the .csv file name This data is a four-column CSV of
%                       position label, then x, y, z coordinate (assumed in cm).
%                       The first five rows should be Nz, Iz, Ar, Al, Cz. The
%                       following rows should be Src1, Src2, Src3 ... then
%                       Det1, Det2, Det3...

% OUTPUTS: ###############################################################
% SD_3D = The SD file containing the 3D information associated with polhemus
% data. Also contains the a 'Landmarks' variable that saves the cranial
% landmarks (Nz Iz Ar Al Cz). Also save out into an .SD file with a name
% matching the input .csv appended with _3D.

% RJC, UCL, April 2020

% UPDATES ################################################################

% TO DO LIST #############################################################
% WRITE THIS FUNCTION!

%Manage Inputs ###########################################################
%#########################################################################
if ~exist('polhemusFilename','var')
    [file,path] = uigetfile('*.csv','Select Polhemus data set (.csv)','MultiSelect','on');
    polhemusFilename = fullfile(path,file);
end

%Target SD output names
outname3D = [polhemusFilename(1:end-4) '_3D.SD'];

%Wavelengths
wavelength1 = SD.Lambda(1);
wavelength2 = SD.Lambda(2);

%#########################################################################
%#########################################################################
%Load data
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%f%f%[^\n\r]';
fileID = fopen(polhemusFilename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
allPos = [dataArray{2} dataArray{3} dataArray{4}]*10; %Polhemus output in cm, convert to mm.%#################################################%#################################################%#################################################%#################################################

%Crop out landmarks and polhemus measurement points
landmarks = allPos(1:5,:);
PolPoints = allPos(6:end,:);

%Translate and Rotate so that Iz is at 0 0 0, Nz is at 0 y 0, Ar and Al have same z
%value and Cz is on top

Iz = landmarks(2,:);
landmarks = landmarks - repmat(Iz,size(landmarks,1),1);
PolPoints = PolPoints - repmat(Iz,size(PolPoints,1),1);
Cz = landmarks(5,:);
%Rotate Cz to be in +ve z;
[azimuth,elevation,r] = cart2sph(Cz(1),Cz(2),Cz(3));
landmarks = rotz(landmarks,90-rad2deg(azimuth),[0 0 0]);
PolPoints = rotz(PolPoints,90-rad2deg(azimuth),[0 0 0]);
landmarks = rotx(landmarks,-rad2deg(elevation),[0 0 0]);
PolPoints = rotx(PolPoints,-rad2deg(elevation),[0 0 0]);

Nz = landmarks(1,:);
Iz = landmarks(2,:);
Cz = landmarks(5,:);
[azimuth,elevation,r] = cart2sph(Nz(1)-Iz(1),Nz(2)-Iz(2),Nz(3)-Iz(3));
landmarks = rotz(landmarks,90-rad2deg(azimuth),[0 0 0]);
landmarks = rotx(landmarks,-rad2deg(elevation),[0 0 0]);
PolPoints = rotz(PolPoints,90-rad2deg(azimuth),[0 0 0]);
PolPoints = rotx(PolPoints,-rad2deg(elevation),[0 0 0]);

%ADD ROTATION FOR Ar Al
%Rotate so Ar and Al have same z
Ar = landmarks(3,:);
Al = landmarks(4,:);
[azimuth,elevation,r] = cart2sph(Ar(1)-Al(1),Ar(2)-Al(2),Ar(3)-Al(3));
landmarks = roty(landmarks,rad2deg(elevation),[0 0 0]);
PolPoints = roty(PolPoints,rad2deg(elevation),[0 0 0]);

%Initial plot to check data;
f1 = figure;
set(f1,'Name','Polhemus Data');
subplot(1,2,1);
plotmesh(landmarks,'g.','MarkerSize', 30);hold on;
landmarkLabels = {'Nz','Iz','Ar','Al','Cz'};
for i = 1:size(landmarks)
    text(landmarks(i,1),landmarks(i,2)+3,landmarks(i,3)+3,landmarkLabels{i});
end
for i = 1:size(PolPoints,1)
    plotmesh(PolPoints(i,:),'m.','MarkerSize', 30);hold on;
    text(PolPoints(i,1),PolPoints(i,2),PolPoints(i,3)+3,num2str(i));
end
axis equal
xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
title('Raw Polhemus Data');


tile_offset = tile - repmat([0 0 offset],size(tile,1),1);
tile_all = [tile; tile_offset];
%Assume there are a full number of tiles polhemus measurements
nTiles = floor(size(PolPoints,1)/3);

%Tile loop;
SrcPos = [];
DetPos = [];
figure(f1);
subplot(1,2,2);
for i = 1:nTiles
    
    p1 = PolPoints(1+(i-1)*3,:);
    p2 = PolPoints(2+(i-1)*3,:);
    p3 = PolPoints(3+(i-1)*3,:);
    PolPoints_tmp = ([p1; p2; p3]);
    midpoint = mean(PolPoints_tmp);
    %use measured points to map known tile optic positions via affine transformation:
    %pfrom = tile(1:3,:);
    %pto = [p1;p2;p3];
    %bsubmat=eye(3);
    %ptnum=size(pfrom,1);
    %amat=zeros(ptnum*3,9);
    %for i=1:ptnum
    %    amat(i*3-2:i*3,:)=kron(bsubmat,pfrom(i,:));
    %end
    %amat=[amat,repmat(bsubmat,ptnum,1)];
    %bvec=pto';
    %bvec=bvec(:);
    %x=amat\bvec;
    %A=reshape(x(1:9),3,3)';
    %B=x(end-2:end);
    %size_in = size(tile);
    %mappedTile = ((A*tile') + repmat(B,1,size_in(1)))';
    
    %define vector perpendicular to plane containing 3 measured positions
    plane = [];
    
    x1 = p1(1);
    y1 = p1(2);
    z1 = p1(3);
    
    x2 = p2(1);
    y2 = p2(2);
    z2 = p2(3);
    
    x3 = p3(1);
    y3 = p3(2);
    z3 = p3(3);
    
    % Calculate equation of the plane containing p1, p2, p3 and normal to
    % plane
    A = y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2);
    B = z1 * (x2 - x3) + z2 * (x3 - x1) + z3 * (x1 - x2);
    C = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    D = -(x1 * (y2 * z3 - y3 * z2) + x2 * (y3 * z1 - y1 * z3) + x3 * (y1 * z2 - y2 * z1));
    
    plane = [A B C D];
    normal = [A B C];
    normnorm = normal./sqrt(sum(normal.^2));
    
    %Use rigid transformation
    H = tile(1:3,:)'*((PolPoints_tmp - repmat(midpoint,3,1)));
    [U,~,V] = svd(H);
    Rot = V*U';
    mappedTile = (Rot*tile')' + repmat(midpoint',1,size(tile,1))';
    
    %Now projection tile in along normal
    projMappedTile = mappedTile - normnorm.*offset;
    
    %Plot measured pos, mapped tile, normal vector
    plot3(mappedTile(1:3,1),mappedTile(1:3,2),mappedTile(1:3,3),'r.','MarkerSize',30);hold on;
    plot3(mappedTile(4:7,1),mappedTile(4:7,2),mappedTile(4:7,3),'b.','MarkerSize',30);
    plot3(midpoint(1),midpoint(2),midpoint(3),'m.','MarkerSize',30);
    quiver3(midpoint(1),midpoint(2),midpoint(3),10*normnorm(1),10*normnorm(2),10*normnorm(3));
    plot3(projMappedTile(1:3,1),projMappedTile(1:3,2),projMappedTile(1:3,3),'r*','MarkerSize',20);hold on;
    plot3(projMappedTile(4:7,1),projMappedTile(4:7,2),projMappedTile(4:7,3),'b*','MarkerSize',20);
    
    SrcPos = [SrcPos ; projMappedTile(1:3,:)];
    DetPos = [DetPos ; projMappedTile(4:7,:)];
end
plotmesh(landmarks,'g.','MarkerSize', 30);hold on;
landmarkLabels = {'Nz','Iz','Ar','Al','Cz'};
for i = 1:size(landmarks)
    text(landmarks(i,1),landmarks(i,2)+3,landmarks(i,3)+3,landmarkLabels{i});
end
axis equal
xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
title('Calculated Projections');

f2 = figure;
set(f2,'Name','Polhemus-derived array');
for i = 1:size(SrcPos,1)
    plot3(SrcPos(i,1),SrcPos(i,2),SrcPos(i,3),'r.','MarkerSize',30);hold on;
    text(SrcPos(i,1),SrcPos(i,2)+3,SrcPos(i,3),['S' num2str(i)],'Color','r');
end
for i = 1:size(DetPos,1)
    plot3(DetPos(i,1),DetPos(i,2),DetPos(i,3),'b.','MarkerSize',30);hold on;
    text(DetPos(i,1),DetPos(i,2)+3,DetPos(i,3),['D' num2str(i)],'Color','b');
end
plotmesh(landmarks,'g.','MarkerSize', 30);hold on;
landmarkLabels = {'Nz','Iz','Ar','Al','Cz'};
for i = 1:size(landmarks)
    text(landmarks(i,1),landmarks(i,2)+3,landmarks(i,3)+3,landmarkLabels{i});
end
axis equal
xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
title('Polhemus-derived array');

SD.SrcPos = SrcPos;
SD.DetPos = DetPos;
SD.nSrcs = size(SrcPos,1);
SD.nDets = size(DetPos,1);
SD.Lambda = [wavelength1 wavelength2];
SD.SpatialUnit = 'mm';
SD.Landmarks = landmarks;
count = 1;
for i = 1:SD.nSrcs
    for j = 1:SD.nDets
        %Insert channel distance limit filter?
        distTmp = sqrt(sum((SD.SrcPos(i,:) - SD.DetPos(j,:)).^2));
        if distTmp <= distLimit
            dists(count) = distTmp; %might be useful later
            SD.MeasList(count,1) = i;
            SD.MeasList(count,2) = j;
            SD.MeasList(count,3) = 0;
            SD.MeasList(count,4) = 1;
            count = count+1;
        end
    end
end
SD.MeasList = [SD.MeasList; SD.MeasList];
SD.MeasList(end/2 + 1:end,4) = 2;

save(outname3D,'SD');

SD_3D = SD;