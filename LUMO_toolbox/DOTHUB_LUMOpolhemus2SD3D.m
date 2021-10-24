 function [SD3D, SD3DFileName] = DOTHUB_LUMOpolhemus2SD3D(posCSVFileName,infantFlag,saveFlag)

%This script takes the polhemus data associated with the measurement
%points on the backs of the tiles and creates a 3D map of source and
%detector positions associated with that array via a plane and projection
%method. It then creates an SD_3D structure assuming all source-detector
%pairs form a channel
%
% ASSUMPTIONS:
% [1] The user has collected the 3 points for an integer number of tiles
% [2] These points are collected in the correct order (anticlockwise from A)
%
% INPUTS: ################################################################
% posCSVFileName =      the .csv file name This data is a four-column CSV of
%                       position label, then x, y, z coordinate (assumed in cm).
%                       The first five rows should be Nz, Iz, Ar, Al, Cz. The
%                       following rows should be SrcA, SrcB, SrcC of each of the N
%                       tiles in turn.
%
% infantFlag =          true if light-guides associated with this array are the infant style
%                       If a vector of length nTiles, take tile-specific value.
%
% saveFlag =            true if resulting SD file should be saved (default true)
%
% OUTPUTS: ###############################################################
% SD_3D = The SD file containing the 3D information associated with polhemus
% data. Also contains the a 'Landmarks' variable that saves the cranial
% landmarks (Nz Iz Ar Al Cz). Also save out into an .SD file with a name
% matching the input .csv appended with _3D.
%
% SD_3DFileName = Filename of the saved SD_3D variable.
%
% RJC, UCL, December 2019
%
% UPDATES ################################################################
% RJC - Dec 2019 - Converted the affine transformation of the tile positions 
%                  to the polhemus measurements to a rigid transformation.
%
% TO DO LIST #############################################################
% Allow multiple file inputs and averaging of polhemus data?

%Manage Inputs ###########################################################
%#########################################################################
if ~exist('posCSVFileName','var')
    [file,path] = uigetfile('*.csv','Select Polhemus data set (.csv)','MultiSelect','on');
    posCSVFileName = fullfile(path,file);
end

mixedFlag = 0;
if ~exist('infantFlag','var')
    answer = questdlg('Which light-guides were used in this array?','Select light guide type','ADULT','INFANT','MIXED','ADULT');
    if strcmp(answer,'INFANT')
        infantFlag = 1;
    elseif strcmp(answer,'ADULT')
        infantFlag = 0;
    elseif  strcmp(answer,'MIXED')
        mixedFlag = 1;
    else
        return
    end
elseif isempty(infantFlag)
    answer = questdlg('Which light-guides were used in this array?','Select light guide type','ADULT','INFANT','MIXED','ADULT');
    if strcmp(answer,'INFANT')
        infantFlag = 1;
    elseif strcmp(answer,'ADULT')
        infantFlag = 0;
    elseif  strcmp(answer,'MIXED')
        mixedFlag = 1;
    else
        return
    end   
end

if ~exist('saveFlag','var')
    saveFlag = 1;
end

%Target SD output names and force full path name
[path,name,ext] = fileparts(posCSVFileName);
if isempty(path)
    path = pwd;
end
if isempty(ext)
    ext = '.csv';
end
posCSVFileName = fullfile(path,[name ext]);
SD3DFileName = fullfile(path,[name '.SD3D']);

%Wavelengths
wavelength1 = 735;
wavelength2 = 850;

%#########################################################################
%#########################################################################
%Load data
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%f%f%[^\n\r]';
fileID = fopen(posCSVFileName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
allPos = [dataArray{2} dataArray{3} dataArray{4}]*10; %Polhemus output in cm, convert to mm.%#################################################%#################################################%#################################################%#################################################

%Crop out landmarks and polhemus measurement points
landmarks = allPos(1:5,:);
PolPoints = allPos(6:end,:);

if mod(size(PolPoints,1),3)~=0
    error('The selected .csv file does not seem to contain a multiple of 3 points as is expected for LUMO');
else
    nTiles = size(PolPoints,1)/3;
end
%Define offsets vector
if mixedFlag
    prompt = 'Enter nTile space-separated values of 1 (Adult) and 0 (Infant)';
    dlgtitle = 'Light guide specification';
    dims = [1 45];
    answer = inputdlg(prompt,dlgtitle,dims);
    answer = str2num(cell2mat(answer));
    if isempty(answer)
        error('No light guide specification provided')
    else
        offset = zeros(nTiles,1);
        offset(answer==1) = 18.75;
        offset(answer==0) = 12.94;
        if any(offset==0)
            error('Error in light guide specification')
        end
    end      
elseif infantFlag==1
    offset = ones(nTiles,1)*12.94; %Infant offset length (mm)
elseif infantFlag==0
    offset = ones(nTiles,1)*18.75; %Adult offset length (mm)
end

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

%Tile optic positions definition, ordered source then detector
%Correct as per May 2019
% Source-detector organization per tile
%  x: Source        s1       d4
%  o: detectors      x       o
%
%              d2 o      o d3    x s3
% (d3 central
%  detector)         x       o
%                   s2       d1
tile(1,:) = [0 10.81 0];
tile(2,:) = [-10.81*cosd(30) -10.81*sind(30) 0];
tile(3,:) = [10.81*cosd(30) -10.81*sind(30) 0];
tile(4,:) = [0 -8.88 0];
tile(5,:) = [-8.88*cosd(30) 8.88*sind(30) 0];
tile(6,:) = [0 0 0];
tile(7,:) = [8.88*cosd(30) 8.88*sind(30) 0];

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
    projMappedTile = mappedTile - normnorm.*offset(i);
    
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
pause(0.5);drawnow

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
pause(0.5);drawnow

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
        %if distTmp <= distLimit
            dists(count) = distTmp; %might be useful later
            SD.MeasList(count,1) = i;
            SD.MeasList(count,2) = j;
            SD.MeasList(count,3) = 0;
            SD.MeasList(count,4) = 1;
            count = count+1;
        %end
    end
end
SD.MeasList = [SD.MeasList; SD.MeasList];
SD.MeasList(end/2 + 1:end,4) = 2;
SD3D = SD;

if saveFlag
    fprintf(['Saving SD3D to ' SD3DFileName ' ...\n']);
    save(SD3DFileName,'SD3D');
end


