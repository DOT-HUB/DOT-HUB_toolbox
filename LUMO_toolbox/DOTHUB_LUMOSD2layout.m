function DOTHUB_LUMOSD2layout(SDinput2D,SDinput3D,groupid,outname)

%This function takes a Homer2-style SD file or variable and outputs a LUMO 
%layout file in .json format.

%############################### INPUTS ###################################

%SDinput2D  =   The path of an .SD file OR the SD variable containing 2D
%               probe layout information (for visualization purposes)
%SDinput3D  =   The path of an .SD file OR the SD variable containing 3D
%               probe layout information (usually phantom-measured for a 
%               given cap). If 3D positions are not available, just parse
%               the 2D SD instead.
%groupid    =   The intended group ID for this layout (in decimal).
%               Defaults to 1

%outname    =   The intended output name for the new .json file. Default is
%               'LUMO_LayoutFile_' datestr(clock) '.JSON'

%############################# Dependencies ###############################
%This script uses the 'loadjson.m' function from jsonlab- 1.5, which is
%available here: http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?jsonlab

% #########################################################################
% RJC, UCL, February 2020
% 
% ############################## Updates ##################################
% #########################################################################

% ############################### TO DO ###################################
% #########################################################################

% Manage variables
if ischar(SDinput2D)
    load(SDinput2D,'-mat');
    SD2D = SD; clear SD;
elseif isstruct(SDinput2D)
    SD2D = SDinput2D;
end

if ischar(SDinput3D)
    load(SDinput3D,'-mat');
    if exist('SD','var')
        SD3D = SD; clear SD;
    end
elseif isstruct(SDinput3D)
    SD3D = SDinput3D;
end

% Define groupid if not specified
if ~exist('groupid','var')
    groupid = 1;
elseif isempty(groupid)
    groupid = 1;
end

% Define Outname
if ~exist('outname','var')
    outname = ['layout_' datestr(clock,'ddmmyyyyHHMMSS') '.json'];
end

%Check both SD variables have same number of srcs and dets;
if SD2D.nSrcs ~= SD3D.nSrcs || SD2D.nDets ~= SD3D.nDets
    error('The selected 2D and 3D SD fiiles have different numbers of docks...');
end

%###########################################################

%Define group ID
jsonstruct.group_uid = groupid;

%edit boundary coordinates 
%[dimensions_2d = Coordinates2d((max_2d - min_2d + [30,30])...) # saving 2d coordinates]
min_2d = min([SD2D.DetPos; SD2D.SrcPos]);
max_2d = max([SD2D.DetPos; SD2D.SrcPos]);

%Subtract minimum, add 15 mm buffer;
buffer = 15; %mm
SD2D.DetPos = SD2D.DetPos - repmat(min_2d,SD2D.nDets,1) + buffer;
SD2D.SrcPos = SD2D.SrcPos - repmat(min_2d,SD2D.nSrcs,1) + buffer;

%Define dimensions_2d
tmp = max([SD2D.DetPos; SD2D.SrcPos]) + buffer;
jsonstruct.dimensions.dimensions_2d.x = tmp(1);
jsonstruct.dimensions.dimensions_2d.y = tmp(2);

if isfield(SD3D,'Landmarks') %ASSUME Nasion, Inion, Ar, Al, Cz
    names = {'Nasion', 'Inion', 'Ar', 'Al', 'Cz'};
    for i = 1:size(SD3D.Landmarks,1)
        jsonstruct.Landmarks{i}.name = names{i};
        jsonstruct.Landmarks{i}.x = SD3D.Landmarks(i,1);
        jsonstruct.Landmarks{i}.y = SD3D.Landmarks(i,2);
        jsonstruct.Landmarks{i}.z = SD3D.Landmarks(i,3);
    end
end

%Put 3D coords in the positive quadrant;
min_3d = min([SD3D.DetPos; SD3D.SrcPos]);
%Subtract minimum, add 15 mm buffer;
buffer = 15; %mm
SD3D.DetPos = SD3D.DetPos - repmat(min_3d,SD3D.nDets,1) + buffer;
SD3D.SrcPos = SD3D.SrcPos - repmat(min_3d,SD3D.nSrcs,1) + buffer;

%Define dimensions_3d
tmp = max([SD3D.DetPos; SD3D.SrcPos]) + buffer;
jsonstruct.dimensions.dimensions_3d.x = tmp(1);
jsonstruct.dimensions.dimensions_3d.y = tmp(2);
jsonstruct.dimensions.dimensions_3d.z = tmp(3);

%Docks structure;
nDocks = SD2D.nDets/4;
src_lab = {'a', 'b', 'c'};
for dock = 1:nDocks
    jsonstruct.docks{1,dock}.dock_id = ['dock_' num2str(dock)];
    
    for d = 1:4
        op = d;
        jsonstruct.docks{1,dock}.optodes{1,op}.optode_id = ['optode_' num2str(op)];
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.x = SD2D.DetPos((dock-1)*4+d,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.y = SD2D.DetPos((dock-1)*4+d,2);
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.x = SD3D.DetPos((dock-1)*4+d,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.y = SD3D.DetPos((dock-1)*4+d,2);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.z = SD3D.DetPos((dock-1)*4+d,3);
    end
    
    for s = 1:3
        op = 4+s;
        jsonstruct.docks{1,dock}.optodes{1,op}.optode_id = ['optode_' src_lab{s}];
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.x = SD2D.SrcPos((dock-1)*3+s,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.y = SD2D.SrcPos((dock-1)*3+s,2);
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.x = SD3D.SrcPos((dock-1)*3+s,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.y = SD3D.SrcPos((dock-1)*3+s,2);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.z = SD3D.SrcPos((dock-1)*3+s,3);
    end
end

jsonStr = jsonencode(jsonstruct);
fid = fopen(outname, 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);









