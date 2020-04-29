function DOTHUB_LUMOSD2layout(SDinput2D,SDinput3D,groupid,outname)

%This function takes a Homer2-style SD file or variable and outputs a LUMO 
%layout file in .json format.

%############################### INPUTS ###################################

%SDinput2D  =   The path of an .SD file OR the SD variable containing 2D
%               probe layout information (for visualization purposes)
%SDinput3D  =   The path of an .SD file OR the SD variable containing 3D
%               probe layout information (usually phantom-measured for a 
%               given cap). If 3D positions are not available, just parse
%               the 2D SD in stead.
%groupid    =   The intended group ID for this layout (in decimal)

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
    SD_2D = SD;
elseif isstruct(SDinput2D)
    SD_2D = SDinput2D;
end

if ischar(SDinput3D)
    load(SDinput3D,'-mat');
    SD_3D = SD;
elseif isstruct(SDinput3D)
    SD_3D = SDinput3D;
end

% Define groupid if not specified
if ~exist('groupid','var')
    groupid = 1;
elseif isempty(groupid)
    groupid = 1;
end

% Define Outname
if ~exist('outname','var')
    outname = ['LUMO_LayoutFile_' datestr(clock) '.JSON'];
end

%Check both SD variables have same number of srcs and dets;
if SD_2D.nSrcs ~= SD_3D.nSrcs || SD_2D.nDets ~= SD_3D.nDets
    error('The selected 2D and 3D SD fiiles have different numbers of docks...');
end

%###########################################################

%Define group ID
jsonstruct.group_uid = groupid;

%edit boundary coordinates 
%[dimensions_2d = Coordinates2d((max_2d - min_2d + [30,30])...) # saving 2d coordinates]
min_2d = min([SD_2D.DetPos; SD_2D.SrcPos]);
max_2d = max([SD_2D.DetPos; SD_2D.SrcPos]);

%Subtract minimum, add 15 mm buffer;
buffer = 15; %mm
SD_2D.DetPos = SD_2D.DetPos - repmat(min_2d,SD.nDets,1) + buffer;
SD_2D.SrcPos = SD_2D.SrcPos - repmat(min_2d,SD.nSrcs,1) + buffer;
tmp = max([SD_2D.DetPos; SD_2D.SrcPos]) + buffer;

%Define dimensions_2d
jsonstruct.dimensions.dimensions_2d.x = tmp(1);
jsonstruct.dimensions.dimensions_2d.y = tmp(2);

%Put 3D coords in the positive quadrant;
min_3d = min([SD_3D.DetPos; SD_3D.SrcPos]);
max_3d = max([SD_3D.DetPos; SD_3D.SrcPos]);

%Subtract minimum, add 15 mm buffer;
buffer = 15; %mm
SD_3D.DetPos = SD_3D.DetPos - repmat(min_3d,SD.nDets,1) + buffer;
SD_3D.SrcPos = SD_3D.SrcPos - repmat(min_3d,SD.nSrcs,1) + buffer;
tmp = max([SD_3D.DetPos; SD_3D.SrcPos]) + buffer;

%Define dimensions_3d
jsonstruct.dimensions.dimensions_3d.x = tmp(1);
jsonstruct.dimensions.dimensions_3d.y = tmp(2);
jsonstruct.dimensions.dimensions_3d.z = tmp(3);

%Docks structure;
nDocks = SD_2D.nDets/4;
src_lab = {'a', 'b', 'c'};
for dock = 1:nDocks
    jsonstruct.docks{1,dock}.dock_id = ['dock_' num2str(dock)];
    
    for d = 1:4
        op = d;
        jsonstruct.docks{1,dock}.optodes{1,op}.optode_id = ['optode_' num2str(op)];
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.x = SD_2D.DetPos((dock-1)*4+d,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.y = SD_2D.DetPos((dock-1)*4+d,2);
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.x = SD_3D.DetPos((dock-1)*4+d,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.y = SD_3D.DetPos((dock-1)*4+d,2);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.z = SD_3D.DetPos((dock-1)*4+d,3);
    end
    
    for s = 1:3
        op = 4+s;
        jsonstruct.docks{1,dock}.optodes{1,op}.optode_id = ['optode_' src_lab{s}];
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.x = SD_2D.SrcPos((dock-1)*3+s,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_2d.y = SD_2D.SrcPos((dock-1)*3+s,2);
        
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.x = SD_3D.SrcPos((dock-1)*3+s,1);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.y = SD_3D.SrcPos((dock-1)*3+s,2);
        jsonstruct.docks{1,dock}.optodes{1,op}.coordinates_3d.z = SD_3D.SrcPos((dock-1)*3+s,3);
    end
end

jsonStr = jsonencode(jsonstruct);
fid = fopen(outname, 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);









