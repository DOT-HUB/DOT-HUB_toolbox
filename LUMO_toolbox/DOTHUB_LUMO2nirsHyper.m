function [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirsHyper(lumoDIRs,layoutFileNames,posCSVFileName,eventsCSVFileName)

% This script is designed to combine multiple LUMO files in to one .nirs
% file. The assumed use case is when two or more independent LUMO arrays are
% applied to a single individual. Therefore, if a .csv polhemus file is
% parsed, it should be for both arrays, measured in the same order as
% listed in lumoDIRs. If it is parsed, the function
% DOTHUB_LUMOPolhemus2SD is called to determine optode positions from the
% polhemus data. The resulting 3D optode positions overwrite the default 3D
% data from the .JSON file.
%
% This function then outputs the .nirs variables and (if saveFlag == 1 or isempty)
% produces a .nirs file with the same name as the .LUMO. Note re-running this function
% will overwrite previous .nirs derived from the parse LUMO file. Note that the .nirs
% file contains both 2D (SD) and 3D (SD_3D) info.
%
% ######################## INPUTS #########################################
%
% lumoDIRs           :   A cell array of strings, each a path to the .LUMO
%                       data directory. If not parsed, a parent directory,
%                       which should contain the .lumo files will be
%                       requested. Note that the order matlab loads the two
%                       will depend on file naming, and must match the
%                       order tiles are listed in the combined .csv file.
%
% layoutFileNames    :   A cell array of strings specifying the the path of the .json layout file
%                       (only required if not present inside .LUMO directories). If layoutFile variable is
%                       parsed it is assued to take precedence over any .json within
%                       the .LUMO directory. If layoutFile is not parsed (or parsed empty)
%                       and it is not present in the .LUMO direcroty, it will be requested.
%
% posCSVFileName     :  The path to a .csv file containing the 3D positioning data
%                       associated with all the selected .LUMO files in a
%                       single .csv within the same coordinate space with
%                       the same landmarks. This data is a four-column CSV of position label, then x, y, z (cm).
%                       The first five rows should be Nz, Iz, Ar, Al, Cz. The
%                       following rows should be SrcA, SrcB, SrcC of each of the N
%                       tiles in turn. If not parse, it is requested.
%
% eventsCSVFileName :   (Optional) The path of a .csv file containing two columns: time of event
%                       (in seconds) and event label. If the eventsCSV is parsed,
%                       it takes precedence over the contents of the .lumo file.
%
% ####################### OUTPUTS #########################################
%
% nirs                      :   the contents of the resulting .nirs file as a structure.
%
% nirsFileName              :   Full path and name of resulting .nirs file
%
% nirsFileName.nirs         :   A .nirs file saved in the same location as the
%                               first entry of the .LUMO variable. This .nirs
%                               combines the parsed .LUMO files and includes SD_3D, which
%                               is the same as SD, but with 3D source and detector
%                               locations.  SD_3D also contains a 'Landmarks' field
%                               that contains the five cranial landmark locations.
%
% ####################### Dependencies ####################################
% This script requires the TOML library for Matlab, which is available here:
% https://github.com/g-s-k/matlab-toml
%
% #########################################################################
% RJC, UCL, April 2020
% TO DO:
% Merge events from multiple files - if CondName is the same, assume events
% can be merged. At the moment I am jut eliminating replicated conditions.
% #########################################################################

% MANAGE VARIABLES  #######################################################
if ~exist('lumoDIRs','var')
    disp('Select parent directory containing .LUMO folders');
    lumoParentDIR = uigetdir(pwd,'Select parent directory containing .LUMO folders...');
    tmp = dir([lumoParentDIR '/*.LUMO']);
    for i = 1:length(tmp)
        lumoDIRs{i} = [tmp(i).folder '/' tmp(i).name];
    end
elseif isempty(lumoDIRs)
    disp('Select parent directory containing .LUMO folders');
    lumoParentDIR = uigetdir(pwd,'Select parent directory containing .LUMO folders...');
    tmp = dir([lumoParentDIR '/*.LUMO']);
    for i = 1:length(tmp)
        lumoDIRs{i} = [tmp(i).folder '/' tmp(i).name];
    end
elseif ~iscell(lumoDIRs)
    error('This function expects multiple lumoDIR inputs as a cell')
end

nFiles = length(lumoDIRs);
for ff = 1:nFiles
    %Ensure full path
    [lumoPath{ff},lumoName{ff},ext] = fileparts(lumoDIRs{ff});
    if ~strcmpi(ext,'.lumo')
        error('Each lumoDIRs input must point to a .LUMO directory');
    end
    if isempty(lumoPath{ff})
        lumoDIRs{ff} = fullfile(pwd,lumoDIRs{ff});
        [lumoPath{ff},lumoName{ff},~] = fileparts(lumoDIRs);
    end
end

if ~exist('layoutFileNames','var') %Not parsed so check exists in .LUMO
    for ff = 1:nFiles
        jsonTmp = dir([lumoDIRs{ff} '/*.json']);
        if isempty(jsonTmp) %Not contained in .LUMO, so load
            [~,tmp,~] = fileparts(lumoDIRs{ff});
            disp(['Select .json layout file associated ' tmp '.LUMO...']);
            [filename, pathname, ~] = uigetfile('*.json',['Select .json layout file associated ' tmp '.LUMO...']);
            layoutFileNames{ff} = [pathname '/' filename];
        else
            layoutFileNames{ff} = [lumoDIRs{ff} '/' jsonTmp(1).name];
        end
    end
elseif isempty(layoutFileNames) %Parsed empty so check exists in .LUMO
    for ff = 1:nFiles
        jsonTmp = dir([lumoDIRs{ff} '/*.json']);
        if isempty(jsonTmp) %Not contained in .LUMO, so load
            [~,tmp,~] = fileparts(lumoDIRs{ff});
            disp(['Select .json layout file associated ' tmp '.LUMO...']);
            [filename, pathname, ~] = uigetfile('*.json',['Select .json layout file associated ' tmp '.LUMO...']);
            layoutFileNames{ff} = [pathname '/' filename];
        else
            layoutFileNames{ff} = [lumoDIRs{ff} '/' jsonTmp(1).name];
        end
    end
end

if ~iscell(layoutFileNames) || length(layoutFileNames) ~= length(lumoDIRs)
    error('Parsed layoutFileName does not match expected size or type')
end

polhemusFlag = 1; %REQUIRED
if ~exist('posCSVFileName','var')
    disp('Load combined Polhemus data...');
    [filename, pathname, ~] = uigetfile('*.csv','Load combined Polhemus data...');
    posCSVFileName = [pathname '/' filename];
elseif isempty(posCSVFileName)
    disp('Load combined Polhemus data...');
    [filename, pathname, ~] = uigetfile('*.csv','Load combined Polhemus data...');
    posCSVFileName = [pathname '/' filename];
end

eventsFileFlag = 0;
if exist('eventsCSVFileName','var')
    if ~isempty(eventsCSVFileName)
        eventsFileFlag = 1;
    end
end

% FILE LOOP ###############################################################
for ff = 1:nFiles
    lumoDIR = lumoDIRs{ff};
    layoutFileName = layoutFileNames{ff};
    
    % LOAD DATA ###############################################################
    %Load toml contents
    disp(['Loading .toml contents of LUMO file ' num2str(ff) '...']);
    metadata = toml.read([lumoDIR '/metadata.toml']);
    recordingdata = toml.read([lumoDIR '/' metadata.file_names.recordingdata_file]); %To be updated to recording.toml
    events = toml.read([lumoDIR '/' metadata.file_names.event_file]);
    if isfield(metadata.file_names,'hardware_file') %Updated file naming structure
        hardware = toml.read([lumoDIR '/' metadata.file_names.hardware_file]);
    elseif ~contains(metadata.file_names.layout_file,'json') %Old file naming structure (layout_file is hardware.toml!)
        hardware = toml.read([lumoDIR '/' metadata.file_names.layout_file]);
    else
        error('No hardware.toml or layout.toml found in LUMO directory');
    end
    
    %Load layout JSON
    layoutData = jsondecode(fileread(layoutFileName));
    
    %Load log file
    logFileID = fopen([lumoDIR '/' metadata.file_names.log_file]);
    logtxt = textscan(logFileID,'%s');
    
    %Load intensity binaries
    nIntFiles = length(metadata.intensity_files);
    for i = 1:nIntFiles
        fname = metadata.intensity_files{1,i}.file_name;
        time_range = metadata.intensity_files{1,i}.time_range;
        fileID = fopen([lumoDIR '/' fname]);
        int(i).magic = fread(fileID,1);
        int(i).major_version = fread(fileID,1);
        int(i).minor_version = fread(fileID,1);
        int(i).patch = fread(fileID,1,'uint8',4);
        int(i).recordings_per_frame = fread(fileID,1,'uint64');
        int(i).nframe = fread(fileID,1,'uint64',24);
        int(i).data = fread(fileID,[int(i).recordings_per_frame int(i).nframe],'float32')';
    end
    
    % #########################################################################
    % COMBINE MULTIPLE INTENSITY FILES and yield .nirs-style 'd' variable
    intDataAll = [];
    for i = 1:length(int)
        intDataAll = [intDataAll; int(i).data];
    end
    % #########################################################################
    
    % Extract useful things for conversion ####################################
    nChan = recordingdata.variables.n_chans;
    chansList = reshape(recordingdata.variables.chans_list,3,nChan)';
    fs = recordingdata.variables.framerate;
    nFrames = recordingdata.variables.number_of_frames;
    nDocks = size(layoutData.docks,1);
    nodes = recordingdata.variables.nodes;
    nNodes = length(nodes);
    
    % Define t ################################################################
    t = (0:1/fs:nFrames/fs - 1/fs)';
    
    % Define SD ###############################################################
    SD.nSrcs = recordingdata.variables.n_srcs;
    SD.nDets = recordingdata.variables.n_dets;
    SD.Lambda = recordingdata.variables.wavelength;
    SD.SpatialUnit = 'mm';
    MLAtmp = recordingdata.variables.chans_list_act'; %temporary measlistact
    
    % Now determine optode positions from 2D information in layout JSON file
    % Also use same loop to extract source power information
    SD.SrcPos = zeros(SD.nSrcs,3);
    SD.DetPos = zeros(SD.nDets,3);
    SD.SrcPowers = zeros(SD.nSrcs*2,2);
    SD.SrcPowers(1:SD.nSrcs,2) = 1;
    SD.SrcPowers(SD.nSrcs+1:end,2) = 2;
    for n = 1:nNodes
        nid = nodes(n);
        for det = 1:4
            SD.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_2d.x;
            SD.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_2d.y;
        end
        for src = 1:3
            SD.SrcPos(src+(n-1)*3,1) = layoutData.docks(nid).optodes(src+4).coordinates_2d.x;
            SD.SrcPos(src+(n-1)*3,2) = layoutData.docks(nid).optodes(src+4).coordinates_2d.y;
            
            SD.SrcPowers(src+(n-1)*3,1) = hardware.Hub.Group.Node{1,n}.Source{1,src*2-1}.Source_power;
            SD.SrcPowers(src+(n-1)*3 + SD.nSrcs,1) = hardware.Hub.Group.Node{1,n}.Source{1,src*2}.Source_power;
        end
    end
    
    %Make SD.MeasList
    SD.MeasList = ones(size(chansList,1),4);
    SD.MeasList(:,1:2) = chansList(:,1:2);
    SD.MeasList(:,4) = chansList(:,3);
    [SD.MeasList, tmpInd] = sortrows(SD.MeasList,[4,1,2]);
    
    % ########## Use sorted ML order to correctly sort data and measlistact
    d = intDataAll(:,tmpInd);
    SD.MeasListAct = MLAtmp(tmpInd);
    % ######################################################################
    
    % Include 3D information for future 
        
    SD3D = SD;
    for n = 1:nNodes
        nid = nodes(n);
        for det = 1:4
            SD3D.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_3d.x;
            SD3D.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_3d.y;
            SD3D.DetPos(det+(n-1)*4,3) = layoutData.docks(nid).optodes(det).coordinates_3d.z;
        end
        for src = 1:3
            SD3D.SrcPos(src+(n-1)*3,1) = layoutData.docks(nid).optodes(src+4).coordinates_3d.x;
            SD3D.SrcPos(src+(n-1)*3,2) = layoutData.docks(nid).optodes(src+4).coordinates_3d.y;
            SD3D.SrcPos(src+(n-1)*3,3) = layoutData.docks(nid).optodes(src+4).coordinates_3d.z;
        end
    end
    
    % Define stim vector, s ###################################################
    if eventsFileFlag %Read events from simple 2-column CSV file.
        delimiter = ',';
        startRow = 1;
        formatSpec = '%f%s%[^\n\r]';
        fileID = fopen(eventsCSVFileName,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        if isempty(dataArray)
            warning('Specified eventsCSV file contains no events!');
            s = zeros(size(t,1),1);
            CondNames = {''};
        else
            
            timeStamp = dataArray{1};
            eventStr = dataArray{2};
            [tmp,~,occuranceInd] = unique(eventStr,'stable'); %find unique events
            CondNames = arrayfun(@cellstr,tmp);
            nCond = size(CondNames,2); %number of conditions
            s = zeros(size(t,1),nCond); %preassign s matrix
            for i = 1:nCond
                timeStampTmp = timeStamp(occuranceInd==i); %Can have multiple entries
                [~,indTmp] = min(abs(repmat(t,1,length(timeStampTmp))-repmat(timeStampTmp,size(t,1),1)));
                s(indTmp,i) = 1;
            end
        end
        
    else
        
        %Now check if there are any serial events
        if isempty(fieldnames(events))
            s = zeros(size(t,1),1); %return empty s vector
            CondNames = {''};
        else
            for i = 1:size(events.events,2)
                timeStamp(i) = str2num(events.events{1,i}.Timestamp)*1e-3; %to seconds
                eventStr{i} = events.events{1,i}.name;
            end
            [tmp,~,occuranceInd] = unique(eventStr,'stable'); %find unique events, maintain order in which they occured
            CondNames = arrayfun(@cellstr,tmp);
            nCond = size(CondNames,2); %number of conditions
            s = zeros(size(t,1),nCond); %preassign s matrix
            for i = 1:nCond
                timeStampTmp = timeStamp(occuranceInd==i); %Can have multiple entries
                [~,indTmp] = min(abs(repmat(t,1,length(timeStampTmp))-repmat(timeStampTmp,size(t,1),1)));
                s(indTmp,i) = 1;
            end
        end
    end
    
    % Define AUX matrix #######################################################
    aux = zeros(length(t),8);
    
    % OUTPUT ##################################################################
    allnirs{ff}.SD = SD;
    allnirs{ff}.SD3D = SD3D;
    allnirs{ff}.aux = aux;
    allnirs{ff}.d = d;
    allnirs{ff}.s = s;
    allnirs{ff}.t = t;
    allnirs{ff}.CondNames = CondNames;
    
    tLength(ff) = length(t); %Assume start point is the same, remove extra datapoints from end
end

%%
mnLength = min(tLength);
dcomb = [];
scomb = [];
SD3Dcomb.SrcPos = [];
SD3Dcomb.DetPos = [];
SD3Dcomb.MeasList = [];
SD3Dcomb.MeasListAct = [];
SD3Dcomb.nSrcs = 0;
SD3Dcomb.nDets = 0;
SD3Dcomb.SrcPowers = [];
CondNamesAll = {};
for ff = 1:nFiles
    dtmp = allnirs{ff}.d;
    dtmp = dtmp(1:mnLength,:);
    dcomb = [dcomb dtmp];

    stmp = allnirs{ff}.s;
    stmp = stmp(1:mnLength,:);
    scomb = [scomb stmp];
    
    SD3Dcomb.SrcPos = [SD3Dcomb.SrcPos; allnirs{ff}.SD3D.SrcPos];
    SD3Dcomb.DetPos = [SD3Dcomb.DetPos; allnirs{ff}.SD3D.DetPos];   
    
    MeasListTmp = allnirs{ff}.SD3D.MeasList;
    MeasListTmp(:,1) = MeasListTmp(:,1)+SD3Dcomb.nSrcs;
    MeasListTmp(:,2) = MeasListTmp(:,2)+SD3Dcomb.nDets;
    SD3Dcomb.MeasList = [SD3Dcomb.MeasList; MeasListTmp];
    SD3Dcomb.MeasListAct = [SD3Dcomb.MeasListAct; allnirs{ff}.SD3D.MeasListAct];
    SD3Dcomb.SrcPowers = [SD3Dcomb.SrcPowers; allnirs{ff}.SD3D.SrcPowers];
    SD3Dcomb.nSrcs = max(SD3Dcomb.MeasList(:,1));
    SD3Dcomb.nDets = max(SD3Dcomb.MeasList(:,2));
    
    CondNamesAll = [CondNamesAll allnirs{ff}.CondNames];
end

%Now sort MeasList and data!
[~,sind] = sort(SD3Dcomb.MeasList(:,4));
for i = 1:4
    SD3Dcomb.MeasList(:,i) = SD3Dcomb.MeasList(sind,i);
end
SD3Dcomb.MeasListAct = SD3Dcomb.MeasListAct(sind);
dcomb = dcomb(:,sind);

if polhemusFlag %If polhemus information is parsed, calculate S-D positions from that file
    fprintf(['Using ' posCSVFileName ' to define optode positions...\n']);
    [SD_POL, SD3DFileName] = DOTHUB_LUMOpolhemus2SD3D(posCSVFileName,[],0); %Set saveflag to 0 because this function assumes all source-detector pairs form channels - not true in this case.
    SD3Dcomb.SrcPos = SD_POL.SrcPos;
    SD3Dcomb.DetPos = SD_POL.DetPos;
    SD3Dcomb.Landmarks = SD_POL.Landmarks;
end

%Straight assignements
tcomb = allnirs{1}.t(1:mnLength);
Auxcomb = allnirs{1}.aux; %Will be zeros anyway
SD3Dcomb.SpatialUnit = allnirs{1}.SD.SpatialUnit;
SD3Dcomb.Lambda = allnirs{1}.SD.Lambda;

%Handle stimuli;
%For now remove conditions that are repeated (so will wipe events assigned
%same condition name for ff>1
[CondNamesComb,IA] = unique(CondNamesAll,'stable');
scomb = scomb(:,IA);

%Now create 2DSD - shift each additional array to below the last
SDcomb = SD3Dcomb;
if isfield(allnirs{1}.SD3D,'Landmarks') %Don't want landmarks in 2D
    SD3Dcomb.Landmarks = allnirs{1}.SD3D.Landmarks;
end
SDcomb.SrcPos = [];
SDcomb.DetPos = [];

for ff = 1:nFiles
    %First centre the array
    tmpSP = allnirs{ff}.SD.SrcPos;
    tmpDP = allnirs{ff}.SD.DetPos;
    tmpALL = [tmpSP; tmpDP];
    tmpSP(:,1) = tmpSP(:,1) - min(tmpALL(:,1)) - range(tmpALL(:,1))/2;
    tmpDP(:,1) = tmpDP(:,1) - min(tmpALL(:,1)) - range(tmpALL(:,1))/2;
    tmpSP(:,2) = tmpSP(:,2) - min(tmpALL(:,2)) - range(tmpALL(:,2))/2;
    tmpDP(:,2) = tmpDP(:,2) - min(tmpALL(:,2)) - range(tmpALL(:,2))/2;
    
    if ff>1 %offset 
        tmpALL = [tmpSP; tmpDP];
        yOffset = min([SDcomb.SrcPos(:,2); SDcomb.DetPos(:,2)]) - range(tmpALL(:,2));
        tmpSP(:,2) = tmpSP(:,2) + yOffset;
        tmpDP(:,2) = tmpDP(:,2) + yOffset;
    end
        
    SDcomb.SrcPos = [SDcomb.SrcPos; tmpSP];
    SDcomb.DetPos = [SDcomb.DetPos; tmpDP];   
end

nirs.SD = SDcomb;
nirs.SD3D = SD3Dcomb;
nirs.aux = Auxcomb;
nirs.d = dcomb;
nirs.s = scomb;
nirs.t = tcomb;
nirs.CondNames = CondNamesComb;


disp(['Saving SD3D file to ' SD3DFileName '...']);
SD3D = SD3Dcomb;
save(SD3DFileName,'SD3D');

outName = lumoName{1};
for ff = 2:nFiles
    outName = [outName '+' lumoName{ff}];
end
nirsFileName = fullfile(lumoPath{1}, [outName '.nirs']);
disp(['Saving file to ' nirsFileName ' ...\n']);
save(nirsFileName,'-struct','nirs','-v7.3');




