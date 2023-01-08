%HELLO GINA
% Hi Rob
% LUMO2NIRS 5 wavelengths
function [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs_5wav(lumoDIR,layoutFileName,posCSVFileName,eventsCSVFileName)

%This script reads a .LUMO file (which is actually just a directory)
%and (if it is not contained within the .LUMO - which it should be)
%also loads the associated .JSON layout from a seperate file.

%If subject-specific registration data (from polhemus or other method) is available
%in .cvs format, this too can be parsd. If it is parsed, the function
%DOTHUB_LUMOPolhemus2SD is called to determine optode positions from the
%polhemus data. The resulting 3D optode positions overwrite the default 3D
%data from the .JSON file.

%This function then outputs the .nirs variables and (if saveFlag == 1 or isempty)
%produces a .nirs file with the same name as the .LUMO. Note re-running this function
%will overwrite previous .nirs derived from the parse LUMO file. Note that the .nirs
%file contains both 2D (SD) and 3D (SD_3D) info.

%######################## INPUTS ##########################################

% lumoDIR :     The path of the .LUMO data directory (if not parsed, requested)

% layoutFileName    :   (Optional) The path of the .json layout file (only required if not
%                       present inside .LUMO directory). If layoutFile variable is
%                       parsed it is assued to take precedence over any .json within
%                       the .LUMO directory. If layoutFile is not parsed (or parsed empty)
%                       and it is not present in the .LUMO direcroty, it will be requested.

% posCSVFileName    :   (Optional) The path of a .csv file containing the 3D positioning data associated
%                       with this LUMO recording. This data is a four-column CSV of
%                       position label, then x, y, z coordinate (assumed in cm).
%                       The first five rows should be Nz, Iz, Ar, Al, Cz. The
%                       following rows should be SrcA, SrcB, SrcC of each of the N
%                       tiles in turn. If polhemusCSV is not parsed, the 3D
%                       information in the .json layout file is assumed.

% eventsCSVFileName :   (Optional) The path of a .csv file containing two columns: time of event
%                       (in seconds) and event label. If the eventsCSV is parsed,
%                       it takes precedence over the contents of the .lumo file.

%######################## OUTPUTS #########################################

% nirs                      :   the contents of the resulting .nirs file as a structure.

% nirsFileName              :   Full path and name of resulting .nirs file

% nirsFileName.nirs         :   A .nirs file saved in the same location as the 
%                               .LUMO variable. Has extra variable: SD_3D, which
%                               is the same as SD, but with 3D source and detector 
%                               locations.  SD_3D also contains a 'Landmarks' field
%                               that contains the five cranial landmark locations.

%######################## Dependencies ####################################
%This script requires the TOML library for Matlab, which is available here:
%https://github.com/g-s-k/matlab-toml

% #########################################################################
% RJC, UCL, Summer 2019
%
% ############################# Updates ###################################
% RJC 20190815 - Updated to correctly re-arrange chan_list and data to match
% homer2 arrangement (sort by wavelength, source, detector)
% RJC 20190829 - Updated to work when events field is empty
% RJC 20191128 - Updated to handle multiple binary files (long datasets)
% RJC 20191202 - Updated to parse source powers into SD.SrcPowers (organized a la Homer2)
% Also note that the plan is to cap the max channel length in the .json layout files to
% 60 mm by default, which will trickle down to this conversion.
% RJC 20200217 - Tidied for GITHUB first commit.
% EGJ 20200923 - Updated to handle zero intensity values
% KET 20210808 - Multiwavelength functionality.
% #########################################################################

% #########################################################################
% TO DO:
% Update event marking options to allow both .csv and automated/serial
% #########################################################################

% MANAGE VARIABLES  #######################################################
if ~exist('lumoDIR','var')
    disp('Select .LUMO directory...');
    lumoDIR = uigetdir(pwd,'Select .LUMO directory');
elseif isempty(lumoDIR)
    disp('Select .LUMO directory...');
    lumoDIR = uigetdir(pwd,'Select .LUMO directory');
elseif ~exist(lumoDIR,'dir')
    disp('Specified directory not found, please select .LUMO directory...');
    lumoDIR = uigetdir(pwd,'Select .LUMO directory');
end

%Ensure full path
[lumoPath,lumoName,ext] = fileparts(lumoDIR);
if ~strcmpi(ext,'.lumo')
    error('The lumoDIR input must point to a .LUMO directory');
end
if isempty(lumoPath)
    lumoDIR = fullfile(pwd,lumoDIR);
    [lumoPath,lumoName,~] = fileparts(lumoDIR);
end

if ~exist('layoutFileName','var') %Not parsed so check exists in .LUMO
    jsonTmp = dir([lumoDIR '/layout.json']);
    if isempty(jsonTmp) %Not contained in .LUMO, so load
        disp('Layout file not parsed or found, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        disp('Using default layout.json within .LUMO directory');
        layoutFileName = [lumoDIR '/' jsonTmp(1).name];
    end
elseif isempty(layoutFileName)
 jsonTmp = dir([lumoDIR '/layout.json']);
    if isempty(jsonTmp) %Not contained in .LUMO, so load
        disp('Layout file not parsed or found, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        disp('Using layout.json within .LUMO directory');
        layoutFileName = [lumoDIR '/' jsonTmp(1).name];
    end
end

polhemusFlag = 0;
if exist('posCSVFileName','var')
    if ~isempty(posCSVFileName)
        polhemusFlag = 1;
    end
end

eventsFileFlag = 0;
if exist('eventsCSVFileName','var')
    if ~isempty(eventsCSVFileName)
        eventsFileFlag = 1;
    end
end

% LOAD DATA  ###############################################################
%Load toml contents
disp('Loading .toml contents...');
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
disp('Loading intensity data');
nIntFiles = length(metadata.intensity_files);
for i = 1:nIntFiles
    disp(['Loading intensity data file ' num2str(i) ' of ' num2str(nIntFiles) '...']);
    fname = metadata.intensity_files{1,i}.file_name;
    %time_range = metadata.intensity_files{1,i}.time_range;
    fileID = fopen([lumoDIR '/' fname]);
    int(i).magic = fread(fileID,1);
    int(i).major_version = fread(fileID,1);
    int(i).minor_version = fread(fileID,1);
    int(i).patch = fread(fileID,1,'uint8',4);
    int(i).recordings_per_frame = fread(fileID,1,'uint64');
    int(i).nframe = fread(fileID,1,'uint64',24); %skip checksum
    int(i).data = fread(fileID,[int(i).recordings_per_frame int(i).nframe],'float32')';
    fclose(fileID);
end

% #########################################################################
% COMBINE MULTIPLE INTENSITY FILES and yield .nirs-style 'd' variable
intDataAll = [];
for i = 1:length(int)
    intDataAll = [intDataAll; int(i).data];
end

% %%%%..GCVL..#### ************************************************
	% show raw data
	printRawData(intDataAll);
% %%%%..GCVL..#### ************************************************



% #########################################################################
% Extract useful things for conversion ####################################
nChan = recordingdata.variables.n_chans;
chansListOrig = reshape(recordingdata.variables.chans_list,3,nChan)';
fs = recordingdata.variables.framerate;
nFrames = recordingdata.variables.number_of_frames;
nDocks = size(layoutData.docks,1);
nodes = recordingdata.variables.nodes;
nNodes = length(nodes);
nSources = length(recordingdata.variables.nodes);
nDets = recordingdata.variables.n_dets;

% Map from LUMO 2-wav chansListOrig to 5-wav version
% DO NOT EDIT ##########################################
% LUMO                          Actual (5 wav)
% Source 1, Wav 1 (735nm)       Source 1, Wav 1 (720nm)
% Source 1, Wav 2 (850nm)       Source 1, Wav 2 (760nm)
% Source 2, Wav 1 (735nm)       Source 1, Wav 3 (800nm)
% Source 2, Wav 2 (850nm)       Source 1, Wav 4 (850nm)
% Source 3, Wav 1 (735nm)       Source 1, Wav 5 (890nm)
% Source 3, Wav 2 (850nm)       Source 1, Wav 6 (NULL)
% etc....
% ######################################################

% Create mapping matrix from LUMO source to 5wav source
%(1,1) = (1,1)
%(1,2) = (1,2)
%(2,1) = (1,3)
%(2,2) = (1,4)
%(3,1) = (1,5)
%(3,2) = (1,6)

for node = 1:nNodes
    for src = 1:3
        for wav = 1:2
            mappingMat{src+3*(node-1),wav} = [node, 2*src + (wav-2)];
        end
    end
end

%Use mapping matrix to create chansList from chansListOrig
chansList = chansListOrig; %preallocate to Orig version
for i = 1:size(chansListOrig,1)
    tmp = chansListOrig(i,[1 3]);
    chansList(i,1) = mappingMat{tmp(1),tmp(2)}(1);
    chansList(i,3) = mappingMat{tmp(1),tmp(2)}(2);
end

% Remove dark measurements and save (6th wavelength)
% %%%%..GCVL..#### ************************************************
darkInd = chansList(:,3)==6;    % chansList has 96 rows and 3 cols    Col-1 rows 1-48=1 49-96=2   Col-2 sequential repeat 1,2,3,4,5,6,7,8   Col-3 8x1, 8x2, 8x3, 8x4, 8x5, 8x6, then repeats
d_dark = intDataAll(:,darkInd);
intDataAll = intDataAll(:,~darkInd);  %select only LED channels  (80 in total)   % %%%%..GCVL..#### ************************************************

	%show dark data
	printdarkchans(chansList, darkInd, d_dark);
% %%%%..GCVL..#### ************************************************
	
% Update chansList
chansList = chansList(~darkInd,:);

% Define t ################################################################
t = (0:1/fs:nFrames/fs - 1/fs)';

% Define SD ###############################################################
% SD.nSrcs = recordingdata.variables.n_srcs;
SD.nSrcs = nSources;
SD.nDets = recordingdata.variables.n_dets;
% SD.Lambda = recordingdata.variables.wavelength;
% HARDCODED WAVELENGTHS
SD.Lambda = [720,760,800,850,890];
SD.SpatialUnit = 'mm';
MLAtmp = recordingdata.variables.chans_list_act'; %temporary saturation list

% Now determine optode positions from 2D information in layout JSON file
% Also use same loop to extract source power information
SD.SrcPos = zeros(SD.nSrcs,3); %12 sources 3D positions
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
    src = 1; %Single source in AV3
    SD.SrcPos(src+(n-1),1) = layoutData.docks(nid).optodes(src+4).coordinates_2d.x;
    SD.SrcPos(src+(n-1),2) = layoutData.docks(nid).optodes(src+4).coordinates_2d.y;

    SD.SrcPowers(src+(n-1),1) = hardware.Hub.Group.Node{1,n}.Source{1,src*2-1}.Source_power;
    SD.SrcPowers(src+(n-1) + SD.nSrcs,1) = hardware.Hub.Group.Node{1,n}.Source{1,src*2}.Source_power;
end

%Make SD.MeasList for the light measurements
SD.MeasList = ones(size(chansList,1),4);
SD.MeasList(:,1:2) = chansList(:,1:2);
SD.MeasList(:,4) = chansList(:,3);
[SD.MeasList, tmpInd] = sortrows(SD.MeasList,[4,1,2]);
SD.MeasListAct = ones(size(SD.MeasList,1),1);

% ########## Use sorted ML order to correctly sort data and measlistact ##
d = intDataAll(:,tmpInd);   %re-order channels per tmpInd   % %%%%..GCVL..#### ************************************************

% %%%%..GCVL..#### ************************************************
	%show nirs array "d"
	printnirsd(d, intDataAll);
% %%%%..GCVL..#### ************************************************

%% flag zero or negative values of d
fprintf('\n Flag zero or negative values of d ');
for rowtest = 1:size(d,1)
	for coltest = 1:size(d,2)
		if d(rowtest, coltest) <= 0
			fprintf('\n  row  col  d    %d %d %f ', rowtest, coltest, d(rowtest,coltest));
			d(rowtest, coltest) = 0.000001;
		end
	end
end	
	
% %%%%..GCVL..#### ************************************************



MLAtmp = MLAtmp(tmpInd);

% ########## Consider channels flagged as saturated ######################
% The default LUMO chan_list_act marks channels that saturate at any point
% in the recording. This can be unfairly conservative as transient motion
% can cause a transient saturation flag. For now, the saturation list is
% saved into SD.MeasListActSat, but the default is set to all active.
SD.MeasListActSat = MLAtmp;
% ########################################################################

% Include 3D information for future
if polhemusFlag %If polhemus information is parsed, calculate S-D positions from that file
    fprintf(['Using ' posCSVFileName ' to define SD3D...\n']);
    [SD_POL, SD3DFileName] = DOTHUB_LUMOpolhemus2SD3D_5(posCSVFileName); %This line saves the .SD3D
    SD_POL.MeasListAct = SD.MeasListAct;
    
    if nNodes<nDocks %Crop SD file accordingly if the number of docks populated in the datafile differs from the number in the SD_3D data
        SD3D = SD; %Define based on SD which contains correct measlist
        for n = 1:nNodes
            nid = nodes(n);
            for det = 1:4
                SD3D.DetPos(det+(n-1)*4,:) = SD.DetPos(det+(nid-1)*4,:);
            end
            src = 1;
            SD3D.SrcPos(src+(n-1),:) = SD.SrcPos(src+(nid-1),:);
            
        end
        SD3D.Landmarks = SD_POL.Landmarks;
        
        %Plot cropped array
        f2 = figure;
        set(f2,'Name','Final (SUBSET) 3D Array Layout');
        for i = 1:size(SD3D.SrcPos,1)
            plot3(SD3D.SrcPos(i,1),SD3D.SrcPos(i,2),SD3D.SrcPos(i,3),'r.','MarkerSize',30);hold on;
            text(SD3D.SrcPos(i,1),SD3D.SrcPos(i,2),SD3D.SrcPos(i,3),['S' num2str(i)],'Color','r');
        end
        for i = 1:size(SD3D.DetPos,1)
            plot3(SD3D.DetPos(i,1),SD3D.DetPos(i,2),SD3D.DetPos(i,3),'b.','MarkerSize',30);hold on;
            text(SD3D.DetPos(i,1),SD3D.DetPos(i,2),SD3D.DetPos(i,3),['D' num2str(i)],'Color','b');
        end
        plotmesh(SD3D.Landmarks,'g.','MarkerSize', 30);hold on;
        landmarkLabels = {'Nz','Iz','Ar','Al','Cz'};
        for i = 1:size(SD3D.Landmarks,1)
            text(SD3D.Landmarks(i,1),SD3D.Landmarks(i,2)+3,SD3D.Landmarks(i,3)+3,landmarkLabels{i});
        end
        axis equal
        xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
        title('Final Subset Array') 
    else
        SD3D = SD;
        
        % Populate SD3D x, y, and z values with SD_POL x, y, and z values
        for src = 1:nSources
            SD3D.SrcPos(src,:) = SD_POL.SrcPos(src,:);
        end
        
        SD3D.DetPos = SD_POL.DetPos;
        
    end
    
else %Assume 3D contents of layout file remains and save as SD3D
    
    SD3D = SD;
    for n = 1:nNodes
        nid = nodes(n);
        for det = 1:4
            SD3D.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_3d.x;
            SD3D.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_3d.y;
            SD3D.DetPos(det+(n-1)*4,3) = layoutData.docks(nid).optodes(det).coordinates_3d.z;
        end
        src = 1;
        SD3D.SrcPos(src+(n-1),1) = layoutData.docks(nid).optodes(src+4).coordinates_3d.x;
        SD3D.SrcPos(src+(n-1),2) = layoutData.docks(nid).optodes(src+4).coordinates_3d.y;
        SD3D.SrcPos(src+(n-1),3) = layoutData.docks(nid).optodes(src+4).coordinates_3d.z;

    end
    if isfield(layoutData,'Landmarks') %Nasion, Inion, Ar, Al, Cz
        for i = 1:size(layoutData.Landmarks,1)
            SD3D.Landmarks(i,1) = layoutData.Landmarks(i).x;
            SD3D.Landmarks(i,2) = layoutData.Landmarks(i).y;
            SD3D.Landmarks(i,3) = layoutData.Landmarks(i).z;
        end
    end
    SD3DFileName = fullfile(lumoPath, [lumoName '_default.SD3D']);
                                                                                %fprintf(['Saving SD3D to ' SD3DFileName ' ...\n']);   %%%%..GCVL..####
    fprintf('\n Saving SD3D to %s   >>> continuing >>> \n', SD3DFileName);
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

% Convert zero intensity values to noise floor estimate ###################
if size(d,1) == 1 % if only one sample
    mnD = d;
else
    mnD = mean(d);
end
dists_3D = DOTHUB_getSDdists(SD3D);
SDS_noise = 70; % SD channels greater than this (mm) are deemed to be in the the noise floor
n_zeros = sum(d==0,'default'); % was 'all'
if n_zeros > 0
    if max(dists_3D) >= SDS_noise
        noisefloorest = mean(mnD(dists_3D>SDS_noise));
        d(d == 0) = noisefloorest;
        disp(['Warning - zero intensity values have been converted to noise floor estimate'])
    else
        d(d == 0) = 1e-6;
        disp(['Warning - zero intensity values have been converted to 1e-6'])
    end
end

% OUTPUT ##################################################################
nirs.SD = SD;
nirs.SD3D = SD3D;
nirs.aux = aux; 
nirs.d = d;
nirs.d_dark = d_dark;
nirs.s = s;
nirs.t = t;
nirs.CondNames = CondNames;

nirsFileName = fullfile(lumoPath, [lumoName '.nirs']);
                                                                    %fprintf(['Saving file to ' nirsFileName ' ...\n']);     %%%%..GCVL..####
fprintf(' Saving file to %s   >>> continuing >>> \n', nirsFileName);
save(nirsFileName,'-struct','nirs','-v7.3');




% %%%%..GCVL..#### ************************************************
fprintf('\n nSources   nDets  %d %d', nSources, nDets);
fprintf('\n length(SD.Lambda) %d \n', length(SD.Lambda));
nWavs = length(SD.Lambda);
% %%%%..GCVL..#### ************************************************

		displaygraph(t, d, nSources, nDets, nWavs);
end		%% needs an end of main function (DOTHUB_LUMO2nirs_5wav) to allow print functions in same file

%%  print functions

	function prtRawData = printRawData(intDataAll)
		sz=size(intDataAll);
		display=sprintf('%d  ',sz);
		fprintf('\n intDataAll size: %s', display);
		fprintf('\n printing 6 rows all columns but remember prints columns first so after 6 values switches to next column \n');
		fprintf('---------------- first 6 1st col ------------------- | ----------------- first 6 2nd col ----------------- | ------------- first 6 3rd col ----------------- |  ---------->> \n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', intDataAll(1:6,:));

		fprintf('\n Transposed intDataAll *******************************');
		fprintf('\n Data Structure .... ');
		fprintf('\n Notation:  Columns  \n            SourceTile s1 & Detectors d1-4 /  SourceTile s2 & Detectors d5-8  /  Wavelengths 720-760-800-850-890-X \n');
		fprintf('   1        2        3        4        5        6        7        8        9       10       11       12       13       14       15       16       17       18       19       20       21       22       23       24       25       26       27       28       29       30       31       32       33       34       35       36       37       38       39       40       41       42       43       44       45       46       47       48     ');
		fprintf('  49       50       51       52       53       54       55       56       57       58       59       60       61       62       63       64       65       66       67       68       69       70       71       72       73       74       75       76       77       78       79       80       81       82       83       84       85       86       87       88       89       90       91       92       93       94       95       96   \n');
		fprintf('s1d1/720 s1d2/720 s1d3/720 s1d4/720 s1d5/720 s1d6/720 s1d7/720 s1d8/720 s1d1/760 s1d2/760 s1d3/760 s1d4/760 s1d5/760 s1d6/760 s1d7/760 s1d8/760 s1d1/800 s1d2/800 s1d3/800 s1d4/800 s1d5/800 s1d6/800 s1d7/800 s1d8/800 s1d1/850 s1d2/850 s1d3/850 s1d4/850 s1d5/850 s1d6/850 s1d7/850 s1d8/850 s1d1/890 s1d2/890 s1d3/890 s1d4/890 s1d5/890 s1d6/890 s1d7/890 s1d8/890  s1d1/X   s1d2/X   s1d3/X   s1d4/X   s1d5/X   s1d6/X   s1d7/X   s1d8/X   ');
		fprintf('s2d1/720 s2d2/720 s2d3/720 s2d4/720 s2d5/720 s2d6/720 s2d7/720 s2d8/720 s2d1/760 s2d2/760 s2d3/760 s2d4/760 s2d5/760 s2d6/760 s2d7/760 s2d8/760 s2d1/800 s2d2/800 s2d3/800 s2d4/800 s2d5/800 s2d6/800 s2d7/800 s2d8/800 s2d1/850 s2d2/850 s2d3/850 s2d4/850 s2d5/850 s2d6/850 s2d7/850 s2d8/850 s2d1/890 s2d2/890 s2d3/890 s2d4/890 s2d5/890 s2d6/890 s2d7/890 s2d8/890  s2d1/X   s2d2/X   s2d3/X   s2d4/X   s2d5/X   s2d6/X   s2d7/X   s2d8/X \n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', intDataAll(1:6,:).');
		fprintf('\n');
	end

	function prtdark = printdarkchans(chansList, darkInd, d_dark)
		sz=size(chansList);         % chansList has 96 rows and 3 cols    Col-1 rows 1-48=1 49-96=2   Col-2 sequential repeat 1,2,3,4,5,6,7,8   Col-3 8x1, 8x2, 8x3, 8x4, 8x5, 8x6, then repeats
		display=sprintf('%d  ',sz);
		fprintf('\n chansList size: %s', display);
		fprintf('\n chansList    96 rows 3 cols    displayed for ease as 96 cols and 3 rows\n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', chansList(:,:));

		sz=size(darkInd);
		display=sprintf('%d  ',sz);
		fprintf('\n darkInd size: %s', display);
		fprintf('\n darkInd  has 96 rows 1 col   and row number has 0 where LED exists and 1 where no LED exists \n');

		sz=size(d_dark);
		display=sprintf('%d  ',sz);
		fprintf('\n d_dark size: %s', display);
		fprintf('\n d_dark   only columns with no LED first 6 rows\n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n', d_dark(1:6,:).');
	end
	
	function prtd = printnirsd(d, intDataAll)
		sz=size(d);
		display=sprintf('%d  ',sz);
		fprintf('\n "d" size: %s', display);
		fprintf('\n printing 6 rows all columns but remember prints columns first so after 6 values switches to next column \n');
		fprintf('---------------- first 6 1st col ------------------- | ----------------- first 6 2nd col ----------------- | ------------- first 6 3rd col ----------------- |  ---------->> \n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', intDataAll(1:6,:));

		fprintf('\n Transposed d  *******************************');
		fprintf('\n Note: d has fewer columns thab intDataAll since dark columns removed');
		fprintf('\n Data Structure .... HAS CHANGED');
		fprintf('\n Notation:  Columns  \n            SourceTile s1 & Detectors d1-4 /  SourceTile s2 & Detectors d5-8  /  Wavelengths 720-760-800-850-890-X \n');
		fprintf('   1        2        3        4        5        6        7        8        9       10       11       12       13       14       15       16       17       18       19       20       21       22       23       24       25       26       27       28       29       30       31       32       33       34       35       36       37       38       39       40     ');
		fprintf('  41       42       43       44       45       46       47       48       49       50       51       52       53       54       55       56       57       58       59       60       61       62       63       64       65       66       67       68       69       70       71       72       73       74       75       76       77       78       79       80   \n');
		fprintf('s1d1/720 s1d2/720 s1d3/720 s1d4/720 s1d5/720 s1d6/720 s1d7/720 s1d8/720 s2d1/720 s2d2/720 s2d3/720 s2d4/720 s2d5/720 s2d6/720 s2d7/720 s2d8/720 s1d1/760 s1d2/760 s1d3/760 s1d4/760 s1d5/760 s1d6/760 s1d7/760 s1d8/760 s2d1/760 s2d2/760 s2d3/760 s2d4/760 s2d5/760 s2d6/760 s2d7/760 s2d8/760 s1d1/800 s1d2/800 s1d3/800 s1d4/800 s1d5/800 s1d6/800 s1d7/800 s1d8/800 ');
		fprintf('s2d1/800 s2d2/800 s2d3/800 s2d4/800 s2d5/800 s2d6/800 s2d7/800 s2d8/800 s1d1/850 s1d2/850 s1d3/850 s1d4/850 s1d5/850 s1d6/850 s1d7/850 s1d8/850 s2d1/850 s2d2/850 s2d3/850 s2d4/850 s2d5/850 s2d6/850 s2d7/850 s2d8/850 s1d1/890 s1d2/890 s1d3/890 s1d4/890 s1d5/890 s1d6/890 s1d7/890 s1d8/890 s2d1/890 s2d2/890 s2d3/890 s2d4/890 s2d5/890 s2d6/890 s2d7/890 s2d8/890 \n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', d(1:6,:).');

	end

%%  display graphs functions

	function showgraph = displaygraph(time, d, nSources, nDets, nWavs)
 		set(groot, 'DefaultFigureWindowStyle','docked');  %docks graphic windows in matlab screen
		%initialise
		rawd_graph = d;
		fprintf('\n rawd_graph \n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',rawd_graph(1:6,:).');
		fprintf('\n');

		s1d1 = 1:16:66;
		s1d2 = 2:16:67;
		s1d3 = 3:16:68;
		s1d4 = 4:16:69;
		s1d5 = 5:16:70;
		s1d6 = 6:16:71;
		s1d7 = 7:16:72;
		s1d8 = 8:16:73;
		s2d1 = 9:16:74;
		s2d2 = 10:16:75;
		s2d3 = 11:16:76;
		s2d4 = 12:16:77;
		s2d5 = 13:16:78;
		s2d6 = 14:16:79;
		s2d7 = 15:16:80;
		s2d8 = 16:16:81;
 
		%%%% Raw Signal Figures	
			figure('Name','Raw Int','NumberTitle','off');
			plot(time, rawd_graph);
			title('Raw Signals I(t)');
			xlabel('Time (s)');

		%%%% Raw Signal By Detector
			figure('Name','Raw Int A','NumberTitle','off');
			t = tiledlayout(2,2);
			t.Title.String = 'Raw Signal I(t)';
			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
			t.Subtitle.FontSize = 10;
				nexttile(t);
				plot(time, rawd_graph(:,s1d1,':'));
				title('s1d1');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s1d3,':'));
				title('s1d3');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s1d2,':'));
				title('s1d2');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s1d4,':'));
				title('s1d4');
				xlabel('Time (s)');

			figure('Name','Raw Int B','NumberTitle','off');
			t = tiledlayout(2,2);
			t.Title.String = 'Raw Signal I(t)';
			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
			t.Subtitle.FontSize = 10;
			nexttile(t);
				plot(time, rawd_graph(:,s1d5,':'));
				title('s1d5');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s1d7,':'));
				title('s1d7');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s1d6,':'));
				title('s1d6');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s1d8,':'));
				title('s1d8');
				xlabel('Time (s)');
 
			figure('Name','Raw Int C','NumberTitle','off');
			t = tiledlayout(2,2);
			t.Title.String = 'Raw Signal I(t)';
			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
			t.Subtitle.FontSize = 10;
				nexttile(t);
				plot(time, rawd_graph(:,s2d1,':'));
				title('s2d1');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s2d3,':'));
				title('s2d3');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s2d2,':'));
				title('s2d2');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s2d4,':'));
				title('s2d4');
				xlabel('Time (s)');
 
			figure('Name','Raw Int D','NumberTitle','off');
			t = tiledlayout(2,2);
			t.Title.String = 'Raw Signal I(t)';
			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
			t.Subtitle.FontSize = 10;
				nexttile(t);
				plot(time, rawd_graph(:,s2d5,':'));
				title('s2d5');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s2d7,':'));
				title('s2d7');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s2d6,':'));
				title('s2d6');
				xlabel('Time (s)');
				nexttile(t);
				plot(time, rawd_graph(:,s2d8,':'));
				title('s2d8');
				xlabel('Time (s)');
 
 
			figure('Name','Raw Int E','NumberTitle','off');
			bar1125 = [rawd_graph(1,s1d1) 0; rawd_graph(1,s1d2) 0; rawd_graph(1,s1d3) 0; rawd_graph(1,s1d4) 0; rawd_graph(1,s2d5) 0; rawd_graph(1,s2d6) 0; rawd_graph(1,s2d7) 0; rawd_graph(1,s2d8) 0];
			xnames1125 = categorical({'s1d1','s1d2','s1d3','s1d4','s2d5','s2d6','s2d7','s2d8'});
			xnames1125 = reordercats(xnames1125,{'s1d1','s1d2','s1d3','s1d4','s2d5','s2d6','s2d7','s2d8'});
			bar(xnames1125, bar1125);
			title('Raw Signal Strength I(0)');

			figure('Name','Raw Int F','NumberTitle','off');
			bar2115 = [rawd_graph(1,s2d1) 0; rawd_graph(1,s2d2) 0; rawd_graph(1,s2d3) 0; rawd_graph(1,s2d4) 0; rawd_graph(1,s1d5) 0; rawd_graph(1,s1d6) 0; rawd_graph(1,s1d7) 0; rawd_graph(1,s1d8) 0];
			xnames2115 = categorical({'s2d1','s2d2','s2d3','s2d4','s1d5','s1d6','s1d7','s1d8'});
			xnames2115 = reordercats(xnames2115,{'s2d1','s2d2','s2d3','s2d4','s1d5','s1d6','s1d7','s1d8'});
			bar(xnames2115, bar2115);
			title('Raw Signal Strength I(0)');

		%%%% Modified Raw Ratios	AND   Raw Attenuation	
				%% look at raw ratios
				%% modify (expand signal *blowup) and add nW (wavelength number) to separate on graph
% 				%% this can be a long loop time and could be faster in main body but put here in figure section so can ignore if wish
% 				blowup = 2;
% 				rows = height(d);
% 				for ii = 1:rows
% 					for nW = 1:nWavs
% 						for nS = 1:nSources
% 							for nD = 1:nDets
% 								jj = (nW-1)*16 + (nS-1)*8 + nD;
% 								ratio(ii,jj) = d(1,jj)/d(ii,jj); 
% 								rawratio(ii,jj) = blowup*(ratio(ii,jj) - 1) + nW;  
% 								log10ratio(ii,jj) = log10(ratio(ii,jj));
% 							end
% 						end
% 					end
% 				end
% 
% 				
% 		fprintf('\n Transposed rawratio 6 rows  ******************************* \n');
% 		fprintf('   1        2        3        4        5        6        7        8        9       10       11       12       13       14       15       16       17       18       19       20       21       22       23       24       25       26       27       28       29       30       31       32       33       34       35       36       37       38       39       40     ');
% 		fprintf('  41       42       43       44       45       46       47       48       49       50       51       52       53       54       55       56       57       58       59       60       61       62       63       64       65       66       67       68       69       70       71       72       73       74       75       76       77       78       79       80   \n');
% 		fprintf('s1d1/720 s1d2/720 s1d3/720 s1d4/720 s1d5/720 s1d6/720 s1d7/720 s1d8/720 s2d1/720 s2d2/720 s2d3/720 s2d4/720 s2d5/720 s2d6/720 s2d7/720 s2d8/720 s1d1/760 s1d2/760 s1d3/760 s1d4/760 s1d5/760 s1d6/760 s1d7/760 s1d8/760 s2d1/760 s2d2/760 s2d3/760 s2d4/760 s2d5/760 s2d6/760 s2d7/760 s2d8/760 s1d1/800 s1d2/800 s1d3/800 s1d4/800 s1d5/800 s1d6/800 s1d7/800 s1d8/800 ');
% 		fprintf('s2d1/800 s2d2/800 s2d3/800 s2d4/800 s2d5/800 s2d6/800 s2d7/800 s2d8/800 s1d1/850 s1d2/850 s1d3/850 s1d4/850 s1d5/850 s1d6/850 s1d7/850 s1d8/850 s2d1/850 s2d2/850 s2d3/850 s2d4/850 s2d5/850 s2d6/850 s2d7/850 s2d8/850 s1d1/890 s1d2/890 s1d3/890 s1d4/890 s1d5/890 s1d6/890 s1d7/890 s1d8/890 s2d1/890 s2d2/890 s2d3/890 s2d4/890 s2d5/890 s2d6/890 s2d7/890 s2d8/890 \n');
% 		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', rawratio(1:6,:).');
% 
% 		fprintf('\n Transposed log10ratio 6 rows  ******************************* \n');
% 		fprintf('   1        2        3        4        5        6        7        8        9       10       11       12       13       14       15       16       17       18       19       20       21       22       23       24       25       26       27       28       29       30       31       32       33       34       35       36       37       38       39       40     ');
% 		fprintf('  41       42       43       44       45       46       47       48       49       50       51       52       53       54       55       56       57       58       59       60       61       62       63       64       65       66       67       68       69       70       71       72       73       74       75       76       77       78       79       80   \n');
% 		fprintf('s1d1/720 s1d2/720 s1d3/720 s1d4/720 s1d5/720 s1d6/720 s1d7/720 s1d8/720 s2d1/720 s2d2/720 s2d3/720 s2d4/720 s2d5/720 s2d6/720 s2d7/720 s2d8/720 s1d1/760 s1d2/760 s1d3/760 s1d4/760 s1d5/760 s1d6/760 s1d7/760 s1d8/760 s2d1/760 s2d2/760 s2d3/760 s2d4/760 s2d5/760 s2d6/760 s2d7/760 s2d8/760 s1d1/800 s1d2/800 s1d3/800 s1d4/800 s1d5/800 s1d6/800 s1d7/800 s1d8/800 ');
% 		fprintf('s2d1/800 s2d2/800 s2d3/800 s2d4/800 s2d5/800 s2d6/800 s2d7/800 s2d8/800 s1d1/850 s1d2/850 s1d3/850 s1d4/850 s1d5/850 s1d6/850 s1d7/850 s1d8/850 s2d1/850 s2d2/850 s2d3/850 s2d4/850 s2d5/850 s2d6/850 s2d7/850 s2d8/850 s1d1/890 s1d2/890 s1d3/890 s1d4/890 s1d5/890 s1d6/890 s1d7/890 s1d8/890 s2d1/890 s2d2/890 s2d3/890 s2d4/890 s2d5/890 s2d6/890 s2d7/890 s2d8/890 \n');
% 		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', log10ratio(1:6,:).');

% 				
% 			figure('Name','Mod-Ratio','NumberTitle','off');
% 			plot(time, rawratio, ':');
% 			title('All Modified Raw Ratios ~I_0/I');
% 
% 			figure('Name','Mod-Ratio A','NumberTitle','off');
% 			sgtitle('Modified Raw Ratio ~I_0/I');
% 			subplot(2,2,1);
% 			plot(time, rawratio(:,s1d1,':'));
% 			title('s1d1');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,3);
% 			plot(time, rawratio(:,s1d2,':'));
% 			title('s1d2');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,2);
% 			plot(time, rawratio(:,s1d3,':'));
% 			title('s1d3');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,4);
% 			plot(time, rawratio(:,s1d4,':'));
% 			title('s1d4');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 
% 			figure('Name','Mod-Ratio B','NumberTitle','off');
% 			sgtitle('Modified Raw Ratio ~I_0/I');
% 			subplot(2,2,1);
% 			plot(time, rawratio(:,s1d5,':'));
% 			title('s1d5');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,3);
% 			plot(time, rawratio(:,s1d6,':'));
% 			title('s1d6');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,2);
% 			plot(time, rawratio(:,s1d7,':'));
% 			title('s1d7');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,4);
% 			plot(time, rawratio(:,s1d8,':'));
% 			title('s1d8');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 
% 			figure('Name','Mod-Ratio C','NumberTitle','off');
% 			sgtitle('Modified Raw Ratio ~I_0/I');
% 			subplot(2,2,1);
% 			plot(time, rawratio(:,s2d1,':'));
% 			title('s2d1');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,3);
% 			plot(time, rawratio(:,s2d2,':'));
% 			title('s2d2');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,2);
% 			plot(time, rawratio(:,s2d3,':'));
% 			title('s2d3');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,4);
% 			plot(time, rawratio(:,s2d4,':'));
% 			title('s2d4');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 
% 			figure('Name','Mod-Ratio D','NumberTitle','off');
% 			sgtitle('Modified Raw Ratio ~I_0/I');
% 			subplot(2,2,1);
% 			plot(time, rawratio(:,s2d5,':'));
% 			title('s2d5');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,3);
% 			plot(time, rawratio(:,s2d6,':'));
% 			title('s2d6');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,2);
% 			plot(time, rawratio(:,s2d7,':'));
% 			title('s2d7');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});
% 			subplot(2,2,4);
% 			plot(time, rawratio(:,s2d8,':'));
% 			title('s2d8');
% 			yticks([1 2 3 4 5]);
% 			yticklabels({'720','760','800','850','890'});

			
% 		%%%%  Raw Attenuation Figures
% 		
% 			figure('Name','Raw Atten','NumberTitle','off');
% 			plot(time, log10ratio, ':');
% 			title('All Raw Attenuations log_1_0(I_0/I)');
% 		
% 		
% 			figure('Name','Raw Atten A','NumberTitle','off');
% 			t = tiledlayout(2,2);
% 			t.Title.String = 'Raw Attenuation log_1_0(I_0/I)';
% 			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
% 			t.Subtitle.FontSize = 10;
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d1,':'));
% 				title('s1d1');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d3,':'));
% 				title('s1d3');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d2,':'));
% 				title('s1d2');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d4,':'));
% 				title('s1d4');
% 				xlabel('Time (s)');
% 			
% 			figure('Name','Raw Atten B','NumberTitle','off');
% 			t = tiledlayout(2,2);
% 			t.Title.String = 'Raw Attenuation log_1_0(I_0/I)';
% 			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
% 			t.Subtitle.FontSize = 10;
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d5,':'));
% 				title('s1d5');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d7,':'));
% 				title('s1d7');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d6,':'));
% 				title('s1d6');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s1d8,':'));
% 				title('s1d8');
% 				xlabel('Time (s)');
% 
% 			figure('Name','Raw Atten C','NumberTitle','off');
% 			t = tiledlayout(2,2);
% 			t.Title.String = 'Raw Attenuation log_1_0(I_0/I)';
% 			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
% 			t.Subtitle.FontSize = 10;
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d1,':'));
% 				title('s2d1');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d3,':'));
% 				title('s2d3');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d2,':'));
% 				title('s2d2');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d4,':'));
% 				title('s2d4');
% 				xlabel('Time (s)');
% 
% 			figure('Name','Raw Atten D','NumberTitle','off');
% 			t = tiledlayout(2,2);
% 			t.Title.String = 'Raw Attenuation log_1_0(I_0/I)';
% 			t.Subtitle.String = 'blu=720 org=760 yel=800 pur=850 grn=890 nm';
% 			t.Subtitle.FontSize = 10;
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d5,':'));
% 				title('s2d5');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d7,':'));
% 				title('s2d7');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d6,':'));
% 				title('s2d6');
% 				xlabel('Time (s)');
% 				nexttile(t);
% 				plot(time, log10ratio(:,s2d8,':'));
% 				title('s2d8');
% 				xlabel('Time (s)');
% 		
 			end
			

	% %%%%..GCVL..#### ************************************************


