function [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUFR2nirs(lufrFileName,layoutFileName,posCSVFileName,eventsCSVFileName)

%This script reads a .LUFR file and also loads the associated .JSON layout from a seperate file.

%If subject-specific registration data (from polhemus or other method) is available
%in .cvs format, this too can be parsed. If it is parsed, the function
%DOTHUB_LUMOPolhemus2SD is called to determine optode positions from the
%polhemus data. The resulting 3D optode positions overwrite the default 3D
%data from the .JSON file.

%This function then outputs the .nirs variables and produces a .nirs file with 
%the same name as the .LUFR. Note re-running this function will overwrite 
%previous .nirs derived from the parsed LUFR file. Note that the .nirs
%file contains both 2D (SD) and 3D (SD_3D) info.

%######################## INPUTS ##########################################

% lufrFileName      :   The filename of the .LUMO data directory (if not parsed, requested)

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

% #########################################################################
% RJC, UCL, October 2021
%
% ############################# Updates ###################################
% #########################################################################

% #########################################################################
% TO DO:
% #########################################################################

% MANAGE VARIABLES  #######################################################
if ~exist('lufrFileName','var')
    disp('Select .lufr file...');
    [lufrFileName,lufrPath] = uigetfile({'*.lufr';'*.LUFR'},'Select .lufr file...');
    lufrFileName = fullfile(lufrPath,lufrFileName);
elseif isempty(lufrFileName)
    disp('Select .lufr file...');
    [lufrFileName,lufrPath] = uigetfile({'*.lufr';'*.LUFR'},'Select .lufr file...');
    lufrFileName = fullfile(lufrPath,lufrFileName);
elseif ~exist(lufrFileName,'file')
    disp('Specified lufr file not found, please select...');
    [lufrFileName,lufrPath] = uigetfile({'*.lufr';'*.LUFR'},'Select .lufr file...');
    lufrFileName = fullfile(lufrPath,lufrFileName);
end

%Ensure full path
[lufrPath,lufrName,ext] = fileparts(lufrFileName);
if ~strcmpi(ext,'.lufr')
    error('The lufrFileName input must point to a .lufr file');
end
if isempty(lufrPath)
    lufrFileName = fullfile(pwd,lufrFileName);
    [lufrPath,lufrName,~] = fileparts(lufrFileName);
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
[infoblks, ... % Free-form information fields
          enum, ...     % JSON format enumeration from the back end
          t_frm, ...    % The time increment of a single data frame (1/fps)
          chdat, ...    % Channel data (channels x frames)
          dkdat, ...    % Dark channel data (often not recorded) (dark channels x frames)
          tmpdat, ...   % Tile internal temperatures (tiles x frames)
          vindat, ...   % Tile input voltages (tiles x frames)
          dcmax, ...    % Detector maximum values for saturation (tils x dets x frames)
          srcpwr, ...   % Source powers (nodes x wavelengths x frames)
          evtim, ...    % Time of each event (ms)
          evstr, ...    % Associated string for each marked event
          t_frmmpu, ... % The time increment of a single MPU frame
          gyrdat, ...   % Gyroscope data (nodes x dim x meas/frame x frame)
          accdat] ...   % Accellerometer data (nodes x dim x meas/frame x frame)
          = loadlufr(lufrFileName);

%Unpacking
fs = 1/t_frm;
fs_mpu = 1/t_frmmpu;
      
% Determine Status of Layout Data
if ~exist('layoutFileName','var') %Not parsed so check exists in .lufr
    if ~exist('enum.groups.index.layout', 'var') %Not contained in .lufr, so load
        disp('Layout file not parsed or found in lufr, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        disp('Using default layout within .lufr directory');
        layoutFileName = [];
    end
elseif isempty(layoutFileName)
    if ~exist('enum.groups.index.layout', 'var') %Not contained in .lufr, so load
        disp('Layout file not parsed or found in lufr, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        disp('Using default layout within .lufr directory');
        layoutFileName = [];
    end
end

% Convert to .nirs format #########################################################################

% Intensity data
d = chdat';

% Temporal data
t = [0:1/fs:(size(chdat, 2) - 1)/fs]'; 

% SD info
nNodes = size(enum.groups.nodes, 1);
nodes = [enum.groups.nodes.id];
SD.Lambda = unique([enum.groups.channels.src_wl]);
SD.nSrcs = nNodes * 3;
SD.nDets = nNodes * 4;

% Get 2D Layout info
if isempty(layoutFileName) %This means no layout file parsed - must use data in lufr (should be there)
    for n = 1:nNodes
    nid = nodes(n);
        for det = 1:4
            SD.DetPos(det+(n-1)*4,1) = enum.groups.index.layout.docks(nid).optodes(det).coordinates_2d.x;
            SD.DetPos(det+(n-1)*4,2) = enum.groups.index.layout.docks(nid).optodes(det).coordinates_2d.y;
            SD.DetPos(det+(n-1)*4,3) = 0;
        end
        for src = 1:3
            SD.SrcPos(src+(n-1)*3,1) = enum.groups.index.layout.docks(nid).optodes(src+4).coordinates_2d.x;
            SD.SrcPos(src+(n-1)*3,2) = enum.groups.index.layout.docks(nid).optodes(src+4).coordinates_2d.y;
            SD.SrcPos(src+(n-1)*3,3) = 0;
        end
    end
else %Load specified layout file
    %Load layout JSON
    layoutData = jsondecode(fileread(layoutFileName));
    nDocks = size(layoutData.docks,1);%Is this info stored in lufr too?
    for n = 1:nNodes
        nid = nodes(n);
        for det = 1:4
            SD.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_2d.x;
            SD.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_2d.y;
            SD.DetPos(det+(n-1)*4,3) = 0;
        end
        for src = 1:3
            SD.SrcPos(src+(n-1)*3,1) = layoutData.docks(nid).optodes(src+4).coordinates_2d.x;
            SD.SrcPos(src+(n-1)*3,2) = layoutData.docks(nid).optodes(src+4).coordinates_2d.y;
            SD.SrcPos(src+(n-1)*3,3) = 0;
        end
    end
    if isfield(layoutData,'Landmarks')
        for i = 1:size(layoutData.Landmarks,1)
            SD.Landmarks(i,1) = layoutData.Landmarks(i).x;
            SD.Landmarks(i,2) = layoutData.Landmarks(i).y;
            SD.Landmarks(i,3) = layoutData.Landmarks(i).z;
        end
    end
end
   
SD.SpatialUnit = 'mm';

%########## Output node-wise peripheral datatypes ###########:

%Acc data
periph.Acc = accdat;
periph.Gyro = gyrdat;
periph.t_mpu = (0:t_frmmpu:(size(accdat, 3) - 1)*t_frmmpu)'; 

% Source power structure - source power logged for each node, wavelength and time point
periph.SrcPowers = srcpwr;

% Max DC voltage on detectors throughout measurement
periph.DCMax = dcmax;

% Temperature of each node
periph.NodeTemp = tmpdat;

% Input voltage
periph.InputVolt = vindat;

% Measurement list
SD.MeasList = ones(SD.nSrcs * SD.nDets * length(SD.Lambda), 4);

% Creating temporary object to change indexing of sources - keep 1:3 instead of 1:6
enum_temp = enum.groups.channels;
change_index = find([enum_temp.src_wl] == 850);
for i = change_index
    enum_temp(i).src_idx = enum_temp(i).src_idx - 3;
end

% Build MeasList array (3rd column unused)
for i = 1:size(chdat, 1)
    SD.MeasList(i, 1) = 3 * enum_temp(i).src_node_id - (3 - enum_temp(i).src_idx - 1);
    SD.MeasList(i, 2) = 4 * enum.groups.channels(i).det_node_id - (4 - enum.groups.channels(i).det_idx - 1); 
    SD.MeasList(i, 4) = find(enum.groups.channels(i).src_wl == SD.Lambda);
end

% Sort MeasList by wavelength, update d accordingly
[SD.MeasList, sort_indx] = sortrows(SD.MeasList, 4);
d = double(d(:,sort_indx));

% Force create MeasListAct - all ones
SD.MeasListAct = ones(size(SD.MeasList,1),1);

% Pad events
s = zeros(size(t));

% Pad aux
aux = zeros(size(t,1),8);

% Sort 3D SD information
if polhemusFlag %If polhemus information is parsed, calculate S-D positions from that file
    fprintf(['Using ' posCSVFileName ' to define SD3D...\n']);
    [SD_POL, SD3DFileName] = DOTHUB_LUMOpolhemus2SD3D(posCSVFileName); %This line saves the .SD3D
    SD_POL.MeasListAct = SD.MeasListAct;
    
    if nNodes<nDocks %Crop SD file accordingly if the number of docks populated in the datafile differs from the number in the SD_3D data
        SD3D = SD; %Define based on SD which contains correct measlist
        for n = 1:nNodes
            nid = nodes(n);
            for det = 1:4
                SD3D.DetPos(det+(n-1)*4,:) = SD_POL.DetPos(det+(nid-1)*4,:);
            end
            for src = 1:3
                SD3D.SrcPos(src+(n-1)*3,:) = SD_POL.SrcPos(src+(nid-1)*3,:);
            end
        end
        SD3D.Landmarks = SD_POL.Landmarks;
        
        %Plot cropped array
        f2 = figure;
        set(f2,'Name','Final (SUBSET) 3D Array Layout');
        for i = 1:size(SD3D.SrcPos,1)
            plot3(SD3D.SrcPos(i,1),SD3D.SrcPos(i,2),SD3D.SrcPos(i,3),'r.','MarkerSize',30);hold on;
            text(SD3D.SrcPos(i,1),SD3D.SrcPos(i,2)+3,SD3D.SrcPos(i,3),['S' num2str(i)],'Color','r');
        end
        for i = 1:size(SD3D.DetPos,1)
            plot3(SD3D.DetPos(i,1),SD3D.DetPos(i,2),SD3D.DetPos(i,3),'b.','MarkerSize',30);hold on;
            text(SD3D.DetPos(i,1),SD3D.DetPos(i,2)+3,SD3D.DetPos(i,3),['D' num2str(i)],'Color','b');
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
        SD3D = SD_POL;
    end
    
else %Assume 3D contents of layout file (or lufr) is to be used and save as SD3D
    
    SD3D = SD;
    if isempty(layoutFileName) %This means no layout file parsed - must use data in lufr (should be there)
        
        for n = 1:nNodes
            nid = nodes(n);
            for det = 1:4
                SD3D.DetPos(det+(n-1)*4,1) = enum.groups.index.layout.docks(nid).optodes(det).coordinates_3d.x;
                SD3D.DetPos(det+(n-1)*4,2) = enum.groups.index.layout.docks(nid).optodes(det).coordinates_3d.y;
                SD3D.DetPos(det+(n-1)*4,3) = enum.groups.index.layout.docks(nid).optodes(det).coordinates_3d.z;
            end
            for src = 1:3
                SD3D.SrcPos(src+(n-1)*3,1) = enum.groups.index.layout.docks(nid).optodes(src+4).coordinates_3d.x;
                SD3D.SrcPos(src+(n-1)*3,2) = enum.groups.index.layout.docks(nid).optodes(src+4).coordinates_3d.y;
                SD3D.SrcPos(src+(n-1)*3,3) = enum.groups.index.layout.docks(nid).optodes(src+4).coordinates_3d.z;
            end
        end
        
    else %use layout file
        
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
        if isfield(layoutData,'Landmarks')
            for i = 1:size(layoutData.Landmarks,1)
                SD3D.Landmarks(i,1) = layoutData.Landmarks(i).x;
                SD3D.Landmarks(i,2) = layoutData.Landmarks(i).y;
                SD3D.Landmarks(i,3) = layoutData.Landmarks(i).z;
            end
        end
    end

    SD3DFileName = fullfile(lufrPath, [lufrName '_default.SD3D']);
    fprintf(['Saving SD3D to ' SD3DFileName ' ...\n']);
    save(SD3DFileName,'SD3D');
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

    %Events hack
    if isempty(evtim)
        disp('No events found.');
        s = zeros(size(t,1),1);
        CondNames = {''};
    else
        %Filter out newlines and data start/stop characters
        evtim_filt = evtim(~strcmpi(evstr,newline) & ~strcmpi(evstr,'$'))*1e-3; %To seconds
        evstr_filt = evstr(~strcmpi(evstr,newline) & ~strcmpi(evstr,'$'));
        evList = unique(evstr_filt,'stable');
        nCond = length(evList);
        s = zeros(length(t),nCond);
        
        for i = 1:length(evtim_filt)
            [~,ind] = min(abs(t - evtim_filt(i)));
            cond = find(strcmpi(evList,evstr_filt(i)));
            s(ind,cond) = 1;
        end
        CondNames = evList; %Check if this is a cell array as is expected.
    end
end

% % Convert zero intensity values to noise floor estimate ###################
% if size(d,1) == 1 % if only one sample
%     mnD = d;
% else
%     mnD = mean(d);
% end
% dists_3D = DOTHUB_getSDdists(SD3D);
% SDS_noise = 70; % SD channels greater than this (mm) are deemed to be in the the noise floor
% n_zeros = sum(d==0,'all');
% if n_zeros > 0
%     if max(dists_3D) >= SDS_noise
%         noisefloorest = mean(mnD(dists_3D>SDS_noise));
%         d(d == 0) = noisefloorest;
%         disp(['Warning - zero intensity values have been converted to noise floor estimate'])
%     else
%         d(d == 0) = 1e-6;
%         disp(['Warning - zero intensity values have been converted to 1e-6'])
%     end
% end

% OUTPUT ##################################################################
nirs.SD = SD;
nirs.SD3D = SD3D;
nirs.periph = periph;
nirs.aux = aux;
nirs.d = d;
nirs.s = s;
nirs.t = t;
nirs.CondNames = CondNames;

nirsFileName = fullfile(lufrPath, [lufrName '.nirs']);
fprintf(['Saving file to ' nirsFileName ' ...\n']);
save(nirsFileName,'-struct','nirs','-v7.3');

end


function [infoblks, ... % Free-form information fields
          enum, ...     % JSON format enumeration from the back end
          t_frm, ...    % The time increment of a single data frame (1/fps)
          chdat, ...    % Channel data (channels x frames)
          dkdat, ...    % Dark channel data (often not recorded) (dark channels x frames)
          tmpdat, ...   % Tile internal temperatures (tiles x frames)
          vindat, ...   % Tile input voltages (tiles x frames)
          dcmax, ...    % Detector maximum values for saturation (tils x dets x frames)
          srcpwr, ...   % Source powers (nodes x wavelengths x frames)
          evtim, ...    % Time of each event (ms)
          evstr, ...    % Associated string for each marked event
          t_frmmpu, ... % The time increment of a single MPU frame
          gyrdat, ...   % Gyroscope data (nodes x dim x meas/frame x frame), units of degrees per second
          accdat] ...   % Accellerometer data (nodes x dim x meas/frame x frame), units of g
          = loadlufr(fn)

%S Powell, Z Kovacsova, Gowerlabs Ltd. 2021

    % Some constants for version 1
    tag_information = 1;
    tag_enumeration = 2;
    tag_frame = 3;
    tag_event = 4;
    
    % Open the file
    if(~isfile(fn))
        error('The lufr file %s cannot be found, get lost', fn);
    end
      
    fid = fopen(fn, 'rb');
    
    % Get file size
    fseek(fid, 0, 'eof');
    filesize = ftell(fid);
    fseek(fid, 0, 'bof');
    
    % Check the file is valid
    filehdr = string(fread(fid, 2, 'int8=>char').');
    if(filehdr ~= "LU")
        error('The lumo frame file does not contain the correct file header');
    end
    
    % Check the file version
    filever = fread(fid, 1, 'uint8=>double');
    if(filever ~= 1)
        error('The lumo frame file is of an unknown version %i', filever);
    end
        
    % Skip over zeros
    fseek(fid, 3, 'cof');
    
    % Get endieness marker
    endmarker = fread(fid, 1, 'uint16=>uint16');
    if(endmarker ~= 0x1234)
        error('The lumo frame file endienness is not supported');
    end
    
    % Read the record counter
    record_count = fread(fid, 1, 'uint32=>double');
    if(record_count == 0)
        warning('Record counter not set, file may be corrupted');
    end
    
    % Now, seek through the file to count the number of frames, storing
    % the offset in an array alongside the length
    rcoffset = [];
    rclength = [];
    rctag = [];
    
    % Also enumerate the number of particular record types
    n_frames = 0;
    n_enums = 0;
    n_events = 0;
    n_infos = 0;
    
    wb = waitbar(ftell(fid)/filesize, 'Scanning records');
    i = 0;
    
    while(~feof(fid))
                
        % Check for a tag
        recordtag = fread(fid, 1, 'uint32=>double');
        if(feof(fid))
            break;
        end
        
        if(recordtag > 4)
           if(feof(fid))
                break;
            end
            error('The lumo frame file contains a bad frame header');
        end
                
        % Get the length
        recordlen = fread(fid, 1, 'uint32=>double');
        
        % Record it
        rcoffset(end+1) = ftell(fid);
        rclength(end+1) = recordlen; 
        rctag(end+1) = recordtag;
        
        % Check that standard frame record lengths are constant
        if(recordtag == tag_frame)

            sizeparam = fread(fid, 9, 'int32=>double');
            if n_frames == 0
                % On the first run, we store this as a reference  
                sizeparam_ref = sizeparam;
                firstframe = false;
            else
                % Subsequently, check no changes have happned
                if ~all(sizeparam(2:end) == sizeparam_ref(2:end))
                    error('Frame data varies, I cannot handle it');
                end

            end
            
            n_frames = n_frames + 1;

        end
        
        if(recordtag == tag_enumeration) 
            n_enums = n_enums + 1;
        end
        
        if(recordtag == tag_information)
            n_infos = n_infos + 1;
        end
        
        if(recordtag == tag_event)
            n_events = n_events + 1;
        end
        
        % Skip to the next
        fseek(fid, rcoffset(end) + recordlen, 'bof');
        
        % Increment counter
        i = i+1;
        if ~mod(i, 100)
            waitbar(ftell(fid)/filesize, wb, sprintf('Scanning frames %d',length(rclength)));
        end
       
    end
    close(wb);
    
    if(n_enums < 1)
        warning('File does not contain an enumeration block');
        enum = [];
    end
    
    if(n_enums > 1)
        warning('File contains more than one enumeration block, only the last will be returned');
    end
            
    if length(rclength) == record_count           
        fprintf('File contains %d records, %d frames, %d events \n', length(rclength), n_frames, n_events);
    else
        warning('Found %d records, %d frames, %d events in the input data file\n', length(rclength), n_frames, n_events);
    end
    
    % Compute size of the output data an allocate
    n_nodes  = sizeparam_ref(3);
    n_schans = sizeparam_ref(4);
    n_dchans = sizeparam_ref(5);
    n_mpu    = sizeparam_ref(6);
    n_spw    = sizeparam_ref(7);
    n_det    = sizeparam_ref(8);
    n_row    = sizeparam_ref(9);
    
    % Get frames per second
    t_frm = sizeparam_ref(2)*1e-6; %in seconds
    fps = 1/(sizeparam_ref(2)*1e-6);
    
    fprintf('Frame rate is %.2f fps\n', fps);
        
    chdat = zeros(n_schans, n_frames, 'single');        % Channel data
    dkdat = zeros(n_dchans, n_frames, 'single');        % Dark channel data
    tmpdat = zeros(n_nodes, n_frames, 'single');        % Temperature 
    vindat = zeros(n_nodes, n_frames, 'single');        % Input voltage
    srcpwr = zeros(n_nodes, n_spw, n_frames, 'single'); % Source powers
    dcmax = zeros(n_nodes, n_det, n_frames, 'single');  % DC max (per det)
    accdat = zeros(n_nodes, 3, n_mpu, n_frames);        % MPU data (concatenated later)
    gyrdat = zeros(n_nodes, 3, n_mpu, n_frames);        % MPU data (concatenated later)
    
    % Allocate information blocks
    infoblks = strings(n_infos, 1);
    
    % Allocate event times and strings
    evtim = zeros(n_events, 1);
    evstr = strings(n_events, 1);
    
    % Loop over the frames and get them data    
    wb = waitbar(0, 'Processing frames');
    
    i_fr = 1;   % Frame count
    i_if = 1;   % Information block count
    i_ev = 1;   % Event count
    
    for i = 1:length(rclength)
        
        % Seek to the frame
         fseek(fid, rcoffset(i), 'bof');

         % Standard measurement frame
         if rctag(i) == tag_frame

             % Skip over the frame dimension data
             fread(fid, 9, 'int32=>double');

             % Get the data
             chdat(:,i_fr) = fread(fid, n_schans, 'single=>single');
             dkdat(:,i_fr) = fread(fid, n_dchans, 'single=>single');

             % Skip over the MPU for now
             mpuframe = fread(fid, n_mpu*9*n_nodes, 'single=>single');
             mpudat_raw = permute(reshape(mpuframe, 9, n_mpu, n_nodes), [3 1 2]);
             
             for j = 1:n_mpu
                 % mpudat_raw : [n_nodes, 9, n_mpu]
                 gyrdat(:, :, j, i_fr) = mpudat_raw(:, 1:3, j);
                 accdat(:, :, j, i_fr) = mpudat_raw(:, 4:6, j);
             end
            
             % Gerrabit more      
             tmpdat(:,i_fr) = fread(fid, n_nodes, 'single=>single');
             vindat(:,i_fr) = fread(fid, n_nodes, 'single=>single');

             % Source powers
             srcpwr(:,:,i_fr) = reshape(fread(fid, n_spw*n_nodes, 'single=>single'), n_spw, n_nodes).';

             % DX maximum per detector on each frame
             dcmax_byrow = reshape(fread(fid, n_det * n_row * n_nodes, 'single=>single'), n_det, n_row, n_nodes);
             dcmax(:,:,i_fr) = squeeze(max(dcmax_byrow, [], 2)).';

             i_fr = i_fr+1;
         end
         
         % Enumeration
         if(rctag(i) == tag_enumeration)
             enumjson = fread(fid, rclength(i), 'char=>char');
             enum = jsondecode(convertCharsToStrings(enumjson));
         end
         
         if(rctag(i) == tag_event)
            evtim(i_ev) = fread(fid, 1, 'uint32=>double');
            evtxt = fread(fid, rclength(i)-4, 'char=>char');
            evstr(i_ev) = convertCharsToStrings(evtxt);
            i_ev = i_ev + 1;
             
         end
         
         if(rctag(i) == tag_information)
             infotxt = fread(fid, rclength(i), 'char=>char');
             infoblks(i_if) = convertCharsToStrings(infotxt);
             i_if = i_if + 1;
         end
     
        if ~mod(i,100)
            waitbar(i/length(rclength), wb, sprintf('Processing record %d / %d', i, length(rclength)));
        end

    end
    
    close(wb);
     
    % Close the file
    fclose(fid);
    
    % Reorder MPU data
    dim = size(gyrdat);
    gyrdat = reshape(gyrdat,[dim(1) dim(2) dim(3)*dim(4)]);
    dim = size(accdat);
    accdat = reshape(accdat,[dim(1) dim(2) dim(3)*dim(4)]);
    t_frmmpu = 10e-3; %Always at 100Hz = 1/10ms
end
