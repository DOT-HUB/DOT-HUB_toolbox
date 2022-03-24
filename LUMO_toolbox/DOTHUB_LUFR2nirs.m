function [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUFR2nirs(lufrFileName,layoutFileName,posCSVFileName,eventsCSVFileName,distLimit)

%This script reads a .LUFR file and also loads the associated .JSON layout from a seperate file.
%
%If subject-specific registration data (from polhemus or other method) is available
%in .cvs format, this too can be parsed. If it is parsed, the function
%DOTHUB_LUMOPolhemus2SD is called to determine optode positions from the
%polhemus data. The resulting 3D optode positions overwrite the default 3D
%data from the .JSON file.
%
%This function then outputs the .nirs variables and produces a .nirs file with
%the same name as the .LUFR. Note re-running this function will overwrite
%previous .nirs derived from the parsed LUFR file. Note that the .nirs
%file contains both 2D (SD) and 3D (SD_3D) info.
%
%######################## INPUTS ##########################################
%
% lufrFileName      :   The filename of the .LUMO data directory (if not parsed, requested)
%
% layoutFileName    :   (Optional) The path of the .json layout file (only required if not
%                       present inside .LUMO directory). If layoutFile variable is
%                       parsed it is assued to take precedence over any .json within
%                       the .LUMO directory. If layoutFile is not parsed (or parsed empty)
%                       and it is not present in the .LUMO direcroty, it will be requested.
%
% posCSVFileName    :   (Optional) The path of a .csv file containing the 3D positioning data associated
%                       with this LUMO recording. This data is a four-column CSV of
%                       position label, then x, y, z coordinate (assumed in cm).
%                       The first five rows should be Nz, Iz, Ar, Al, Cz. The
%                       following rows should be SrcA, SrcB, SrcC of each of the N
%                       tiles in turn. If polhemusCSV is not parsed, the 3D
%                       information in the .json layout file is assumed.
%
% eventsCSVFileName :   (Optional) The path of a .csv file containing two columns: time of event
%                       (in seconds) and event label. If the eventsCSV is parsed,
%                       it takes precedence over the contents of the .lumo file.
%
% distLimit         :   (Optional) The maximum mm length of a S-D distance for
%                       that channel to parsed in to the .nirs. This limit
%                       is applied to the default layout file information
%                       so that the resulting data contains the same
%                       channels, irreleevant of their subject-specific
%                       separation as defined in the polhemus data. Default
%                       is 0 (all channels parsed).
%
%######################## OUTPUTS #########################################
%
% nirs                      :   the contents of the resulting .nirs file as a structure.
%
% nirsFileName              :   Full path and name of resulting .nirs file
%
% nirsFileName.nirs         :   A .nirs file saved in the same location as the
%                               .LUMO variable. Has extra variable: SD_3D, which
%                               is the same as SD, but with 3D source and detector
%                               locations.  SD_3D also contains a 'Landmarks' field
%                               that contains the five cranial landmark locations.
%
% #########################################################################
% RJC, UCL, October 2021
%
% ############################# Updates ###################################
% #########################################################################
%
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

maxDist = 0;
if exist('distLimit','var')
    if ~isempty(distLimit)
        maxDist = distLimit;
    end
end

% Determine Status of Layout Data ###############################################################
if ~exist('layoutFileName','var') %Not parsed so check exists in .lufr
    if ~exist('enum.groups.index.layout', 'var') %Not contained in .lufr, so load
        disp('Layout file not parsed or found in lufr, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        error('Parsing a layout .json is currently required')
        %disp('Using default layout within .lufr directory (WARNING: UNTESTED)');
        %layoutFileName = [];
    end
elseif isempty(layoutFileName)
    if ~exist('enum.groups.index.layout', 'var') %Not contained in .lufr, so load
        disp('Layout file not parsed or found in lufr, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
        layoutData = jsondecode(fileread(layoutFileName));
    else
        error('Parsing a layout .json is currently required')
        %disp('Using default layout within .lufr directory (WARNING: UNTESTED)');
        %layoutFileName = [];
    end
end

% Extract Default Layout Optode Positions and create chnfilter #############################
if ~isempty(layoutFileName) && maxDist>0
    layoutData = jsondecode(fileread(layoutFileName));
    nDocks = size(layoutData.docks,1);
    %Extract 3D optode positions
    %[ src_node_id src_opt_id det_node_id det_opt_id]
    chfilter = [];
    for n = 1:nDocks
        for src = 1:3
            SrcPos(n,src,1) = layoutData.docks(n).optodes(src+4).coordinates_3d.x;
            SrcPos(n,src,2) = layoutData.docks(n).optodes(src+4).coordinates_3d.y;
            SrcPos(n,src,3) = layoutData.docks(n).optodes(src+4).coordinates_3d.z;
        end
        for det = 1:4
            DetPos(n,det,1) = layoutData.docks(n).optodes(det).coordinates_3d.x;
            DetPos(n,det,2) = layoutData.docks(n).optodes(det).coordinates_3d.y;
            DetPos(n,det,3) = layoutData.docks(n).optodes(det).coordinates_3d.z;
        end
    end
    for sn = 1:nDocks
        for src = 1:3
            for dn = 1:nDocks
                for det = 1:4
                    %horrific nested loop, but my brain can't process a better
                    %solution right now
                    if vecnorm(SrcPos(sn,src,:) - DetPos(dn,det,:)) <= maxDist
                        chfilter = [chfilter; sn src dn det];
                    end
                end
            end
        end
    end
end

% LOAD LUFR DATA  ###############################################################
[infoblks, enum, tchdat, chdat, satflag, ...                    % Primary data
          tmpdat, vindat, srcpwr, ...                           % Environmental
          evtim, evstr, ...                                     % Event markers
          tmpudat, gyrdat, accdat, ...                          % Motion
          chperm] ...                                           % Channel permutation
          = loadlufr(lufrFileName,'chfilter',chfilter);         % Variable options

% Convert to .nirs format #########################################################################

% Intensity data
d = chdat'; clear chdat

% Temporal data
t = (0:tchdat:(size(d, 1)-1)*tchdat)';

% Define SD (2D) ###############################################################
nNodes = size(enum.groups.nodes, 1);
nodes = [enum.groups.nodes.id];
SD.Lambda = unique([enum.groups.channels.src_wl]);
SD.nSrcs = nNodes * 3;
SD.nDets = nNodes * 4;
SD.SpatialUnit = 'mm';
enumperm = enum.groups(1).channels(chperm);
for i = 1:size(d,2)
    SD.MeasList(i,1) = 3*(enumperm(i).src_node_idx) + (enumperm(i).src_optode_idx - 3);
    SD.MeasList(i,2) = 4*(enumperm(i).det_node_idx) + enumperm(i).det_idx + 1;
    SD.MeasList(i,3) = 1;
    SD.MeasList(i,4) = find(SD.Lambda == enumperm(i).src_wl);
end

[~, sortInd] = sort(SD.MeasList(:,4));
SD.MeasList = SD.MeasList(sortInd,:);
d = d(:,sortInd);
satflag = satflag(sortInd,:)';

% Force create MeasListAct - all ones
SD.MeasListAct = ones(size(SD.MeasList,1),1);

% Define SD (3D) ###############################################################
% Using either this or (if available) polhemus data.
if polhemusFlag %If polhemus information is parsed, calculate S-D positions from that file
    fprintf(['Using ' posCSVFileName ' to define SD3D...\n']);
    [SD_POL, SD3DFileName] = DOTHUB_LUMOpolhemus2SD3D(posCSVFileName); %This line saves the .SD3D
    
    SD_POL.MeasList = SD.MeasList;
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

%########## Output node-wise peripheral datatypes ###########:

%Acc data
periph.Acc = accdat;
periph.Gyro = gyrdat;
periph.t_mpu = (0:tmpudat:(size(accdat, 3) - 1)*tmpudat)';

% Source power structure - source power logged for each node, wavelength and time point
periph.SrcPowers = srcpwr;

% Saturation flag
periph.SatFlag = satflag;

% Temperature of each node
periph.NodeTemp = tmpdat;

% Input voltage
periph.InputVolt = vindat;

% Pad events
s = zeros(size(t));

% Pad aux
aux = zeros(size(t,1),8);

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
            cond = find(strcmp(evList,evstr_filt(i)));
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


% % Full Measurement List
% SD.MeasList = ones(SD.nSrcs * SD.nDets * length(SD.Lambda), 4);
% 
% % Creating temporary object to change indexing of sources - keep 1:3 instead of 1:6
% enum_temp = enum.groups.channels;
% change_index = find([enum_temp.src_wl] == 850);
% for i = change_index
%     enum_temp(i).src_idx = enum_temp(i).src_idx - 3;
% end

% Build MeasList array (3rd column unused)
% for i = 1:size(chdat, 1)
%     SD.MeasList(i, 1) = 3*(enum_temp(i).src_node_idx) + (enum_temp(i).src_idx + 1);
%     SD.MeasList(i, 2) = 4*(enum_temp(i).det_node_idx) + (enum_temp(i).det_idx + 1);
%     SD.MeasList(i, 4) = find((SD.Lambda == enum_temp(i).src_wl));
% end

% % Sort MeasList by wavelength, update d and satflag accordingly
% [SD.MeasList, sort_indx] = sortrows(SD.MeasList, 4);
% d = double(d(:,sort_indx));
% periph.SatFlag = periph.SatFlag(:,sort_indx);