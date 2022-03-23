function [infoblks, ... % Free-form information fields
          enum, ...     % JSON format enumeration from the back end
          tchdat, ...   % The time increment of a single data frame (1/fps)
          chdat, ...    % Channel data (channels x frames)
          satflag, ...  % Saturation flag (channels x frames)
          tmpdat, ...   % Tile internal temperatures (tiles x frames)
          vindat, ...   % Tile input voltages (tiles x frames)
          srcpwr, ...   % Source powers (nodes x wavelengths x frames)
          evtim, ...    % Time of each event (ms)
          evstr, ...    % Associated string for each marked event
          tmpudat, ...  % The time increment of a single MPU frame
          gyrdat, ...   % Gyroscope data (nodes x dim x meas/frame x frame), units of degrees per second
          accdat] ...   % Accellerometer data (nodes x dim x meas/frame x frame), units of g
          = loadlufr(fn, varargin)  % Specify the group index (defaults to zero)

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
    if(filever ~= 1 && filever ~=2)
        error('The lumo frame file is of an unknown version %i', filever);
    end
    
    % On version two, we will need to set the particular group index from 
    % which we wish to extract data
    if(filever > 1)
        if(nargin > 1)
            groupidx_sel = varargin{1};
        else
            warning('No group index selected for data, defaulting to the first');
            groupidx_sel = 0;
        end
    fprintf('Group index %d selected \n', groupidx_sel);
    else
        groupidx_sel = 0;
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
            
            % On file versions beyond 1, a group index is inserted into
            % frame records. If the group index is not that which has been
            % selected, we skip over the data
            if(filever > 1)
                group_idx = fread(fid, 1, 'int32=>doble');
                if(group_idx ~= groupidx_sel)
                    fseek(fid, rcoffset(end) + recordlen, 'bof');
                    continue;
                end
            end

            % On file versions beyond 1, an error counter is inserted between the frame
            % index and the remaining metadata.
            if(filever > 1)
                sizeparam = fread(fid, 10, 'int32=>double');
                sizeparam = sizeparam([1 3 4 5 6 7 8 9 10]);
            else
                sizeparam = fread(fid, 9, 'int32=>double');
            end
            
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
    tchdat = sizeparam_ref(2)*1e-6; %in seconds
    fps = 1/(sizeparam_ref(2)*1e-6);
    
    fprintf('Frame rate is %.2f fps\n', fps);
        
    chdat = zeros(n_schans, n_frames, 'single');        % Channel data
    dkdat = zeros(n_dchans, n_frames, 'single');        % Dark channel data
    tmpdat = zeros(n_nodes, n_frames, 'single');        % Temperature 
    vindat = zeros(n_nodes, n_frames, 'single');        % Input voltage
    srcpwr = zeros(n_nodes, n_spw, n_frames, 'single'); % Source powers
    detmax = zeros(n_nodes, n_row, n_det, n_frames, 'single'); % Detector maximum
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
             
            
            % Skip the record if it is not from the desired group
            if(filever > 1)
                group_idx = fread(fid, 1, 'int32=>doble');
                if(group_idx ~= groupidx_sel)
                    continue;
                end
                
               % Skip frame dimension data
               fread(fid, 10, 'int32=>double');
            else

             % Skip over the frame dimension data
             fread(fid, 9, 'int32=>double');
             
            end

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
             detmax_i = reshape(fread(fid, n_det * n_row * n_nodes, 'single=>single'), n_det, n_row, n_nodes);
             detmax(:,:,:,i_fr) = permute(detmax_i, [3,2,1]);
             

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
    
    % Fix enmeration acquistion data in some versions
    if(isfield(enum, 'group'))
        
        % Merge fields into correct structure. Note that versions with this
        % erroneous field naming did not have hyperscanning support, so the
        % group index can be fixed.
        for i = 1:length(enum.groups(1).channels)
            enum.groups(1).channels(i).acq_row = enum.group.channels(i).acq_row;
            enum.groups(1).channels(i).acq_offset_us = enum.group.channels(i).acq_offset_us;
        end
        
        enum = rmfield(enum, 'group');
    
    end
    
    % Build the saturation flag matrix
    det_sat_limit = 95;

    satflag = false(n_schans, n_frames);

    % Loop over the frames and get them data    
    wb = waitbar(0, 'Building saturation flag matrix');
    
%     for j = 1:n_frames
%         i_row_idx = [enum.groups(groupidx_sel+1).channels.acq_row] + 1;
%         i_det_idx = [enum.groups(groupidx_sel+1).channels.det_optode_idx] + 1;
%         i_node_idx = [enum.groups(groupidx_sel+1).channels.det_node_idx] + 1;
%         idx = sub2ind(size(detmax), i_node_idx, i_row_idx, i_det_idx, j*ones(size(i_node_idx)));
%         satflag(:,j) = detmax(idx) > det_sat_limit;
%         if ~mod(j,100)
%             waitbar(j/n_frames, wb, sprintf('Building saturation flag matrix %d / %d', j, n_frames));
%         end
%         
%     end
        
    close(wb);
    
    % Reorder MPU data
    dim = size(gyrdat);
    gyrdat = reshape(gyrdat,[dim(1) dim(2) dim(3)*dim(4)]);
    dim = size(accdat);
    accdat = reshape(accdat,[dim(1) dim(2) dim(3)*dim(4)]);
    tmpudat = 10e-3; %Always at 100Hz = 1/10ms
end
