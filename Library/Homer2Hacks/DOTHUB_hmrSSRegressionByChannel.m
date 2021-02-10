%Homer2 Software License Agreement 
%Copyright © 2015, David Boas, Jay Dubb, Ted Huppert
%All rights reserved. 

function y_reg = DOTHUB_hmrSSRegressionByChannel(y, SD, rhoSD_ssThresh, flagSSmethod)

%This function regresses a SS channel from each long channel for the entire
%data period, rather than ony within the HRF GLM process. This should
%allow data that does not have enough trials for a GLM approach to take
%advantage of the SS concept.  ONLY WORKS FOR 2 wavelength DOD data.
%
%INPUTS
% y - data (dc or 2-wav DOD)
% SD - source detector stucture (units should be consistent with rhoSD_ssThresh)
% rhoSD_ssThresh - max distance for a short separation measurement. Set =0
%          if you do not want to regress the short separation measurements.
% flagSSmethod - 0 if short separation regression is performed with the nearest
%               short separation channel.
%            1 if performed with the short separation channel with the
%               greatest correlation.
%            2 if performed with average of all short separation channels.
%            3 if performed with the average of the short separation
%              channels closer ( < rhoSD_ssThresh) to the long
%              source-detector separation channel
%            4 if performed with the average of the short separation
%              channels closer ( < rhoSD_ssThresh) to the location of the
%              source and detector of a long Src-Det separation channel
%
%OUTPUTS
%y_reg - the channel-wise-regressed version of y;
%
%####################################################################
% Robert J Cooper, University College London, April 2016
% Amended by EEVR, UCL, 2020
%####################################################################
%
%####################################################################
%Change log
% 15/06/2017 Updated to accept DOD data, and simplify.
% 19/04/2020 Correction: the case for flagSSmethod==1  (use SS with highest 
%            correlation) was regressing the SS channels only for the 1st 
%            wavelength and not the 2nd.
% 20/04/2020 (EEVR, UCL) Updated to regress two new types of regressors:
%             A) Regress the mean of the SS channels closer to the long SD
%                separation channel. Select: flagSSmethod = 3
%             B) Regress the mean of the SS channels that are closer to the
%                source and detector sides of a long SD channels: Select:
%                flagSSmethod = 4
%####################################################################

%Debug variables
% rhoSD_ssThresh = 26;
% flagSSmethod = 1;
% load('Finger_Tapping_1.nirs','-mat');
% tmpOD = hmrIntensity2OD(d);
% procResult.dc = hmrOD2Conc(tmpOD,SD,[6 6]);
% y = procResult.dc;
%########################

ml = SD.MeasList;
mlAct = [SD.MeasListAct(1:end/2); SD.MeasListAct(1:end/2)];
lstPerWav = find(ml(:,4)==1);
nChan = length(ml);
nChanPerWav = length(lstPerWav);
lst = 1:nChan;

rhoSDPerWav = zeros(nChanPerWav,1);
posMPerWav = zeros(nChanPerWav,3);

for iML = 1:nChanPerWav
    rhoSDPerWav(iML) = sum((SD.SrcPos(ml(lstPerWav(iML),1),:) - SD.DetPos(ml(lstPerWav(iML),2),:)).^2).^0.5;
    posMPerWav(iML,:) = (SD.SrcPos(ml(lstPerWav(iML),1),:) + SD.DetPos(ml(lstPerWav(iML),2),:)) / 2;
end

rhoSD = [rhoSDPerWav; rhoSDPerWav];
posM = [posMPerWav; posMPerWav];

lstSS = lst(find(rhoSD<=rhoSD_ssThresh & mlAct(lst)==1)); %#ok<*FNDSB>
lstSSPerWav = lstPerWav(find(rhoSDPerWav<=rhoSD_ssThresh & mlAct(lstPerWav)==1));

if length(size(y)) == 2 %DOD DATA
    %Populate regressors matrix
    ssRegressors = nan(size(y,1),nChan);
    
    %Check we have some SS channels
    if isempty(lstSS)
        fprintf('There are no short separation channels in this probe\n');
        y_reg = y;
        return
    end
    
    if flagSSmethod==0  % use nearest SS
        for iML = 1:nChan
            rho = sum((ones(length(lstSS),1)*posM(iML,:) - posM(lstSS,:)).^2,2).^0.5;
            [~,ii] = min(rho);
            
            %populate regressors (same for HbO and HbR)
            ssRegressors(:,iML) = y(:,lstSS(ii));
        end
        
    elseif flagSSmethod==1 % use SS with highest correlation
                
        tmpy = y;
        tmpy = (tmpy-ones(length(tmpy),1)*mean(tmpy,1))./(ones(length(tmpy),1)*std(tmpy,[],1));
        cc(:,:) = tmpy'*tmpy / length(tmpy);
        
        clear tmpy
        
        % find short separation channel with highest correlation (separate by wavelength)
        for iML = 1:nChanPerWav
            [~,ii] = max(cc(iML,lstSS(1:end/2)));
            ssRegressors(:,iML) = y(:,lstSS(ii));
        end
        for iML = nChanPerWav+1:nChan
            [~,ii] = max(cc(iML,lstSS(end/2+1:end)));
            ssRegressors(:,iML) = y(:,lstSS(ii+numel(lstSS)/2));% Corrected! now pointing to 2nd half of the array
        end
        
    elseif flagSSmethod==2 % use average of all active SS as regressor (separate by wavelength)
        
        mnSS1 = squeeze(mean(y(:,lstSS(1:end/2)),2));
        mnSS2 = squeeze(mean(y(:,lstSS(end/2+1:end)),2));
        tmp1 = repmat(mnSS1,1,nChanPerWav);
        tmp2 = repmat(mnSS2,1,nChanPerWav); 
        ssRegressors = [tmp1 tmp2];
       
    elseif flagSSmethod==3 % use local average of many SS around channel as regressor (separate by wavelength)

        for iML = 1:nChanPerWav
            rho = sum((ones(length(lstSS(1:end/2)),1)*posM(iML,:) - posM(lstSS(1:end/2),:)).^2,2).^0.5;        
            ii = find(rho <= rhoSD_ssThresh);
            if ~isempty(ii)
                 ssRegressors(:,iML) = mean(y(:,lstSS(ii)),2);
            end
        end            

        for iML = nChanPerWav+1:nChan
            rho = sum((ones(length(lstSS(end/2+1:end)),1)*posM(iML,:) - posM(lstSS(end/2+1:end),:)).^2,2).^0.5;        
            ii = find(rho <= rhoSD_ssThresh);
            if ~isempty(ii)
                 ssRegressors(:,iML) = mean(y(:,lstSS(ii+numel(lstSS)/2)),2);
            end
        end        

    elseif flagSSmethod==4 % use local average of many SS around source and detector as regressor (separate by wavelength)

        for iML = 1:nChanPerWav
            
            % SOURCE SIDE------------------------------------------------            
            Src = SD.SrcPos(ml(lst(iML),1),:);
            rho = sum((ones(length(lstSS(1:end/2)),1)*Src - posM(lstSS(1:end/2),:)).^2,2).^0.5;        
            ii = find(rho <= rhoSD_ssThresh);
            
            % DETECTOR SIDE------------------------------------------------
            Det = SD.DetPos(ml(lst(iML),2),:);
            rho = sum((ones(length(lstSS(1:end/2)),1)*Det - posM(lstSS(1:end/2),:)).^2,2).^0.5;        
            jj = find(rho <= rhoSD_ssThresh);
            
            if ~isempty(ii) && isempty(jj)
                ssRegressors(:,iML) = mean(y(:,lstSS(ii)),2);
            elseif isempty(ii)&& ~isempty(jj)
                ssRegressors(:,iML) = mean(y(:,lstSS(jj)),2);
            elseif ~isempty(ii)&& ~isempty(jj)
                kk = unique([ii; jj]);% Delete repeated channels
                ssRegressors(:,iML) = mean(y(:,lstSS(kk(:))),2);
            end
        end
        
        for iML = nChanPerWav+1:nChan
            
            % SOURCE SIDE------------------------------------------------            
            Src = SD.SrcPos(ml(lst(iML),1),:);
            rho = sum((ones(length(lstSS(end/2+1:end)),1)*Src - posM(lstSS(end/2+1:end),:)).^2,2).^0.5;        
            ii = find(rho <= rhoSD_ssThresh);
            
            % DETECTOR SIDE------------------------------------------------
            Det = SD.DetPos(ml(lst(iML),2),:);
            rho = sum((ones(length(lstSS(end/2+1:end)),1)*Det - posM(lstSS(end/2+1:end),:)).^2,2).^0.5;        
            jj = find(rho <= rhoSD_ssThresh);
            
            if ~isempty(ii) && isempty(jj)
                ssRegressors(:,iML) = mean(y(:,lstSS(ii+numel(lstSS)/2)),2);
            elseif isempty(ii)&& ~isempty(jj)
                ssRegressors(:,iML) = mean(y(:,lstSS(jj+numel(lstSS)/2)),2);                
            elseif ~isempty(ii)&& ~isempty(jj)
                kk = unique([ii; jj]);% Delete repeated channels
                ssRegressors(:,iML) = mean(y(:,lstSS(kk(:)+numel(lstSS)/2)),2);
            end            
            
        end
        
    end
    
    %Now loop around channels and perform regression;
    y_reg = nan(size(y));
    
    for i = 1:nChan %Loop over channels
        
            lChan = y(:,i);
            ssChan = ssRegressors(:,i);
            
            %We are saying that lChan = B*ssChan + residual
            %B ~ inv(ssChan)*lChan;
            %So residual = lChan - B*ssChan
            if ~any(isnan(ssChan))
                B = pinv(ssChan)*lChan;
                tmpReg = lChan - B*ssChan;
                y_reg(:,i) = tmpReg;   
            else
                y_reg(:,i) = lChan;
            end
    end
    
elseif length(size(y)) == 3 %DC DATA
    %Populate regressors matrix
    ssRegressors = nan(size(y,1),2,nChanPerWav);
    
    %Check we have some SS channels
    if isempty(lstSSPerWav)
        fprintf('There are no short separation channels in this probe\n');
        y_reg = y;
        return
    end
    
    if flagSSmethod==0  % use nearest SS
        for iML = 1:nChanPerWav
            rho = sum((ones(length(lstSSPerWav),1)*posMPerWav(iML,:) - posMPerWav(lstSSPerWav,:)).^2,2).^0.5;
            [~,ii] = min(rho);
            
            %populate regressors (same for HbO and HbR)
            ssRegressors(:,:,iML) = y(:,[1 2],lstSSPerWav(ii));
        end
        
    elseif flagSSmethod==1 % use SS with highest correlation
                
        % HbO
        dc = squeeze(y(:,1,:));
        dc = (dc-ones(length(dc),1)*mean(dc,1))./(ones(length(dc),1)*std(dc,[],1));
        cc(:,:,1) = dc'*dc / length(dc);
        
        % HbR
        dc = squeeze(y(:,2,:));
        dc = (dc-ones(length(dc),1)*mean(dc,1))./(ones(length(dc),1)*std(dc,[],1));
        cc(:,:,2) = dc'*dc / length(dc);
        
        clear dc
        
        % find short separation channel with highest correlation
        for iML = 1:size(cc,1)
            % HbO
            [~,ii] = max(cc(iML,lstSSPerWav,1));
            ssRegressors(:,1,iML) = y(:,1,lstSSPerWav(ii));
            % HbR
            [~,ii] = max(cc(iML,lstSSPerWav,2));
            ssRegressors(:,2,iML) = y(:,2,lstSSPerWav(ii));
        end
        
    elseif flagSSmethod==2 % use average of all active SS as regressor
        
        mnSS = squeeze(mean(y(:,[1 2],lstSSPerWav),3));
        tmp = repmat(mnSS,1,size(y,3));
        ssRegressors = reshape(tmp,[size(y,1),2,size(y,3)]);

    elseif flagSSmethod==3 % use local average of many SS around channel as regressor (separate by wavelength)

        for iML = 1:nChanPerWav
            rho = sum((ones(length(lstSSPerWav),1)*posM(iML,:) - posM(lstSSPerWav,:)).^2,2).^0.5;        
            ii = find(rho <= rhoSD_ssThresh);
            if ~isempty(ii)
                % HbO
                 ssRegressors(:,1,iML) = mean(squeeze(y(:,1,lstSSPerWav(ii))),2);
                % HbR
                 ssRegressors(:,2,iML) = mean(squeeze(y(:,2,lstSSPerWav(ii))),2);                
            end
        end

    elseif flagSSmethod==4 % use local average of many SS around source and detector as regressor (separate by wavelength)
        
        dists = DOTHUB_getSDdists(SD);
        dists = [dists dists];
        for iML = 1:nChanPerWav

            % SOURCE SIDE------------------------------------------------            
            Src = SD.SrcPos(ml(lst(iML),1),:);
            rho = sum((ones(length(lstSSPerWav),1)*Src - posM(lstSSPerWav,:)).^2,2).^0.5;        
            ii = find(rho <= rhoSD_ssThresh);
            
            % DETECTOR SIDE------------------------------------------------
            Det = SD.DetPos(ml(lst(iML),2),:);
            rho = sum((ones(length(lstSSPerWav),1)*Det - posM(lstSSPerWav,:)).^2,2).^0.5;        
            jj = find(rho <= rhoSD_ssThresh);
            
            if ~isempty(ii) && isempty(jj)
              % HbO
                ssRegressors(:,1,iML) = mean(squeeze(y(:,1,lstSS(ii))),2);
              % HbR
                ssRegressors(:,2,iML) = mean(squeeze(y(:,2,lstSS(ii))),2);
            elseif isempty(ii) && ~isempty(jj)
              % HbO
                ssRegressors(:,1,iML) = mean(squeeze(y(:,1,lstSS(jj))),2);
              % HbR
                ssRegressors(:,2,iML) = mean(squeeze(y(:,2,lstSS(jj))),2);
            elseif ~isempty(ii) && ~isempty(jj)
              kk = unique([ii(:); jj(:)]); % Delete repeated channels
              % HbR
                ssRegressors(:,1,iML) = mean(squeeze(y(:,1,lstSS(kk(:)))),2);
              % HbO
                ssRegressors(:,2,iML) = mean(squeeze(y(:,2,lstSS(kk(:)))),2);
            end

        end    
        
    end
    
    %Now loop around channels and perform regression;
    y_reg = zeros(size(y));
    
    for i = 1:size(y,3) %Loop over channels
        for j = 1:2 %Loop over haem
            
            lChan = squeeze(y(:,j,i));
            ssChan = squeeze(ssRegressors(:,j,i));
            
            %We are saying that lChan = B*ssChan + residual
            %B ~ inv(ssChan)*lChan;
            %So residual = lChan - B*ssChan
            if ~any(isnan(ssChan))
                B = pinv(ssChan)*lChan;
                tmpReg = lChan - B*ssChan;
                y_reg(:,j,i) = tmpReg;
            else
                disp(['No local short channels found for channel ' num2str(i)]);
                y_reg(:,j,i) = lChan;
            end
        end
        tmpHbT = squeeze(y_reg(:,1,i)) + squeeze(y_reg(:,2,i));
        y_reg(:,3,i) = tmpHbT;
    end
end

%Debug output
%###########
%procResult.dc = y_reg;
%save('Finger_Tapping_reg.nirs','s','t','d','aux','SD','procInput','procResult','userdata','tIncMan');



