function SD = DOTHUB_balanceMeasListAct(SD)

%This function ensures every wavelength of a given channel is marked
%inactive if any wavelength of that channel is marked inactive

%RJC UCL, April 2020 ############################################################
nWavs = length(SD.Lambda);
tmp = reshape(SD.MeasListAct,length(SD.MeasListAct)/nWavs,nWavs);
tmp2 = all(tmp')';
SD.MeasListAct = repmat(tmp2,nWavs,1);