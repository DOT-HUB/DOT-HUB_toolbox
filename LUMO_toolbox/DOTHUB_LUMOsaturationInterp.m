function [dclean, SD] = DOTHUB_LUMOsaturationInterp(d,SD,satflag,maxSatTime)

%This function takes the raw intensity data (d) and the channel and time
%wise saturation flag availiable within LUMO lufr data to interpolated
%across the periods of saturation for each channel. If a period of
%saturation exceeds maxTime, that channel is deactivated in SD.MeasListAct
%
% ####################### INPUTS ##########################################
%
% d                     :   default homer style intensity matrix
%
% SD                    :   .nirs SD variable
%
% satflag               :   saturation flag matrix where == 1 means
%                           saturation. Same dimensions as d.
%
% maxSatTime            :   max duration of any single saturation period before
%                           channel is discarded. Default 20s?
%
% maxSatBurden          :   max proportion of recording a channel can be
%                           saturated before being discarded (0:1)

% ####################### OUTPUTS #########################################

% dclean                :   cleaned homer style intensity matrix

% SD                    :   updated .nirs SD variable

% ####################### Dependencies ####################################

% #########################################################################
% Written by RJC, April 2022
%
% ############################# Updates ###################################
% #########################################################################

%Input handling

nChan = size(d,2);
nDP = size(d,1);
for i = 1:nChan
    if any(satflag(:,i)) %There is saturation
        if sum(satflag(:,i))/nDP > maxSatBurden
            SD.MeasListAct(i) = 0;
        else

end
