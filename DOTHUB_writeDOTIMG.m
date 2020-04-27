function [dotimg, dotimgFileName] = DOTHUB_writeDOTIMG(dotimgFileName,logData,hbo,hbr,mua,timebase,saveFlag)

% This script creates a DOT images file (.dotimg)

% ####################### INPUTS ##########################################

% dotimgFileName       :  The desired path &/ filename for the .dotimg file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet.
                        
%                       [pathstr, name, ~] = fileparts(prepro.fileName);
%                       ds = datestr(now,'yyyymmDDHHMMSS');
%                       dotimgFileName = fullfile(pathstr,[name '.dotimg']);
%                       logData(1,:) = {'Created on: ', ds};
%                       logData(2,:) = {'Associated prepro file: ', prepro.fileName};
%                       logData(3,:) = {'Associated invjac file: ', invjac.fileName};
%                       logData(4,:) = {'reconMethod: ', varInputs.reconMethod};
%                       logData(5,:) = {'regMethod: ', varInputs.regMethod};
%                       logData(6,:) = {'hyperParameter: ', varInputs.regMethod};

% logData           :  logData is a cell array of strings containing useful
%                      info as per snippet above. Parse empty to ignore.
 
% hbo               :  Structure containing hbo images (.vol, .gm), of
%                      dimensions nFrames x nNodes. Note that both can be
%                      empty if not requested in DOTHUB_reconstruction

% hbr               :  Structure containing hbr images (.vol, .gm), of
%                      dimensions nFrames x nNodes. Note that both can be
%                      empty if not requested in DOTHUB_reconstruction

% mua               :  Cell of length nWavs with structures (.vol, .gm) where
%                      images of dimensions nFrames x nNodes are stored. 
%                      Note that one or both can beempty if not requested 
%                      in DOTHUB_reconstruction

% timebase          :  A time vector with the same length as the first
%                      dimension of each image set.

% saveFlag          :  True if .dotimg file is to be saved (not always the
%                      case e.g. for online reconstruction

% ####################### OUTPUTS #########################################

% dotimg               :  Structure containing all relevant data inputs

% dotimgFileName       :  The full path of the resulting .dotimg file

% dotimgFileName.dotimg   :  File containing contents of dotimg structure

% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

dotimg.logData = logData;
dotimg.hbo = hbo;
dotimg.hbr = hbr;
dotimg.mua = mua;
dotimg.timebase = timebase;

%Save .dotimg file ###########################################################
%Create filename ##########################################################
[pathstr, name, ext] = fileparts(dotimgFileName);
if isempty(ext) || ~strcmpi(ext,'.dotimg')
    ext = '.dotimg';
end
if isempty(pathstr)
    pathstr = pwd;
end
dotimgFileName = fullfile(pathstr,[name ext]);
dotimg.fileName = dotimgFileName; %including the fileName within the structure 
%is very useful for tracking and naming things derived further downstream.

%Save .dotimg file ###########################################################
if saveFlag
    if exist(dotimgFileName,'file')
        warning([name ext ' will be overwritten...']);
    end
    save(dotimgFileName,'-struct','dotimg');
    fprintf('################### Writing .dotimgimg file #######################\n');
    fprintf(['.dotimg data file saved as ' dotimgFileName '\n']);
end

