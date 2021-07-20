function [dotimg, dotimgFileName] = DOTHUB_writeDOTIMG(dotimgFileName,logData,hbo,hbr,mua,tImg,saveFlag)

% This script creates a DOT images file (.dotimg)
%
% ####################### INPUTS ##########################################
%
% dotimgFileName       :  The desired path &/ filename for the .dotimg file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet.
%
%                       [pathstr, name, ~] = fileparts(prepro.fileName);
%                       ds = datestr(now,'yyyymmDDHHMMSS');
%                       dotimgFileName = fullfile(pathstr,[name '.dotimg']);
%                       logData(1,:) = {'Created on: ', ds};
%                       logData(2,:) = {'Associated prepro file: ', prepro.fileName};
%                       logData(3,:) = {'Associated invjac file: ', invjac.fileName};
%                       logData(4,:) = {'reconMethod: ', varInputs.reconMethod};
%                       logData(5,:) = {'regMethod: ', varInputs.regMethod};
%                       logData(6,:) = {'hyperParameter: ', varInputs.regMethod};
%
% logData           :  logData is a cell array of strings containing useful
%                      info as per snippet above. Parse empty to ignore.
%
% hbo               :  Structure containing hbo images (.vol, .gm), of
%                      dimensions nFrames x nNodes. Note that both can be
%                      empty if not requested in DOTHUB_reconstruction
%
% hbr               :  Structure containing hbr images (.vol, .gm), of
%                      dimensions nFrames x nNodes. Note that both can be
%                      empty if not requested in DOTHUB_reconstruction
%
% mua               :  Cell of length nWavs with structures (.vol, .gm) where
%                      images of dimensions nFrames x nNodes are stored.
%                      Note that one or both can beempty if not requested
%                      in DOTHUB_reconstruction
%
% tIMG              :  A time vector with the same length as the first
%                      dimension of each image set (prepro.tDOD)
%
% saveFlag          :  True if .dotimg file is to be saved (not always the
%                      case e.g. for online reconstruction
%
% ####################### OUTPUTS #########################################
%
% dotimg               :  Structure containing all relevant data inputs
%
% dotimgFileName       :  The full path of the resulting .dotimg file
%
% .dotimg              :  (If saveFlag) A file containing a structure of:
%                       dotimg.logData      - as defined above
%                       dotimg.hbo          - as defined above
%                       dotimg.hbr          - as defined above
%                       dotimg.mua          - as defined above
%                       dotimg.tImg         - as defined above
%                       dotimg.fileName  	- the path of the saved dotimg file
%
% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################
%
% MANAGE VARIABLES
% #########################################################################

if isempty(logData)
    logData = {};
    warning('logData is empty: this might make it harder to keep track of your data...');
end

dotimg.logData = logData;
dotimg.hbo = hbo;
dotimg.hbr = hbr;
dotimg.mua = mua;
dotimg.tImg = tImg;

%Save .dotimg file ###########################################################
%Create filename ##########################################################
%Save .dotimg file ###########################################################
if saveFlag
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
    if exist(dotimgFileName,'file')
        warning([name ext ' will be overwritten...']);
    end
    save(dotimgFileName,'-struct','dotimg','-v7.3');
    fprintf('################### Writing .dotimgimg file #######################\n');
    fprintf(['.dotimg data file saved as ' dotimgFileName '\n']);
end

