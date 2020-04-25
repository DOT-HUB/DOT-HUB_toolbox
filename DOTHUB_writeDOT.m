function [dot, dotFileName] = DOTHUB_writeDOT(dotFileName,logData,hbo,hbr,mua,timebase,saveFlag)

% This script creates a DOT images file (.dot)

% ####################### INPUTS ##########################################

% dotFileName       :  The desired path &/ filename for the .dot file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet.
                        
%                       [pathstr, name, ~] = fileparts(prepro.fileName);
%                       ds = datestr(now,'yyyymmDDHHMMSS');
%                       dotFileName = fullfile(pathstr,[name '_' ds '.dot']);
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

% saveFlag          :  True if .dot file is to be saved (not always the
%                      case e.g. for online reconstruction

% ####################### OUTPUTS #########################################

% dot               :  Structure containing all relevant data inputs

% dotFileName       :  The full path of the resulting .dot file

% dotFileName.dot   :  File containing contents of dot structure

% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

dot.logData = logData;
dot.hbo = hbo;
dot.hbr = hbr;
dot.mua = mua;
dot.timebase = timebase;

%Save .dot file ###########################################################
%Create filename ##########################################################
[pathstr, name, ext] = fileparts(dotFileName);
if isempty(ext) || ~strcmpi(ext,'.dot')
    ext = '.dot';
end
if isempty(pathstr)
    pathstr = pwd;
end
dotFileName = fullfile(pathstr,[name ext]);
dot.fileName = dotFileName; %including the fileName within the structure 
%is very useful for tracking and naming things derived further downstream.

%Save .dot file ###########################################################
if saveFlag
    save(dotFileName,'-struct','dot');
    fprintf(['.dot data file saved as ' dotFileName '\n']);
end

