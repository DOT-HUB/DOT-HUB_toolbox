function DOTHUB_LUMOplotMPU(nirsFileName,varargin)

% This function creates a matrix plot of log intensity
%
%######################## INPUTS ##########################################
%
% nirsFileName      = nirs file pathname or structure
% varargin          =  optional input pairs:
%                       'hAxes' - optional axis handle

%######################## OUTPUTS #########################################
%
%Outputs are figures...
%
%######################## Dependencies ####################################
%This script requires other functions in the DOTHUB function library
%
% #########################################################################
% RJC, UCL, Jan 2022
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES  ##################################################
if ~exist('nirsFileName','var')
    disp('Select .nirs file...');
    [file,path] = uigetfile('*.nirs','Select .nirs file');
    nirsFileName = fullfile(path,file);
end

varInputs = inputParser;
addParameter(varInputs,'hAxes','',@ishandle);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;
if isempty(varInputs.hAxes)
    varInputs.hAxes = gca;
end

%Load or rename data
if ischar(nirsFileName)
    nirs = load(nirsFileName,'-mat');
else
    nirs = nirsFileName;
end

%Plot by tile
nTiles = size(nirs.periph.Acc,1);
cols = hsv(nTiles);
figure;
p1 = subplot(3,1,1);
for i = 1:nTiles
    plot(nirs.periph.t_mpu,squeeze(nirs.periph.Acc(i,:,:))','color',cols(i,:));
    hold on
end
xlabel('Time (s)')
ylabel('g');
title('Accelerometers')
p2 = subplot(3,1,2);
for i = 1:nTiles
    plot(nirs.periph.t_mpu,squeeze(nirs.periph.Gyro(i,:,:))','color',cols(i,:));
    hold on
end
xlabel('Time (s)')
ylabel('deg.s^-^1');
title('Gyroscopes')
p3 = subplot(3,1,3);
plot(nirs.t,nirs.d);
xlabel('Time (s)')
ylabel('Intensity (arb.)');
title('Intensity data')

linkaxes([p1 p2 p3],'x')





