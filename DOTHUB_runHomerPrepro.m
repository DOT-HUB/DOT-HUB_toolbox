function [prepro, preproFileName] = DOTHUB_runHomerPrepro(nirsFileName,cfgFileName,variableToReconstruct)

% Runs Homer processing stream defined in cfgFileName on nirs file
% nirsFileName, and outputs relevant results into a prepro file. The output
% variable to be reconstructed is specified by 'variableToReconstruct'
% this script was hacked out of Homer2 functions, written by Jay Dubb at
% MGH.  See Homer2 https://www.nitrc.org/projects/homer2/
% ####################### INPUTS ##########################################
%
% nirsFileName              : Path to .nirs file.
%
% cfgFileName               : Path to .cfg file.
%
% variableToReconstruct     : 'full' or 'HRF' (default 'HRF').  If 'full',
%                              the dod output from the .cfg file processing
%                              stream (i.e. the full time course after
%                              filter/motion correction or whatever) is 
%                              saved into prepro.dod (which is the variable 
%                              actually reconstructed downstream in the
%                              DOTHUB pipeline).  If 'HRF' is selected, the
%                              dcAvg variable is converted back to dod and
%                              this is saved into prepro.dod. dcAvg and
%                              dcStd are also saved 
%
% ####################### OUTPUTS #########################################
%
% prepro            : The prepro stucture
%
% preproFileName    : The full pathname of the resulting prepro file
%
% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020

fprintf('################# Running DOTHUB_runHomerPrepro ##################\n');

% MANAGE VARIABLES
% #########################################################################
[pathstr, name, ext] = fileparts(nirsFileName);
if isempty(ext) || ~strcmpi(ext,'.nirs')
    ext = '.nirs';
end
if isempty(pathstr) %Enforce full pathname
    pathstr = pwd;
end
nirsFileName = fullfile(pathstr,[name ext]);

[pathstr, name, ext] = fileparts(cfgFileName);
if isempty(ext) || ~strcmpi(ext,'.cfg')
    ext = '.cfg';
end
if isempty(pathstr) %Enforce full pathname
    pathstr = pwd;
end
cfgFileName = fullfile(pathstr,[name ext]);

if exist('variableToReconstruct','var')
    ind = find(strcmpi(variableToReconstruct,{'full','hrf'}));
    if isempty(ind)
        error('The variableToReconstruct input must be ''full'' or ''hrf''.');
    end
else
    variableToReconstruct = 'HRF';
end
if strcmpi(variableToReconstruct,'full')
    fullFlag = 1;
else
    fullFlag = 0;
end    
   
% LOAD DATA
% #########################################################################
% Load .nirs
hmr = load(nirsFileName,'-mat');

% Overwrite hmr.SD with hmr.SD3D to match downstream Homer2 function calls
SD2D = hmr.SD;
hmr.SD = hmr.SD3D;

% Load .cfg
fid = fopen(cfgFileName,'r');
[procInput, ~] = parseProcessOpt(fid);
fclose(fid);

% RUN STREAM
% #########################################################################
paramOut = {};
fcallList = {};
hwait = waitbar(0, 'Processing...' );
flagNoRun = 0;
for iFunc = 1:procInput.procFunc.nFunc    
    wb = waitbar( iFunc/procInput.procFunc.nFunc, hwait, sprintf('Processing... %s',procInput.procFunc.funcName{iFunc}) );
    % Extract input arguments from hmr
    argIn = parseProcessFuncArgsIn(procInput.procFunc.funcArgIn{iFunc});
    for ii = 1:length(argIn)
        if ~exist(argIn{ii},'var')
            if isfield(hmr,argIn{ii})
                eval(sprintf('%s = hmr.%s;',argIn{ii},argIn{ii}));
            else
                if strcmpi(argIn{ii},'tIncMan')
                    tIncMan = ones(length(hmr.t),1); %RJC addition
                else
                    eval(sprintf('%s = [];',argIn{ii}));  % if variable doesn't exist and not in hmr then make it empty DAB 11/8/11
                end
            end
        end
    end

    % parse input parameters
    p = [];
    sargin = '';
    sarginVal = '';
    for iP = 1:procInput.procFunc.nFuncParam(iFunc)
        if ~procInput.procFunc.nFuncParamVar(iFunc)
            p{iP} = procInput.procFunc.funcParamVal{iFunc}{iP};
        else
            p{iP}.name = procInput.procFunc.funcParam{iFunc}{iP};
            p{iP}.val = procInput.procFunc.funcParamVal{iFunc}{iP};
        end
        if length(procInput.procFunc.funcArgIn{iFunc})==1 & iP==1
            sargin = sprintf('%sp{%d}',sargin,iP);
            if isnumeric(p{iP})
                if length(p{iP})==1
                    sarginVal = sprintf('%s%s',sarginVal,num2str(p{iP}));
                else
                    sarginVal = sprintf('%s[%s]',sarginVal,num2str(p{iP}));
                end
            elseif ~isstruct(p{iP})
                sarginVal = sprintf('%s,%s',sarginVal,p{iP});
            else
                sarginVal = sprintf('%s,[XXX]',sarginVal);
            end
        else
            sargin = sprintf('%s,p{%d}',sargin,iP);
            if isnumeric(p{iP})
                if length(p{iP})==1
                    sarginVal = sprintf('%s,%s',sarginVal,num2str(p{iP}));
                else
                    sarginVal = sprintf('%s,[%s]',sarginVal,num2str(p{iP}));
                end
            elseif ~isstruct(p{iP})
                sarginVal = sprintf('%s,%s',sarginVal,p{iP});
            else
                sarginVal = sprintf('%s,[XXX]',sarginVal);
            end
        end
    end
    
    % set up output format
    sargout = procInput.procFunc.funcArgOut{iFunc};
    for ii=1:length(procInput.procFunc.funcArgOut{iFunc})
        if sargout(ii)=='#'
            sargout(ii) = ' ';
        end
    end
    
    % call function
    fcall = sprintf( '%s = %s%s%s);', sargout, ...
        procInput.procFunc.funcName{iFunc}, ...
        procInput.procFunc.funcArgIn{iFunc}, sargin );
    if flagNoRun==0
        try 
            eval( fcall );
        catch ME
	        msg = sprintf('Function %s generated ERROR at line %d: %s', procInput.procFunc.funcName{iFunc}, ME.stack(1).line, ME.message);
            menu(msg,'OK');
            close(hwait);
            assert(logical(0), msg);
        end
    end
    fcallList{end+1} = sprintf( '%s = %s%s%s);', sargout, ...
        procInput.procFunc.funcName{iFunc}, ...
        procInput.procFunc.funcArgIn{iFunc}, sarginVal );
    
    % parse output parameters
    foos = procInput.procFunc.funcArgOut{iFunc};
    % remove '[', ']', and ','
    for ii=1:length(foos)
        if foos(ii)=='[' | foos(ii)==']' | foos(ii)==',' | foos(ii)=='#'
            foos(ii) = ' ';
        end
    end
    % get parameters for Output to hmr.procResult
    lst = strfind(foos,' ');
    lst = [0 lst length(foos)+1];
    param = [];
    for ii=1:length(lst)-1
        foo2 = foos(lst(ii)+1:lst(ii+1)-1);
        lst2 = strmatch( foo2, paramOut, 'exact' );
        idx = strfind(foo2,'foo');
        if isempty(lst2) & (isempty(idx) || idx>1) & ~isempty(foo2)
            paramOut{end+1} = foo2;
        end
    end
end
delete(wb);

% Output variables for reconstruction ###################################
%Force MeasListAct to be the same across wavelengths
SD3D = DOTHUB_balanceMeasListAct(hmr.SD3D);

if fullFlag %Full timecourse is to be reconstructed
    dodRecon = dod;
    tRecon = t;
    if ~exist(dcAvg,'var') %User may have calculated HRFs but want to recon 
                           %full timecourse. Still save HRF is exists, else 
                           %parse empty (so I only need one writePREPRO
                           %line
        dcAvg = [];
    end
    if ~exist(dcAvgStd,'var')
        dcAvgStd = [];
    end
    if ~exist(tHRF,'var')
        tHRF = [];
    end
else %HRF is to be reconstructed
    %Find the DPF values used in cfg
    funcInd = find(strcmpi(procInput.procFunc.funcName,'hmrOD2Conc'));
    DPFs = procInput.procFunc.funcParamVal{funcInd}{:};
    dodRecon = DOTHUB_hmrConc2OD(dcAvg,SD3D,DPFs);
    tRecon = tHRF;
    
    if isfield(hmr,'CondNames')
        condNames = hmr.CondNames;
    else
        condNames = {};
    end
end

% USE CODE SNIPPET FROM DOTHUB_writePREPRO to define filename and logData
[pathstr, name, ~] = fileparts(nirsFileName);
ds = datestr(now,'yyyymmDDHHMMSS');
preproFileName = fullfile(pathstr,[name '.prepro']);
logData(1,:) = {'Created on: '; ds};
logData(2,:) = {'Derived from data: ', nirsFileName};
logData(3,:) = {'Pre-processed using:', cfgFileName};

%(preproFileName,logData,dod,tDOD,SD3D,s,dcAvg,dcAvgStd,tHRF)
%Convert to uM
dcAvg = 1e6*dcAvg;
dcAvgStd = 1e6*dcAvgStd;
[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,dodRecon,tRecon,SD3D,hmr.s,dcAvg,dcAvgStd,tHRF,condNames,SD2D);

