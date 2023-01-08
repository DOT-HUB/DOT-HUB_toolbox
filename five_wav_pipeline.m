% DOT-HUB toolbox script to run five-wavelength data

% based off RJC 2020 UCL
% modified by GCVL 2022 Cambridge

%% Specify paths of pre-defined elements (.LUMO, atlas .mshs, Homer2 preprocessing .cfg file).

[filepath,~,~] = fileparts(mfilename('fullpath'));
if ~exist('LUMODirName','var')
    disp('Select .LUMO directory...');
    LUMODirName = uigetdir(pwd,'Select .LUMO directory');
elseif isempty(LUMODirName)
    disp('Select .LUMO directory...');
    LUMODirName = uigetdir(pwd,'Select .LUMO directory');
elseif ~exist(LUMODirName,'dir')
    disp('Specified directory not found, please select .LUMO directory...');
    LUMODirName = uigetdir(pwd,'Select .LUMO directory');
end
origMeshFileName = [filepath '/ExampleMeshes/AdultMNI152.mshs'];
%cfgFileName = [filepath '/ExampleData/Example1/preproPipelineExample1.cfg'];

%%
addpath('LUMO_toolbox');

%% Covert .LUMO to .nirs
[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs_5wav(LUMODirName);
%s[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs (LUMODirName);

%% Run data quality checks - this produces multiple figures, so comment out for speed.
% DOTHUB_dataQualityCheck(nirsFileName);
% disp('Examine data quality figures, press any key to continue');
% pause 

%% Run Homer2 pre-processing pipeline using .cfg file. Alternatively you can run line by line (as per commented below).

%[prepro, preproFileName] = DOTHUB_runHomerPrepro(nirsFileName,cfgFileName);

%%%%%Equivalent line-by-line Homer2 calls and prepro write:
 dod = hmrIntensity2OD(nirs.d);
 SD3D = enPruneChannels(nirs.d,nirs.SD3D,ones(size(nirs.t)),[0 1e11],12,[0 100],0); 

 %Force MeasListAct to be the same across wavelengths
 SD3D = DOTHUB_balanceMeasListAct(SD3D);

%Set SD2D
 SD2D = nirs.SD; 
 SD2D.MeasListAct = SD3D.MeasListAct;
 
%Bandpass filter and convert to Concentration
 dod = hmrBandpassFilt(dod,nirs.t,0,0.5);
 %%
 ppf = 6;
 ppf_OD2conc = repmat(ppf,1,length(SD2D.Lambda));
 %dc = hmrOD2Conc(dod,SD3D,[6 6]);
 dc = DOTHUB_hmrOD2Conc5wavnew(dod,SD3D,ppf_OD2conc);
 dc = dc*1e6; %Homer works in Molar by default, we use uMolar.
 
 % %%%%..GCVL..#### ************************************************ 
%%  plot concentrations
%%  Concentration Change Figures
	time = nirs.t;
% dc  (datapoints,chromophores,detectors)  chromophores 1=HbO, 2=HbR, 3=HbT, 4=CCO 
% detectors 1->16  (==channels 1->16)    1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16
% 										s1d1  s1d2  s1d3  s1d4  s1d5  s1d6  s1d7  s1d8  s2d1  s2d2  s2d3  s2d4  s2d5  s2d6  s2d7  s2d8
	figure('Name','Conc A','NumberTitle','off');
	t = tiledlayout(2,2);
	t.Title.String = 'Chromophore Concentration Change with Time';
	t.Subtitle.FontSize = 10;
	newcolors = {'red','blue','green'};
	colororder(newcolors);
	nexttile(t);
		plot(time, dc(:,[1 2 4],1));
		title('s1d1');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],3));
		title('s1d3');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],2));
		title('s1d2');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],4));
		title('s1d4');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
 
	figure('Name','Conc B','NumberTitle','off');
	t = tiledlayout(2,2);
	t.Title.String = 'Chromophore Concentration Change with Time';
	t.Subtitle.FontSize = 10;
	newcolors = {'red','blue','green'};
	colororder(newcolors);
	nexttile(t);
		plot(time, dc(:,[1 2 4],5));
		title('s1d5');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],7));
		title('s1d7');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],6));
		title('s1d6');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],8));
		title('s1d8');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');

	figure('Name','Conc C','NumberTitle','off');
	t = tiledlayout(2,2);
	t.Title.String = 'Chromophore Concentration Change with Time';
	t.Subtitle.FontSize = 10;
	newcolors = {'red','blue','green'};
	colororder(newcolors);
	nexttile(t);
		plot(time, dc(:,[1 2 4],9));
		title('s2d1');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],11));
		title('s2d3');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],10));
		title('s2d2');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],12));
		title('s2d4');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');

	figure('Name','Conc D','NumberTitle','off');
	t = tiledlayout(2,2);
	t.Title.String = 'Chromophore Concentration Change with Time';
	t.Subtitle.FontSize = 10;
	newcolors = {'red','blue','green'};
	colororder(newcolors);
	nexttile(t);
		plot(time, dc(:,[1 2 4],13));
		title('s2d5');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],15));
		title('s2d7');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],14));
		title('s2d6');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
	nexttile(t);
		plot(time, dc(:,[1 2 4],16));
		title('s2d8');
		xlabel('Time (s)');
		ylabel('\Delta Concentration (microM)');
		legend({'HbO','HbR','CCO'},'Location','Best');
		
% %%%%..GCVL..#### ************************************************  
