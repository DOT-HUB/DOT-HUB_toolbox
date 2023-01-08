% dc = hmrOD2Conc( dod, SD, ppf )
%
% UI NAME:
% OD_to_Conc
%
% dc = hmrOD2Conc( dod, SD, ppf )
% Convert OD to concentrations
%
% INPUTS:
% dod: the change in OD (#time points x #channels)
% SD:  the SD structure. A spatial unit of mm is assumed, but if
%      SD.SpatialUnit = 'cm' then cm will be used.
% ppf: partial pathlength factors for each wavelength. If there are 2
%      wavelengths of data, then this is a vector ot 2 elements.
%      Typical value is ~6 for each wavelength if the absorption change is 
%      uniform over the volume of tissue measured. To approximate the
%      partial volume effect of a small localized absorption change within
%      an adult human head, this value could be as small as 0.1.
%
% OUTPUTS:
% dc: the concentration data (#time points x 3 x #SD pairs
%     3 concentrations are returned (HbO, HbR, HbT)
%

function dc = DOTHUB_hmrOD2Conc5wavnew( dod, SD, ppf )

nWav = length(SD.Lambda);
ml = SD.MeasList;
% %%%%..GCVL..#### ************************************************ 
    sz=size(ml);
    display=sprintf('%d  ',sz);
    fprintf('\n ml size: %s\n', display);
		fprintf('\n ml   80 rows   4 cols  displayed for ease as 80 cols and 4 rows \n');
		fprintf('   1        2        3        4        5        6        7        8        9       10       11       12       13       14       15       16       17       18       19       20       21       22       23       24       25       26       27       28       29       30       31       32       33       34       35       36       37       38       39       40     ');
		fprintf('  41       42       43       44       45       46       47       48       49       50       51       52       53       54       55       56       57       58       59       60       61       62       63       64       65       66       67       68       69       70       71       72       73       74       75       76       77       78       79       80   \n');
		fprintf('s1d1/720 s1d2/720 s1d3/720 s1d4/720 s1d5/720 s1d6/720 s1d7/720 s1d8/720 s2d1/720 s2d2/720 s2d3/720 s2d4/720 s2d5/720 s2d6/720 s2d7/720 s2d8/720 s1d1/760 s1d2/760 s1d3/760 s1d4/760 s1d5/760 s1d6/760 s1d7/760 s1d8/760 s2d1/760 s2d2/760 s2d3/760 s2d4/760 s2d5/760 s2d6/760 s2d7/760 s2d8/760 s1d1/800 s1d2/800 s1d3/800 s1d4/800 s1d5/800 s1d6/800 s1d7/800 s1d8/800 ');
		fprintf('s2d1/800 s2d2/800 s2d3/800 s2d4/800 s2d5/800 s2d6/800 s2d7/800 s2d8/800 s1d1/850 s1d2/850 s1d3/850 s1d4/850 s1d5/850 s1d6/850 s1d7/850 s1d8/850 s2d1/850 s2d2/850 s2d3/850 s2d4/850 s2d5/850 s2d6/850 s2d7/850 s2d8/850 s1d1/890 s1d2/890 s1d3/890 s1d4/890 s1d5/890 s1d6/890 s1d7/890 s1d8/890 s2d1/890 s2d2/890 s2d3/890 s2d4/890 s2d5/890 s2d6/890 s2d7/890 s2d8/890 \n');
		fprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n', ml(:,:));

% %%%%..GCVL..#### ************************************************ 

if length(ppf)~=nWav
    errordlg('The length of PPF must match the number of wavelengths in SD.Lambda');
    dc = zeros(size(dod,1),3,length(find(ml(:,4)==1)));
    return
end

nTpts = size(dod,1);


% %%%%..GCVL..#### ************************************************ 
e = GetExtinctions(SD.Lambda);       %%%%%%%..GCVL..######   note have not used GetExtinctions5wav since they look the same!
%%%  HbO and HbR in GetExtinctions as OD/M/cm and as extinction coefficients (ie. no x2.303) but leave GetExtinctions x2.303 (ie. as absorption coefficients)
%%%  CCO (AA3) in GetExtinctions as OD/milli-Mole/mm and as absorption coefficient (ie. x2.303) and leaves GetExtinctions without change

%Fifth column is AA3=CCO 
		fprintf('\n  e from from GetExtinctions \n    HbO         HbR         HbT    X    CCO \n');
		fprintf(' %f %f %f %f %f \n', e.');

%adjust CCO from OD/milli-Mole/mm to OD/M/cm
e(:,5) = e(:,5)*1000*10;

%%  Overide - Convolution values of e - hand entered    HbO  HbR  HbT=0  X=0  CCO
e = [950.98 3396.78 0 0 4551.48; 1322.15 3421.01 0 0 4418.97; 1925.43 2218.47 0 0 5044.69; 2501.14 1863.10 0 0 5231.78; 2884.31 1919.97 0 0 4724.56];   % Gemma-1 Convolution absorption (ie.x2.303)
%	e = [952.98 3388.63 0 0 4544.91; 1322.28 3418.43 0 0 4417.04; 1924.24 2218.31 0 0 5040.22; 2499.08 1863.50 0 0 5228.00; 2879.21 1920.81 0 0 4719.24];   % Gemma-2 Convolution absorption (ie.x2.303)
%	e = [1130.97 3378.00 0 0 4548.30; 1427.85 3402.92 0 0 4420.35; 1925.80 2214.36 0 0 5044.02; 2388.82 1855.05 0 0 5231.93; 2683.17 1871.94 0 0 4722.78];  % Homer   Convolution absorption (ie.x2.303)


%check Spatial Unit
		fprintf('\n SD.SpatialUnit: %s \n', SD.SpatialUnit);

% 1 is HbO, 2 is HbR and 5 is CCO
if ~isfield(SD,'SpatialUnit')
    e = e(:,[1 2 5]) / 10; % convert from /cm to /mm
elseif strcmpi(SD.SpatialUnit,'mm')
    e = e(:,[1 2 5]) / 10; % convert from /cm to /mm
elseif strcmpi(SD.SpatialUnit,'cm')
    e = e(:,[1 2 5]) ;
end

%  wavelengths used
		fprintf('\n Wavelengths used (nm): %d %d %d %d %d \n', SD.Lambda');
		fprintf(' note: no convolution over spectrum \n', display);

%  check values of extinction coefficients and units 
	    sz=size(e);
		display=sprintf('%d  ',sz);
		fprintf('\n e size: %s \n', display);
		temparry = [SD.Lambda' e];
		fprintf('\n Values of e (absorption coefficients)   "includes 2.303 factor" \n');   % absorption = 2.303 extinction 
		fprintf(' Lambda      HbO        HbR      CCO \n');
		fprintf('  %d      %f %f %f \n', temparry.');
 		fprintf('\n Example values @800 in OD/M/cm   .... ext/abs  \n HbO 865/1992   HbR 840/1934   CCO 2261/5207 \n');

%left inverse (ie. rows > columns) defined as  inverse[A(transposed) * A] *A(transposed)
einv = inv( e'*e )*e';
		fprintf('\n Inverse \n');
	    sz=size(einv);
		display=sprintf('%d  ',sz);
		fprintf(' einv size: %s \n', display);
		fprintf(' %f %f %f %f %f \n', einv.');

% %%%%..GCVL..#### ************************************************ 

% %%%%..GCVL..#### ************************************************ 
		fprintf('\n Detector (channel) grouping by wavelength \n');
lst = find( ml(:,4)==1 );
for idx=1:length(lst)
    idx1 = lst(idx);
    idx2 = find( ml(:,4)>1 & ml(:,1)==ml(idx1,1) & ml(:,2)==ml(idx1,2) );
		fprintf('\n idx %d   idx1 %d   idx2 %d %d %d %d', idx, idx1, idx2);
    rho = norm(SD.SrcPos(ml(idx1,1),:)-SD.DetPos(ml(idx1,2),:));
% dc  (datapoints,chromophores,detectors)  chromophores 1=HbO, 2=HbR, 3=HbT, 4 =CCO
    dc(:,:,idx) = ( einv * (dod(:,[idx1 idx2'])./(ones(nTpts,1)*rho*ppf))' )';
end
% %%%%..GCVL..#### ************************************************ 

	dc(:,4,:) = dc(:,3,:); %Put CCO in fourth column
	dc(:,3,:) = dc(:,1,:) + dc(:,2,:); %Total Haem

% %%%%..GCVL..#### ************************************************ 
    sz=size(dc);
    display=sprintf('%d  ',sz);
    fprintf('\n\n dc size: %s\n', display);
% %%%%..GCVL..#### ************************************************ 

