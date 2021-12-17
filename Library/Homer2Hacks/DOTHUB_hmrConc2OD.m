
% Convert Concentrations to OD for reconstruction purposes
%
% INPUTS:
% dc:   the Average concentration data (#time points x 3 x #SD pairs)
%       or (#time points x 3 x #SD pairs x condition) with 3 concentrations 
%       as (HbO, HbR, HbT). UNITS ARE MOLAR as this is a homer2 function
% SD:   the SD structure
% ppf:  partial pathlength factors for each wavelength. If there are 2
%       wavelengths of data, then this is a vector ot 2 elements.
%       Typical value is ~6 for each wavelength if the absorption change is 
%       uniform over the volume of tissue measured. To approximate the
%       partial volume effect of a small localized absorption change within
%       an adult human head, this value could be as small as 0.1.
%
% OUTPUTS:
% dod:  the OD data (#time points x #channels) or (#time points x #channels x condition)
%
% RJC, UCL

function dod = DOTHUB_hmrConc2OD(dc, SD, ppf)

nWav = length(SD.Lambda);
ml = SD.MeasList;

if length(ppf)~=nWav
    errordlg('The length of PPF must match the number of wavelengths in SD.Lambda');
    dod = [];
    return
end

nTpts = size(dc,1);

e = GetExtinctions( SD.Lambda );
e = e(:,1:2) / 10; % convert from /cm to /mm. This requires input to be in MOLAR
%einv = inv( e'*e )*e';

lst = find( ml(:,4)==1 );

if length(size(dc))==3
    %one condition or non-HRF data
    dod = nan(size(dc,1),size(dc,3)*2);
    
    for idx=1:length(lst)
        idx1 = lst(idx);
        idx2 = find( ml(:,4)>1 & ml(:,1)==ml(idx1,1) & ml(:,2)==ml(idx1,2) );
        rho = norm(SD.SrcPos(ml(idx1,1),:)-SD.DetPos(ml(idx1,2),:));
        dod(:,[idx1 idx2']) = (e * dc(:,1:2,idx)')' .* (ones(nTpts,1)*rho*ppf);
    end
    
elseif length(size(dc))==4
    %multiple conditions
    ncond = size(dc,4);    
    dod = nan(size(dc,1),size(dc,3)*2,ncond);
   
    for iCond = 1:ncond
        for idx=1:length(lst)
            idx1 = lst(idx);
            idx2 = find( ml(:,4)>1 & ml(:,1)==ml(idx1,1) & ml(:,2)==ml(idx1,2) );
            rho = norm(SD.SrcPos(ml(idx1,1),:)-SD.DetPos(ml(idx1,2),:));
            dod(:,[idx1 idx2'],iCond) = (e * dc(:,1:2,idx,iCond)')' .* (ones(nTpts,1)*rho*ppf);
        end
    end
end
