function dists = DOTHUB_getSDdists(SD)

%This function takes an SD file and spits out the S-D distances for every
%channel.

%RJC UCL, Jan 2020 ############################################################

for i = 1:sum(SD.MeasList(:,4)==1)  
    dists(i) = sqrt(sum((SD.SrcPos(SD.MeasList(i,1),:) - SD.DetPos(SD.MeasList(i,2),:)).^2));
end