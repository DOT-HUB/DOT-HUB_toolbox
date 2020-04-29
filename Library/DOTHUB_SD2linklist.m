function linklist = DOTHUB_SD2linklist(SD)

%This function takes the SD.MeasList variable and creates the toast++
%equivalent.

%RJC, UCL April 2020
nwav = length(SD.Lambda);
tmp = SD.MeasList(1:end/nwav,1:2);
nsrc = length(unique(tmp(:,1)));

%Find max number of detectors for a given source
i = 1;
count = 1;
while i <= length(tmp)
    tmp2 = sum(tmp(:,1) == tmp(i));
    ndet(count) = tmp2;
    i = i + tmp2;
    count = count+1;
end

max_det = max(ndet);

linklist = -1*ones(nsrc,max_det);

count = 1;
for i = 1:nsrc
    for j = 1:ndet(i)
        linklist(i,j) = tmp(count,2) - 1; %base 0
        count = count+1;
    end
end

