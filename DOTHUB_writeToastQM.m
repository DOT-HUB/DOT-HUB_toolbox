function DOTHUB_writeToastQM(qmfilename,source_pos,det_pos,linklist)

% Function to create and save a qmfile for toast, from a pos and linklist
% Function possibly written by Adam Gibson, UCL?
% RJC April 2020

% Create qm file
qmfile = fopen(qmfilename,'W');
disp ([' - writing to qm file ', qmfilename]);

% write header
dim=size(det_pos,2);
fprintf(qmfile,'QM file %sD\n', num2str(dim));
fprintf(qmfile,'Dimension %s\n', num2str(dim));

% write source list
NSource = size(source_pos,1);

fprintf(qmfile,'\nSourceList %s\n', num2str(NSource));
for i = 1 : NSource
  for j = 1 : dim
    fprintf(qmfile,'%s ', num2str(source_pos(i,j)));
  end
  fprintf(qmfile,'\n', num2str(source_pos(i,j)));
end

% write measurement list
Ndet = size(det_pos,1);
fprintf(qmfile,'\nMeasurementList %s\n', num2str(Ndet));
for i = 1 : Ndet
  for j = 1 : dim
    fprintf(qmfile,'%s ', num2str(det_pos(i,j)));
  end
  fprintf(qmfile,'\n', num2str(det_pos(i,j)));
end

% write link list
ctr=1:NSource;
fprintf(qmfile,'\nLinkList\n');
for i = 1 : NSource
  list=linklist(i,:);
  list=list(list>-1);
  fprintf(qmfile,'%s: %s\n',int2str(length(list)), int2str(list));
end

fclose(qmfile);



