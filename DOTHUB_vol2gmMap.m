function vol2gm = DOTHUB_vol2gmMap(VolNodes,GMNodes,radius,saveFlag)

%This function performs the mapping from a tetrahedral volume mesh to the
%associated GM mesh.  Output is in the form of a sparse transformation matrix with
%dimensions NxM where M is the number of tetrahedral mesh nodes and N is
%the number of GM surface nodes;

%Defualt radius = 3 mm
if ~exist('radius','var');
    radius = 3;
end

count = 1;
h = waitbar(0,'Calculating vol2gm transform...');
for n = 1:length(GMNodes);
    waitbar(n/length(GMNodes));
    p = GMNodes(n,:);
    ind = find(VolNodes(:,1) < (p(1)+radius) & VolNodes(:,2) < (p(2)+radius) & VolNodes(:,3) < (p(3)+radius) ...
        & VolNodes(:,1) > (p(1)-radius) & VolNodes(:,2) > (p(2)-radius) & VolNodes(:,3) > (p(3)-radius));
    nind = length(ind);
    if nind == 0; 
        err_message = sprintf('GM node %d has zero Head nodes within radius.  Check meshes are aligned or increase radius...',n);
        error(err_message)
        return
    end
    i(count:count+nind-1) = n;
    j(count:count+nind-1) = ind;
    s(count:count+nind-1) = 1/nind;
    count = count+nind;
end
delete(h)

vol2gm = sparse(i,j,s,length(GMNodes),length(VolNodes));

if ~exist('saveFlag','var');
    saveFlag = 1;
end

if saveFlag==1;
    [mapping_filename, mapping_pathname] = uiputfile('*.mat','Save vol2gm mapping as...');
    out_full = [mapping_pathname, mapping_filename ];
    save(out_full,'vol2gm', 'radius');
end
