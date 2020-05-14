function mesh = DOTHUB_ply2Mesh(plyFileName)

if ~exist('plyFileName','var')
    [filename, pathname] = uigetfile('*.ply','Select .PLY file...');
    plyFileName = fullfile(pathname, filename);
end

tmp = pcread(plyFileName);
mesh.CData = tmp.Color;
[mesh.node,mesh.face] = read_ply(plyFileName);

