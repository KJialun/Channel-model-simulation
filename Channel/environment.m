%% Set environment geometry and material properties.


%% Geometric properties of wall objects:
object_geometry(1,1:6)=[[0,100],[70,100],[0,70]];
object_geometry(2,1:6)=[[0,90],[30,90],[0,30]];
object_geometry(3,1:6)=[[0,35],[30,35],[0,30]];
object_geometry(4,1:6)=[[0,25],[30,25],[0,30]];
object_geometry(5,1:6)=[[0,10],[30,10],[0,30]];
object_geometry(6,1:6)=[[80,75],[110,75],[80,110]];
object_geometry(7,1:6)=[[80,65],[110,65],[80,110]];
object_geometry(8,1:6)=[[80,35],[110,35],[80,110]];
object_geometry(9,1:6)=[[80,25],[110,25],[80,110]];
object_geometry(10,1:6)=[[80,10],[110,10],[80,110]];
object_geometry(11,1:6)=[[70,130],[70,100],[100,130]];
object_geometry(12,1:6)=[[30,90],[30,35],[35,90]];
object_geometry(13,1:6)=[[30,25],[30,10],[10,25]];
object_geometry(14,1:6)=[[80,130],[80,75],[75,130]];
object_geometry(15,1:6)=[[80,65],[80,35],[35,65]];
object_geometry(16,1:6)=[[80,25],[80,10],[10,25]];
object_geometry(17,1:6)=[[0,0],[110,0],[0,110]];
[object_number,~] = size(object_geometry);

%% Geometric properties of building:
building_geometry(1,1:4) = [[0,130],[70,100]];
building_geometry(2,1:4) = [[0,90],[30,35]];
building_geometry(3,1:4) = [[0,25],[30,10]];
building_geometry(4,1:4) = [[80,130],[110,75]];
building_geometry(5,1:4) = [[80,65],[110,35]];
building_geometry(6,1:4) = [[80,25],[110,10]];
% building_geometry(7,1:4) = [[0,0],[110,0]];

[building_number,~] = size(building_geometry);








