%%
% This is a matlab script to import the void area and plot them.
% Note, must be in the IndividualVoid results folder.

close all;
clear all;

%%
% Import data

meshNoGhostInfRaw = importdata('Mesh/NoGhosts/InfiniteVT/voidArea.dat');
meshNoGhostInfData = meshNoGhostInfRaw.data;

meshNoGhostFinRaw = importdata('Mesh/NoGhosts/FiniteVT/voidArea.dat');
meshNoGhostFinData = meshNoGhostFinRaw.data;

meshGhostRaw = importdata('Mesh/Ghosts/Void/voidArea.dat');
meshGhostData = meshGhostRaw.data;

%%%

NodeDefaultRaw = importdata('Node/DefaultCutOff/Post-Void/voidArea.dat');
NodeDefaultData = NodeDefaultRaw.data;

NodeLargeRaw = importdata('Node/LargeCutoff/Post-Void/voidArea.dat');
NodeLargeData = NodeLargeRaw.data;

NodeDmallRaw = importdata('Node/SmallCutoff/Post-Void/voidArea.dat');
NodeSmallData = NodeDmallRaw.data;

%%%

VertexJaggedRaw = importdata('Vertex/Jagged/voidArea.dat');
VertexJaggedData = VertexJaggedRaw.data;

VertexSmoothRaw = importdata('Vertex/Smooth/voidArea.dat');
VertexSmoothData = VertexSmoothRaw.data;

VertexCurvedRaw = importdata('Vertex/Curved/voidArea.dat');
VertexCurvedData = VertexCurvedRaw.data;

%%
% Clean up and organise data

timeData = meshNoGhostInfData(:,1);

meshNoGhostInfArea = meshNoGhostInfData(:,2);
meshNoGhostFinArea = meshNoGhostFinData(:,3);
meshGhostArea = meshGhostData(:,3);

NodeDefaultArea = NodeDefaultData(:,2);
NodeLargeArea = NodeLargeData(:,2);
NodeSmallArea = NodeSmallData(:,2);

VertexJaggedArea = VertexJaggedData(:,3);
VertexSmoothArea = VertexSmoothData(:,3);
VertexCurvedArea = VertexCurvedData(:,3);

%%
% Plot the void area

mLineWidth = 2.5;

hold on;
plot(timeData,VertexJaggedArea,'-','linewidth',1.5)
plot(timeData,VertexSmoothArea,'-','linewidth',1.5)
plot(timeData,VertexCurvedArea,'-','linewidth',1.5)

plot(timeData,NodeDefaultArea,':','linewidth',2)
plot(timeData,NodeLargeArea,':','linewidth',2)
plot(timeData,NodeSmallArea,':','linewidth',2)

plot(timeData,meshGhostArea,'--','linewidth',1.5)
plot(timeData,meshNoGhostInfArea,'--','linewidth',1.5)
plot(timeData,meshNoGhostFinArea,'--','linewidth',1.5)
hold off;

legend('Vertex Jagged','Vertex Smooth','Vertex Curved','Node Default','Node Large','Node Small','Mesh Ghosts','Mesh Infinite','Mesh Finite','fontsize',12,'interpreter','Latex')

xlabel('Time (seconds)','fontsize',14,'interpreter','Latex')
ylabel('Void Area (CD^2)','fontsize',14,'interpreter','Latex')
title('Void Area with Domain scaling of 0.8','fontsize',16,'interpreter','Latex')
