%%
% This is a matlab script to import the void area and plot them.
% Note, must be in the IndividualVoid results folder.

close all;
clear all;

%%
% Import data

meshNoGhostRaw = importdata('Mesh/NoGhosts/Void/voidArea.dat');
meshNoGhostData = meshNoGhostRaw.data;

meshGhostRaw = importdata('Mesh/Ghosts/Void/voidArea.dat');
meshGhostData = meshGhostRaw.data;

NodeRaw = importdata('Node/DefaultCutOff/voidArea.dat');
NodeData = NodeRaw.data;

VertexJaggedRaw = importdata('Vertex/Jagged/voidArea.dat');
VertexJaggedData = VertexJaggedRaw.data;

VertexSmoothRaw = importdata('Vertex/Smooth/voidArea.dat');
VertexSmoothData = VertexSmoothRaw.data;

%%
% Clean up and organise data

timeData = meshNoGhostData(:,1);

meshNoGhostVoidArea = meshNoGhostData(:,2);
meshGhostVoidArea = meshGhostData(:,2);
NodeVoidArea = NodeData(:,2);
VertexJaggedVoidArea = VertexJaggedData(:,2);
VertexSmoothVoidArea = VertexSmoothData(:,2);

%%
% Plot the void area

mLineWidth = 2.5;

hold on;

plot(timeData,meshNoGhostVoidArea,'linewidth',mLineWidth)
plot(timeData,meshGhostVoidArea,'linewidth',mLineWidth)
plot(timeData,NodeVoidArea,'--','linewidth',mLineWidth)
plot(timeData,VertexJaggedVoidArea,':','linewidth',mLineWidth)
plot(timeData,VertexSmoothVoidArea,':','linewidth',mLineWidth)

hold off;
legend('Mesh','Mesh Ghosts','Node','Vertex Jagged','Vertex Smooth','fontsize',14,'interpreter','Latex')
xlabel('Time (seconds)','fontsize',14,'interpreter','Latex')
ylabel('Void Area (CD^2)','fontsize',14,'interpreter','Latex')
title('Void Area with Domain scaling of 0.8','fontsize',16,'interpreter','Latex')
