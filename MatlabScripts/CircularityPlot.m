VertexCircularityJaggedRaw = importdata('GrowingMonolayer/Vertex/Jagged/CicularityCalc.dat');
VertexCircularitySmoothRaw = importdata('GrowingMonolayer/Vertex/Smooth/CicularityCalc.dat');
VertexCircularityCurvedRaw = importdata('GrowingMonolayer/Vertex/Curved/CicularityCalc.dat');
NodeCircularityDefaultRaw = importdata('GrowingMonolayer/Node/Default/CicularityCalc.dat');
NodeCircularityLargeRaw = importdata('GrowingMonolayer/Node/LargeCutoff/CicularityCalc.dat');
NodeCircularitySmallRaw = importdata('GrowingMonolayer/Node/SmallCutoff/CicularityCalc.dat');
MeshCircularityGhostRaw = importdata('GrowingMonolayer/Mesh/Ghosts/CicularityCalc.dat');
MeshCircularityInfRaw = importdata('GrowingMonolayer/Mesh/NoGhosts/InfiniteVT/CicularityCalc.dat');
MeshCircularityFinRaw = importdata('GrowingMonolayer/Mesh/NoGhosts/FiniteVT/CicularityCalc.dat');

%%

VertexCircularityJaggedData = VertexCircularityJaggedRaw.data;
VertexCircularitySmoothData = VertexCircularitySmoothRaw.data;
VertexCircularityCurvedData = VertexCircularityCurvedRaw.data;

NodeCircularityDefaultData = NodeCircularityDefaultRaw.data;
NodeCircularityLargeData = NodeCircularityLargeRaw.data;
NodeCircularitySmallData = NodeCircularitySmallRaw.data;

MeshCircularityGhostData = MeshCircularityGhostRaw.data;
MeshCircularityInfData = MeshCircularityInfRaw.data;
MeshCircularityFinData = MeshCircularityFinRaw.data;

%%

time = VertexCircularityJaggedData(:,1);

Vertex_J_A = VertexCircularityJaggedData(:,2);
Vertex_J_P = VertexCircularityJaggedData(:,3);
Vertex_J_C = VertexCircularityJaggedData(:,4);

Vertex_S_A = VertexCircularitySmoothData(:,2);
Vertex_S_P = VertexCircularitySmoothData(:,3);
Vertex_S_C = VertexCircularitySmoothData(:,4);

Vertex_C_A = VertexCircularityCurvedData(:,2);
Vertex_C_P = VertexCircularityCurvedData(:,3);
Vertex_C_C = VertexCircularityCurvedData(:,4);

Node_D_A = NodeCircularityDefaultData(:,2);
Node_D_P = NodeCircularityDefaultData(:,3);
Node_D_C = NodeCircularityDefaultData(:,4);

Node_L_A = NodeCircularityLargeData(:,2);
Node_L_P = NodeCircularityLargeData(:,3);
Node_L_C = NodeCircularityLargeData(:,4);

Node_S_A = NodeCircularitySmallData(:,2);
Node_S_P = NodeCircularitySmallData(:,3);
Node_S_C = NodeCircularitySmallData(:,4);

Mesh_G_A = MeshCircularityGhostData(:,2);
Mesh_G_P = MeshCircularityGhostData(:,3);
Mesh_G_C = MeshCircularityGhostData(:,4);

Mesh_I_A = MeshCircularityInfData(:,2);
Mesh_I_P = MeshCircularityInfData(:,3);
Mesh_I_C = MeshCircularityInfData(:,4);

Mesh_F_A = MeshCircularityFinData(:,2);
Mesh_F_P = MeshCircularityFinData(:,3);
Mesh_F_C = MeshCircularityFinData(:,4);


%%

figure;
subplot(1,3,1)
hold on;
plot(time,Vertex_J_A,'-','linewidth',1.5)
plot(time,Vertex_S_A,'-','linewidth',1.5)
plot(time,Vertex_C_A,'-','linewidth',1.5)

plot(time,Node_D_A,':','linewidth',1.5)
plot(time,Node_L_A,':','linewidth',1.5)
plot(time,Node_S_A,':','linewidth',1.5)

plot(time,Mesh_G_A,'--','linewidth',1.5)
plot(time,Mesh_I_A,'--','linewidth',1.5)
plot(time,Mesh_F_A,'--','linewidth',1.5)
hold off
xlabel('time')
ylabel('Tissue Area')

subplot(1,3,2)
hold on;
plot(time,Vertex_J_P,'-','linewidth',1.5)
plot(time,Vertex_S_P,'-','linewidth',1.5)
plot(time,Vertex_C_P,'-','linewidth',1.5)

plot(time,Node_D_P,':','linewidth',1.5)
plot(time,Node_L_P,':','linewidth',1.5)
plot(time,Node_S_P,':','linewidth',1.5)

plot(time,Mesh_G_P,'--','linewidth',1.5)
plot(time,Mesh_I_P,'--','linewidth',1.5)
plot(time,Mesh_F_P,'--','linewidth',1.5)
hold off
xlabel('time')
ylabel('Tissue Perimiter')

subplot(1,3,3)
hold on;
plot(time,Vertex_J_C,'-','linewidth',1.5)
plot(time,Vertex_S_C,'-','linewidth',1.5)
plot(time,Vertex_C_C,'-','linewidth',1.5)

plot(time,Node_D_C,':','linewidth',1.5)
plot(time,Node_L_C,':','linewidth',1.5)
plot(time,Node_S_C,':','linewidth',1.5)

plot(time,Mesh_G_C,'--','linewidth',1.5)
plot(time,Mesh_I_C,'--','linewidth',1.5)
plot(time,Mesh_F_C,'--','linewidth',1.5)
hold off
xlabel('time')
ylabel('Tissue Circularity')

legend('Vertex Jagged','Vertex Smooth','Vertex Curved','Node Default','Node Large','Node Small','Mesh Ghosts','Mesh Infinite','Mesh Finite')