master = 'GrowingMonolayer_25hrs';

VertexCircularityJaggedRaw = importdata(append(master, '/Vertex/Jagged/CicularityCalc.dat'));
VertexCircularitySmoothRaw = importdata(append(master, '/Vertex/Smooth/CicularityCalc.dat'));
VertexCircularityCurvedRaw = importdata(append(master, '/Vertex/Curved/CicularityCalc.dat'));
NodeCircularityDefaultRaw = importdata(append(master, '/Node/Default/CicularityCalc.dat'));
NodeCircularityLargeRaw = importdata(append(master, '/Node/LargeCutoff/CicularityCalc.dat'));
NodeCircularitySmallRaw = importdata(append(master, '/Node/SmallCutoff/CicularityCalc.dat'));
MeshCircularityGhostRaw = importdata(append(master, '/Mesh/Ghosts/CicularityCalc.dat'));
MeshCircularityInfRaw = importdata(append(master, '/Mesh/NoGhosts/InfiniteVT/CicularityCalc.dat'));
MeshCircularityFinRaw = importdata(append(master, '/Mesh/NoGhosts/FiniteVT/CicularityCalc.dat'));

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
subplot(2,2,1)
hold on;
plot(time,Vertex_J_A,'-','linewidth',1.5)
plot(time,Vertex_S_A,'-','linewidth',1.5)
plot(time,Vertex_C_A,'-','linewidth',1.5)

plot(time,Node_D_A,':','linewidth',2)
plot(time,Node_L_A,':','linewidth',2)
plot(time,Node_S_A,':','linewidth',2)

plot(time,Mesh_G_A,'--','linewidth',1.5)
plot(time,Mesh_I_A,'--','linewidth',1.5)
plot(time,Mesh_F_A,'--','linewidth',1.5)
hold off
xlabel('time','fontsize',12,'interpreter','Latex')
ylabel('Tissue Area','fontsize',12,'interpreter','Latex')

subplot(2,2,2)
hold on;
plot(time,Vertex_J_P,'-','linewidth',1.5)
plot(time,Vertex_S_P,'-','linewidth',1.5)
plot(time,Vertex_C_P,'-','linewidth',1.5)

plot(time,Node_D_P,':','linewidth',2)
plot(time,Node_L_P,':','linewidth',2)
plot(time,Node_S_P,':','linewidth',2)

plot(time,Mesh_G_P,'--','linewidth',1.5)
plot(time,Mesh_I_P,'--','linewidth',1.5)
plot(time,Mesh_F_P,'--','linewidth',1.5)
hold off
xlabel('time','fontsize',12,'interpreter','Latex')
ylabel('Tissue Perimiter','fontsize',12,'interpreter','Latex')

subplot(2,2,3)
hold on;
plot(time,Vertex_J_C,'-','linewidth',1.5)
plot(time,Vertex_S_C,'-','linewidth',1.5)
plot(time,Vertex_C_C,'-','linewidth',1.5)

plot(time,Node_D_C,':','linewidth',2)
plot(time,Node_L_C,':','linewidth',2)
plot(time,Node_S_C,':','linewidth',2)

plot(time,Mesh_G_C,'--','linewidth',1.5)
plot(time,Mesh_I_C,'--','linewidth',1.5)
plot(time,Mesh_F_C,'--','linewidth',1.5)
hold off
xlabel('time','fontsize',12,'interpreter','Latex')
ylabel('Tissue Circularity','fontsize',12,'interpreter','Latex')

legend('Vertex Jagged','Vertex Smooth','Vertex Curved','Node Default','Node Large','Node Small','Mesh Ghosts','Mesh Infinite','Mesh Finite','fontsize',12,'interpreter','Latex')