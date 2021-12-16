

% Code to visualise initial and contours as void closes.
% Note, must be in the IndividualVoid results folder.


addpath /Users/domenicgermano/Workspace/ChasteDom/anim/matlab

how_many_contours = 10;
contour_step = 40;

tissue_line_width = 1.0;
tissue_line_width_node = 0.75;

contour_line_width = 1.4;
contour_node_width = 8;

color_map = GetColourMap(how_many_contours);

%%

figure;

subplot(3,3,6)

ContourData = LoadNonConstantLengthData("Node/SmallCutOff/CircularityContour.dat");


ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

max_x = max(data_x);
min_x = min(data_x);
max_y = max(data_y);
min_y = min(data_y);

max_xy = max(max_x,max_y);
min_xy = max(min_x,min_y);

axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
title('Node Small','FontSize',16,'interpreter','Latex')


%%
subplot(3,3,1)

ContourData = LoadNonConstantLengthData("Vertex/Jagged/CircularityContour.dat");


t_end = size(ContourData,2);

ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

axis(1.1*[min_xy max_xy min_xy max_xy])
pbaspect([1 1 1])
% axis square

% set(gca,'visible','off')
title('Vertex Jagged','FontSize',16,'interpreter','Latex')
%%

subplot(3,3,2)
ContourData = LoadNonConstantLengthData("Vertex/Smooth/CircularityContour.dat");


t_end = size(ContourData,2);

ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end
axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square

% set(gca,'visible','off')
title('Vertex Smooth','FontSize',16,'interpreter','Latex')

%%


subplot(3,3,3)

ContourData = LoadNonConstantLengthData("Vertex/Curved/CircularityContour.dat");

t_end = size(ContourData,2);

ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
title('Vertex Curved','FontSize',16,'interpreter','Latex')

%%

subplot(3,3,4)


ContourData = LoadNonConstantLengthData("Node/Default/CircularityContour.dat");

t_end = size(ContourData,2);

ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
title('Node Default','FontSize',16,'interpreter','Latex')


%%

subplot(3,3,5)

ContourData = LoadNonConstantLengthData("Node/LargeCutOff/CircularityContour.dat");


ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
title('Node Large','FontSize',16,'interpreter','Latex')
%%
% subplot(3,3,6)
% 
% ContourData = LoadNonConstantLengthData("Node/SmallCutOff/CircularityContour.dat");
% 
% 
% ci=0;
% for i = 1:contour_step:(how_many_contours*contour_step)
%     ci = ci + 1;
%     data_1 = cell2mat(ContourData(1,i));
% 
%     data_1(1) = [];
% 
%     data_x = data_1(1:2:end);
%     data_y = data_1(2:2:end);
%     
%     hold on
%     for j=1:2:(length(data_x)-1)
%         plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
%     end
% %     scatter(data_x,data_y,10,'k.')
% 
% end
% 
% % pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
% title('Node Small','FontSize',16,'interpreter','Latex')
% % 

%%

subplot(3,3,7)

ContourData = LoadNonConstantLengthData("Mesh/Ghosts/CircularityContour.dat");

t_end = size(ContourData,2);
ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

axis(1.1*[min_xy max_xy min_xy max_xy])
pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
title('Mesh Ghosts','FontSize',16,'interpreter','Latex')

%%

subplot(3,3,8)
ContourData = LoadNonConstantLengthData("Mesh/NoGhosts/InfiniteVT/CircularityContour.dat");


ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end

axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
title('Mesh Infinite VT','FontSize',16,'interpreter','Latex')

%%

subplot(3,3,9)
ContourData = LoadNonConstantLengthData("Mesh/NoGhosts/FiniteVT/CircularityContour.dat");

t_end = size(ContourData,2);

ci=0;
for i = 1:contour_step:(how_many_contours*contour_step)
    ci = ci + 1;
    data_1 = cell2mat(ContourData(1,i));

    data_1(1) = [];

    data_x = data_1(1:2:end);
    data_y = data_1(2:2:end);
    
    hold on
    for j=1:2:(length(data_x)-1)
        plot([data_x(j) data_x(j+1)],[data_y(j) data_y(j+1)],'k-','linewidth',contour_line_width,'color',color_map(ci,:))
    end
%     scatter(data_x,data_y,10,'k.')

end
axis(1.1*[min_xy max_xy min_xy max_xy])

pbaspect([1 1 1])
% axis square
% set(gca,'visible','off')
% set(findall(gca, 'type', 'text'), 'visible', 'on')

title('Mesh Finite VT','FontSize',16,'interpreter','Latex')


%%
function color_map = GetColourMap(how_many_contours)

    % https://colordesigner.io/gradient-generator
    % using RGB colours:
    % #173f5f -> #ed553b

    if(how_many_contours == 20)
    color_map = 1/255*[23, 63, 95; 34, 64, 93; 46, 65, 91; 57, 66, 89; 68, 68, 87; 79, 69, 86; 91, 70, 84; 102, 71, 82; 113, 72, 80; 124, 73, 78; 136, 75, 76; 147, 76, 74; 158, 77, 72; 169, 78, 70; 181, 79, 68; 192, 80, 67; 203, 82, 65; 214, 83, 63; 226, 84, 61; 237, 85, 59];
    
    elseif(how_many_contours == 19)
    color_map = 1/255*[23, 63, 95; 35, 64, 93; 47, 65, 91; 59, 67, 89; 71, 68, 87; 82, 69, 85; 94, 70, 83; 106, 72, 81; 118, 73, 79; 130, 74, 77; 142, 75, 75; 154, 76, 73; 166, 78, 71; 178, 79, 69; 189, 80, 67; 201, 81, 65; 213, 83, 63; 225, 84, 61; 237, 85, 59];

    elseif(how_many_contours == 18) 
    color_map = 1/255*[23, 63, 95; 36, 64, 93; 48, 66, 91; 61, 67, 89; 73, 68, 87; 86, 69, 84; 99, 71, 82; 111, 72, 80; 124, 73, 78; 136, 75, 76; 149, 76, 74; 161, 77, 72; 174, 79, 70; 187, 80, 67; 199, 81, 65; 212, 82, 63; 224, 84, 61; 237, 85, 59];
    
    elseif(how_many_contours == 17)
    color_map = 1/255*[23, 63, 95; 36, 64, 93; 50, 66, 91; 63, 67, 88; 77, 69, 86; 90, 70, 84; 103, 71, 82; 117, 73, 79; 130, 74, 77; 143, 75, 75; 157, 77, 73; 170, 78, 70; 184, 80, 68; 197, 81, 66; 210, 82, 64; 224, 84, 61; 237, 85, 59];

    elseif(how_many_contours == 16)
    color_map = 1/255*[23, 63, 95; 37, 64, 93; 52, 66, 90; 66, 67, 88; 80, 69, 85; 94, 70, 83; 109, 72, 81; 123, 73, 78; 137, 75, 76; 151, 76, 73; 166, 78, 71; 180, 79, 69; 194, 81, 66; 208, 82, 64; 223, 84, 61; 237, 85, 59];
    
    elseif(how_many_contours == 15)
    color_map = 1/255*[23, 63, 95; 38, 65, 92; 54, 66, 90; 69, 68, 87; 84, 69, 85; 99, 71, 82; 115, 72, 80; 130, 74, 77; 145, 76, 74; 161, 77, 72; 176, 79, 69; 191, 80, 67; 206, 82, 64; 222, 83, 62; 237, 85, 59];
        
    elseif(how_many_contours == 14)
    color_map = 1/255*[23, 63, 95; 39, 65, 92; 56, 66, 89; 72, 68, 87; 89, 70, 84; 105, 71, 81; 122, 73, 78; 138, 75, 76; 155, 77, 73; 171, 78, 70; 188, 80, 67; 204, 82, 65; 221, 83, 62; 237, 85, 59];
        
    elseif(how_many_contours == 13)
    color_map = 1/255*[23, 63, 95; 41, 65, 92; 59, 67, 89; 77, 69, 86; 94, 70, 83; 112, 72, 80; 130, 74, 77; 148, 76, 74; 166, 78, 71; 184, 80, 68; 201, 81, 65; 219, 83, 62; 237, 85, 59];

    elseif(how_many_contours == 12)
    color_map = 1/255*[23, 63, 95; 42, 65, 92; 62, 67, 88; 81, 69, 85; 101, 71, 82; 120, 73, 79; 140, 75, 75; 159, 77, 72; 179, 79, 69; 198, 81, 66; 218, 83, 62; 237, 85, 59];

    elseif(how_many_contours == 11)
    color_map = 1/255*[23, 63, 95; 44, 65, 91; 66, 67, 88; 87, 70, 84; 109, 72, 81; 130, 74, 77; 151, 76, 73; 173, 78, 70; 194, 81, 66; 216, 83, 63; 237, 85, 59];

    elseif(how_many_contours == 10)
    color_map = 1/255*[23, 63, 95; 47, 65, 91; 71, 68, 87; 94, 70, 83; 118, 73, 79; 142, 75, 75; 166, 78, 71; 189, 80, 67; 213, 83, 63; 237, 85, 59];
    
    elseif(how_many_contours == 9)
    color_map = 1/255*[23, 63, 95; 50, 66, 91; 77, 69, 86; 103, 71, 82; 130, 74, 77; 157, 77, 73; 184, 80, 68; 210, 82, 64; 237, 85, 59];
        
    elseif(how_many_contours == 8)
    color_map = 1/255*[23, 63, 95; 54, 66, 90; 84, 69, 85; 115, 72, 80; 145, 76, 74; 176, 79, 69; 206, 82, 64; 237, 85, 59];
    
    elseif(how_many_contours == 7)
    color_map = 1/255*[23, 63, 95; 59, 67, 89; 94, 70, 83; 130, 74, 77; 166, 78, 71; 201, 81, 65; 237, 85, 59];
        
    elseif(how_many_contours == 6)
    color_map = 1/255*[23, 63, 95; 66, 67, 88; 109, 72, 81;	151, 76, 73; 194, 81, 66; 237, 85, 59];

    elseif(how_many_contours == 5)
    color_map = 1/255*[23, 63, 95; 77, 69, 86; 130, 74, 77; 184, 80, 68; 237, 85, 59];
    
    elseif(how_many_contours == 4)
    color_map = 1/255*[23, 63, 95; 94, 70, 83; 166, 78, 71; 237, 85, 59];
        
    elseif(how_many_contours == 3)
    color_map = 1/255*[23, 63, 95; 130, 74, 77; 237, 85, 59];
    
    elseif(how_many_contours == 2)
    color_map = 1/255*[23, 63, 95; 237, 85, 59];
    
    end
    
end


