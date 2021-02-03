
clear all
close all

A = readtable('R_0_Surface_Plot_Data_for_Revision.csv');
A_sorted = sortrows(A, 'R_0');


A_mat = table2array(A_sorted);
A_mat;
R_0 = A_mat(:,1);
b_a = A_mat(:,2);
b_p = A_mat(:,3);
class(A)
class(R_0)
%plot3(A)
% plot3(b_a, b_p, R_0, '.')
%surf(b_a, b_p, R_0)

% scatter3(b_a, b_p, R_0, 'r', 'filled')

 color_vec = linspace(1,3, numel(R_0));
% numel(R_0)
% C = repmat([1],numel(R_0),1);
% C
% c = C(:)
 size_big = 100;
% scatter3(b_a, b_p, R_0,size_big,color_vec, 'filled')
% xlabel('b_p');
% ylabel('b_a');
% zlabel('R_0');
% view(40,35)
% saveas(gcf,'R_0_Surface_Plot_Matlab.png');


B = readtable('R_0_NGM_Surface_Plot_Data_for_Revision.csv');
B_sorted = sortrows(B, 'R_0_NGM');
B;
B_mat = table2array(B_sorted);
B_mat;
R_0_NGM = B_mat(:,1);
b_a_B = B_mat(:,2);
b_p_B = B_mat(:,3);

color_vec_B = linspace(1,3, numel(R_0_NGM));
%scatter3(b_a_B, b_p_B, R_0_NGM,size_big,color_vec_B, 'filled')
%xlabel('b_p');
%ylabel('b_a');
%zlabel('R_0_NGM');
%view(40,35)
%createfigure(b_a_B, b_p_B, R_0_NGM, size_big, color_vec_B)
%createfigure1(b_a_B, b_p_B, R_0_NGM, size_big, color_vec_B)

figure1 = figure;


%%% PANEL A %%%%%

subplot(1,3,1);
axes1 = subplot(1,3,1);
%axes1.FontSize = 32; 

scatter(b_p, b_a, size_big, R_0,'MarkerFaceColor','flat','MarkerEdgeColor','none');
%cb = colorbar; % creates the colorbar on side
%colorTitleHandle = get(cb,'Title');
%colorBarPosition = get(cb,'Position');

%colorBarPosition = cb.Position
%cpbA_l = colorBarPosition(1)
%cpbA_b = colorBarPosition(2)
%cpbA_w = colorBarPosition(3)
%cpbA_h = colorBarPosition(4)

%CPB_scaling = 0.5
%new_CPB_height = cpbA_h*CPB_scaling

% lbwh_3a = axes1.Position;
% fig3a_l = lbwh_3a(1);
% fig3a_b = lbwh_3a(2);
% fig3a_w = lbwh_3a(3);
% fig3a_h = lbwh_3a(4);
% 
% left_gap = 0.01;
% bottom_offset = 0.27;
% left_pos = fig3a_l + fig3a_w +left_gap;
% cb_a_bottom_pos = bottom_offset + fig3a_b;
%pos2 = [0.2757 0.1100 0.0138 0.4075];
%pos2 = [0.35 0.3775 0.0138 0.2775];
% cb_width = 0.0138;
% cb_height = 0.2775;
% pos_a_cb = [left_pos cb_a_bottom_pos  cb_width cb_height];

cb = colorbar();
colorTitleHandle = get(cb,'Title');

titleString_A = 'R_{0}';
set(colorTitleHandle ,'String',titleString_A, 'FontSize',24);

%colorBarPosition = get(cb,'Position');

%set(colorBarPosition, [0.2757 0.1100 0.0138 0.4075])
axis square;
% Create xlabel
xlabel('b_p','FontSize',24);

% Create ylabel
ylabel('b_a','FontSize',24);

% Create zlabel
zlabel('R_{0}','FontSize',24);



%view(axes1,[-130.8 15]);
%view(axes1,[134.4 8.59999999999994]);
grid(axes1,'on');
t = title('A)', 'Units', 'normalized', 'Position', [-0.2, 1.215, 0], 'FontSize',28 ); % Set Title with correct Position


%%% PANEL B %%%%

%subplot(1,3,2);
% gap_from_a_cb = 0.075;
% plot_b_l = left_pos + cb_width + gap_from_a_cb;
% fig_3b_bottom_offset = 0.005
% fig3_b_set = fig3a_b + fig_3b_bottom_offset;
% pos_plot_B = [plot_b_l fig3_b_set  fig3a_w fig3a_h];

% Create axes
%axes2 = subplot(1,3,2);
% axes2 = subplot('Position',pos_plot_B);

%hold(axes1,'on');
%axes2.FontSize = 32; 

% Create scatter3
axes2 = subplot(1,3,2);
scatter(b_a_B, b_p_B, size_big, R_0_NGM, 'MarkerFaceColor','flat','MarkerEdgeColor','none');

% left_pos_b_cb = plot_b_l + fig3a_w +left_gap;
% pos_b_cb = [left_pos_b_cb cb_a_bottom_pos  cb_width cb_height];

cb = colorbar();

colorTitleHandle = get(cb,'Title');
titleString_B = 'R_{0_{NGM}}';
set(colorTitleHandle ,'String',titleString_B, 'FontSize',24);
axis square;

% Create xlabel
xlabel('b_a','FontSize',24);

% Create ylabel
ylabel('b_p','FontSize',24);

% Create zlabel
zlabel('R_{0_{NGM}}','FontSize',24);



% view(axes2,[-130.8 15]);
%view(axes2,[134.4 8.59999999999994]);
grid(axes2,'on');
hold on;
t = title('B)', 'Units', 'normalized', 'Position', [-0.2, 1.220, 0], 'FontSize',28 ); % Set Title with correct Position

% lbwh_3b = axes2.Position;
% fig3b_l = lbwh_3b(1);
% fig3b_b = lbwh_3b(2);
% fig3b_w = lbwh_3b(3);
% fig3b_h = lbwh_3b(4);
% fig3a_l;
% fig3a_b;
% fig3a_w;
% fig3a_h;

% wh_ratio = fig3b_w/fig3b_h;
% gap_w = fig3b_l - fig3a_l - fig3a_w - left_gap - cb_width;
% expected_left_pos_3c = fig3b_l + fig3b_w + gap_w + cb_width;

%%% PANEL C %%%%%

%subplot(1,3,3, [expected_left_pos_3c,fig3b_b, fig3b_w, fig3b_h]);
% scaling_3c = .70;
% fig_3c_height = scaling_3c*fig3b_h;
% fig_3c_width = scaling_3c*fig3b_w;
% fig3c_b = fig3b_b + .1175;
% pos_plot_3c = [expected_left_pos_3c,fig3c_b, fig_3c_width, fig_3c_height];
% axes3 = subplot('Position',pos_plot_3c);
%axes3 = subplot(1,3,3);
%axes3.FontSize = 32;
axes3 = subplot(1,3,3);
C_data_table = readtable('product_data_for_Fig_3C.csv');
C_sorted = sortrows(C_data_table, 'b_a');
C_data_table_R_0_NGM = C_sorted.R_0_NGM;
C_data_table_R_0 = C_sorted.R_0;
C_data_table_b_a = C_sorted.b_a;
class(C_data_table_b_a)
%C_data_R_0_NGM = table2array(C_data_table_R_0_NGM);
%C_data_R_0 = table2array(C_data_table_R_0);
%C_data_b_a = table2array(C_data_table_b_a);

color_vec_C = linspace(1,3, numel(C_data_table_b_a));

scatter(C_data_table_R_0, C_data_table_R_0_NGM, size_big, C_data_table_b_a, 'MarkerFaceColor','flat','MarkerEdgeColor','none')
t = title('C)', 'Units', 'normalized', 'Position', [-0.2, 1.220, 0], 'FontSize',28 ); % Set Title with correct Position

% 
% lbwh_3c = axes3.Position;
% fig3c_l = lbwh_3c(1);
% fig3c_b = lbwh_3c(2);
% fig3c_w = lbwh_3c(3);
% fig3c_h = lbwh_3c(4);
% 
% left_pos_c_cb = fig3c_l + fig3c_w +left_gap;
% pos_c_cb = [left_pos_c_cb cb_a_bottom_pos  cb_width cb_height];

cb = colorbar();

%cb = colorbar; % creates the colorbar on side
colorTitleHandle = get(cb,'Title');
titleString_C = 'b_a';
set(colorTitleHandle ,'String',titleString_C, 'FontSize',24);

caxis([0,1]); % low end is 0, high end is 1

 axis square;
%pbaspect([.50 .50 .50])
%title('C)') 
% Create xlabel
xlabel('R_{0}','FontSize',24);

% Create ylabel
ylabel('R_{0_{NGM}}','FontSize',24);


%saveas(gcf,'Figure_3_Revised_Matlab');

%saveas(gcf,'Figure_3_Revised_Matlab','epsc');


