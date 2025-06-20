clear all; clc; close all;

figure;
M = 100;
fid = fopen('SOL_ERROR', 'rb');
data = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error= data(M+1);
t = 1:M+1;
loglog(t, data, 'b-', 'LineWidth', 2);
xlabel("time step");
ylabel("error e_t");
grid on;



figure;
fid = fopen('SOL_ERROR_IM_3_1', 'rb');
data1 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error1= data1(M+1);
loglog(t, data1,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_1");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");



% figure;
% fid = fopen('SOL_ERROR_IM_3_2', 'rb');
% data2 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error2= data2(M+1);
% loglog(t, data,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_2");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_3', 'rb');
% data3 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error3= data3(M+1);
% loglog(t, data,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_3");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_4', 'rb');
% data4 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error4= data4(M+1);
% loglog(t, data4,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_4");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_5', 'rb');
% data5 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error5= data5(M+1);
% loglog(t, data5,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_5");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_6', 'rb');
% data6 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error6= data6(M+1);
% loglog(t, data6,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_6");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_7', 'rb');
% data7 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error7= data7(M+1);
% loglog(t, data7,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_7");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_8', 'rb');
% data8 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error8= data8(M+1);
% loglog(t, data8,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_8");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_9', 'rb');
% data9 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error9= data9(M+1);
% loglog(t, data9,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_9");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");
% 
% figure;
% fid = fopen('SOL_ERROR_IM_3_10', 'rb');
% data10 = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
% fclose(fid);
% error10= data10(M+1);
% loglog(t, data10,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_10");
% xlabel("time step");
% ylabel("error e_t");
% grid on;
% legend("show");

figure;
hold on;
loglog(t, data1,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_1");
% loglog(t, data2,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_2");
% loglog(t, data3,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_3");
% loglog(t, data4,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_4");
% loglog(t, data5,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_5");
% loglog(t, data6,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_6");
% loglog(t, data7,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_7");
% loglog(t, data8,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_8");
% loglog(t, data9,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_9");
% loglog(t, data10,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_10");

xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");
