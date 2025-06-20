clear all; clc; close all;

figure;
M = 16000;
fid = fopen('SOL_ERROR_3_TEMPORAL_1', 'rb');
data = fread(fid, [M+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error= data(M+1);
t = 1:M+1;
loglog(t, data, 'b-', 'LineWidth', 2);
xlabel("time step");
ylabel("error e_t");
grid on;


figure;
M1 = 16000;
fid = fopen('SOL_ERROR_3_TEMPORAL_1', 'rb');
data1 = fread(fid, [M1+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error1= data1(M1+1);
t1 = 1:M1+1;
loglog(t1, data1,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_1");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");


figure;
M2 = 18000;
fid = fopen('SOL_ERROR_3_TEMPORAL_2', 'rb');
data2 = fread(fid, [M2+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error2= data2(M2+1);
t2 = 1:M2+1;
loglog(t2, data2,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_2");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M3 = 20000;
fid = fopen('SOL_ERROR_3_TEMPORAL_3', 'rb');
data3 = fread(fid, [M3+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error3= data3(M3+1);
t3 = 1:M3+1;
loglog(t3, data3,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_3");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M4 = 22000;
fid = fopen('SOL_ERROR_3_TEMPORAL_4', 'rb');
data4 = fread(fid, [M4+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error4= data4(M4+1);
t4 = 1:M4+1;
loglog(t4, data4,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_4");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M5 = 24000;
fid = fopen('SOL_ERROR_3_TEMPORAL_5', 'rb');
data5 = fread(fid, [M5+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error5= data5(M5+1);
t5 = 1:M5+1;
loglog(t5, data5,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_5");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M6 = 26000;
fid = fopen('SOL_ERROR_3_TEMPORAL_6', 'rb');
data6 = fread(fid, [M6+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error6= data6(M6+1);
t6 = 1:M6+1;
loglog(t6, data6,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_6");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M7 = 28000;
fid = fopen('SOL_ERROR_3_TEMPORAL_7', 'rb');
data7 = fread(fid, [M7+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error7= data7(M7+1);
t7 = 1:M7+1;
loglog(t7, data7,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_7");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M8 = 30000;
fid = fopen('SOL_ERROR_3_TEMPORAL_8', 'rb');
data8 = fread(fid, [M8+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error8= data8(M8+1);
t8 = 1:M8+1;
loglog(t8, data8,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_8");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M9 = 32000;
fid = fopen('SOL_ERROR_3_TEMPORAL_9', 'rb');
data9 = fread(fid, [M9+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error9= data9(M9+1);
t9 = 1:M9+1;
loglog(t9, data9,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_9");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
M10 = 34000;
fid = fopen('SOL_ERROR_3_TEMPORAL_10', 'rb');
data10 = fread(fid, [M10+1, 1], 'double');  % 按列读取为向量
fclose(fid);
error10= data10(M10+1);
t10 = 1:M10+1;
loglog(t10, data10,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_10");
xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");

figure;
hold on;
loglog(t1, data1,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_1");
loglog(t2, data2,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_2");
loglog(t3, data3,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_3");
loglog(t4, data4,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_4");
loglog(t5, data5,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_5");
loglog(t6, data6,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_6");
loglog(t7, data7,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_7");
loglog(t8, data8,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_8");
loglog(t9, data9,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_9");
loglog(t10, data10,  'LineWidth', 2, "DisplayName", "SOL\_ERROR\_10");

xlabel("time step");
ylabel("error e_t");
grid on;
legend("show");
