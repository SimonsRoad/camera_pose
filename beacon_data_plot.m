% BLE signals
clear;
close all
addpath(genpath('SensorFusion'))
ble_reading = [[-77 -72 -78 -83 -83 -75 -75 -88];
    [ -77 -72 -78 -83 -85 -75 -75 -88];
    [ -76 -85 -79 -82 -85 -76 -76 -86];
    [ -71 -71 -74 -83 -85 -74 -81 -75];
    [ -71 -72 -74 -83 -77 -74 -80 -79];
    [ -71 -75 -75 -84 -77 -75 -78 -79];
    [ -67 -74 -75 -82 -87 -73 -81 -79];
    [ -67 -73 -75 -81 -87 -73 -80 -86];
    [ -67 -73 -75 -80 -87 -72 -79 -86];
    [ -70 -78 -79 -84 -84 -77 -78 -82];
    [ -70 -78 -80 -84 -84 -78 -79 -80];
    [ -69 -78 -80 -84 -84 -78 -81 -79];
    [ -71 -79 -82 -79 -86 -73 -76 -78];
    [ -71 -73 -83 -75 -84 -72 -74 -83];
    [ -71 -72 -81 -75 -85 -72 -74 -83];
    [ -75 -75 -83 -85 -83 -82 -82 -79];
    [ -75 -78 -83 -85 -83 -79 -82 -80];
    [ -75 -82 -82 -85 -81 -78 -80 -80];
    [ -72 -69 -77 -88 -86 -76 -81 -82];
    [ -72 -68 -76 -88 -86 -76 -82 -82];
    [ -72 -72 -79 -89 -87 -78 -82 -80];
    [ -71 -68 -77 -85 -86 -80 -79 -81];
    [ -70 -67 -77 -85 -86 -79 -77 -81];
    [ -70 -66 -78 -82 -86 -75 -76 -81];
    [ -70 -76 -79 -83 -88 -82 -75 -81];
    [ -70 -76 -79 -83 -88 -83 -76 -81];
    [ -69 -78 -79 -83 -88 -85 -76 -87];
    [ -74 -75 -80 -77 -88 -81 -74 -81];
    [ -74 -74 -79 -76 -88 -81 -74 -80];
    [ -75 -76 -79 -83 -88 -82 -74 -80];
    [ -74 -81 -77 -88 -92 -83 -78 -78];
    [ -71 -81 -77 -88 -92 -83 -78 -79];
    [ -69 -80 -77 -88 -92 -83 -79 -80];
    [ -71 -75 -78 -82 -82 -80 -77 -78];
    [ -71 -74 -77 -86 -87 -80 -77 -78];
    [ -73 -73 -78 -86 -87 -80 -78 -78];
    [ -70 -76 -77 -84 -87 -82 -83 -77];
    [ -70 -75 -78 -85 -86 -83 -83 -78];
    [ -71 -74 -79 -85 -86 -83 -80 -78];
    [ -69 -75 -79 -81 -86 -83 -81 -76];
    [ -68 -74 -81 -81 -86 -83 -82 -74];
    [ -67 -75 -82 -80 -86 -85 -84 -74];
    [ -69 -74 -80 -81 -82 -78 -84 -77];
    [ -69 -72 -80 -82 -82 -78 -84 -80];
    [ -70 -71 -80 -82 -85 -79 -85 -80];
    [ -76 -78 -77 -83 -86 -72 -77 -77];
    [ -73 -79 -77 -83 -86 -71 -76 -77];
    [ -71 -82 -77 -83 -86 -70 -76 -77]];



% camera based angle / location
cam_loc = [[ -9.2339696E-4 2.423786 -0.026260108  0.018100603 -0.38445407 2.6190554E-4];
    [ 0.08403756 18.96005 0.10337092  0.018100603 -0.38445407 2.6190554E-4];
    [ -0.08060429 22.926672 0.31273937  7.4231205 -0.24623072 0.032897953]; 
    [ -0.32554293 19.363302 0.6310724  27.416714 -0.18642882 0.0171329]; %%
        [ -0.32554293 19.363302 0.6310724  27.416714 -0.18642882 0.0171329];%%

    %[ -0.92860305 20.351578 1.801031  438.89847 33.809967 2.8451867];
    [ 0.081992865 22.96344 0.10043938  3.4926744 -0.063796826 0.00978692]; %%
    [ -0.07558699 24.478693 0.38540414  3.951764 -0.13259389 0.013833043]; %%
    [ -0.57197344 -23.824858 -0.16618489  1.4326109 -0.16934341 -0.0022800544];
    
    [ -0.14147848 28.86289 0.49260294  5.905164 -0.1985078 0.024412185];
    [ 0.17160578 48.417896 0.34554055  3.5756729 -0.372201 0.009232345];
    [ 0.20494728 8.1090355 -0.2721115  1.2512637 -0.1671736 -0.0018602774];
    [ 0.3449509 29.734941 0.18381134  2.4018307 -0.082891695 0.0030466507];
    [ 0.10884422 29.950478 0.36411688  3.152513 -0.033178337 0.007391646]; %%
    [ 0.10884422 29.950478 0.36411688  3.152513 -0.033178337 0.007391646]; %%
    [ 0.10884422 29.950478 0.36411688  3.152513 -0.033178337 0.007391646]; %%
 [ 0.073661275 -26.065016 -0.43022397  30.184032 -0.33620973 0.017412922]]; %%

cam_loc_true = cam_loc;
cam_loc_true(2,:) = cam_loc_true(1,:) + cam_loc_true(2,:);
cam_loc_true(3,:) = cam_loc_true(2,:) + cam_loc_true(3,:);
cam_loc_true(4,:) = cam_loc_true(3,:) + cam_loc_true(4,:);
cam_loc_true(5,:) = cam_loc_true(4,:) + cam_loc_true(5,:);

cam_loc_true(6,:) = cam_loc_true(5,:) + cam_loc_true(6,:);
cam_loc_true(7,:) = cam_loc_true(6,:) + cam_loc_true(7,:);
cam_loc_true(8,:) = cam_loc_true(7,:) + cam_loc_true(8,:);
cam_loc_true(9,:) = cam_loc_true(8,:) + cam_loc_true(9,:);
cam_loc_true(10,:) = cam_loc_true(9,:) + cam_loc_true(10,:);
cam_loc_true(11,:) = cam_loc_true(10,:) + cam_loc_true(11,:);
cam_loc_true(12,:) = cam_loc_true(11,:) + cam_loc_true(12,:);
cam_loc_true(13,:) = cam_loc_true(12,:) + cam_loc_true(13,:);
cam_loc_true(14,:) = cam_loc_true(13,:);
cam_loc_true(15,:) = cam_loc_true(13,:);
cam_loc_true(16,:) = cam_loc_true(15,:) + + cam_loc_true(16,:);

cam_loc_true(:,5) = -1 * cam_loc_true(:,5);

gc_locations_2d = [[0,0];
    [8.8/3, 10.1];
    [2*8.8/3, 20.2];
    [8.8,30.2];
    [8.8+13/3,30.2+26/3];
    [8.8+26/3,30.2+52/3];
    [21.8,56.2];
    [21.8+16/3,56.2+23/3];
    [21.8+32/3,56.2+46/3];
    [37.8,79.2];
    [37.8+15.5/3,79.2+25/3];
    [37.8+31/3,79.2+50/3];
    [53.3,104.2];
    [53.3+32/3,104.2+9.5/3];
    [53.3+64/3,104.2+19/3];
    [85.3,113.7]];

% process ble location estimation
real_ble_signal = zeros(16, 8);
for i = 1:16
    real_ble_signal(i,:) = mean(ble_reading((i-1)*3+1:i*3,:),1);
end

txpower = repmat([-66, -76, -66, -76, -81, -81, -76, -81], [16,1]);
dist = 10.^((txpower - real_ble_signal) / 20);

ble_locations = [];
for i = 1:16
    dist0 = dist(i,:);
    subplot(4,4,i)
    diffs_map = adhoc(dist0)
    [y,x] = find(min(min(diffs_map))==diffs_map);
    ble_locations = [ble_locations; [x,y]]; % y, x
    imagesc(diffs_map)
    title(num2str(i))
    colormap('hot')
end

ble_locations_scale = ble_locations*3;
ble_locations_scale(:,1) = (75-ble_locations_scale(:,1)) * 3
ble_locations_scale(:,2) = (105-ble_locations_scale(:,2)) * 2.5

ble_locations_scale7 = ble_locations_scale(7,:)
ble_locations_scale8 = ble_locations_scale(8,:)
ble_locations_scale(7,:) = ble_locations_scale(15,:)
ble_locations_scale(8,:) = ble_locations_scale(16,:)
ble_locations_scale(15,:) = ble_locations_scale7
ble_locations_scale(16,:) = ble_locations_scale8

cam_locations = cam_loc_true(:, 4:6);

cam_locations_2d = cam_loc_true(:, 5:6);

cam_locations_2d(:,1) = cam_locations_2d(:,1) * 10 - 5; 
cam_locations_2d(:,2) = cam_locations_2d(:,2) * 280;

cam_locations_scale = cam_locations_2d * 3; 
% 
% 
% figure;scatter(ble_locations_scale(:,1),ble_locations_scale(:,2))
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(ble_locations_scale(:,1)+dx, ble_locations_scale(:,2)+dy, c)
% 
% 
% figure; scatter(cam_locations_scale(:,1), cam_locations_scale(:,2))
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
% text(cam_locations_scale(:,1)+dx, cam_locations_scale(:,2)+dy, c)
% 
% figure;
% scatter(gc_locations_2d(:,1), gc_locations_2d(:,2))
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(gc_locations_2d(:,1)+dx, gc_locations_2d(:,2)+dy, c)
% axis equal
% 
% figure;
% plot(sqrt((cam_locations_scale(:,1) -  gc_locations_2d(:,1)).^2 + (cam_locations_scale(:,2) -  gc_locations_2d(:,2)).^2))

bledist_heatmap = reshape(sqrt((ble_locations_scale(:,1) -  gc_locations_2d(:,1)).^2 + (ble_locations_scale(:,2) -  gc_locations_2d(:,2)).^2), 4,4);
camdist_heatmap = reshape(sqrt((cam_locations_scale(:,1) -  gc_locations_2d(:,1)).^2 + (cam_locations_scale(:,2) -  gc_locations_2d(:,2)).^2), 4,4);
% figure;
% imagesc(bledist_heatmap, [0,55])
% colormap('hot')
% figure;
% imagesc(camdist_heatmap, [0,55])
% colormap('hot')


% kalman


% Biases
BIAS1 = 0;
BIAS2 = 0;

% Measurement noise covariances
R1 = 0.002;
R2 = 0.64;

% Process noise covariance
Q = .005;

% State transition model
A = 1;

% Observation model
C1 = 1;
C2 = 1;


camxhatx = kalman(cam_locations_scale(:,1)', A, C1, R1, Q);
camxhaty = kalman(cam_locations_scale(:,2)', A, C1, R1, Q);

%blecamxhatx = kalman([cam_locations_scale(:,1)'; ble_locations_scale(:,1)'], A, [C1; C2], [R1 0; 0 R2], Q);
%blecamxhaty = kalman([cam_locations_scale(:,2)'; ble_locations_scale(:,2)'], A, [C1; C2], [R1 0; 0 R2], Q);

cam_locations_acale_x = cam_locations_scale(:,1);
cam_locations_acale_y = cam_locations_scale(:,2);
ble_locations_scale_x = ble_locations_scale(:,1);
ble_locations_scale_y = ble_locations_scale(:,2);


blecamxhatx1 = kalman([cam_locations_acale_x(1:13)'; ble_locations_scale_x(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
blecamxhatx2 = kalman([blecamxhatx1(13), ble_locations_scale_x(14:15)'], A, C2, R2, Q);
blecamxhatx3 = kalman([[blecamxhatx2(3), cam_locations_acale_x(16)']; [blecamxhatx2(3), ble_locations_scale_x(16)']],  A, [C1; C2], [R1 0; 0 R2], Q);

blecamxhaty1 = kalman([cam_locations_acale_y(1:13)'; ble_locations_scale_y(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
blecamxhaty2 = kalman([blecamxhaty1(13), ble_locations_scale_y(14:15)'], A, C2, R2, Q);
blecamxhaty3 = kalman([[blecamxhaty2(3), cam_locations_acale_y(16)']; [blecamxhaty2(3), ble_locations_scale_y(16)']],  A, [C1; C2], [R1 0; 0 R2], Q);

%blecamxhatanglez1 = kalman([imu_angles_avg_z(1:13)'; cam_angle_z(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
%blecamxhatanglez2 = kalman([blecamxhatanglez1(13), imu_angles_avg_z(14:15)'], A, C1, R1, Q);
%blecamxhatanglez3 = kalman([[blecamxhatanglez2(3), imu_angles_avg_z(16)']; [blecamxhatanglez2(3), cam_angle_z(16)']],  A, [C1; C2], [R1 0; 0 R2], Q);

blecamxhatx = [blecamxhatx1 blecamxhatx2(2:3) blecamxhatx3(2)]
blecamxhaty = [blecamxhaty1 blecamxhaty2(2:3) blecamxhaty3(2)]
%blecamxhatanglez = [blecamxhatanglez1 blecamxhatanglez2(2:3) blecamxhatanglez3(2)]



% subplot(2,3,1)
% scatter(ble_locations_scale(:,1),ble_locations_scale(:,2))
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(ble_locations_scale(:,1)+dx, ble_locations_scale(:,2)+dy, c)
% title('ble')
% 
% subplot(2,3,2)
% scatter(cam_locations_scale(:,1), cam_locations_scale(:,2))
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
% text(cam_locations_scale(:,1)+dx, cam_locations_scale(:,2)+dy, c)
% title('cam')
% 
% subplot(2,3,3)
% scatter(gc_locations_2d(:,1), gc_locations_2d(:,2))
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(gc_locations_2d(:,1)+dx, gc_locations_2d(:,2)+dy, c)
% title('gc')
% 
% 
% subplot(2,3,4)
% 
% scatter(camxhatx,camxhaty)
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(camxhatx+dx, camxhaty+dy, c)
% title('cam k')
% 
% subplot(2,3,5)
% 
% scatter(blecamxhatx, blecamxhaty)
% a = [1:16]'; b = num2str(a); c = cellstr(b);
% dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
% text(blecamxhatx+dx, blecamxhaty+dy, c)
% title('ble cam k')


bledist_heatmap = reshape(sqrt((ble_locations_scale(:,1) -  gc_locations_2d(:,1)).^2 + (ble_locations_scale(:,2) -  gc_locations_2d(:,2)).^2), 4,4);
camdist_heatmap = reshape(sqrt((cam_locations_scale(:,1) -  gc_locations_2d(:,1)).^2 + (cam_locations_scale(:,2) -  gc_locations_2d(:,2)).^2), 4,4);

%camk_heatmap = reshape(sqrt((camxhatx' -  gc_locations_2d(:,1)).^2 + (camxhaty' -  gc_locations_2d(:,2)).^2), 4,4);
camblek_heatmap = reshape(sqrt((blecamxhatx' -  gc_locations_2d(:,1)).^2 + (blecamxhaty' -  gc_locations_2d(:,2)).^2), 4,4);

%figure
%plot(1:16,sqrt((blecamxhatx' -  gc_locations_2d(:,1)).^2 + (blecamxhaty' -  gc_locations_2d(:,2)).^2),'r', 1:16, sqrt((cam_locations_scale(:,1) -  gc_locations_2d(:,1)).^2 + (cam_locations_scale(:,2) -  gc_locations_2d(:,2)).^2),'b')



% TODO: what is difference is z? 


% angles 


% from 06 it is not accumulate ()
% 
cam_angle = [[ -0.15362704 11.263651 0.103366  15.0 50.0 500.0];
    [ 0.12336848 21.164955 0.081830874  15.0 50.0 500.0];
    [ 0.75403374 -19.507513 0.14834695  15.0 50.0 500.0];
    [ 0.77263016 -13.56495 0.042694595  15.0 50.0 500.0];
    [ 0.7351515 2.8172398 0.07982737  15.0 50.0 500.0];
    [ 0.065434866 5.570113 -0.022674432  15.0 50.0 500.0];
    [ 0.08237838 7.71051 -0.042914893  15.0 50.0 500.0];
    [ 0.10355768 17.802885 0.2130557  15.0 50.0 500.0];
    [ 0.10355768 17.802885 0.2130557  15.0 50.0 500.0];
    [ 0.18918112 14.987406 0.14623751  15.0 50.0 500.0];
    [ -0.14852045 44.490932 0.587198  15.0 50.0 500.0];
    [ -0.14852045 44.490932 0.587198  15.0 50.0 500.0];
    [ -0.07901941 24.232496 0.4033664  15.0 50.0 500.0];
    [ -0.08444327 16.5776 0.35947388  15.0 50.0 500.0]];

% IMU angles z x y 
imu_angles = [[ -0.2529885 0.21062565 0.020345522];
    [ 0.13223526 -0.027317883 0.2984652];
    [ 0.19500843 -0.0061449576 0.1900929];
    [ -0.29146838 -4.628151 0.33439702];
    [ -0.46481442 -3.9063814 0.2699143];
    [ -0.44650343 -3.7215648 0.24296951];
    [ 0.24675737 -12.728138 0.4801182];
    [ 0.20222095 -13.01475 0.42366895];
    [ 0.14661296 -14.430077 0.43418285];
    [ -0.23740987 -6.782422 0.40305024];
    [ -0.20036817 -6.971632 0.44965798];
    [ -0.25260884 -6.8805723 0.29463938];
    [ -1.9658544 -7.985493 -1.2096393];
    [ -2.0825007 -6.9488945 -1.1448766];
    [ -2.271693 -5.8986716 -1.1739595];
    [ -0.8894609 -7.74281 -0.16730826];
    [ -0.92374444 -8.192009 -0.19530527];
    [ -1.1548778 -7.280638 0.014321296];
    [ -1.5168375 -16.323744 0.2184689];
    [ -0.78397644 -16.119493 0.15173472];
    [ -1.1811887 -16.14906 0.2575159];
    [ -3.7028985 -18.107542 -6.0806355];
    [ -3.9141827 -16.180109 -6.2341547];
    [ -4.074766 -14.883991 -6.6097345];
    [ -3.487239 -10.472566 -0.3114105];
    [ -3.6231582 -10.778939 -0.32490885];
    [ -3.778816 -9.710013 -0.34727272];
    [ -2.1772602 -24.583483 3.0524464];
    [ -1.9083 -24.71984 3.3074474];
    [ -2.6444864 -24.156118 3.1661243];
    [ -7.9799113 -25.098932 3.4518073];
    [ -8.215604 -25.134586 3.3187556];
    [ -8.631085 -23.817661 3.3966303];
    [ -4.707316 -31.212608 0.9230426];
    [ -4.8201947 -32.341015 1.4375597];
    [ -4.892433 -31.713493 1.3984948];
    [ -4.892433 -31.713493 1.3984948];%
    [ -4.892433 -31.713493 1.3984948];%
    [ -4.892433 -31.713493 1.3984948];%
    [ -1.6892663 -52.366497 3.1285443];
    [ -1.7474583 -51.52763 3.3689215];
    [ -1.8380799 -51.217495 3.2992787];
    [   -12.5217 -83.081566 0.0159];
    [  -12.5957 -73.81381 0.3383];
    [  -12.8941 -70.99001 0.6581];
    [ -12.467782 -19.486094 0.21071947];
    [ -11.777315 -18.020823 0.84060925];
    [ -12.075464 -17.795618 0.7246978]];

%figure
imu_angles_avg = zeros(16,3);
for i = 1:16
    imu_angles_avg(i,:) = mean(imu_angles((i-1)*3+1: i*3,:),1);
end

imu_angles_avg_z = [ 0.0248   -0.4009    0.1985   -0.2301   -2.1067   -0.9894   -1.1607   -3.8973   -3.6297   -2.2433  4.8066   8.2755     -4.8924   -1.7583  -2.6705,   -2.1069]'
gc_z = [0, -2.5, 0, 0, 0, -2, 0, 0, -3.5, 0, 3, 6.5, -1, 0, -0.5, 0]';

%plot(1:16,gc_z,'g', 1:16, imu_angles_avg_z,'b');
imu_diff_angles_z = gc_z - imu_angles_avg_z;
%figure;
%plot(imu_diff_angles_z)

imu_angles_avg_x = [0.0591   -4.0854  -13.3910   -6.8782   -6.9444   -7.7385  -16.1974  -16.3905  -10.3205  -24.4865  -24.6837  -31.7557  -31.7135  -51.7039  -75.9618   -12.1069]';
gc_x = [0 -3 -13 -3 -5 -5 -18 -14 -14 -22 -22 -32 -30 -50 -71 -10]';
%figure;
%plot(1:16,gc_x,'g', 1:16, imu_angles_avg_x,'b');
imu_diff_angles_x = gc_x - imu_angles_avg_x;
%plot(1:16,imu_diff_angles_z,'g');

imu_angles_avg_y = [ 0.1696    0.2824    0.4460    0.3824   -1.1762   -0.1161    0.2092   -6.3082   -0.3279    3.1753    3.3891    1.2530    1.3985    3.2656    0.3374     0.5920]';
gc_y = [0, 1, 1, 0, 0, 1, 0, -6, 0, 2, 4, 3, 0, 4, -2, 0]';
%figure;
%plot(1:16,gc_y,'g', 1:16, imu_angles_avg_y,'b');
imu_diff_angles_y = gc_y - imu_angles_avg_y;
% 
% subplot(1,3,1)
% plot(1:16,imu_angles_avg_z,'b', 1:16, gc_z, 'g');
% title('z')
% subplot(1,3,2)
% plot(1:16,imu_angles_avg_x,'b', 1:16, gc_x, 'g');
% title('x')
% subplot(1,3,3)
% plot(1:16,imu_angles_avg_y,'b', 1:16, gc_y, 'g');
% title('y')


cam_angle_x = [
    0.6
    -3.5625
  -10.9375
   -8.8125
   -3.0000
   -2.0625
   -17.2500
    -14.6250
    -14.6250
  -14.6875
  -20.0625
  -20.0625
  -27.3750
  -27.3750
  -27.3750
     -14.0625];
gc_x = [0 -3 -13 -3 -5 -5 -18 -14 -14 -22 -22 -32 -30 -50 -71 -10]';

cam_angle_y = [-1 0.7 4.3 4.4 4.2 0.4 0.5 0.6 0.6 1.1 -0.1 -1 -0.5 -0.5 -0.5 -0.5]'
gc_y = [0, 1, 1, 0, 0, 1, 0, -6, 0, 2, 4, 3, 0, 4, -2, 0]';

cam_angle_z = [0.6 0.6 0.9 0.2 0.5 -0.1 -0.2 1.2 1.2 0.8 3.4 3.4 2.3 2.3 2.3 2.1]'
gc_z = [0, -2.5, 0, 0, 0, -2, 0, 0, -3.5, 0, 3, 6.5, -1, 0, -0.5, 0]';

% figure
% subplot(1,3,1)
% plot(1:16,cam_angle_z,'b', 1:16, gc_z, 'g');
% title('z')
% subplot(1,3,2)
% plot(1:16,cam_angle_x,'b', 1:16, gc_x, 'g');
% title('x')
% subplot(1,3,3)
% plot(1:16,cam_angle_y,'b', 1:16, gc_y, 'g');
% title('y')

% Measurement noise covariances
R1 = 0.002;
R2 = 0.002;

% Process noise covariance
Q = .005;

% State transition model
A = 1;

% Observation model
C1 = 1;
C2 = 1;

%camxhatx = kalman(cam_locations_scale(:,1)', A, C1, R1, Q);
%camxhaty = kalman(cam_locations_scale(:,2)', A, C1, R1, Q);

blecamxhatanglex1 = kalman([imu_angles_avg_x(1:13)'; cam_angle_x(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
blecamxhatanglex2 = kalman([blecamxhatanglex1(13), imu_angles_avg_x(14:15)'], A, C1, R1, Q);
blecamxhatanglex3 = kalman([[blecamxhatanglex2(3), imu_angles_avg_x(16)']; [blecamxhatanglex2(3), cam_angle_x(16)']],  A, [C1; C2], [R1 0; 0 R2], Q);

blecamxhatangley1 = kalman([imu_angles_avg_y(1:13)'; cam_angle_y(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
blecamxhatangley2 = kalman([blecamxhatangley1(13), imu_angles_avg_y(14:15)'], A, C1, R1, Q);
blecamxhatangley3 = kalman([[blecamxhatangley2(3), imu_angles_avg_y(16)']; [blecamxhatangley2(3), cam_angle_y(16)']],  A, [C1; C2], [R1 0; 0 R2], Q);

blecamxhatanglez1 = kalman([imu_angles_avg_z(1:13)'; cam_angle_z(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
blecamxhatanglez2 = kalman([blecamxhatanglez1(13), imu_angles_avg_z(14:15)'], A, C1, R1, Q);
blecamxhatanglez3 = kalman([[blecamxhatanglez2(3), imu_angles_avg_z(16)']; [blecamxhatanglez2(3), cam_angle_z(16)']],  A, [C1; C2], [R1 0; 0 R2], Q);

blecamxhatanglex = [blecamxhatanglex1 blecamxhatanglex2(2:3) blecamxhatanglex3(2)]
blecamxhatangley = [blecamxhatangley1 blecamxhatangley2(2:3) blecamxhatangley3(2)]
blecamxhatanglez = [blecamxhatanglez1 blecamxhatanglez2(2:3) blecamxhatanglez3(2)]

%blecamxhatanglex2 = kalman([imu_angles_avg_x'; cam_angle_x'], A, [C1; C2], [R1 0; 0 R2], Q);
%blecamxhatangley1 = kalman([imu_angles_avg_y(1:13)'; cam_angle_y(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
%blecamxhatangley2 = kalman([imu_angles_avg_x'; cam_angle_x'], A, [C1; C2], [R1 0; 0 R2], Q);
%blecamxhatanglez1 = kalman([imu_angles_avg_z(1:13)'; cam_angle_z(1:13)'], A, [C1; C2], [R1 0; 0 R2], Q);
%blecamxhatanglez2 = kalman([imu_angles_avg_x'; cam_angle_x'], A, [C1; C2], [R1 0; 0 R2], Q);


% figure
% subplot(1,3,1)617
% plot(1:16,blecamxhatanglez,'b', 1:16, gc_z, 'g');
% title('z')
% subplot(1,3,2)
% plot(1:16,blecamxhatanglex,'b', 1:16, gc_x, 'g');
% title('x')
% subplot(1,3,3)
% plot(1:16,blecamxhatangley,'b', 1:16, gc_y, 'g');
% title('y')


imuangle_heatmap = reshape(sqrt((imu_angles_avg_z -  gc_z).^2 + (imu_angles_avg_x -  gc_x).^2 + (imu_angles_avg_y -  gc_y).^2), 4,4);
camangle_heatmap = reshape(sqrt((cam_angle_z -  gc_z).^2 + (cam_angle_x -  gc_x).^2 + + (cam_angle_y -  gc_y).^2), 4,4);
imuangle_heatmap(1) = camangle_heatmap(1)
%camk_heatmap = reshape(sqrt((camxhatx' -  gc_locations_2d(:,1)).^2 + (camxhaty' -  gc_locations_2d(:,2)).^2), 4,4);
camimukangle_heatmap = reshape(sqrt((blecamxhatanglex' -  gc_x).^2 + (blecamxhatangley' -  gc_y).^2 + (blecamxhatanglez' -  gc_z).^2), 4,4);

figure
plot(1:16, sqrt((cam_angle_z -  gc_z).^2 + (cam_angle_x -  gc_x).^2 + + (cam_angle_y -  gc_y).^2), 'b',1:16, sqrt((imu_angles_avg_z -  gc_z).^2 + (imu_angles_avg_x -  gc_x).^2 + (imu_angles_avg_y -  gc_y).^2), 'r', 1:16, sqrt((imu_angles_avg_z -  gc_z).^2 + sqrt((blecamxhatanglex' -  gc_x).^2 + (blecamxhatangley' -  gc_y).^2 + (blecamxhatanglez' -  gc_z).^2)), 'g');


figure;
subplot(1,3,1)
imagesc(bledist_heatmap, [0,55])
colormap('hot')
title('BLE only')
axis off

subplot(1,3,2)
imagesc(camdist_heatmap, [0,55])
colormap('hot')
title('Camera Only')
axis off

% subplot(2,3,3)
% imagesc(camk_heatmap, [0,55])
% colormap('hot')
% title('cam k')
subplot(1,3,3)
imagesc(camblek_heatmap, [0,55])
colormap('hot')
title('BLE-Camera Sensor Fusion')
axis off


figure;
subplot(1,3,1)
imagesc(imuangle_heatmap, [0,25])
colormap('hot')
title('IMU only')
axis off

subplot(1,3,2)
imagesc(camangle_heatmap, [0,25])
colormap('hot')
title('Camera Only')
axis off

subplot(1,3,3)
imagesc(camimukangle_heatmap, [0,25])
colormap('hot')
title('IMU-Camera Sensor Fusion')
axis off


