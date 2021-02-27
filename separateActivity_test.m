clear
close all

addpath IMUdata
load Test5_walk_square_slope20_nostop
% load Test3_walk_speed40_slope0.mat
% load testDataAll.mat
% load ratio

%%
% lin acc*0.002394*2
% gyro*0.07*2
% angle*0.0125

% IMU data
% Linear acc x y z
% Gyro x y z
% Angle z
% x forward y up

GAIT = KneelingData.GAIT.Data;

walkIMUraw = double(KneelingData.IMU_s.Data);

walkIMU = zeros(size(walkIMUraw,1),size(walkIMUraw,2));

for i = 1:7
    walkIMU(:,i*7-6:i*7-4) = walkIMUraw(:,i*7-6:i*7-4)*0.002394*2;
    walkIMU(:,i*7-3:i*7-1) = walkIMUraw(:,i*7-3:i*7-1)*0.07*2;
    walkIMU(:,i*7) = walkIMUraw(:,i*7)*0.0125;
end


figure(51)
clf
for i = 1:7
    subplot(7,1,i)
    plot(walkIMU(:,42+i))

end

backIMU = walkIMU(:,end-6:end);

walkIMU(~any(walkIMU,2), :) = [];

writematrix(walkIMU,'IMU_test.csv') 

% writematrix(flatGroundIMU,'flatGroundIMU_test.csv') 
% writematrix(walkUpIMU,'walkUpIMU_test.csv') 
% writematrix(walkAcrossIMU,'walkAcrossIMU_test.csv')
% writematrix(walkDownIMU,'walkDownIMU_test.csv') 

