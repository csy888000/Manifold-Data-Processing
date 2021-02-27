function IMUtest1 = getIMUdata(IMU_FILENAME)

% clear
% close all
% 
% dbstop if error
% 
% addpath Vicon_Data

%%
% LShoe RShoe LShank RShank LThigh RThigh LowerBack
% lin acc*0.002394*2
% gyro*0.07*2
% angle*0.0125

% IMU data
% Linear acc x y z
% Gyro x y z
% Angle z
% x forward y up


% test1.IMU_FILENAME = '210129_223336_test11.xls';
% IMU_data = readmatrix(test1.IMU_FILENAME);

IMU_data = readmatrix(IMU_FILENAME);

frameIMU = IMU_data(:,1);

periodIMU = IMU_data(:,3);
angleIMU = IMU_data(:,4:6);
gyroIMU = IMU_data(:,7:9)*180/pi;
linAccIMU = IMU_data(:,10:12)*9.81;
magnetIMU = IMU_data(:,13:15);
angleIMU(angleIMU(:,1)<0,1) = angleIMU(angleIMU(:,1)<0,1)+360;
angleIMU(angleIMU(:,3)>75,3) = angleIMU(angleIMU(:,3)>75,3)-360;
linAccIMU(linAccIMU>40) = NaN;
linAccIMU(linAccIMU<-40) = NaN;
linAccIMU = fillgaps(linAccIMU);

magnetIMU(magnetIMU>500) = NaN;
magnetIMU(magnetIMU<-500) = NaN;
magnetIMU = fillgaps(magnetIMU);

gyroIMU(gyroIMU>600) = NaN;
gyroIMU(gyroIMU<-600) = NaN;
gyroIMU = fillgaps(gyroIMU);

% sevenIMUraw = IMU_data(:,16:64);
% sevenIMU = zeros(size(sevenIMUraw,1),size(sevenIMUraw,2));
% sevenIMU = sevenIMUraw;
% for i = 1:7
%     sevenIMU(:,i*7-6) = sevenIMUraw(:,i*7-6)/0.0125/2*0.002394*2;
%     sevenIMU(:,i*7-3) = sevenIMUraw(:,i*7-3)/0.002394/2*0.07*2;
%     sevenIMU(:,i*7) = sevenIMUraw(:,i*7)/0.07/2*0.0125;
% end


% sevenIMU(sevenIMU==0) = NaN; 

markerIMU = IMU_data(:,16:end);

IMUtest1.angleIMU = angleIMU;
IMUtest1.gyroIMU = gyroIMU;
IMUtest1.linAccIMU = linAccIMU;
IMUtest1.magnetIMU = magnetIMU;
IMUtest1.periodIMU = periodIMU;
% IMUtest1.sevenIMU = sevenIMU;
IMUtest1.markerIMU = markerIMU;



%%
figure(100)
clf
title('single IMU')
for i = 1:3
    subplot(9,1,i)
    hold on
    plot(linAccIMU(:,i),'b')
end
for i = 1:3
    subplot(9,1,i+3)
    hold on
    plot(gyroIMU(:,i),'b')
end
for i = 1:3
    subplot(9,1,i+6)
    hold on
    plot(angleIMU(:,i),'b')
end






