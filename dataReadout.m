clear
close all

dbstop if error

addpath IMU_data
addpath Vicon_Data


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% treadmill speed 1.25m/s
test1.VICON_FILENAME = 'Subj01_Test11.csv';
test1.testObj = getVicondata(test1.VICON_FILENAME);
test1.testNum = 1;
range = 5178:1:5451;
test1.range = range;
disp('finish: 1')


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


test1.IMU_FILENAME = '210129_223336_test11.xls';
IMU_data = readmatrix(test1.IMU_FILENAME);



frameIMU = IMU_data(:,1);

periodIMU = IMU_data(:,3);
angleIMU = IMU_data(:,4:6);
gyroIMU = IMU_data(:,7:9)*180/pi;
linAccIMU = IMU_data(:,10:12)*9.81;
magnetIMU = IMU_data(:,13:15);


sevenIMUraw = IMU_data(:,16:64);
% sevenIMU = zeros(size(sevenIMUraw,1),size(sevenIMUraw,2));
sevenIMU = sevenIMUraw;
for i = 1:7
    sevenIMU(:,i*7-6) = sevenIMUraw(:,i*7-6)/0.0125/2*0.002394*2;
    sevenIMU(:,i*7-3) = sevenIMUraw(:,i*7-3)/0.002394/2*0.07*2;
    sevenIMU(:,i*7) = sevenIMUraw(:,i*7)/0.07/2*0.0125;
%     sevenIMU(:,i*7-6:i*7-4) = sevenIMUraw(:,i*7-6:i*7-4)*0.002394*2;
%     sevenIMU(:,i*7-3:i*7-1) = sevenIMUraw(:,i*7-3:i*7-1)*0.07*2;
%     sevenIMU(:,i*7) = sevenIMUraw(:,i*7)*0.0125;
end


sevenIMU(sevenIMU==0) = NaN;


markerIMU = IMU_data(:,65:112);


LShankAcc = sevenIMUraw(:,15:17);
LShankGyro = sevenIMUraw(:,18:20);
LShankAngle = sevenIMUraw(:,21);

LShankAcc(LShankAcc==0) = NaN;
LShankAcc(LShankAcc>50) = NaN;
LShankAcc(LShankAcc<-50) = NaN;

LShankGyro(LShankGyro==0) = NaN;
LShankGyro(LShankGyro>400) = NaN;
LShankGyro(LShankGyro<-400) = NaN;

LShankAngle(LShankAngle==0) = NaN;
LShankAngle(LShankAngle>200) = NaN;
LShankAngle(LShankAngle<-200) = NaN;

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


% figure(2)
% clf
% title('seven IMU')
% for i = 1:7
%     subplot(7,1,i)
%     plot(sevenIMU(:,14+i))
% end
% 
% figure(3)
% clf
% title('seven IMU')
% for i = 1:3
%     subplot(7,1,i)
%     hold on
%     plot(LShankAcc(:,i),'b')
% end
% for i = 1:3
%     subplot(7,1,i+3)
%     hold on
%     plot(LShankGyro(:,i),'b')
% end
% for i = 1
%     subplot(7,1,i+6)
%     hold on
%     plot(LShankAngle(:,i),'b')
% end
% 
% 
% %%
% figure(11)
% clf 
% hold on
% 
% subplot(3,1,1)
% hold on
% plot(1:length(linAccIMU),linAccIMU(:,1),'b')
% plot((1:length(linAccIMU))+0,-LShankAcc(:,2),'r')
% 
% subplot(3,1,2)
% hold on
% plot(1:length(linAccIMU),linAccIMU(:,2),'b')
% plot((1:length(linAccIMU))+50,LShankAcc(:,1),'r')
% 
% subplot(3,1,3)
% hold on
% plot(1:length(linAccIMU),linAccIMU(:,3),'b')
% plot((1:length(linAccIMU))+50,LShankAcc(:,3),'r')
% 
% 
% 
% figure(12)
% % clf 
% hold on
% 
% subplot(3,1,1)
% hold on
% plot(1:length(linAccIMU),gyroIMU(:,1),'b')
% plot((1:length(linAccIMU))+0,-LShankGyro(:,2),'r')
% 
% subplot(3,1,2)
% hold on
% plot(1:length(linAccIMU),gyroIMU(:,2),'b')
% plot((1:length(linAccIMU))-60,LShankGyro(:,1),'r')
% 
% subplot(3,1,3)
% hold on
% plot(1:length(linAccIMU),gyroIMU(:,3),'b')
% plot((1:length(linAccIMU))-70,LShankGyro(:,3),'r')







