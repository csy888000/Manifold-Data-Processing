clear
close all

addpath IMU_data
load Test4_walk_square_slope20_stopatcorner
% load Test3_walk_speed40_slope0.mat
load testDataAll.mat
load ratio

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



%%

jointAngle = cell(1,length(testDataAll));
for i = 1:length(testDataAll)
    range = testDataAll{i}.range;
    
    hipLeftSeg = testDataAll{i}.testObj.angleData{7}(range);
    hipRightSeg = testDataAll{i}.testObj.angleData{8}(range);
    kneeLeftSeg = testDataAll{i}.testObj.angleData{9}(range);
    kneeRightSeg = testDataAll{i}.testObj.angleData{10}(range);
    ankleLeftSeg = testDataAll{i}.testObj.angleData{11}(range);
    ankleRightSeg = testDataAll{i}.testObj.angleData{12}(range);

    jointAngle{i} = [hipLeftSeg, hipRightSeg, kneeLeftSeg, kneeRightSeg, ankleLeftSeg, ankleRightSeg];
end








%%
rangeVicon = testDataAll{2}.range;
timeVicon = rangeVicon./100;
frameIMUraw = rangeVicon*ratio(1)+ratio(2);
startEndIMU = [round(frameIMUraw(1)) round(frameIMUraw(end))];
frameIMU = startEndIMU(1):startEndIMU(end);
timeIMU = (frameIMU-ratio(2))./ratio(1)./100;

figure(60)
subplot(3,1,1)
hold on
plot(timeVicon,jointAngle{2}(:,4));
plot(timeIMU,GAIT(frameIMU,2)*-60)
plot(timeIMU,GAIT(frameIMU,1)*-60)
ylabel('Knee Joint Angle')
subplot(3,1,2)
plot(timeIMU,walkIMU(frameIMU,2*7-2))
ylabel('Right Shank Linear Acc x')
subplot(3,1,3)
plot(timeIMU,walkIMU(frameIMU,7*7-2))
ylabel('Back Linear Acc x')


%%
imuNum = 2;
figure(61)
clf
subplot(2,1,1)
hold on
plot(timeIMU,walkIMU(frameIMU,imuNum*7-6))

plot(timeIMU,walkIMU(frameIMU,imuNum*7-5))

plot(timeIMU,walkIMU(frameIMU,imuNum*7-4))
legend('x','y','z')
ylabel('lin acc')

subplot(2,1,2)
hold on
plot(timeIMU,walkIMU(frameIMU,imuNum*7-3))

plot(timeIMU,walkIMU(frameIMU,imuNum*7-2))

plot(timeIMU,walkIMU(frameIMU,imuNum*7-1))
ylabel('gyro')


%%
figure(62)
clf
subplot(3,1,1)
hold on
plot(timeIMU,walkIMU(frameIMU,2*7-3)*pi/180)
subplot(3,1,2)
plot(timeIMU,walkIMU(frameIMU,2*7-2)*pi/180)
subplot(3,1,3)
plot(timeIMU,walkIMU(frameIMU,2*7-1)*pi/180)




%%
bodyAngle = cell(1,length(testDataAll));
for i = 1:length(testDataAll)
    
    thetaBody = testDataAll{i}.testObj.angleData{13}(1:end);

    bodyAngle{i} = [thetaBody];
end

rangeViconAll = 1:length(bodyAngle{4});
timeViconAll = rangeViconAll./100;
frameIMUrawAll = rangeViconAll*ratio(1)+ratio(2);
startEndIMUAll = [round(frameIMUrawAll(1)) round(frameIMUrawAll(end))];
frameIMUAll = startEndIMUAll(1):startEndIMUAll(end);

timeIMUAll = (frameIMUAll-ratio(2))./ratio(1)./100;


%%
figure(70)
clf
hold on
plot(timeViconAll,bodyAngle{4});
% plot(timeIMUAll,GAIT(frameIMUAll,2)*-60)
% plot(timeIMUAll,GAIT(frameIMUAll,1)*-60)
plot(timeIMUAll,backIMU(frameIMUAll,2))
ylabel('Knee Joint Angle')



theta = bodyAngle{4};
state = ones(length(theta),1)*1000;
flatGround = zeros(length(theta),1);
walkUp = zeros(length(theta),1);
walkAcross = zeros(length(theta),1);
walkDown = zeros(length(theta),1);
for i = 1:length(theta)
    th = theta(i);

    if th > 45 && th < 135
        disp('flat ground')
        state(i) = 90;
    elseif th > -45 && th < 45
        disp('walk up')
        state(i) = 0;
    elseif th > -135 && th < -45
        disp('walk across')
        state(i) = -90;
    elseif th > 135 || th < -135
        disp('walk down')
        state(i) = 180;
    end
end

plot(timeViconAll, state)

%%
frameFlatGround = unique(round(timeViconAll(state == 90)*100*ratio(1)+ratio(2)));
frameWalkup = unique(round(timeViconAll(state == 0)*100*ratio(1)+ratio(2)));
frameWalkAcross = unique(round(timeViconAll(state == -90)*100*ratio(1)+ratio(2)));
frameWalkDown = unique(round(timeViconAll(state == 180)*100*ratio(1)+ratio(2)));


frameWalk = zeros(1,length(frameFlatGround)+length(frameWalkup)+length(frameWalkAcross)+length(frameWalkDown));
frameWalk(frameFlatGround) = 90;
frameWalk(frameWalkup) = 0;
frameWalk(frameWalkAcross) = -90;
frameWalk(frameWalkDown) = 180;

figure(72)
plot(frameWalk)


flatGroundIMU = walkIMU(frameFlatGround,:);
walkUpIMU = walkIMU(frameWalkup,:);
walkAcrossIMU = walkIMU(frameWalkAcross,:);
walkDownIMU = walkIMU(frameWalkDown,:);

writematrix(flatGroundIMU,'flatGroundIMU.csv') 
writematrix(walkUpIMU,'walkUpIMU.csv') 
writematrix(walkAcrossIMU,'walkAcrossIMU.csv')
writematrix(walkDownIMU,'walkDownIMU.csv') 

