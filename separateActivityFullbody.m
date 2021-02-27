clear
close all

dbstop if error

addpath IMU_data
addpath Vicon_Data

%%

load testDataAll.mat
load IMUDataAll.mat
load ratio.mat

testStruct = cell(1,length(testDataAll));
for i = 1:length(testDataAll)
    testStruct{i} = cell2struct([struct2cell(testDataAll{i});struct2cell(IMUDataAll{i})],[fieldnames(testDataAll{i});fieldnames(IMUDataAll{i})]);
end

for i = 1:length(testStruct)
    eval(['test',num2str(i), '= testStruct{i};']);
end



%%

jointAngle = cell(1,length(testStruct));
for i = 1:length(testStruct)
%     range = testStruct{i}.range;

    viconMarker = testStruct{i}.testObj.viconData(:,1:3);
    
    hipLeftSeg = testStruct{i}.testObj.angleJointData{1}(:);
    hipRightSeg = testStruct{i}.testObj.angleJointData{2}(:);
    kneeLeftSeg = testStruct{i}.testObj.angleJointData{3}(:);
    kneeRightSeg = testStruct{i}.testObj.angleJointData{4}(:);
    ankleLeftSeg = testStruct{i}.testObj.angleJointData{5}(:);
    ankleRightSeg = testStruct{i}.testObj.angleJointData{6}(:);
    
    
    shoulderLeftSeg = testStruct{i}.testObj.angleJointData{8}(:);
    shoulderRightSeg = testStruct{i}.testObj.angleJointData{9}(:);
    elbowLeftSeg = testStruct{i}.testObj.angleJointData{10}(:);
    elbowRightSeg = testStruct{i}.testObj.angleJointData{11}(:);
    handLeftSeg = testStruct{i}.testObj.angleJointData{12}(:);
    handRightSeg = testStruct{i}.testObj.angleJointData{13}(:);

    jointAngle{i} = [hipLeftSeg, hipRightSeg, kneeLeftSeg, kneeRightSeg, ankleLeftSeg, ankleRightSeg,...
        shoulderLeftSeg, shoulderRightSeg, elbowLeftSeg, elbowRightSeg, handLeftSeg, handRightSeg];
    markerOne{i} = viconMarker;
end

imuData = cell(1,length(testStruct));
for i = 1:length(testStruct)
    
    angleIMU = testStruct{i}.IMUobj.angleIMU;
    gyroIMU = testStruct{i}.IMUobj.gyroIMU;
    linAccIMU = testStruct{i}.IMUobj.linAccIMU;
    magnetIMU = testStruct{i}.IMUobj.magnetIMU;
    markerIMU = testStruct{i}.IMUobj.markerIMU(:,1:3);
    imuData{i} = [angleIMU gyroIMU linAccIMU magnetIMU markerIMU];
end


bodyAngle = cell(1,length(testStruct));
for i = 1:length(testStruct)
    
    thetaBody = testStruct{i}.testObj.angleJointData{7}(1:end);

    bodyAngle{i} = [thetaBody];
end




for num = 1:length(testStruct)
    %%

    ratio = ratioRound{num};

    frameViconRaw = 1:length(bodyAngle{num});
    timeViconAll = frameViconRaw./100;

    frameIMURaw = 1:length(imuData{num});

    frameViconAlign = frameViconRaw*ratio(1)+ratio(2);
    startEndViconAlign = [round(frameViconAlign(1)) round(frameViconAlign(end))];

    frameEND = min(frameIMURaw(end), startEndViconAlign(end));

    frameIMU = startEndViconAlign(1):frameEND;



    frameVicon = 1:((frameEND-ratio(2))/ratio(1));

    timeIMU = frameIMU./100;
    timeVicon = frameVicon./100;

    testStruct{num}.frameIMU = frameIMU;
    testStruct{num}.frameVicon = frameVicon;
    
    
    
    %%
    figure(70)
    clf
    hold on
    plot(timeViconAll,bodyAngle{num});
    % plot(timeIMUAll,imuAngle{num}(frameIMUAll,2))
    ylabel('Knee Joint Angle')



    theta = bodyAngle{num};
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
    stateIMU = state(frameVicon);
    vMarker = markerOne{num}(frameVicon, 1:3);

    frameWalkup = unique(round(timeViconAll(stateIMU == 0)*100));
    frameWalkDown = unique(round(timeViconAll(stateIMU == 180)*100));
%     if rem(num,2)
%         frameFlatGround = unique(round(timeViconAll(stateIMU == 90)*100));
%         frameWalkAcross = unique(round(timeViconAll(stateIMU == -90)*100));     
%     else
%         frameFlatGround = unique(round(timeViconAll(stateIMU == -90)*100));
%         frameWalkAcross = unique(round(timeViconAll(stateIMU == 90)*100));   
%     end

    frameFlatGround = unique(round(timeViconAll(stateIMU == 90)*100));
    frameWalkAcross = unique(round(timeViconAll(stateIMU == -90)*100));   

    stateWalk = zeros(1,length(frameFlatGround)+length(frameWalkup)+length(frameWalkAcross)+length(frameWalkDown));
    stateWalk(frameFlatGround) = 90;
    stateWalk(frameWalkup) = 0;
    stateWalk(frameWalkAcross) = -90;
    stateWalk(frameWalkDown) = 180;
    

    labelWalk = zeros(length(stateWalk),1);
    labelWalk(stateWalk == 90) = 1;
    labelWalk(stateWalk == 0) = 2;
    labelWalk(stateWalk == -90) = 3;
    labelWalk(stateWalk == 180) = 4;
    
    rangeVicon = testStruct{num}.rangeVicon;
    labelWalk(rangeVicon(length(rangeVicon)):end) = 3;
    
    testStruct{num}.stateWalk = stateWalk;
    testStruct{num}.labelWalk = labelWalk';
    
    
    figure(71)
    clf
    hold on
    plot(stateWalk)
    plot(imuData{num}(frameIMU,1))



    figure(72)
    clf
    hold on
    % plot(frameWalk)
    plot(imuData{num}(frameIMU,13))
    plot(vMarker(:,1))

    %%
    imuDataSeg = imuData{num}(frameIMU,:);
    flatGroundIMU = imuDataSeg(frameFlatGround,:);
    walkUpIMU = imuDataSeg(frameWalkup,:);
    walkAcrossIMU = imuDataSeg(frameWalkAcross,:);
    walkDownIMU = imuDataSeg(frameWalkDown,:);

    figure(73)
    clf
    subplot(4,1,1)
    hold on
    plot(flatGroundIMU(:,13))
    plot(flatGroundIMU(:,14))
    subplot(4,1,2)
    hold on
    plot(walkUpIMU(:,13))
    plot(walkUpIMU(:,14))
    subplot(4,1,3)
    hold on
    plot(walkAcrossIMU(:,13))
    plot(walkAcrossIMU(:,14))
    subplot(4,1,4)
    hold on
    plot(walkDownIMU(:,13))
    plot(walkDownIMU(:,14))

    figure(74)
    clf
    subplot(4,1,1)
    hold on
    plot(flatGroundIMU(:,1))
    plot(flatGroundIMU(:,4))
    subplot(4,1,2)
    hold on
    plot(walkUpIMU(:,1))
    plot(walkUpIMU(:,4))
    subplot(4,1,3)
    hold on
    plot(walkAcrossIMU(:,1))
    plot(walkAcrossIMU(:,4))
    subplot(4,1,4)
    hold on
    plot(walkDownIMU(:,1))
    plot(walkDownIMU(:,4))
    
    
    
    %%
    walkUpIMUsave{num} = [walkUpIMU(:,1:12) ones(length(walkUpIMU),1)*2];
    walkDownIMUsave{num} = [walkDownIMU(:,1:12) ones(length(walkDownIMU),1)*4];
    
    flatGroundIMUsave{num} = [flatGroundIMU(:,1:12) ones(length(flatGroundIMU),1)*3];
    walkAcrossIMUsave{num} = [walkAcrossIMU(:,1:12) ones(length(walkAcrossIMU),1)*1];

    
    IMUdataSave{num} = [flatGroundIMUsave{num};walkUpIMUsave{num};walkAcrossIMUsave{num};walkDownIMUsave{num}];
%     if rem(num,2)
%         walkAcrossIMUsave{num} = [walkAcrossIMU(:,1:12) ones(length(walkAcrossIMU),1)*3];
%     else
%         walkAcrossIMUsave{num} = [walkAcrossIMU(:,1:12) ones(length(walkAcrossIMU),1)*4];
%     end
    IMUdataSave_test{num} = [imuData{num}(frameIMU,1:12) labelWalk jointAngle{num}(frameVicon,1:end)];
    labelAll{num} = labelWalk;
    
    
end


disp('Saving...')
writematrix(IMUdataSave_test{1},'Subj04_IMU_walk_5deg_test1.csv')
writematrix(IMUdataSave_test{2},'Subj04_IMU_walk_5deg_test2.csv')
writematrix(IMUdataSave_test{3},'Subj04_IMU_walk_10deg_test3.csv')
writematrix(IMUdataSave_test{4},'Subj04_IMU_walk_10deg_test4.csv')
writematrix(IMUdataSave_test{5},'Subj04_IMU_walk_15deg_test5.csv')
writematrix(IMUdataSave_test{6},'Subj04_IMU_walk_15deg_test6.csv')
writematrix(IMUdataSave_test{7},'Subj04_IMU_walk_20deg_test7.csv')
writematrix(IMUdataSave_test{8},'Subj04_IMU_walk_20deg_test8.csv')




save('testStruct.mat','testStruct')


disp('Save done')










