clear
close all

dbstop if error

addpath IMU_data
addpath Vicon_Data


load testStruct.mat
%%

jointAngle = cell(1,length(testStruct));
imuData = cell(1,length(testStruct));
labelWalk = cell(1,length(testStruct));
stateWalk = cell(1,length(testStruct));
for i = 1:length(testStruct)
    rangeVicon = testStruct{i}.rangeVicon;
    frameVicon = testStruct{i}.frameVicon;
    frameIMU = testStruct{i}.frameIMU;

    hipLeftSeg = testStruct{i}.testObj.angleJointData{1}(frameVicon);
    hipRightSeg = testStruct{i}.testObj.angleJointData{2}(frameVicon);
    kneeLeftSeg = testStruct{i}.testObj.angleJointData{3}(frameVicon);
    kneeRightSeg = testStruct{i}.testObj.angleJointData{4}(frameVicon);
    ankleLeftSeg = testStruct{i}.testObj.angleJointData{5}(frameVicon);
    ankleRightSeg = testStruct{i}.testObj.angleJointData{6}(frameVicon);

    shoulderLeftSeg = testStruct{i}.testObj.angleJointData{8}(frameVicon);
    shoulderRightSeg = testStruct{i}.testObj.angleJointData{9}(frameVicon);
    elbowLeftSeg = testStruct{i}.testObj.angleJointData{10}(frameVicon);
    elbowRightSeg = testStruct{i}.testObj.angleJointData{11}(frameVicon);
    handLeftSeg = testStruct{i}.testObj.angleJointData{12}(frameVicon);
    handRightSeg = testStruct{i}.testObj.angleJointData{13}(frameVicon);
    
    
%     jointAngle{i} = [hipLeftSeg, hipRightSeg, kneeLeftSeg, kneeRightSeg, ankleLeftSeg, ankleRightSeg];
    jointAngle{i} = [hipLeftSeg, hipRightSeg, kneeLeftSeg, kneeRightSeg, ankleLeftSeg, ankleRightSeg,...
        shoulderLeftSeg, shoulderRightSeg, elbowLeftSeg, elbowRightSeg, handLeftSeg, handRightSeg];
    
    angleIMU = testStruct{i}.IMUobj.angleIMU(frameIMU,:);
    gyroIMU = testStruct{i}.IMUobj.gyroIMU(frameIMU,:);
    linAccIMU = testStruct{i}.IMUobj.linAccIMU(frameIMU,:);
    magnetIMU = testStruct{i}.IMUobj.magnetIMU(frameIMU,:);
    imuData{i} = [angleIMU gyroIMU linAccIMU magnetIMU];
    
    labelWalk{i} = testStruct{i}.labelWalk;
    stateWalk{i} = testStruct{i}.stateWalk;
    
    
    heelLeftSeg{i} = testStruct{i}.testObj.viconData(frameVicon,94:96);
    heelRightSeg{i} = testStruct{i}.testObj.viconData(frameVicon,112:114);
% end


%%
    num = i;
% for num = 7%1:length(testStruct)
%     num = 7;
    % rangeStep = 1:23450;
    rangeStep = intersect(rangeVicon, frameVicon);

    heelLeft = heelLeftSeg{num};
    heelRight = heelRightSeg{num};

    heelLeftVel = diff(heelLeft);
    heelRightVel = diff(heelRight);

    % filter data, get pos vel acc
    [B,A] = butter(4,20/(100/2),'low');
    heelLeft_filt = filtfilt(B,A,heelLeft);
    heelRight_filt = filtfilt(B,A,heelRight);

    heelLeftVel_filt = diff(heelLeft_filt);
    heelRightVel_filt = diff(heelRight_filt);

    heelLeftAcc_filt = diff(heelLeftVel_filt);
    heelRightAcc_filt = diff(heelRightVel_filt);


    % get heel peaks
    [heelLeftPeaks, heelLeftLocs] = findpeaks(heelLeft(rangeStep,3),'MinPeakProminence',30);
    [heelRightPeaks, heelRightLocs] = findpeaks(heelRight(rangeStep,3),'MinPeakProminence',30);

    zeroCrossingInx = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    heelLeftZerosAll = zeroCrossingInx(heelLeftVel_filt(rangeStep,3));
    heelRightZerosAll = zeroCrossingInx(heelRightVel_filt(rangeStep,3));


    % get heel strike
    heelLeftZerosInx = zeros(length(heelLeftLocs)-1,1);
    heelLeftZeros = zeros(length(heelLeftLocs)-1,1);
    for j = 1:length(heelLeftZerosInx)
        heelLeftZerosInx(j) = find(heelLeftZerosAll > heelLeftLocs(j), 1, 'first');
        heelLeftZeros(j) = heelLeftZerosAll(heelLeftZerosInx(j))+1;
    end

    heelRightZerosInx = zeros(length(heelRightLocs)-1,1);
    heelRightZeros = zeros(length(heelRightLocs)-1,1);
    for j = 1:length(heelRightZerosInx)
        heelRightZerosInx(j) = find(heelRightZerosAll > heelRightLocs(j), 1, 'first');
        heelRightZeros(j) = heelRightZerosAll(heelRightZerosInx(j))+1;
    end


    % plot data
    figure(201)
    clf
    subplot(3,1,1)
    hold on
    plot(heelLeftAcc_filt(rangeStep,3))
    plot(heelRightAcc_filt(rangeStep,3))
    subplot(3,1,2)
    hold on
    plot(heelLeftVel_filt(rangeStep,3))
    plot(heelRightVel_filt(rangeStep,3))
    subplot(3,1,3)
    hold on
    p1 = plot(heelLeft(rangeStep,3),'b');
    p2 = plot(heelRight(rangeStep,3),'r');
    plot(heelLeftLocs, heelLeftPeaks,'ko')
    plot(heelRightLocs, heelRightPeaks,'cx')
    plot(heelLeftZeros, heelLeft(rangeStep(1)+heelLeftZeros-1,3),'md')
    plot(heelRightZeros, heelRight(rangeStep(1)+heelRightZeros-1,3),'gd')

    legend([p1, p2], 'Left','Right')


    %%
    labelSeg = labelWalk{num}(rangeStep);

    figure(202)
    clf
    hold on
    plot(labelWalk{num}(rangeStep))
    plot(heelLeftZeros, labelSeg(heelLeftZeros),'md')
    plot(heelRightZeros, labelSeg(heelRightZeros),'gd')


    heelStrike = sort([heelLeftZeros;heelRightZeros]);


    labelChange = abs(diff(labelSeg));
    labelChange(labelChange>1) = 1;
    labelMark = find(labelChange == 1);

    plot(labelMark, labelSeg(labelMark),'ko')

    heelStrikeSave = [];
    heelStrikeValueSave = cell(1,1);
    for j = 1:length(labelMark)
        if j == 1
            segmentLabel{j} = labelSeg(1:labelMark(j));
            segmentLimit = 1:labelMark(j);
            heelStrikeInx = find(heelStrike < labelMark(j));
        else
            segmentLabel{j} = labelSeg(labelMark(j-1)+1:labelMark(j));
            segmentLimit = labelMark(j-1)+1:labelMark(j);
            heelStrikeInx = find(heelStrike < labelMark(j) & heelStrike > labelMark(j-1)+1);
        end

        heelStrikeValue = heelStrike(heelStrikeInx);
        heelStrikeDiff = diff(heelStrikeValue);
        heelStrikeDiffInx = find(heelStrikeDiff < 85);
        heelStrikeDiffInxDiff = diff(heelStrikeDiffInx);
        heelRemove = find(heelStrikeDiffInxDiff>1);
        heelStrikeDiffInx(heelRemove) = [];

        if isempty(heelStrikeInx)
            heelStrikeValueSave{j} = [];
        else
            heelStrikeValueSave{j} = heelStrikeValue(heelStrikeDiffInx(1):heelStrikeDiffInx(end)+1);
        end
        heelStrikeSave = [heelStrikeSave;heelStrikeValueSave{j}];
    end

    % plot(heelStrikeSave, labelSeg(heelStrikeSave),'rx')

    heelLeftSave = intersect(heelStrikeSave, heelLeftZeros);
    heelRightSave = intersect(heelStrikeSave, heelRightZeros);

    plot(heelLeftSave, labelSeg(heelLeftSave),'rx')
    plot(heelRightSave, labelSeg(heelRightSave),'bx')
    legend('','Left','Right','','Left','Right')


    %%
    jointTh = jointAngle{num};
    imuTrial = imuData{num};


    walkSeg = cell(1,length(heelStrikeValueSave));
    imuSeg = cell(1,length(heelStrikeValueSave));
    label_segment = cell(1,length(heelStrikeValueSave));    
        
    for j = 1:length(heelStrikeValueSave)
        heelLeftValueSave{j} = intersect(heelStrikeValueSave{j}, heelLeftSave);
        heelRightValueSave{j} = intersect(heelStrikeValueSave{j}, heelRightSave);

        walkHeel = heelStrikeValueSave{j};
        if isempty(walkHeel)
            walkSeg{j} = [];
            imuSeg{j} = [];
            label_segment{j} = [];
        else
            walkSeg{j} = jointTh(walkHeel(1):walkHeel(end),:);
            imuSeg{j} = imuTrial(walkHeel(1):walkHeel(end),:);
            label_segment{j} = labelWalk{num}(walkHeel(1):walkHeel(end));
        end

    end
    
    
    walkSeg = walkSeg(~cellfun(@isempty, walkSeg));
    imuSeg = imuSeg(~cellfun(@isempty, imuSeg));
    label_segment = label_segment(~cellfun(@isempty, label_segment));
    
    heelLeftValueSave = heelLeftValueSave(~cellfun(@isempty, heelLeftValueSave));
    heelRightValueSave = heelRightValueSave(~cellfun(@isempty, heelRightValueSave));

    figure(203)
    clf
    hold on
    plot(walkSeg{5}(:,3))


    walkAll{num}.jointTh = jointTh;
    walkAll{num}.imuTrial = imuTrial;
    walkAll{num}.heelLeftValueSave = heelLeftValueSave;
    walkAll{num}.heelRightValueSave = heelRightValueSave;
    walkAll{num}.walkSeg = walkSeg;
    walkAll{num}.imuSeg = imuSeg;
    walkAll{num}.label_segment = label_segment;

end


%%

disp('Saving...')
save('MatSave/walkAll.mat','walkAll')
disp('Save done.')




%%




% figure(301)
% subplot(3,1,1)
% plot(labelWalk{num}(rangeVicon))
% subplot(3,1,2)
% plot(jointAngle{num}(rangeVicon,6))
% subplot(3,1,3)
% plot(imuData{num}(rangeVicon,3))




%%


















