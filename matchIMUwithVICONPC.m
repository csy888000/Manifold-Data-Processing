close all
clear

dbstop if error

addpath IMU_data
addpath Vicon_Data

%% Read data

load testDataAll.mat
load IMUDataAll.mat

%%
testStruct = cell(1,length(testDataAll));
for i = 1:length(testDataAll)
    testStruct{i} = cell2struct([struct2cell(testDataAll{i});struct2cell(IMUDataAll{i})],[fieldnames(testDataAll{i});fieldnames(IMUDataAll{i})]);
end

for i = 1:length(testDataAll)
    eval(['test',num2str(i), '= testStruct{i};']);
end



%% Match NUC data with VICON PC data
ratioRound = cell(1,length(testDataAll));
for i = 1:length(testDataAll)
    viconFromPC = testStruct{i}.testObj.viconData(:,109:end);
    viconFrame = 1:size(viconFromPC,1);

    VICON_Data_From_NUC = testStruct{i}.IMUobj.markerIMU;
    imuFrame = 1:size(VICON_Data_From_NUC,1);

    
    VICON_Data_From_NUC_round=roundn(VICON_Data_From_NUC,-3);
    viconFromPC_round=roundn(viconFromPC,-3);

    [dataDual, idPC, idNUC] = intersect(viconFromPC_round(:,1:3), VICON_Data_From_NUC_round(:,1:3), 'rows');
    ratio = polyfit(idPC,idNUC,1); 

    ratioRound{i} = round(ratio);
    
end

save ratio ratioRound

%% Plot motion data
figure(1)
clf
hold on
plot(imuFrame./100, VICON_Data_From_NUC(:,8),'b')
plot(((1:length(viconFrame))*ratio(1)+ratio(2))./100, viconFromPC(:,8),'r')
legend('Motion from NUC','Motion from VICON PC')

% figure(2)
% clf
% hold on
% plot(TimeSimulink, RowanIMU(:,21)*0.0125,'b')
% plot(((1:length(viconFrame))*ratio(1)+ratio(2))./50, thetaZShankLeft,'r')
% legend('IMU Angle','Vicon Angle')



