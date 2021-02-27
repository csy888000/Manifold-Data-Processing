function testObj = getVicondata(viconFilename)

viconPCData_FILL = importVICONdataFullbody(viconFilename);
viconFromPC = viconPCData_FILL(1:end-1,3:end);


%% Get time
timeVicon = viconPCData_FILL(1:end-1,1)./100;

%% VICON data preprocess
[viconFilledFromPC, ViconObj] = fillViconGapsFullbody(viconFromPC);

%% Get joint angles
[angleSubj, angleJointSubj] = getViconAngles(ViconObj);

testObj.viconData = viconFilledFromPC;
testObj.angleData = angleSubj;
testObj.angleJointData = angleJointSubj;
testObj.timeVicon = timeVicon;


end