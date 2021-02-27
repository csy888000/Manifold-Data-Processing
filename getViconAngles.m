function [angleSubj, angleJointSubj] = getViconAngles(ViconObj)

LFHD = ViconObj.LFHD;
RFHD = ViconObj.RFHD;
LBHD = ViconObj.LBHD;
RBHD = ViconObj.RBHD;

C7 = ViconObj.C7;
T10 = ViconObj.T10;
CLAV = ViconObj.CLAV;
STRN = ViconObj.STRN;
RBAK = ViconObj.RBAK;

LSHO = ViconObj.LSHO;
LUPA = ViconObj.LUPA;
LELB = ViconObj.LELB;
LFRM = ViconObj.LFRM;
LWRA = ViconObj.LWRA;
LWRB = ViconObj.LWRB;
LFIN = ViconObj.LFIN;

RSHO = ViconObj.RSHO;
RUPA = ViconObj.RUPA;
RELB = ViconObj.RELB;
RFRM = ViconObj.RFRM;
RWRA = ViconObj.RWRA;
RWRB = ViconObj.RWRB;
RFIN = ViconObj.RFIN;


LASI = ViconObj.LASI;
RASI = ViconObj.RASI;
LPSI = ViconObj.LPSI;
RPSI = ViconObj.RPSI;

LTHI = ViconObj.LTHI;
LKNE = ViconObj.LKNE;
LTIB = ViconObj.LTIB;
LANK = ViconObj.LANK;
LHEE = ViconObj.LHEE;
LTOE = ViconObj.LTOE;

RTHI = ViconObj.RTHI;
RKNE = ViconObj.RKNE;
RTIB = ViconObj.RTIB;
RANK = ViconObj.RANK;
RHEE = ViconObj.RHEE;
RTOE = ViconObj.RTOE;

thetaBody = NaN(length(LASI),1);
thetaZShoeLeft = NaN(length(LASI),1);
thetaZShoeRight = NaN(length(LASI),1);
thetaZShankLeft = NaN(length(LASI),1);
thetaZShankRight = NaN(length(LASI),1);
thetaZThighLeft = NaN(length(LASI),1);
thetaZThighRight = NaN(length(LASI),1);

thetaZHandLeft = NaN(length(LASI),1);
thetaZHandRight = NaN(length(LASI),1);
thetaZLowerArmLeft = NaN(length(LASI),1);
thetaZLowerArmRight = NaN(length(LASI),1);
thetaZUpperArmLeft = NaN(length(LASI),1);
thetaZUpperArmRight = NaN(length(LASI),1);

showFig = false;

if showFig
    figure(1)
    hold on
    axis equal
end

for i = 1:1:length(LASI) %8128:1:8259 
    frontMid = 0.5*(LASI(i,:)+RASI(i,:));
    rearMid = 0.5*(LPSI(i,:)+RPSI(i,:));
    
    if showFig
        LL1 = line([LASI(i,1) LPSI(i,1)],[LASI(i,2) LPSI(i,2)]);
        LL2 = line([LPSI(i,1) RPSI(i,1)],[LPSI(i,2) RPSI(i,2)]);
        LL3 = line([RPSI(i,1) RASI(i,1)],[RPSI(i,2) RASI(i,2)]);
        LL4 = line([RASI(i,1) LASI(i,1)],[RASI(i,2) LASI(i,2)]);
        set(LL1, 'Color', 'r')
        set(LL2, 'Color', 'b')
        set(LL3, 'Color', 'g')
        set(LL4, 'Color', 'm')

        LLmid = line([frontMid(1) rearMid(1)],[frontMid(2) rearMid(2)]);
        set(LLmid, 'Color', 'k','LineWidth',2)
        plot(frontMid(1),frontMid(2),'ko')
    end
    
    theta = atan2(frontMid(2)-rearMid(2),frontMid(1)-rearMid(1))*180/pi;
    [~, rotMat3] = rotAngle(-theta);
    thetaBody(i) = theta;

    
    LSHOrot = rotMat3*LSHO(i,:)';
    LELBrot = rotMat3*LELB(i,:)';
    LWRArot = rotMat3*LWRA(i,:)';
    LWRBrot = rotMat3*LWRB(i,:)';
    LFINrot = rotMat3*LFIN(i,:)';

    RSHOrot = rotMat3*RSHO(i,:)';
    RELBrot = rotMat3*RELB(i,:)';
    RWRArot = rotMat3*RWRA(i,:)';
    RWRBrot = rotMat3*RWRB(i,:)';
    RFINrot = rotMat3*RFIN(i,:)';

    LWRTmid = 0.5*(LWRArot+LWRBrot);
    RWRTmid = 0.5*(RWRArot+RWRBrot);    
    
    LASIrot = rotMat3*LASI(i,:)';
    RASIrot = rotMat3*RASI(i,:)';
    LPSIrot = rotMat3*LPSI(i,:)';
    RPSIrot = rotMat3*RPSI(i,:)';
    
    LHIPmid = 0.5*(LASIrot+LPSIrot);
    RHIPmid = 0.5*(RASIrot+RPSIrot);
    
    LTHIrot = rotMat3*LTHI(i,:)';
    LKNErot = rotMat3*LKNE(i,:)';
    LANKrot = rotMat3*LANK(i,:)';
    LHEErot = rotMat3*LHEE(i,:)';
    LTOErot = rotMat3*LTOE(i,:)';
    
    RTHIrot = rotMat3*RTHI(i,:)';
    RKNErot = rotMat3*RKNE(i,:)';
    RANKrot = rotMat3*RANK(i,:)';
    RHEErot = rotMat3*RHEE(i,:)';
    RTOErot = rotMat3*RTOE(i,:)';
    
    
    
    
    thetaZHandLeft(i) = atan2(LWRTmid(3)-LFINrot(3),LWRTmid(1)-LFINrot(1))*180/pi-90;
    thetaZHandRight(i) = atan2(RWRTmid(3)-RFINrot(3),RWRTmid(1)-RFINrot(1))*180/pi-90;     
    
    thetaZLowerArmLeft(i) = atan2(LELBrot(3)-LWRTmid(3),LELBrot(1)-LWRTmid(1))*180/pi-90;
    thetaZLowerArmRight(i) = atan2(RELBrot(3)-RWRTmid(3),RELBrot(1)-RWRTmid(1))*180/pi-90;

    thetaZUpperArmLeft(i) = atan2(LSHOrot(3)-LELBrot(3),LSHOrot(1)-LELBrot(1))*180/pi-90;
    thetaZUpperArmRight(i) = atan2(RSHOrot(3)-RELBrot(3),RSHOrot(1)-RELBrot(1))*180/pi-90;    
    
  
    
    
    thetaZShoeLeft(i) = atan2(LTOErot(3)-LHEErot(3),LTOErot(1)-LHEErot(1))*180/pi;
    thetaZShoeRight(i) = atan2(RTOErot(3)-RHEErot(3),RTOErot(1)-RHEErot(1))*180/pi;

    thetaZShankLeft(i) = atan2(LKNErot(3)-LANKrot(3),LKNErot(1)-LANKrot(1))*180/pi-90;
    thetaZShankRight(i) = atan2(RKNErot(3)-RANKrot(3),RKNErot(1)-RANKrot(1))*180/pi-90;

    thetaZThighLeft(i) = atan2(LHIPmid(3)-LKNErot(3),LHIPmid(1)-LKNErot(1))*180/pi-90;
    thetaZThighRight(i) = atan2(RHIPmid(3)-RKNErot(3),RHIPmid(1)-RKNErot(1))*180/pi-90;
    
end



figure(2)
thetaBodyPlot = thetaBody(thetaBody ~= 0);
plot(thetaBodyPlot)


thetaZHandLeft = thetaZHandLeft(~isnan(thetaZHandLeft));
thetaZHandRight = thetaZHandRight(~isnan(thetaZHandRight));
thetaZLowerArmLeft = thetaZLowerArmLeft(~isnan(thetaZLowerArmLeft));
thetaZLowerArmRight = thetaZLowerArmRight(~isnan(thetaZLowerArmRight));
thetaZUpperArmLeft = thetaZUpperArmLeft(~isnan(thetaZUpperArmLeft));
thetaZUpperArmRight = thetaZUpperArmRight(~isnan(thetaZUpperArmRight));


thetaZShoeLeft = thetaZShoeLeft(~isnan(thetaZShoeLeft));
thetaZShoeRight = thetaZShoeRight(~isnan(thetaZShoeRight));
thetaZShankLeft = thetaZShankLeft(~isnan(thetaZShankLeft));
thetaZShankRight = thetaZShankRight(~isnan(thetaZShankRight));
thetaZThighLeft = thetaZThighLeft(~isnan(thetaZThighLeft));
thetaZThighRight = thetaZThighRight(~isnan(thetaZThighRight));


shoulderAngleLeft = thetaZUpperArmLeft;
shoulderAngleRight = thetaZUpperArmRight;

elbowAngleLeft = thetaZLowerArmLeft - thetaZUpperArmLeft;
elbowAngleRight = thetaZLowerArmRight - thetaZUpperArmRight;

handAngleLeft = thetaZHandLeft - thetaZLowerArmLeft;
handAngleRight = thetaZHandRight - thetaZLowerArmRight;

hipAngleLeft = thetaZThighLeft;
hipAngleRight = thetaZThighRight;

kneeAngleLeft = thetaZShankLeft - thetaZThighLeft;
kneeAngleRight = thetaZShankRight - thetaZThighRight;

ankleAngleLeft = thetaZShoeLeft - thetaZShankLeft;
ankleAngleRight = thetaZShoeRight - thetaZShankRight;


% angleSubj = {thetaZShoeLeft, thetaZShoeRight, thetaZShankLeft, thetaZShankRight, thetaZThighLeft, thetaZThighRight,...
%     hipAngleLeft, hipAngleRight, kneeAngleLeft, kneeAngleRight, ankleAngleLeft, ankleAngleRight, thetaBody};

angleSubj = {thetaZShoeLeft, thetaZShoeRight, thetaZShankLeft, thetaZShankRight, thetaZThighLeft, thetaZThighRight,...
    thetaZLowerArmLeft, thetaZLowerArmRight, thetaZUpperArmLeft, thetaZUpperArmRight};

angleJointSubj = {hipAngleLeft, hipAngleRight, kneeAngleLeft, kneeAngleRight, ankleAngleLeft, ankleAngleRight,...
    thetaBody,...
    shoulderAngleLeft, shoulderAngleRight, elbowAngleLeft, elbowAngleRight, handAngleLeft, handAngleRight};
    


angleLabel = {'Shoe','Shank','Thigh','Hip','Knee','Ankle'};




end




function [rotMat2, rotMat3] = rotAngle(th)

rotMat2 = [cosd(th) -sind(th); sind(th) cosd(th)];
rotMat3 = [cosd(th) -sind(th) 0; sind(th) cosd(th) 0; 0 0 1];
end