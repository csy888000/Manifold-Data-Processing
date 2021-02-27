clear
close all

addpath IMU_data
addpath gpdm gpdm/netlab gpdm/gpdm gpdm/gplvm gpdm/util

format long

global USE_GAMMA_PRIOR  % gamma prior for dynamics, only works with RBF kernel
global GAMMA_ALPHA % defines shape of the gamma prior
global USE_LAWRENCE % fix dynamics HPs, as Lawrence suggested (use with thetad = [0.2 0.01 1e6];) 
global FIX_HP % fix all HPs
global MARGINAL_W % marginalize over W while learning X
global MARGINAL_DW % marginalize over scale in dynamics while learning X
global LEARN_SCALE % use different scales for different output dimensions
global REMOVE_REDUNDANT_SCALE % let W absorb the overall scale of reconstruction
global W_VARIANCE % kappa^2 in the paper, not really the variance though
global M_CONST % M value in Jack's master's thesis
global BALANCE % Constant in front of dynamics term, set to D/q for the B-GPDM
global SUBSET_SIZE % Number of data to select for EM, set -1 for all data. 
global USE_OLD_MISSING_DATA

M_CONST = 1; 
REMOVE_REDUNDANT_SCALE = 1;
LEARN_SCALE = 1; 
MARGINAL_W = 0; 
MARGINAL_DW = 0; 
W_VARIANCE = 1e6; 
FIX_HP = 0; 
USE_GAMMA_PRIOR = 0; 
GAMMA_ALPHA = [5 10 2.5]; 
USE_LAWRENCE = 0;
BALANCE = 1;
SUBSET_SIZE = -1; 


opt = foptions;
opt(1) = 1;
opt(9) = 0;
if MARGINAL_W == 1
    opt(14) = 100; % total number of iterations
    extItr = 1; 
else
    opt(14) = 10; % rescaling every 10 iterations
    extItr = 100; % do extItr*opt(14) iterations in total
end  

% modelType(1) : input of dynamics
%   0 => [x_t, x_{t-1}]
%   1 => [x_t, x_t - x_{t-1}]
%   2 => [x_t]
% modelType(2) : output of dynamics 
%   0 => x_{t+1} 
%   1 => x_{t+1} - x_t
% modelType(3) : kernel type
%   0 => RBF kernel with weighted dimensions, use with input 0 or 1
%   1 => RBF kernel 
%   2 => Linear kernel
%   3 => weighted Linear kernel + RBF kernel with weighted dimensions, use with
%   input 0 or 1
%   4 => weighted linear kernel
%   5 => linear + RBF

% Learn single walker model from lin+rbf kernel.
modelType = [2 0 5]; 
% [Y, initY, varY, segments] = loadMocapData({['07_01.amc']}, [1], [2],[260]);

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % treadmill speed 1.25m/s
% test1.viconFilename = 'Subj_walk_v1_Test1.csv';
% test1.testObj = getVicondata(test1.viconFilename);
% test1.testNum = 1;
% range = 5178:1:5451;
% test1.range = range;
% disp('finish: 1')
% 
% % testDataAll = {test1};
% % disp('Saving...')
% % save('testDataAll.mat','testDataAll')
% % disp('Save done.')
% 




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



jointAngle = cell(1,length(testDataAll));
for i = 1:length(testDataAll)
    range = testDataAll{i}.range;
    
    % thetaZShoeLeft = testDataAll{1}.testObj.angleData{1}(range);
    % thetaZShoeRight = testDataAll{1}.testObj.angleData{2}(range);
    % thetaZShankLeft = testDataAll{1}.testObj.angleData{3}(range);
    % thetaZShankRight = testDataAll{1}.testObj.angleData{4}(range);
    % thetaZThighLeft = testDataAll{1}.testObj.angleData{5}(range);
    % thetaZThighRight = testDataAll{1}.testObj.angleData{6}(range);

    hipLeftSeg = testStruct{i}.testObj.angleJointData{1}(:);
    hipRightSeg = testStruct{i}.testObj.angleJointData{2}(:);
    kneeLeftSeg = testStruct{i}.testObj.angleJointData{3}(:);
    kneeRightSeg = testStruct{i}.testObj.angleJointData{4}(:);
    ankleLeftSeg = testStruct{i}.testObj.angleJointData{5}(:);
    ankleRightSeg = testStruct{i}.testObj.angleJointData{6}(:);


    % jointAngle = [thetaZShoeLeft, thetaZShoeRight, thetaZShankLeft, thetaZShankRight, thetaZThighLeft, thetaZThighRight,...
    %     hipLeftSeg, hipRightSeg, kneeLeftSeg, kneeRightSeg, ankleLeftSeg, ankleRightSeg];

    jointAngle{i} = [hipLeftSeg, hipRightSeg, kneeLeftSeg, kneeRightSeg, ankleLeftSeg, ankleRightSeg];
end

figure(21)
clf
colorX = {'b','k','r','g','m','c','y','k'};
for i = 1:length(jointAngle)
    for j = 1:6
        subplot(3,2,j)
        hold on
        plot(jointAngle{i}(:,j),colorX{i})
    end
end




%%
walkX = cell(1,length(jointAngle));
walkX_pred = cell(1,length(jointAngle));
for j = 1:length(jointAngle)
    Y = jointAngle{j};
    segments = 1;

    missing = [];
    N = size(Y, 1); D = size(Y, 2);
    q = 3; % dimensionality of latent space

    % PCA
    X = zeros(N, q);
    refY = Y; meanData = mean(Y);
    Y = Y - repmat(meanData, N, 1);
    [v, u] = pcaGP(Y);
    [coeff,score,latent] = pca(Y);
    v(v<0)=0;
    X = Y*u(:, 1:q)*diag(1./sqrt(v(1:q)));

    % initialize hyperparameters
    theta = [1 1 exp(1)];
    thetad = [0.9 1 0.1 exp(1)];
    w = ones(D,1);

    [X,theta,thetad,w] = gpdmfitFull(X, Y, w, segments, theta, thetad, opt, ... 
         extItr, modelType, missing);

    % save example_model X Y w theta thetad modelType N D q meanData  ...
    % segments initY varY missing refY;

    % save example_model X Y w theta thetad modelType N D q meanData  ...
    % segments missing refY;

    % generate samples from learned model

    % load example_model
    [K,invK] = computeKernel(X, theta);
    [Xin,Xout] = priorIO(X, segments, modelType);
    [Kd,invKd] = computePriorKernel(Xin, thetad, modelType(3));
    simSteps = 256;
    % starts at ened of training sequence;
    simStart = [X(segments(1)+1,:) X(end,:)]; %  inputs 2 points in case using 2nd order model
    [X_pred,XRand_pred] = simulatedynamics(X, segments, thetad, invKd, simSteps, simStart, modelType);

    % uncomment if want to generate new samples

    % hmcopt = foptions;      % Default options vector.
    % hmcopt(1) = 1;			% Switch on diagnostics.
    % hmcopt(7) = 100;	    	% Number of steps in leapfrog trajectory.
    % hmcopt(9) = 0; 
    % hmcopt(14) = 60;		% Number of Monte Carlo samples returned. 
    % hmcopt(15) = 20;		% Number of samples omitted at start of chain.
    % hmcopt(18) = 0.01;  	% leapfrog step size
    % X_samples = sampledynamics('example_model', X_pred, 'samples', hmcopt);
    
    walku{j} = u;
    walkv{j} = v;
    walkY{j} = Y;
    walkX{j} = X;
    walkX_pred{j} = X_pred;
    
end

figure(31)
clf
colorX = {'b','k','r','g','m','c'};
for i = 1:length(jointAngle)
    for j = 1:3
        subplot(3,1,j)
        hold on
        plot(walkX{i}(:,j),colorX{i})
    end
end


save walkLatent walkX walkX_pred



%%
% 
% load walkSamples
% load walkX
% load walkX_pred


load walkLatent
segments = 1;

figure(10)
clf;
hold on;

colorX = {'b','k','r','g','m','c'};
colorXpred = {'b','k','r','g','m','c'};
shapeU = {'d','s','x','+','*','^'};

walk_pred_rot = walkX_pred;
for k = [2 4 5 6]%1:length(jointAngle)
    
    X = walkX{k};
    X_pred = walkX_pred{k};
    % for n=1:4:size(X_samples,2)
    %     X_samples{n}(:,[1 2 3]) = X_samples{n}(:,[3 1 2]);
    %     X_samples{n}(:,3) = -X_samples{n}(:,3);
    %     X_samples{n}(:,1) = -X_samples{n}(:,1);
    % end
    X(:,[1 2 3]) = X(:,[3 1 2]);
    X(:,3) = -X(:,3);
    
    X_pred(:,[1 2 3]) = X_pred(:,[3 1 2]);
    X_pred(:,3) = -X_pred(:,3);
    if k == 2
        X(:,1) = -X(:,1);
        X_pred(:,1) = -X_pred(:,1);

    end
    if k == 3
        X(:,1) = X(:,1);
        X_pred(:,1) = X_pred(:,1);
        X(:,2) = X(:,2);
        X_pred(:,2) = X_pred(:,2);
    end

    if k == 1
        X = X*rotAngle(10,-10,0);
        X_pred = X_pred*rotAngle(10,-10,0);
    end
    
    if k == 3
        X = X*rotAngle(5,10,0);
        X_pred = X_pred*rotAngle(-5,10,0);
    end
    
    if k == 4
        X = X*rotAngle(10,130,-50);
        X_pred = X_pred*rotAngle(10,130,-50);
    end

    
    if k == 5
        X = X*rotAngle(-5,80,0);
        X_pred = X_pred*rotAngle(-5,80,0);
    end
    
    if k == 6
        X = X*rotAngle(-5,30,0);
        X_pred = X_pred*rotAngle(-5,30,0);
    end
    
    % for n=1:4:size(X_samples,2)
    %     plotseries(X_samples{n}, [1], 'g');
    % end

%     plotseries(X, segments, colorX{k});
    walk_pred_rot{k} = X_pred;
    plotseries(X_pred, [1], colorXpred{k});
    
    hold on
%     plot3(0,0,0,'d','Color',colorX{i})
%     plot3(walku{k}(1:3,1),walku{k}(1:3,2),walku{k}(1:3,3),shapeU{k},'Color',colorX{k})
%     for j = 1%1:3
%         line([0 walku{k}(j,1)],[0 walku{k}(j,2)],[0 walku{k}(j,3)])
%     end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    view(9,8)

end

text(-1.5,-2, -1.5,'Blue: speed 1.25m/s  0 degree')
text(-1.5,-2, -1.6,'Black: speed 1.5m/s  0 degree')
text(-1.5,-2, -1.7,'Red: speed 1.75m/s  0 degree')
text(-1.5,-2, -1.8,'Green: Walking up slope 20 degree')
text(-1.5,-2, -1.9,'Megenta: Walking across slope 20 degree')
text(-1.5,-2, -2,'Cyan: Walking down slope 20 degree')


%%
REF_CURVE = cell(1,length(walk_pred_rot));

figure(11)
clf
for i = 1:6 %[2 4 5 6]
    predCurve = walk_pred_rot{i};
    firstHalf = predCurve(1:128,:);
    secondHalf = predCurve(129:end,:);
    meanPredCurve = (firstHalf+secondHalf)./2;
    REF_CURVE{i} = meanPredCurve;
    for j = 1:3
        subplot(3,1,j)
        hold on
        plot(meanPredCurve(:,j))
        
    end
end   

for i = 1:6
    REF_CURVE_FLAT(:,i*4-3:i*4-1)  = REF_CURVE{i};
    if i < 4
        REF_CURVE_FLAT(:,i*4) = ones(length(REF_CURVE{i}),1);
    else
        REF_CURVE_FLAT(:,i*4) = ones(length(REF_CURVE{i}),1)*(i-2);
    end
end

writematrix(REF_CURVE_FLAT,'REF_CURVE_FLAT.csv') 


%%

function rotMat3 = rotAngle(x, y, z)

rotMatX = [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];
rotMatY = [cosd(y) 0 -sind(y); 0 1 0; sind(y) 0 cosd(y)];
rotMatZ = [cosd(z) sind(z) 0; -sind(z) cosd(z) 0; 0 0 1];

rotMat3 = rotMatX*rotMatY*rotMatZ;
end


