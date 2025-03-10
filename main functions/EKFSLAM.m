function [EstTraject, GroundTruth, AIE, wp, lm, movingLMidx, METHOD, removeCheck] = EKFSLAM(varargin)
% EKFSLAM - Extended Kalman Filter Simultaneous Localization and Mapping (EKF-SLAM) function.
%
% DESCRIPTION:
%   This function estimates the trajectory of a robot's path using the EKF-SLAM algorithm and returns various relevant parameters.
% 
% SYNTAX:
%   [EstTraject, GroundTruth, RMSE, wp, lm, movingLMidx, METHOD, removeCheck] = EKFSLAM(wp, lm, movingLMidx, 'method', 'pro', 'draw', 'on')
%
% INPUTS:
%   wp          - Waypoints coordinates. If not provided, waypoints are randomly generated using the Traveling Salesman Problem solution.
%   lm          - Coordinates of the landmarks. If not provided, landmarks are placed randomly based on the waypoints' field of view.
%   movingLMidx - Indices of the dynamic landmarks.
%   'method'    - SLAM method selection: 'pro' (proposed) or 'conv' (conventional). Default is conventional.
%   'draw'      - Visualization option: 'on' (default) or 'off' to visualize the trajectory.
%   'seed'      - Fixed rng for same results.
%   'fixNdynLM' - Fixing the number of moving landmarks.
%   'rMethod'   - Removing method (A and B).
%
% OUTPUTS:
%   EstTraject  - Estimated trajectory of the robot path.
%   GroundTruth - Provides the ground truth for the EKF-SLAM, useful as a reference for the estimated trajectory.
%   RMSE        - Root Mean Square Error between EstTraject and GroundTruth.
%   wp          - Coordinates of the waypoints.
%   lm          - Coordinates of the landmarks.
%   movingLMidx - Indices of the dynamic landmarks.
%   METHOD      - Indicates the method used: 'proposed' (for removing moving landmarks) or 'conventional'.
%   removeCheck - Returns '1' if any landmark was removed by the proposed method, '0' otherwise.
%
% NOTES:
%   1. The function can be called with no inputs using EKFSLAM().
%   2. To design custom maps with waypoints, landmarks, and moving landmarks, use the mapDesigner.m function and then load the generated data into this function.
%
% EXAMPLE:
%   [EstTraject, GroundTruth] = EKFSLAM(wp, lm, movingLMidx, 'method', 'pro', 'draw', 'on');
%
% Author: Masoud DORVASH (masoud.dorvash@gmail.com) - Ali ESLAMIAN (alieslamian@gmail.com)
% Date: 09/11/2023

count = 0;
METHOD = 'pro';
draw_option = 'on';
seedFlag = false;
fixNdynLM_flag = false;
rMethodA = true;

% Check if there are additional input arguments and parse them
if ~isempty(varargin)
    for i = 1:length(varargin)

        if ischar(varargin{i}) && strcmpi(varargin{i}, 'draw')
            draw_option = varargin{i + 1};
        end

        if ischar(varargin{i}) && strcmpi(varargin{i}, 'method')
            METHOD = varargin{i + 1};
        end

        if ischar(varargin{i}) && strcmpi(varargin{i}, 'seed')
            seedFlag = true;
            seed = cell2mat(varargin(i+1));
        end

        if ischar(varargin{i}) && strcmpi(varargin{i}, 'fixNdynLM')
            fixNdynLM_flag = true;
            NMovLM = cell2mat(varargin(i+1));
        end

        if ischar(varargin{i}) && strcmpi(varargin{i}, 'rMethod')
            if ischar(varargin{i+1}) && strcmpi(varargin{i+1}, 'B')
                rMethodA = false;
            end
        end

    end

end

UsePropsedEKFSLAM = strcmpi(METHOD,'pro');
UseConvEKFSLAM = strcmpi(METHOD,'conv');

for k = 1:length(varargin)
    if isempty(varargin{k}) || isnumeric(varargin{k})
        count = count + 1;
    end
end

if count < 3

    if ~fixNdynLM_flag

        prompt = {'Number of moving landmarks (0 means stationary environmet):'};
        Title = 'Condition';
        dims = [1 60];
        NMovLM = inputdlg(prompt, Title, dims);
        NMovLM = str2double(NMovLM{1});

    end

    movingLMidx     = [];
    varMin          = 0;
    varMax          = 10;
    nWP             = 20;
    nLM             = 10;

    p               = RandPath(varMin,varMax);
    [xwp, ywp]      = adjustPoints(p.WayPoints(1,:), p.WayPoints(2,:), nWP);
    wp              = [xwp; ywp];
    thetawp         = zeros(1,length(xwp));

    [lm, sDev]       = RandLandmark(p,nLM);

else

    wp = varargin{1};
    lm = varargin{2};
    movingLMidx = varargin{3};
    xwp = wp(1,:);
    ywp = wp(2,:);
    thetawp         = zeros(1,size(wp,1));
    nLM = length(lm);
    DX = diff(xwp);
    DY = diff(ywp);
    p.Length = sum(sqrt(DX.^2 + DY.^2));
    sDev    = ceil(p.Length/size(wp,2));
    NMovLM = numel(movingLMidx);
    if NMovLM == 0
        MovLMFlag = 0;
    else
        MovLMFlag = 1;
    end

end

%%  Proposed Dynamic EKF-SLAM

if seedFlag == true
    rng(seed)
end

if UsePropsedEKFSLAM

    %%  EKF Initializing

    Xtrue           = [xwp(1), ywp(1), thetawp(1)]';
    xLand           = lm(1,:);
    yLand           = lm(2,:);

    % - Initialize EKF estimate of state and covariance
    Xest            = Xtrue;
    Pest            = zeros(3, 3);

    % - Define robot's state noise (assumed to be Gaussian)
    sigmaX          = 0.01;
    sigmaY          = 0.01;
    sigmaTheta      = 0.01;

    R               = diag([sigmaX^2 sigmaY^2 sigmaTheta^2]);

    % - Landmark measurement noise
    sigmaM          = 0.02;

    % - Define noise covariance for landmarks
    R_lm            = sigmaM*eye(2);

    % - Dynamic landmarks initialization
    vLM             = 6;

    % - Define time step
    dt              = 0.1;

    % - Create cell array to store true and estimated states over time
    GroundTruth     = cell(1,round(2.5 * p.Length / dt));
    EstTraject      = cell(1,round(2.5 * p.Length / dt));
    LandMark        = cell(1,round(2.5 * p.Length / dt));
    GroundTruth{1}  = Xtrue;
    EstTraject{1}   = Xest;
    dist            = zeros(nLM,1);
    counterLM       = zeros(size(lm,2));
    mStep           = 2;
    dirLM           = zeros(size(lm,2),2);

    % - Dynamic landmarks initialization
    sigmaD          = 1;
    mat2            = [inf,inf,inf];

    % - Moving Landmark Indices
    if NMovLM > length(lm)
        MovLMFlag = 1;
        numbers = randperm(10);
        movingLMidx = numbers(1:length(lm));
        disp('Moving landmarks are more than total number of landmarks')
        disp('All landmarks are assumed to be dynamic')
    elseif NMovLM == 0 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 0;
    elseif NMovLM == 1 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 1;
        movingLMidx = randi([1, size(lm,2)]);
    elseif NMovLM == 2 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 1;
        LMDistance = pdist(lm');
        [~, idxMaxDis] = max(LMDistance);
        x1Idx = 1;
        while idxMaxDis > (nLM-x1Idx)
            idxMaxDis = idxMaxDis - (nLM-x1Idx);
            x1Idx = x1Idx + 1;
        end
        x2Idx           = x1Idx + idxMaxDis;
        movingLMidx     = [x1Idx x2Idx];
    elseif NMovLM >= 3 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 1;
        numbers = randperm(NMovLM);
        movingLMidx = numbers(1:NMovLM);
    end

    counterStep     = 1;
    removeCheck     = '0';

    % - Moving landmarks path
    [centerGTx, centerGTy] = centroid(polyshape(xwp,ywp));
    distLand = [(xLand-centerGTx); (yLand-centerGTy)];
    [~,idxMove] = max(abs(distLand), [],1);

    selIdx = sub2ind(size(dirLM),1:size(lm,2),idxMove);
    dirLM(selIdx) = -sign(distLand(selIdx));

    %%  Main Loop

    for i = 1:length(xwp)

        dx = xwp(i) - Xtrue(1);
        dy = ywp(i) - Xtrue(2);

        while sqrt(dx^2 + dy^2) > 0.01

            Xest = Xest(1:3);
            Pest = Pest(1:3,1:3);
            R = R(1:3,1:3);

            idxP = zeros(nLM,1);
            LandOld = [xLand; yLand];
            DistOld = dist;

            % - Check distance to each landmark
            for j = 1:nLM

                % - Moving landmark path
                if any(j == movingLMidx)

                    inCheck = inpolygon(xLand(j),yLand(j),xwp,ywp);

                    if counterLM(j) < mStep

                        xLand(j) = xLand(j)+ vLM * dt * dirLM(j,1);
                        yLand(j) = yLand(j)+ vLM * dt * dirLM(j,2);
                        counterLM(j) = counterLM(j) + 1;

                    else
                        counterLM(j) = 0;
                        if inCheck
                            thetaLM = atan2(dirLM(j,2),dirLM(j,1)) + randi([1 3]) * pi/2 + pi;
                            dirLM(j,:) = [cos(thetaLM); sin(thetaLM)];
                        else
                            distLand = [(xLand(j)-centerGTx); (yLand(j)-centerGTy)];
                            [~,idxMove] = max(abs(distLand), [],1);
                            dirLM(j,:) = zeros(size(lm,1),1);
                            dirLM(j,idxMove) = -sign(distLand(idxMove));
                        end
                        xLand(j) = xLand(j)+ vLM * dt * dirLM(j,1);
                        yLand(j) = yLand(j)+ vLM * dt * dirLM(j,2);
                        counterLM(j) = counterLM(j) + 1;
                    end

                end

                LandMark{counterStep}(j,:) = [xLand(j) yLand(j)];

                dist(j) = sqrt((Xtrue(1) - xLand(j))^2 + (Xtrue(2) - yLand(j))^2) + sigmaM*randn;

                if rMethodA

                    if dist(j) <= sDev

                        idxP(j) = 1;

                        if sum(idxP) >= 3

                            idx = find(idxP);

                            zCoor = [idx, xLand(idx)'+sigmaM*randn(size(idx)), yLand(idx)'+sigmaM*randn(size(idx))];
                            Distance = pdist([zCoor(:,2),zCoor(:,3)])';
                            combinations = nchoosek(zCoor(:,1), 2);
                            mat1 = [combinations,Distance];
                            [mismatch] = detectMLM(mat1, mat2,sigmaD);
                            mat2 = mat1;
                            [~,index] = removedynamic(mismatch,zCoor);

                            if ~isempty(index)
                                removeCheck = '1';
                            end

                            % - Add the landmark to the state vector and covariance
                            if j ~= mismatch
                                XestP = Xest;
                                Xest = [XestP; sigmaM*randn+xLand(j); sigmaM*randn+yLand(j)];
                                Pest = blkdiag(Pest, 0.01*eye(2));
                                R = blkdiag(R, R_lm);
                            end

                        end

                    else

                        % - Remove dynamic landmarks
                        if dist(j) <= sDev && DistOld(j) <= sDev && counterStep >= 2
                            Angle = angleBetweenPoints([LandOld(:,j)],[GroundTruth{counterStep-1}(1) GroundTruth{counterStep-1}(2)], ...
                                [GroundTruth{counterStep}(1) GroundTruth{counterStep}(2)]);
                            estDistance = estDist(DistOld(j), v*dt, Angle);
                            if abs(estDistance - dist(j)) >= threshold
                                removeCheck = '1';
                                continue
                            end
                        end

                        % - Add the landmark to the state vector and covariance
                        if dist(j) <= sDev
                            XestP = Xest;
                            Xest = [XestP; sigmaM*randn+xLand(j); sigmaM*randn+yLand(j)];
                            Pest = blkdiag(Pest, 0.01*eye(2));
                            R = blkdiag(R, R_lm);
                        end

                    end

                end

            end

            thetaGoal = atan2(dy, dx);
            omega = wrapToPi(thetaGoal - Xtrue(3)) / dt;
            v = min(sqrt(dx^2 + dy^2) / dt, 0.5);
            u = [v; omega];

            Xtrue = Xtrue + [v*cos(Xtrue(3))*dt; v*sin(Xtrue(3))*dt; omega*dt] + chol(R(1:3,1:3))*randn(3,1)*sqrt(dt);

            % - Predict the state and covariance
            Xest(1:3) = Xest(1:3) + [u(1)*cos(Xest(3))*dt; u(1)*sin(Xest(3))*dt; u(2)*dt];
            F = [1 0 -u(1)*sin(Xest(3))*dt; 0 1 u(1)*cos(Xest(3))*dt; 0 0 1];
            F = blkdiag(F, eye(size(Xest,1)-3));
            Pest = F * Pest * F' + R;

            % - Correction Step
            for j = 4:2:length(Xest)
                dx = Xest(j) - Xest(1);
                dy = Xest(j+1) - Xest(2);
                q = dx^2 + dy^2;
                z = [sqrt(q); atan2(dy, dx) - Xest(3)];

                % Jacobian of the measurement model
                H = [-dx/sqrt(q), -dy/sqrt(q), 0, dx/sqrt(q), dy/sqrt(q), zeros(1, length(Xest) - 5);
                    dy/q, -dx/q, -1, -dy/q, dx/q, zeros(1, length(Xest) - 5)];

                K = Pest * H' / (H * Pest * H' + R_lm);
                Xest = Xest + K * (z - [sqrt(q); atan2(dy, dx) - Xest(3)]);
                Pest = (eye(size(K,1)) - K * H) * Pest;
            end

            Xest(3) = wrapToPi(Xest(3));

            counterStep = counterStep + 1;

            GroundTruth{counterStep} = Xtrue;
            EstTraject{counterStep} = Xest;

            dx = xwp(i) - Xtrue(1);
            dy = ywp(i) - Xtrue(2);

        end

    end

    %%  Post Processing

    % - Removing zeros (empty cells)
    GroundTruth(counterStep+1:end) = [];
    EstTraject(counterStep+1:end) = [];

    % - Convert true histogram into a matrix
    GroundTruth = [GroundTruth{1:end}];

    % - Select the first 3 elements of each cell (plot)
    EstTrajTrunc    = cellfun(@(x) x(1:min(3, end)), EstTraject, 'UniformOutput', false);

    % - Ensure that each cell contains exactly 3 elements (plot)
    EstTrajPadded   = cellfun(@(x) [x, zeros(1, max(0, 3-length(x)))], EstTrajTrunc, 'UniformOutput', false);

    % Convert the estimated histogram into a matrix (plot)
    EstTraject         = cell2mat(EstTrajPadded);


    %%  Results

    AIE = calculate_RMSE(EstTraject(1:2,:),GroundTruth(1:2,:));
    
    % - Plot Results
    if strcmpi(draw_option, 'on')

        ScreenSize  = get(0,'ScreenSize');
        figure('Name','Map','NumberTitle','off',...
            'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2.5)]);
        movegui('center')
        plot(GroundTruth(1, :), GroundTruth(2, :), 'LineWidth', 1.5), hold on
        plot(EstTraject(1, :), EstTraject(2, :), 'LineWidth', 1.5);

        scatter(lm(1,:),lm(2,:),70,'MarkerEdgeColor',[0 .5 .5], ...
            'MarkerFaceColor',[0 .7 .7], ...
            'LineWidth',2)
        if MovLMFlag ~= 0
            scatter(lm(1,movingLMidx),lm(2,movingLMidx),70,'MarkerEdgeColor',[0 .5 .5], ...
                'MarkerFaceColor','y', ...
                'LineWidth',2)
        end
        xvec = {GroundTruth(1, :), EstTraject(1, :), lm(1,:)};
        yvec = {ywp, GroundTruth(2, :), EstTraject(2, :), lm(2,:)};
        xmin = min(cellfun(@min, xvec));
        xmax = max(cellfun(@max, xvec));
        ymin = min(cellfun(@min, yvec));
        ymax = max(cellfun(@max, yvec));
        xlim([floor(xmin)-1, ceil(xmax)+1])
        ylim([floor(ymin)-1, ceil(ymax)+1])
        xlabel('$x$ (m)','Interpreter','latex','FontSize',12')
        ylabel('$y$ (m)','Interpreter','latex','FontSize',12')
        if MovLMFlag ~= 0
            legend('Ground Truth', 'Estimated Path', 'Stationary Landmarks', 'Moving Landmarks' ,'Interpreter','latex','Location','northoutside','FontSize',11)
        else
            legend('Ground Truth', 'Estimated Path', 'Stationary Landmarks', 'Interpreter','latex','Location','northoutside','FontSize',11)
        end
        grid on

    end

end

%%  Static EKF-SLAM

if seedFlag == true
    rng(seed)
end

if UseConvEKFSLAM

    %%  EKF Initializing

    Xtrue           = [xwp(1), ywp(1), thetawp(1)]';
    xLand           = lm(1,:);
    yLand           = lm(2,:);

    % - Initialize EKF estimate of state and covariance
    Xest            = Xtrue;
    Pest            = zeros(3, 3);

    % - Define robot's state noise (assumed to be Gaussian)
    sigmaX          = 0.01;
    sigmaY          = 0.01;
    sigmaTheta      = 0.01;

    R  = diag([sigmaX^2 sigmaY^2 sigmaTheta^2]);

    % - Landmark measurement noise
    sigmaM          = 0.02;

    % - Define noise covariance for landmarks
    R_lm            = sigmaM*eye(2);

    % - Dynamic landmarks initialization
    vLM             = 6;

    % - Define time step
    dt              = 0.1;

    % - Create cell array to store true and estimated states over time
    GroundTruth     = cell(1,round(2.5 * p.Length / dt));
    EstTraject      = cell(1,round(2.5 * p.Length / dt));
    GroundTruth{1}  = Xtrue;
    EstTraject{1}   = Xest;
    counterLM       = zeros(size(lm,2));
    mStep           = 2;
    dirLM           = zeros(size(lm,2),2);

    % - Moving Landmark Indices
    if NMovLM > length(lm)
        MovLMFlag = 1;
        numbers = randperm(10);
        movingLMidx = numbers(1:length(lm));
        disp('Moving landmarks are more than total number of landmarks')
        disp('All landmarks are assumed to be dynamic')
    elseif NMovLM == 0 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 0;
    elseif NMovLM == 1 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 1;
        movingLMidx = randi([1, size(lm,2)]);
    elseif NMovLM == 2 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 1;
        LMDistance = pdist(lm');
        [~, idxMaxDis] = max(LMDistance);
        x1Idx = 1;
        while idxMaxDis > (nLM-x1Idx)
            idxMaxDis = idxMaxDis - (nLM-x1Idx);
            x1Idx = x1Idx + 1;
        end
        x2Idx           = x1Idx + idxMaxDis;
        movingLMidx     = [x1Idx x2Idx];
    elseif NMovLM >= 3 && isempty(movingLMidx) && NMovLM <= length(lm)
        MovLMFlag = 1;
        numbers = randperm(NMovLM);
        movingLMidx = numbers(1:NMovLM);
    end

    counterStep     = 1;
    removeCheck     = '0';


    % - Moving landmarks path
    [centerGTx, centerGTy] = centroid(polyshape(xwp,ywp));
    distLand = [(xLand-centerGTx); (yLand-centerGTy)];
    [~,idxMove] = max(abs(distLand), [],1);

    selIdx = sub2ind(size(dirLM),1:size(lm,2),idxMove);
    dirLM(selIdx) = -sign(distLand(selIdx));

    %%  Main Loop

    for i = 1:length(xwp)

        dx = xwp(i) - Xtrue(1);
        dy = ywp(i) - Xtrue(2);

        while sqrt(dx^2 + dy^2) > 0.01

            Xest = Xest(1:3);
            Pest = Pest(1:3,1:3);
            R = R(1:3,1:3);

            % - Check distance to each landmark




            for j = 1:nLM

                % - Moving landmark path
                if any(j == movingLMidx)

                    inCheck = inpolygon(xLand(j),yLand(j),xwp,ywp);

                    if counterLM(j) < mStep

                        xLand(j) = xLand(j)+ vLM * dt * dirLM(j,1);
                        yLand(j) = yLand(j)+ vLM * dt * dirLM(j,2);
                        counterLM(j) = counterLM(j) + 1;

                    else
                        counterLM(j) = 0;
                        if inCheck
                            thetaLM = atan2(dirLM(j,2),dirLM(j,1)) + randi([1 3]) * pi/2 + pi;
                            dirLM(j,:) = [cos(thetaLM); sin(thetaLM)];
                        else
                            distLand = [(xLand(j)-centerGTx); (yLand(j)-centerGTy)];
                            [~,idxMove] = max(abs(distLand), [],1);
                            dirLM(j,:) = zeros(size(lm,1),1);
                            dirLM(j,idxMove) = -sign(distLand(idxMove));
                        end
                        xLand(j) = xLand(j)+ vLM * dt * dirLM(j,1);
                        yLand(j) = yLand(j)+ vLM * dt * dirLM(j,2);
                        counterLM(j) = counterLM(j) + 1;
                    end

                end

                dist(j) = sqrt((Xtrue(1) - xLand(j))^2 + (Xtrue(2) - yLand(j))^2) + sigmaM*randn;

                if dist <= sDev

                    % - Add the landmark to the state vector and covariance

                    XestP = Xest;
                    Xest = [XestP; sigmaM*randn+xLand(j); sigmaM*randn+yLand(j)];
                    Pest = blkdiag(Pest, 0.01*eye(2));
                    R = blkdiag(R, R_lm);

                end

            end

            thetaGoal = atan2(dy, dx);
            omega = wrapToPi(thetaGoal - Xtrue(3)) / dt;
            v = min(sqrt(dx^2 + dy^2) / dt, 0.5);
            u = [v; omega];

            Xtrue = Xtrue + [v*cos(Xtrue(3))*dt; v*sin(Xtrue(3))*dt; omega*dt] + chol(R(1:3,1:3))*randn(3,1)*sqrt(dt);

            % - Predict the state and covariance
            Xest(1:3) = Xest(1:3) + [u(1)*cos(Xest(3))*dt; u(1)*sin(Xest(3))*dt; u(2)*dt];
            F = [1 0 -u(1)*sin(Xest(3))*dt; 0 1 u(1)*cos(Xest(3))*dt; 0 0 1];
            F = blkdiag(F, eye(size(Xest,1)-3));
            Pest = F * Pest * F' + R;

            % - Correction Step
            for j = 4:2:length(Xest)
                dx = Xest(j) - Xest(1);
                dy = Xest(j+1) - Xest(2);
                q = dx^2 + dy^2;
                z = [sqrt(q); atan2(dy, dx) - Xest(3)];

                % Jacobian of the measurement model
                H = [-dx/sqrt(q), -dy/sqrt(q), 0, dx/sqrt(q), dy/sqrt(q), zeros(1, length(Xest) - 5);
                    dy/q, -dx/q, -1, -dy/q, dx/q, zeros(1, length(Xest) - 5)];

                K = Pest * H' / (H * Pest * H' + R_lm);
                Xest = Xest + K * (z - [sqrt(q); atan2(dy, dx) - Xest(3)]);
                Pest = (eye(size(K,1)) - K * H) * Pest;
            end

            Xest(3) = wrapToPi(Xest(3));

            counterStep = counterStep + 1;

            GroundTruth{counterStep} = Xtrue;
            EstTraject{counterStep} = Xest;

            dx = xwp(i) - Xtrue(1);
            dy = ywp(i) - Xtrue(2);

        end

    end

    %%  Post Processing

    % - Removing zeros (empty cells)
    GroundTruth(counterStep+1:end) = [];
    EstTraject(counterStep+1:end) = [];

    % - Convert true histogram into a matrix
    GroundTruth = [GroundTruth{1:end}];

    % - Select the first 3 elements of each cell (plot)
    EstTrajTrunc    = cellfun(@(x) x(1:min(3, end)), EstTraject, 'UniformOutput', false);

    % - Ensure that each cell contains exactly 3 elements (plot)
    EstTrajPadded   = cellfun(@(x) [x, zeros(1, max(0, 3-length(x)))], EstTrajTrunc, 'UniformOutput', false);

    % Convert the estimated histogram into a matrix (plot)
    EstTraject         = cell2mat(EstTrajPadded);


    %%  Results

    AIE = calculate_RMSE(EstTraject(1:2,:),GroundTruth(1:2,:));

    % - Plot Results
    if strcmpi(draw_option, 'on')

        ScreenSize  = get(0,'ScreenSize');
        figure('Name','Map','NumberTitle','off',...
            'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2.5)]);
        movegui('center')
        plot(GroundTruth(1, :), GroundTruth(2, :), 'LineWidth', 1.5), hold on
        plot(EstTraject(1, :), EstTraject(2, :), 'LineWidth', 1.5);

        scatter(lm(1,:),lm(2,:),70,'MarkerEdgeColor',[0 .5 .5], ...
            'MarkerFaceColor',[0 .7 .7], ...
            'LineWidth',2)
        if MovLMFlag ~= 0
            scatter(lm(1,movingLMidx),lm(2,movingLMidx),70,'MarkerEdgeColor',[0 .5 .5], ...
                'MarkerFaceColor','y', ...
                'LineWidth',2)
        end
        xvec = {GroundTruth(1, :), EstTraject(1, :), lm(1,:)};
        yvec = {ywp, GroundTruth(2, :), EstTraject(2, :), lm(2,:)};
        xmin = min(cellfun(@min, xvec));
        xmax = max(cellfun(@max, xvec));
        ymin = min(cellfun(@min, yvec));
        ymax = max(cellfun(@max, yvec));
        xlim([floor(xmin)-1, ceil(xmax)+1])
        ylim([floor(ymin)-1, ceil(ymax)+1])
        xlabel('$x$ (m)','Interpreter','latex','FontSize',12')
        ylabel('$y$ (m)','Interpreter','latex','FontSize',12')
        if MovLMFlag ~= 0
            legend('Ground Truth', 'Estimated Path', 'Stationary Landmarks', 'Moving Landmarks' ,'Interpreter','latex','Location','northoutside','FontSize',11)
        else
            legend('Ground Truth', 'Estimated Path', 'Stationary Landmarks', 'Interpreter','latex','Location','northoutside','FontSize',11)
        end
        grid on
    end

end

end