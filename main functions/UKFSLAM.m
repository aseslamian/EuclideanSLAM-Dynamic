function [EstTraject, GroundTruth, AIE, wp, lm, movingLMidx, METHOD, removeCheck] = UKFSLAM(varargin)
% UKFSLAM - Unscented Kalman Filter Simultaneous Localization and Mapping (UKF-SLAM) function.
%
% DESCRIPTION:
%   This function estimates the trajectory of a robot's path using the UKF-SLAM algorithm and returns various relevant parameters.
%
% SYNTAX:
%   [EstTraject, GroundTruth, RMSE, wp, lm, movingLMidx, METHOD, removeCheck] = UKFSLAM(wp, lm, movingLMidx, 'method', 'pro', 'draw', 'on', 'rtm')
%
% INPUTS:
%   wp          - Waypoints coordinates. If not provided, waypoints are randomly generated using the Traveling Salesman Problem solution.
%   lm          - Coordinates of the landmarks. If not provided, landmarks are placed randomly based on the waypoints' field of view.
%   movingLMidx - Indices of the dynamic landmarks.
%   'method'    - SLAM method selection: 'pro' (proposed) or 'conv' (conventional). Default is conventional.
%   'draw'      - Visualization option: 'on' (default) or 'off' to visualize the trajectory.
%   'rtm'       - Option for real-time mapping. When selected, the function will plot the robot's movement in real-time along with the landmarks, ground truth, and the robot's point of view for landmark detection.
%   'seed'      - Fixed rng for same results.
%   'fixNdynLM' - Fixing the number of moving landmarks.
%   'rMethod'   - Removing method (A and B).
%
% OUTPUTS:
%   EstTraject  - Estimated trajectory of the robot path.
%   GroundTruth - Provides the ground truth for the UKF-SLAM, useful as a reference for the estimated trajectory.
%   RMSE        - Root Mean Square Error between EstTraject and GroundTruth.
%   wp          - Coordinates of the waypoints.
%   lm          - Coordinates of the landmarks.
%   movingLMidx - Indices of the dynamic landmarks.
%   METHOD      - Indicates the method used: 'proposed' (for removing moving landmarks) or 'conventional'.
%   removeCheck - Returns '1' if any landmark was removed by the proposed method, '0' otherwise.
%
% NOTES:
%   1. The function can be called with no inputs using UKFSLAM().
%   2. To design custom maps with waypoints, landmarks, and moving landmarks, use the mapDesigner.m function and then load the generated data into this function.
%
% EXAMPLE:
%   [EstTraject, GroundTruth] = UKFSLAM(wp, lm, movingLMidx, 'method', 'pro', 'draw', 'on');
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

UsePropsedUKFSLAM = strcmpi(METHOD,'pro');
UseConvUKFSLAM = strcmpi(METHOD,'conv');

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
    p.Length = 0;
    for i = 1:length(xwp)-1
        distance = sqrt((xwp(i+1) - xwp(i))^2 + (ywp(i+1) - ywp(i))^2);
        p.Length = p.Length + distance;
    end
    sDev    = max(ceil(p.Length/20),2.5);
    NMovLM = numel(movingLMidx);
    if NMovLM == 0
        MovLMFlag = 0;
    else
        MovLMFlag = 1;
    end

end


%%  Proposed Dynamic UKF-SLAM

if seedFlag == true
    rng(seed)
end

if UsePropsedUKFSLAM

    %%  UKF Initializing

    Xtrue           = [xwp(1), ywp(1), thetawp(1)]';
    xLand           = lm(1,:);
    yLand           = lm(2,:);

    % - Initialize UKF estimate of state and covariance
    Xest            = Xtrue;
    Pest            = zeros(3, 3);
    alpha           = 1e-2;                                             % Spread of sigma points
    beta            = 2;                                                % Optimal for Gaussian distribution

    % - Define robot's state noise (assumed to be Gaussian and using wheel encoders)
    sigmaX          = 0.01;
    sigmaY          = 0.01;
    sigmaTheta      = 0.01;

    R               = diag([sigmaX^2 sigmaY^2 sigmaTheta^2]);

    % - Landmark measurement noise
    sigmaM          = 0.02;
    R_lm            = sigmaM*eye(2);

    % - Dynamic landmarks initialization
    vLM             = 6;

    % - Define time step
    dt              = 0.1;

    % - Moving landmark threshold
    % The noise in robot states measurement, represented by sigmaX, follows a
    % Rayleigh distribution, while the sensor's distance measuring noise is
    % Gaussian. 99% of these distributions fall withing the threshold region.
    threshold       = 4*sqrt((4-pi)/2)*sigmaX + 3*sigmaM;

    % - Create array to hold data
    GroundTruth     = cell(1,round(2.5 * p.Length / dt));
    EstTraject      = cell(1,round(2.5 * p.Length / dt));
    LandMark        = cell(1,round(2.5 * p.Length / dt));
    GroundTruth{1}  = Xtrue;
    EstTraject{1}   = Xest;
    dist            = zeros(nLM,1);
    counterStep     = 1;
    removeCheck     = '0';
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

    % - Moving landmarks path
    [centerGTx, centerGTy] = centroid(polyshape(xwp,ywp));
    distLand = [(xLand-centerGTx); (yLand-centerGTy)];
    [~,idxMove] = max(abs(distLand), [],1);

    selIdx = sub2ind(size(dirLM),1:size(lm,2),idxMove);
    dirLM(selIdx) = -sign(distLand(selIdx));

    %%  Main Loop

    for i = 1:length(xwp)

        dx_wp = xwp(i) - Xtrue(1);
        dy_wp = ywp(i) - Xtrue(2);

        while sqrt(dx_wp^2 + dy_wp^2) > 0.01

            Xest = Xest(1:3);
            Pest = Pest(1:3,1:3);
            R = R(1:3,1:3);
            idxP = zeros(nLM,1);

            % Update kappa and lambda for changing state size
            kappa = 3 - length(Xest);
            lambda = alpha^2 * (length(Xest) + kappa) - length(Xest);

            thetaGoal = atan2(dy_wp, dx_wp);
            omega = wrapToPi(thetaGoal - Xtrue(3)) / dt;
            v = min(sqrt(dx_wp^2 + dy_wp^2) / dt, 0.5);
            u = [v; omega];

            LandOld = [xLand; yLand];
            DistOld = dist;

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

            % Generate weights
            Wmean = [lambda / (length(Xest) + lambda); ones(2 * length(Xest), 1) * 1/(2*(length(Xest) + lambda))];
            Wcov = Wmean;
            Wcov(1) = Wcov(1) + (1 - alpha^2 + beta);

            % Generate Sigma Points
            L = covariance_root(Pest);
            Xsig = [Xest, Xest + sqrt(length(Xest) + lambda)*L, Xest - sqrt(length(Xest) + lambda)*L];

            % Propagate Sigma Points
            for j = 1:size(Xsig, 2)
                theta = Xsig(3, j);
                Xsig(1:3, j) = Xsig(1:3, j) + [u(1)*cos(theta)*dt; u(1)*sin(theta)*dt; u(2)*dt];
            end

            % Recover Mean and Covariance
            Xest = sum(Wmean' .* Xsig, 2);
            Pest = zeros(size(Pest));
            for j = 1:size(Xsig, 2)
                Pest = Pest + Wcov(j) * (Xsig(:, j) - Xest) * (Xsig(:, j) - Xest)';
            end
            Pest = Pest + R;

            % UKF Correction Step
            for j = 4:2:length(Xest)
                % Generate Sigma Points for landmark
                L = covariance_root(Pest);
                Xsig = [Xest, Xest + sqrt(length(Xest) + lambda)*L, Xest - sqrt(length(Xest) + lambda)*L];

                % Observation Model
                zSig = zeros(2, size(Xsig, 2));
                for k = 1:size(Xsig, 2)
                    dx = Xsig(j, k) - Xsig(1, k);
                    dy = Xsig(j+1, k) - Xsig(2, k);
                    zSig(:, k) = [sqrt(dx^2 + dy^2); atan2(dy, dx) - Xsig(3, k)];
                end

                zPred = sum(Wmean' .* zSig, 2);
                S = zeros(2, 2);
                for k = 1:size(zSig, 2)
                    S = S + Wcov(k) * (zSig(:, k) - zPred) * (zSig(:, k) - zPred)';
                end
                S = S + R_lm;

                % Cross covariance
                C = zeros(length(Xest), 2);
                for k = 1:size(Xsig, 2)
                    C = C + Wcov(k) * (Xsig(:, k) - Xest) * (zSig(:, k) - zPred)';
                end

                % Kalman Gain
                K = C / S;

                % Update state mean and covariance
                dx = Xest(j) - Xest(1);
                dy = Xest(j+1) - Xest(2);
                zTrue = [sqrt(dx^2 + dy^2); atan2(dy, dx) - Xest(3)];
                Xest = Xest + K * (zTrue - zPred);
                Pest = Pest - K * S * K';
            end

            % Simulate True System Dynamics
            Xtrue = Xtrue + [v*cos(Xtrue(3))*dt; v*sin(Xtrue(3))*dt; omega*dt] + diag(R(1:3,1:3));

            % Update GroundTruth and EstTraject
            counterStep = counterStep + 1;
            GroundTruth{counterStep} = Xtrue;
            EstTraject{counterStep} = Xest;

            dx_wp = xwp(i) - Xtrue(1);
            dy_wp = ywp(i) - Xtrue(2);

        end

    end

    %%  Post Processing

    % - Saving compelete data
    EstTrajectFull = EstTraject;

    % - Removing zeros (empty cells)
    GroundTruth(counterStep:end) = [];
    EstTraject(counterStep:end) = [];
    LandMark(counterStep:end) = [];

    % - Convert true histogram into a matrix
    GroundTruth     = [GroundTruth{1:end}];

    % - Select the first 3 elements of each cell (plot)
    EstTrajTrunc    = cellfun(@(x) x(1:min(3, end)), EstTraject, 'UniformOutput', false);

    % - Ensure that each cell contains exactly 3 elements (plot)
    EstTrajPadded   = cellfun(@(x) [x, zeros(1, max(0, 3-length(x)))], EstTrajTrunc, 'UniformOutput', false);

    % Convert the estimated histogram into a matrix (plot)
    EstTraject      = cell2mat(EstTrajPadded);


    %%  Results

    AIE = calculate_RMSE(EstTraject(1:2,:),GroundTruth(1:2,:));

    % - Plot Results
    if strcmpi(draw_option, 'on')

        ScreenSize  = get(0,'ScreenSize');
        figure('Name','Map','NumberTitle','off',...
            'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2.5)]);
        movegui('center')
        plot(GroundTruth(1, :), GroundTruth(2, :), 'LineWidth', 2), hold on
        plot(EstTraject(1, :), EstTraject(2, :), 'LineWidth', 2);

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

    if any(strcmp('rtm', varargin))
        figure('Name','Map','NumberTitle','off',...
            'Position', [0 0 600 600]);
        movegui('center')
        xlim([min(min(min(GroundTruth(1,:)),min(GroundTruth(2,:))),min(EstTraject(1,:)))-1.05*sDev ...
            max(max(max(GroundTruth(1,:)),max(GroundTruth(2,:))),max(EstTraject(2,:)))+1.05*sDev])
        ylim([min(min(min(GroundTruth(1,:)),min(GroundTruth(2,:))),min(EstTraject(1,:)))-1.05*sDev ...
            max(max(max(GroundTruth(1,:)),max(GroundTruth(2,:))),max(EstTraject(2,:)))+1.05*sDev])
        xlabel('$x$ (m)','Interpreter','latex','FontSize',12')
        ylabel('$y$ (m)','Interpreter','latex','FontSize',12')
        title('Real Time Map','Interpreter','latex','FontSize',14')
        grid on
        hold on;

        plot(GroundTruth(1,:), GroundTruth(2,:), 'LineWidth', 2)
        h = plot(EstTraject(1, 1), EstTraject(2, 1),'LineWidth',2,'color', '#D95319'); % Initial point

        triangleSize = 0.34;
        [xs, ys] = computeTriangle(EstTraject(:,1), EstTraject(:,1), triangleSize);
        triangle = patch(xs, ys, [0.8500 0.3250 0.0980]);

        for i = 2:counterStep-1
            LandmarkData = EstTrajectFull{i}(4:end);
            set(h, 'XData', EstTraject(1, 1:i), 'YData', EstTraject(2, 1:i));
            c = drawCircle(EstTraject(1,i),EstTraject(2,i),sDev,'FaceAlpha', 0.2, 'LineWidth', 1.5);

            hs = scatter(LandMark{i}(:,1), LandMark{i}(:,2),75,'+','markeredgecolor','r','LineWidth',2);


            [xs, ys] = computeTriangle(EstTraject(:,i-1), EstTraject(:,i), triangleSize);
            set(triangle, 'XData', xs, 'YData', ys);

            hl = zeros(size(LandmarkData,1)/2,1);
            if  ~isempty(LandmarkData)
                for j = 1:length(LandmarkData)/2
                    hl(j) = plot([EstTraject(1, i) LandmarkData(2*j-1)], [EstTraject(2, i) LandmarkData(2*j)],'LineWidth',1.5,'color','c');
                end

            end

            drawnow
            pause(dt)
            delete(hs)
            delete(c)
            if ~isempty(hl)
                delete(hl(1:j))
            end
        end

        drawCircle(EstTraject(1,end),EstTraject(2,end),sDev,'FaceAlpha', 0.2, 'LineWidth', 1.5);
        scatter(LandMark{end}(:,1), LandMark{end}(:,2),75,'+','markeredgecolor','r','LineWidth',2)
        for j = 1:length(LandmarkData)/2
            hl(j) = plot([EstTraject(1, end) LandmarkData(2*j-1)], [EstTraject(2, end) LandmarkData(2*j)],'LineWidth',1.5,'color','c');
        end
        hold off;
    end

end

%%  Static UKF-SLAM

if seedFlag == true
    rng(seed)
end

if UseConvUKFSLAM

    Xtrue           = [xwp(1), ywp(1), thetawp(1)]';
    xLand           = lm(1,:);
    yLand           = lm(2,:);

    % - Initialize UKF estimate of state and covariance
    Xest            = Xtrue;
    Pest            = zeros(3, 3);
    alpha           = 1e-2;                                             % Spread of sigma points
    beta            = 2;                                                % Optimal for Gaussian distribution

    % - Define robot's state noise (assumed to be Gaussian and using wheel encoders)
    sigmaX          = 0.01;
    sigmaY          = 0.01;
    sigmaTheta      = 0.01;

    R               = diag([sigmaX^2 sigmaY^2 sigmaTheta^2]);

    % - Landmark measurement noise
    sigmaM          = 0.02;
    R_lm            = sigmaM*eye(2);

    % - Dynamic landmarks initialization
    vLM             = 6;

    % - Define time step
    dt              = 0.1;

    % - Create array to hold data
    GroundTruth     = cell(1,round(2.5 * p.Length / dt));
    EstTraject      = cell(1,round(2.5 * p.Length / dt));
    LandMark        = cell(1,round(2.5 * p.Length / dt));
    GroundTruth{1}  = Xtrue;
    EstTraject{1}   = Xest;
    dist            = zeros(nLM,1);
    counterStep     = 1;
    removeCheck     = '0';
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

    % - Moving landmarks path
    [centerGTx, centerGTy] = centroid(polyshape(xwp,ywp));
    distLand = [(xLand-centerGTx); (yLand-centerGTy)];
    [~,idxMove] = max(abs(distLand), [],1);

    selIdx = sub2ind(size(dirLM),1:size(lm,2),idxMove);
    dirLM(selIdx) = -sign(distLand(selIdx));

    %%  Main Loop

    for i = 1:length(xwp)

        dx_wp = xwp(i) - Xtrue(1);
        dy_wp = ywp(i) - Xtrue(2);

        while sqrt(dx_wp^2 + dy_wp^2) > 0.01

            Xest = Xest(1:3);
            Pest = Pest(1:3,1:3);
            R = R(1:3,1:3);

            % Update kappa and lambda for changing state size
            kappa = 3 - length(Xest);
            lambda = alpha^2 * (length(Xest) + kappa) - length(Xest);

            thetaGoal = atan2(dy_wp, dx_wp);
            omega = wrapToPi(thetaGoal - Xtrue(3)) / dt;
            v = min(sqrt(dx_wp^2 + dy_wp^2) / dt, 0.5);
            u = [v; omega];

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

                % - Add the landmark to the state vector and covariance
                if dist(j) <= sDev
                    XestP = Xest;
                    Xest = [XestP; sigmaM*randn+xLand(j); sigmaM*randn+yLand(j)];
                    Pest = blkdiag(Pest, 0.01*eye(2));
                    R = blkdiag(R, R_lm);
                end
            end

            % Generate weights
            Wmean = [lambda / (length(Xest) + lambda); ones(2 * length(Xest), 1) * 1/(2*(length(Xest) + lambda))];
            Wcov = Wmean;
            Wcov(1) = Wcov(1) + (1 - alpha^2 + beta);

            % Generate Sigma Points
            L = covariance_root(Pest);
            Xsig = [Xest, Xest + sqrt(length(Xest) + lambda)*L, Xest - sqrt(length(Xest) + lambda)*L];

            % Propagate Sigma Points
            for j = 1:size(Xsig, 2)
                theta = Xsig(3, j);
                Xsig(1:3, j) = Xsig(1:3, j) + [u(1)*cos(theta)*dt; u(1)*sin(theta)*dt; u(2)*dt];
            end

            % Recover Mean and Covariance
            Xest = sum(Wmean' .* Xsig, 2);
            Pest = zeros(size(Pest));
            for j = 1:size(Xsig, 2)
                Pest = Pest + Wcov(j) * (Xsig(:, j) - Xest) * (Xsig(:, j) - Xest)';
            end
            Pest = Pest + R;

            % UKF Correction Step
            for j = 4:2:length(Xest)
                % Generate Sigma Points for landmark
                L = covariance_root(Pest);
                Xsig = [Xest, Xest + sqrt(length(Xest) + lambda)*L, Xest - sqrt(length(Xest) + lambda)*L];

                % Observation Model
                zSig = zeros(2, size(Xsig, 2));
                for k = 1:size(Xsig, 2)
                    dx = Xsig(j, k) - Xsig(1, k);
                    dy = Xsig(j+1, k) - Xsig(2, k);
                    zSig(:, k) = [sqrt(dx^2 + dy^2); atan2(dy, dx) - Xsig(3, k)];
                end

                zPred = sum(Wmean' .* zSig, 2);
                S = zeros(2, 2);
                for k = 1:size(zSig, 2)
                    S = S + Wcov(k) * (zSig(:, k) - zPred) * (zSig(:, k) - zPred)';
                end
                S = S + R_lm;

                % Cross covariance
                C = zeros(length(Xest), 2);
                for k = 1:size(Xsig, 2)
                    C = C + Wcov(k) * (Xsig(:, k) - Xest) * (zSig(:, k) - zPred)';
                end

                % Kalman Gain
                K = C / S;

                % Update state mean and covariance
                dx = Xest(j) - Xest(1);
                dy = Xest(j+1) - Xest(2);
                zTrue = [sqrt(dx^2 + dy^2); atan2(dy, dx) - Xest(3)];
                Xest = Xest + K * (zTrue - zPred);
                Pest = Pest - K * S * K';
            end

            % Simulate True System Dynamics
            Xtrue = Xtrue + [v*cos(Xtrue(3))*dt; v*sin(Xtrue(3))*dt; omega*dt] + diag(R(1:3,1:3));

            % Update GroundTruth and EstTraject
            counterStep = counterStep + 1;
            GroundTruth{counterStep} = Xtrue;
            EstTraject{counterStep} = Xest;

            dx_wp = xwp(i) - Xtrue(1);
            dy_wp = ywp(i) - Xtrue(2);

        end

    end

    %%  Post Processing

    % - Saving compelete data
    EstTrajectFull = EstTraject;

    % - Removing zeros (empty cells)
    GroundTruth(counterStep:end) = [];
    EstTraject(counterStep:end) = [];
    LandMark(counterStep:end) = [];


    % - Convert true histogram into a matrix
    GroundTruth     = [GroundTruth{1:end}];

    % - Select the first 3 elements of each cell (plot)
    EstTrajTrunc    = cellfun(@(x) x(1:min(3, end)), EstTraject, 'UniformOutput', false);

    % - Ensure that each cell contains exactly 3 elements (plot)
    EstTrajPadded   = cellfun(@(x) [x, zeros(1, max(0, 3-length(x)))], EstTrajTrunc, 'UniformOutput', false);

    % Convert the estimated histogram into a matrix (plot)
    EstTraject      = cell2mat(EstTrajPadded);


    %%  Results

    AIE = calculate_RMSE(EstTraject(1:2,:),GroundTruth(1:2,:));

    % - Plot Results
    if strcmpi(draw_option, 'on')

        ScreenSize  = get(0,'ScreenSize');
        figure('Name','Map','NumberTitle','off',...
            'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2.5)]);
        movegui('center')
        plot(GroundTruth(1, :), GroundTruth(2, :), 'LineWidth', 2), hold on
        plot(EstTraject(1, :), EstTraject(2, :), 'LineWidth', 2);

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

    if any(strcmp('rtm', varargin))
        figure('Name','Map','NumberTitle','off',...
            'Position', [0 0 600 600]);
        movegui('center')
        xlim([min(min(min(GroundTruth(1,:)),min(GroundTruth(2,:))),min(EstTraject(1,:)))-1.05*sDev ...
            max(max(max(GroundTruth(1,:)),max(GroundTruth(2,:))),max(EstTraject(2,:)))+1.05*sDev])
        ylim([min(min(min(GroundTruth(1,:)),min(GroundTruth(2,:))),min(EstTraject(1,:)))-1.05*sDev ...
            max(max(max(GroundTruth(1,:)),max(GroundTruth(2,:))),max(EstTraject(2,:)))+1.05*sDev])
        xlabel('$x$ (m)','Interpreter','latex','FontSize',12')
        ylabel('$y$ (m)','Interpreter','latex','FontSize',12')
        title('Real Time Map','Interpreter','latex','FontSize',14')
        grid on
        hold on;

        plot(GroundTruth(1,:), GroundTruth(2,:), 'LineWidth', 2)
        h = plot(EstTraject(1, 1), EstTraject(2, 1),'LineWidth',2,'color', '#D95319'); % Initial point

        triangleSize = 0.34;
        [xs, ys] = computeTriangle(EstTraject(:,1), EstTraject(:,1), triangleSize);
        triangle = patch(xs, ys, [0.8500 0.3250 0.0980]);

        for i = 2:counterStep-1
            LandmarkData = EstTrajectFull{i}(4:end);
            set(h, 'XData', EstTraject(1, 1:i), 'YData', EstTraject(2, 1:i));
            c = drawCircle(EstTraject(1,i),EstTraject(2,i),sDev,'FaceAlpha', 0.2, 'LineWidth', 1.5);

            hs = scatter(LandMark{i}(:,1), LandMark{i}(:,2),75,'+','markeredgecolor','r','LineWidth',2);


            [xs, ys] = computeTriangle(EstTraject(:,i-1), EstTraject(:,i), triangleSize);
            set(triangle, 'XData', xs, 'YData', ys);

            hl = zeros(size(LandmarkData,1)/2,1);
            if  ~isempty(LandmarkData)
                for j = 1:length(LandmarkData)/2
                    hl(j) = plot([EstTraject(1, i) LandmarkData(2*j-1)], [EstTraject(2, i) LandmarkData(2*j)],'LineWidth',1.5,'color','c');
                end

            end

            drawnow
            pause(dt)
            delete(hs)
            delete(c)
            if ~isempty(hl)
                delete(hl(1:j))
            end
        end

        drawCircle(EstTraject(1,end),EstTraject(2,end),sDev,'FaceAlpha', 0.2, 'LineWidth', 1.5);
        scatter(LandMark{end}(:,1), LandMark{end}(:,2),75,'+','markeredgecolor','r','LineWidth',2)
        for j = 1:length(LandmarkData)/2
            hl(j) = plot([EstTraject(1, end) LandmarkData(2*j-1)], [EstTraject(2, end) LandmarkData(2*j)],'LineWidth',1.5,'color','c');
        end
        hold off;
    end

end

end