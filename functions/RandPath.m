function p = RandPath(Min,Max)

N = 10;

try

    x = randi([Min,Max],[1 , N]);
    y = randi([Min,Max],[1 , N]);
    D = zeros(N,N);

    for i = 1:N-1

        for j = i+1:N

            D(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);

            D(j,i) = D(i,j);

        end

    end

    model.n = N;
    model.x = x;
    model.y = y;
    model.D = D;

    CostFunction = @(tour) TourLength(tour,model);

    nVar = model.n;

    %% ACO Parameters

    MaxIt = 100;      % Maximum Number of Iterations

    nAnt = 50;        % Number of Ants (Population Size)

    Q = 1;

    tau0 = 10*Q/(nVar*mean(model.D(:)));	% Initial Phromone

    alpha = 1;        % Phromone Exponential Weight
    beta = 1;         % Heuristic Exponential Weight

    rho = 0.05;       % Evaporation Rate


    %% Initialization

    eta = 1./model.D;             % Heuristic Information Matrix

    tau = tau0*ones(nVar,nVar);   % Phromone Matrix

    BestCost = zeros(MaxIt,1);    % Array to Hold Best Cost Values

    % Empty Ant
    empty_ant.Tour = [];
    empty_ant.Cost = [];

    % Ant Colony Matrix
    ant = repmat(empty_ant,nAnt,1);

    % Best Ant
    BestAnt.Cost = inf;


    %% ACO Main Loop

    for it = 1:MaxIt

        % Move Ants
        for k = 1:nAnt

            ant(k).Tour = randi([1 nVar]);

            for l = 2:nVar

                i = ant(k).Tour(end);

                P = tau(i,:).^alpha.*eta(i,:).^beta;

                P(ant(k).Tour) = 0;

                P = P/sum(P);

                j = RouletteWheelSelection(P);

                ant(k).Tour = [ant(k).Tour j];

            end

            ant(k).Cost = CostFunction(ant(k).Tour);

            if ant(k).Cost < BestAnt.Cost
                BestAnt = ant(k);
            end

        end

        % Update Phromones
        for k=1:nAnt

            tour = ant(k).Tour;
            tourP = tour;
            tour = [tourP tour(1)];

            for l = 1:nVar

                i = tour(l);
                j = tour(l+1);
                tau(i,j) = tau(i,j)+Q/ant(k).Cost;

            end

        end

        % Evaporation
        tau = (1-rho)*tau;

        % Plot Solution
        if (it >= 20 && abs(sum(BestCost(it-19:it)) - 20*BestCost(it)) <=1)
            break
        end

    end

    wp(1,:) = model.x(BestAnt.Tour);
    wp(2,:) = model.y(BestAnt.Tour);

    p.WayPoints = [wp wp(:,1)];
    p.Length = BestAnt.Cost;

catch ME
    if (strcmp(ME.identifier,'MATLAB:badsubscript'))
        p = RandPath(Min,Max);
    else
        rethrow(ME);
    end
end

end