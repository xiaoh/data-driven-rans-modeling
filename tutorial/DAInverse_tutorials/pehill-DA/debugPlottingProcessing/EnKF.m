%% Debug EnKF

% full observation
debugDataPath = '../debugData/';

X_1 = load([debugDataPath, 'X_1.0.txt']);
HX_1 = X_1(1:end-4, :);

[NstateAug, Ns] = size(X_1);
Nstate = length(HX_1);

tr = ones(Nstate, 1);
ze = zeros(4, Nstate);

H = diag(tr); H = [H; ze];
Hpy = load([debugDataPath, 'H_1.0']);
deltaH = H - Hpy';

sigma = 1e-7;
R = sigma^2 * diag(tr);

% HX = H' * X_1;
% diffH = HX - HX_1;

Xmean_1Vec = mean(X_1, 2);
Xmean_1 = [];
for i = 1:Ns
    Xmean_1 = horzcat(Xmean_1, Xmean_1Vec);
end

P_1 = (1/(Ns-1))*(X_1 - Xmean_1)*(X_1 - Xmean_1)';

K_1 = P_1 * H / (H' * P_1 * H + R); 

K_1Py = load([debugDataPath, 'kalmanGain_1.0']);

diff = K_1 - K_1Py;
max(diff(:));

disp(['max delta H = ', max(deltaH(:))]);