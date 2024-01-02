clear;
clc;
close all;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Note that if bus index is not numbered in order, some matpower internal functions need to be called
%casename = 'case4gs';
casename = 'case14';
%casename = 'case24_ieee_rts';
basecase = loadcase(casename);
basecase.branch(:, TAP) = 0; % for static analysis; this is just slightly adjust BR_X
mpopt = mpoption('verbose', 0, 'out.all', 0, 'model', 'DC');
baseresult = runpf(basecase, mpopt);
bus = baseresult.bus;
branch = baseresult.branch;
gen = baseresult.gen;
[ref, pv, pq] = bustypes(bus, gen);

nline = length(baseresult.branch(:,1));
nbus = length(baseresult.bus(:, 1));
gbus = gen(:, GEN_BUS); 

%% This is manually calculated active power injection at each bus (Pg - Pd)
all_p = 0 - baseresult.bus(:, PD);
all_p(gbus) = baseresult.gen(:, PG) + all_p(gbus);

%% study about makeBdc; can be comment off 
% Bbus = -imag(Ybus)
[Ybus, Yf, Yt] = makeYbus(baseresult.baseMVA, baseresult.bus, baseresult.branch);
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseresult.baseMVA, baseresult.bus, baseresult.branch);
 
%% Calculate H matrix; assuming fully redundant z
deviceI = [branch(:, F_BUS); branch(:, T_BUS); bus(:, BUS_I)];
b = 1.0 ./ branch(:, BR_X);
branches = branch(:, [F_BUS, T_BUS]);
H = zeros(2 * nline + nbus, nbus);
for loop1 = 1 : length(deviceI)
    if loop1 <= nline
        H(loop1, branches(loop1, 1)) = b(loop1);
        H(loop1, branches(loop1, 2)) = -b(loop1);
    elseif (loop1 > nline ) & (loop1 <= 2*nline)
        H(loop1, branches(loop1-nline, 2)) = b(loop1-nline);
        H(loop1, branches(loop1-nline, 1)) = -b(loop1-nline);
    elseif loop1 > 2*nline
        curBus = bus(loop1 - 2*nline, BUS_I);
        %fneigh = find(mpccase.branch(:, F_BUS) == curBus);
        %H(loop1, :) = H(loop1, :) + sum(H(fneigh, :), 1);
        %tneigh = find(mpccase.branch(:, T_BUS) == curBus);
        %H(loop1, :) = H(loop1, :) + sum(H(tneigh+nline, :), 1);
        neigh = find(deviceI(1:2*nline) == curBus);
        H(loop1, :) = H(loop1, :) + sum(H(neigh, :), 1);
    end
end

%% use baseresult to generate the true x and true z
x_true = bus(:, VA) / 180 * pi; % change from degree to rad
% fully redundance z value is, PF, PT, PG; using this order, dont change
z_true = [branch(:, PF); branch(:, PT); all_p]/baseresult.baseMVA;

% verification z_true equal to H*x_true or not; 1e-15 means 0
% z_true - H*x_true

%% in case we want to less redudant measurement; folloing matpower interface

idx_zPF = transpose(1:nline);
%idx_zPT = transpose(1:nline);
idx_zPT = [];
idx_zPG = transpose(1:nbus);
%idx_zPG = [];

H_red = H([idx_zPF; (nline+idx_zPT) ; (nline+nline+idx_zPG)], :);
z_true_red = z_true([idx_zPF; (nline+idx_zPT) ; (nline+nline+idx_zPG)]);

% verify again to see the index are not messedup
% z_true_red - H_red * x_true

% further reduce H_red to remove the column corresponding to ref bus
H_red(:, ref) = [];
x_true_red = x_true;
x_true_red(ref) = [];

% verify again to see the index are not messedup
% z_true_red - H_red * x_true_red


%% following steps to generate noisy measurements
NOISE = 1;
sigma_PF = 0.001;
sigma_PT = [];
% sigma_PG = [];
sigma_PG = 0.001;

sigma_vector = [
    sigma_PF*ones(size(idx_zPF, 1), 1)
    sigma_PT*ones(size(idx_zPT, 1), 1)
    sigma_PG*ones(size(idx_zPG, 1), 1)
    ]; 
sigma_square = sigma_vector.^2;
R_inv = diag(1./sigma_square); % or W matrix

% z = z_true_red + NOISE * normrnd(0, sigma_vector);
z = z_true_red + NOISE * rand(size(sigma_vector)).* sigma_vector; % to remove the chance that random noise generate bad data 

%% perform weight least square error DC state estimation

gain = transpose(H_red) * R_inv * H_red;
x_est = inv(gain) * transpose(H_red) * R_inv * z;
z_est = H_red * x_est;

error_sqrsum = sum((z - z_est).^2./sigma_square)

%% bad data detection
alpha = 0.95;
freedom = length(z_est) - length(x_est);
baddata_t = chi2inv(alpha, freedom)

if error_sqrsum > baddata_t
    fprintf('bad data detected!');
end


