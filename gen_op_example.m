clear;
close all;
clc;

%% add path of matpower if necessary
% addpath(genpath('/home/hugo/Dropbox/Experiments/matpower5.0'));


%mpccase = case5;
%casefile = 'case9';
casefile = 'case24_ieee_rts';
%casefile = 'case24';
%casefile = 'case30';
%casefile = 'case39';
%casefile = 'caseRTS96i';
%casefile = 'case118';
%casefile = case2746wop;
%casefile = 'case2737sop-1';
%casefile = 'case2737sop-3';

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

rng(23, 'simdTwister');


sample_size = 100;

for loop1 = 1 : sample_size
    % need to always load case here: so changes are not accumulated
    mpccase = loadcase(casefile);
    mpccase.branch(:, TAP) = 0; % simplify transformer set up
    nline = length(mpccase.branch(:,1));
    nbus = length(mpccase.bus(:, 1));

    % as an example change load randomly
    % in practice, we should use a benchmark 
    pd_change_ratio = rand(nbus, 1) * 0.2 - 0.1; % (-0.1, 0.1)
    mpccase.bus(:, PD) = mpccase.bus(:, PD)./(1 + pd_change_ratio);
    qd_change_ratio = rand(nbus, 1) * 0.2 - 0.1;
    mpccase.bus(:, QD) = mpccase.bus(:, QD)./(1 + pd_change_ratio);

    % perform an appopriate power flow analysis to obtain true measurements
    %  or measurements without noise
    % here DC power flow analysis as an example
    mpopt = mpoption('model', 'DC');
    baseresult = runpf(mpccase, mpopt);

end