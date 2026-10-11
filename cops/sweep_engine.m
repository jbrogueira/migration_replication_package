%%This script implements the DHA to compute the migration and welfare effects of the 2004 EU enlargement
%% Initialize environment
% Set environment
%% SWEEP-ENGINE HEADER EDIT: 'clear all / close all / clc' removed — this copy is
%% run inside the sweep_phi.m function workspace, which carries the injected
%% value-track state (INJECT_V) and the SWEEP_MODE flag; clearing would wipe them.
%% Everything below is an exact copy of COPS_shock_europe_skills.m apart from the
%% two marked SWEEP-ENGINE insertions.

endyear=2014;
tic

% Parameters
N = 17;
J = 2;
%% TIME_HORIZON env override (default 200). For T>200 the loaded baseline driver
%% series are padded with their stationary t=200 slice below (after the load).
th_env = getenv('TIME_HORIZON');
if isempty(th_env); time = 200; else; time = round(str2double(th_env)); end
beta = 0.97;
rho=4;

% Initial Conditions (relative time differences in values u_hat^(1/nu))
V_eu15_h = ones(J*N,time);
V_nms8_h = ones(J*N,time);
V_oth_h = ones(J*N,time);
V_eu15_h(:,1)=0;
V_nms8_h(:,1)=0;
V_oth_h(:,1)=0;
V_eu15_l = ones(J*N,time);
V_nms8_l = ones(J*N,time);
V_oth_l = ones(J*N,time);
V_eu15_l(:,1)=0;
V_nms8_l(:,1)=0;
V_oth_l(:,1)=0;

% Approach 2: split high-skill value functions at the target country into
% std and pref tracks. Used at the c-columns (states i, i2 = 16, 33);
% outside-c values are read from the unsplit V_*_h. Populated by the
% backward sweep in Step 4 (plan sec 2.5). Allocated here matching the
% existing V_*_h shape.
V_eu15_h_std  = ones(J*N,time);
V_eu15_h_pref = ones(J*N,time);
V_nms8_h_std  = ones(J*N,time);
V_nms8_h_pref = ones(J*N,time);
V_oth_h_std   = ones(J*N,time);
V_oth_h_pref  = ones(J*N,time);
V_eu15_h_std(:,1)=0;  V_eu15_h_pref(:,1)=0;
V_nms8_h_std(:,1)=0;  V_nms8_h_pref(:,1)=0;
V_oth_h_std(:,1)=0;   V_oth_h_pref(:,1)=0;

%%%%Loading Baseline 
load('Baseline_europe_Factual_skills.mat','series_wages_h','series_wages_l','series_mu_eu15_h','series_mu_nms8_h','series_mu_oth_h','series_mu_eu15_l','series_mu_nms8_l','series_mu_oth_l',...
    'L0_eu15_h','L0_nms8_h','L0_oth_h','L0_eu15_l','L0_nms8_l','L0_oth_l','series_Lnhat_h','series_Lnhat_l','series_Lnhat','series_pi','series_xbilat', 'series_sh_h', 'series_sh_l','series_sh_hat_l','series_sh_hat_h','series_Ldyn_nms8_h','series_Ldyn_nms8_l','series_Ldyn_eu15_h','series_Ldyn_eu15_l','series_Ldyn_oth_h','series_Ldyn_oth_l')

%% Horizon extension: pad loaded baseline driver series from 200 to `time` with
%% their terminal (t=200) slice. Justified only because the baseline is stationary
%% at t=200 (verified: wage hat flat to ~1e-6, stock drift ~1e-3/decade). Gated so
%% the default time=200 run is byte-identical (loop body never executes).
if time > 200
    pad_names = {'series_wages_h','series_wages_l', ...
        'series_mu_eu15_h','series_mu_nms8_h','series_mu_oth_h', ...
        'series_mu_eu15_l','series_mu_nms8_l','series_mu_oth_l', ...
        'series_Lnhat_h','series_Lnhat_l','series_Lnhat', ...
        'series_pi','series_xbilat','series_sh_h','series_sh_l', ...
        'series_sh_hat_l','series_sh_hat_h', ...
        'series_Ldyn_nms8_h','series_Ldyn_nms8_l', ...
        'series_Ldyn_eu15_h','series_Ldyn_eu15_l', ...
        'series_Ldyn_oth_h','series_Ldyn_oth_l'};
    for kk = 1:numel(pad_names)
        A = eval(pad_names{kk});
        nd = ndims(A); n = size(A, nd);
        if n < time
            if nd == 2
                A = [A, repmat(A(:, end), 1, time - n)];
            else
                A = cat(3, A, repmat(A(:, :, end), 1, 1, time - n));
            end
            eval([pad_names{kk} ' = A;']);
        end
    end
    clear A nd n kk pad_names;
    fprintf('[horizon] padded baseline series 200 -> time=%d (terminal-slice hold)\n', time);
end

%%distribution of households
L0_eu15_h=L0_eu15_h;
L0_eu15_l=L0_eu15_l;
L0_nms8_h=L0_nms8_h;
L0_nms8_l=L0_nms8_l;
L0_oth_h=L0_oth_h;
L0_oth_l=L0_oth_l;

Ln_hat0_h=series_Lnhat_h(1:N,:);
Ln_hat0_l=series_Lnhat_l(1:N,:);
Ln_hat0=series_Lnhat(1:N,:);
sh_hat_h0=series_sh_hat_h;
sh_hat_l0=series_sh_hat_l;

Din_baseline=series_pi;
wages0_h=series_wages_h;
wages0_l=series_wages_l;
%% CONVENTION FIX (2026-07-02, docs/cops_hat_convention_fix_plan.md): series_wages_*
%% are the baseline SEQUENTIAL wage changes w(t)/w(t-1) (entries converge to 1.0000),
%% not levels vs slot 1. The budget-clearing closure needs the baseline wage LEVEL
%% vs the base year (the HSV schedule's calibration anchor): cumulative products.
W_base_h_lev = cumprod(series_wages_h, 2);   % baseline wage level vs slot 1, per country
W_base_l_lev = cumprod(series_wages_l, 2);

%%%Changes in migration costs define as exp(change in taum)^(-1/nu)
load('data_xbilat_tariffs_tau_skills.mat', 'tau_mig_hat_eu15', 'tau_mig_hat_nms8');
tau_mig_hat_eu15_data=tau_mig_hat_eu15(:,:,2:end);
tau_mig_hat_nms8_data=tau_mig_hat_nms8(:,:,2:end);

tau_mig_hat_eu15=ones(34,34,time);
tau_mig_hat_nms8=ones(34,34,time);
tau_mig_hat_oth=ones(34,34,time);

%% Approach 2 + R2-light + temporary pref: parameters used downstream.
%% lambda_r2: per-period rate at which inelig outsiders regain eligibility.
%% k_temp_pref: number of years a worker stays in the pref pool at c after
%%   arrival; modelled as Poisson decay at rate 1/k_temp_pref per period
%%   (transition pref -> std at c).
%% A0b (with budget clearing): R2-light on (lambda_r2=0.1) and temp pref
%% on (k_temp_pref=10). With tau_P_pref=tau_h_base, the std pool gets
%% replenished by (i) inelig workers returning to c, (ii) pref workers
%% expiring to std. L_std stays bounded above zero. Tests whether budget
%% clearing converges cleanly when L_std is well-behaved.
lambda_r2    = 0.1;
k_temp_pref  = 10;
%% SWEEP-ENGINE ADDITION: K_TEMP_PREF env override (faster-expiry stepping
%% stone — smaller k caps the preferential pool earlier, shortens the
%% transition). Default unchanged (k=10).
ktp_env = getenv('K_TEMP_PREF');
if ~isempty(ktp_env); k_temp_pref = str2double(ktp_env); end
expire_rate  = 1 / k_temp_pref;

% tau_mig_hat_nms8(:,:,1:12)=tau_mig_hat_nms8_data;  % DISABLED: isolate preferential tax shock

shock_val = 1;  % Approach 2 disables the Approach 1 migration-cost shock; the
                % preferential rate enters via the after-tax wedge in Step 7.
%% TARGET_CTRY env override: target-country index in the 18-country ordering
%% (PT=16 default; ES=8, DE=5). Sector-2 (non-employment) index is i+17.
%% The native correction below assumes the target is an EU-15 member.
tc_env = getenv('TARGET_CTRY');
if isempty(tc_env); i = 16; else; i = round(str2double(tc_env)); end
i2 = i + 17;

% Skill-specific tau_mig_hat: create separate h-skill and l-skill matrices.
% High-skill: shock at the PT employment column (sector 1 = employed; the
% preferential tax is an income tax, so only moves into employment at the
% target are affected).
tau_mig_hat_eu15_h = tau_mig_hat_eu15;
tau_mig_hat_eu15_l = tau_mig_hat_eu15;
tau_mig_hat_nms8_h = tau_mig_hat_nms8;
tau_mig_hat_nms8_l = tau_mig_hat_nms8;
tau_mig_hat_oth_h  = tau_mig_hat_oth;
tau_mig_hat_oth_l  = tau_mig_hat_oth;

% One-time permanent shock at t=1. The DHA propagates it forward.
tau_mig_hat_eu15_h(:, i, 1) = shock_val;
tau_mig_hat_nms8_h(:, i, 1) = shock_val;
tau_mig_hat_oth_h(:, i, 1)  = shock_val;
% Low-skill: no shock (tau_mig_hat_*_l remain ones)

% Native correction: PT natives sit inside the EU-15 bloc and cannot be
% distinguished from other EU-15 nationals. For EU-15 workers already at PT
% (employed or non-employed), reset the cost of transitioning into PT
% employment so the shock does not lower the cost of remaining / becoming
% employed in PT for target-country natives. NMS-8 and ROW workers at PT
% are necessarily immigrants, so their cells are left shocked.
tau_mig_hat_eu15_h([i i2], i, 1) = 1;

%%%tarifffs
    tau_hat=ones(N+1,N+1,time);

    load('Baseline_europe_data_skills.mat', 'series_tau');
series_tau_data=series_tau;
    for t=1:endyear-2001
        series_tau(:,:,t)=series_tau_data(:,:,t);
    end
    for t=endyear-2001+1:time
        series_tau(:,:,t)=series_tau_data(:,:,endyear-2001);
    end

    tau_hat=ones(N+1,N+1,time);
    % tau_hat(:,:,2:end)=series_tau(:,:,2:end)./series_tau(:,:,1:end-1);  % DISABLED: isolate preferential tax shock
 

% Create mu_tilde for the first period (mu from the baseline economy)
mu_baseline_eu15_h=series_mu_eu15_h;
mu_baseline_nms8_h=series_mu_nms8_h;
mu_baseline_oth_h=series_mu_oth_h;
mu_baseline_eu15_l=series_mu_eu15_l;
mu_baseline_nms8_l=series_mu_nms8_l;
mu_baseline_oth_l=series_mu_oth_l;
mu00_eu15_h=mu_baseline_eu15_h(:,:,1); 
mu00_nms8_h=mu_baseline_nms8_h(:,:,1);
mu00_oth_h=mu_baseline_oth_h(:,:,1);
mu00_eu15_l=mu_baseline_eu15_l(:,:,1); 
mu00_nms8_l=mu_baseline_nms8_l(:,:,1);
mu00_oth_l=mu_baseline_oth_l(:,:,1);

% Target country indices in the 34-location vector (J*N = 34)
% Set above via TARGET_CTRY (default PT: sector 1 = 16, sector 2 = 16+17 = 33)

toldyn = 0.0001;
max_iter = 1E+20;
outer_max_iter = 500;  % outer V iter cap (diagnostic safety)
om_env = getenv('OUTER_MAX');
if ~isempty(om_env); outer_max_iter = round(str2double(om_env)); end
%% AGGLOM_A: agglomeration elasticity 'a' override (default 1.0 = committed).
%% Diagnostic: a=0 turns off the Ahat=(Ln_hat)^a productivity feedback to test
%% whether the early-slot eigenvalue>1 (outer Ymax floor) is the agglomeration loop.
aa_env = getenv('AGGLOM_A');
if isempty(aa_env); agglom_a = 1.0; else; agglom_a = str2double(aa_env); end
%% REGAIN_PRICED: re-eligibility pricing extension (docs/cops_regain_pricing_plan.md).
%% When 1, the ineligible (std) track's off-c destination columns price the regain
%% of eligibility abroad at the stock rate lambda_r2 (blend in the Yaux build below).
%% Default 0/unset = current engine: regain priced as immediate at off-c
%% continuations; the guarded blocks below do not execute.
rp_env = getenv('REGAIN_PRICED');
if isempty(rp_env); regain_priced = false; else; regain_priced = (str2double(rp_env) ~= 0); end
%% TAU_P_PREF: override the preferential rate (default 0.235). Set = tau_h_base
%% (0.24) for the neutral-policy audit: split fully operational but economically
%% inert, so re-aggregating the pools must reproduce unsplit COPS (identity test).
tpp_env = getenv('TAU_P_PREF');
if isempty(tpp_env); tau_P_pref_override = []; else; tau_P_pref_override = str2double(tpp_env); end
%% TAU_H_BASE / TAU_L_BASE: override the baseline average rates of the HSV
%% schedule (defaults 0.24 / 0.18). Unset = engine defaults bit-for-bit.
thb_env = getenv('TAU_H_BASE');
if isempty(thb_env); tau_h_base_override = []; else; tau_h_base_override = str2double(thb_env); end
tlb_env = getenv('TAU_L_BASE');
if isempty(tlb_env); tau_l_base_override = []; else; tau_l_base_override = str2double(tlb_env); end
Ymax = 1;
iter = 1;

%% Route 1 diagnostic (2026-05-28): per-outer-iter capture of the budget-block
%% amplifiers (closure denominator, lambda_hat, f_std, wage hats) to locate the
%% iter~26 spike under the corrected lambda-rescale closure. Read-only: writes
%% logs, does not alter the iteration. Full per-slot tables for iters in
%% [diag_iter_lo, diag_iter_hi]; a one-line summary every iter in route1_summary.log.
diag_iter_lo = 1;
diag_iter_hi = 10000;  % dump full per-slot table every iter; run hangs at the kink anyway
%% Config knobs read from the OS env (getenv survives the top-of-file `clear all`).
%% DAMP_SPLIT overrides the V_std/V_pref damping; DIAG_DIR sets the log directory.
ds_env = getenv('DAMP_SPLIT');
if isempty(ds_env); alpha_damp_split_cfg = 0.1; else; alpha_damp_split_cfg = str2double(ds_env); end
dd_env = getenv('DIAG_DIR');
if isempty(dd_env); diag_dir = 'iter_diag'; else; diag_dir = dd_env; end
%% LAB_DAMP: relaxation on the PT labor stocks across outer iters (1.0 = off =
%% committed behavior). Ldyn = f(V) is undamped today; it compounds over the
%% 200-period forward sim and runs away at mid-transition slots. Damp the stock
%% path: Ldyn_used = alpha_L*Ldyn_raw + (1-alpha_L)*Ldyn_used_prev. At a fixed
%% point raw=prev, so the equilibrium is unchanged (pure algorithm change).
ld_env = getenv('LAB_DAMP');
if isempty(ld_env); alpha_L = 1.0; else; alpha_L = str2double(ld_env); end
have_prev_L = false;
%% Reformulated closure (2026-05-28): feasibility clamp + optional lambda damping.
%% The raw closed-form lambda = numer/denom has the endogenous PT tax base in its
%% denominator, so a transient base dip sends lambda -> +/-inf and the wedge
%% detonates (the pole). CLAMP_TAU keeps tau_h_std, tau_l in [tau_floor, tau_max]
%% (default [0, 0.99]), bounding lambda and removing the pole. LAMBDA_DAMP relaxes
%% lambda across outer iters (clamp first, then damp). G stays fixed. At a feasible
%% fixed point the clamp does not bind and damping has settled, so lambda_applied =
%% lambda_raw and the budget clears exactly. Defaults reproduce committed behavior
%% (CLAMP_TAU off, LAMBDA_DAMP=1.0 off).
ct_env = getenv('CLAMP_TAU');
clamp_tau = ~isempty(ct_env) && (strcmp(ct_env,'1') || strcmpi(ct_env,'on') || strcmpi(ct_env,'true'));
tau_max_clamp   = 0.99;   % tau_h_std, tau_l <= this (no >100% / negative wedge)
tau_floor_clamp = 0.0;    % tau_h_std, tau_l >= this (no subsidy)
ll_env = getenv('LAMBDA_DAMP');
if isempty(ll_env); alpha_lambda = 1.0; else; alpha_lambda = str2double(ll_env); end
%% BUDGET_OFF: switch off the budget-clearing closure (B-iii rule). When set, the
%% standard-rate wedges are forced f_std=f_l=1 (standard h and l rates held at
%% baseline) and the period budget gap is left as an unrebated diagnostic; the
%% preferential cut (f_pref_const) stays on. This removes the lambda<->V Laffer
%% feedback that makes the lambda-rescale closure diverge.
bo_env = getenv('BUDGET_OFF');
budget_off = ~isempty(bo_env) && (strcmp(bo_env,'1') || strcmpi(bo_env,'on') || strcmpi(bo_env,'true'));
%% PREF_SCHED: operational-neutrality switch — the pref pool is taxed on the baseline
%% HSV schedule like everyone else (wedge = f_std, base inside A_hat, no flat-rate
%% revenue term). This is the exact eps=0 anchor of the closure-ON system under the
%% same-date baseline convention (PATH_SOLVER_SPEC.md sec 9). Default off.
ps_env = getenv('PREF_SCHED');
pref_sched = ~isempty(ps_env) && (strcmp(ps_env,'1') || strcmpi(ps_env,'on') || strcmpi(ps_env,'true'));
%% SWEEP-ENGINE ADDITION (continuation): PREF_MIX in [0,1] interpolates the
%% preferential pool's treatment between the schedule (0 = the exact neutral
%% anchor; equals PREF_SCHED) and the flat preferential rate (1 = policy).
%% Wedge and budget are interpolated consistently: a fraction (1-m) of the
%% pool's base rides the schedule (joins A_hat, scales with lambda), a
%% fraction m pays the flat rate (revenue term in B_hat).
pm_env = getenv('PREF_MIX');
if isempty(pm_env); pref_mix = 1.0; else; pref_mix = str2double(pm_env); end
if pref_sched; pref_mix = 0.0; end
%% Baseline-schedule treatment (author decision 2026-07-02): the HSV schedule is
%% INDEXED to the baseline wage path — its reference wage moves with the baseline
%% wage, so along the baseline the average rates stay at their calibrated values
%% (tau_h_base, tau_l_base) and rates respond only to the policy-induced wage
%% deviation What_*. The preferential advantage is then constant,
%% (1-tau_P)/(1-tau_h_base), matching the paper's static comparison. The unindexed
%% alternative (schedule fixed in wage units: average rates rise with the baseline
%% wage level, which roughly doubles over the data years in all 17 countries and
%% carries a nominal/dollar valuation component) is available for comparison via
%% TAU_INDEXED=0. Default: indexed.
ti_env = getenv('TAU_INDEXED');
tau_indexed = isempty(ti_env) || ~(strcmp(ti_env,'0') || strcmpi(ti_env,'off') || strcmpi(ti_env,'false'));
lambda_applied_prev = ones(1, time);  % persists across outer iters (NOT re-init in the loop)
if ~exist(diag_dir, 'dir'); mkdir(diag_dir); end
fprintf('[route1] alpha_damp_split=%.3f  alpha_L=%.3f  clamp_tau=%d  alpha_lambda=%.3f  budget_off=%d  time=%d  diag_dir=%s\n', ...
        alpha_damp_split_cfg, alpha_L, clamp_tau, alpha_lambda, budget_off, time, diag_dir);
route1_fid = fopen(fullfile(diag_dir, 'route1_summary.log'), 'w');
fprintf(route1_fid, '# Route 1 per-iter summary. Extremes across slots t=2..time-1.\n');
fprintf(route1_fid, '# denom = h_factor*A_hat + l_factor*C_hat (closure denominator; ->0 blows up lambda_hat).\n');
fprintf(route1_fid, '# Ytrk: track driving Ymax at Yarg slot (0=unsplit 1=std 2=pref).\n');
fprintf(route1_fid, '# n_clamp: # slots where the tau-clamp bound. max_resid: max|budget residual| (0=cleared).\n');
fprintf(route1_fid, '%5s %13s %13s %5s %13s %5s %13s %5s %12s %5s %12s %5s %5s %5s %6s %13s %5s\n', ...
    'iter','Ymax','min_denom','@t','max|lam|','@t','max_fstd','@t','min_tauh','@t','max|wh|','@t','Yarg','Ytrk','nclmp','max_resid','@t');

%% SWEEP-ENGINE INSERTION 1 of 2 (sweep_phi.m interface; this file is otherwise an
%% exact copy of COPS_shock_europe_skills.m — regenerate by cp + these two blocks).
%% Inject the value-track state supplied by sweep_phi in place of the default
%% all-ones start. Twelve tracks; slot-1/terminal conventions are the caller's.
if exist('INJECT_V', 'var') && isstruct(INJECT_V)
    V_eu15_h = INJECT_V.V_eu15_h;  V_nms8_h = INJECT_V.V_nms8_h;  V_oth_h = INJECT_V.V_oth_h;
    V_eu15_l = INJECT_V.V_eu15_l;  V_nms8_l = INJECT_V.V_nms8_l;  V_oth_l = INJECT_V.V_oth_l;
    V_eu15_h_std  = INJECT_V.V_eu15_h_std;   V_nms8_h_std  = INJECT_V.V_nms8_h_std;   V_oth_h_std  = INJECT_V.V_oth_h_std;
    V_eu15_h_pref = INJECT_V.V_eu15_h_pref;  V_nms8_h_pref = INJECT_V.V_nms8_h_pref;  V_oth_h_pref = INJECT_V.V_oth_h_pref;
end

while (iter <= outer_max_iter) && (Ymax > toldyn)
  
    for t = 1:time
        YY_eu15_h(:,:,t) = reshape(V_eu15_h(:, t), 1, J*N);
        Yaux0_eu15_h(:, :, t) = reshape(YY_eu15_h(:, :, t), J*N, 1);
    end
    Yaux_eu15_h = permute(Yaux0_eu15_h, [2, 1, 3]);
    Yaux1_eu15_h = zeros(J*N, J*N, time);
    for t = 1:time
        Yaux1_eu15_h(:, :, t) = repmat(Yaux_eu15_h(:, :, t), J*N, 1);
    end
    
       for t = 1:time
         YY_eu15_l(:,:,t) = reshape(V_eu15_l(:, t), 1, J*N);
        Yaux0_eu15_l(:, :, t) = reshape(YY_eu15_l(:, :, t), J*N, 1);
    end
    Yaux_eu15_l = permute(Yaux0_eu15_l, [2, 1, 3]);
    Yaux1_eu15_l = zeros(J*N, J*N, time);
    for t = 1:time
        Yaux1_eu15_l(:, :, t) = repmat(Yaux_eu15_l(:, :, t), J*N, 1);
    end
    
          for t = 1:time
        YY_nms8_h(:,:,t) = reshape(V_nms8_h(:, t), 1, J*N);
        Yaux0_nms8_h(:, :, t) = reshape(YY_nms8_h(:, :, t), J*N, 1);
    end
    Yaux_nms8_h = permute(Yaux0_nms8_h, [2, 1, 3]);
    Yaux1_nms8_h = zeros(J*N, J*N, time);
    for t = 1:time
        Yaux1_nms8_h(:, :, t) = repmat(Yaux_nms8_h(:, :, t), J*N, 1);
    end
    
       for t = 1:time
        YY_nms8_l(:,:,t) = reshape(V_nms8_l(:, t), 1, J*N);
        Yaux0_nms8_l(:, :, t) = reshape(YY_nms8_l(:, :, t), J*N, 1);
    end
    Yaux_nms8_l = permute(Yaux0_nms8_l, [2, 1, 3]);
    Yaux1_nms8_l = zeros(J*N, J*N, time);
    for t = 1:time
        Yaux1_nms8_l(:, :, t) = repmat(Yaux_nms8_l(:, :, t), J*N, 1);
    end
       
    for t = 1:time
       YY_oth_h(:,:,t) = reshape(V_oth_h(:, t), 1, J*N);
       Yaux0_oth_h(:, :, t) = reshape(YY_oth_h(:, :, t), J*N, 1);
    end
    Yaux_oth_h = permute(Yaux0_oth_h, [2, 1, 3]);
    Yaux1_oth_h = zeros(J*N, J*N, time);
    for t = 1:time
        Yaux1_oth_h(:, :, t) = repmat(Yaux_oth_h(:, :, t), J*N, 1);
    end
    
     for t = 1:time
        YY_oth_l(:,:,t) = reshape(V_oth_l(:, t), 1, J*N);
        Yaux0_oth_l(:, :, t) = reshape(YY_oth_l(:, :, t), J*N, 1);
    end
    Yaux_oth_l = permute(Yaux0_oth_l, [2, 1, 3]);
    Yaux1_oth_l = zeros(J*N, J*N, time);
    for t = 1:time
        Yaux1_oth_l(:, :, t) = repmat(Yaux_oth_l(:, :, t),J*N, 1);
    end

    %% Approach 2: build std/pref Y-matrices for high-skill blocs (plan sec 2.5).
    %% Y_outsider = Y_pref by construction (newcomers to c become preferential),
    %% so only two distinct vectors per bloc per period: Y_std and Y_pref.
    %% Outside-c columns inherit from Y_h; the c-columns ([i, i2] = [16, 33])
    %% are overridden with V_h_std / V_h_pref values at those c-states.
    %% Temporary pref policy: V_pref destination at c reflects expected residency.
    %% Each period, fraction expire_rate of pref workers at c transition to std.
    %% So V_eff_pref(c) = (1 - expire_rate) * V_pref(c) + expire_rate * V_std(c).
    %% V_std unchanged (std workers don't expire).
    Yaux1_eu15_h_std  = Yaux1_eu15_h;
    Yaux1_eu15_h_pref = Yaux1_eu15_h;
    Yaux1_nms8_h_std  = Yaux1_nms8_h;
    Yaux1_nms8_h_pref = Yaux1_nms8_h;
    Yaux1_oth_h_std   = Yaux1_oth_h;
    Yaux1_oth_h_pref  = Yaux1_oth_h;
    if regain_priced; offc_cols = setdiff(1:J*N, [i i2]); end
    for t = 1:time
        V_eff_eu15 = (1 - expire_rate) * V_eu15_h_pref([i i2], t) + expire_rate * V_eu15_h_std([i i2], t);
        V_eff_nms8 = (1 - expire_rate) * V_nms8_h_pref([i i2], t) + expire_rate * V_nms8_h_std([i i2], t);
        V_eff_oth  = (1 - expire_rate) * V_oth_h_pref([i i2], t)  + expire_rate * V_oth_h_std([i i2], t);
        %% Unsplit Yaux1_*_h: per plan sec 2.6 (mu^outsider = mu^pref under R3),
        %% outsiders crossing into c become pref. Override c-columns with V_eff
        %% so that V_h converges to a coherent outsider value (and matches the
        %% pref-track c-column override below). Without this, V_h carries a
        %% legacy V_h(c) at c-destinations that has no policy meaning and feeds
        %% noise into the convergence check via the unsplit Y_h backward sweep.
        Yaux1_eu15_h(:, [i i2], t) = repmat(V_eff_eu15', J*N, 1);
        Yaux1_nms8_h(:, [i i2], t) = repmat(V_eff_nms8', J*N, 1);
        Yaux1_oth_h(:, [i i2], t)  = repmat(V_eff_oth',  J*N, 1);
        Yaux1_eu15_h_std(:, [i i2], t)  = repmat(V_eu15_h_std([i i2], t)', J*N, 1);
        Yaux1_eu15_h_pref(:, [i i2], t) = repmat(V_eff_eu15', J*N, 1);
        Yaux1_nms8_h_std(:, [i i2], t)  = repmat(V_nms8_h_std([i i2], t)', J*N, 1);
        Yaux1_nms8_h_pref(:, [i i2], t) = repmat(V_eff_nms8', J*N, 1);
        Yaux1_oth_h_std(:, [i i2], t)   = repmat(V_oth_h_std([i i2], t)', J*N, 1);
        Yaux1_oth_h_pref(:, [i i2], t)  = repmat(V_eff_oth', J*N, 1);
        %% Re-eligibility pricing extension (docs/cops_regain_pricing_plan.md):
        %% price the regain of eligibility abroad at the stock rate lambda_r2 —
        %% the SAME rate the ineligible-stock dynamics use — instead of as
        %% immediate. Ineligible (std) track, off-c destination columns only:
        %%   Y_std(:, j', t) = lambda_r2 * Y_h(j', t) + (1 - lambda_r2) * V_std(j', t)
        %% for all j' not in [i i2]: the worker chooses destination j'; at the
        %% start of next period regain resolves with probability lambda_r2
        %% (continuation = eligible value) or not (continuation = own ineligible
        %% value). Arithmetic blend on the value hats, per period, mirroring the
        %% expiry blend V_eff at the c-columns. The c-columns keep the V_std
        %% override above (no regain at the adopting country). The unsplit
        %% Yaux1_*_h off-c columns read here are untouched by the c-column
        %% overrides above, so the blend reads the eligible-outsider values.
        if regain_priced
            Yaux1_eu15_h_std(:, offc_cols, t) = lambda_r2 * Yaux1_eu15_h(:, offc_cols, t) ...
                + (1 - lambda_r2) * repmat(V_eu15_h_std(offc_cols, t)', J*N, 1);
            Yaux1_nms8_h_std(:, offc_cols, t) = lambda_r2 * Yaux1_nms8_h(:, offc_cols, t) ...
                + (1 - lambda_r2) * repmat(V_nms8_h_std(offc_cols, t)', J*N, 1);
            Yaux1_oth_h_std(:, offc_cols, t)  = lambda_r2 * Yaux1_oth_h(:, offc_cols, t) ...
                + (1 - lambda_r2) * repmat(V_oth_h_std(offc_cols, t)', J*N, 1);
        end
    end

    %% Step 2 (create path of migration matrices)
    % Calculate mu1_tilde
    mu1_tilde_eu15_h = mu_baseline_eu15_h(:, :, 2).*(Yaux1_eu15_h(:, :, 2).^beta);
    mu_eu15_h = zeros(J*N, J*N, time); 
    mu_eu15_h(:, :, 1) = mu1_tilde_eu15_h;
    
    mu1_tilde_eu15_l = mu_baseline_eu15_l(:, :, 2).*(Yaux1_eu15_l(:, :, 2).^beta);
    mu_eu15_l = zeros(J*N, J*N, time); 
    mu_eu15_l(:, :, 1) = mu1_tilde_eu15_l;
    
    mu1_tilde_nms8_h = mu_baseline_nms8_h(:, :, 2).*(Yaux1_nms8_h(:, :, 2).^beta);
    mu_nms8_h = zeros(J*N, J*N, time); 
    mu_nms8_h(:, :, 1) = mu1_tilde_nms8_h;
    
    mu1_tilde_nms8_l = mu_baseline_nms8_l(:, :, 2).*(Yaux1_nms8_l(:, :, 2).^beta);
    mu_nms8_l = zeros(J*N, J*N, time); 
    mu_nms8_l(:, :, 1) = mu1_tilde_nms8_l;
    
    mu1_tilde_oth_h = mu_baseline_oth_h(:, :, 2).*(Yaux1_oth_h(:, :, 2).^beta);
    mu_oth_h = zeros(J*N, J*N, time); 
    mu_oth_h(:, :, 1) = mu1_tilde_oth_h;
    
    mu1_tilde_oth_l = mu_baseline_oth_l(:, :, 2).*(Yaux1_oth_l(:, :, 2).^beta);
    mu_oth_l = zeros(J*N, J*N, time);
    mu_oth_l(:, :, 1) = mu1_tilde_oth_l;

    %% Approach 2: split high-skill migration matrices into std and pref tracks.
    %% At slot 1 (pre-shock baseline): V_std = V_pref = V_h, so the std and pref
    %% tracks coincide with the unsplit migration matrix. They diverge from
    %% slot 2 onwards as V_*_h_std and V_*_h_pref pick up the after-tax wedge.
    mu1_tilde_eu15_h_std  = mu_baseline_eu15_h(:, :, 2).*(Yaux1_eu15_h_std(:, :, 2).^beta);
    mu1_tilde_eu15_h_pref = mu_baseline_eu15_h(:, :, 2).*(Yaux1_eu15_h_pref(:, :, 2).^beta);
    mu_eu15_h_std  = zeros(J*N, J*N, time);
    mu_eu15_h_pref = zeros(J*N, J*N, time);
    mu_eu15_h_std(:, :, 1)  = mu1_tilde_eu15_h_std;
    mu_eu15_h_pref(:, :, 1) = mu1_tilde_eu15_h_pref;

    mu1_tilde_nms8_h_std  = mu_baseline_nms8_h(:, :, 2).*(Yaux1_nms8_h_std(:, :, 2).^beta);
    mu1_tilde_nms8_h_pref = mu_baseline_nms8_h(:, :, 2).*(Yaux1_nms8_h_pref(:, :, 2).^beta);
    mu_nms8_h_std  = zeros(J*N, J*N, time);
    mu_nms8_h_pref = zeros(J*N, J*N, time);
    mu_nms8_h_std(:, :, 1)  = mu1_tilde_nms8_h_std;
    mu_nms8_h_pref(:, :, 1) = mu1_tilde_nms8_h_pref;

    mu1_tilde_oth_h_std  = mu_baseline_oth_h(:, :, 2).*(Yaux1_oth_h_std(:, :, 2).^beta);
    mu1_tilde_oth_h_pref = mu_baseline_oth_h(:, :, 2).*(Yaux1_oth_h_pref(:, :, 2).^beta);
    mu_oth_h_std  = zeros(J*N, J*N, time);
    mu_oth_h_pref = zeros(J*N, J*N, time);
    mu_oth_h_std(:, :, 1)  = mu1_tilde_oth_h_std;
    mu_oth_h_pref(:, :, 1) = mu1_tilde_oth_h_pref;

    % For t = 1
    num_eu15_h = tau_mig_hat_eu15_h(:,:,1).*mu_eu15_h(:, :, 1).*(Yaux1_eu15_h(:, :, 3).^beta);
    num_eu15_h(isnan(num_eu15_h)) = 0;
    den_eu15_h = sum(num_eu15_h')';
    den_eu15_h = den_eu15_h*ones(1, J*N); 
    mu_eu15_h(:, :, 2) = num_eu15_h./den_eu15_h;
    
    num_eu15_l = tau_mig_hat_eu15_l(:,:,1).*mu_eu15_l(:, :, 1).*(Yaux1_eu15_l(:, :, 3).^beta);
    num_eu15_l(isnan(num_eu15_l)) = 0;
    den_eu15_l = sum(num_eu15_l')';
    den_eu15_l = den_eu15_l*ones(1, J*N);  
    mu_eu15_l(:, :, 2) = num_eu15_l./den_eu15_l;
    
    
    num_nms8_h = tau_mig_hat_nms8_h(:,:,1).*mu_nms8_h(:, :, 1).*(Yaux1_nms8_h(:, :, 3).^beta);
    num_nms8_h(isnan(num_nms8_h)) = 0;
    den_nms8_h = sum(num_nms8_h')';
    den_nms8_h = den_nms8_h*ones(1, J*N); 
    mu_nms8_h(:, :, 2) = num_nms8_h./den_nms8_h;
    
     num_nms8_l = tau_mig_hat_nms8_l(:,:,1).*mu_nms8_l(:, :, 1).*(Yaux1_nms8_l(:, :, 3).^beta);
    num_nms8_l(isnan(num_nms8_l)) = 0;
    den_nms8_l = sum(num_nms8_l')';
    den_nms8_l = den_nms8_l*ones(1, J*N); 
    mu_nms8_l(:, :, 2) = num_nms8_l./den_nms8_l;
    
    num_oth_h = tau_mig_hat_oth_h(:,:,1).*mu_oth_h(:, :, 1).*(Yaux1_oth_h(:, :, 3).^beta);
    num_oth_h(isnan(num_oth_h)) = 0;
    den_oth_h = sum(num_oth_h')';
    den_oth_h = den_oth_h*ones(1, J*N); 
    mu_oth_h(:, :, 2) = num_oth_h./den_oth_h;
    
     num_oth_l = tau_mig_hat_oth_l(:,:,1).*mu_oth_l(:, :, 1).*(Yaux1_oth_l(:, :, 3).^beta);
    num_oth_l(isnan(num_oth_l)) = 0;
    den_oth_l = sum(num_oth_l')';
    den_oth_l = den_oth_l*ones(1, J*N);
    mu_oth_l(:, :, 2) = num_oth_l./den_oth_l;

    %% Approach 2: first-period std/pref migration matrix update for high-skill blocs.
    % EU-15
    num_eu15_h_std  = tau_mig_hat_eu15_h(:,:,1).*mu_eu15_h_std(:, :, 1).*(Yaux1_eu15_h_std(:, :, 3).^beta);
    num_eu15_h_std(isnan(num_eu15_h_std)) = 0;
    den_eu15_h_std = sum(num_eu15_h_std')';
    den_eu15_h_std = den_eu15_h_std*ones(1, J*N);
    mu_eu15_h_std(:, :, 2)  = num_eu15_h_std./den_eu15_h_std;

    num_eu15_h_pref = tau_mig_hat_eu15_h(:,:,1).*mu_eu15_h_pref(:, :, 1).*(Yaux1_eu15_h_pref(:, :, 3).^beta);
    num_eu15_h_pref(isnan(num_eu15_h_pref)) = 0;
    den_eu15_h_pref = sum(num_eu15_h_pref')';
    den_eu15_h_pref = den_eu15_h_pref*ones(1, J*N);
    mu_eu15_h_pref(:, :, 2) = num_eu15_h_pref./den_eu15_h_pref;

    % NMS-8
    num_nms8_h_std  = tau_mig_hat_nms8_h(:,:,1).*mu_nms8_h_std(:, :, 1).*(Yaux1_nms8_h_std(:, :, 3).^beta);
    num_nms8_h_std(isnan(num_nms8_h_std)) = 0;
    den_nms8_h_std = sum(num_nms8_h_std')';
    den_nms8_h_std = den_nms8_h_std*ones(1, J*N);
    mu_nms8_h_std(:, :, 2)  = num_nms8_h_std./den_nms8_h_std;

    num_nms8_h_pref = tau_mig_hat_nms8_h(:,:,1).*mu_nms8_h_pref(:, :, 1).*(Yaux1_nms8_h_pref(:, :, 3).^beta);
    num_nms8_h_pref(isnan(num_nms8_h_pref)) = 0;
    den_nms8_h_pref = sum(num_nms8_h_pref')';
    den_nms8_h_pref = den_nms8_h_pref*ones(1, J*N);
    mu_nms8_h_pref(:, :, 2) = num_nms8_h_pref./den_nms8_h_pref;

    % ROW (Other)
    num_oth_h_std  = tau_mig_hat_oth_h(:,:,1).*mu_oth_h_std(:, :, 1).*(Yaux1_oth_h_std(:, :, 3).^beta);
    num_oth_h_std(isnan(num_oth_h_std)) = 0;
    den_oth_h_std = sum(num_oth_h_std')';
    den_oth_h_std = den_oth_h_std*ones(1, J*N);
    mu_oth_h_std(:, :, 2)  = num_oth_h_std./den_oth_h_std;

    num_oth_h_pref = tau_mig_hat_oth_h(:,:,1).*mu_oth_h_pref(:, :, 1).*(Yaux1_oth_h_pref(:, :, 3).^beta);
    num_oth_h_pref(isnan(num_oth_h_pref)) = 0;
    den_oth_h_pref = sum(num_oth_h_pref')';
    den_oth_h_pref = den_oth_h_pref*ones(1, J*N);
    mu_oth_h_pref(:, :, 2) = num_oth_h_pref./den_oth_h_pref;

    % And now for all the other time periods
    for t = 2:time-2
        num = tau_mig_hat_eu15_h(:,:,t).*(mu_baseline_eu15_h(:, :, t+1)./mu_baseline_eu15_h(:, :, t)).*mu_eu15_h(:, :, t).*...
            (Yaux1_eu15_h(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')';
        den = den*ones(1, J*N);  
        mu_eu15_h(:, :, t+1) = num./den;
    end
    
     for t = 2:time-2
        num = tau_mig_hat_eu15_l(:,:,t).*(mu_baseline_eu15_l(:, :, t+1)./mu_baseline_eu15_l(:, :, t)).*mu_eu15_l(:, :, t).*...
            (Yaux1_eu15_l(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')';
        den = den*ones(1, J*N); 
        mu_eu15_l(:, :, t+1) = num./den;
    end
    
    for t = 2:time-2
        num = tau_mig_hat_nms8_h(:,:,t).*(mu_baseline_nms8_h(:, :, t+1)./mu_baseline_nms8_h(:, :, t)).*mu_nms8_h(:, :, t).*...
            (Yaux1_nms8_h(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')';
        den = den*ones(1, J*N); 
        mu_nms8_h(:, :, t+1) = num./den;
    end
    
    for t = 2:time-2
        num = tau_mig_hat_nms8_l(:,:,t).*(mu_baseline_nms8_l(:, :, t+1)./mu_baseline_nms8_l(:, :, t)).*mu_nms8_l(:, :, t).*...
            (Yaux1_nms8_l(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')';
        den = den*ones(1, J*N); 
        mu_nms8_l(:, :, t+1) = num./den;
    end
    
    
    for t = 2:time-2
        num = tau_mig_hat_oth_h(:,:,t).*(mu_baseline_oth_h(:, :, t+1)./mu_baseline_oth_h(:, :, t)).*mu_oth_h(:, :, t).*...
            (Yaux1_oth_h(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')';
        den = den*ones(1, J*N); 
        mu_oth_h(:, :, t+1) = num./den;
    end
    
     for t = 2:time-2
        num = tau_mig_hat_oth_l(:,:,t).*(mu_baseline_oth_l(:, :, t+1)./mu_baseline_oth_l(:, :, t)).*mu_oth_l(:, :, t).*...
            (Yaux1_oth_l(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')';
        den = den*ones(1, J*N);
        mu_oth_l(:, :, t+1) = num./den;
    end

    %% Approach 2: forward-period std/pref migration matrix update for high-skill blocs.
    for t = 2:time-2
        % EU-15 std
        num = tau_mig_hat_eu15_h(:,:,t).*(mu_baseline_eu15_h(:, :, t+1)./mu_baseline_eu15_h(:, :, t)).*mu_eu15_h_std(:, :, t).*...
            (Yaux1_eu15_h_std(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')'; den = den*ones(1, J*N);
        mu_eu15_h_std(:, :, t+1)  = num./den;
        % EU-15 pref
        num = tau_mig_hat_eu15_h(:,:,t).*(mu_baseline_eu15_h(:, :, t+1)./mu_baseline_eu15_h(:, :, t)).*mu_eu15_h_pref(:, :, t).*...
            (Yaux1_eu15_h_pref(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')'; den = den*ones(1, J*N);
        mu_eu15_h_pref(:, :, t+1) = num./den;

        % NMS-8 std
        num = tau_mig_hat_nms8_h(:,:,t).*(mu_baseline_nms8_h(:, :, t+1)./mu_baseline_nms8_h(:, :, t)).*mu_nms8_h_std(:, :, t).*...
            (Yaux1_nms8_h_std(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')'; den = den*ones(1, J*N);
        mu_nms8_h_std(:, :, t+1)  = num./den;
        % NMS-8 pref
        num = tau_mig_hat_nms8_h(:,:,t).*(mu_baseline_nms8_h(:, :, t+1)./mu_baseline_nms8_h(:, :, t)).*mu_nms8_h_pref(:, :, t).*...
            (Yaux1_nms8_h_pref(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')'; den = den*ones(1, J*N);
        mu_nms8_h_pref(:, :, t+1) = num./den;

        % ROW (Other) std
        num = tau_mig_hat_oth_h(:,:,t).*(mu_baseline_oth_h(:, :, t+1)./mu_baseline_oth_h(:, :, t)).*mu_oth_h_std(:, :, t).*...
            (Yaux1_oth_h_std(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')'; den = den*ones(1, J*N);
        mu_oth_h_std(:, :, t+1)  = num./den;
        % ROW (Other) pref
        num = tau_mig_hat_oth_h(:,:,t).*(mu_baseline_oth_h(:, :, t+1)./mu_baseline_oth_h(:, :, t)).*mu_oth_h_pref(:, :, t).*...
            (Yaux1_oth_h_pref(:, :, t+2).^beta);
        num(isnan(num)) = 0;
        den = sum(num')'; den = den*ones(1, J*N);
        mu_oth_h_pref(:, :, t+1) = num./den;
    end

    %% Step 3 (create path of labor)
    L00 = reshape(L0_eu15_h, 1, J*N); 
    L00_aux0 = reshape(L00, J*N, 1);
    L00_aux1 = repmat(L00_aux0, 1, J*N);
    L00_aux1 = mu00_eu15_h.*L00_aux1;
    L1 = sum(L00_aux1)';
    L1_eu15_h = reshape(L1, 1, J*N);
    
    L00 = reshape(L0_eu15_l, 1, J*N); 
    L00_aux0 = reshape(L00, J*N, 1);
    L00_aux1 = repmat(L00_aux0, 1, J*N);
    L00_aux1 = mu00_eu15_l.*L00_aux1;
    L1 = sum(L00_aux1)';
    L1_eu15_l = reshape(L1, 1, J*N);
    
    L00 = reshape(L0_nms8_h, 1, J*N); 
    L00_aux0 = reshape(L00, J*N, 1);
    L00_aux1 = repmat(L00_aux0, 1, J*N);
    L00_aux1 = mu00_nms8_h.*L00_aux1;
    L1 = sum(L00_aux1)';
    L1_nms8_h = reshape(L1, 1, J*N);
    
    L00 = reshape(L0_nms8_l, 1, J*N); 
    L00_aux0 = reshape(L00, J*N, 1);
    L00_aux1 = repmat(L00_aux0, 1, J*N);
    L00_aux1 = mu00_nms8_l.*L00_aux1;
    L1 = sum(L00_aux1)';
    L1_nms8_l = reshape(L1, 1, J*N);
    
    L00 = reshape(L0_oth_h, 1, J*N); 
    L00_aux0 = reshape(L00, J*N, 1);
    L00_aux1 = repmat(L00_aux0, 1, J*N);
    L00_aux1 = mu00_oth_h.*L00_aux1;
    L1 = sum(L00_aux1)';
    L1_oth_h = reshape(L1, 1, J*N);
    
    L00 = reshape(L0_oth_l, 1, J*N);
    L00_aux0 = reshape(L00, J*N, 1);
    L00_aux1 = repmat(L00_aux0, 1, J*N);
    L00_aux1 = mu00_oth_l.*L00_aux1;
    L1 = sum(L00_aux1)';
    L1_oth_l = reshape(L1, 1, J*N);

    %% Approach 2: split L1 into std and pref pools at PT for high-skill blocs (plan sec 2.7).
    %% At slot 1 (model t=0): L_std = L0 at c-states (all PT residents are std); L_pref = 0.
    %% At slot 2 (model t=1): L_std = forward via mu00 from c-state sources; L_pref = forward
    %% from outside-c sources (newcomers, who become pref under the R3 rule). The transition
    %% uses the baseline mu00; the policy first affects mu(:,:,2) and so reaches Ldyn_*_h at
    %% slot 3 onwards.
    %% EU-15
    L00_eu15_h_src_std = zeros(J*N, 1);  L00_eu15_h_src_std([i i2])  = reshape(L0_eu15_h, J*N, 1)([i i2]);
    L00_eu15_h_src_pref = reshape(L0_eu15_h, J*N, 1);  L00_eu15_h_src_pref([i i2]) = 0;
    L1_eu15_h_std  = sum(mu00_eu15_h .* repmat(L00_eu15_h_src_std,  1, J*N))';
    L1_eu15_h_pref = sum(mu00_eu15_h .* repmat(L00_eu15_h_src_pref, 1, J*N))';
    L1_eu15_h_std  = reshape(L1_eu15_h_std,  1, J*N);
    L1_eu15_h_pref = reshape(L1_eu15_h_pref, 1, J*N);
    %% NMS-8
    L00_nms8_h_src_std = zeros(J*N, 1);  L00_nms8_h_src_std([i i2])  = reshape(L0_nms8_h, J*N, 1)([i i2]);
    L00_nms8_h_src_pref = reshape(L0_nms8_h, J*N, 1);  L00_nms8_h_src_pref([i i2]) = 0;
    L1_nms8_h_std  = sum(mu00_nms8_h .* repmat(L00_nms8_h_src_std,  1, J*N))';
    L1_nms8_h_pref = sum(mu00_nms8_h .* repmat(L00_nms8_h_src_pref, 1, J*N))';
    L1_nms8_h_std  = reshape(L1_nms8_h_std,  1, J*N);
    L1_nms8_h_pref = reshape(L1_nms8_h_pref, 1, J*N);
    %% ROW (oth)
    L00_oth_h_src_std = zeros(J*N, 1);  L00_oth_h_src_std([i i2])  = reshape(L0_oth_h, J*N, 1)([i i2]);
    L00_oth_h_src_pref = reshape(L0_oth_h, J*N, 1);  L00_oth_h_src_pref([i i2]) = 0;
    L1_oth_h_std  = sum(mu00_oth_h .* repmat(L00_oth_h_src_std,  1, J*N))';
    L1_oth_h_pref = sum(mu00_oth_h .* repmat(L00_oth_h_src_pref, 1, J*N))';
    L1_oth_h_std  = reshape(L1_oth_h_std,  1, J*N);
    L1_oth_h_pref = reshape(L1_oth_h_pref, 1, J*N);

    Ldyn_eu15_h = zeros(1, J*N, time);
    Ldyn_eu15_h(:, :, 2) = L1_eu15_h;
    Ldyn_eu15_h(:, :, 1) = reshape(L0_eu15_h, 1, J*N);
    
    Ldyn_eu15_l = zeros(1, J*N, time); 
    Ldyn_eu15_l(:, :, 2) = L1_eu15_l;
    Ldyn_eu15_l(:, :, 1) = reshape(L0_eu15_l, 1, J*N);
    
    Ldyn_nms8_h = zeros(1, J*N, time); 
    Ldyn_nms8_h(:, :, 2) = L1_nms8_h;
    Ldyn_nms8_h(:, :, 1) = reshape(L0_nms8_h, 1, J*N);
    
    Ldyn_nms8_l = zeros(1, J*N, time); 
    Ldyn_nms8_l(:, :, 2) = L1_nms8_l;
    Ldyn_nms8_l(:, :, 1) = reshape(L0_nms8_l, 1, J*N);
    
    Ldyn_oth_h = zeros(1, J*N, time); 
    Ldyn_oth_h(:, :, 2) = L1_oth_h;
    Ldyn_oth_h(:, :, 1) = reshape(L0_oth_h, 1, J*N);
    
    Ldyn_oth_l = zeros(1, J*N, time);
    Ldyn_oth_l(:, :, 2) = L1_oth_l;
    Ldyn_oth_l(:, :, 1) = reshape(L0_oth_l, 1, J*N);

    %% ============================================================
    %% Approach 2: split high-skill at target country c into std/pref pools.
    %% At array slot 1 (model t = 0, pre-shock): the entire baseline h-skill
    %% stock at PT is in the std pool; the pref pool is empty by construction.
    %% Slot 2+ is populated by the std/pref population update (plan sec 2.7,
    %% implemented in Step 5). Outside c (states != [i, i2]) the std/pref
    %% arrays remain zero; the unsplit Ldyn_*_h arrays carry outside-c stocks.
    %% ============================================================
    c_states = [i, i2];  % PT employment (16) and non-employment (33)

    Ldyn_eu15_h_std  = zeros(1, J*N, time);
    Ldyn_eu15_h_pref = zeros(1, J*N, time);
    Ldyn_eu15_h_std(:, c_states, 1)  = Ldyn_eu15_h(:, c_states, 1);
    Ldyn_eu15_h_std(:, c_states, 2)  = L1_eu15_h_std(:, c_states);
    Ldyn_eu15_h_pref(:, c_states, 2) = L1_eu15_h_pref(:, c_states);

    Ldyn_nms8_h_std  = zeros(1, J*N, time);
    Ldyn_nms8_h_pref = zeros(1, J*N, time);
    Ldyn_nms8_h_std(:, c_states, 1)  = Ldyn_nms8_h(:, c_states, 1);
    Ldyn_nms8_h_std(:, c_states, 2)  = L1_nms8_h_std(:, c_states);
    Ldyn_nms8_h_pref(:, c_states, 2) = L1_nms8_h_pref(:, c_states);

    Ldyn_oth_h_std  = zeros(1, J*N, time);
    Ldyn_oth_h_pref = zeros(1, J*N, time);
    Ldyn_oth_h_std(:, c_states, 1)  = Ldyn_oth_h(:, c_states, 1);
    Ldyn_oth_h_std(:, c_states, 2)  = L1_oth_h_std(:, c_states);
    Ldyn_oth_h_pref(:, c_states, 2) = L1_oth_h_pref(:, c_states);

    %% R2-light: ineligible-outside stocks (plan sec 2.13). Workers leaving
    %% c on the std rate become ineligible while abroad; eligibility decays
    %% geometrically at rate lambda_r2 per period; on re-entering c they
    %% re-join the std pool. Slot 1: no inelig stock (everyone abroad is elig).
    % lambda_r2 and k_temp_pref defined at top of script.
    Ldyn_eu15_h_inelig = zeros(1, J*N, time);
    Ldyn_nms8_h_inelig = zeros(1, J*N, time);
    Ldyn_oth_h_inelig  = zeros(1, J*N, time);

    %% Approach 2 + R2-light: high-skill forward update (plan sec 2.7, 2.13).
    %% Four sources: c-std, c-pref, elig-outside, inelig-outside.
    %% mu_*_h_std applies to c-std and inelig-outside (both see V_std at c).
    %% mu_*_h_pref applies to c-pref and elig-outside (both see V_pref at c).
    %% Std pool at c gets own stayers + inelig re-entries.
    %% Pref pool at c gets own stayers + elig arrivals.
    %% Inelig at outside is replenished by std-pool leavers and itself, then
    %% decays geometrically at rate lambda_r2 per period.
    min_Lelig_sweep = Inf;  % H3-3 diagnostic: min eligible-outside mass this sweep
    for t = 2:time-1
        L_t_std     = zeros(J*N, 1);  L_t_std([i i2])     = Ldyn_eu15_h_std(:,  [i i2], t)';
        L_t_pref    = zeros(J*N, 1);  L_t_pref([i i2])    = Ldyn_eu15_h_pref(:, [i i2], t)';
        L_t_total_outside = reshape(Ldyn_eu15_h(:, :, t), J*N, 1);  L_t_total_outside([i i2]) = 0;
        L_t_inelig  = reshape(Ldyn_eu15_h_inelig(:, :, t), J*N, 1);
        L_t_elig    = L_t_total_outside - L_t_inelig;
        min_Lelig_sweep = min(min_Lelig_sweep, min(L_t_elig));  % H3-3
        flow_std     = sum(mu_eu15_h_std(:,  :, t) .* repmat(L_t_std,    1, J*N))';
        flow_pref    = sum(mu_eu15_h_pref(:, :, t) .* repmat(L_t_pref,   1, J*N))';
        flow_elig    = sum(mu_eu15_h_pref(:, :, t) .* repmat(L_t_elig,   1, J*N))';
        flow_inelig  = sum(mu_eu15_h_std(:,  :, t) .* repmat(L_t_inelig, 1, J*N))';
        Ldyn_eu15_h_std(:,  [i i2], t+1) = flow_std([i i2])'  + flow_inelig([i i2])';
        Ldyn_eu15_h_pref(:, [i i2], t+1) = flow_pref([i i2])' + flow_elig([i i2])';
        % Temporary pref: fraction expire_rate of pref pool at c transitions to std.
        expire_eu15 = expire_rate * Ldyn_eu15_h_pref(:, [i i2], t+1);
        Ldyn_eu15_h_std(:,  [i i2], t+1) = Ldyn_eu15_h_std(:,  [i i2], t+1) + expire_eu15;
        Ldyn_eu15_h_pref(:, [i i2], t+1) = Ldyn_eu15_h_pref(:, [i i2], t+1) - expire_eu15;
        % Inelig at outside-c: new std-leavers enter at full mass; pre-existing
        % inelig stock decays by (1-lambda_r2). Earlier code applied (1-lambda)
        % to both, mis-classifying ~lambda% of fresh std-leavers as elig each
        % period — a R2-light timing bug that biased pref/std pool composition.
        flow_std_offc    = flow_std;     flow_std_offc([i i2])    = 0;
        flow_inelig_offc = flow_inelig;  flow_inelig_offc([i i2]) = 0;
        Ldyn_eu15_h_inelig(:, :, t+1) = (flow_std_offc + (1 - lambda_r2) * flow_inelig_offc)';
        Ldyn_eu15_h(:, :, t+1) = (flow_std + flow_pref + flow_elig + flow_inelig)';
    end
    Ldyn_eu15_h(:, :, time) = 0;
    
     for t = 2:time-1
        aux0 = reshape(Ldyn_eu15_l(:, :, t), J*N,1); 
        aux1 = repmat(aux0, 1, J*N);
        aux1 = mu_eu15_l(:, :, t).*aux1;
        aux2 = sum(aux1)';
       Ldyn_eu15_l(:, :, t+1) = reshape(aux2, 1, J*N);
    end
    Ldyn_eu15_l(:, :, time) = 0;
    
    
    for t = 2:time-1
        L_t_std     = zeros(J*N, 1);  L_t_std([i i2])     = Ldyn_nms8_h_std(:,  [i i2], t)';
        L_t_pref    = zeros(J*N, 1);  L_t_pref([i i2])    = Ldyn_nms8_h_pref(:, [i i2], t)';
        L_t_total_outside = reshape(Ldyn_nms8_h(:, :, t), J*N, 1);  L_t_total_outside([i i2]) = 0;
        L_t_inelig  = reshape(Ldyn_nms8_h_inelig(:, :, t), J*N, 1);
        L_t_elig    = L_t_total_outside - L_t_inelig;
        min_Lelig_sweep = min(min_Lelig_sweep, min(L_t_elig));  % H3-3
        flow_std     = sum(mu_nms8_h_std(:,  :, t) .* repmat(L_t_std,    1, J*N))';
        flow_pref    = sum(mu_nms8_h_pref(:, :, t) .* repmat(L_t_pref,   1, J*N))';
        flow_elig    = sum(mu_nms8_h_pref(:, :, t) .* repmat(L_t_elig,   1, J*N))';
        flow_inelig  = sum(mu_nms8_h_std(:,  :, t) .* repmat(L_t_inelig, 1, J*N))';
        Ldyn_nms8_h_std(:,  [i i2], t+1) = flow_std([i i2])'  + flow_inelig([i i2])';
        Ldyn_nms8_h_pref(:, [i i2], t+1) = flow_pref([i i2])' + flow_elig([i i2])';
        expire_nms8 = expire_rate * Ldyn_nms8_h_pref(:, [i i2], t+1);
        Ldyn_nms8_h_std(:,  [i i2], t+1) = Ldyn_nms8_h_std(:,  [i i2], t+1) + expire_nms8;
        Ldyn_nms8_h_pref(:, [i i2], t+1) = Ldyn_nms8_h_pref(:, [i i2], t+1) - expire_nms8;
        % Inelig decay: fresh std-leavers enter full mass; pre-existing inelig
        % decays by (1-lambda_r2). See note in EU-15 block.
        flow_std_offc    = flow_std;     flow_std_offc([i i2])    = 0;
        flow_inelig_offc = flow_inelig;  flow_inelig_offc([i i2]) = 0;
        Ldyn_nms8_h_inelig(:, :, t+1) = (flow_std_offc + (1 - lambda_r2) * flow_inelig_offc)';
        Ldyn_nms8_h(:, :, t+1) = (flow_std + flow_pref + flow_elig + flow_inelig)';
    end
    Ldyn_nms8_h(:, :, time) = 0;
    
     for t = 2:time-1
        aux0 = reshape(Ldyn_nms8_l(:, :, t), J*N,1); 
        aux1 = repmat(aux0, 1, J*N);
        aux1 = mu_nms8_l(:, :, t).*aux1;
        aux2 = sum(aux1)';
        Ldyn_nms8_l(:, :, t+1) = reshape(aux2, 1, J*N);
    end
    Ldyn_nms8_l(:, :, time) = 0;
    
    
    for t = 2:time-1
        L_t_std     = zeros(J*N, 1);  L_t_std([i i2])     = Ldyn_oth_h_std(:,  [i i2], t)';
        L_t_pref    = zeros(J*N, 1);  L_t_pref([i i2])    = Ldyn_oth_h_pref(:, [i i2], t)';
        L_t_total_outside = reshape(Ldyn_oth_h(:, :, t), J*N, 1);  L_t_total_outside([i i2]) = 0;
        L_t_inelig  = reshape(Ldyn_oth_h_inelig(:, :, t), J*N, 1);
        L_t_elig    = L_t_total_outside - L_t_inelig;
        min_Lelig_sweep = min(min_Lelig_sweep, min(L_t_elig));  % H3-3
        flow_std     = sum(mu_oth_h_std(:,  :, t) .* repmat(L_t_std,    1, J*N))';
        flow_pref    = sum(mu_oth_h_pref(:, :, t) .* repmat(L_t_pref,   1, J*N))';
        flow_elig    = sum(mu_oth_h_pref(:, :, t) .* repmat(L_t_elig,   1, J*N))';
        flow_inelig  = sum(mu_oth_h_std(:,  :, t) .* repmat(L_t_inelig, 1, J*N))';
        Ldyn_oth_h_std(:,  [i i2], t+1) = flow_std([i i2])'  + flow_inelig([i i2])';
        Ldyn_oth_h_pref(:, [i i2], t+1) = flow_pref([i i2])' + flow_elig([i i2])';
        expire_oth = expire_rate * Ldyn_oth_h_pref(:, [i i2], t+1);
        Ldyn_oth_h_std(:,  [i i2], t+1) = Ldyn_oth_h_std(:,  [i i2], t+1) + expire_oth;
        Ldyn_oth_h_pref(:, [i i2], t+1) = Ldyn_oth_h_pref(:, [i i2], t+1) - expire_oth;
        % Inelig at outside-c: new std-leavers enter at full mass; pre-existing
        % inelig stock decays by (1-lambda_r2). Earlier code applied (1-lambda)
        % to both, mis-classifying ~lambda% of fresh std-leavers as elig.
        flow_std_offc    = flow_std;     flow_std_offc([i i2])    = 0;
        flow_inelig_offc = flow_inelig;  flow_inelig_offc([i i2]) = 0;
        Ldyn_oth_h_inelig(:, :, t+1) = (flow_std_offc + (1 - lambda_r2) * flow_inelig_offc)';
        Ldyn_oth_h(:, :, t+1) = (flow_std + flow_pref + flow_elig + flow_inelig)';
    end
    Ldyn_oth_h(:, :, time) = 0;
    
    for t = 2:time-1
         aux0 = reshape(Ldyn_oth_l(:, :, t), J*N,1); 
        aux1 = repmat(aux0, 1, J*N);
        aux1 = mu_oth_l(:, :, t).*aux1;
        aux2 = sum(aux1)';
        Ldyn_oth_l(:, :, t+1) = reshape(aux2, 1, J*N);
    end
    Ldyn_oth_l(:, :, time) = 0;

    %% Route 1 follow-up (2026-05-28): relaxation on the PT labor stocks.
    %% Blend this iter's forward-simulated stock path with the previous iter's
    %% (uniform alpha_L over all bloc primitives preserves every linear relation
    %% among them; the combined sums below inherit the damping). Next iter's
    %% forward sim rebuilds from L0, so this does not contaminate the simulation.
    if (alpha_L < 1.0) && have_prev_L
        Ldyn_eu15_h = alpha_L*Ldyn_eu15_h + (1-alpha_L)*Ldyn_eu15_h_prev;
        Ldyn_nms8_h = alpha_L*Ldyn_nms8_h + (1-alpha_L)*Ldyn_nms8_h_prev;
        Ldyn_oth_h  = alpha_L*Ldyn_oth_h  + (1-alpha_L)*Ldyn_oth_h_prev;
        Ldyn_eu15_l = alpha_L*Ldyn_eu15_l + (1-alpha_L)*Ldyn_eu15_l_prev;
        Ldyn_nms8_l = alpha_L*Ldyn_nms8_l + (1-alpha_L)*Ldyn_nms8_l_prev;
        Ldyn_oth_l  = alpha_L*Ldyn_oth_l  + (1-alpha_L)*Ldyn_oth_l_prev;
        Ldyn_eu15_h_std  = alpha_L*Ldyn_eu15_h_std  + (1-alpha_L)*Ldyn_eu15_h_std_prev;
        Ldyn_eu15_h_pref = alpha_L*Ldyn_eu15_h_pref + (1-alpha_L)*Ldyn_eu15_h_pref_prev;
        Ldyn_nms8_h_std  = alpha_L*Ldyn_nms8_h_std  + (1-alpha_L)*Ldyn_nms8_h_std_prev;
        Ldyn_nms8_h_pref = alpha_L*Ldyn_nms8_h_pref + (1-alpha_L)*Ldyn_nms8_h_pref_prev;
        Ldyn_oth_h_std   = alpha_L*Ldyn_oth_h_std   + (1-alpha_L)*Ldyn_oth_h_std_prev;
        Ldyn_oth_h_pref  = alpha_L*Ldyn_oth_h_pref  + (1-alpha_L)*Ldyn_oth_h_pref_prev;
    end
    if (alpha_L < 1.0)
        Ldyn_eu15_h_prev = Ldyn_eu15_h;  Ldyn_nms8_h_prev = Ldyn_nms8_h;  Ldyn_oth_h_prev = Ldyn_oth_h;
        Ldyn_eu15_l_prev = Ldyn_eu15_l;  Ldyn_nms8_l_prev = Ldyn_nms8_l;  Ldyn_oth_l_prev = Ldyn_oth_l;
        Ldyn_eu15_h_std_prev  = Ldyn_eu15_h_std;   Ldyn_eu15_h_pref_prev = Ldyn_eu15_h_pref;
        Ldyn_nms8_h_std_prev  = Ldyn_nms8_h_std;   Ldyn_nms8_h_pref_prev = Ldyn_nms8_h_pref;
        Ldyn_oth_h_std_prev   = Ldyn_oth_h_std;    Ldyn_oth_h_pref_prev  = Ldyn_oth_h_pref;
        have_prev_L = true;
    end

    Ldyn_h=Ldyn_eu15_h+Ldyn_nms8_h+Ldyn_oth_h;
    Ldyn_l=Ldyn_eu15_l+Ldyn_nms8_l+Ldyn_oth_l;
    Ldyn=Ldyn_h+Ldyn_l;
        
    %% Step 4 i(find temporary eq.)
     load('Base_year2002_skill.mat', 'VARn00', 'VALn00','Din00', 'xbilat00', 'Sn00', 'io', 'theta',...
         'B', 'N', 'maxit', 'tol', 'vfactor', 'VALS0', 'VALL0','sh_h00','sh_l00');  
      
    vfactor = -0.05;
    tol = 1E-10;  % original COPS value, paired with the sum-of-squares metric in
    %% SWEEP-ENGINE ADDITION (spectrum hygiene): TEMPEQ_TOL env override — finite-
    %% difference Jacobian-vector products need the inner solver at 1e-12 so that
    %% eigenvalue error bars are ~1e-3, not ~0.06 (A1 test battery). Default 1e-10.
    tt_env = getenv('TEMPEQ_TOL');
    if ~isempty(tt_env); tol = str2double(tt_env); end
                  % tempeq_skill.m. 1E-6 under that metric allows per-entry wage
                  % errors ~1e-3 (it was calibrated for the per-entry max-abs
                  % metric used 2026-05-05..06-08) and displaces the neutral-run
                  % fixed point: hats=1 stops being a temporary-equilibrium
                  % solution at the non-data-pinned slots (>=14), which the slot-2
                  % surprise equation amplifies by 1/(1-beta)~33 into CEV ~ -0.1%.
    max_iter = 1E+20;
    %% SWEEP-ENGINE ADDITION (hang guard): TEMPEQ_MAXIT env caps the inner
    %% temporary-equilibrium iterations. The original cap is effectively
    %% unbounded, so a probe evaluation pushed into a region where the
    %% tatonnement cannot reach tolerance hangs forever (observed 2026-06-12).
    %% A capped probe returns a slightly imprecise temporary equilibrium,
    %% which the outer solver tolerates. Default: unchanged (unbounded).
    tm_env = getenv('TEMPEQ_MAXIT');
    if ~isempty(tm_env); max_iter = round(str2double(tm_env)); end

    % Set initial conditions
    VALn0 = VALn00;
    VARn0 = VARn00;
    sh_h=sh_h00;
    sh_l=sh_l00;
    Sn = Sn00;
    pi = Din00;
    xbilat = xbilat00;
    
    realwages_h = ones(1, N, time);
    realwages_l = ones(1, N, time);
    realwages_h_std  = ones(1, N, time);  % Approach 2: after-tax wedge for std pool at PT sector 1
    realwages_h_pref = ones(1, N, time);  % Approach 2: after-tax wedge for pref pool at PT sector 1
    aux_nonemp=ones(N,time); %non employment income with public goods

    %% Approach 2: budget-clearing closure parameters (paper Sec 4 lambda-rescale).
    %% Closure follows paper main.tex Sec 4: rescale Benabou-HSV lambda at PT to
    %% clear the budget. Average rate at wage w under HSV is 1 - lambda*w^(-tau_prog).
    %% Three groups: std h-skill, pref h-skill (fixed at tau_P_pref), l-skill.
    %% At PT, lambda adjusts to clear PT's budget at the baseline G level.
    %% Other countries: rates frozen at baseline (plan sec 2.17).
    tau_h_base = 0.24;  % main.tex Table 5 native row (avg rate at calibrated w_h)
    tau_l_base = 0.18;  % main.tex Table 5 native row (avg rate at calibrated w_l)
    if ~isempty(tau_h_base_override); tau_h_base = tau_h_base_override; end
    if ~isempty(tau_l_base_override); tau_l_base = tau_l_base_override; end
    if isempty(tau_P_pref_override); tau_P_pref = 0.235; else; tau_P_pref = tau_P_pref_override; end  % 0.235=0.5pp cut; TAU_P_PREF env override (neutral audit: set =tau_h_base)
    tau_progressivity = 0.216;  % paper Table 1: Bénabou-HSV progressivity (Holter et al.)
    rho_c      = VALL0(i) / VALS0(i);          % ~3.93 at PT (slot-1 dollar wage bill ratio)
    L_h_base_c = L0_eu15_h(i) + L0_nms8_h(i) + L0_oth_h(i);
    L_l_base_c = L0_eu15_l(i) + L0_nms8_l(i) + L0_oth_l(i);
    %% lambda_hat_path: lambda_PT(t) / lambda_PT_baseline. Closure variable.
    %% f_std_path, f_l_path: multiplicative wedges on realwages at PT under
    %% lambda-rescale. tau_h_t_path, tau_l_t_path: implied average rates.
    lambda_hat_path = ones(1, time);
    tau_h_t_path = tau_h_base * ones(1, time);
    tau_l_t_path = tau_l_base * ones(1, time);
    f_std_path   = ones(1, time);
    f_l_path     = ones(1, time);
    f_pref_path  = ones(1, time);
    %% CONVENTION FIX: f_*_path store the LEVEL wedges F(t) (after-tax wedge on the
    %% counterfactual-ratio wage). What realwages receives is their SEQUENTIAL change
    %% F(t)/F(t-1), logged in mult_*_path; What_*_path log the counterfactual/baseline
    %% wage level ratios at PT.
    mult_std_path  = ones(1, time);
    mult_pref_path = ones(1, time);
    mult_l_path    = ones(1, time);
    What_h_path    = ones(1, time);
    What_l_path    = ones(1, time);
    %% Slot-1 pref wedge (reference). The applied wedge is the same-date
    %% f_pref_t = (1-tau_P_pref)/(1-tau_bar_h_t), computed in the closure block:
    %% the flat preferential rate is compared to the drifting baseline schedule
    %% rate (convention resolved 2026-06-10, PATH_SOLVER_SPEC.md sec 9).
    f_pref_const = (1 - tau_P_pref) / (1 - tau_h_base);
    G_level_hat  = tau_h_base + tau_l_base * rho_c;  % baseline G in W_h_bar units
    G_actual_path = G_level_hat * ones(1, time);
    T_hat_path    = zeros(1, time);  % retained as zero under lambda-rescale (diagnostic).
    alpha_l_path  = ones(1, time);
    %% Shadow lump-sum tax: T_hat that WOULD clear budget if LS closure were active.
    %% Diagnostic only under lambda-rescale; provides comparison to the deprecated LS path.
    tau_h_t_candidate_path = tau_h_base * ones(1, time);
    alpha_tau = 1.0;  % unused under lambda-rescale closure (kept for compatibility)
    %% Per-period diagnostic dumps (filled inside the t loop, used at the end).
    diag_w_h_hat   = zeros(1, time);
    diag_w_l_hat   = zeros(1, time);
    diag_L_h_std   = zeros(1, time);
    diag_L_h_pref  = zeros(1, time);
    diag_L_l_c     = zeros(1, time);
    diag_L_h_base  = zeros(1, time);
    diag_L_l_base  = zeros(1, time);
    diag_W_h_base  = zeros(1, time);
    diag_W_l_base  = zeros(1, time);
    diag_rho_t     = zeros(1, time);
    diag_G_base_t  = zeros(1, time);
    diag_L_h_std_neemp  = zeros(1, time);
    diag_L_h_pref_neemp = zeros(1, time);
    diag_L_l_c_neemp    = zeros(1, time);
    %% Std-pool outflow rate at PT from origin states i (employed) and i2 (non-employed).
    %% Per-period probability mass leaving the std pool at PT for any non-PT destination state.
    %% Plan §2.12: outflow_std_b(t) = sum_{d not in [i, i2]} mu_b_h_std(s, d, t) for s in {i, i2}.
    %% Compared against the same row sum of the baseline mu (no-policy outflow rate).
    diag_outflow_std_eu15_emp    = zeros(1, time);
    diag_outflow_std_nms8_emp    = zeros(1, time);
    diag_outflow_std_oth_emp     = zeros(1, time);
    diag_outflow_std_eu15_neemp  = zeros(1, time);
    diag_outflow_std_nms8_neemp  = zeros(1, time);
    diag_outflow_std_oth_neemp   = zeros(1, time);
    diag_outflow_base_emp        = zeros(1, time);
    diag_outflow_base_neemp      = zeros(1, time);

    %% Route 1: closure-internal quantities (not stored elsewhere) for spike diagnosis.
    diag_h_factor   = zeros(1, time);
    diag_l_factor   = zeros(1, time);
    diag_denom      = zeros(1, time);  % h_factor*A_hat + l_factor*C_hat
    diag_A_hat      = zeros(1, time);
    diag_B_hat      = zeros(1, time);
    diag_C_hat      = zeros(1, time);
    diag_L_std_hat  = zeros(1, time);
    diag_L_pref_hat = zeros(1, time);
    diag_L_l_hat    = zeros(1, time);
    %% Reformulated closure diagnostics: raw (unclamped) lambda, clamp-binding flag,
    %% and the budget residual at the applied lambda (0 = budget cleared).
    diag_lambda_raw    = ones(1, time);
    diag_clamp_binding = zeros(1, time);
    diag_budget_resid  = zeros(1, time);

    Ltemp_h = permute(Ldyn_h, [2, 1, 3]);
    Ltemp_l = permute(Ldyn_l, [2, 1, 3]);
    Ltemp = permute(Ldyn, [2, 1, 3]);

    %% CONVENTION FIX: counterfactual/baseline wage LEVEL ratios at PT, anchored at
    %% w'(1)=w(1) (the counterfactual coincides with the baseline before the
    %% announcement); accumulated inside the loop from the sequential tempeq ratios.
    %% F_*_prev hold the previous slot's level wedges for the sequential application.
    What_h_lev = 1;  What_l_lev = 1;
    F_std_prev = 1;  F_pref_prev = 1;  F_l_prev = 1;

        for t = 1:  time-2;

        % Shocks
        kappa_hat = ones(N, N); % relative change in trade costs
        kappa_hat=1./tau_hat(:,:,t); %shock to tariffs
        %taup=series_tau(:,:,t+2);
        taup=series_tau(:,:,t+1);  % matches Preferential_tax timing (no tariff changes)
        
        Snp = 0*Sn;
        
        % Initialize vectors of factor prices (w,r) and good prices (p)
        om_s = ones(N,1); % initial guess for factor prices
        om_l = ones(N,1);
          
        % Initialize vectors of factor/good prices (w, r, p)
        Ln_hat = ones(N,1);
       Ln_hat(1:17,:)= (Ltemp(1:N-1, :, t+1)+Ltemp(N:2*N-2, :, t+1))./(Ltemp(1:N-1, :, t)+Ltemp(N:2*N-2, :, t)); 
         
        Ln_hat_h = ones(N,1);
        Ln_hat_h(1:17,:)= Ltemp_h(1:N-1, :, t+1)./Ltemp_h(1:N-1, :, t); 
        
        Ln_hat_l = ones(N,1);
       Ln_hat_l(1:17,:)= Ltemp_l(1:N-1, :, t+1)./Ltemp_l(1:N-1, :, t); 
        
        pi_tilde0 = Din_baseline(:, :, t);
        pi_tilde1 = Din_baseline(:, :, t+1);
        w0_h = wages0_h(:, t+1);
        w0_l = wages0_l(:, t+1);
        
        Ln_hat00_h=ones(N,1);
        Ln_hat00_l=ones(N,1);
        Ln_hat00=ones(N,1);
        Ln_hat00_h(1:17,1) = Ln_hat0_h(:,  t+1); 
        Ln_hat00_l(1:17,1) = Ln_hat0_l(:,  t+1);
        Ln_hat00(1:17,1) =   Ln_hat0(:,  t+1);
        
        
        sh_hat_h00=series_sh_hat_h(:,t+1);
        sh_hat_l00=series_sh_hat_l(:,t+1);
        
        sh_l0=series_sh_l(:,t);
        sh_h0=series_sh_h(:,t);
        
        om0_s = (w0_h.^(1./(1-B))).*((Ln_hat00_h./sh_hat_h00).^(B./(1-B)));
        om0_l = (w0_l.^(1./(1-B))).*((Ln_hat00_l./sh_hat_l00).^(B./(1-B)));
        om_s=om0_s;
        om_l=om0_l;
      
        a=agglom_a; %agglomeration elasticity (default 1; AGGLOM_A env override)
        Ahat= (Ln_hat./Ln_hat00).^a;
       
        [om_s, om_l, wf0s, wf0l, VARnp, VALLp, VALSp, VALnp,Phat, rf0, p_hat, pip, Xp, Snp,...
    xbilatp,sh_hp,sh_lp] = tempeq_skill( VARn0,  pi, Snp,...
    kappa_hat, io, theta, B, J,... 
    N, max_iter, tol, vfactor, pi_tilde1, pi_tilde0, om_s, om_l, om0_s, om0_l, taup, sh_h, sh_l, sh_h0, sh_l0,rho, Ln_hat_h, Ln_hat_l, VALS0, VALL0, Ahat);
 
  VALn0=VALnp;
  VARn0=VARnp;
  Sn=Snp;
  Din=pip;
  pi=Din;
  xbilat=xbilatp;  
  
  VALn0=VALnp;
  VARn0=VARnp;
  VALS0 = VALSp;
  VALL0 = VALLp;
  sh_h=sh_hp;
  sh_l=sh_lp;
  Sn=Snp;
  Din=pip;
  pi=Din;
  xbilat=xbilatp;

           %% Approach 2: lump-sum tax on residents (all but pref pool).
           %% Proportional rates fixed at baseline (tau_h_t = tau_h_base, tau_l = tau_l_base).
           %% A per-capita lump-sum T(t) is collected from std h-skill employed and
           %% l-skill employed at PT to clear the budget G_cf(t) = G_base(t). T is in
           %% units of W_h_base(t). Pref pool exempt. Non-employed exempt (their flow
           %% utility is home production, untaxed).
           %%
           %% Budget identity (in W_h_base*L_h_base units):
           %%   tau_h_base*w_h_hat*L_std_hat + tau_P*w_h_hat*L_pref_hat
           %%   + tau_l*rho_t*w_l_hat*L_l_hat + T_hat*(L_std_hat + q_t*L_l_hat)
           %%   = G_level_hat_t = tau_h_base + tau_l*rho_t
           %% with q_t = L_l_base(t)/L_h_base(t).
           %%
           %% Effect on V via realwages: std worker disposable income hat becomes
           %%   w_h_real_hat - T_hat/(1 - tau_h_base);
           %% l-skill worker disposable income hat (at PT) becomes
           %%   w_l_real_hat - T_hat/((1 - tau_l_base)*alpha_dollar_t),
           %% where alpha_dollar_t = W_l_base_dollar(t)/W_h_base_dollar(t) is the
           %% same-t baseline DOLLAR wage ratio l/h. In hat decomposition,
           %% alpha_dollar_t = (rho_c/q_slot1) * (w_l_hat_baseline/w_h_hat_baseline)
           %% with rho_c = VALL0(i)/VALS0(i) the slot-1 wage-bill ratio and
           %% q_slot1 = L_l_base_c/L_h_base_c the slot-1 labor ratio. The earlier
           %% buggy version used wage hats only, missing the rho_c/q_slot1 factor.
           %% CONVENTION FIX: wf0s/w0_h is the SEQUENTIAL counterfactual/baseline
           %% wage ratio (both numerator and denominator are period-over-period
           %% changes). The budget bases need the same-date LEVEL ratio w'(t)/w(t):
           %% the cumulative product, anchored at w'(1)=w(1).
           What_h_lev   = What_h_lev * (wf0s(i) / w0_h(i));
           What_l_lev   = What_l_lev * (wf0l(i) / w0_l(i));
           w_h_hat_c    = What_h_lev;
           w_l_hat_c    = What_l_lev;
           L_std_h_c    = Ldyn_eu15_h_std(1, i, t+1) + Ldyn_nms8_h_std(1, i, t+1) + Ldyn_oth_h_std(1, i, t+1);
           L_pref_h_c   = Ldyn_eu15_h_pref(1, i, t+1) + Ldyn_nms8_h_pref(1, i, t+1) + Ldyn_oth_h_pref(1, i, t+1);
           L_l_c        = Ldyn_eu15_l(1, i, t+1) + Ldyn_nms8_l(1, i, t+1) + Ldyn_oth_l(1, i, t+1);
           L_h_base_t   = series_Ldyn_eu15_h(1, i, t+1) + series_Ldyn_nms8_h(1, i, t+1) + series_Ldyn_oth_h(1, i, t+1);
           L_l_base_t   = series_Ldyn_eu15_l(1, i, t+1) + series_Ldyn_nms8_l(1, i, t+1) + series_Ldyn_oth_l(1, i, t+1);
           L_std_h_hat  = L_std_h_c  / L_h_base_t;
           L_pref_h_hat = L_pref_h_c / L_h_base_t;
           L_l_hat      = L_l_c      / L_l_base_t;
           %% CONVENTION FIX: baseline wage LEVEL vs the base year (cumulative
           %% product of the sequential series), not the raw sequential entry.
           W_h_base_t   = W_base_h_lev(i, t+1);
           W_l_base_t   = W_base_l_lev(i, t+1);
           %% Schedule reference wages (Treatment 2, author decision 2026-07-02):
           %% under indexation the schedule's reference wage moves with the
           %% baseline wage, so the SCHEDULE terms (tau_bar, h_factor, l_factor)
           %% see W_*_sched = 1 — baseline rates constant at calibration, rates
           %% moved only by the policy deviation What_*. The wage-bill ratio
           %% rho_c_t and alpha_l_t keep the actual baseline drift (W_*_base_t):
           %% they are tax-base sizes, not schedule objects.
           if tau_indexed
               W_h_sched_t = 1;           W_l_sched_t = 1;
           else
               W_h_sched_t = W_h_base_t;  W_l_sched_t = W_l_base_t;
           end
           q_slot1      = L_l_base_c / L_h_base_c;
           q_t          = L_l_base_t / L_h_base_t;
           %% rho_c_t: same-t baseline DOLLAR wage-bill ratio low/high at PT.
           %% rho_c (slot 1) = VALL0(i)/VALS0(i) is the slot-1 dollar wage bill
           %% ratio. W_*_base_t (above) are baseline wage LEVELS vs slot 1, so the
           %% same-t correction adds the labor-ratio drift and the wage-level ratio.
           %% (The raw series_wages_* entries are sequential changes — 2026-07-02 fix.)
           rho_c_t      = rho_c * (q_t / q_slot1) * (W_l_base_t / W_h_base_t);
           %% alpha_l_t: same-t baseline DOLLAR wage ratio l/h. Slot-1 dollar
           %% wage ratio = rho_c / q_slot1 (~0.474 at PT). Same-t correction
           %% multiplies by the wage-hat ratio.
           alpha_l_t    = (rho_c / q_slot1) * (W_l_base_t / W_h_base_t);
           %% Baseline rates: the HSV schedule at the schedule-reference wage.
           %% Indexed (default, Treatment 2): W_sched = 1 — baseline rates constant
           %% at the calibrated tau_*_base; the counterfactual deviation is the
           %% policy alone. Unindexed (TAU_INDEXED=0): W_sched = the drifting
           %% baseline wage level — rates rise with baseline wage growth (the
           %% superseded 2026-06-10 design).
           tau_bar_h_t  = 1 - (1 - tau_h_base) * W_h_sched_t^(-tau_progressivity);
           tau_bar_l_t  = 1 - (1 - tau_l_base) * W_l_sched_t^(-tau_progressivity);
           G_level_hat_t = tau_bar_h_t + tau_bar_l_t * rho_c_t;
           %% Lambda-rescale closure (paper Sec 4). Under Benabou-HSV with single
           %% lambda_PT clearing PT's budget: avg_tax(w) = 1 - lambda*w^(-tau_prog).
           %% Std h-skill: tau_h_std = 1 - lambda_hat*(1-tau_h_base)*w_h_hat^(-tau_prog)
           %% l-skill:     tau_l    = 1 - lambda_hat*(1-tau_l_base)*w_l_hat^(-tau_prog)
           %% Pref h-skill is exempt: stays at tau_P_pref. Closed-form lambda_hat:
           %%   lambda_hat = (A + B + C - G_target) / (h_factor*A + l_factor*C)
           %% with A = w_h*L_std, B = tau_P*w_h*L_pref, C = rho_t*w_l*L_l;
           %% h_factor = (1-tau_h_base)*w_h_hat^(-tau_prog),
           %% l_factor = (1-tau_l_base)*w_l_hat^(-tau_prog).
           %% At slot 1 (w_hat=1, L_std=1, L_pref=0, L_l=1): lambda_hat = 1 by construction.
           %% Continuation interpolation (PREF_MIX = m): fraction (1-m) of the
           %% pref base on the schedule (inside A_hat), fraction m at the flat
           %% rate (revenue in B_hat). m=0: neutral anchor; m=1: policy.
           A_hat    = w_h_hat_c * (L_std_h_hat + (1 - pref_mix) * L_pref_h_hat);
           B_hat    = pref_mix * tau_P_pref * w_h_hat_c * L_pref_h_hat;
           C_hat        = rho_c_t * w_l_hat_c * L_l_hat;
           %% h_factor / l_factor: HSV mapping from lambda_hat to avg rate.
           %% Derivation: tau_h_t = 1 - lambda_hat * (1-tau_h_base) * [w_cf/w_ref]^(-tau_prog)
           %% where w_cf/w_ref = w_h_hat_c * W_h_sched_t: the counterfactual wage
           %% relative to the schedule's reference wage. Indexed default: the
           %% reference moves with the baseline wage, so w_cf/w_ref = w_h_hat_c
           %% (the policy deviation alone). All factors are cumulative level
           %% ratios (2026-07-02 convention fix).
           h_factor     = (1 - tau_h_base) * (w_h_hat_c * W_h_sched_t)^(-tau_progressivity);
           l_factor     = (1 - tau_l_base) * (w_l_hat_c * W_l_sched_t)^(-tau_progressivity);
           closure_numer_t = A_hat + B_hat + C_hat - G_level_hat_t;
           closure_denom_t = h_factor * A_hat + l_factor * C_hat;
           %% --- Reformulated closure: feasibility clamp + optional lambda damping ---
           %% tau = 1 - lambda*factor in [tau_floor, tau_max]  <=>
           %%   lambda in [(1-tau_max)/factor, (1-tau_floor)/factor], for each of h, l.
           %% Intersect the two skills' admissible lambda ranges (single-lambda closure).
           if exist('INJECT_LAMBDA', 'var')
               lambda_raw_t = NaN;
               lambda_hat_t = INJECT_LAMBDA(t+1);
               lambda_applied_prev(t+1) = lambda_hat_t;
               clamp_binding_t = 0;
           else
               lambda_raw_t = closure_numer_t / closure_denom_t;
               if clamp_tau
                   lam_lo = max((1 - tau_max_clamp)   / h_factor, (1 - tau_max_clamp)   / l_factor);
                   lam_hi = min((1 - tau_floor_clamp) / h_factor, (1 - tau_floor_clamp) / l_factor);
                   lambda_clamped_t = min(max(lambda_raw_t, lam_lo), lam_hi);
                   clamp_binding_t  = double(lambda_clamped_t ~= lambda_raw_t);
               else
                   lambda_clamped_t = lambda_raw_t;
                   clamp_binding_t  = 0;
               end
               %% Damp across outer iters (clamp first, then damp). alpha_lambda=1 => no damping.
               lambda_hat_t = alpha_lambda * lambda_clamped_t ...
                            + (1 - alpha_lambda) * lambda_applied_prev(t+1);
               lambda_applied_prev(t+1) = lambda_hat_t;
           end
           %% SWEEP-ENGINE ADDITION (augmented solver): when sweep_phi supplies a
           %% lambda path (INJECT_LAMBDA), use it for the wedges instead of the
           %% per-period closed form; the budget then enters the outer solver as
           %% explicit rows denom*lambda - numer = 0 (multiplied form).
           tau_h_std_t  = 1 - lambda_hat_t * h_factor;
           tau_l_t      = 1 - lambda_hat_t * l_factor;
           f_std_t      = (1 - tau_h_std_t) / (1 - tau_bar_h_t);  % = lambda_hat * w_h_hat^(-tau_prog)
           f_l_t        = (1 - tau_l_t)     / (1 - tau_bar_l_t);  % = lambda_hat * w_l_hat^(-tau_prog)
           if budget_off
               %% B-iii: hold standard h and l rates at the baseline schedule (no
               %% closure wedge). The pref cut stays on. Budget gap left unrebated.
               f_std_t = 1.0;  f_l_t = 1.0;
               tau_h_std_t = tau_bar_h_t;  tau_l_t = tau_bar_l_t;  lambda_hat_t = 1.0;
           end
           %% Pref wedge: interpolated between the schedule wedge (f_std_t) and
           %% the flat preferential rate vs the same-date baseline schedule rate.
           f_pref_t = (1 - pref_mix) * f_std_t ...
                    + pref_mix * (1 - tau_P_pref) / (1 - tau_bar_h_t);
           %% Budget residual at the APPLIED lambda (W_h_base units). 0 when lambda is
           %% interior and undamped; nonzero flags an unfunded slot (clamp bound or mid-damp).
           budget_resid_t = closure_numer_t - lambda_hat_t * closure_denom_t;
           %% --- end reformulated closure ---
           lambda_hat_path(t+1) = lambda_hat_t;
           tau_h_t_path(t+1) = tau_h_std_t;
           tau_l_t_path(t+1) = tau_l_t;
           f_std_path(t+1)   = f_std_t;
           f_l_path(t+1)     = f_l_t;
           f_pref_path(t+1)  = f_pref_t;
           What_h_path(t+1)  = What_h_lev;
           What_l_path(t+1)  = What_l_lev;
           G_actual_path(t+1) = G_level_hat_t;  % clears exactly under lambda-rescale
           %% Shadow LS T_hat (diagnostic): what LS would have been at this state.
           %% Not used in iteration under lambda-rescale; logged for comparison.
           T_denom_diag = L_std_h_hat + q_t * L_l_hat;
           G_prop_LS    = tau_h_base * w_h_hat_c * L_std_h_hat ...
                        + tau_P_pref * w_h_hat_c * L_pref_h_hat ...
                        + tau_l_base * rho_c_t * w_l_hat_c * L_l_hat;
           T_hat_path(t+1)   = (G_level_hat_t - G_prop_LS) / T_denom_diag;
           %% Shadow proportional rate (tau_h_t-only closure): diagnostic comparison.
           tau_h_t_candidate = (tau_h_base + tau_l_base * rho_c_t ...
                              - tau_P_pref * w_h_hat_c * L_pref_h_hat ...
                              - tau_l_base * rho_c_t * w_l_hat_c * L_l_hat) ...
                             / (w_h_hat_c * L_std_h_hat);
           tau_h_t_candidate_path(t+1) = tau_h_t_candidate;
           alpha_l_path(t+1) = alpha_l_t;
           %% Route 1: closure-internal amplifiers for spike diagnosis.
           diag_h_factor(t+1)   = h_factor;
           diag_l_factor(t+1)   = l_factor;
           diag_denom(t+1)      = closure_denom_t;
           diag_A_hat(t+1)      = A_hat;
           diag_B_hat(t+1)      = B_hat;
           diag_C_hat(t+1)      = C_hat;
           diag_L_std_hat(t+1)  = L_std_h_hat;
           diag_L_pref_hat(t+1) = L_pref_h_hat;
           diag_L_l_hat(t+1)    = L_l_hat;
           diag_lambda_raw(t+1)    = lambda_raw_t;
           diag_clamp_binding(t+1) = clamp_binding_t;
           diag_budget_resid(t+1)  = budget_resid_t;
           %% Per-period diagnostic state (PT, sector 1 emp / sector 2 non-emp).
           diag_w_h_hat(t+1)   = w_h_hat_c;
           diag_w_l_hat(t+1)   = w_l_hat_c;
           diag_L_h_std(t+1)   = L_std_h_c;
           diag_L_h_pref(t+1)  = L_pref_h_c;
           diag_L_l_c(t+1)     = L_l_c;
           diag_L_h_base(t+1)  = L_h_base_t;
           diag_L_l_base(t+1)  = L_l_base_t;
           diag_W_h_base(t+1)  = W_h_base_t;
           diag_W_l_base(t+1)  = W_l_base_t;
           diag_rho_t(t+1)     = rho_c_t;
           diag_G_base_t(t+1)  = G_level_hat_t;
           diag_L_h_std_neemp(t+1)  = Ldyn_eu15_h_std(1, i2, t+1)  + Ldyn_nms8_h_std(1, i2, t+1)  + Ldyn_oth_h_std(1, i2, t+1);
           diag_L_h_pref_neemp(t+1) = Ldyn_eu15_h_pref(1, i2, t+1) + Ldyn_nms8_h_pref(1, i2, t+1) + Ldyn_oth_h_pref(1, i2, t+1);
           diag_L_l_c_neemp(t+1)    = Ldyn_eu15_l(1, i2, t+1) + Ldyn_nms8_l(1, i2, t+1) + Ldyn_oth_l(1, i2, t+1);
           %% Std-pool outflow rate at PT: probability mass leaving the std pool at PT
           %% for any non-PT destination state. Origins: i (employed) and i2 (non-employed).
           %% mu is row-stochastic by construction (sum_d mu(s,d) = 1), so outflow from
           %% PT = 1 - prob(stay-emp) - prob(switch-to-neemp). Exact and avoids vector
           %% indexing that triggered an Octave 3D-tensor indexing bug.
           diag_outflow_std_eu15_emp(t+1)   = 1 - mu_eu15_h_std(i,  i, t+1) - mu_eu15_h_std(i,  i2, t+1);
           diag_outflow_std_nms8_emp(t+1)   = 1 - mu_nms8_h_std(i,  i, t+1) - mu_nms8_h_std(i,  i2, t+1);
           diag_outflow_std_oth_emp(t+1)    = 1 - mu_oth_h_std(i,   i, t+1) - mu_oth_h_std(i,   i2, t+1);
           diag_outflow_std_eu15_neemp(t+1) = 1 - mu_eu15_h_std(i2, i, t+1) - mu_eu15_h_std(i2, i2, t+1);
           diag_outflow_std_nms8_neemp(t+1) = 1 - mu_nms8_h_std(i2, i, t+1) - mu_nms8_h_std(i2, i2, t+1);
           diag_outflow_std_oth_neemp(t+1)  = 1 - mu_oth_h_std(i2,  i, t+1) - mu_oth_h_std(i2,  i2, t+1);
           %% Baseline reference: EU-15 baseline (97.8% of pre-policy h-skill at PT per audit T3).
           diag_outflow_base_emp(t+1)       = 1 - mu_baseline_eu15_h(i,  i, t+1) - mu_baseline_eu15_h(i,  i2, t+1);
           diag_outflow_base_neemp(t+1)     = 1 - mu_baseline_eu15_h(i2, i, t+1) - mu_baseline_eu15_h(i2, i2, t+1);

           %Adjust real wages with price index Phat)
         %% PG variant (kept for reference; deactivated by alpha=0 below)
         %alpha=[0.208491979; 0.257003806;0.157763309;0.228906625;0.215113269;...
         %       0.308133528;0.181491956;0.171517815;0.251488209;0.163053985;...
         %       0.230000277;0.203164869;0.207104508;0.209422095;0.183216559;...
         %       0.185904936;0.208868624;0.186909573]; %share of public goods
         alpha=0;  % COPS baseline: no public-goods congestion. Approach 2 fiscal channel only.
            
          realwages_h(:, :, t+1) = (((wf0s./w0_h)./Phat).^(1-alpha)).*((Ln_hat./Ln_hat00).^(-alpha));
          realwages_l(:, :, t+1) = (((wf0l./w0_l)./Phat).^(1-alpha)).*((Ln_hat./Ln_hat00).^(-alpha));
          aux_nonemp(:,t+1)=((Ln_hat./Ln_hat00).^(-alpha));

          %% Approach 2 (lambda-rescale closure): apply HSV-implied wedges at PT.
          %% f_std_t / f_pref_t / f_l_t are LEVEL wedges (after-tax wedge on the
          %% counterfactual-ratio wage). realwages carries the TIME-DIFFERENCE of
          %% the counterfactual/baseline consumption ratio, so the wedges enter as
          %% their SEQUENTIAL changes F(t+1)/F(t): the full wedge lands once at the
          %% announcement slot and the multiplier returns to ~1 as the paths settle
          %% — a one-off permanent tax change, not a compounding one.
          %% (CONVENTION FIX 2026-07-02, docs/cops_hat_convention_fix_plan.md;
          %% the previous code multiplied the full level wedge at every slot.)
          %% l-skill outside PT and h-skill outside PT: unchanged (paper Sec 4 +
          %% plan sec 2.17 hold rates outside PT at baseline).
          %% Non-employed (sector 2 = i2): no labour income, no wedge.
          mult_std_t  = f_std_t  / F_std_prev;
          mult_pref_t = f_pref_t / F_pref_prev;
          mult_l_t    = f_l_t    / F_l_prev;
          F_std_prev  = f_std_t;  F_pref_prev = f_pref_t;  F_l_prev = f_l_t;
          mult_std_path(t+1)  = mult_std_t;
          mult_pref_path(t+1) = mult_pref_t;
          mult_l_path(t+1)    = mult_l_t;
          realwages_h_std(:, :, t+1)  = realwages_h(:, :, t+1);
          realwages_h_pref(:, :, t+1) = realwages_h(:, :, t+1);
          realwages_h_std(:, i, t+1)  = realwages_h(:, i, t+1) * mult_std_t;
          realwages_h_pref(:, i, t+1) = realwages_h(:, i, t+1) * mult_pref_t;
          realwages_l(:, i, t+1)      = realwages_l(:, i, t+1) * mult_l_t;

    end

   N=17;
    realwagesaux_h= realwages_h;
    realwagesaux_l= realwages_l;
    %% SWEEP-ENGINE DIAGNOSTIC (welfare-sign experiment, 2026-06-14): WAGE_INJECT
    %% multiplies the UNSPLIT high-skill real wage at PT-employed (state i=16) by
    %% a constant factor at all periods — a pure exogenous real-wage change on the
    %% ORIGINAL COPS channel (realwagesaux_h -> rwage_h -> unsplit V_eu15_h),
    %% untouched by the std/pref split. Run in neutral config to read the sign of
    %% V_eu15_h(16,2): factor>1 (wage gain) giving V_hat<1 means a wage gain maps
    %% to a gain (so the f_pref>1 wedge that yields V_hat>1 is backwards); V_hat>1
    %% means the inherited channel maps a wage gain to V_hat>1. Default off.
    wi_env = getenv('WAGE_INJECT');
    if ~isempty(wi_env)
        realwagesaux_h(:, i, :) = realwagesaux_h(:, i, :) * str2double(wi_env);
    end
    
    % Check value functions
    sigma = 1./0.62; 
    
    rwage_h = ones(N+1, 1, time);    
    for t = 1:time
        rwage_h(:, :, t) = permute(realwagesaux_h(:, :, t), [2, 1, 3]);
    end
    
   rwage_l = ones(N+1, 1, time);    
    for t = 1:time
        rwage_l(:, :, t) = permute(realwagesaux_l(:, :, t), [2, 1, 3]);
    end
    
    %correction non emp public goods
    aux_nonemp2=ones(N,1,time);%non employment with publi
    for t = 1:time
        aux_nonemp2(:,:,t) =aux_nonemp(1:N,t) ;
    end
    rwage_h=[rwage_h(1:17,:,:); (aux_nonemp2(1:N,:,:))]; %non employment with public goods
    rwage_l=[rwage_l(1:17,:,:); (aux_nonemp2(1:N,:,:))]; %non employment with public goods

    %% Approach 2: build rwage_h_std and rwage_h_pref the same way, sourced from
    %% the std/pref-adjusted realwages arrays. The std array carries the LS
    %% subtraction; the pref array carries the multiplicative pref tax wedge.
    rwage_h_std  = ones(N+1, 1, time);
    rwage_h_pref = ones(N+1, 1, time);
    for t = 1:time
        rwage_h_std(:,  :, t) = permute(realwages_h_std(:,  :, t), [2, 1, 3]);
        rwage_h_pref(:, :, t) = permute(realwages_h_pref(:, :, t), [2, 1, 3]);
    end
    rwage_h_std  = [rwage_h_std(1:17,:,:);  (aux_nonemp2(1:N,:,:))];
    rwage_h_pref = [rwage_h_pref(1:17,:,:); (aux_nonemp2(1:N,:,:))];

    % Iterate first period 
    % Deal with t = 1 (surprise in change of fundamentals)
    t1 = 2;
    maxdixH1 = 10;
    Y0 = V_eu15_h(:, t1);
    Hvf = Y0'; % V1(theta)-V0/(V1-V0)
    Hvf2 = repmat(Hvf,J*N, 1); 
    mu1_tilde_eu15_h = mu_baseline_eu15_h(:, :, 2).*(Yaux1_eu15_h(:, :, t1).^beta);  %% H3-1 fix 2026-06-10: c-columns carry V_eff, consistent with forward mu (line ~376) and the std/pref slot-2 blocks
    auxshock = rwage_h(:, :, t1).^(1/sigma);
    auxvf4 = auxshock';
    auxw = reshape(auxvf4, J*N, 1);
    auxvf1 = mu1_tilde_eu15_h.*(Yaux1_eu15_h(:, :, t1+1).^beta);
    Y_aux = sum(auxvf1')';
    Y1 = Y_aux.*auxw;
    Y0_eu15_h = Y1;

    %% Approach 2: Y0_eu15_h_std and Y0_eu15_h_pref using std/pref V, Yaux1, and rwage.
    %% Build mu1_tilde from Yaux1_*_h_{std|pref}(:,:,t1) (hybrid: V_h offshore,
    %% V_std/V_pref at c-columns) to match the forward mu recursion at line ~289.
    %% Earlier code broadcast V_std/V_pref across ALL destination columns, which
    %% mismatched offshore once V_h diverged from V_std/V_pref and acted as a
    %% persistent slot-2 boundary forcing term.
    %% rwage_h_std / rwage_h_pref already encode the std LS subtraction and pref
    %% multiplicative wedge respectively.
    % std:
    mu1_tilde_eu15_h_std_init  = mu_baseline_eu15_h(:, :, 2).*(Yaux1_eu15_h_std(:, :, t1).^beta);
    auxshock_std  = rwage_h_std(:, :, t1).^(1/sigma);
    auxw_std = reshape(auxshock_std', J*N, 1);
    auxvf1_std = mu1_tilde_eu15_h_std_init.*(Yaux1_eu15_h_std(:, :, t1+1).^beta);
    Y0_eu15_h_std = sum(auxvf1_std')'.*auxw_std;
    % pref:
    mu1_tilde_eu15_h_pref_init = mu_baseline_eu15_h(:, :, 2).*(Yaux1_eu15_h_pref(:, :, t1).^beta);
    auxshock_pref = rwage_h_pref(:, :, t1).^(1/sigma);
    auxw_pref = reshape(auxshock_pref', J*N, 1);
    auxvf1_pref = mu1_tilde_eu15_h_pref_init.*(Yaux1_eu15_h_pref(:, :, t1+1).^beta);
    Y0_eu15_h_pref = sum(auxvf1_pref')'.*auxw_pref;

    t1 = 2;
    maxdixH1 = 10;
    Y0 = V_eu15_l(:, t1);
    Hvf = Y0'; % V1(theta)-V0/(V1-V0)
    Hvf2 = repmat(Hvf,J*N, 1); 
    mu1_tilde_eu15_l = mu_baseline_eu15_l(:, :, 2).*(Hvf2.^beta);
    auxshock = rwage_l(:, :, t1).^(1/sigma);
    auxvf4 = auxshock';
    auxw = reshape(auxvf4, J*N, 1); 
    auxvf1 = mu1_tilde_eu15_l.*(Yaux1_eu15_l(:, :, t1+1).^beta);
    Y_aux = sum(auxvf1')';
    Y1 = Y_aux.*auxw;
    Y0_eu15_l = Y1;
    
    t1 = 2;
    maxdixH1 = 10;
    Y0 = V_nms8_h(:, t1);
    Hvf = Y0'; % V1(theta)-V0/(V1-V0)
      Hvf2 = repmat(Hvf,J*N, 1); 
    mu1_tilde_nms8_h = mu_baseline_nms8_h(:, :, 2).*(Yaux1_nms8_h(:, :, t1).^beta);  %% H3-1 fix 2026-06-10: c-columns carry V_eff, consistent with forward mu (line ~376) and the std/pref slot-2 blocks
    auxshock = rwage_h(:, :, t1).^(1/sigma);
    auxvf4 = auxshock';
    auxw = reshape(auxvf4, J*N, 1);
    auxvf1 = mu1_tilde_nms8_h.*(Yaux1_nms8_h(:, :, t1+1).^beta);
    Y_aux = sum(auxvf1')';
    Y1 = Y_aux.*auxw;
    Y0_nms8_h = Y1;

    %% Approach 2: Y0_nms8_h_std and Y0_nms8_h_pref using Yaux1 hybrid.
    mu1_tilde_nms8_h_std_init  = mu_baseline_nms8_h(:, :, 2).*(Yaux1_nms8_h_std(:, :, t1).^beta);
    auxvf1_std = mu1_tilde_nms8_h_std_init.*(Yaux1_nms8_h_std(:, :, t1+1).^beta);
    Y0_nms8_h_std = sum(auxvf1_std')'.*auxw_std;
    mu1_tilde_nms8_h_pref_init = mu_baseline_nms8_h(:, :, 2).*(Yaux1_nms8_h_pref(:, :, t1).^beta);
    auxvf1_pref = mu1_tilde_nms8_h_pref_init.*(Yaux1_nms8_h_pref(:, :, t1+1).^beta);
    Y0_nms8_h_pref = sum(auxvf1_pref')'.*auxw_pref;

    t1 = 2;
    maxdixH1 = 10;
    Y0 = V_nms8_l(:, t1);
    Hvf = Y0'; % V1(theta)-V0/(V1-V0)
    Hvf2 = repmat(Hvf,J*N, 1); 
    mu1_tilde_nms8_l = mu_baseline_nms8_l(:, :, 2).*(Hvf2.^beta);
    auxshock = rwage_l(:, :, t1).^(1/sigma);
    auxvf4 = auxshock';
    auxw = reshape(auxvf4, J*N, 1); 
    auxvf1 = mu1_tilde_nms8_l.*(Yaux1_nms8_l(:, :, t1+1).^beta);
    Y_aux = sum(auxvf1')';
    Y1 = Y_aux.*auxw;
    Y0_nms8_l = Y1;
     
     
   t1 = 2;
    maxdixH1 = 10;
    Y0 = V_oth_h(:, t1);
    Hvf = Y0'; % V1(theta)-V0/(V1-V0)
     Hvf2 = repmat(Hvf,J*N, 1); 
    mu1_tilde_oth_h = mu_baseline_oth_h(:, :, 2).*(Yaux1_oth_h(:, :, t1).^beta);  %% H3-1 fix 2026-06-10: c-columns carry V_eff, consistent with forward mu (line ~376) and the std/pref slot-2 blocks
    auxshock = rwage_h(:, :, t1).^(1/sigma);
    auxvf4 = auxshock';
    auxw = reshape(auxvf4, J*N, 1);
    auxvf1 = mu1_tilde_oth_h.*(Yaux1_oth_h(:, :, t1+1).^beta);
    Y_aux = sum(auxvf1')';
    Y1 = Y_aux.*auxw;
    Y0_oth_h = Y1;

    %% Approach 2: Y0_oth_h_std and Y0_oth_h_pref using Yaux1_*_std/pref hybrid
    %% (V_h offshore, V_std/V_pref at c) to match the forward mu recursion.
    mu1_tilde_oth_h_std_init  = mu_baseline_oth_h(:, :, 2).*(Yaux1_oth_h_std(:, :, t1).^beta);
    auxvf1_std = mu1_tilde_oth_h_std_init.*(Yaux1_oth_h_std(:, :, t1+1).^beta);
    Y0_oth_h_std = sum(auxvf1_std')'.*auxw_std;
    mu1_tilde_oth_h_pref_init = mu_baseline_oth_h(:, :, 2).*(Yaux1_oth_h_pref(:, :, t1).^beta);
    auxvf1_pref = mu1_tilde_oth_h_pref_init.*(Yaux1_oth_h_pref(:, :, t1+1).^beta);
    Y0_oth_h_pref = sum(auxvf1_pref')'.*auxw_pref;

     t1 = 2;
    maxdixH1 = 10;
    Y0 = V_oth_l(:, t1);
    Hvf = Y0'; % V1(theta)-V0/(V1-V0)
    Hvf2 = repmat(Hvf,J*N, 1); 
    mu1_tilde_oth_l = mu_baseline_oth_l(:, :, 2).*(Hvf2.^beta);
    auxshock = rwage_l(:, :, t1).^(1/sigma);
    auxvf4 = auxshock';
    auxw = reshape(auxvf4, J*N, 1); 
    auxvf1 = mu1_tilde_oth_l.*(Yaux1_oth_l(:, :, t1+1).^beta);
    Y_aux = sum(auxvf1')';
    Y1 = Y_aux.*auxw;
    Y0_oth_l = Y1;

    %%%%%Update wages, et al.
    rw = rwage_h.^(1/sigma);
    rw = permute(rw, [2, 1, 3]);

    
    rw_aux = zeros(J*N, 1, time);
    rw_aux1 = zeros(J*N, J*N, time);
    for t = 1:time
        rw_aux(:, :, t) = reshape(rw(:, :, t), J*N, 1);
        rw_aux1_h(:, :, t) = repmat(rw_aux(:, :, t), 1, J*N);
    end

    %% Approach 2: rw_aux1_h_std and rw_aux1_h_pref are built from rwage_h_std
    %% and rwage_h_pref (constructed earlier with the non-emp concat), which
    %% already encode both the multiplicative pref tax wedge for pref and the
    %% ADDITIVE lump-sum subtraction for std at PT sector 1.
    rw_std  = permute(rwage_h_std.^(1/sigma),  [2, 1, 3]);
    rw_pref = permute(rwage_h_pref.^(1/sigma), [2, 1, 3]);
    rw_aux1_h_std  = zeros(J*N, J*N, time);
    rw_aux1_h_pref = zeros(J*N, J*N, time);
    for t = 1:time
        rw_aux_std  = reshape(rw_std(:,  :, t), J*N, 1);
        rw_aux_pref = reshape(rw_pref(:, :, t), J*N, 1);
        rw_aux1_h_std(:,  :, t) = repmat(rw_aux_std,  1, J*N);
        rw_aux1_h_pref(:, :, t) = repmat(rw_aux_pref, 1, J*N);
    end

    rw = rwage_l.^(1/sigma);
    rw = permute(rw, [2, 1, 3]);

    rw_aux = zeros(J*N, 1, time);
    rw_aux1 = zeros(J*N, J*N, time);
    for t = 1:time
        rw_aux(:, :, t) = reshape(rw(:, :, t), J*N, 1);
        rw_aux1_l(:, :, t) = repmat(rw_aux(:, :, t), 1, J*N);
    end

    %% Step 5  (deliver new path of Y_{t+1}) 
    lambda = zeros(J*N, J*N, time); % NOT productivity
        % is [what/Phat]^{1/nu} * mu
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_eu15_h(:,:,t).*(mu_baseline_eu15_h(:, :, t+1)./mu_baseline_eu15_h(:, :, t)).*...
            mu_eu15_h(:, : ,t).*rw_aux1_h(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_eu15_h(:, :, t+1).^beta);
    end

    Y_eu15_h = zeros(J*N, 1, time);
    Y_eu15_h(:, :, 2) = Y0_eu15_h;
    for t = (t1 + 1):time
        Y_eu15_h(:, :, t) = sum(num(:,:,t)')';
    end
    Y_eu15_h(:, :, time) = 1;

    %% Approach 2: Y_eu15_h_std backward update (uses mu_eu15_h_std, rw_aux1_h_std, Yaux1_eu15_h_std).
    lambda = zeros(J*N, J*N, time);
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_eu15_h(:,:,t).*(mu_baseline_eu15_h(:, :, t+1)./mu_baseline_eu15_h(:, :, t)).*...
            mu_eu15_h_std(:, :, t).*rw_aux1_h_std(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_eu15_h_std(:, :, t+1).^beta);
    end
    Y_eu15_h_std = zeros(J*N, 1, time);
    Y_eu15_h_std(:, :, 2) = Y0_eu15_h_std;
    for t = (t1 + 1):time
        Y_eu15_h_std(:, :, t) = sum(num(:,:,t)')';
    end
    Y_eu15_h_std(:, :, time) = 1;

    %% Approach 2: Y_eu15_h_pref backward update.
    lambda = zeros(J*N, J*N, time);
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_eu15_h(:,:,t).*(mu_baseline_eu15_h(:, :, t+1)./mu_baseline_eu15_h(:, :, t)).*...
            mu_eu15_h_pref(:, :, t).*rw_aux1_h_pref(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_eu15_h_pref(:, :, t+1).^beta);
    end
    Y_eu15_h_pref = zeros(J*N, 1, time);
    Y_eu15_h_pref(:, :, 2) = Y0_eu15_h_pref;
    for t = (t1 + 1):time
        Y_eu15_h_pref(:, :, t) = sum(num(:,:,t)')';
    end
    Y_eu15_h_pref(:, :, time) = 1;

    lambda = zeros(J*N, J*N, time); % NOT productivity
        % is [what/Phat]^{1/nu} * mu
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_eu15_l(:,:,t).*(mu_baseline_eu15_l(:, :, t+1)./mu_baseline_eu15_l(:, :, t)).*...
            mu_eu15_l(:, : ,t).*rw_aux1_l(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_eu15_l(:, :, t+1).^beta);
    end
    
    Y_eu15_l = zeros(J*N, 1, time); 
    Y_eu15_l(:, :, 2) = Y0_eu15_l;
    for t = (t1 + 1):time
        Y_eu15_l(:, :, t) = sum(num(:,:,t)')';
    end    
    Y_eu15_l(:, :, time) = 1;
    

    lambda = zeros(J*N, J*N, time); % NOT productivity
        % is [what/Phat]^{1/nu} * mu
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_nms8_h(:,:,t).*(mu_baseline_nms8_h(:, :, t+1)./mu_baseline_nms8_h(:, :, t)).*...
            mu_nms8_h(:, : ,t).*rw_aux1_h(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_nms8_h(:, :, t+1).^beta);
    end
    
    Y_nms8_h = zeros(J*N, 1, time);
    Y_nms8_h(:, :, 2) = Y0_nms8_h;
    for t = (t1 + 1):time
        Y_nms8_h(:, :, t) = sum(num(:,:,t)')';
    end
    Y_nms8_h(:, :, time) = 1;

    %% Approach 2: Y_nms8_h_std and Y_nms8_h_pref backward updates.
    lambda = zeros(J*N, J*N, time);
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_nms8_h(:,:,t).*(mu_baseline_nms8_h(:, :, t+1)./mu_baseline_nms8_h(:, :, t)).*...
            mu_nms8_h_std(:, :, t).*rw_aux1_h_std(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_nms8_h_std(:, :, t+1).^beta);
    end
    Y_nms8_h_std = zeros(J*N, 1, time);
    Y_nms8_h_std(:, :, 2) = Y0_nms8_h_std;
    for t = (t1 + 1):time
        Y_nms8_h_std(:, :, t) = sum(num(:,:,t)')';
    end
    Y_nms8_h_std(:, :, time) = 1;

    lambda = zeros(J*N, J*N, time);
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_nms8_h(:,:,t).*(mu_baseline_nms8_h(:, :, t+1)./mu_baseline_nms8_h(:, :, t)).*...
            mu_nms8_h_pref(:, :, t).*rw_aux1_h_pref(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_nms8_h_pref(:, :, t+1).^beta);
    end
    Y_nms8_h_pref = zeros(J*N, 1, time);
    Y_nms8_h_pref(:, :, 2) = Y0_nms8_h_pref;
    for t = (t1 + 1):time
        Y_nms8_h_pref(:, :, t) = sum(num(:,:,t)')';
    end
    Y_nms8_h_pref(:, :, time) = 1;

    lambda = zeros(J*N, J*N, time); % NOT productivity
    % is [what/Phat]^{1/nu} * mu
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_nms8_l(:,:,t).*(mu_baseline_nms8_l(:, :, t+1)./mu_baseline_nms8_l(:, :, t)).*...
            mu_nms8_l(:, : ,t).*rw_aux1_l(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_nms8_l(:, :, t+1).^beta);
    end

    Y_nms8_l = zeros(J*N, 1, time); 
    Y_nms8_l(:, :, 2) = Y0_nms8_l;
    for t = (t1 + 1):time
        Y_nms8_l(:, :, t) = sum(num(:,:,t)')';
    end    
    Y_nms8_l(:, :, time) = 1;

    lambda = zeros(J*N, J*N, time); % NOT productivity
        % is [what/Phat]^{1/nu} * mu
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_oth_h(:,:,t).*(mu_baseline_oth_h(:, :, t+1)./mu_baseline_oth_h(:, :, t)).*...
            mu_oth_h(:, : ,t).*rw_aux1_h(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_oth_h(:, :, t+1).^beta);
    end
    
    Y_oth_h = zeros(J*N, 1, time);
    Y_oth_h(:, :, 2) = Y0_oth_h;
    for t = (t1 + 1):time
        Y_oth_h(:, :, t) = sum(num(:,:,t)')';
    end
    Y_oth_h(:, :, time) = 1;

    %% Approach 2: Y_oth_h_std and Y_oth_h_pref backward updates.
    lambda = zeros(J*N, J*N, time);
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_oth_h(:,:,t).*(mu_baseline_oth_h(:, :, t+1)./mu_baseline_oth_h(:, :, t)).*...
            mu_oth_h_std(:, :, t).*rw_aux1_h_std(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_oth_h_std(:, :, t+1).^beta);
    end
    Y_oth_h_std = zeros(J*N, 1, time);
    Y_oth_h_std(:, :, 2) = Y0_oth_h_std;
    for t = (t1 + 1):time
        Y_oth_h_std(:, :, t) = sum(num(:,:,t)')';
    end
    Y_oth_h_std(:, :, time) = 1;

    lambda = zeros(J*N, J*N, time);
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_oth_h(:,:,t).*(mu_baseline_oth_h(:, :, t+1)./mu_baseline_oth_h(:, :, t)).*...
            mu_oth_h_pref(:, :, t).*rw_aux1_h_pref(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_oth_h_pref(:, :, t+1).^beta);
    end
    Y_oth_h_pref = zeros(J*N, 1, time);
    Y_oth_h_pref(:, :, 2) = Y0_oth_h_pref;
    for t = (t1 + 1):time
        Y_oth_h_pref(:, :, t) = sum(num(:,:,t)')';
    end
    Y_oth_h_pref(:, :, time) = 1;

    lambda = zeros(J*N, J*N, time); % NOT productivity
        % is [what/Phat]^{1/nu} * mu
    for t = 1:time-1
        lambda(:, :, t+1) = tau_mig_hat_oth_l(:,:,t).*(mu_baseline_oth_l(:, :, t+1)./mu_baseline_oth_l(:, :, t)).*...
            mu_oth_l(:, : ,t).*rw_aux1_l(:, :, t+1);
        lambda(isnan(lambda)) = 0;
    end
    
    num = zeros(size(lambda));
    for t = 1:time-1
        num(:, :, t) = lambda(:, :, t).*(Yaux1_oth_l(:, :, t+1).^beta);
    end
    
    Y_oth_l = zeros(J*N, 1, time); 
    Y_oth_l(:, :, 2) = Y0_oth_l;
    for t = (t1 + 1):time
        Y_oth_l(:, :, t) = sum(num(:,:,t)')';
    end    
    Y_oth_l(:, :, time) = 1;

    Y=[Y_eu15_h;Y_nms8_h; Y_oth_h; Y_eu15_l;Y_nms8_l; Y_oth_l];
    
    %% Step 6 (store new set of initial conditions)
    % Store new values
     Ynew = zeros(6*N*J, time); 
    for t = 1:time
        Ynew(:, t) = [Y_eu15_h(:, :, t);Y_nms8_h(:, :, t); Y_oth_h(:,:,t);Y_eu15_l(:, :, t);Y_nms8_l(:, :, t); Y_oth_l(:,:,t)];
    end
    Ynew(:, time) = 1;
    V_eu15_h(:, time) = 1;
    V_nms8_h(:, time) = 1;
    V_oth_h(:, time) = 1;
    V_eu15_l(:, time) = 1;
    V_nms8_l(:, time) = 1;
    V_oth_l(:, time) = 1;
    V=[V_eu15_h;V_nms8_h; V_oth_h;V_eu15_l;V_nms8_l; V_oth_l];
    
    %% Step 7 (check and update guess)
    % Excess function — extended for Approach 2 to include V_h_std and V_h_pref
    % so the iteration runs until the std/pref tracks converge, not just the
    % unsplit V_h (which is unaffected when shock_val = 1).
    %% Load-bearing metric (2026-06-10): the split-track rows other than the PT
    %% c-rows [i i2] are outputs only — they re-enter nothing in the map (Yaux1
    %% copies the unsplit track everywhere and overrides only the c-columns).
    %% After the H3-1 fix the unsplit h-bloc c-rows are likewise inert (read only
    %% by legacy welfare diagnostics). The stopping metric Ydiff therefore reads
    %% (a) the unsplit stack excluding the h-bloc c-rows, and (b) the split tracks
    %% at the c-rows only. The full-stack metric is kept as Ydiff_full (logged,
    %% not used for stopping).
    Ydiff = zeros(time, 1);
    Ydiff_full = zeros(time, 1);
    d_unsplit_path = zeros(time, 1);  % Route 1: per-slot V-change by track
    d_std_path     = zeros(time, 1);
    d_pref_path    = zeros(time, 1);
    Y_h_std_2d  = reshape([Y_eu15_h_std;  Y_nms8_h_std;  Y_oth_h_std],  3*J*N, time);
    Y_h_pref_2d = reshape([Y_eu15_h_pref; Y_nms8_h_pref; Y_oth_h_pref], 3*J*N, time);
    V_h_std_2d  = [V_eu15_h_std;  V_nms8_h_std;  V_oth_h_std];
    V_h_pref_2d = [V_eu15_h_pref; V_nms8_h_pref; V_oth_h_pref];
    lb_split_rows   = [i i2, J*N+i J*N+i2, 2*J*N+i 2*J*N+i2];  % c-rows per h bloc in the 3-bloc stacks
    ghost_unsplit   = [i i2, J*N+i J*N+i2, 2*J*N+i 2*J*N+i2];  % h-bloc c-rows in the 6-bloc stack (first 3 blocs)
    lb_unsplit_rows = setdiff(1:6*J*N, ghost_unsplit);
    for t = 2:time
        d_unsplit = max(abs(Y(lb_unsplit_rows, t) - V(lb_unsplit_rows, t)));
        d_std     = max(abs(Y_h_std_2d(lb_split_rows, t)  - V_h_std_2d(lb_split_rows, t)));
        d_pref    = max(abs(Y_h_pref_2d(lb_split_rows, t) - V_h_pref_2d(lb_split_rows, t)));
        d_unsplit_path(t) = d_unsplit;
        d_std_path(t)     = d_std;
        d_pref_path(t)    = d_pref;
        Ydiff(t, 1) = max([d_unsplit, d_std, d_pref]);
        Ydiff_full(t, 1) = max([max(abs(Y(:, t) - V(:, t))), ...
                                max(abs(Y_h_std_2d(:, t)  - V_h_std_2d(:, t))), ...
                                max(abs(Y_h_pref_2d(:, t) - V_h_pref_2d(:, t)))]);
    end
    % Update guess
    V = 0.5*(Ynew + V);
    V_eu15_h=V(1:J*N,:);
    V_nms8_h=V(J*N+1:4*N,:);
    V_oth_h=V(4*N+1:6*N,:);
    V_eu15_l=V(6*N+1:8*N,:);
    V_nms8_l=V(8*N+1:10*N,:);
    V_oth_l=V(10*N+1:12*N,:);

    %% Approach 2: V damping for std/pref tracks at high-skill blocs (plan sec 2.10).
    V_eu15_h_std(:, time)  = 1;  V_eu15_h_pref(:, time) = 1;
    V_nms8_h_std(:, time)  = 1;  V_nms8_h_pref(:, time) = 1;
    V_oth_h_std(:, time)   = 1;  V_oth_h_pref(:, time)  = 1;
    Ynew_h_std  = zeros(3*J*N, time);
    Ynew_h_pref = zeros(3*J*N, time);
    for t = 1:time
        Ynew_h_std(:, t)  = [Y_eu15_h_std(:, :, t);  Y_nms8_h_std(:, :, t);  Y_oth_h_std(:, :, t)];
        Ynew_h_pref(:, t) = [Y_eu15_h_pref(:, :, t); Y_nms8_h_pref(:, :, t); Y_oth_h_pref(:, :, t)];
    end
    Ynew_h_std(:, time)  = 1;
    Ynew_h_pref(:, time) = 1;
    V_h_std  = [V_eu15_h_std;  V_nms8_h_std;  V_oth_h_std];
    V_h_pref = [V_eu15_h_pref; V_nms8_h_pref; V_oth_h_pref];
    %% D-1 asymmetric damping (2026-05-16): V_std and V_pref tracks where the
    %% fiscal feedback lives are damped at alpha=0.1 (90% old + 10% new), while
    %% the unsplit V (V_h, V_l) keeps the standard 0.5/0.5 mix at line 1483.
    %% In linear regime the LS Picard map has eigenvalue ~1.02 on V_std/V_pref;
    %% damping alpha gives effective eigenvalue 0.1*1.02 + 0.9 = 0.92 -> stable.
    %% Route 1 (2026-05-28): value taken from alpha_damp_split_cfg (env DAMP_SPLIT,
    %% default 0.1 = as-committed) so the 0.1-vs-0.5 configs can be compared.
    alpha_damp_split = alpha_damp_split_cfg;
    V_h_std  = alpha_damp_split*Ynew_h_std  + (1-alpha_damp_split)*V_h_std;
    V_h_pref = alpha_damp_split*Ynew_h_pref + (1-alpha_damp_split)*V_h_pref;
    V_eu15_h_std = V_h_std(1:J*N, :);
    V_nms8_h_std = V_h_std(J*N+1:4*N, :);
    V_oth_h_std  = V_h_std(4*N+1:6*N, :);
    V_eu15_h_pref = V_h_pref(1:J*N, :);
    V_nms8_h_pref = V_h_pref(J*N+1:4*N, :);
    V_oth_h_pref  = V_h_pref(4*N+1:6*N, :);

    Ymax = max(Ydiff); %#ok<NOPTS>
    fprintf('iter=%d  Ymax=%.6e\n', iter, Ymax); fflush(stdout);

    %% --- Route 1 diagnostic dump (does not affect iteration) ---
    %% Budget block fills slot t+1 for t=1..time-2, i.e. slots 2..time-1.
    slots = 2:(time-1);
    [mn_denom, k1] = min(diag_denom(slots));
    [mx_lam,   k2] = max(abs(lambda_hat_path(slots)));
    [mx_fstd,  k3] = max(f_std_path(slots));
    [mn_tauh,  k4] = min(tau_h_t_path(slots));
    [mx_wh,    k5] = max(abs(diag_w_h_hat(slots)));
    [mx_resid, k6] = max(abs(diag_budget_resid(slots)));
    n_clamp   = sum(diag_clamp_binding(slots));   % # slots where the tau-clamp bound this iter
    [~, Yarg] = max(Ydiff);
    Ytrk = 0;  % 0=unsplit 1=std 2=pref
    if Yarg >= 2
        [~, trk_i] = max([d_unsplit_path(Yarg), d_std_path(Yarg), d_pref_path(Yarg)]);
        Ytrk = trk_i - 1;
    end
    fprintf(route1_fid, '%5d %13.5e %13.5e %5d %13.5e %5d %13.5e %5d %12.6f %5d %12.6f %5d %5d %5d %6d %13.5e %5d\n', ...
        iter, Ymax, mn_denom, slots(k1), mx_lam, slots(k2), mx_fstd, slots(k3), ...
        mn_tauh, slots(k4), mx_wh, slots(k5), Yarg, Ytrk, n_clamp, mx_resid, slots(k6));
    %% H3 diagnostics (Phase A, 2026-06-10). H3-1: slot-2 consistency gap between the
    %% unsplit V_h and the expiry-weighted V_eff at the PT rows (raw-V vs V_eff
    %% inconsistency at the slot-2 surprise block). H3-2: state index driving Ymax
    %% within the dominant track. H3-3: min eligible-outside mass this sweep.
    e_h3 = expire_rate;
    dVeff_h3 = max(abs(V_eu15_h([i i2],2) - ((1-e_h3)*V_eu15_h_pref([i i2],2) + e_h3*V_eu15_h_std([i i2],2))));
    dVeff_h3 = max(dVeff_h3, max(abs(V_nms8_h([i i2],2) - ((1-e_h3)*V_nms8_h_pref([i i2],2) + e_h3*V_nms8_h_std([i i2],2)))));
    dVeff_h3 = max(dVeff_h3, max(abs(V_oth_h([i i2],2)  - ((1-e_h3)*V_oth_h_pref([i i2],2)  + e_h3*V_oth_h_std([i i2],2)))));
    if Ytrk == 0
        [~, h3row] = max(abs(Y(:, Yarg) - V(:, Yarg)));
    elseif Ytrk == 1
        [~, h3row] = max(abs(Y_h_std_2d(:, Yarg) - V_h_std_2d(:, Yarg)));
    else
        [~, h3row] = max(abs(Y_h_pref_2d(:, Yarg) - V_h_pref_2d(:, Yarg)));
    end
    h3state = mod(h3row - 1, J*N) + 1;
    fprintf(route1_fid, '  h3: dVeff2=%.4e  argmax_state=%d (stack-row %d)  min_Lelig=%.4e  Ymax_full=%.4e\n', ...
        dVeff_h3, h3state, h3row, min_Lelig_sweep, max(Ydiff_full));
    fflush(route1_fid);
    if (iter >= diag_iter_lo) && (iter <= diag_iter_hi)
        itf = fopen(fullfile(diag_dir, sprintf('route1_iter_%02d.log', iter)), 'w');
        fprintf(itf, '# iter=%d Ymax=%.6e Yarg_slot=%d dom_track=%d (0=unsplit 1=std 2=pref)\n', ...
            iter, Ymax, Yarg, Ytrk);
        fprintf(itf, '%5s %11s %11s %13s %13s %13s %13s %13s %13s %13s %13s %11s %11s %11s %5s %13s %13s\n', ...
            'slot','w_h_hat','w_l_hat','lambda','lam_raw','denom', ...
            'A_hat','C_hat','G_target','numer','resid','tau_h_std','f_std','f_l','clmp','Ydiff','B_hat');
        for s = slots
            numer_s = (diag_A_hat(s)+diag_B_hat(s)+diag_C_hat(s)) - diag_G_base_t(s);  % A+B+C - G_target
            fprintf(itf, '%5d %11.6f %11.6f %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %11.6f %11.6f %11.6f %5d %13.5e %13.5e\n', ...
                s, diag_w_h_hat(s), diag_w_l_hat(s), lambda_hat_path(s), diag_lambda_raw(s), diag_denom(s), ...
                diag_A_hat(s), diag_C_hat(s), diag_G_base_t(s), numer_s, diag_budget_resid(s), ...
                tau_h_t_path(s), f_std_path(s), f_l_path(s), diag_clamp_binding(s), Ydiff(s), diag_B_hat(s));
        end
        fclose(itf);
    end
    %% --- end Route 1 diagnostic dump ---

    iter = iter + 1;
end
if exist('route1_fid', 'var'); fclose(route1_fid); end

%% SWEEP-ENGINE INSERTION 3 of 3 (welfare dump, 2026-06-14): SAVE_WELFARE gated.
%% Writes the equilibrium objects needed to (a) evaluate the COPS online-appendix
%% household-welfare formula directly (cev_direct.m) and (b) reconstruct the
%% value-hat V_hat from these objects as a cross-check (reconstruct_vhat.m logic).
%% No effect on the iteration; default-off (env unset -> skip). Runs before the
%% sweep:done stop so all forward-block objects (mu, realwages, tau) are in scope.
sw_env = getenv('SAVE_WELFARE');
if ~isempty(sw_env)
    wel = struct();
    wel.i = i; wel.i2 = i2; wel.N = N; wel.J = J; wel.time = time;
    wel.beta = beta; wel.sigma = sigma; wel.expire_rate = expire_rate;
    %% consumption hats (after-tax real-wage hats) at PT, country-level 1xNcx time
    wel.rw_pref = realwages_h_pref; wel.rw_std = realwages_h_std; wel.rw_l = realwages_l;
    wel.rw_h = realwagesaux_h;  %% unsplit high-skill real-wage hat (for the V_h track)
    %% non-employed consumption hat (public-goods/lump-sum adjusted), N x time
    wel.aux_nonemp = aux_nonemp;
    %% counterfactual migration matrices (34x34xtime) for the split + unsplit tracks
    wel.mu_eu15_h_pref = mu_eu15_h_pref; wel.mu_nms8_h_pref = mu_nms8_h_pref; wel.mu_oth_h_pref = mu_oth_h_pref;
    wel.mu_eu15_h_std  = mu_eu15_h_std;  wel.mu_nms8_h_std  = mu_nms8_h_std;  wel.mu_oth_h_std  = mu_oth_h_std;
    wel.mu_eu15_h = mu_eu15_h; wel.mu_nms8_h = mu_nms8_h; wel.mu_oth_h = mu_oth_h;
    wel.mu_eu15_l = mu_eu15_l; wel.mu_nms8_l = mu_nms8_l; wel.mu_oth_l = mu_oth_l;
    %% baseline migration matrices
    wel.mub_eu15_h = mu_baseline_eu15_h; wel.mub_nms8_h = mu_baseline_nms8_h; wel.mub_oth_h = mu_baseline_oth_h;
    wel.mub_eu15_l = mu_baseline_eu15_l; wel.mub_nms8_l = mu_baseline_nms8_l; wel.mub_oth_l = mu_baseline_oth_l;
    %% migration-cost hats (== 1 under Approach 2, shock_val=1; saved for generality)
    wel.tau_h = tau_mig_hat_eu15_h; wel.tau_nms8_h = tau_mig_hat_nms8_h; wel.tau_oth_h = tau_mig_hat_oth_h;
    wel.tau_eu15_l = tau_mig_hat_eu15_l; wel.tau_nms8_l = tau_mig_hat_nms8_l; wel.tau_oth_l = tau_mig_hat_oth_l;
    %% solved value-hat tracks (input state), for reconstruction cross-check
    wel.V_eu15_h_pref = V_eu15_h_pref; wel.V_nms8_h_pref = V_nms8_h_pref; wel.V_oth_h_pref = V_oth_h_pref;
    wel.V_eu15_h_std  = V_eu15_h_std;  wel.V_nms8_h_std  = V_nms8_h_std;  wel.V_oth_h_std  = V_oth_h_std;
    wel.V_eu15_l = V_eu15_l; wel.V_nms8_l = V_nms8_l; wel.V_oth_l = V_oth_l;
    %% engine backward outputs (pre-damping Ynew), ground truth for reconstruction debug
    wel.Yeng_eu15_h_pref = squeeze(Y_eu15_h_pref); wel.Yeng_nms8_h_pref = squeeze(Y_nms8_h_pref); wel.Yeng_oth_h_pref = squeeze(Y_oth_h_pref);
    wel.Yeng_eu15_h_std  = squeeze(Y_eu15_h_std);  wel.Yeng_nms8_h_std  = squeeze(Y_nms8_h_std);  wel.Yeng_oth_h_std  = squeeze(Y_oth_h_std);
    wel.Yeng_eu15_l = squeeze(Y_eu15_l); wel.Yeng_nms8_l = squeeze(Y_nms8_l); wel.Yeng_oth_l = squeeze(Y_oth_l);
    wel.Yeng_eu15_h = squeeze(Y_eu15_h); wel.Yeng_nms8_h = squeeze(Y_nms8_h); wel.Yeng_oth_h = squeeze(Y_oth_h);
    save('-v7', sw_env, 'wel');
    fprintf('[SAVE_WELFARE] dumped PT welfare objects to %s\n', sw_env); fflush(stdout);
end

%% SWEEP-ENGINE INSERTION 2 of 2: in sweep mode, stop here — the welfare /
%% reporting tail below is not needed per sweep evaluation. The undamped sweep
%% output lives in Ynew (stacked unsplit) and Ynew_h_std / Ynew_h_pref; the
%% raised flag is caught by sweep_phi.m, whose workspace retains all variables.
if exist('SWEEP_MODE', 'var')
    error('sweep:done', 'sweep complete');
end


%Computing welfare and migration outcomes
%
% NOTE (Approach 2): the unsplit V_*_h is iterated using the pre-tax wage hat
% rwage_h (line 940; no LS subtraction, no pref wedge) and the unsplit
% migration matrix mu_*_h. The Welfare/Welfare_results below therefore mix
% pre-tax wages with post-policy mu — they are the legacy COPS object, not
% the policy welfare under the preferential-tax scheme. The Approach 2
% policy welfare is the std/pref split CEV table built below from
% V_*_h_std and V_*_h_pref (see lines starting "WELFARE: std/pref split").

Welfare=[-100*(1-beta).*sigma.*log(V_eu15_h(:,2)) -100*(1-beta).*sigma.*log(V_eu15_l(:,2)) -100*(1-beta).*sigma.*log(V_nms8_h(:,2)) -100*(1-beta).*sigma.*log(V_nms8_l(:,2))];
share_eu15_h=L0_eu15_h./sum(L0_eu15_h);
share_eu15_l=L0_eu15_l./sum(L0_eu15_l);
share_nms8_h=L0_nms8_h./sum(L0_nms8_h);
share_nms8_l=L0_nms8_l./sum(L0_nms8_l);
shares=[share_eu15_h share_eu15_l share_nms8_h share_nms8_l ];
Welfare_agg=sum(Welfare.*shares)';

Welfare=[-100*(1-beta).*sigma.*log(V_eu15_h(:,2)) -100*(1-beta).*sigma.*log(V_eu15_l(:,2)) -100*(1-beta).*sigma.*log(V_nms8_h(:,2)) -100*(1-beta).*sigma.*log(V_nms8_l(:,2))];

L0=[L0_eu15_h L0_eu15_l L0_nms8_h L0_nms8_l];
W=L0.*Welfare;

share_eu15_h= L0_eu15_h(1)+L0_eu15_h(2)+L0_eu15_h(5)+L0_eu15_h(6)+L0_eu15_h(8)+L0_eu15_h(9)+L0_eu15_h(10)+L0_eu15_h(12)+L0_eu15_h(16)+L0_eu15_h(17)+....
              L0_nms8_h(1)+ L0_nms8_h(2)+L0_nms8_h(5)+L0_nms8_h(6)+L0_nms8_h(8)+L0_nms8_h(9)+L0_nms8_h(10)+L0_nms8_h(12)+L0_nms8_h(16)+L0_nms8_h(17)+...
              L0_eu15_h(1+17)+L0_eu15_h(2+17)+L0_eu15_h(5+17)+L0_eu15_h(6+17)+L0_eu15_h(8+17)+L0_eu15_h(9+17)+L0_eu15_h(10+17)+L0_eu15_h(12+17)+L0_eu15_h(16+17)+L0_eu15_h(17+17)+.... 
              L0_nms8_h(1+17)+ L0_nms8_h(2+17)+L0_nms8_h(5+17)+L0_nms8_h(6+17)+L0_nms8_h(8+17)+L0_nms8_h(9+17)+L0_nms8_h(10+17)+L0_nms8_h(12+17)+L0_nms8_h(16+17)+L0_nms8_h(17+17);


share_eu15_l=L0_eu15_l(1)+ L0_eu15_l(2)+L0_eu15_l(5)+L0_eu15_l(6)+L0_eu15_l(8)+L0_eu15_l(9)+L0_eu15_l(10)+L0_eu15_l(12)+L0_eu15_l(16)+L0_eu15_l(17)+....
             L0_nms8_l(1)+ L0_nms8_l(2)+L0_nms8_l(5)+L0_nms8_l(6)+L0_nms8_l(8)+L0_nms8_l(9)+L0_nms8_l(10)+L0_nms8_l(12)+L0_nms8_l(16)+L0_nms8_l(17)+...
             L0_eu15_l(1+17)+ L0_eu15_l(2+17)+L0_eu15_l(5+17)+L0_eu15_l(6+17)+L0_eu15_l(8+17)+L0_eu15_l(9+17)+L0_eu15_l(10+17)+L0_eu15_l(12+17)+L0_eu15_l(16+17)+L0_eu15_l(17+17)+....
             L0_nms8_l(1+17)+ L0_nms8_l(2+17)+L0_nms8_l(5+17)+L0_nms8_l(6+17)+L0_nms8_l(8+17)+L0_nms8_l(9+17)+L0_nms8_l(10+17)+L0_nms8_l(12+17)+L0_nms8_l(16+17)+L0_nms8_l(17+17);

share_nms8_h=L0_eu15_h(3)+ L0_eu15_h(4)+ L0_eu15_h(7)+L0_eu15_h(11)+L0_eu15_h(13)+L0_eu15_h(14)+L0_eu15_h(15)+....
             L0_nms8_h(3)+ L0_nms8_h(4)+ L0_nms8_h(7)+L0_nms8_h(11)+L0_nms8_h(13)+L0_nms8_h(14)+L0_nms8_h(15)+...
             L0_eu15_h(3+17)+ L0_eu15_h(4+17)+ L0_eu15_h(7+17)+L0_eu15_h(11+17)+L0_eu15_h(13+17)+L0_eu15_h(14+17)+L0_eu15_h(15+17)+....
             L0_nms8_h(3+17)+ L0_nms8_h(4+17)+ L0_nms8_h(7+17)+L0_nms8_h(11+17)+L0_nms8_h(13+17)+L0_nms8_h(14+17)+L0_nms8_h(15+17);

share_nms8_l=L0_eu15_l(3)+ L0_eu15_l(4)+ L0_eu15_l(7)+L0_eu15_l(11)+L0_eu15_l(13)+L0_eu15_l(14)+L0_eu15_l(15)+....
             L0_nms8_l(3)+ L0_nms8_l(4)+ L0_nms8_l(7)+L0_nms8_l(11)+L0_nms8_l(13)+L0_nms8_l(14)+L0_nms8_l(15)+...
             L0_eu15_l(3+17)+ L0_eu15_l(4+17)+ L0_eu15_l(7+17)+L0_eu15_l(11+17)+L0_eu15_l(13+17)+L0_eu15_l(14+17)+L0_eu15_l(15+17)+....
             L0_nms8_l(3+17)+ L0_nms8_l(4+17)+ L0_nms8_l(7+17)+L0_nms8_l(11+17)+L0_nms8_l(13+17)+L0_nms8_l(14+17)+L0_nms8_l(15+17);

Welfare_eu15_h=(W(1,1)+W(2,1)+W(5,1)+W(6,1)+W(8,1)+W(9,1)+W(10,1)+W(12,1)+W(16,1)+W(17,1)+...
                W(1,3)+W(2,3)+W(5,3)+W(6,3)+W(8,3)+W(9,3)+W(10,3)+W(12,3)+W(16,3)+W(17,3)+...
                W(1+17,1)+W(2+17,1)+W(5+17,1)+W(6+17,1)+W(8+17,1)+W(9+17,1)+W(10+17,1)+W(12+17,1)+W(16+17,1)+W(17+17,1)+...
                W(1+17,3)+W(2+17,3)+W(5+17,3)+W(6+17,3)+W(8+17,3)+W(9+17,3)+W(10+17,3)+W(12+17,3)+W(16+17,3)+W(17+17,3))./share_eu15_h;

Welfare_eu15_l=(W(1,2)+W(2,2)+W(5,2)+W(6,2)+W(8,2)+W(9,2)+W(10,2)+W(12,2)+W(16,2)+W(17,2)+...
                W(1,4)+W(2,4)+W(5,4)+W(6,4)+W(8,4)+W(9,4)+W(10,4)+W(12,4)+W(16,4)+W(17,4)+...
                W(1+17,2)+W(2+17,2)+W(5+17,2)+W(6+17,2)+W(8+17,2)+W(9+17,2)+W(10+17,2)+W(12+17,2)+W(16+17,2)+W(17+17,2)+...
                W(1+17,4)+W(2+17,4)+W(5+17,4)+W(6+17,4)+W(8+17,4)+W(9+17,4)+W(10+17,4)+W(12+17,4)+W(16+17,4)+W(17+17,4))./share_eu15_l;

Welfare_nms8_h=(W(3,1)+W(4,1)+W(7,1)+W(11,1)+W(13,1)+W(14,1)+W(15,1)+...
                W(3,3)+W(4,3)+W(7,3)+W(11,3)+W(13,3)+W(14,3)+W(15,3)+...
                W(3+17,1)+W(4+17,1)+W(7+17,1)+W(11+17,1)+W(13+17,1)+W(14+17,1)+W(15+17,1)+...
                W(3+17,3)+W(4+17,3)+W(7+17,3)+W(11+17,3)+W(13+17,3)+W(14+17,3)+W(15+17,3))./share_nms8_h;

Welfare_nms8_l=(W(3,2)+W(4,2)+W(7,2)+W(11,2)+W(13,2)+W(14,2)+W(15,2)+...
                W(3,4)+W(4,4)+W(7,4)+W(11,4)+W(13,4)+W(14,4)+W(15,4)+...
                W(3+17,2)+W(4+17,2)+W(7+17,2)+W(11+17,2)+W(13+17,2)+W(14+17,2)+W(15+17,2)+...
                W(3+17,4)+W(4+17,4)+W(7+17,4)+W(11+17,4)+W(13+17,4)+W(14+17,4)+W(15+17,4))./share_nms8_l;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

Welfare_eu15=(Welfare_eu15_h*share_eu15_h+Welfare_eu15_l*share_eu15_l)./(share_eu15_h+share_eu15_l);
Welfare_nms8=(Welfare_nms8_h*share_nms8_h+Welfare_nms8_l*share_nms8_l)./(share_nms8_h+share_nms8_l);

Welfare_europe=(Welfare_eu15*(share_eu15_h+share_eu15_l)+Welfare_nms8*(share_nms8_h+share_nms8_l))./(share_eu15_h+share_eu15_l+share_nms8_h+share_nms8_l);

Welfare_results=[Welfare_eu15_h;Welfare_eu15_l;Welfare_nms8_h;Welfare_nms8_l;Welfare_eu15;Welfare_nms8;Welfare_europe]

% ==== WELFARE: std/pref split for target country ====
% --- (a) PT h-skill incumbents (std pool): CEV using V_std at array slot 2 ---
Welfare_PT_eu15_h_std_s1 = -100*(1-beta)*sigma*log(V_eu15_h_std(i,  2));
Welfare_PT_eu15_h_std_s2 = -100*(1-beta)*sigma*log(V_eu15_h_std(i2, 2));
Welfare_PT_nms8_h_std_s1 = -100*(1-beta)*sigma*log(V_nms8_h_std(i,  2));
Welfare_PT_nms8_h_std_s2 = -100*(1-beta)*sigma*log(V_nms8_h_std(i2, 2));
Welfare_PT_oth_h_std_s1  = -100*(1-beta)*sigma*log(V_oth_h_std(i,   2));
Welfare_PT_oth_h_std_s2  = -100*(1-beta)*sigma*log(V_oth_h_std(i2,  2));

% --- (b) PT h-skill arrivals (pref pool): CEV using V_pref at array slot 2 ---
Welfare_PT_eu15_h_pref_s1 = -100*(1-beta)*sigma*log(V_eu15_h_pref(i,  2));
Welfare_PT_eu15_h_pref_s2 = -100*(1-beta)*sigma*log(V_eu15_h_pref(i2, 2));
Welfare_PT_nms8_h_pref_s1 = -100*(1-beta)*sigma*log(V_nms8_h_pref(i,  2));
Welfare_PT_nms8_h_pref_s2 = -100*(1-beta)*sigma*log(V_nms8_h_pref(i2, 2));
Welfare_PT_oth_h_pref_s1  = -100*(1-beta)*sigma*log(V_oth_h_pref(i,   2));
Welfare_PT_oth_h_pref_s2  = -100*(1-beta)*sigma*log(V_oth_h_pref(i2,  2));

% --- (c) PT l-skill ---
Welfare_PT_eu15_l_s1 = -100*(1-beta)*sigma*log(V_eu15_l(i,  2));
Welfare_PT_eu15_l_s2 = -100*(1-beta)*sigma*log(V_eu15_l(i2, 2));
Welfare_PT_nms8_l_s1 = -100*(1-beta)*sigma*log(V_nms8_l(i,  2));
Welfare_PT_nms8_l_s2 = -100*(1-beta)*sigma*log(V_nms8_l(i2, 2));

% --- (d) Aggregate welfare for PT workers ---
L0_PT_eu15_h = L0_eu15_h(i)  + L0_eu15_h(i2);
L0_PT_nms8_h = L0_nms8_h(i)  + L0_nms8_h(i2);
L0_PT_oth_h  = L0_oth_h(i)   + L0_oth_h(i2);
L0_PT_l      = L0_eu15_l(i)  + L0_eu15_l(i2) + ...
               L0_nms8_l(i)  + L0_nms8_l(i2) + ...
               L0_oth_l(i)   + L0_oth_l(i2);
L0_PT_h      = L0_PT_eu15_h + L0_PT_nms8_h + L0_PT_oth_h;

Welfare_PT_h_std = ( L0_eu15_h(i)  * Welfare_PT_eu15_h_std_s1 + ...
                     L0_eu15_h(i2) * Welfare_PT_eu15_h_std_s2 + ...
                     L0_nms8_h(i)  * Welfare_PT_nms8_h_std_s1 + ...
                     L0_nms8_h(i2) * Welfare_PT_nms8_h_std_s2 + ...
                     L0_oth_h(i)   * Welfare_PT_oth_h_std_s1 + ...
                     L0_oth_h(i2)  * Welfare_PT_oth_h_std_s2 ) ./ L0_PT_h;

Welfare_PT_l     = ( L0_eu15_l(i)  * Welfare_PT_eu15_l_s1 + ...
                     L0_eu15_l(i2) * Welfare_PT_eu15_l_s2 + ...
                     L0_nms8_l(i)  * Welfare_PT_nms8_l_s1 + ...
                     L0_nms8_l(i2) * Welfare_PT_nms8_l_s2 + ...
                     L0_oth_l(i)   * (-100*(1-beta)*sigma*log(V_oth_l(i,  2))) + ...
                     L0_oth_l(i2)  * (-100*(1-beta)*sigma*log(V_oth_l(i2, 2))) ) ./ L0_PT_l;

Welfare_PT_agg   = (Welfare_PT_h_std * L0_PT_h + Welfare_PT_l * L0_PT_l) ./ ...
                   (L0_PT_h + L0_PT_l);

% --- Display ---
fprintf('\n=== Welfare: PT workers (Direct preferential tax) ===\n');
fprintf('PT h-skill incumbents (std rate)  : sector 1 = %.4f, sector 2 = %.4f\n', ...
    Welfare_PT_eu15_h_std_s1, Welfare_PT_eu15_h_std_s2);
fprintf('PT h-skill option value (pref rate): sector 1 = %.4f, sector 2 = %.4f\n', ...
    Welfare_PT_eu15_h_pref_s1, Welfare_PT_eu15_h_pref_s2);
fprintf('PT l-skill (standard rate)         : sector 1 = %.4f, sector 2 = %.4f\n', ...
    Welfare_PT_eu15_l_s1, Welfare_PT_eu15_l_s2);
fprintf('PT h-skill per capita (std, agg)   : %.4f\n', Welfare_PT_h_std);
fprintf('PT l-skill per capita (std, agg)   : %.4f\n', Welfare_PT_l);
fprintf('PT aggregate welfare               : %.4f\n', Welfare_PT_agg);

Welfare_PT_results = [Welfare_PT_h_std; Welfare_PT_l; Welfare_PT_agg] %#ok<NOPTS>

%% ============================================================================
%% Aggregate CEV table by skill x pool x employment status, weighted by slot-2
%% pool sizes. l-skill has no pool split (always faces standard rate).
%% ============================================================================
Welfare_PT_oth_l_s1 = -100*(1-beta)*sigma*log(V_oth_l(i,  2));
Welfare_PT_oth_l_s2 = -100*(1-beta)*sigma*log(V_oth_l(i2, 2));

L_PT_h_std_emp_2  = Ldyn_eu15_h_std(1, i, 2)  + Ldyn_nms8_h_std(1, i, 2)  + Ldyn_oth_h_std(1, i, 2);
L_PT_h_pref_emp_2 = Ldyn_eu15_h_pref(1, i, 2) + Ldyn_nms8_h_pref(1, i, 2) + Ldyn_oth_h_pref(1, i, 2);
L_PT_h_std_neemp_2  = Ldyn_eu15_h_std(1, i2, 2)  + Ldyn_nms8_h_std(1, i2, 2)  + Ldyn_oth_h_std(1, i2, 2);
L_PT_h_pref_neemp_2 = Ldyn_eu15_h_pref(1, i2, 2) + Ldyn_nms8_h_pref(1, i2, 2) + Ldyn_oth_h_pref(1, i2, 2);
L_PT_l_emp_2   = Ldyn_eu15_l(1, i, 2)  + Ldyn_nms8_l(1, i, 2)  + Ldyn_oth_l(1, i, 2);
L_PT_l_neemp_2 = Ldyn_eu15_l(1, i2, 2) + Ldyn_nms8_l(1, i2, 2) + Ldyn_oth_l(1, i2, 2);

CEV_h_pref_emp = (Ldyn_eu15_h_pref(1, i, 2) * Welfare_PT_eu15_h_pref_s1 + ...
                  Ldyn_nms8_h_pref(1, i, 2) * Welfare_PT_nms8_h_pref_s1 + ...
                  Ldyn_oth_h_pref(1, i, 2)  * Welfare_PT_oth_h_pref_s1) / max(L_PT_h_pref_emp_2, eps);
CEV_h_pref_neemp = (Ldyn_eu15_h_pref(1, i2, 2) * Welfare_PT_eu15_h_pref_s2 + ...
                    Ldyn_nms8_h_pref(1, i2, 2) * Welfare_PT_nms8_h_pref_s2 + ...
                    Ldyn_oth_h_pref(1, i2, 2)  * Welfare_PT_oth_h_pref_s2) / max(L_PT_h_pref_neemp_2, eps);
CEV_h_std_emp = (Ldyn_eu15_h_std(1, i, 2) * Welfare_PT_eu15_h_std_s1 + ...
                 Ldyn_nms8_h_std(1, i, 2) * Welfare_PT_nms8_h_std_s1 + ...
                 Ldyn_oth_h_std(1, i, 2)  * Welfare_PT_oth_h_std_s1) / max(L_PT_h_std_emp_2, eps);
CEV_h_std_neemp = (Ldyn_eu15_h_std(1, i2, 2) * Welfare_PT_eu15_h_std_s2 + ...
                   Ldyn_nms8_h_std(1, i2, 2) * Welfare_PT_nms8_h_std_s2 + ...
                   Ldyn_oth_h_std(1, i2, 2)  * Welfare_PT_oth_h_std_s2) / max(L_PT_h_std_neemp_2, eps);
CEV_l_emp = (Ldyn_eu15_l(1, i, 2) * Welfare_PT_eu15_l_s1 + ...
             Ldyn_nms8_l(1, i, 2) * Welfare_PT_nms8_l_s1 + ...
             Ldyn_oth_l(1, i, 2)  * Welfare_PT_oth_l_s1) / max(L_PT_l_emp_2, eps);
CEV_l_neemp = (Ldyn_eu15_l(1, i2, 2) * Welfare_PT_eu15_l_s2 + ...
               Ldyn_nms8_l(1, i2, 2) * Welfare_PT_nms8_l_s2 + ...
               Ldyn_oth_l(1, i2, 2)  * Welfare_PT_oth_l_s2) / max(L_PT_l_neemp_2, eps);

L_PT_h_total_2 = L_PT_h_std_emp_2 + L_PT_h_pref_emp_2 + L_PT_h_std_neemp_2 + L_PT_h_pref_neemp_2;
L_PT_l_total_2 = L_PT_l_emp_2 + L_PT_l_neemp_2;
CEV_h_agg = (L_PT_h_std_emp_2 * CEV_h_std_emp + L_PT_h_pref_emp_2 * CEV_h_pref_emp + ...
             L_PT_h_std_neemp_2 * CEV_h_std_neemp + L_PT_h_pref_neemp_2 * CEV_h_pref_neemp) / max(L_PT_h_total_2, eps);
CEV_l_agg = (L_PT_l_emp_2 * CEV_l_emp + L_PT_l_neemp_2 * CEV_l_neemp) / max(L_PT_l_total_2, eps);
CEV_PT_agg = (L_PT_h_total_2 * CEV_h_agg + L_PT_l_total_2 * CEV_l_agg) / max(L_PT_h_total_2 + L_PT_l_total_2, eps);

fprintf('\n=== PT CEV summary (%% consumption equivalent at slot 2) ===\n');
fprintf('%-30s %12s %12s\n', 'Group', 'Employed', 'Non-empl');
fprintf('%-30s %12.4f %12.4f\n', 'h-skill (pref rate)', CEV_h_pref_emp, CEV_h_pref_neemp);
fprintf('%-30s %12.4f %12.4f\n', 'h-skill (std rate)',  CEV_h_std_emp,  CEV_h_std_neemp);
fprintf('%-30s %12.4f %12.4f\n', 'l-skill (std rate)',  CEV_l_emp,      CEV_l_neemp);
fprintf('%-30s %12.4f\n', 'h-skill aggregate', CEV_h_agg);
fprintf('%-30s %12.4f\n', 'l-skill aggregate', CEV_l_agg);
fprintf('%-30s %12.4f\n', 'PT aggregate', CEV_PT_agg);
fprintf('Slot-2 pool sizes (h std emp, h pref emp, h std non, h pref non, l emp, l non):\n');
fprintf('  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e\n', L_PT_h_std_emp_2, L_PT_h_pref_emp_2, L_PT_h_std_neemp_2, L_PT_h_pref_neemp_2, L_PT_l_emp_2, L_PT_l_neemp_2);

%% ============================================================================
%% Per-period diagnostic dump to pt_diagnostic.log. Includes: tax rates,
%% wage hats, labor stocks (counter and same-t baseline), G actual vs G base,
%% value functions at PT (sector 1 emp + sector 2 non-emp) for std/pref/l.
%% ============================================================================
diag_t_grid = [2:31, 51, 76, 101, 126, 151, 175, 199];
if time > 200
    %% extend the dump grid into the padded horizon so the tail is visible
    diag_t_grid = [diag_t_grid, round(linspace(225, time-1, 12))];
end
diag_t_grid = unique(diag_t_grid(diag_t_grid <= time));
diag_fid = fopen('pt_diagnostic.log', 'w');
fprintf(diag_fid, '# PT per-period diagnostic dump. tau_P_pref=%.4f, tau_h_base=%.4f, tau_l_base=%.4f.\n', tau_P_pref, tau_h_base, tau_l_base);
fprintf(diag_fid, '# Slot 1 = pre-shock (t=0). Slot 2 = first counterfactual period (t=1).\n');
fprintf(diag_fid, '# Hats are counter/baseline at the same slot.\n\n');
fprintf(diag_fid, '## Tax rates and budget (G in units of baseline W_h*L_h at PT)\n');
fprintf(diag_fid, '# tau_h_cand = shadow proportional rate (closed-form, plan sec 2.8). Diagnostic only under LS.\n');
fprintf(diag_fid, '%5s %12s %12s %12s %12s %12s %12s %12s %12s\n', 'slot', 'tau_h_t', 'tau_h_cand', 'f_std', 'T_hat', 'rho_t', 'G_actual', 'G_baseline', 'G_diff');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %12.6f %12.6f %12.6f %12.6e %12.6f %12.6f %12.6f %12.6e\n', ...
        k, tau_h_t_path(k), tau_h_t_candidate_path(k), f_std_path(k), T_hat_path(k), diag_rho_t(k), G_actual_path(k), diag_G_base_t(k), G_actual_path(k) - diag_G_base_t(k));
end

%% Full-path closure-variable trajectory (every slot, not just diag_t_grid).
%% Under lambda-rescale closure (paper Sec 4): lambda_hat = lambda_PT(t)/lambda_PT_baseline
%% is the live instrument; rates fall when lambda_hat>1 (surplus), rise when <1 (deficit).
%% f_std and f_l are the resulting flow-utility wedges at PT sector 1.
%% tau_h_t and tau_l_t are the implied average tax rates.
%% T_hat and tau_h_cand are diagnostics from the deprecated LS / pure-tau_h closures.
fprintf(diag_fid, '\n## Closure variable lambda_hat + implied wedges (every t)\n');
fprintf(diag_fid, '# lambda_hat: paper Sec 4 closure variable. f_std/f_l: wedges on PT sector-1 realwages.\n');
fprintf(diag_fid, '# tau_h_t/tau_l_t: implied avg rates. T_hat/tau_h_cand: diagnostic comparisons.\n');
fprintf(diag_fid, '%5s %12s %12s %12s %12s %12s %14s %12s\n', ...
    'slot', 'lambda_hat', 'tau_h_t', 'tau_l_t', 'f_std', 'f_l', 'T_hat_diag', 'tau_h_cand');
for k = 2:(time-1)
    fprintf(diag_fid, '%5d %12.6f %12.6f %12.6f %12.6f %12.6f %14.6e %12.6f\n', ...
        k, lambda_hat_path(k), tau_h_t_path(k), tau_l_t_path(k), f_std_path(k), f_l_path(k), T_hat_path(k), tau_h_t_candidate_path(k));
end

fprintf(diag_fid, '\n## Wage hats at PT (counter/baseline, same-t)\n');
fprintf(diag_fid, '%5s %12s %12s %12s %12s\n', 'slot', 'w_h_hat', 'w_l_hat', 'W_h_base', 'W_l_base');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %12.6f %12.6f %12.6f %12.6f\n', k, diag_w_h_hat(k), diag_w_l_hat(k), diag_W_h_base(k), diag_W_l_base(k));
end

fprintf(diag_fid, '\n## Labor stocks at target country (i), employed (sector 1, i=16). Counter = std + pref + (= L_h total).\n');
fprintf(diag_fid, '%5s %14s %14s %14s %14s %14s %14s\n', 'slot', 'L_h_std', 'L_h_pref', 'L_h_total_cf', 'L_h_baseline', 'L_l_cf', 'L_l_baseline');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e\n', ...
        k, diag_L_h_std(k), diag_L_h_pref(k), diag_L_h_std(k) + diag_L_h_pref(k), diag_L_h_base(k), diag_L_l_c(k), diag_L_l_base(k));
end

fprintf(diag_fid, '\n## Labor stocks at target country (i), non-employed (sector 2, i2=33).\n');
fprintf(diag_fid, '%5s %14s %14s %14s %14s\n', 'slot', 'L_h_std_neemp', 'L_h_pref_neemp', 'L_h_neemp_cf', 'L_l_neemp_cf');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %14.4e %14.4e %14.4e %14.4e\n', ...
        k, diag_L_h_std_neemp(k), diag_L_h_pref_neemp(k), diag_L_h_std_neemp(k) + diag_L_h_pref_neemp(k), diag_L_l_c_neemp(k));
end

fprintf(diag_fid, '\n## Value functions V_hat at target country (i), employed (i=16). EU-15 bloc (representative).\n');
fprintf(diag_fid, '%5s %14s %14s %14s\n', 'slot', 'V_eu15_h_std', 'V_eu15_h_pref', 'V_eu15_l');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %14.6f %14.6f %14.6f\n', k, V_eu15_h_std(i, k), V_eu15_h_pref(i, k), V_eu15_l(i, k));
end

fprintf(diag_fid, '\n## Value functions V_hat at target country (i), non-employed (i2=33). EU-15 bloc.\n');
fprintf(diag_fid, '%5s %14s %14s %14s\n', 'slot', 'V_eu15_h_std', 'V_eu15_h_pref', 'V_eu15_l');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %14.6f %14.6f %14.6f\n', k, V_eu15_h_std(i2, k), V_eu15_h_pref(i2, k), V_eu15_l(i2, k));
end

fprintf(diag_fid, '\n## V_hat at target country, employed (i), all blocs std/pref.\n');
fprintf(diag_fid, '%5s %12s %12s %12s %12s %12s %12s\n', 'slot', 'eu15_std', 'eu15_pref', 'nms8_std', 'nms8_pref', 'oth_std', 'oth_pref');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', k, ...
        V_eu15_h_std(i, k), V_eu15_h_pref(i, k), ...
        V_nms8_h_std(i, k), V_nms8_h_pref(i, k), ...
        V_oth_h_std(i, k),  V_oth_h_pref(i, k));
end

%% Std-pool outflow rate at PT: per-period probability mass leaving the std pool at PT
%% for any non-PT destination state. Plan §2.12 diagnostic.
%% Compared against baseline EU-15 outflow rate (mu_baseline_eu15_h same row sum).
fprintf(diag_fid, '\n## Std-pool outflow rate at target country, employed origin (i=16). Sum over d not in [i, i2] of mu_b_h_std(i, d, t).\n');
fprintf(diag_fid, '%5s %14s %14s %14s %14s\n', 'slot', 'eu15_std', 'nms8_std', 'oth_std', 'base_eu15');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %14.6e %14.6e %14.6e %14.6e\n', k, ...
        diag_outflow_std_eu15_emp(k), diag_outflow_std_nms8_emp(k), ...
        diag_outflow_std_oth_emp(k), diag_outflow_base_emp(k));
end

fprintf(diag_fid, '\n## Std-pool outflow rate at target country, non-employed origin (i2=33). Sum over d not in [i, i2] of mu_b_h_std(i2, d, t).\n');
fprintf(diag_fid, '%5s %14s %14s %14s %14s\n', 'slot', 'eu15_std', 'nms8_std', 'oth_std', 'base_eu15');
for k = diag_t_grid
    fprintf(diag_fid, '%5d %14.6e %14.6e %14.6e %14.6e\n', k, ...
        diag_outflow_std_eu15_neemp(k), diag_outflow_std_nms8_neemp(k), ...
        diag_outflow_std_oth_neemp(k), diag_outflow_base_neemp(k));
end

fclose(diag_fid);
fprintf('Per-period diagnostic dump written to pt_diagnostic.log\n');

L_region_nms8_h=zeros(time,1,1);

for t = 1:time
        L_stock_nms8(:, :, t) = permute(Ldyn_nms8_h(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(Ldyn_eu15_h(:, :, t), [2, 1, 3]);
        num_L_region_nms8_h(t,1,1)=L_stock_nms8(1,:,t)+L_stock_nms8(2,:,t)+L_stock_nms8(5,:,t)+L_stock_nms8(6,:,t)...
            +L_stock_nms8(8,:,t)+L_stock_nms8(9,:,t)+L_stock_nms8(10,:,t)+L_stock_nms8(12,:,t)+L_stock_nms8(16,:,t)...
            +L_stock_nms8(17,:,t)+...
            L_stock_nms8(1+17,:,t)+L_stock_nms8(2+17,:,t)+L_stock_nms8(5+17,:,t)+L_stock_nms8(6+17,:,t)...
            +L_stock_nms8(8+17,:,t)+L_stock_nms8(9+17,:,t)+L_stock_nms8(10+17,:,t)+L_stock_nms8(12+17,:,t)+L_stock_nms8(16+17,:,t)...
            +L_stock_nms8(17+17,:,t);
        L_region_nms8_h(t,1,1)=num_L_region_nms8_h(t,1,1);
        L_region_nms8_h(t,1,1)=L_region_nms8_h(t,1,1);
       
end
    
L_region_nms8_baseline_h=zeros(time,1,1);

for t = 1:time
        L_stock_nms8(:, :, t) = permute(series_Ldyn_nms8_h(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(series_Ldyn_eu15_h(:, :, t), [2, 1, 3]);
       num_L_region_nms8_h(t,1,1)=L_stock_nms8(1,:,t)+L_stock_nms8(2,:,t)+L_stock_nms8(5,:,t)+L_stock_nms8(6,:,t)...
            +L_stock_nms8(8,:,t)+L_stock_nms8(9,:,t)+L_stock_nms8(10,:,t)+L_stock_nms8(12,:,t)+L_stock_nms8(16,:,t)...
            +L_stock_nms8(17,:,t)+...
            L_stock_nms8(1+17,:,t)+L_stock_nms8(2+17,:,t)+L_stock_nms8(5+17,:,t)+L_stock_nms8(6+17,:,t)...
            +L_stock_nms8(8+17,:,t)+L_stock_nms8(9+17,:,t)+L_stock_nms8(10+17,:,t)+L_stock_nms8(12+17,:,t)+L_stock_nms8(16+17,:,t)...
            +L_stock_nms8(17+17,:,t);
        L_region_nms8_baseline_h(t,1,1)=num_L_region_nms8_h(t,1,1);
        L_region_nms8_baseline_h(t,1,1)=L_region_nms8_baseline_h(t,1,1);
       
end

L_region_nms8_l=zeros(time,1,1);
 
for t = 1:time
        L_stock_nms8(:, :, t) = permute(Ldyn_nms8_l(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(Ldyn_eu15_l(:, :, t), [2, 1, 3]);
        num_L_region_nms8_l(t,1,1)=L_stock_nms8(1,:,t)+L_stock_nms8(2,:,t)+L_stock_nms8(5,:,t)+L_stock_nms8(6,:,t)...
            +L_stock_nms8(8,:,t)+L_stock_nms8(9,:,t)+L_stock_nms8(10,:,t)+L_stock_nms8(12,:,t)+L_stock_nms8(16,:,t)...
            +L_stock_nms8(17,:,t)+...
            L_stock_nms8(1+17,:,t)+L_stock_nms8(2+17,:,t)+L_stock_nms8(5+17,:,t)+L_stock_nms8(6+17,:,t)...
            +L_stock_nms8(8+17,:,t)+L_stock_nms8(9+17,:,t)+L_stock_nms8(10+17,:,t)+L_stock_nms8(12+17,:,t)+L_stock_nms8(16+17,:,t)...
            +L_stock_nms8(17+17,:,t);
        L_region_nms8_l(t,1,1)=num_L_region_nms8_l(t,1,1);
        L_region_nms8_l(t,1,1)=L_region_nms8_l(t,1,1);
       
end
    
L_region_nms8_baseline_l=zeros(time,1,1);
 
for t = 1:time
        L_stock_nms8(:, :, t) = permute(series_Ldyn_nms8_l(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(series_Ldyn_eu15_l(:, :, t), [2, 1, 3]);
       num_L_region_nms8_l(t,1,1)=L_stock_nms8(1,:,t)+L_stock_nms8(2,:,t)+L_stock_nms8(5,:,t)+L_stock_nms8(6,:,t)...
            +L_stock_nms8(8,:,t)+L_stock_nms8(9,:,t)+L_stock_nms8(10,:,t)+L_stock_nms8(12,:,t)+L_stock_nms8(16,:,t)...
            +L_stock_nms8(17,:,t)+...
            L_stock_nms8(1+17,:,t)+L_stock_nms8(2+17,:,t)+L_stock_nms8(5+17,:,t)+L_stock_nms8(6+17,:,t)...
            +L_stock_nms8(8+17,:,t)+L_stock_nms8(9+17,:,t)+L_stock_nms8(10+17,:,t)+L_stock_nms8(12+17,:,t)+L_stock_nms8(16+17,:,t)...
            +L_stock_nms8(17+17,:,t);
        L_region_nms8_baseline_l(t,1,1)=num_L_region_nms8_l(t,1,1);
        L_region_nms8_baseline_l(t,1,1)=L_region_nms8_baseline_l(t,1,1);
       
end

L_region_eu15_h=zeros(time,1,1);
 
for t = 1:time
        L_stock_eu15(:, :, t) = permute(Ldyn_eu15_h(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(Ldyn_eu15_h(:, :, t), [2, 1, 3]);
        num_h_region_eu15_h(t,1,1)=L_stock_eu15(1,:,t)+L_stock_eu15(2,:,t)+L_stock_eu15(5,:,t)+L_stock_eu15(6,:,t)...
            +L_stock_eu15(8,:,t)+L_stock_eu15(9,:,t)+L_stock_eu15(10,:,t)+L_stock_eu15(12,:,t)+L_stock_eu15(16,:,t)...
            +L_stock_eu15(17,:,t)+...
             L_stock_eu15(1+17,:,t)+L_stock_eu15(2+17,:,t)+L_stock_eu15(5+17,:,t)+L_stock_eu15(6+17,:,t)...
            +L_stock_eu15(8+17,:,t)+L_stock_eu15(9+17,:,t)+L_stock_eu15(10+17,:,t)+L_stock_eu15(12+17,:,t)+L_stock_eu15(16+17,:,t)...
            +L_stock_eu15(17+17,:,t);
        L_region_eu15_h(t,1,1)=num_h_region_eu15_h(t,1,1);
        L_region_eu15_h(t,1,1)=L_region_eu15_h(t,1,1);
       
end
    
L_region_eu15_baseline_h=zeros(time,1,1);
 
for t = 1:time
        L_stock_eu15(:, :, t) = permute(series_Ldyn_eu15_h(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(series_Ldyn_eu15_h(:, :, t), [2, 1, 3]);
       num_h_region_eu15_h(t,1,1)=L_stock_eu15(1,:,t)+L_stock_eu15(2,:,t)+L_stock_eu15(5,:,t)+L_stock_eu15(6,:,t)...
            +L_stock_eu15(8,:,t)+L_stock_eu15(9,:,t)+L_stock_eu15(10,:,t)+L_stock_eu15(12,:,t)+L_stock_eu15(16,:,t)...
            +L_stock_eu15(17,:,t)+...
            L_stock_eu15(1+17,:,t)+L_stock_eu15(2+17,:,t)+L_stock_eu15(5+17,:,t)+L_stock_eu15(6+17,:,t)...
            +L_stock_eu15(8+17,:,t)+L_stock_eu15(9+17,:,t)+L_stock_eu15(10+17,:,t)+L_stock_eu15(12+17,:,t)+L_stock_eu15(16+17,:,t)...
            +L_stock_eu15(17+17,:,t);
        L_region_eu15_baseline_h(t,1,1)=num_h_region_eu15_h(t,1,1);
        L_region_eu15_baseline_h(t,1,1)=L_region_eu15_baseline_h(t,1,1);
       
end


L_region_eu15_l=zeros(time,1,1);
 
for t = 1:time
        L_stock_eu15(:, :, t) = permute(Ldyn_eu15_l(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(Ldyn_eu15_l(:, :, t), [2, 1, 3]);
        num_L_region_eu15_l(t,1,1)=L_stock_eu15(1,:,t)+L_stock_eu15(2,:,t)+L_stock_eu15(5,:,t)+L_stock_eu15(6,:,t)...
            +L_stock_eu15(8,:,t)+L_stock_eu15(9,:,t)+L_stock_eu15(10,:,t)+L_stock_eu15(12,:,t)+L_stock_eu15(16,:,t)...
            +L_stock_eu15(17,:,t)+...
            L_stock_eu15(1+17,:,t)+L_stock_eu15(2+17,:,t)+L_stock_eu15(5+17,:,t)+L_stock_eu15(6+17,:,t)...
            +L_stock_eu15(8+17,:,t)+L_stock_eu15(9+17,:,t)+L_stock_eu15(10+17,:,t)+L_stock_eu15(12+17,:,t)+L_stock_eu15(16+17,:,t)...
            +L_stock_eu15(17+17,:,t);
        L_region_eu15_l(t,1,1)=num_L_region_eu15_l(t,1,1);
        L_region_eu15_l(t,1,1)=L_region_eu15_l(t,1,1);
       
end
    
L_region_eu15_baseline_l=zeros(time,1,1);
 
for t = 1:time
        L_stock_eu15(:, :, t) = permute(series_Ldyn_eu15_l(:, :, t), [2, 1, 3]);
        L_stock_eu15(:, :, t) = permute(series_Ldyn_eu15_l(:, :, t), [2, 1, 3]);
       num_L_region_eu15_l(t,1,1)=L_stock_eu15(1,:,t)+L_stock_eu15(2,:,t)+L_stock_eu15(5,:,t)+L_stock_eu15(6,:,t)...
            +L_stock_eu15(8,:,t)+L_stock_eu15(9,:,t)+L_stock_eu15(10,:,t)+L_stock_eu15(12,:,t)+L_stock_eu15(16,:,t)...
            +L_stock_eu15(17,:,t)+...
            L_stock_eu15(1+17,:,t)+L_stock_eu15(2+17,:,t)+L_stock_eu15(5+17,:,t)+L_stock_eu15(6+17,:,t)...
            +L_stock_eu15(8+17,:,t)+L_stock_eu15(9+17,:,t)+L_stock_eu15(10+17,:,t)+L_stock_eu15(12+17,:,t)+L_stock_eu15(16+17,:,t)...
            +L_stock_eu15(17+17,:,t);
        L_region_eu15_baseline_l(t,1,1)=num_L_region_eu15_l(t,1,1);
        L_region_eu15_baseline_l(t,1,1)=L_region_eu15_baseline_l(t,1,1);
       
    end

L_region_nms8_baseline=L_region_nms8_baseline_h+L_region_nms8_baseline_l;
L_region_nms8=L_region_nms8_h+L_region_nms8_l;

L_region_eu15_baseline=L_region_eu15_baseline_h+L_region_eu15_baseline_l;
L_region_eu15=L_region_eu15_h+L_region_eu15_l;

L_region_baseline=L_region_nms8_baseline./(L_region_nms8_baseline+L_region_eu15_baseline);
L_region=L_region_nms8./(L_region_nms8+L_region_eu15);
L_tot0=100*L_region_baseline;
L_tot1=100*L_region;
L_tot=100*L_region_baseline - 100*L_region;
L_tot=[L_tot(6);L_tot(14);L_tot(198)]

%High skills
L_region_nms8_baseline=L_region_nms8_baseline_h;
L_region_nms8=L_region_nms8_h;
L_region_eu15_baseline=L_region_eu15_baseline_h;
L_region_eu15=L_region_eu15_h;
L_region_baseline=L_region_nms8_baseline./(L_region_nms8_baseline+L_region_eu15_baseline);
L_region=L_region_nms8./(L_region_nms8+L_region_eu15);
L_h0=100*L_region_baseline;
L_h1=100*L_region;
L_h=100*L_region_baseline - 100*L_region;
L_h=[L_h(6);L_h(14);L_h(198)]

%Low skills
L_region_nms8_baseline=L_region_nms8_baseline_l;
L_region_nms8=L_region_nms8_l;
L_region_eu15_baseline=L_region_eu15_baseline_l;
L_region_eu15=L_region_eu15_l;
L_region_baseline=L_region_nms8_baseline./(L_region_nms8_baseline+L_region_eu15_baseline);
L_region=L_region_nms8./(L_region_nms8+L_region_eu15);
L_l=100*L_region_baseline - 100*L_region;
L_l0=100*L_region_baseline;
L_l1=100*L_region;
L_l=[L_l(6);L_l(14);L_l(198)]

Counterfactual_SS=[Ldyn_eu15_h(:,:,199)' Ldyn_eu15_l(:,:,199)' Ldyn_nms8_h(:,:,199)' Ldyn_nms8_l(:,:,199)' ];
Counterfactual_14=[Ldyn_eu15_h(:,:,14)' Ldyn_eu15_l(:,:,14)' Ldyn_nms8_h(:,:,14)' Ldyn_nms8_l(:,:,14)' ];
Counterfactual_6=[Ldyn_eu15_h(:,:,6)' Ldyn_eu15_l(:,:,6)' Ldyn_nms8_h(:,:,6)' Ldyn_nms8_l(:,:,6)' ];

Baseline_SS=[series_Ldyn_eu15_h(:,:,199)' series_Ldyn_eu15_l(:,:,199)' series_Ldyn_nms8_h(:,:,199)' series_Ldyn_nms8_l(:,:,199)' ];
Baseline_14=[series_Ldyn_eu15_h(:,:,14)' series_Ldyn_eu15_l(:,:,14)' series_Ldyn_nms8_h(:,:,14)' series_Ldyn_nms8_l(:,:,14)' ];
Baseline_6=[series_Ldyn_eu15_h(:,:,6)' series_Ldyn_eu15_l(:,:,6)' series_Ldyn_nms8_h(:,:,6)' series_Ldyn_nms8_l(:,:,6)' ];

Emp_abs_h=100*(series_Ldyn_eu15_h(:,:,199)'+series_Ldyn_nms8_h(:,:,199)'-Ldyn_eu15_h(:,:,199)'-Ldyn_nms8_h(:,:,199)'+series_Ldyn_oth_h(:,:,199)'-Ldyn_oth_h(:,:,199)');
Emp_abs_l=100*(series_Ldyn_eu15_l(:,:,199)'+series_Ldyn_nms8_l(:,:,199)'-Ldyn_eu15_l(:,:,199)'-Ldyn_nms8_l(:,:,199)'+series_Ldyn_oth_l(:,:,199)'-Ldyn_oth_l(:,:,199)');
Emp_abs_agg=Emp_abs_h+Emp_abs_l;


% Table 7
disp(['_________________________________________________________________'])
disp(['                         Section 5.3.1                                   '])
disp(['Welfare effects of trade and migration policy with public goods          '])
disp(['EU-15 High Skill'' = ' num2str(round(Welfare_eu15_h*100)/100) '% '])
disp(['EU-15 Low Skill '' = ' num2str(round(Welfare_eu15_l*100)/100) '% '])
disp(['EU-15 Aggregate '' = ' num2str(round(Welfare_eu15*100)/100) '% '])
disp(['NMS High Skill  '' = ' num2str(round(Welfare_nms8_h*100)/100) '% '])
disp(['NMS Low Skill   '' = ' num2str(round(Welfare_nms8_l*100)/100) '% '])
disp(['NMS Aggregate   '' = ' num2str(round(Welfare_nms8*100)/100) '% '])
disp(['Europe          '' = ' num2str(round(Welfare_europe*100)/100) '% '])
disp(['_________________________________________________________________'])

% Table 5
disp(['_________________________________________________________________'])
disp(['                           Section 5.3.1                                       '])
disp(['Migration effects with public goods         '])
disp(['         All NMS    Low-Skill    High-Skill    '])
disp(['2002'' = ' '   0         0           0      ' ])
disp(['2007'' = ' num2str(round(L_tot(1)*100)/100) '       ' num2str(round(L_l(1)*100)/100) '       ' num2str(round(L_h(1)*100)/100) ])
disp(['2015'' = ' num2str(round(L_tot(2)*100)/100) '       ' num2str(round(L_l(2)*100)/100) '       ' num2str(round(L_h(2)*100)/100)])
disp(['S.S.'' = ' num2str(round(L_tot(3)*100)/100) '       ' num2str(round(L_l(3)*100)/100) '       ' num2str(round(L_h(3)*100)/100)])
disp(['_________________________________________________________________'])

%% Audit save (gated by SAVE_AUDIT env): dump arrays for the split audit Tiers 0-2.
if ~isempty(getenv('SAVE_AUDIT'))
    mu_slots = [2, max(2, round(time/2)), time-2];
    aud = struct();
    aud.tau_P_pref=tau_P_pref; aud.tau_h_base=tau_h_base; aud.tau_l_base=tau_l_base;
    aud.f_pref_const=f_pref_const; aud.f_pref_path=f_pref_path; aud.f_std_path=f_std_path; aud.f_l_path=f_l_path;
    aud.k_temp_pref=k_temp_pref; aud.lambda_r2=lambda_r2; aud.expire_rate=expire_rate;
    aud.time=time; aud.i=i; aud.i2=i2; aud.N=N; aud.J=J; aud.Ymax=Ymax; aud.toldyn=toldyn;
    aud.Ldyn_eu15_h=Ldyn_eu15_h; aud.Ldyn_nms8_h=Ldyn_nms8_h; aud.Ldyn_oth_h=Ldyn_oth_h;
    aud.Ldyn_eu15_h_std=Ldyn_eu15_h_std; aud.Ldyn_nms8_h_std=Ldyn_nms8_h_std; aud.Ldyn_oth_h_std=Ldyn_oth_h_std;
    aud.Ldyn_eu15_h_pref=Ldyn_eu15_h_pref; aud.Ldyn_nms8_h_pref=Ldyn_nms8_h_pref; aud.Ldyn_oth_h_pref=Ldyn_oth_h_pref;
    aud.Ldyn_eu15_h_inelig=Ldyn_eu15_h_inelig; aud.Ldyn_nms8_h_inelig=Ldyn_nms8_h_inelig; aud.Ldyn_oth_h_inelig=Ldyn_oth_h_inelig;
    aud.V_eu15_h=V_eu15_h; aud.V_nms8_h=V_nms8_h; aud.V_oth_h=V_oth_h; aud.V_eu15_l=V_eu15_l;
    aud.V_eu15_h_std=V_eu15_h_std; aud.V_nms8_h_std=V_nms8_h_std; aud.V_oth_h_std=V_oth_h_std;
    aud.V_eu15_h_pref=V_eu15_h_pref; aud.V_nms8_h_pref=V_nms8_h_pref; aud.V_oth_h_pref=V_oth_h_pref;
    aud.mu_slots=mu_slots;
    aud.mu_eu15_h_std=mu_eu15_h_std(:,:,mu_slots); aud.mu_eu15_h_pref=mu_eu15_h_pref(:,:,mu_slots); aud.mu_eu15_h=mu_eu15_h(:,:,mu_slots);
    aud.realwages_h=realwages_h; aud.realwages_h_std=realwages_h_std; aud.realwages_h_pref=realwages_h_pref;
    aud.CEV_h_std_emp=CEV_h_std_emp; aud.CEV_h_pref_emp=CEV_h_pref_emp;
    aud.CEV_h_std_neemp=CEV_h_std_neemp; aud.CEV_h_pref_neemp=CEV_h_pref_neemp;
    aud.CEV_l_emp=CEV_l_emp; aud.CEV_l_neemp=CEV_l_neemp;
    save('-v7','audit_neutral.mat','aud');
    fprintf('[audit] saved audit_neutral.mat (Ymax=%.3e)\n', Ymax);
end

toc
