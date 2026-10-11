function run_cops(cc)
%% Rerun of the canonical exercise at the baseline average rates of the static
%% calibration sigma = 1.323, theta = 3.202, rho = 3, phi = 0.72 (2026-10-11):
%% native high-skill 0.2436, native low-skill 0.1668
%% (baseline_sig_1p323_taul_0.216 equilibrium). Previous rerun (recalib_v5):
%% 0.2372 / 0.1598; engine defaults 0.24 / 0.18. All other
%% settings as R4/R5: indexed baseline schedule, closed budget at the target
%% country, k_temp_pref = 10, tau_P homotopy 0.2127 -> 0.19 -> 0.16 -> 0.13.
%% cc = 'pt' (16), 'es' (8) or 'de' (5).
%%
%% Run from cops/runs (outputs and diagnostic folders land there):
%%   octave --no-gui --eval "addpath('..'); run_cops('pt')"
%%
%% Stages: zero-policy sanity-check run (must reproduce the baseline; hard
%% gate), open-budget solve at 0.2127 from the V = 1 seed, closed-budget
%% Anderson homotopy seeded from the open solution, target-country stock dump
%% per step, per-period value-based CEV at tau_P = 0.13.
  switch cc
      case 'pt'; ic = 16;
      case 'es'; ic = 8;
      case 'de'; ic = 5;
      otherwise; error('v6:ctry', 'unknown country tag %s (use pt|es|de)', cc);
  end
  %% read by sweep_engine in every solve and in the stock dumps
  setenv('TAU_H_BASE', '0.2436');
  setenv('TAU_L_BASE', '0.1668');
  pfx = sprintf('v6_%s', cc);
  tracks = {'V_eu15_h','V_nms8_h','V_oth_h','V_eu15_l','V_nms8_l','V_oth_l', ...
            'V_eu15_h_std','V_nms8_h_std','V_oth_h_std', ...
            'V_eu15_h_pref','V_nms8_h_pref','V_oth_h_pref'};
  time = 200; nr = 34;
  V1 = struct();
  for q = 1:numel(tracks); V1.(tracks{q}) = ones(nr, time); end
  np = struct('tol_F', 5e-6, 'maxit', 400, 'm', 15, 'beta_mix', 0.5);
  slots_show = [2 10 50 100 150 199];

  %% ---- gate: zero-policy sanity-check run ----
  dd = sprintf('iter_diag_%s_zero', pfx); mkdir(dd);
  opts = struct('tempeq_tol', 1e-12, 'tau_p_pref', 0.26, 'pref_mix', 0, ...
                'target_ctry', ic, 'diag_dir', dd);
  printf('===== %s zero-policy gate =====\n', pfx); fflush(stdout);
  t0 = tic;
  [Vz, infoz] = anderson_solve(V1, opts, np);
  s2 = 2:time-1;
  devL = max(abs(infoz.dg.lambda_hat_path(s2) - 1));
  devV = max(cellfun(@(tr) max(max(abs(Vz.(tr)(:,s2) - 1))), tracks));
  printf('===== %s zero-policy: |F|=%.3e iters=%d elapsed=%.1fs  max|lambda-1|=%.3e max|V-1|=%.3e =====\n', ...
      pfx, infoz.F_final, infoz.iters, toc(t0), devL, devV); fflush(stdout);
  Vcur = Vz; %#ok<NASGU>
  save('-v7', sprintf('%s_zero_policy.mat', pfx), 'Vcur', 'infoz');
  setenv('PREF_MIX', '0');
  Stz = dump_pt_stocks(sprintf('%s_zero_policy.mat', pfx), sprintf('%s_zero', pfx), false, 0.26, ic);
  unsetenv('PREF_MIX');
  devS = max(abs([Stz.h_ratio(2:end-1), Stz.l_ratio(2:end-1)] - 1));
  printf('%s zero-policy stock check: max|ratio-1|=%.3e\n', pfx, devS); fflush(stdout);
  if devL > 1e-3 || devV > 1e-3 || devS > 1e-3
      error('v6:gate', 'zero-policy gate FAIL at %s: lambda %.3e, V %.3e, stocks %.3e', cc, devL, devV, devS);
  end

  %% ---- open-budget solve at 0.2127 (R1 protocol) ----
  dd = sprintf('iter_diag_%s_open', pfx); mkdir(dd);
  opts = struct('tempeq_tol', 1e-12, 'tau_p_pref', 0.2127, 'budget_off', true, ...
                'target_ctry', ic, 'diag_dir', dd);
  printf('===== %s open solve tau_P=0.2127 =====\n', pfx); fflush(stdout);
  t0 = tic;
  [Vsol, info] = anderson_solve(V1, opts, np);
  printf('===== %s open DONE: |F|=%.3e iters=%d elapsed=%.1fs =====\n', ...
      pfx, info.F_final, info.iters, toc(t0)); fflush(stdout);
  Vcur = Vsol; %#ok<NASGU>
  save('-v7', sprintf('%s_open.mat', pfx), 'Vcur', 'info');

  %% ---- closed-budget homotopy 0.2127 -> 0.19 -> 0.16 -> 0.13 ----
  steps = [0.2127, 0.19, 0.16, 0.13];
  seedfile = sprintf('%s_open.mat', pfx);
  for k = 1:numel(steps)
      tp = steps(k);
      tag = strrep(sprintf('%.4g', tp), '.', 'p');
      s = load(seedfile);
      if isfield(s, 'Vcur'); V0 = s.Vcur; else; V0 = s.Vsol; end
      dd = sprintf('iter_diag_%s_%s', pfx, tag); mkdir(dd);
      opts = struct('tempeq_tol', 1e-12, 'tau_p_pref', tp, 'target_ctry', ic, 'diag_dir', dd);
      printf('===== %s closed step tau_P=%.4f (seed %s) =====\n', pfx, tp, seedfile); fflush(stdout);
      t0 = tic;
      [Vsol, info] = anderson_solve(V0, opts, np);
      printf('===== %s tau_P=%.4f DONE: |F|=%.3e iters=%d elapsed=%.1fs =====\n', ...
          pfx, tp, info.F_final, info.iters, toc(t0));
      %% save FIRST — a print error must not lose a solved step
      Vcur = Vsol; %#ok<NASGU>
      fname = sprintf('%s_taup_%s.mat', pfx, tag);
      save('-v7', fname, 'Vcur', 'info');
      dg = info.dg;
      printf('lambda in [%.4f, %.4f]  max|budget resid|=%.2e\n', ...
          min(dg.lambda_hat_path(2:end-1)), max(dg.lambda_hat_path(2:end-1)), ...
          max(abs(dg.budget_resid)));
      printf('%5s %8s %10s %8s %8s %8s %8s\n', ...
          'slot', 'lambda', 'tau_h_std', 'tau_l', 'f_pref', 'What_h', 'What_l');
      for q = slots_show
          printf('%5d %8.4f %10.4f %8.4f %8.4f %8.4f %8.4f\n', q, ...
              dg.lambda_hat_path(q), dg.tau_h_t_path(q), dg.tau_l_t_path(q), ...
              dg.f_pref_path(q), dg.What_h_path(q), dg.What_l_path(q));
      end
      fflush(stdout);
      dump_pt_stocks(fname, sprintf('%s_%s', pfx, tag), false, tp, ic);
      seedfile = fname;
  end

  %% ---- per-period value-based CEV at the target rate ----
  printf('===== %s CEV paths at tau_P=0.13 =====\n', pfx); fflush(stdout);
  compute_cev_paths(sprintf('%s_taup_0p13.mat', pfx), sprintf('%s_taup13', pfx), ic);
  printf('===== %s DONE =====\n', pfx); fflush(stdout);
end
