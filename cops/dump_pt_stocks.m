function St = dump_pt_stocks(solfile, label, budget_off, taup, target_ctry)
%% dump_pt_stocks — run ONE sweep at a solved value path and extract target-country
%% worker stocks by skill and sub-type along the transition. Standalone diagnostic:
%% mirrors sweep_phi's environment setup, runs sweep_engine in this function's
%% workspace, and reads the Ldyn_* arrays. Does NOT edit sweep_engine. Default
%% target PT: states 16 = employed, 33 = non-employed; "total" sums both.
%%   solfile    : .mat carrying the value struct (Vcur / Vsol / Vout)
%%   label      : tag for the diag dir
%%   budget_off : true = open budget; false = closed budget (engine uses closed-form lambda)
%%   taup       : preferential rate
%%   target_ctry: target-country index in the 18-country ordering (default 16 = PT)
  if nargin < 5; target_ctry = 16; end
  s = load(solfile);
  if isfield(s,'Vcur'); V=s.Vcur; elseif isfield(s,'Vsol'); V=s.Vsol; else; V=s.Vout; end

  setenv('OUTER_MAX','1');
  if target_ctry == 16; unsetenv('TARGET_CTRY'); else; setenv('TARGET_CTRY', num2str(target_ctry)); end
  if budget_off; setenv('BUDGET_OFF','1'); else; unsetenv('BUDGET_OFF'); end
  unsetenv('PREF_SCHED');
  setenv('TAU_P_PREF', num2str(taup,'%.10g'));
  dd = sprintf('iter_diag_stocks_%s', label); mkdir(dd); setenv('DIAG_DIR', dd);
  setenv('TEMPEQ_TOL','1e-12'); setenv('TEMPEQ_MAXIT','20000');
  %% snapshot bookkeeping BEFORE the run: sweep_engine reuses names like `taup`,
  %% so these locals would be clobbered by run('sweep_engine.m') otherwise.
  bk_label = label; bk_bo = budget_off; bk_taup = taup;
  INJECT_V = V;     %#ok<NASGU>  read by sweep_engine insertion 1
  SWEEP_MODE = 1;   %#ok<NASGU>  engine raises sweep:done after the forward/backward sweep
  try
      run('sweep_engine.m');
      error('dump:noStop','engine ran to completion without sweep:done');
  catch err
      if ~strcmp(err.identifier,'sweep:done'); rethrow(err); end
  end

  ie = target_ctry; ine = target_ctry + 17;
  %% counterfactual stocks at PT (sum over origin blocs eu15/nms8/oth), emp+nonemp
  St.h_std   = ptvec(Ldyn_eu15_h_std  + Ldyn_nms8_h_std  + Ldyn_oth_h_std,  ie, ine);
  St.h_pref  = ptvec(Ldyn_eu15_h_pref + Ldyn_nms8_h_pref + Ldyn_oth_h_pref, ie, ine);
  St.h_unspl = ptvec(Ldyn_eu15_h      + Ldyn_nms8_h      + Ldyn_oth_h,      ie, ine);
  St.l_tot   = ptvec(Ldyn_eu15_l      + Ldyn_nms8_l      + Ldyn_oth_l,      ie, ine);
  St.h_tot   = St.h_std + St.h_pref;
  %% baseline stocks at PT (no std/pref split in the baseline)
  St.h_base  = ptvec(series_Ldyn_eu15_h + series_Ldyn_nms8_h + series_Ldyn_oth_h, ie, ine);
  St.l_base  = ptvec(series_Ldyn_eu15_l + series_Ldyn_nms8_l + series_Ldyn_oth_l, ie, ine);
  %% shares and cf/baseline ratios
  St.pref_share_of_h = St.h_pref ./ St.h_tot;          % preferential share of PT high-skill
  St.h_share_of_pt   = St.h_tot ./ (St.h_tot + St.l_tot);
  St.h_ratio = St.h_tot ./ St.h_base;                  % cf / baseline, high-skill
  St.l_ratio = St.l_tot ./ St.l_base;                  % cf / baseline, low-skill
  St.label = bk_label; St.budget_off = bk_bo; St.taup = bk_taup;
  if exist('tau_P_pref','var'); St.taup_used = tau_P_pref(1); else; St.taup_used = NaN; end
  St.time = numel(St.l_tot);

  save('-v7', sprintf('pt_stocks_%s.mat', bk_label), 'St');
  snaps = [1 2 10 50 100 150 St.time-1];
  printf('\n=== PT stocks: %s (budget_off=%d, tau_P arg=%.4f, engine tau_P_pref=%.4f) — emp+nonemp, blocs summed ===\n', bk_label, bk_bo, bk_taup, St.taup_used);
  printf(' t     L_low      L_high     h_std      h_pref   pref/h%%  h-share%%  Lh/base  Ll/base\n');
  for t = snaps
      printf('%4d  %9.3e %9.3e %9.3e %9.3e  %6.2f   %6.2f   %6.4f   %6.4f\n', ...
          t, St.l_tot(t), St.h_tot(t), St.h_std(t), St.h_pref(t), ...
          100*St.pref_share_of_h(t), 100*St.h_share_of_pt(t), St.h_ratio(t), St.l_ratio(t));
  end
  fflush(stdout);
end

function y = ptvec(X, ie, ine)
%% extract PT (employed + non-employed) total as a 1 x time vector from a stock array
  Xs = squeeze(X);
  if size(Xs,1) ~= 34 && size(Xs,2) == 34; Xs = Xs.'; end
  y = Xs(ie,:) + Xs(ine,:);
end
