function out = compute_cev_paths(solfile, label, target_ctry)
%% compute_cev_paths — per-period CEV from value functions alone.
%% The solved value tracks are hats of TIME-CHANGES (counterfactual/baseline);
%% the reform is unanticipated at period 2, so period-1 values coincide across
%% paths and the same-date value LEVEL ratio telescopes:
%%   R(d,t) = V'(d,t)/V_base(d,t) = prod_{s=2..t} Vhat(d,s),
%% comparing a worker in state d at period t on the reform path with the SAME
%% demographic group (origin bloc, skill, location/employment state) under the
%% baseline (standard tax). CEV in consumption-equivalent %:
%%   CEV(d,t) = 100*(1-beta)*nu*log R(d,t)      (nu = code sigma = 1.6129)
%% At t=2 this is the announcement-period direct measure (= -codeW of
%% cev_direct). No consumption/migration decomposition is used, so the
%% preferential group's baseline-counterpart attribution problem does not
%% arise: the baseline object is the unsplit high-skill value of the same
%% origin bloc, which is exactly the denominator the recursion's hats carry.
%% Terminal caveat: Vhat(T)=1 is imposed, so paths flatten mechanically near T
%% (truncation error ~ beta^(T-t); negligible for t <= ~150).
  if nargin < 2; label = solfile; end
  if nargin < 3; target_ctry = 16; end
  s = load(solfile);
  if isfield(s, 'Vcur'); V = s.Vcur; else; V = s.Vsol; end
  beta = 0.97; nu = 1/0.62;
  blocs = {'eu15', 'nms8', 'oth'};
  kinds = {'h_pref', 'h_std', 'h', 'l'};   % pref, standard-schedule h, unsplit h, low-skill
  states = [target_ctry, target_ctry + 17]; % target-country employed, non-employed (default PT: 16, 33)
  T = size(V.V_eu15_h, 2);
  show = [2 5 10 25 50 100 150 199];
  out = struct('label', label, 'beta', beta, 'nu', nu, 'states', states, 'T', T);
  printf('=== CEV paths from value functions alone — %s ===\n', label);
  printf('CEV(d,t) = 100*(1-beta)*nu*log prod_{s=2..t} Vhat(d,s); vs same demographic group under baseline.\n');
  for ki = 1:numel(kinds)
      for st = 1:2
          d = states(st);
          printf('\n%-7s %s (state %d):  t=  ', kinds{ki}, merge_state(st), d);
          printf('%8d ', show); printf('\n');
          for b = 1:numel(blocs)
              tr = sprintf('V_%s_%s', blocs{b}, kinds{ki});
              vh = V.(tr)(d, :);
              R = [1, cumprod(vh(2:T))];       % R(1)=1; R(t) = prod_{s=2..t}
              cev = 100 * (1 - beta) * nu * log(R);
              out.(sprintf('cev_%s_%s_%d', blocs{b}, kinds{ki}, d)) = cev;
              printf('%-16s ', blocs{b}); printf('%8.3f ', cev(show)); printf('\n');
          end
      end
  end
  save('-v7', sprintf('cev_paths_%s.mat', strrep(label, '.', 'p')), 'out');
  printf('\nsaved cev_paths_%s.mat\n', strrep(label, '.', 'p'));
end

function s = merge_state(st)
  if st == 1; s = 'employed'; else; s = 'non-employed'; end
end
