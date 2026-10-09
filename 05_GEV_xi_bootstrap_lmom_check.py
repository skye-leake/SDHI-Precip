#!/usr/bin/env python3
"""
gev_xi_bootstrap.py -- compute step, with selectable estimator.

Same as before, plus --estimator {lmom,mle}.

  lmom  Hosking & Wallis L-moment estimator. Closed form, ~20 us per fit.
        The estimator Papalexiou & Koutsoyiannis used, and the one
        Hosking et al. (1985) show outperforms MLE for n < 100. Whole
        run takes seconds.

  mle   scipy genextreme.fit. ~50 ms per fit, ~2500x slower.

Recommended: run lmom as primary. Failing to constrain xi with the
estimator most favourable to free estimation is the stronger result,
and it keeps "short-record estimates are biased toward zero" true --
that holds for L-moments but NOT for MLE, which is biased upward at
n = 15.

Sign convention: Hosking's k = scipy's c = -xi, so xi = -k throughout,
matching the previous script. Validated by recovering known shapes
from large samples (see validate_estimator() below -- run it once).

     nohup python3 05_GEV_xi_bootstrap_lmom_check.py --base /home/scratch/idf_curves --tag v1 --estimator lmom --outdir gev_xi_lmom --workers 8 > gev_xi_lmom.log 2>&1 &

"""

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from itertools import product

import numpy as np
import pandas as pd
import xarray as xr
from scipy.special import gamma as _gamma
from scipy.stats import genextreme

CITIES = ['albany', 'amarillo', 'grand_junction', 'minneapolis',
          'nashville', 'phoenix', 'seattle', 'tallahassee']
DURATIONS = ['15min', '24hr']
SCENARIO = "HIST"
YEARS = "1990-2005"
CELL_STRIDE = 5
N_BOOT = 1000
N_MIN_YEARS = 10
MASTER_SEED = 20260726
BOOT_CONVERGE_MIN = 0.90


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S} pid{os.getpid():>6}] {msg}",
          flush=True)


# ---------------------------------------------------------------- estimators
def fit_xi_lmom(x):
    """Hosking & Wallis L-moment estimator for the GEV shape parameter.
    Returns xi in the climate convention (xi > 0 = heavy tail)."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n < 4:
        return None
    j = np.arange(1, n + 1)
    b0 = x.mean()
    b1 = np.sum((j - 1) / (n - 1) * x) / n
    b2 = np.sum((j - 1) * (j - 2) / ((n - 1) * (n - 2)) * x) / n
    l2 = 2 * b1 - b0
    l3 = 6 * b2 - 6 * b1 + b0
    if not np.isfinite(l2) or l2 <= 0:
        return None
    t3 = l3 / l2
    c = 2.0 / (3.0 + t3) - np.log(2.0) / np.log(3.0)
    k = 7.8590 * c + 2.9554 * c * c          # Hosking's k = scipy's c
    return float(-k) if np.isfinite(k) else None


def fit_xi_mle(x):
    """scipy maximum-likelihood GEV fit. No plausibility filtering."""
    try:
        c, loc, scale = genextreme.fit(x)
        xi = -c
        if not np.isfinite([xi, loc, scale]).all() or scale <= 0:
            return None
        return float(xi)
    except Exception:
        return None


ESTIMATORS = {'lmom': fit_xi_lmom, 'mle': fit_xi_mle}


def validate_estimator(fit, rng=None):
    """Recover known shapes from large samples. Run once before trusting
    a new estimator -- catches sign-convention errors, which are the
    likely failure mode here."""
    rng = rng or np.random.default_rng(0)
    print("estimator validation (n = 20000):")
    ok = True
    for xi_true in (-0.2, 0.0, 0.114, 0.3):
        x = genextreme.rvs(-xi_true, loc=50, scale=14, size=20000,
                           random_state=rng)
        est = fit(x)
        err = abs(est - xi_true)
        flag = "OK " if err < 0.03 else "FAIL"
        ok &= err < 0.03
        print(f"   {flag}  true {xi_true:+.3f}  estimated {est:+.4f}")
    if not ok:
        raise SystemExit("estimator failed validation -- check sign convention")
    print("   validation passed\n")


# ---------------------------------------------------------------- bootstrap
def bootstrap_xi(x, n_boot, rng, fit):
    out = np.full(n_boot, np.nan)
    for b in range(n_boot):
        v = fit(rng.choice(x, size=x.size, replace=True))
        if v is not None:
            out[b] = v
    return out


def task_paths(outdir, city, dur):
    stem = os.path.join(outdir, f"{city}__{dur}")
    return stem + ".csv", stem + "_boot.npz"


def run_task(payload):
    city, dur, seed_seq, cfg = payload
    fit = ESTIMATORS[cfg['estimator']]
    csv_path, npz_path = task_paths(cfg['outdir'], city, dur)
    t0 = time.time()
    rng = np.random.default_rng(seed_seq)

    nc = os.path.join(cfg['base'], SCENARIO, city,
                      f"combined_annual_maxima_{SCENARIO}_{YEARS}_{city}"
                      f"_{cfg['tag']}.nc")
    log(f"{city}/{dur}: START [{cfg['estimator']}]  {nc}")
    if not os.path.exists(nc):
        log(f"{city}/{dur}: ERROR file not found")
        return {'city': city, 'duration': dur, 'status': 'missing_file'}

    ds = xr.open_dataset(nc)
    try:
        lat_dim = ('south_north' if 'south_north' in ds.sizes
                   else [d for d in ds.sizes if d != 'year'][0])
        lon_dim = ('west_east' if 'west_east' in ds.sizes
                   else [d for d in ds.sizes if d not in ('year', lat_dim)][0])
        var = f"max_intensity_{dur}"
        if var not in ds.data_vars:
            log(f"{city}/{dur}: ERROR variable {var} absent")
            return {'city': city, 'duration': dur, 'status': 'missing_var'}

        pts = list(product(range(0, ds.sizes[lat_dim], CELL_STRIDE),
                           range(0, ds.sizes[lon_dim], CELL_STRIDE)))
        log(f"{city}/{dur}: grid {ds.sizes[lat_dim]}x{ds.sizes[lon_dim]}, "
            f"{len(pts)} cells at stride {CELL_STRIDE}")

        rows, draws = [], np.full((len(pts), cfg['n_boot']), np.nan)

        for k, (i, j) in enumerate(pts):
            am = ds[var].isel({lat_dim: i, lon_dim: j}).values
            x = am[~np.isnan(am)]
            row = {'city': city, 'duration': dur, 'i': i, 'j': j,
                   'n_years': int(x.size), 'estimator': cfg['estimator']}

            if x.size < N_MIN_YEARS:
                row['outcome'] = 'too_few_years'
            else:
                xi_hat = fit(x)
                if xi_hat is None:
                    row['outcome'] = 'fit_failed'
                else:
                    b = bootstrap_xi(x, cfg['n_boot'], rng, fit)
                    draws[k, :] = b
                    n_ok = int(np.isfinite(b).sum())
                    row.update({'xi': xi_hat, 'n_boot_ok': n_ok})
                    if n_ok < BOOT_CONVERGE_MIN * cfg['n_boot']:
                        row['outcome'] = 'bootstrap_failed'
                    else:
                        lo, hi = np.nanpercentile(b, [2.5, 97.5])
                        row.update({'outcome': 'ok', 'xi_lo': float(lo),
                                    'xi_hi': float(hi)})
            rows.append(row)
            if (k + 1) % cfg['log_every'] == 0 or (k + 1) == len(pts):
                log(f"{city}/{dur}: cell {k+1}/{len(pts)}")
    finally:
        ds.close()

    df = pd.DataFrame(rows)
    tmp = csv_path + ".tmp"
    df.to_csv(tmp, index=False)
    os.replace(tmp, csv_path)
    tmp = npz_path + ".tmp.npz"
    np.savez_compressed(tmp, draws=draws, ij=np.array(pts, dtype=int))
    os.replace(tmp, npz_path)

    counts = df.outcome.value_counts().to_dict()
    log(f"{city}/{dur}: DONE in {time.time()-t0:.1f}s  {counts}")
    return {'city': city, 'duration': dur, 'status': 'ok',
            'seconds': time.time() - t0, **counts}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--estimator', choices=list(ESTIMATORS), default='lmom')
    ap.add_argument('--outdir', default=None,
                    help='default: gev_xi_<estimator>')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--log-every', type=int, default=5)
    ap.add_argument('--cities', nargs='*', default=CITIES)
    ap.add_argument('--durations', nargs='*', default=DURATIONS)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--skip-validation', action='store_true')
    args = ap.parse_args()

    outdir = args.outdir or f"gev_xi_{args.estimator}"
    os.makedirs(outdir, exist_ok=True)

    if not args.skip_validation:
        validate_estimator(ESTIMATORS[args.estimator])

    cfg = {'base': args.base, 'tag': args.tag, 'outdir': outdir,
           'n_boot': args.n_boot, 'log_every': args.log_every,
           'estimator': args.estimator}

    tasks = [(c, d) for c in sorted(args.cities) for d in args.durations]
    children = np.random.SeedSequence(MASTER_SEED).spawn(len(tasks))

    todo = []
    for (city, dur), seed in zip(tasks, children):
        csv_path, _ = task_paths(outdir, city, dur)
        if os.path.exists(csv_path) and not args.force:
            log(f"{city}/{dur}: SKIP (already done)")
            continue
        todo.append((city, dur, seed, cfg))

    log(f"estimator={args.estimator} | seed {MASTER_SEED} | {len(todo)}/"
        f"{len(tasks)} tasks | {args.workers} workers | {args.n_boot} resamples")
    if not todo:
        log("nothing to do")
        return

    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_task, p): (p[0], p[1]) for p in todo}
        for n, fut in enumerate(as_completed(futures), 1):
            city, dur = futures[fut]
            try:
                results.append(fut.result())
            except Exception as e:
                log(f"{city}/{dur}: FAILED {type(e).__name__}: {e}")
                results.append({'city': city, 'duration': dur,
                                'status': 'exception'})
            log(f"--- {n}/{len(todo)} complete ({time.time()-t0:.0f}s) ---")

    log(f"finished in {time.time()-t0:.0f}s")
    log(pd.DataFrame(results).to_string(index=False))

    frames, missing = [], []
    for city, dur in tasks:
        csv_path, _ = task_paths(outdir, city, dur)
        (frames.append(pd.read_csv(csv_path)) if os.path.exists(csv_path)
         else missing.append(f"{city}/{dur}"))
    if missing:
        log(f"WARNING incomplete, missing: {', '.join(missing)}")
    if frames:
        combined = os.path.join(outdir,
                                f"gev_xi_cells_{args.tag}_{args.estimator}.csv")
        pd.concat(frames, ignore_index=True).to_csv(combined, index=False)
        log(f"combined -> {combined}")


if __name__ == '__main__':
    sys.exit(main())