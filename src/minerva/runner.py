"""
Run MINERVA on a single target or on a catalog of targets.

run(config_path) reads the configuration and runs run_single_target and/or
run_catalog. For each target, the Monte Carlo simulations are run by
simulation.run_simulation; the per-target statistics by summary.print_execution_summary.
"""
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

from . import utils as minerva_utils
from .context import RunContext
from .logs import ToggleableLogger
from .simulation import run_simulation, setup_seasons_logic
from .summary import print_execution_summary

# Header of the per-target results file (one line per simulation)
TEXTHEAD = 'Target Name, Vmag, Planetary Mass [Mj], Planetary Radius [Rj], Planetary Period [d], Eccentricity, Periastron [deg], Age [Gyr], Tag, Prot [d], logRHK, Kact [m/s], Kpl [m/s], std(RV) [m/s], sig_ph [m/s], std(noise_rv) [m/s], Nobs, Niter, Kfit [m/s], Kfit_err [m/s], SNRfinal, dBIC, jitter_act [m/s], simult_flag, Kstar_fit [m/s], Pstar_fit [d], status_flag, sig_best, bootstrap_flag\n'


def run_single_target(var):
    """Monte Carlo simulations for the target defined in single_target_params."""
    temp_params = minerva_utils.load_single_target_params(var.config)
    target_params = temp_params
    outseries = 'plot_series_' + target_params['targetname']
    log_filename = var.OutDir + f"Log_{temp_params['targetname']}_simulation.txt"

    # Log file: from here on, print() writes to the log only
    sys.stdout.write(f"Single target: {temp_params['targetname']}\nLog: {log_filename}\n")

    logger = ToggleableLogger(log_filename, to_file=True, to_screen=False)
    sys.stdout = logger

    print(f"--- LOG {temp_params['targetname']} ---\nDate: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Progress bar (on screen, not in the log)
    iterator = tqdm(range(var.numero_iterazioni),
                    desc=f"Simulations {temp_params['targetname']}",
                    unit="iter", file=sys.__stdout__)

    # SNR thresholds, in increasing order
    if isinstance(var.sig_best, (list, tuple, np.ndarray)):
        thresholds = sorted(list(var.sig_best))
    else:
        thresholds = [float(var.sig_best)]

    with PdfPages(var.FigDir + outseries + '.pdf') as pdf:

        # One random generator for the whole target (reproducible results)
        rng = np.random.default_rng(var.analysis["random_seed"])
        # NaN = simulation not run (run_simulation returned None): excluded from the statistics
        var.nobs_history = {th: np.full(var.numero_iterazioni, np.nan) for th in thresholds}
        var.status_flags = {th: np.full(var.numero_iterazioni, np.nan) for th in thresholds}

        debug_pdf_filename = var.FigDir + f"DEBUG_{target_params['targetname']}.pdf"
        debug_pdf = PdfPages(debug_pdf_filename) if var.options["export_debug_pdf"] else None
        log_Rhk = None   # defined even if every simulation fails
        for niter in iterator:
            target_params = minerva_utils.load_single_target_params(var.config)

            print(f"\n{'-'*50}\nSIMULATION {niter+1}/{var.numero_iterazioni}\n")
            setup_seasons_logic(var, target_params)
            minerva_utils.print_params(target_params)

            # Per-target results file
            var.outputfilename = var.OutDir + 'results_' + target_params['targetname'] + '.txt'
            if niter == 0:
                with open(var.outputfilename, 'w') as f:   # new file at every run
                    f.write(TEXTHEAD)

            # One Monte Carlo simulation
            result = run_simulation(var, target_params, niter, pdf, debug_pdf, rng)
            if result is None: continue   # simulation not run (invalid star)
            log_Rhk, nobs_dict, status_dict = result

            # Number of observations and status for each SNR threshold
            for th in thresholds:
                var.nobs_history[th][niter] = nobs_dict[th]
                var.status_flags[th][niter] = status_dict[th]

            if niter == var.numero_iterazioni - 1:
                print('Target completed')

        # Statistics over all the simulations of the target
        print_execution_summary(
                target_params['targetname'], log_Rhk, var.nobs_history, var.status_flags, var.outputfilename, thresholds,
                is_bootstrapped=var.options["bootstrap"]
            )
        sys.__stdout__.write("\nTarget completed!\n")

        if debug_pdf is not None:
            debug_pdf.close()


def _simulate_target(var, i, thresholds, pdf, show_progress=True):
    """All the Monte Carlo simulations of the i-th catalog target.

    Writes the per-target results file and returns (target name, statistics per threshold,
    failures per threshold, number of simulations). Used by both the sequential and the
    parallel catalog runs, so that the two give identical results.
    """
    # One random generator per target, with a different seed for each target
    rng = np.random.default_rng(var.analysis["random_seed"] + i)
    # NaN = simulation not run: excluded from the statistics
    var.nobs_history = {th: np.full(var.numero_iterazioni, np.nan) for th in thresholds}

    # Target parameters, loaded once per target (each simulation works on a copy)
    target_params = minerva_utils.load_target_params_from_catalog(var.catalog, i)
    var.status_flags = {th: np.full(var.numero_iterazioni, np.nan) for th in thresholds}
    print(f"\n{'#'*60}\n# TARGET {i+1}/{len(var.catalog)}: {target_params['targetname']}\n{'#'*60}\n")
    # Per-target results file
    var.outputfilename = var.OutDir + 'results_' + target_params['targetname'] + '.txt'
    with open(var.outputfilename, 'w') as f:   # new file at every run
        f.write(TEXTHEAD)

    debug_pdf_filename = var.FigDir + f"DEBUG_{target_params['targetname']}.pdf"
    debug_pdf = PdfPages(debug_pdf_filename) if var.options["export_debug_pdf"] else None
    # Inner progress bar: disappears when done (leave=False), shown on screen despite the logger
    inner_iterator = tqdm(range(var.numero_iterazioni),
                          desc=f"Simulations {target_params['targetname']}",
                          unit="it",
                          leave=False,
                          file=sys.__stdout__,
                          disable=not show_progress)
    log_Rhk = None   # defined even if every simulation fails
    for niter in inner_iterator:
        print(f"   > Simulation {niter+1}/{var.numero_iterazioni}...")

        setup_seasons_logic(var, target_params)
        minerva_utils.print_params(target_params)

        # One Monte Carlo simulation
        result = run_simulation(var, target_params, niter, pdf, debug_pdf, rng)
        if result is None: continue

        log_Rhk, nobs_dict, status_dict = result

        # Number of observations and status for each SNR threshold
        for th in thresholds:
            var.nobs_history[th][niter] = nobs_dict[th]
            var.status_flags[th][niter] = status_dict[th]

        if niter == var.numero_iterazioni - 1:
            print(f"Target completed: {target_params['targetname']}")

    # Statistics of the target: log, screen and per-target results file
    results_stats, fails_tot = print_execution_summary(
        target_params['targetname'], log_Rhk, var.nobs_history, var.status_flags, var.outputfilename, thresholds,
        is_bootstrapped=var.options["bootstrap"]
    )

    if debug_pdf is not None:
        debug_pdf.close()

    return target_params['targetname'], results_stats, fails_tot, niter + 1


def _catalog_line(targetname, results_stats, fails_tot, thresholds, n_sim):
    """One line of the catalog summary file, with every SNR threshold."""
    cat_str = f"{targetname}, "
    for th in thresholds:
        m_all, p16_all, p84_all, m_succ, p16_succ, p84_succ = results_stats[th]
        cat_str += f"{m_all:.1f}, {p16_all:.1f}, {p84_all:.1f}, {m_succ:.1f}, {p16_succ:.1f}, {p84_succ:.1f}, "
    # Failures per threshold, written as "[n1, n2]"
    fails_str = ", ".join(str(fails_tot[th]) for th in thresholds)
    cat_str += f"{n_sim}, \"[{fails_str}]\" \n"
    return cat_str


# --- parallel catalog runs: one RunContext per worker process ---------------------
_WORKER_VAR = None


def _init_worker(config_path, overrides):
    """Initialise a worker process: configuration, reference data and catalog."""
    global _WORKER_VAR
    import matplotlib
    matplotlib.use('Agg')   # workers never show figures on screen
    var = RunContext(config_path, overrides)
    var.namecatalog = var.config['namecatalog']
    var.catalog = minerva_utils.import_file_catalog(var.resolve(var.namecatalog), InDir='')
    _WORKER_VAR = var


def _worker_target(i, thresholds):
    """Run one catalog target in a worker, with its own log and RV-series PDF."""
    var = _WORKER_VAR
    name = minerva_utils.load_target_params_from_catalog(var.catalog, i)['targetname']
    logger = ToggleableLogger(var.OutDir + f"Log_{name}.txt", to_file=True, to_screen=False)
    sys.stdout = logger
    try:
        print(f"--- LOG {name} ---\nDate: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        with PdfPages(var.FigDir + f"plot_series_{name}.pdf") as pdf:
            out = _simulate_target(var, i, thresholds, pdf, show_progress=False)
    finally:
        logger.restore_stdout()
    return i, out


def _n_workers(option):
    """Number of worker processes from run_options.n_workers (an integer or 'auto')."""
    if isinstance(option, str) and option.lower() == 'auto':
        return max(1, os.cpu_count() or 1)
    return max(1, int(option))


def run_catalog(var):
    """Monte Carlo simulations for every target of the catalog.

    With run_options.n_workers > 1 the targets run in parallel processes. Every target has
    its own random seed (random_seed + i), so the results are identical to a sequential run;
    each target then writes its own log (Log_<target>.txt) and RV-series PDF.
    """
    var.namecatalog = var.config['namecatalog']
    outseries = 'plot_series_' + var.config['namecatalog']
    var.catalog = minerva_utils.import_file_catalog(var.resolve(var.namecatalog), InDir='')
    n_workers = min(_n_workers(var.options["n_workers"]), len(var.catalog))

    # Log file of the catalog: from here on, print() writes to the log only
    catalog_name_clean = var.namecatalog.split('.')[0]
    log_filename = var.OutDir + f"Log_Catalog_{catalog_name_clean}.txt"

    sys.stdout.write(f"Catalog: {len(var.catalog)} stars.\nLog: {log_filename}\n")
    logger = ToggleableLogger(log_filename, to_file=True, to_screen=False)
    sys.stdout = logger

    print(f"--- LOG CATALOG {catalog_name_clean} ---\nDate: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Summary file of the catalog (one line per target)
    var.outcatalog = var.OutDir + 'results_' + catalog_name_clean + '.txt'

    # SNR thresholds, in increasing order
    if isinstance(var.sig_best, (list, tuple, np.ndarray)):
        thresholds = sorted(list(var.sig_best))
    else:
        thresholds = [float(var.sig_best)]

    with open(var.outcatalog, 'w') as f:  # 'w': a new summary file at every run
        header_str = "Target Name, "

        for th in thresholds:
            header_str += f"Med_ALL_SNR{th}, 16th_ALL_SNR{th}, 84th_ALL_SNR{th}, Med_SUCC_SNR{th}, 16th_SUCC_SNR{th}, 84th_SUCC_SNR{th}, "
        header_str += "Niter, Fail_Tot\n"
        f.write(header_str)

    if n_workers <= 1:
        # Sequential run: one log and one RV-series PDF for the whole catalog
        with PdfPages(var.FigDir + outseries + '.pdf') as pdf:
            catalog_iterator = tqdm(range(len(var.catalog)),
                                    desc="Catalog progress",
                                    unit="star", file=sys.__stdout__)
            for i in catalog_iterator:
                out = _simulate_target(var, i, thresholds, pdf)
                with open(var.outcatalog, 'a') as f:
                    f.write(_catalog_line(*out[:3], thresholds, out[3]))
    else:
        # Parallel run: targets distributed over worker processes, summary written in catalog order
        print(f"Parallel run with {n_workers} worker processes: one log per target (Log_<target>.txt)")
        results = {}
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx,
                                 initializer=_init_worker, initargs=(var.config_path, var.overrides)) as pool:
            futures = [pool.submit(_worker_target, i, thresholds) for i in range(len(var.catalog))]
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc=f"Catalog progress ({n_workers} workers)", unit="star", file=sys.__stdout__):
                i, out = fut.result()
                results[i] = out
                print(f"Target {i+1}/{len(var.catalog)} completed: {out[0]}")
        with open(var.outcatalog, 'a') as f:
            for i in range(len(var.catalog)):
                f.write(_catalog_line(*results[i][:3], thresholds, results[i][3]))


def run(config_path, overrides=None):
    """Run MINERVA with the given configuration file (overrides: see RunContext)."""
    var = RunContext(config_path, overrides)
    try:
        if var.config['single_target']:
            run_single_target(var)
        if var.config['read_catalog']:
            run_catalog(var)
    finally:
        if isinstance(sys.stdout, ToggleableLogger):
            sys.stdout.restore_stdout()
    return var
