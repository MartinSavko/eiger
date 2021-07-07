#!/usr/bin/env ccp4-python
"""Module to run SIMBAD"""

__author__ = "Adam Simpkin, and Felix Simkovic"
__contributing_authors__ = "Jens Thomas, and Ronan Keegan"
__credits__ = "Daniel Rigden, William Shepard, Martin Savko, Charles Ballard, Villi Uski, and Andrey Lebedev"
__date__ = "18 Feb 2019"
__email__ = "hlasimpk@liv.ac.uk"
__version__ = "0.1"

import argparse
import os
import sys

from pyjob.misc import StopWatch
from pyjob.factory import TASK_PLATFORMS

import simbad.command_line
import simbad.exit

logger = None


def _argparse_core_options(p):
    """Add core options to an already existing parser"""
    sg = p.add_argument_group('Basic options')
    sg.add_argument('-ccp4_jobid', type=int,
                    help='Set the CCP4 job id - only needed when running from the CCP4 GUI')
    sg.add_argument('-ccp4i2_xml', help=argparse.SUPPRESS)
    sg.add_argument('-chunk_size', default=0, type=int,
                    help='Max jobs to submit at any given time')
    sg.add_argument('-debug_lvl', type=str, default='info', choices=['info', 'debug', 'warning', 'error', 'critical'],
                    help='The console verbosity level')
    sg.add_argument('-name', type=str, default="simbad",
                    help='The identifier for each job [simbad]')
    sg.add_argument('-output_pdb', type=str,
                    help='Path to the output PDB for the best result')
    sg.add_argument('-output_mtz', type=str,
                    help='Path to the output MTZ for the best result')
    sg.add_argument('-run_dir', type=str, default=".",
                    help='Directory where the SIMBAD work directory will be created')
    sg.add_argument('-results_to_display', type=int, default=10,
                    help='The number of results to display in the GUI')

    # TODO: Update this location to the correct path for the temp dir
    sg.add_argument('-tmp_dir', type=str, default='/nfs/data',
                    help='Directory in which to put temporary files from SIMBAD')

    sg.add_argument('-work_dir', type=str,
                    help='Path to the directory where SIMBAD will run (will be created if it doesn\'t exist)')
    sg.add_argument('-webserver_uri',
                    help='URI of the webserver directory - also indicates we are running as a webserver')
    sg.add_argument('-rvapi_document', help=argparse.SUPPRESS)
    sg.add_argument('-tab_prefix', type=str, default="", help=argparse.SUPPRESS)
    sg.add_argument('--cleanup', default=False,
                    action="store_true", help="Delete all data not reported by the GUI")
    sg.add_argument('--display_gui', default=False,
                    action="store_true", help="Show the SIMBAD GUI")
    sg.add_argument('--process_all', default=False,
                    action="store_true", help="Trial all search models")
    sg.add_argument('--skip_mr', default=False,
                    action="store_true", help="Skip Molecular replacement step")
    sg.add_argument('--version', action='version', version='SIMBAD v{0}'.format(simbad.version.__version__),
                    help='Print the SIMBAD version')

def _argparse_job_submission_options(p):
    """Add the options for submission to a cluster queuing system"""
    sg = p.add_argument_group('Cluster queue submission options')
    sg.add_argument('-nproc', type=int, default=20,
                    help="Number of processors. For local, serial runs the jobs will be split across nproc "
                         "processors. For cluster submission, this should be the number of processors on a node.")
    sg.add_argument('-submit_nproc', type=int, default=1,
                    help="For cluster submission, the number of processors to use on head node when creating "
                         "submission scripts")
    sg.add_argument('-submit_qtype', type=str, default='local', choices=TASK_PLATFORMS.keys(),
                    help='The job submission queue type')
    sg.add_argument('-submit_queue', type=str, default=None,
                    help='The queue to submit to on the cluster.')


def simbad_argparse():
    """Create the argparse options"""
    p = argparse.ArgumentParser(
        description="SIMBAD: Sequence Independent Molecular replacement Based on Available Database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    _argparse_core_options(p)
    _argparse_job_submission_options(p)
    simbad.command_line._argparse_contaminant_options(p)
    simbad.command_line._argparse_lattice_options(p)
    simbad.command_line._argparse_rot_options(p)
    simbad.command_line._argparse_mr_options(p)
    simbad.command_line._argparse_mtz_options(p)
    p.add_argument('mtz', help="The path to the input mtz file")
    return p


def main():
    """Main SIMBAD routine"""
    args = simbad_argparse().parse_args()

    args.work_dir = simbad.command_line.get_work_dir(
        args.run_dir, work_dir=args.work_dir, ccp4_jobid=args.ccp4_jobid, ccp4i2_xml=args.ccp4i2_xml
    )

    log_file = os.path.join(args.work_dir, 'simbad.log')
    debug_log_file = os.path.join(args.work_dir, 'debug.log')
    log_class = simbad.command_line.LogController()
    log_class.add_console(level=args.debug_lvl)
    log_class.add_logfile(log_file, level="info", format="%(message)s")
    log_class.add_logfile(debug_log_file, level="notset",
                          format="%(asctime)s\t%(name)s [%(lineno)d]\t%(levelname)s\t%(message)s")
    global logger
    logger = log_class.get_logger()

    if not os.path.isfile(args.amore_exe):
        raise OSError("amore executable not found")

    simbad.command_line.print_header()
    logger.info("Running in directory: %s\n", args.work_dir)

    stopwatch = StopWatch()
    stopwatch.start()

    end_of_cycle, solution_found, all_results = False, False, {}
    while not (solution_found or end_of_cycle):
        # =====================================================================================
        # Perform the lattice search
        solution_found = simbad.command_line._simbad_lattice_search(args)
        logger.info("Lattice search completed in %d days, %d hours, %d minutes, and %d seconds",
                    *stopwatch.lap.time_pretty)

        if solution_found and not args.process_all:
            logger.info(
                "Lucky you! SIMBAD worked its charm and found a lattice match for you.")
            continue
        elif solution_found and args.process_all:
            logger.info(
                "SIMBAD thinks it has found a solution however process_all is set, continuing to contaminant search")
        else:
            logger.info("No results found - lattice search was unsuccessful")

        if args.output_pdb and args.output_mtz:
            csv = os.path.join(args.work_dir, 'latt/lattice_mr.csv')
            all_results['latt'] = simbad.util.result_by_score_from_csv(csv, 'final_r_free', ascending=True)

        # =====================================================================================
        # Perform the contaminant search
        solution_found = simbad.command_line._simbad_contaminant_search(args)
        logger.info("Contaminant search completed in %d days, %d hours, %d minutes, and %d seconds",
                    *stopwatch.lap.time_pretty)

        if solution_found:
            logger.info(
                "Check you out, crystallizing contaminants! But don't worry, SIMBAD figured it out and found a solution.")
            continue
        else:
            logger.info(
                "No results found - contaminant search was unsuccessful")

        if args.output_pdb and args.output_mtz:
            csv = os.path.join(args.work_dir, 'cont/cont_mr.csv')
            all_results['cont'] = simbad.util.result_by_score_from_csv(csv, 'final_r_free', ascending=True)

        # =====================================================================================
        # Make sure we only run the loop once for now
        end_of_cycle = True

    if len(all_results) >= 1:
        sorted_results = sorted(all_results.iteritems(), key=lambda (k, v): (v[1], k))
        result = sorted_results[0][1]
        simbad.util.output_files(args.work_dir, result, args.output_pdb, args.output_mtz)

    stopwatch.stop()
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds",
                *stopwatch.time_pretty)

    log_class.close()


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.NOTSET)
    try:
        main()
    except Exception:
        simbad.exit.exit_error(*sys.exc_info())
