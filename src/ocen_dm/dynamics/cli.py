"""Command-line entry points; long evolutions require an explicit run command."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser("nbody", help="prepare and run controlled dynamical experiments")
    actions = parser.add_subparsers(dest="nbody_action", required=True)
    prepare = actions.add_parser("prepare", help="validate and sample joint equilibrium; save orbit and inputs")
    prepare.add_argument("config", type=Path)
    prepare.add_argument("--label", required=True)
    prepare.add_argument("--output-root", type=Path, default=None)
    prepare.set_defaults(func=command_prepare)
    run = actions.add_parser("run", help="evolve a prepared model and compute diagnostics")
    run.add_argument("directory", type=Path)
    run.add_argument("--engine", choices=("live", "frozen"), default="live")
    run.add_argument("--label", default=None, help="fresh evolution label within the prepared directory")
    run.add_argument("--nemo-env", type=Path, default=None)
    run.set_defaults(func=command_run)
    analyse = actions.add_parser("analyse", help="recompute diagnostics for a completed evolution")
    analyse.add_argument("directory", type=Path)
    analyse.add_argument("--nemo-env", type=Path, default=None)
    analyse.set_defaults(func=command_analyse)
    plot = actions.add_parser("plot", help="compare enclosed mass, shell density and bound mass")
    plot.add_argument("directories", type=Path, nargs="+")
    plot.add_argument("--output", type=Path, required=True)
    plot.set_defaults(func=command_plot)


def command_prepare(args):
    from ..paths import results_dir
    from .config import Experiment, safe_label
    from .experiment import prepare
    root = args.output_root if args.output_root is not None else results_dir()/"dynamics"
    print(prepare(Experiment.load(args.config), root/safe_label(args.label)))
    return 0


def command_run(args):
    from .experiment import run
    print(run(args.directory, args.engine, args.label, args.nemo_env))
    return 0


def command_analyse(args):
    from .experiment import analyse
    print(analyse(args.directory, args.nemo_env))
    return 0


def command_plot(args):
    from .plotting import plot_evolutions
    print(plot_evolutions(args.directories, args.output))
    return 0
