#!/usr/bin/env python3
import sys
import re
import argparse
import subprocess

from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument(
    "--help", help="Display help message.", action="store_true")
parser.add_argument(
    "positional", action="append",
    nargs="?", metavar="POS",
    help="additional arguments not in slurm parser group to pass to sbatch")

# A subset of SLURM-specific arguments
slurm_parser = parser.add_argument_group("slurm-specific arguments")
slurm_parser.add_argument(
    "-a", "--array", help="job array index values")
slurm_parser.add_argument(
    "--begin", help="defer job until HH:MM MM/DD/YY")
slurm_parser.add_argument(
    "-c", "--cpus-per-task", help="number of cpus required per task")
slurm_parser.add_argument(
    "-d", "--dependency",
    help="defer job until condition on jobid is satisfied")
slurm_parser.add_argument(
    "-D", "--workdir", help="set working directory for batch script")
slurm_parser.add_argument(
    "-e", "--error", help="file for batch script's standard error")
slurm_parser.add_argument(
    "-J", "--job-name", help="name of job")
slurm_parser.add_argument(
    "--mail-type", help="notify on state change: BEGIN, END, FAIL or ALL")
slurm_parser.add_argument(
    "--mail-user", help="who to send email notification for job state changes")
slurm_parser.add_argument(
    "-n", "--ntasks", help="number of tasks to run")
slurm_parser.add_argument(
    "-N", "--nodes", help="number of nodes on which to run (N = min[-max])")
slurm_parser.add_argument(
    "-o", "--output", help="file for batch script's standard output")
slurm_parser.add_argument(
    "-p", "--partition", help="partition requested")
slurm_parser.add_argument(
    "-Q", "--quiet", help="quiet mode (suppress informational messages)")
slurm_parser.add_argument(
    "-t", "--time", help="time limit")
slurm_parser.add_argument(
    "--wrap", help="wrap command string in a sh script and submit")
slurm_parser.add_argument(
    "-C", "--constraint", help="specify a list of constraints")
slurm_parser.add_argument(
    "--mem", help="minimum amount of real memory")
slurm_parser.add_argument(
    "-x", "--exclude", help="Explicitly exclude certain nodes from the resources granted to the job")

args = parser.parse_args()

if args.help:
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

extras = ""
if args.positional:
    for m in args.positional:
        if m is not None:
            extras = extras + " " + m

arg_dict = dict(args.__dict__)


# Process resources
if "resources" in job_properties:
    resources = job_properties["resources"]
    if arg_dict["time"] is None:
        if "runtime" in resources:
            arg_dict["time"] = resources["runtime"] * 60
        elif "walltime" in resources:
            arg_dict["time"] = resources["walltime"] * 60
        elif "time" in resources:
            arg_dict["time"] = resources["time"] * 60
        else:
            arg_dict['time'] = "2:00:00"
    if "mem" in resources and arg_dict["mem"] is None:
        arg_dict["mem"] = resources["mem"] * 1000
    if "cores" in resources:
       arg_dict["cpus-per-task"] = resources["cores"]

# Threads
# defaults to 1 and will overwrite cores from resources 
# if not careful here
if "threads" in job_properties:
    if job_properties["threads"] >1:
        arg_dict["cpus-per-task"] = job_properties["threads"]
    
#job name
arg_dict["job-name"] = job_properties["rule"]


opt_keys = ["array", "begin", "cpus-per-task",
            "depedency", "workdir", "error", "job-name", "mail_type",
            "mail_user", "ntasks", "nodes", "output", "partition",
            "quiet", "time", "wrap", "constraint", "mem", "exclude"]


if 'partition' in job_properties['cluster']:
    #override the default behavior and allow specification of partition when required.
    arg_dict['partition'] = job_properties['cluster']['partition'] #
else:
    if arg_dict["partition"] is None:
        arg_dict["partition"] = "euan,owners,normal"
        arg_dict["exclude"] = "sh02-13n15"

opts = ""
for k, v in arg_dict.items():
    if k not in opt_keys:
        continue
    if v is not None:
        opts += " --{} \"{}\" ".format(k, v)

if arg_dict["wrap"] is not None:
    cmd = "sbatch {opts}".format(opts=opts)
else:
    cmd = "sbatch {opts} {extras}".format(opts=opts, extras=extras)

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise
