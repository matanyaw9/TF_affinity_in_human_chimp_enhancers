import argparse
import importlib
import json
import time
import traceback
import ast
import os
import subprocess
import utils_matanya

JOBS_OUT_DIR = "/home/labs/davidgo/matanyaw/jobs_outputs"
JOBS_ERR_DIR = "/home/labs/davidgo/matanyaw/jobs_errors"

# todo debud and comlete this function to have automated job sending
def run_func_as_job(func, *args, kwargs=None, job_name="run_func_as_job",
                   jobs_outputs_dir=JOBS_OUT_DIR, 
                   jobs_errors_dir=JOBS_ERR_DIR, memory=2000,
                   print_info=True):
    """
    Run a function as a job.
    args:
    func: function, the function to run
    *args: arguments to pass to the function
    job_name: str, the name of the job
    jobs_outputs_dir: str, path to the jobs outputs directory           
    jobs_errors_dir: str, path to the jobs errors directory
    """
    JOB_RUNNER_SCRIPT = "/home/labs/davidgo/matanyaw/backup/job_runner.py"
    os.makedirs(jobs_outputs_dir, exist_ok=True)
    os.makedirs(jobs_errors_dir, exist_ok=True)

    kwargs = kwargs or {}

    job_command = f"python {JOB_RUNNER_SCRIPT} " \
                f"--module {func.__module__} " \
                f"--function {func.__name__} " \
                f"--args '{json.dumps(args)}' " \
                f"--kwargs '{json.dumps(kwargs)}'"

    
    command = [
        "bsub",
        "-q", "short", 
        "-R", f"rusage[mem={memory}]",
        "-J", job_name,
        "-o", os.path.join(jobs_outputs_dir, f"{job_name}.out"),
        "-e", os.path.join(jobs_errors_dir, f"{job_name}.err"),
        job_command
    ]
    
    if print_info:
        print("🔧 Submitting command:\n", ' '.join(command))

    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running the job: {e.stderr}")
        return None
    

def run_any_func(module_name, function_name, args=None, kwargs=None):
    try:
        print(f"🧪 Running: {module_name}.{function_name}")
        mod = importlib.import_module(module_name)
        func = getattr(mod, function_name)

        args = args or []
        kwargs = kwargs or {}

        start_time = time.time()
        result = func(*args, **kwargs)
        duration = time.time() - start_time

        print(f"✅ Job completed in {duration:.2f} seconds.")
        return result

    except Exception as e:
        print(f"❌ Job failed: {e}")
        traceback.print_exc()
        return None

def parse_args():
    parser = argparse.ArgumentParser(description="Dynamic Job Runner")
    parser.add_argument("--module", required=True, help="Module name (e.g. utils_matanya)")
    parser.add_argument("--function", required=True, help="Function name (e.g. create_tf_to_loci_files)")
    parser.add_argument("--args", type=str, help="Positional args as JSON list")
    parser.add_argument("--kwargs", type=str, help="Keyword args as JSON dict")
    return parser.parse_args()

if __name__ == "__main__":
    # set working directory to the script's directory
    import os
    import sys
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    sys.path.append(script_dir)


    args = parse_args()
    arg_list = ast.literal_eval(args.args) if args.args else []
    kwarg_dict = ast.literal_eval(args.kwargs) if args.kwargs else {}

    run_any_func(args.module, args.function, arg_list, kwarg_dict)
