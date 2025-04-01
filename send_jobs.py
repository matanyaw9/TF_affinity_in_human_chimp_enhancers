import subprocess
import os
import pandas as pd 

DATA_DIR = "/home/labs/davidgo/matanyaw/data/"
MPRA_FILE = "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/comparative_analysis_combined/humanMPRA_with_seq_final2.csv"
PBM_FILE = "/home/labs/davidgo/matanyaw/data/pbm_8mer_aggregated_data.csv"
JOBS_OUT_DIR = "/home/labs/davidgo/matanyaw/jobs_outputs"
JOBS_ERR_DIR = "/home/labs/davidgo/matanyaw/jobs_errors"


def run_single_job(mpra_csv_file, pbm_csv_file, output_dir, window_size, start_index, end_index,
                   job_name="process_TF_to_loci_files", 
                   jobs_outputs_dir=JOBS_OUT_DIR, 
                   jobs_errors_dir=JOBS_ERR_DIR):
    """
    Run a single job to create the TF to loci files.
    args:
    mpra_csv_file: str, path to the MPRA CSV file
    pbm_csv_file: str, path to the PBM CSV file
    output_dir: str, path to the output directory
    window_size: int, the window size to use (for finding sequences in the PBM file)
    start_index: int, the start index to use (for the MPRA file, which can be split into multiple jobs)
    end_index: int, the end index to use (for the MPRA file, which can be split into multiple jobs)
    job_name: str, the name of the job
    jobs_outputs_dir: str, path to the jobs outputs directory
    jobs_errors_dir: str, path to the jobs errors directory
    """
    
    if not os.path.exists(jobs_outputs_dir):
        os.makedirs(jobs_outputs_dir)
    if not os.path.exists(jobs_errors_dir):
        os.makedirs(jobs_errors_dir)

    TF_2_LOCI_SCRIPT_FILE = "/home/labs/davidgo/matanyaw/backup/create_TF_to_loci_files.py"
    
    command = [
        "bsub",  # Executable
        "-q", "short", 
        "-R", "rusage[mem=2000]",   # Resource requirements in Mb
        "-J", job_name,
        "-o", jobs_outputs_dir,
        "-e", jobs_errors_dir,
        "python", TF_2_LOCI_SCRIPT_FILE,  # The Python script to run
        "--mpra_file", mpra_csv_file,
        "--pbm_file", pbm_csv_file,
        "--output_dir", output_dir,
        "--window_size", str(window_size),
        "--start_index", str(start_index),
        "--end_index", str(end_index),
    ]

    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running the job: {e.stderr}")
        return None
    
if __name__ == "__main__":
    
    job_name = "job_test_2"
    num_jobs = 1

    OUTPUT_DIR_NAME = "job_test_2"
    WINDOW_SIZE = 8
    start_index = 0
    end_index = 5

    OUTPUT_DIR = os.path.join(DATA_DIR, OUTPUT_DIR_NAME)
    
    
    run_single_job(mpra_csv_file=MPRA_FILE,
            pbm_csv_file=PBM_FILE,
            output_dir=OUTPUT_DIR,
            window_size=WINDOW_SIZE,
            start_index=start_index,
            end_index=end_index, 
            job_name=job_name)


#     # Now let's send the jobs
#     mpra_df = pd.read_csv(MPRA_FILE)
#     mpra_df = mpra_df[mpra_df['differential_activity'] == True]
#     mpra_df.reset_index(drop=True, inplace=True)
#     num_jobs = 150
#     chunk_size = mpra_df.shape[0] // num_jobs
#     for i in range(num_jobs):
#         start_index = i * chunk_size
#         end_index = min((i + 1) * chunk_size, mpra_df.shape[0]) 
#         job_name = f"job_{i}"
#         run_single_job(mpra_file=MPRA_FILE,
#                 pbm_file=PBM_FILE,
#                 output_dir=OUTPUT_DIR,
#                 window_size=WINDOW_SIZE,
#                 start_index=start_index,
#                 end_index=end_index, 
#                 job_name=job_name)
#         print(f"Sent job {i+1}/{num_jobs} for lines: {start_index} to {end_index}")

# print(f"All {num_jobs} Jobs are Sent!")