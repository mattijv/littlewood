from pathlib import Path
from argparse import Action, ArgumentParser

def get_batch_template(job_name, project_id, time_limit, cpus, memory_limit, N, subdivisions, buckets):
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}-bucket-%BUCKET%
#SBATCH --account=Project_{project_id}
#SBATCH --partition=small
#SBATCH --time={time_limit}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --nodes=1
#SBATCH --mem={memory_limit}

srun /projappl/project_2007026/littlewood/lw -N{N} -j{cpus} -s{subdivisions} -B{buckets} -b%BUCKET%
"""

project_id = "2007026"

parser = ArgumentParser(description = 'CLI utility to create CSC batch files.')

parser.add_argument("-N", required = True, type = int)
parser.add_argument("-j", "--job-name")
parser.add_argument("-t", "--time-limit", default = "00:10:00")
parser.add_argument("-c", "--cpus", type = int, default = 10)
parser.add_argument("-m", "--memory-limit", default = "1G")
parser.add_argument("-s", "--subdivisions", type = int, default = 0)
parser.add_argument("-b", "--buckets", type = int, default = 1)

args = parser.parse_args()

job_name = args.job_name
if job_name is None:
    job_name = f"run-n-{args.N}"

# Ensure batch directory exists
Path(f"batch/{job_name}").mkdir(parents = True, exist_ok = True)

run_script = """#!/bin/bash
"""

batch_template = get_batch_template(job_name, project_id, args.time_limit, args.cpus, args.memory_limit, args.N, args.subdivisions, args.buckets)
for b in range(1, args.buckets + 1):
    file_name = f"batch/{job_name}/{job_name}-bucket-{b}.sh"
    run_script += f"sbatch {job_name}-bucket-{b}.sh\n"
    content = batch_template.replace("%BUCKET%", str(b))
    with open(file_name, "w") as batch_file:
        batch_file.write(content)

# Write the run script
run_script_file_name = f"batch/{job_name}/{job_name}-run.sh"
with open(run_script_file_name, "w") as run_script_file:
    run_script_file.write(run_script)
