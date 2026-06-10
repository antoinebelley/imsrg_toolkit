"""Submit per-sample pipelines as chained SLURM job arrays.

Instead of submitting 4 jobs per sample (imsrg -> diag -> density -> expvals,
~4*N jobs per isotope), each stage becomes a single job array with one task
per sample. Consecutive stages are linked with --dependency=aftercorr:<jobid>
so task i of a stage starts only once task i of the previous stage has
completed successfully (see https://slurm.schedmd.com/job_array.html).
This keeps the scheduler load at one job per stage regardless of N.
"""

import os
from pathlib import Path
from subprocess import run, PIPE

# SubMIT slurm.conf: MaxArraySize = 1001, i.e. valid task ids are 0..1000
MAX_ARRAY_SIZE = 1001


class JobArrayStage():
  def __init__(self, name, header, workdir=None):
    self.name = name
    self.header = header      # '#SBATCH ...' lines (no shebang), resources for ONE task
    self.workdir = workdir    # directory each task cds into before running its scripts
    self.sample_ids = []
    self.scripts = []         # list of lists: scripts run sequentially within one task

  def add_task(self, sample_id, scripts):
    if isinstance(scripts, (str, Path)):
      scripts = [scripts]
    self.sample_ids.append(sample_id)
    self.scripts.append([str(s) for s in scripts])


class JobArrayChain():
  def __init__(self, name, script_directory, submit_cmd='sbatch'):
    self.name = name
    self.script_directory = str(script_directory)
    self.submit_cmd = submit_cmd
    self.stages = []
    Path(self.script_directory).mkdir(parents=True, exist_ok=True)

  def new_stage(self, name, header, workdir=None):
    stage = JobArrayStage(name, header, workdir=workdir)
    self.stages.append(stage)
    return stage

  def write_stage_files(self, stage):
    fn_manifest = f"{self.script_directory}/{self.name}_{stage.name}.manifest"
    with open(fn_manifest, "w") as f:
      for sample_id, scripts in zip(stage.sample_ids, stage.scripts):
        f.write("\t".join([str(sample_id)] + scripts) + "\n")

    s = "#!/bin/bash\n"
    s += stage.header.rstrip() + "\n\n"
    s += f'MANIFEST="{fn_manifest}"\n'
    s += 'LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")"\n'
    s += 'if [ -z "$LINE" ]; then\n'
    s += '    echo "No manifest entry for task $SLURM_ARRAY_TASK_ID in $MANIFEST" >&2\n'
    s += '    exit 1\n'
    s += 'fi\n'
    s += 'IFS=$\'\\t\' read -r -a FIELDS <<< "$LINE"\n'
    s += 'SAMPLEID="${FIELDS[0]}"\n'
    s += f'echo "stage={stage.name} SampleID=$SAMPLEID array=$SLURM_ARRAY_JOB_ID task=$SLURM_ARRAY_TASK_ID"\n'
    if stage.workdir is not None:
      s += f'cd "{stage.workdir}" || exit 1\n'
    s += 'for SCRIPT in "${FIELDS[@]:1}"; do\n'
    s += '    "$SCRIPT"\n'
    s += '    RC=$?\n'
    s += '    if [ $RC -ne 0 ]; then\n'
    s += '        echo "ERROR: $SCRIPT exited with code $RC for SampleID=$SAMPLEID" >&2\n'
    s += '        exit $RC\n'
    s += '    fi\n'
    s += 'done\n'

    fn_script = f"{self.script_directory}/{self.name}_{stage.name}_array.sh"
    with open(fn_script, "w") as f:
      f.write(s)
    os.chmod(fn_script, 0o755)
    return fn_script

  def submit(self, max_concurrent=None, dry_run=False, verbose=True):
    """Submit one job array per stage, chained with aftercorr dependencies.

    max_concurrent limits how many tasks of each array run at once
    (sbatch --array=0-N%max_concurrent). Returns {stage_name: jobid}.
    """
    ntasks = {len(stage.sample_ids) for stage in self.stages}
    if len(ntasks) != 1:
      raise ValueError(f"All stages must have the same number of tasks, got {ntasks}")
    n = ntasks.pop()
    if n == 0:
      raise ValueError("No tasks added to the job array chain")
    if n > MAX_ARRAY_SIZE:
      raise ValueError(f"{n} tasks exceeds MaxArraySize={MAX_ARRAY_SIZE}; split the submission")

    array_spec = f"0-{n - 1}"
    if max_concurrent is not None:
      array_spec += f"%{max_concurrent}"

    jobids = {}
    previous_jobid = None
    for stage in self.stages:
      fn_script = self.write_stage_files(stage)
      cmd = [self.submit_cmd, '--parsable', f'--array={array_spec}']
      if previous_jobid is not None:
        cmd += [f'--dependency=aftercorr:{previous_jobid}', '--kill-on-invalid-dep=yes']
      cmd.append(fn_script)
      if dry_run:
        print(f"[dry run] {' '.join(cmd)}")
        previous_jobid = f"<{stage.name}_jobid>"
        continue
      jobid = run(cmd, stdout=PIPE, text=True, check=True).stdout.rstrip()
      jobids[stage.name] = jobid
      if verbose:
        print(f"Submitted {stage.name} array ({n} tasks) with jobid {jobid}")
      previous_jobid = jobid
    return jobids
