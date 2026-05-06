---
name: abacustest-submit
description: "Submit jobs to Bohrium cloud platform or supercomputers. Use when: user wants to run tasks on remote computing resources by specifying directories, Docker images, and commands."
metadata: { "openclaw": { "emoji": "🚀", "requires": { "pip": ["abacustest"], "env": ["BOHRIUM_USERNAME", "BOHRIUM_PASSWORD", "BOHRIUM_PROJECT_ID"] } } }
---

# abacustest Submit

A flexible job submission tool for running tasks on cloud platforms (Bohrium) or supercomputers. 

**How it works:**
1. Specify directories containing your task files
2. Specify a Docker image with required software
3. Specify the command to execute
4. Submit and let abacustest handle the rest (upload, run, download results)

While commonly used for ABACUS calculations, this tool is **generic** and can run any task that can be containerized.

---

## When to Use

✅ **Use this skill**:
- "Run my calculations on Bohrium cloud"
- "Submit multiple jobs to a remote cluster"
- "Execute tasks in a specific Docker image on the cloud"
- "Track job status and download results automatically"
- "Run high-throughput computations on remote resources"

❌ **Do not use**:
- Prepare input files → Use `abacustest-prepare` or `abacustest-prepare-inputs`
- Extract/analyze results → Use `abacustest-extract-dft-results` or custom scripts
- Local execution → Run commands directly on your machine

---

## Core Concept

abacustest-submit is **not limited to ABACUS**. It's a generic job submission framework:

```
┌─────────────────────────────────────────────────────────┐
│  What You Provide:                                      │
│  1. example: Directories with your files                │
│  2. image: Docker image with required software          │
│  3. command: Command to execute in each directory       │
├─────────────────────────────────────────────────────────┤
│  What abacustest Does:                                  │
│  1. Uploads directories to cloud                        │
│  2. Runs command in specified image                     │
│  3. Downloads results when complete                     │
│  4. Tracks job status                                   │
└─────────────────────────────────────────────────────────┘
```

**Example use cases:**
- DFT calculations (ABACUS, VASP, QE, CP2K, etc.)
- Machine learning training jobs
- Data processing pipelines
- Molecular dynamics simulations
- Any containerized workload

---

## Prerequisites

### 1. Bohrium Account Credentials

**Option A: Environment Variables (Recommended)**
```bash
export BOHRIUM_USERNAME="your_bohrium_username"
export BOHRIUM_PASSWORD="your_bohrium_password"
export BOHRIUM_PROJECT_ID="your_bohrium_project_id"
```

**Option B: In param.json**
```json
{
    "config": {
        "bohrium_username": "your_username",
        "bohrium_password": "your_password",
        "bohrium_project_id": "your_project_id"
    },
    "run_dft": [...]
}
```

### 2. Example Directories

Each **example** is a folder containing files for one job:
```
project/
├── 000000/
│   ├── input.dat
│   ├── config.yaml
│   └── script.py
├── 000001/
│   └── ...
└── 000002/
    └── ...
```

⚠️ **Important**: Use relative paths within directories, not absolute paths.

---

## Configuration File (param.json)

### Basic Structure

```json
{
    "bohrium_group_name": "my-jobs",
    "save_path": "results",
    "max_parallel": 100,
    "run_dft": [
        {
            "ifrun": true,
            "image": "registry.dp.tech/dptech/abacus:LTSv3.10.1",
            "example": ["00[0-2]", "00[3-5]"],
            "group_size": 1,
            "bohrium": {
                "scass_type": "c32_m64_cpu",
                "job_type": "container",
                "platform": "ali"
            },
            "command": "python script.py > log",
            "extra_files": []
        }
    ]
}
```

### Parameter Reference

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `save_path` | string | Directory to save results | `"results"` |
| `max_parallel` | int | Maximum parallel jobs | `100` |
| `bohrium_group_name` | string | Bohrium workflow group name | `"abacustest"` |

### run_dft Array

Each dictionary in `run_dft` defines a job batch:

| Field | Type | Description | Default |
|-------|------|-------------|---------|
| `ifrun` | bool | Skip if `false` | `true` |
| `sub_save_path` | string | Subdirectory for results | `save_path` |
| `image` | string | Docker image name | - |
| `example` | list | Job folder names (supports glob patterns) | - |
| `group_size` | int | Jobs per machine | `1` |
| `bohrium` | dict | Bohrium configuration | - |
| `command` | string | Execution command | - |
| `extra_files` | list | Additional files to upload | `[]` |

### Bohrium Configuration

```json
"bohrium": {
    "scass_type": "c32_m64_cpu",
    "job_type": "container",
    "platform": "ali"
}
```

| Field | Description | Examples |
|-------|-------------|----------|
| `scass_type` | Machine type | `c32_m64_cpu`, `c32_m128_cpu` |
| `job_type` | Job type | `container` |
| `platform` | Cloud platform | `ali`, `paratera` |

⚠️ **Important**: For the Ali platform, their machines are typically dual-threaded, so the actual number of physical cores is halved. For example, a c32_m64 instance has 16 physical cores. We usually use 16 when writing MPI commands.

### Alternative: Dispatcher (Supercomputers)

For non-Bohrium platforms (traditional HPC), use `dispatcher`:

```json
"dispatcher": {
    "machine_dict": {
        "remote_root": "/home/username/work_path",
        "remote_profile": {
            "hostname": "xxx.xx.xxx.xxx",
            "username": "Username",
            "password": "password",
            "port": 22
        }
    },
    "resources_dict": {
        "number_node": 1,
        "cpu_per_node": 8,
        "gpu_per_node": 1,
        "queue_name": "Normal"
    }
}
```

### Upload Custom Packages

If remote image lacks required Python packages:

```json
"upload_packages": ["my_package1", "/home/user/my_package2"]
```

These packages will be uploaded and installed on the remote platform before job execution.

---

## Commands

### Submit Jobs

```bash
abacustest submit -p param.json
```

**Run in background** (recommended for long jobs):
```bash
abacustest submit -p param.json &
```

**With prepare step** (prepare inputs then submit):
```bash
# param.json contains both "prepare" and "run_dft" sections
abacustest submit -p param.json
```

### Track Job Status

After submission, you'll see:
```
Workflow has been submitted (ID: abacustest-kbrb2, UID: f14c5d95-655c-47b4-a709-c9a8138a40cf)
Workflow link: https://workflows.deepmodeling.com/workflows/argo/abacustest-kbrb2
```

Visit the URL to track job status in your browser.

### Download Results

**Using Job ID:**
```bash
abacustest download -p param.json abacustest-kbrb2
```

**Using UID** (if ID expired):
```bash
abacustest download -p param.json f14c5d95-655c-47b4-a709-c9a8138a40cf
```

---

## Examples

### Example 1: Simple Job Submission

**param.json:**
```json
{
    "save_path": "results",
    "run_dft": [
        {
            "image": "registry.dp.tech/dptech/abacus:LTSv3.10.1",
            "example": ["job1", "job2", "job3"],
            "bohrium": {
                "scass_type": "c32_m64_cpu",
                "job_type": "container",
                "platform": "ali"
            },
            "command": "OMP_BUN_THREADS=1 mpirun -np 16 abacus | tee out.log"
        }
    ]
}
```

**Submit:**
```bash
abacustest submit -p param.json &
```

### Example 2: Using Glob Patterns

**param.json:**
```json
{
    "save_path": "high_throughput",
    "max_parallel": 50,
    "run_dft": [
        {
            "image": "registry.dp.tech/dptech/abacus-stable:LTSv3.10",
            "example": ["00[0-9]", "01[0-9]"],
            "group_size": 2,
            "bohrium": {
                "scass_type": "c32_m64_cpu",
                "job_type": "container",
                "platform": "ali"
            },
            "command": "OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log",
            "extra_files": ["analyze.py", "config.json"]
        }
    ]
}
```

### Example 3: Prepare + Submit in One Step

**param.json:**
```json
{
    "prepare": {
        "strus": ["structure1.cif", "structure2.cif"],
        "stru_format": "cif",
        "input_template": "INPUT",
        "pp_path": "/path/to/pp"
    },
    "save_path": "results",
    "run_dft": [
        {
            "image": "registry.dp.tech/dptech/abacus-stable:LTSv3.10",
            "bohrium": {
                "scass_type": "c8_m16_cpu",
                "job_type": "container",
                "platform": "ali"
            },
            "command": "mpirun -np 4 abacus > log"
        }
    ]
}
```

**Note**: When `example` is not defined in `run_dft`, the generated inputs from `prepare` will be used automatically.

### Example 4: Supercomputer Submission (via Dispatcher)

**param.json:**
```json
{
    "save_path": "hpc_results",
    "run_dft": [
        {
            "example": ["sim1", "sim2"],
            "dispatcher": {
                "machine_dict": {
                    "remote_root": "/home/user/work",
                    "remote_profile": {
                        "hostname": "hpc.example.com",
                        "username": "user",
                        "password": "pass",
                        "port": 22
                    }
                },
                "resources_dict": {
                    "number_node": 1,
                    "cpu_per_node": 16,
                    "queue_name": "batch"
                }
            },
            "command": "./run_simulation.sh"
        }
    ]
}
```

### Example 5: Machine Learning Training

**param.json:**
```json
{
    "save_path": "ml-results",
    "run_dft": [
        {
            "image": "pytorch/pytorch:2.0-cuda11.7",
            "example": ["exp1", "exp2", "exp3"],
            "bohrium": {
                "scass_type": "gpu_v100_1x",
                "job_type": "container",
                "platform": "ali"
            },
            "command": "python train.py --config config.yaml > train.log",
            "extra_files": ["train.py", "config.yaml", "requirements.txt"]
        }
    ]
}
```

---

## Output Files

After job completion and download:

```
results/
├── 000000/
│   ├── output/
│   │   └── ...
│   ├── log
│   └── time.json
├── 000001/
│   └── ...
├── 000002/
│   └── ...
└── info.json              # Job metadata
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Run submit in background (`&`) | Jobs can take hours to days |
| Use glob patterns for many jobs | `["00[0-9]"]` is cleaner than listing 10 folders |
| Set `max_parallel` appropriately | Avoid overwhelming computing resources |
| Save job ID/UID | Needed for downloading results later |
| Use `extra_files` for scripts | Upload analysis scripts with jobs |
| Check Bohrium credits before submit | Ensure sufficient computing credits |
| Test with one job first | Verify command and image work before scaling |
| Use `group_size > 1` for small jobs | Improves resource utilization |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Authentication failed | Verify `BOHRIUM_*` environment variables or config credentials |
| Job stuck in pending | Check `scass_type` availability; try different machine type |
| Results not downloading | Use UID instead of ID if job is old |
| Absolute path errors | Use relative paths within example directories |
| Image not found | Verify image name; check `ABBREVIATION` in abacustest |
| Job failed immediately | Check `log` or `out.log` for error messages |
| Missing dependencies | Use `upload_packages` or include in custom image |
| Permission denied | Ensure files in example directories are readable |

---

## Related Skills

- **Input preparation**: [`abacustest-prepare`](../abacustest-prepare/SKILL.md), [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- **Result extraction**: [`abacustest-extract-dft-results`](../abacustest-extract-dft-results/SKILL.md)
- **Specialized models**: [`abacustest-models`](../abacustest-models/SKILL.md)

---

## External Resources

- **Bohrium Platform**: https://bohrium.dp.tech
- **dpdispatcher Docs**: https://docs.deepmodeling.com/projects/dpdispatcher/
- **abacustest GitHub**: https://github.com/pxlxingliang/abacus-test
