Execute `abacustest collectdata -h` to see the help message.
Execute `abacustest collectdata -p abacus-pw-abacus.json -j *` to collect data from all jobs in the current directory.
Execute `abacustest collectdata -p abacus-pw-qe.json -j * -t 1` to collect the QE data from all jobs in the current directory.
Execute `abacustest collectdata --outparam -t 0` to check the metric name of abacus job, and `-t 1` to check the metric name of QE job, and `-t 2` to check the metric name of VASP job.
