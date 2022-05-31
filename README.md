# QeX

Library for supproting following experiments with [sequence_parser](https://github.com/qipe-nlab/sequence_parser).

## Supported experimetns
  - [Randomized benchmarking](https://arxiv.org/abs/0707.0963)
  - [Direct fidelity estimation](https://arxiv.org/abs/1104.4695)
  - [Variational quantum gate optimization](https://arxiv.org/abs/1810.12745)
  - [Variational quantum eigensolver](https://arxiv.org/abs/1304.3061)
  - [Subspace-search variational quantum eigensolver](https://arxiv.org/abs/1810.09434)

## Install
Please clone and use the package directly

## Example
```python
from QeX.driver import ExpBase, Circuit
from QeX.util.group import CliffordGroup
from QeX.experiments import RandomizedBenchmarking

# 1. Create circuit
exp = ExpBase(qubit_name_list, cross_name_list)
cir = Circuit(exp)

# 2. Declare the experiment
rb = RandomizedBenchmarking(
  circuit       = cir,
  group         = CliffordGroup(1),
  sequence_list = [(0,10,1000), (2,10,1000), (4,10,1000)],
  seed          = 0,
  interleaved   = None
)

# 3. Define take_data function to meet with your envirionment
def take_data(job_table):
  pass
  
# 4. Execute experiments
rb.execute(take_data)

# 5. Make & Show the report
rb.make_data_table()
rb.make_report()
rb.show_report()
```

## To Do
  - Add setup.py and requirement.txt
  - Add more detail condition on dataset (projector dataset id, clique index)
  - Generalize the driver
