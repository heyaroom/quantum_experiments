# QeX

Library for supproting following experiments with [sequence_parser](https://github.com/qipe-nlab/sequence_parser).

## Supported experimetns
  - [Randomized benchmarking](https://arxiv.org/abs/0707.0963)
  - [Quantum tomograpy & Direct fidelity estimation](https://arxiv.org/abs/1104.4695)
  - [Variational quantum gate optimization](https://arxiv.org/abs/1810.12745)
  - [Variational quantum eigensolver](https://arxiv.org/abs/1304.3061)
  - [Subspace-search variational quantum eigensolver](https://arxiv.org/abs/1810.09434)

## Example
```python
from sequence_parser import Circuit
from QeX.util.group import CliffordGroup
from QeX.experiments import RandomizedBenchmarking

# 1. Import Circuit from "sequence_parser"
cir = Circuit(backend)

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
