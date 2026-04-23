# Implementation of Subtrajectory-clustering

### Requirements
- OpenMP

### Install

Install via `pip install .`

Main reference: https://arxiv.org/abs/2308.14865

#### Running test-suite experiments
```
git submodule update --init --recursive
cd discrete-subtrajectory-clustering-test-suite
python main.py run compile
python main.py run socg_table_athens_centers
python main.py plot socg_table_athens_centers
```

##### Subtrajectory Clustering on test-suite data
```
python experiments/athens.py
python experiments/drifter.py
python experiments/unid.py
python experiments/table3.py
python experiments/drifterTable.py
python experiments/unidTable2.py
```
##### Motion Segmentation
```
python experiments/cmu.py
```

##### Lower Bounds via greedy independent set
```
python experiments/drifters.py
```

#### Scaling behaviour experiment
Switch to `drifterMeasurements' branch
