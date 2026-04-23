import json
from pathlib import Path
import re

import pandas as pd

from genericTable import parse_test_suite_tables, create_tables_from_data, join_tables_by_c2, format_tables

# Paths (relative to this script)
script_dir = Path(__file__).parent
latex_path = (script_dir /
    "../discrete-subtrajectory-clustering-test-suite/data/plots/socg_table_unid_centers.tex"
).resolve()
data_path = (script_dir / "../unid_result.json").resolve()

# Read and parse the test suite output
tables_old = parse_test_suite_tables(latex_path.read_text())

# Read the experiment data and construct the rows to be added
#
# Data format: data[c2][l][delta] = {
#   size, filteredSize, cost, mem_MB, time_s
# }
tables_new = create_tables_from_data(json.loads(data_path.read_text()), 362)

# Concatenate the tables by c2, format and output them
for table_str in format_tables(join_tables_by_c2(tables_old + tables_new)):
    print(table_str)
    print()
