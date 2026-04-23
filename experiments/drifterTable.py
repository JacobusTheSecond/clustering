import json
from pathlib import Path
import re

import pandas as pd

from genericTable import parse_test_suite_tables, create_tables_from_data, format_table

# Paths (relative to this script)
script_dir = Path(__file__).parent
latex_path = (script_dir /
    "../discrete-subtrajectory-clustering-test-suite/data/plots/socg_table_drifter_centers.tex"
).resolve()
data_path = (script_dir / "../drifter_result.json").resolve()

# Read and parse the test suite output
tables_old = parse_test_suite_tables(latex_path.read_text())

# Read the experiment data and construct the rows to be added
#
# Data format: data[c2][l][delta] = {
#   size, filteredSize, cost, mem_MB, time_s
# }
tables_new = create_tables_from_data(json.loads(data_path.read_text()), 2011)

# Concatenate the tables and format them
table_str = format_table(pd.concat(tables_old + tables_new))

# Print the resulting table
print(table_str)
print()
