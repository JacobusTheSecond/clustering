from collections import defaultdict
import math
import pandas as pd
import re

def parse_latex_table_row(row: str) -> list[str]:
    return [
        value.strip()
        for value in row.strip().rstrip("\\").split("&")
    ]

def parse_test_suite_tables(latex: str) -> list[pd.DataFrame]:
    tabular_blocks = re.findall(
        r'\\begin\{tabular\}.*?\\end\{tabular\}',
        latex,
        flags=re.DOTALL
    )

    result = []

    for block in tabular_blocks:
        lines = block.splitlines()

        top_idx = next(i for i, l in enumerate(lines) if r"\toprule" in l)
        mid_idx = next(i for i, l in enumerate(lines) if r"\midrule" in l)
        bottom_idx = next(i for i, l in enumerate(lines) if r"\bottomrule" in l)

        if top_idx + 2 != mid_idx or mid_idx + 1 > bottom_idx:
            raise ValueError("Invalid tabular block: incorrect toprule, midrule and bottomrule")

        column_names = parse_latex_table_row(lines[top_idx + 1])

        rows = [
            parse_latex_table_row(row)
            for row in lines[mid_idx + 1:bottom_idx]
        ]

        for row in rows:
            if len(row) != len(column_names):
                raise ValueError("Invalid tabular block: inconsistent number of columns")

        # Transposing rows (see https://stackoverflow.com/a/6473724)
        columns = list(map(list, zip(*rows)))

        name_values = columns[column_names.index("Name")]
        seconds_values = columns[column_names.index("Seconds")]
        gbytes_values = columns[column_names.index("GBytes")]
        c2_vec_values = columns[column_names.index("$c2$")]
        clusters_values = columns[column_names.index("Clusters")]
        score_values = columns[column_names.index("kCenters")]

        def parse_optional_value_func(base_parse_func):
            return (lambda value: None if value == "---" else base_parse_func(value))

        def parse_c2_vec(c2_vec: str):
            return tuple([
                float(value.strip())
                for value in c2_vec.strip("()").split(",")
            ])

        curr_result = pd.DataFrame({
            "name": name_values,
            "seconds": list(map(parse_optional_value_func(float), seconds_values)),
            "gbytes": list(map(parse_optional_value_func(float), gbytes_values)),
            "c2_vec": list(map(parse_optional_value_func(parse_c2_vec), c2_vec_values)),
            "clusters": list(map(parse_optional_value_func(float), clusters_values)),
            "score": list(map(parse_optional_value_func(float), score_values)),
        })

        result.append(curr_result)

    return result

def create_tables_from_data(data: dict, c3_value: int) -> list[pd.DataFrame]:
    def items_sorted_numerically(data: dict):
        return sorted(map(
            lambda pair: (float(pair[0]), pair[1]),
            data.items(),
        ))

    result = []
    
    for c2_value, c2_data in items_sorted_numerically(data):
        curr_result = pd.DataFrame(
            {
                "name": pd.Series(dtype="str"),
                "seconds": pd.Series(dtype="float64"),
                "gbytes": pd.Series(dtype="float64"),
                "c2_vec": pd.Series(dtype="object"),
                "clusters": pd.Series(dtype="float64"),
                "score": pd.Series(dtype="float64"),
            }
        )

        c2_vec = (1, c2_value, c3_value)

        for l, c2_l_data in items_sorted_numerically(c2_data):
            for delta, info in items_sorted_numerically(c2_l_data):
                name = "our"
                seconds = float(info["time_s"])
                gbytes = float(info["mem_MB"] / 1024)
                clusters = float(info["filteredSize"])
                score = float(info["cost"])

                curr_result.loc[len(curr_result)] = [
                    name,
                    seconds,
                    gbytes,
                    c2_vec + (int(l), float(delta)),
                    clusters,
                    score,
                ]

        result.append(curr_result)

    return result

def join_tables_by_c2(tables: list[pd.DataFrame]) -> list[pd.DataFrame]:
    c2_to_tables = defaultdict(list)

    for table in tables:
        table_c2_values = set(map(
            lambda c2_vec: c2_vec[1],
            filter(lambda value: not value is None, table["c2_vec"]),
        ))

        # if len(table_c2_values) != 1:
        #     raise ValueError("Cannot determine the value of c2 from a table")

        # c2_value = next(iter(table_c2_values))

        # For some test suite result tables, c2 is set incorrectly.
        # This is a dirty workaround

        if len(table_c2_values) == 0:
            raise ValueError("Cannot determine the value of c2 from a table")

        c2_value = min(table_c2_values)

        c2_to_tables[c2_value].append(table)

    return [
        pd.concat(tables, ignore_index=True)
        for c2, tables in sorted(c2_to_tables.items())
    ]

def format_float(value: float):
    return str(int(value) if value.is_integer() else value)

def format_name(name: str, c2_vec: tuple):
    if name == "rightstep":
        return "PSC"

    if name == "basic-length":
        length = str(int(c2_vec[3])) if not c2_vec is None else "?"
        return f"SC-$l$-{length}"

    if name == "basic-size":
        size = str(int(c2_vec[4])) if not c2_vec is None else "?"
        return f"SC-$m$-{size}"

    if name == "java":
        return "MapConstruct"

    if name == "our":
        l = str(c2_vec[3]) if not c2_vec is None else "?"
        delta = format_float(c2_vec[4]) if not c2_vec is None else "?"

        return f"Our ($\\ell = {l}, \\Delta = {delta}$)"

    return name

def format_table(table: pd.DataFrame) -> str:
    TABLE_HEADER = r"""
\begin{tabular}{@{}lrrrrr@{}}
\toprule
Name & Seconds & GBytes & $(c_2, c_3)$ & Clusters & Score \\
\midrule
""".strip()
    TABLE_FOOTER = r"""
\bottomrule
\end{tabular}
""".strip()

    rows = []
    for i in range(len(table)):
        name = table["name"].iloc[i]
        seconds = table["seconds"].iloc[i]
        gbytes = table["gbytes"].iloc[i]
        c2_vec = table["c2_vec"].iloc[i]
        clusters = table["clusters"].iloc[i]
        score = table["score"].iloc[i]

        pretty_name = format_name(name, c2_vec)

        seconds_formatted = f"{seconds:.2f}" if not math.isnan(seconds) else "---"
        gbytes_formatted = f"{gbytes:.2f}" if not math.isnan(gbytes) else "---"

        c2_value_formatted = format_float(c2_vec[1]) if not c2_vec is None else "---"
        c3_value_formatted = format_float(c2_vec[2]) if not c2_vec is None else "---"

        clusters_formatted = f"{int(clusters)}" if not math.isnan(clusters) else "---"
        score_formatted = f"{score:.2f}" if not math.isnan(score) else "---"

        rows.append(
            f"{pretty_name} & {seconds_formatted} & "
            f"{gbytes_formatted} & "
            f"({c2_value_formatted}, {c3_value_formatted}) & "
            f"{clusters_formatted} & {score_formatted}\\\\"
        )

    return "\n".join([TABLE_HEADER] + rows + [TABLE_FOOTER])
    
def format_tables(tables: list[pd.DataFrame]) -> list[str]:
    return [format_table(table) for table in tables]
