#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export projected hierarchical blockmodel partitions from graph-tool results.

This script reads one or more pickled ``NestedBlockState`` results, projects the
partition at each hierarchy level back onto the original graph, and writes a CSV
with one row per original node and one column per hierarchy level.

Examples
--------
python export_hierarchy_partitions.py \
    ./Output/res_HSBM-withDC-MCMC.pkl

python export_hierarchy_partitions.py \
    ./Output/res_HSBM-withDC-MCMC.pkl \
    ./Output/res_HSBM-woutDC-MCMC.pkl \
    --output-dir ./Output/partitions \
    --node-prop name
"""

from __future__ import annotations

import argparse
import csv
import pickle
from pathlib import Path
from typing import List, Optional

import numpy as np
from graph_tool.all import Graph, GraphView


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read graph-tool NestedBlockState result files and export the "
            "partition at every hierarchy level, projected onto the original "
            "graph, as CSV."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more pickled graph-tool NestedBlockState result files.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory where CSV files will be written. Defaults to the input "
            "file directory."
        ),
    )
    parser.add_argument(
        "--node-prop",
        default=None,
        help=(
            "Vertex property to use as node identifier in the CSV. If omitted, "
            "the script tries 'name', 'label', 'id', and 'vertex_id' in that order; "
            "otherwise it falls back to the integer vertex index."
        ),
    )
    parser.add_argument(
        "--no-node-label",
        action="store_true",
        help="Do not include a node label column; export only vertex_index and level columns.",
    )
    parser.add_argument(
        "--drop-trivial-top",
        action="store_true",
        help="Drop the highest hierarchy level if it contains a single block.",
    )
    return parser.parse_args()


def load_state(path: Path):
    with path.open("rb") as fh:
        state = pickle.load(fh)
    if not hasattr(state, "project_partition") or not hasattr(state, "get_levels"):
        raise TypeError(
            f"{path} does not look like a graph-tool NestedBlockState result."
        )
    return state


def get_base_graph(state) -> Graph:
    g = state.g
    if isinstance(g, GraphView):
        return Graph(g, prune=False)
    return g


def choose_node_property_name(g: Graph, requested: Optional[str]) -> Optional[str]:
    available = list(g.vp.keys())
    if requested is not None:
        if requested not in g.vp:
            raise KeyError(
                f"Requested node property '{requested}' was not found. "
                f"Available vertex properties: {available}"
            )
        return requested

    for candidate in ("name", "label", "id", "vertex_id"):
        if candidate in g.vp:
            return candidate
    return None


def vertex_labels(g: Graph, prop_name: Optional[str]) -> List[str]:
    if prop_name is None:
        return [str(int(v)) for v in g.vertices()]
    vp = g.vp[prop_name]
    return [str(vp[v]) for v in g.vertices()]


def projected_partition_array(state, level: int) -> np.ndarray:
    """Return membership at hierarchy level `level`, projected onto base graph."""
    projected = state.project_partition(level, 0)
    if hasattr(projected, "fa"):
        return np.asarray(projected.fa, dtype=np.int64).copy()
    return np.asarray(projected, dtype=np.int64).copy()


def level_indices_to_export(state, drop_trivial_top: bool) -> List[int]:
    levels = state.get_levels()
    idx = list(range(len(levels)))
    if drop_trivial_top and idx:
        top = idx[-1]
        if levels[top].get_nonempty_B() == 1:
            idx = idx[:-1]
    return idx


def build_rows(state, node_prop_name: Optional[str], include_node_label: bool, drop_trivial_top: bool):
    g = get_base_graph(state)
    labels = vertex_labels(g, node_prop_name)
    n = g.num_vertices()

    level_ids = level_indices_to_export(state, drop_trivial_top)
    partitions = {
        f"level_{level}": projected_partition_array(state, level)
        for level in level_ids
    }

    for key, arr in partitions.items():
        if len(arr) != n:
            raise ValueError(
                f"Projected partition '{key}' has length {len(arr)} but the graph has {n} vertices."
            )

    fieldnames = ["vertex_index"]
    if include_node_label:
        fieldnames.append("node")
    fieldnames.extend(partitions.keys())

    rows = []
    for i in range(n):
        row = {"vertex_index": i}
        if include_node_label:
            row["node"] = labels[i]
        for col, arr in partitions.items():
            row[col] = int(arr[i])
        rows.append(row)

    return fieldnames, rows


def output_path_for(input_path: Path, output_dir: Optional[Path]) -> Path:
    parent = output_dir if output_dir is not None else input_path.parent
    parent.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    return parent / f"{stem}_projected_partitions.csv"


def export_one(input_path: Path, output_dir: Optional[Path], node_prop: Optional[str], no_node_label: bool, drop_trivial_top: bool) -> Path:
    state = load_state(input_path)
    g = get_base_graph(state)
    node_prop_name = choose_node_property_name(g, node_prop)

    fieldnames, rows = build_rows(
        state=state,
        node_prop_name=node_prop_name,
        include_node_label=not no_node_label,
        drop_trivial_top=drop_trivial_top,
    )

    out_path = output_path_for(input_path, output_dir)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return out_path


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else None

    for raw_path in args.inputs:
        input_path = Path(raw_path).expanduser().resolve()
        out_path = export_one(
            input_path=input_path,
            output_dir=output_dir,
            node_prop=args.node_prop,
            no_node_label=args.no_node_label,
            drop_trivial_top=args.drop_trivial_top,
        )
        print(f"Wrote {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
