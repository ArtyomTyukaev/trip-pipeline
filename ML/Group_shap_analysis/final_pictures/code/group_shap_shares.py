"""Folding `data/group_shap.tsv` into GROUP shares, once, for every figure that needs it.

    from group_shap_shares import load_group_shap, mean_group_share, spread

This is the `group_shap` counterpart of `group_ablation/plotting/group_shares.py`, and it exists
for one reason: **the fold is one step shorter here, and pretending otherwise would be wrong.**

  * `summary_shap.tsv` of the ablation carries one row per one-hot LEVEL, so that module has to sum
    levels into features and features into groups before anything may be ranked.
  * `group_shap.tsv` carries **one row per group already**. `gshap` is the exact group Shapley
    value with the group as ONE player - `mean_rows |phi_g|`, then averaged over the model's
    outputs (`settings/custom_functions/group_shap.py::_collapse`). There is nothing to sum up. A
    reader who copies the level-folding from the other module will silently sum a group into
    itself.

## Why a share and not `gshap` itself

`gshap` is in the units of the model's own margin, and those units are not shared between models: a
regression model on `log10(1 + RPM)` and a three-class splicing classifier on logits do not have a
common scale, and even two introns of the same dataset differ by the spread of their own targets.
Measured on this table, the per-model total runs from 0.081 to 1.82 - a factor of 22 - so plotting
`gshap` directly would rank the models, not the groups.

Dividing by the model's own total removes exactly that. `share` answers "what fraction of this
model's attention does the group hold", which is the same question the ablation figures ask of
their own `share` column, so the two are comparable as compositions.

**What a share still cannot say** is amount. Both budgets sum to 1, so a group holding more of the
splicing model's budget is not thereby better predicted - `../../CLAUDE.md` §6.2 records that the
splicing readout's own two pools agree at only r ~ 0.42, which caps what any splicing model can
reach. Read every figure built on this module as composition.

## The denominator of the average over experiments

`N_EXPERIMENTS`, imported from the ablation module rather than counted off this table - the same
rule and for the same reason. Here it happens to be safe either way: verified, all 72 models are
present and every one of them carries every group of its dataset, so nothing is missing. It is
imported anyway so that the two analyses cannot drift apart on the one number they must share.
"""

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GROUP_SHAP = os.path.dirname(HERE)
ANALYSIS = os.path.dirname(GROUP_SHAP)
PROJECT = os.path.dirname(ANALYSIS)

# The ablation's own plotting folder, for the palette, the typography and the member counts. It
# goes on the path FIRST because `common.py` there is what registers Arial, and every figure of
# this project has to be set in it.
ABLATION_PLOTTING = os.path.join(ANALYSIS, "group_ablation", "plotting")
sys.path.insert(0, ABLATION_PLOTTING)

from common import (DATASETS, FS, FS_SMALL, GRID, GROUP_COLORS, HAIR, INK,     # noqa: E402,F401
                    INK_SOFT, MM, MUTED, SURFACE, TABLE_TYPES, TASKS, journal)
from group_shares import (N_EXPERIMENTS, N_MEMBERS, members_label, save)      # noqa: E402,F401

DATA = os.environ.get("GROUP_SHAP_DATA") or os.path.join(GROUP_SHAP, "data")
FIGURES = os.path.join(GROUP_SHAP, "figures")

MODEL_KEY = ["dataset", "task", "table_type", "experiment"]


def load_group_shap(path=None):
    """`data/group_shap.tsv` with a `share` column added: `gshap` over the model's own total.

    `gshap` is a mean of absolute values, so it is non-negative by construction and the total is
    never zero - asserted, because a zero total would make every share of that model NaN and the
    figures would draw the hole as if it were data.
    """
    frame = pd.read_csv(path or os.path.join(DATA, "group_shap.tsv"), sep="\t")
    if (frame["gshap"] < 0).any():
        raise ValueError("gshap is a mean of |phi| and cannot be negative")
    total = frame.groupby(MODEL_KEY)["gshap"].transform("sum")
    if (total <= 0).any():
        raise ValueError("a model's group Shapley values are all zero")
    frame["share"] = frame["gshap"] / total
    return frame


def per_experiment_group_share(frame, dataset):
    """(task, table_type, experiment) x group for one dataset. Every row sums to 1.

    Done **per dataset**, never once over the whole table. `nmd` exists only in `single_combined`,
    and one unstack over all three would invent a zero-attention `nmd` column for the two intron
    datasets - a group that was never in those models would then be plotted on the axis and counted
    in every tally. Absent and zero are different things; keeping the datasets apart is what keeps
    absent absent.

    The SAME asymmetry exists one level down, between the two table types, and it is why the NaN
    check below is per table type rather than over the whole frame. `intron_annotations` owns
    `DistToIntron5`, `DistToIntron3`, `RelIntronPos` and `MANE_intron_length`, and all four are on
    `settings/columns_to_drop_in_intergenic.txt` - so the group owns no column in the intergenic
    table, was never a player there, and has no value to report. That leaves a NaN in its column
    for every intergenic row, and it is a true statement, not a hole.

    Within one table type the check is unchanged and still strict: a group missing from SOME models
    of one type is a defect, and this is what catches it. Callers must not fill these NaN with zero
    - `mean_group_share` keeps them with `min_count=1` for exactly that reason.
    """
    part = frame[frame["dataset"] == dataset]
    if part.empty:
        raise KeyError("no rows for dataset %r" % dataset)
    wide = (part.groupby(["task", "table_type", "experiment", "group"])["share"]
            .sum().unstack("group"))

    for table_type, block in wide.groupby(level="table_type"):
        present = [c for c in block.columns if block[c].notna().any()]
        if block[present].isna().any().any():
            raise ValueError("dataset %r / %s does not carry every group in every model"
                             % (dataset, table_type))

    total = wide.sum(axis=1).to_numpy()
    if not np.allclose(total, 1.0, atol=1e-9):
        raise AssertionError("a model's group shares sum to %r"
                             % total[np.argmax(np.abs(total - 1.0))])
    return wide


def mean_group_share(frame, dataset):
    """(task, table_type) x group, averaged over the dataset's experiments.

    `min_count=1` is load-bearing. Pandas sums an all-NaN group to 0.0 by default, which would turn
    "this table type never had this group" into "this group holds no attention" - the one confusion
    this whole module is built to avoid. With it, a group absent from every experiment of a
    (task, table_type) stays NaN and the figures can draw it as absent.
    """
    per = per_experiment_group_share(frame, dataset)
    return per.groupby(level=[0, 1]).sum(min_count=1) / N_EXPERIMENTS[dataset]


def spread(frame, dataset):
    """(task, table_type, experiment) x group - the per-experiment shares, for the whiskers.

    Same object as `per_experiment_group_share`; named for what the figures use it for, so the two
    call sites read as the two different questions they are.
    """
    return per_experiment_group_share(frame, dataset)


# --- compatibility layer -------------------------------------------------------------------
# Two agents independently wrote a loader for this table, `gshap_shares.py` and this file, with
# identical arithmetic - verified numerically, max |difference| 0.000e+00 on every dataset and both
# functions. The duplicate is gone; these two names are what its callers used, kept so the figures
# it drew can be reproduced without editing them further.

load_gshap = load_group_shap


def ink_on(hex_colour):
    """White or ink for a label set inside a fill, by the fill's luminance."""
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return SURFACE if (0.299 * r + 0.587 * g + 0.114 * b) < 150 else INK
