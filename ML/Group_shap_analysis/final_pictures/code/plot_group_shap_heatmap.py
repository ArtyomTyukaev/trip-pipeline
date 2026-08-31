"""The published vertical SHAP heatmap, with GROUPS on the rows instead of features.

    python analysis/group_shap/code/plots/plot_group_shap_heatmap.py

Writes twelve figures, one per (dataset x task x table_type), into `--outdir`
(`analysis/group_shap/pictures/` by default):

    group_shap_heatmap_<dataset>_<task>_<table_type>.png   (+ .pdf)

## What it is

The figure family being copied is
`trip-pipeline/ML/single_combined_pos_features/figures/*_shap_*.png`: a tall heatmap with the
features down the side ordered by decreasing importance, the experiments across the top, a white
-> dark red scale on raw `mean |SHAP|`, the row labels coloured by functional group, the title
`XGBoost | <dataset> | <table_type>` over a one-line subtitle that states the scale and the
maximum, and a colourbar labelled `mean |SHAP|` on the right.

The only thing changed here is what a row is. There, a row is one model column - forty to a
hundred and thirty of them, and the figure is a metre of thread. Here a row is a whole feature
group, so the same picture is one row per group - sixteen or seventeen since 2026-08-25,
thirteen or fourteen before it - and can be read at a glance. The
columns are still the experiments and the colour still means the same quantity.

## It is the project's own drawing function, not a copy of it

`plot_shap_heatmap_vertical` from `settings/custom_functions/shap_heatmaps_vertical.py` is called
unchanged - the same function `single_feature_ablation/plotting/plot_shap_matrix.py` uses, and the
one that drew the reference figures. Nothing about the style is re-implemented, so nothing about it
can drift.

Its input contract is an `agg_df` with `experiment`, `table_type` and one `SHAP_<name>` column per
row-to-be. That contract says nothing about `<name>` having to be a model column, so
`build_group_agg` builds the same frame with one `SHAP_<group>` column per group and hands over
`gshap`. Everything downstream - the descending sort, the Reds scale, the coloured row labels, the
header, the colourbar - is that function's own code operating on a dozen-odd rows instead of a
hundred.

Three ordinary arguments of that function are given non-default values, and each is a consequence
of the row count, not a change of style:

  * `cell_h_in`, because its default row height is `height / n_rows` = 11 in / 13 rows, which
    would stretch a dozen-odd groups over the same metre the figure exists to avoid;
  * `cell_w_in` from the widest experiment name, so `5ss_noP` fits over its own column in
    `single_combined` (2 columns) without the intron panels paying for it. Since the introns are
    drawn under their gene names (`analysis/experiment_labels.py`) that is no longer enough on its
    own - eight columns of `PACRGL Control` want 14 in of grid - so a panel whose names do not fit
    turns them 90 degrees and keeps the floor width instead. `single_combined` still does not;
  * `drop_invisible=False` and `white_pct=0`. The default drops any row with no cell above the
    20th percentile, which on a 13-row matrix would silently delete two or three groups and make
    the panels hold different groups. `white_pct=0` puts the white threshold at exactly 0, which
    is also what the reference figures ended up with (`thr=0` in their subtitle), so white means
    "exactly zero" and nothing else is hidden.

The function writes its own filename, so it is called with `save=False, show=False` - which makes
it RETURN the figure - and this script writes the PNG and the PDF under the names above. That is
the same route `plot_shap_matrix.py` takes and the reason the file must not be edited.

## The one cell worth knowing about before reading the figures

`repeats` has a group Shapley of **exactly zero** in `single_combined / classification /
intragenic / 5ss_noP`: the group was a player in the game, the booster simply never split on
`RepeatGroup`. It is a real answer, not a hole, and it is drawn as one - white, inside the grid,
against the light grey this function reserves for a group the model never had. No cell of these
twelve figures is grey: within one panel every model carries every group of that panel, which is
what the check in `build_group_agg` asserts. Across panels the group sets differ - the intergenic
ones are one row shorter, having no `intron_annotations` - so an intragenic panel and its
intergenic partner are not row-for-row comparable.

The zero cannot mislead the sort either. Rows are ordered by the group's mean over the panel's
experiments, so `repeats` is placed by both its columns - `CPA_noP` is non-zero - and not by the
zero alone.

## What the colour is, and what it is not

`gshap` is the exact group Shapley value with the group as ONE player, all 2^k
coalitions enumerated (`../PLAN.md`). It is a `mean |phi|`-type quantity, non-negative, in the
units of that model's own margin - so it drops straight onto this colour scale, and a cell is
comparable with its neighbours in the same panel. It is NOT comparable across panels: a regressor
on `log10(1 + RPM)` and a three-class splicing classifier on logits have no common unit, and the
per-model total runs over a factor of 22 across this table. Each figure carries its own `max=` in
the subtitle for that reason. The share-based figures of this folder are what compares models.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# group_shap_shares puts group_ablation/plotting on sys.path and imports `common` from it, which
# is what registers Arial (this cluster has neither Arial nor Helvetica, so `font.family="Arial"`
# alone would fall through to DejaVu Sans without a word). Import it before matplotlib is used.
from group_shap_shares import (DATASETS, FIGURES, GROUP_COLORS, TABLE_TYPES,      # noqa: E402
                               TASKS, load_group_shap)
from common import FONT                                                           # noqa: E402

import matplotlib.pyplot as plt                                                   # noqa: E402
from matplotlib.colors import to_hex                                              # noqa: E402

# Файл лежит на два уровня ниже корня `group_shap` (code/plots/), поэтому корень
# берётся ДВУМЯ dirname, а не одним. Каждый скрипт считает его сам, так что при
# следующем переезде править придётся здесь же.
GROUP_SHAP = os.path.dirname(os.path.dirname(HERE))
ANALYSIS = os.path.dirname(GROUP_SHAP)
PROJECT = os.path.dirname(ANALYSIS)
sys.path.insert(0, os.path.join(PROJECT, "settings", "custom_functions"))
sys.path.insert(0, ANALYSIS)

from shap_heatmaps import _text_width_in                                          # noqa: E402
from shap_heatmaps_vertical import plot_shap_heatmap_vertical                     # noqa: E402
from experiment_labels import label, order                                        # noqa: E402
from group_labels import group_label, table_label, dataset_label, task_label

# The rotated-label geometry is not re-implemented here. `plot_shap_matrix.py` measured the two
# candidate rotations against this very drawing function and settled on 90 degrees, and its
# `label_room` / `expand_canvas` pair grows the canvas without moving a single cell that is
# already drawn. Importing them is what keeps the two figure families from drifting apart; the
# module runs nothing at import beyond its own constants.
sys.path.insert(0, os.path.join(ANALYSIS, "single_feature_ablation", "plotting"))
from plot_shap_matrix import expand_canvas, label_room, rotate_x_labels           # noqa: E402

# `common` sets savefig.bbox="tight" for the axes-based figures of that folder. This one is laid
# out in absolute inches by the plotting function and saved by `_save` with no bbox argument, so
# the rcParam would trim a canvas that is already exactly the right size - the very thing that
# function's comment says not to do. Restore the default for this script only.
plt.rcParams["savefig.bbox"] = None

MODEL = "XGBoost"


def build_group_agg(frame, dataset, task, table_type):
    """The `agg_df` the vertical heatmap wants, with one `SHAP_<group>` column per GROUP.

    One row per experiment, value `gshap`. `build_shap_matrix` groups by experiment and takes the
    mean, so one row per key is one value per cell and no aggregation happens.

    Nothing is summed. `group_shap.tsv` already carries one row per group - the group Shapley
    value with the group as a single player - so folding levels into features into groups, which
    is what the ablation's summary tables need, would here sum a group into itself.

    ONE TABLE TYPE AT A TIME, and that is load-bearing rather than tidiness. The two types do not
    carry the same groups: `intron_annotations` owns `DistToIntron5`, `DistToIntron3`,
    `RelIntronPos` and `MANE_intron_length`, all four of which sit on
    `settings/columns_to_drop_in_intergenic.txt`, so the group has no column in the intergenic
    table and was never a player there. Pivoting both types together therefore manufactured a
    column that is NaN for every intergenic row, and the check below - a good check, guarding the
    real difference between "absent" and "zero" - fired on a legitimate absence. Filtering first
    keeps the check pointed at what it was written for: a group missing from SOME models of one
    table type, which is a defect.
    """
    part = frame[(frame["dataset"] == dataset) & (frame["task"] == task)
                 & (frame["table_type"] == table_type)]
    if part.empty:
        raise KeyError("no rows for %s / %s / %s" % (dataset, task, table_type))

    wide = (part.pivot_table(index=["experiment", "table_type"], columns="group",
                             values="gshap", aggfunc="mean")
            .reset_index())
    wide.columns.name = None

    duplicated = part.duplicated(["experiment", "table_type", "group"]).sum()
    if duplicated:
        raise ValueError("%s / %s / %s: %d duplicated (experiment, table_type, group) rows"
                         % (dataset, task, table_type, duplicated))
    groups = [c for c in wide.columns if c not in ("experiment", "table_type")]
    if wide[groups].isna().any().any():
        raise ValueError("%s / %s / %s: a model is missing a group its table type carries; absent "
                         "and zero are different things and this figure cannot tell them apart"
                         % (dataset, task, table_type))

    wide = wide.rename(columns={g: "SHAP_%s" % g for g in groups})
    return wide, groups


# Widest the grid of cells may get before the names are turned on end. The rule is on the WHOLE
# grid, not on one label, because that is the thing that goes wrong: a name only costs too much
# when there are enough columns paying for it. Measured at the 15 pt default:
#
#   single_combined, 2 columns, `5ss_noP`         1.01 in per column ->  2.02 in, horizontal
#   introns, 8 columns, gene names                1.74 in per column -> 13.91 in, rotated
#   introns, 8 columns, letter codes              0.59 in per column ->  4.72 in, horizontal
#
# So 5.0 in sits in a wide gap, and the third line is the point of putting it there: revert the
# labels to the codes and the panels go back to horizontal on their own.
MAX_HEAT_WIDTH_IN = 5.0

# Narrowest a column is ever drawn, and what a rotated panel uses: its names cost it no width.
CELL_FLOOR_IN = 0.58


def cell_width(labels, font_size, floor=CELL_FLOOR_IN, pad=0.14):
    """Inches per column: whatever the widest DRAWN name needs, never below `floor`.

    Horizontal names sit above their own column, so a column narrower than its name overprints its
    neighbour; `5ss_noP` needs about 0.6 in at 15 pt. Once the names are rotated the column owes
    them no width at all and the floor is what is used - `needs_rotation` is what decides which of
    the two regimes a panel is in, and this function is only ever asked for the horizontal one.
    """
    return max(floor, _text_width_in(labels, font_size) + pad)


def needs_rotation(labels, font_size):
    """True when the horizontal names would make the grid of cells wider than it is worth.

    Measured, not assumed: the same `_text_width_in` the column width itself is computed from.
    Deciding per panel rather than once for the script is what leaves `single_combined` exactly as
    it was published - two columns of horizontal names cost 2 in and nothing is gained by turning
    them - while the eight intron columns, which would cost 14, are turned.
    """
    return len(labels) * cell_width(labels, font_size) > MAX_HEAT_WIDTH_IN


def draw(frame, dataset, task, table_type, outdir, args):
    agg, groups = build_group_agg(frame, dataset, task, table_type)

    # Two lists, and keeping them apart is the whole trick. `experiments` are the intron CODES -
    # the values in `agg["experiment"]`, the columns of the matrix the drawing function builds,
    # and what `plot_shap_heatmap_vertical` validates `experiment_order` against
    # (shap_heatmaps_vertical.py:95-97 raises on a name it cannot find). `drawn` are the gene
    # names, and they are substituted onto the axis afterwards, once the matrix is built. Handing
    # gene names in as the order would fail that check; renaming the data would break the join to
    # every table underneath.
    experiments = order(agg["experiment"].unique())
    drawn = [label(e) for e in experiments]
    rotate = needs_rotation(drawn, args.font)

    fig = plot_shap_heatmap_vertical(
        agg_df=agg,
        # The function builds its header from `dataset_name` and `table_type` only, so the task
        # has to travel inside the dataset name or the four panels of one dataset would carry the
        # same title. Title reads `XGBoost | intron_exp | regression | intragenic`.
        dataset_name="%s | %s" % (dataset, task),
        # СЫРОЕ имя, не подпись: функция фильтрует им матрицу (shap_heatmaps_vertical.py:88),
        # а не только строит заголовок. Подпись подставляется в заголовок ниже, после отрисовки.
        table_type=table_type,
        model=MODEL,
        log_scale=False,
        cell_h_in=args.cell_h,
        # Rotated names owe the column no width, so the floor is what a rotated panel gets.
        cell_w_in=(CELL_FLOOR_IN if rotate else cell_width(drawn, args.font)),
        grid_lw=args.grid_lw,
        x_tick_font_size=args.font,
        y_tick_font_size=args.font,
        title_font_size=args.title_font,
        white_pct=0.0,            # thr = 0 exactly: only a true zero is white
        drop_invisible=False,     # every panel holds every group of its dataset AND table type
        experiment_order=experiments,
        save=False, show=False,
        groups_to_features={g: [g] for g in GROUP_COLORS},
        group_colors=GROUP_COLORS)

    if fig is None:
        raise RuntimeError("%s / %s / %s: the plotting function drew nothing"
                           % (dataset, task, table_type))

    # Подзаголовок `linear scale, white_pct=0%, thr=0, max=0.168` снимается: это настройки
    # отрисовки, они одинаковы во всех двенадцати панелях, кроме `max`, а `max` и так виден на
    # шкале цвета справа.
    #
    # Убирается ЗДЕСЬ, до `label_room`, а не в конце: тот считает высоту шапки по самому нижнему
    # из `fig.texts`, и снятие подзаголовка после расчёта оставило бы поле под текст, которого
    # уже нет. Функция отрисовки в списке неприкасаемых, поэтому строку не убрать в источнике -
    # снимается готовый объект, найденный по формату, с проверкой на единственность.
    subtitles = [t for t in fig.texts if "white_pct=" in t.get_text()]
    if len(subtitles) != 1:
        raise AssertionError("%s / %s / %s: ожидался ровно один подзаголовок, найдено %d"
                             % (dataset, task, table_type, len(subtitles)))
    subtitles[0].remove()

    ax = max((a for a in fig.axes if a.images),
             key=lambda a: a.get_position().width * a.get_position().height)
    rows = [t.get_text() for t in ax.get_yticklabels()]
    if sorted(rows) != sorted(groups):
        raise AssertionError("%s / %s / %s: drew %d rows for %d groups"
                             % (dataset, task, table_type, len(rows), len(groups)))

    # Имена групп подставляются на ось ПОСЛЕ проверки выше, позиция в позицию. `rows` прочитаны с
    # самой оси, поэтому сдвиг здесь невозможен, а ключи под рисунком остаются прежними.
    # The codes drawn by the function are replaced by the gene names, position for position. The
    # assertion is what makes that safe: it is the only thing standing between a relabelled axis
    # and a mislabelled one, and a silent off-by-one here would put every column under the wrong
    # gene while the figure looked perfectly well-formed.
    drawn_codes = [t.get_text() for t in ax.get_xticklabels()]
    if drawn_codes != experiments:
        raise AssertionError("%s / %s / %s: asked for columns %s and the function drew %s"
                             % (dataset, task, table_type, experiments, drawn_codes))
    # Только если они действительно различаются. С 2026-08-31 переименовываются и две библиотеки
    # `single_combined` (`5ss_noP` -> `Regular`, `CPA_noP` -> `PCPA`), так что панелей без
    # переименования больше не осталось, но условие оставлено: оно защищает от лишней перерисовки.
    if drawn != drawn_codes:
        ax.set_xticklabels(drawn, fontsize=args.font)

    if rotate:
        rotate_x_labels(ax)

    # ПОРЯДОК ЗДЕСЬ ОБЯЗАТЕЛЕН, и каждый шаг стоит там, где стоит, по своей причине.
    #
    # 1. `check_colors` первым, потому что он сверяет подписи строк с КЛЮЧАМИ палитры групп: после
    #    переименования он упал бы на верном рисунке.
    # 2. Переименование строк - позиция в позицию, `rows` прочитаны с самой оси выше.
    # 3. Расширение холста ПОСЛЕ переименования, а не только при повороте колонок. Имена Игоря
    #    заметно шире ключей (`Distance to intron ends` против `intron_annotations`), левое поле
    #    было рассчитано под ключи, и подписи уезжали за край.
    # 4. `check_labels` последним, на окончательном состоянии. Раньше он стоял до переименования и
    #    честно проверял КОРОТКИЕ подписи, поэтому обрезку длинных пропускал - проверка должна
    #    видеть то, что уйдёт в файл.
    check_colors(ax, "%s/%s/%s" % (dataset, task, table_type))
    ax.set_yticklabels([group_label(r) for r in rows], fontsize=args.font)
    # Заголовок функция строит как `model | dataset_name | table_type` из СЫРОГО имени. Меняем его
    # текст после отрисовки, найдя по точному совпадению: если функция когда-нибудь изменит формат,
    # совпадения не будет и это упадёт, а не подменит молча что-то другое.
    want = "%s | %s | %s | %s" % (MODEL, dataset, task, table_type)
    titles = [t for t in fig.texts if t.get_text() == want]
    if len(titles) != 1:
        raise AssertionError("%s / %s / %s: ожидался ровно один заголовок %r, найдено %d"
                             % (dataset, task, table_type, want, len(titles)))
    # Без `MODEL`: модель во всех двенадцати панелях одна и та же, слово `XGBoost` в заголовке
    # ничего не различает. Искомая строка `want` выше собирается ВМЕСТЕ с ним, потому что функция
    # рисует именно так - убирается он только из того, что остаётся на картинке.
    # Две строки: что предсказывается - сверху, условия прогона - под ним. Одной строкой это
    # `Single-exon expression prediction | Regression | Inside MANE genes`, то есть шире панели.
    titles[0].set_text("%s\n%s | %s"
                       % (dataset_label(dataset), task_label(task), table_label(table_type)))

    # Поля считаются ПОСЛЕ подмены заголовка. Он стал двухстрочным, то есть выше, и `label_room`
    # меряет высоту шапки по фактическим `fig.texts` - посчитай раньше, и вторая строка вырастет
    # в уже занятое место, ровно на повёрнутые подписи колонок.
    extra_top, extra_right, extra_left = label_room(fig, ax)
    extra_left = max(extra_left, row_label_room(fig, ax))
    expand_canvas(fig, extra_top, extra_right, extra_left)

    check_labels(fig, ax, "%s/%s/%s" % (dataset, task, table_type))


    name = "group_shap_heatmap_%s_%s_%s" % (dataset, task, table_type)
    os.makedirs(outdir, exist_ok=True)
    for ext, dpi in (("png", args.dpi), ("pdf", None)):
        fig.savefig(os.path.join(outdir, "%s.%s" % (name, ext)), dpi=dpi, facecolor="white")
    w, h = fig.get_figwidth(), fig.get_figheight()
    print("  %s.png  %.2f x %.2f in = %d x %d px at %d dpi, %d rows x %d columns, "
          "row label %.1f px when embedded at %d px wide"
          % (name, w, h, round(w * args.dpi), round(h * args.dpi), args.dpi,
             len(rows), len(experiments), args.font / (w * 72.0) * args.embed_px,
             args.embed_px))
    print("    rows: %s" % ", ".join(rows))
    plt.close(fig)
    return name


def check_colors(ax, label):
    """Every row label is drawn in its group's colour from `analysis/feature_groups.py`.

    The colouring is done inside the plotting function, by `_x_label_colors`, which falls back to
    black for anything it cannot place - so a group the palette does not know would come out black
    and look deliberate. This reads the colour back off the drawn label instead of trusting the
    dict that went in, and it is the check that keeps these figures comparable with the rest of
    the group figures, which take the same palette.
    """
    wrong = []
    for tick in ax.get_yticklabels():
        group = tick.get_text()
        want = GROUP_COLORS.get(group)
        got = to_hex(tick.get_color())
        if want is None or got.lower() != want.lower():
            wrong.append((group, got, want))
    if wrong:
        raise AssertionError("%s: %d row label(s) not in the group palette, first %s"
                             % (label, len(wrong), wrong[:3]))


def row_label_room(fig, ax, pad=0.08):
    """На сколько дюймов холст должен вырасти СЛЕВА, чтобы подписи строк поместились.

    `label_room` из `plot_shap_matrix.py` меряет только подписи КОЛОНОК - он написан под поворот
    оси X, и строки в нём не учтены вовсе. Пока строки назывались ключами (`intron_annotations`),
    поле под них хватало и без измерения; с именами Игоря (`Distance to intron ends`) перестало.

    Меряется по отрисованной фигуре, как и всё остальное здесь, поэтому держится при любом шрифте
    и любой длине имени.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [lbl.get_window_extent(renderer) for lbl in ax.get_yticklabels() if lbl.get_text()]
    if not boxes:
        return 0.0
    return max(0.0, pad - min(b.x0 for b in boxes) / fig.dpi)


def check_labels(fig, ax, label):
    """No two tick labels of either axis may overlap, and none may leave the canvas.

    Cheap, and the only thing that can catch the failure this geometry is exposed to: the row
    labels are drawn at a fixed font size into a margin sized from the same font, and the column
    labels sit horizontally over columns whose width is computed, not measured.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for axis, ticks in (("y", ax.get_yticklabels()), ("x", ax.get_xticklabels())):
        boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in ticks if t.get_text()]
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i][1], boxes[j][1]
                if a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1:
                    raise AssertionError("%s: %s labels %r and %r overlap"
                                         % (label, axis, boxes[i][0], boxes[j][0]))
        # Заголовок против подписей. `check_labels` смотрел только подписи друг с другом и с краем
        # холста, поэтому наезд ДВУХСТРОЧНОГО заголовка на повёрнутые имена колонок он пропускал -
        # а это ровно то, что здесь один раз и случилось.
        for text in fig.texts:
            tb = text.get_window_extent(renderer)
            for name, b in boxes:
                if tb.x0 < b.x1 and b.x0 < tb.x1 and tb.y0 < b.y1 and b.y0 < tb.y1:
                    raise AssertionError("%s: заголовок %r наезжает на подпись %r"
                                         % (label, text.get_text().replace("\n", " / "), name))

        page = fig.get_window_extent(renderer)
        clipped = [t for t, b in boxes
                   if b.x0 < page.x0 - 0.5 or b.x1 > page.x1 + 0.5
                   or b.y0 < page.y0 - 0.5 or b.y1 > page.y1 + 0.5]
        if clipped:
            raise AssertionError("%s: %s labels clipped by the canvas, first %s"
                                 % (label, axis, clipped[:3]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=FIGURES)
    ap.add_argument("--datasets", nargs="+", default=DATASETS, choices=DATASETS)
    ap.add_argument("--cell-h", type=float, default=0.40, dest="cell_h",
                    help="inches per group row")
    ap.add_argument("--font", type=int, default=15, help="tick label size in points")
    ap.add_argument("--title-font", type=int, default=18, dest="title_font")
    ap.add_argument("--grid-lw", type=float, default=1.0, dest="grid_lw",
                    help="width of the white cell separators, in points")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--embed-px", type=int, default=230, dest="embed_px",
                    help="width the figure is expected to be embedded at, for the legibility "
                         "line printed per figure. Reporting only, nothing is drawn from it")
    args = ap.parse_args(argv)

    print("font in use: %s" % FONT)
    if FONT != "Arial":
        raise RuntimeError("expected Arial, matplotlib resolved %r" % FONT)

    frame = load_group_shap()
    if frame["gshap"].isna().any():
        raise ValueError("gshap carries NaN; this figure cannot draw an unknown")

    zeros = frame[frame["gshap"] == 0]
    print("exact zeros in gshap: %d row(s)%s"
          % (len(zeros), (" - " + ", ".join("%s/%s/%s/%s %s" % tuple(r) for r in
                                            zeros[["dataset", "task", "table_type", "experiment",
                                                   "group"]].to_numpy())) if len(zeros) else ""))

    written = []
    for dataset in args.datasets:
        for task in TASKS:
            for table_type in TABLE_TYPES:
                written.append(draw(frame, dataset, task, table_type, args.outdir, args))
    print("\n%d figures written to %s" % (len(written), args.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
