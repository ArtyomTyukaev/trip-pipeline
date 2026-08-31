"""
Do the models work at all — one metric per task, TRAIN against TEST, all 36 models.

    python analysis/group_shap/plotting/plot_model_quality.py

Пишет в `--outdir` (по умолчанию `analysis/group_shap/figures/`):

    model_quality_regression.png      коэффициент Спирмена
    model_quality_classification.png  ROC AUC

## Одна метрика на задачу, и это осознанное сужение

Воркер пишет два десятка метрик, но читать их все на одной странице нельзя — они в разных
единицах и не делят ось. Здесь по одной:

* **регрессия — Спирмен.** Ранговая корреляция, а не `R²`: у `intron_splicing` таргет доля от 0
  до 1 с потолком воспроизводимости около 0.42 (два пула одного интрона согласуются на `r ≈ 0.42`),
  и `R²` на таком таргете читается хуже ранга.
* **классификация — ROC AUC.** У `single_combined` и `intron_exp` задача бинарная, у
  `intron_splicing` — трёхклассовая, поэтому берётся `TEST_roc_auc_weighted`, взвешенный по
  размеру классов. **На оси написано просто «ROC AUC»**: для бинарной задачи weighted и есть
  обычный AUC, так что подпись верна для всех трёх датасетов, а уточнение только мешало бы.

## Откуда данные

`analysis/group_shap/data/raw/*_results.tsv` — выход прогона 2026-08-30. Отдельного «обычного»
обучения не нужно: прогон группового SHAP обучает модели **тем же кодом**, что и любой другой, и
лишь дописывает колонки `GSHAP_`/`BSHAP_` сверху. Метрики в этих файлах — обычные метрики полной
модели на всех признаках, раунд 0.

Это же и причина не читать сюда `group_ablation/<датасет>/data/summary_metrics.tsv`: та таблица
осталась от прогона на **старом** определении intergenic и старых порогах сплайсинга.

## Что нарисовано

Одна дорожка — одна модель, то есть пара (эксперимент, тип таблицы). Дорожки собраны в шесть
блоков «датасет × тип таблицы», цвет блока зелёный для intragenic и фиолетовый для intergenic —
постоянная пара проекта. Внутри дорожки гантель: **полый кружок — обучающий пул**, оптимистичный
конец, **закрашенный — отложенные 15 % баркодов**, честный. Отрезок между ними это переобучение,
и его длина напечатана справа как Δ.

Два конца подписаны **над точками первой дорожки** — `train` и `test`, с выноской вниз к маркеру.
Подпись стоит один раз на всю картинку: она объясняет условные обозначения, а не дорожку, и
повторённая 36 раз была бы шумом. Именно поэтому нижней легенды нет — она дублировала бы ровно эти
два слова.

**`test` на картинке — это отложенные 15 % баркодов.** В колонках источника конец зовётся
`HOLDOUT_Spearman` у регрессии и `TEST_roc_auc_weighted` у классификации; имена колонок разошлись
исторически, величина одна и та же, и читателю показывается одно слово.

Ось фиксирована и одинакова во всех блоках, иначе длина отрезка перестала бы быть сравнимой, а
она и есть предмет картинки. У ROC AUC нарисован пол 0.5 — это не граница шкалы, а уровень
незнания.

Интроны стоят в постоянном порядке проекта и подписаны именами генов
(`analysis/experiment_labels.py`); `single_combined` не интроны и проходят через тот же порядок
неизменными.

## Две вещи, которые здесь легко сломать

**Заголовок блока выровнен по ИЗМЕРЕННОМУ левому краю имён генов, а не по краю оси.** Имена генов
это тик-лейблы, они живут в поле слева от оси, поэтому заголовок, поставленный в `axes fraction`
от нуля, выглядит съехавшим вправо. Здесь он ставится в координатах ФИГУРЫ по фактическому
`get_window_extent()` тик-лейблов после `canvas.draw()` — то есть выравнивание пересчитывается под
любую ширину и любой шрифт, а не подобрано числом.

**`savefig.bbox = "tight"` в `common.py` делает заголовок картинки частью её ширины.** Длинный
`suptitle` растягивает холст независимо от `figsize`, и попытка сузить картинку, не тронув
заголовок, ничего не даёт. Поэтому заголовок здесь короткий, а условные обозначения ушли к точкам.
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import matplotlib.pyplot as plt                                            # noqa: E402
from matplotlib.transforms import blended_transform_factory                # noqa: E402
from group_shap_shares import (ANALYSIS, DATASETS, FIGURES, FS, FS_SMALL,  # noqa: E402
                               GRID, HAIR, INK, INK_SOFT, MM, MUTED, SURFACE, save)

sys.path.insert(0, ANALYSIS)
from experiment_labels import label, order                                 # noqa: E402

GROUP_SHAP = os.path.dirname(HERE)
RAW = os.path.join(GROUP_SHAP, "data", "raw")

# Зелёный intragenic / фиолетовый intergenic — постоянная пара проекта.
TABLE_COLOUR = {"intragenic": "#008300", "intergenic": "#4a3aa7"}

SPEC = {
    "regression": dict(train="TRAIN_full_Spearman", test="HOLDOUT_Spearman",
                       axis="Spearman correlation", lo=0.0, hi=1.0, floor=None,
                       ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0]),
    # «ROC AUC» без уточнения: у бинарных датасетов weighted совпадает с обычным AUC, поэтому
    # подпись верна для всех трёх, а слово weighted на оси только отвлекало бы.
    "classification": dict(train="TRAIN_roc_auc_weighted", test="TEST_roc_auc_weighted",
                           axis="ROC AUC", lo=0.45, hi=1.0, floor=0.5,
                           ticks=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
}

# Было 183, и `bbox="tight"` отдавал 175 мм готового холста. 142 даёт около 145 — сужение на 17 %.
# Укоротить пришлось и заголовок: пока он был длинным, он сам держал ширину, и `figsize` ни на что
# не влиял.
WIDTH_MM = 108         # две картинки встают рядом в одну строку документа
LANE_H = 0.155         # дюймов на дорожку
BLOCK_GAP = 1.5        # между intragenic и intergenic ОДНОГО датасета: строка под тип таблицы
DATASET_GAP = 3.0      # между датасетами: строка под тип таблицы плюс строка под имя датасета
HEADROOM = 3.9         # дорожек сверху: train / test стоят ВЫШЕ имени первого датасета

# Подписи train / test подняты над строкой с именем датасета. Горизонтально они не пересекаются
# (имя по центру оси, а обе подписи у правого края — первая дорожка это `single_combined`, самая
# сильная модель), но вертикально пересекались бы, и выноска вниз к маркеру это скрыть не может.
TRAIN_OFFSET, TEST_OFFSET = 3.3, 2.5

# Чисто белый фон. Проектный `SURFACE` слегка кремовый (#fcfcfb) и задан глобально в rcParams
# `common.py`; здесь он переопределяется ЛОКАЛЬНО, чтобы не трогать остальные картинки проекта.
# Тем же белым заливается полый маркер, иначе на белом фоне он читался бы как чуть серый.
WHITE = "#ffffff"

# Тонкие марки: на 108 мм прежние 1.1 pt и s=15/17 давали жирную гантель, которая забивала собой
# и без того короткий отрезок. Марка должна показывать положение, а не занимать его.
LINE_W = 0.65          # толщина отрезка train↔test
DOT_S = 7.0            # площадь маркера, pt^2
DOT_EDGE_W = 0.65      # обводка полого маркера
LEFT_IN = 0.94         # дюймов слева: держит имена генов, типы таблиц и имя датасета
                       # (при 0.86 `single_combined` заходил на график на 7.8 px — поймано
                       # проверкой в конце draw(), а не глазами)
HEAD_INSET_IN = 0.04   # отступ заголовков блока от левого края холста
RIGHT_IN = 0.30        # дюймов справа под колонку Δ
TOP_IN = 0.26          # дюймов сверху под заголовок картинки
BOTTOM_IN = 0.34       # дюймов снизу под подпись оси


def load(task):
    frames = []
    for path in sorted(glob.glob(os.path.join(RAW, "*_%s_results.tsv" % task))):
        dataset = os.path.basename(path).replace("_%s_results.tsv" % task, "")
        frames.append(pd.read_csv(path, sep="\t").assign(dataset=dataset))
    if not frames:
        raise SystemExit("нет %s/*_%s_results.tsv — сначала прогон group SHAP" % (RAW, task))
    return pd.concat(frames, ignore_index=True)


def lanes(frame):
    """Список дорожек в порядке рисования: (датасет, тип таблицы, эксперимент, строка)."""
    out = []
    for dataset in DATASETS:
        for table_type in ("intragenic", "intergenic"):
            part = frame[(frame.dataset == dataset) & (frame.table_type == table_type)]
            if part.empty:
                continue
            rows = part.set_index("experiment")
            block = [(dataset, table_type, exp, rows.loc[exp]) for exp in order(rows.index)]
            out.append(block)
    return out


def label_left_edge(fig, ax):
    """Левый край подписей оси y, в долях ширины ФИГУРЫ.

    Меряется, а не подбирается: имена генов это тик-лейблы, их ширина зависит от шрифта и от
    самого длинного имени (`PACRGL Control`), поэтому любая зашитая константа разъезжается при
    первой же смене ширины картинки или набора экспериментов.
    """
    fig.canvas.draw()
    inverse = fig.transFigure.inverted()
    return min(inverse.transform((tick.get_window_extent().x0, 0))[0]
               for tick in ax.get_yticklabels())


def draw(task, frame, outdir):
    spec = SPEC[task]
    blocks = lanes(frame)
    n_lanes = sum(len(b) for b in blocks)

    # Два разных зазора. Имя датасета пишется ОДИН раз на пару таблиц, поэтому перед вторым блоком
    # датасета нужна строка только под тип таблицы, а перед первым — ещё и под имя.
    y, cursor, marks, dataset_marks, previous = [], 0.0, [], [], None
    for block in blocks:
        dataset, table_type = block[0][0], block[0][1]
        if previous is not None:
            cursor += DATASET_GAP if dataset != previous else BLOCK_GAP
        start = cursor
        for _ in block:
            y.append(cursor); cursor += 1.0
        marks.append((table_type, start))
        if dataset != previous:
            dataset_marks.append((dataset, start))
        previous = dataset
    bottom = cursor - 1.0 + 0.5                         # низ оси: полдорожки под последней

    span = bottom + HEADROOM                            # высота оси в дорожках
    width = WIDTH_MM * MM
    height = LANE_H * span + TOP_IN + BOTTOM_IN
    fig = plt.figure(figsize=(width, height))
    fig.patch.set_facecolor(WHITE)
    ax = fig.add_axes([LEFT_IN / width, BOTTOM_IN / height,
                       1 - (LEFT_IN + RIGHT_IN) / width, 1 - (TOP_IN + BOTTOM_IN) / height])
    ax.set_facecolor(WHITE)

    flat = [row for block in blocks for row in block]
    for (dataset, table_type, exp, row), yy in zip(flat, y):
        tr, te = float(row[spec["train"]]), float(row[spec["test"]])
        colour = TABLE_COLOUR[table_type]
        ax.plot([te, tr], [yy, yy], color=colour, linewidth=LINE_W, zorder=2,
                solid_capstyle="round")
        ax.scatter([tr], [yy], s=DOT_S, facecolor=WHITE, edgecolor=colour,
                   linewidth=DOT_EDGE_W, zorder=3)
        ax.scatter([te], [yy], s=DOT_S + 1.5, color=colour, zorder=4)
        ax.annotate("%.3f" % (tr - te), xy=(1.0, yy), xycoords=("axes fraction", "data"),
                    xytext=(4, 0), textcoords="offset points", va="center", ha="left",
                    fontsize=FS_SMALL - 1, color=INK_SOFT, annotation_clip=False)

    ax.set_yticks(y)
    ax.set_yticklabels([label(exp) for _d, _t, exp, _r in flat], fontsize=FS_SMALL)
    ax.set_ylim(bottom, -HEADROOM)
    ax.set_xlim(spec["lo"], spec["hi"])
    ax.set_xticks(spec["ticks"])
    ax.set_xlabel(spec["axis"], fontsize=FS + 2)
    if spec["floor"] is not None:
        ax.axvline(spec["floor"], color=MUTED, linewidth=HAIR, zorder=1)

    # Условные обозначения — над двумя концами ПЕРВОЙ дорожки, с выноской вниз к маркеру. Один раз
    # на картинку: они объясняют, что такое полый и закрашенный кружок, а не эту конкретную модель.
    # Два разных отступа по вертикали, потому что на первой дорожке концы близки и на одной строке
    # слова слиплись бы.
    first = flat[0][3]
    for offset, value, name in ((TRAIN_OFFSET, float(first[spec["train"]]), "train"),
                                (TEST_OFFSET, float(first[spec["test"]]), "test")):
        x = float(np.clip(value, spec["lo"], spec["hi"]))
        edge = (x - spec["lo"]) / (spec["hi"] - spec["lo"])
        # У самого края оси подпись прижимается к нему, а не центрируется: иначе она уезжает в
        # колонку Δ и накрывает её заголовок.
        align = "right" if edge > 0.85 else ("left" if edge < 0.15 else "center")
        ax.annotate(name, xy=(x, y[0] - offset), ha=align, va="bottom",
                    fontsize=FS_SMALL, color=INK_SOFT, zorder=8)
        ax.plot([x, x], [y[0] - offset, y[0] - 0.28], color=MUTED, linewidth=0.5, zorder=3)

    ax.annotate("Δ", xy=(1.0, y[0] - TEST_OFFSET), xycoords=("axes fraction", "data"),
                xytext=(4, 0), textcoords="offset points", ha="left", va="bottom",
                fontsize=FS_SMALL, color=INK_SOFT, annotation_clip=False)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_linewidth(HAIR); ax.spines[side].set_color(INK)
    ax.tick_params(direction="out", width=HAIR, colors=INK, pad=2)
    # Засечка на каждый интрон: 36 дорожек подряд, и без них строка отыскивается только по
    # вертикальному положению текста. Левая линия оси включена ради них же — засечки, растущие из
    # пустоты, читаются как артефакт.
    ax.tick_params(axis="y", length=2.0, labelsize=FS_SMALL)
    # Подписи по X крупнее прочего: это единственная шкала картинки, по ней снимают значения.
    ax.tick_params(axis="x", length=2.4, labelsize=FS + 0.5)
    ax.grid(axis="x", color=GRID, linewidth=HAIR)
    ax.set_axisbelow(True)

    # Заголовок блока — над первой дорожкой блока, выровнен по левому краю имён генов, и **в две
    # строки**: датасет, под ним тип таблицы. Одной строкой `single_combined · intragenic` это
    # 28 знаков, что вдвое шире левого поля, поэтому текст уходил в область графика и пересекал
    # сетку и линию пола 0.5. Каждая половина порознь — около 15 знаков, то есть ровно ширина
    # самого длинного имени гена (`PACRGL Control`), под которое поле и рассчитано.
    # Заголовки ставятся на ФИКСИРОВАННЫЙ отступ от края холста, а не по левому краю имён генов.
    # Имена генов это тик-лейблы, они выровнены по ПРАВОМУ краю к оси, поэтому их левый край едет
    # вместе с осью: расширение `LEFT_IN` сдвигало и якорь, и наезд не уменьшался ни на пиксель.
    # Вдобавок `single_combined` длиннее самого длинного имени гена, так что при таком якоре он не
    # поместился бы ни при какой ширине поля. `label_left_edge()` оставлена — она печатается в лог
    # как справочная величина.
    left = HEAD_INSET_IN / width
    at_labels = blended_transform_factory(fig.transFigure, ax.transData)
    heads = []
    for table_type, start in marks:
        heads.append(ax.annotate(
            table_type, xy=(left, start - 0.55), xycoords=at_labels, ha="left", va="bottom",
            fontsize=FS_SMALL, color=TABLE_COLOUR[table_type], annotation_clip=False))

    # Имя датасета — ОДИН раз на пару таблиц, чёрным жирным по центру оси. Раньше оно повторялось в
    # каждом из шести заголовков, и повтор нёс ноль информации: пара intragenic/intergenic всегда
    # принадлежит одному датасету. Центрирование по оси, а не по полю, отделяет его от цветных
    # подписей типа таблицы, которые остаются слева.
    # Имя датасета — ОДИН раз на пару таблиц, обычным начертанием, в том же левом поле и по тому
    # же левому краю, что имена генов и типы таблиц. Повторять его в каждом из шести заголовков
    # незачем: пара intragenic/intergenic всегда принадлежит одному датасету.
    #
    # Жирным оно не пишется, и это не выбор оформления: зарегистрирован единственный
    # `context/other/Arial.ttf` весом 400, matplotlib на запрос bold молча отдаёт тот же Regular,
    # а настоящего Arial Bold и метрически совместимого Liberation Sans Bold на машине нет.
    for dataset, start in dataset_marks:
        heads.append(ax.annotate(
            dataset, xy=(left, start - 1.50), xycoords=at_labels, ha="left", va="bottom",
            fontsize=FS, color=INK, annotation_clip=False))

    # Утверждение, а не надежда: заголовок обязан кончаться левее оси. Ширина текста зависит от
    # шрифта и от длины имени датасета, поэтому «на глаз в этот раз поместилось» здесь ничего не
    # гарантирует — при добавлении датасета с длинным именем это должно падать, а не разъезжаться.
    fig.canvas.draw()
    axes_left = ax.get_window_extent().x0
    worst = max(h.get_window_extent().x1 for h in heads)
    if worst > axes_left:
        raise SystemExit("заголовок блока заходит на область графика на %.1f px: увеличьте LEFT_IN"
                         % (worst - axes_left))

    # «Quality of the models», а не «Models quality»: во множественном числе существительное в роли
    # определения по-английски не ставится, а картинка показывает 36 моделей сразу. Имя файла
    # остаётся `model_quality_*` — на него ссылаются документ и `final_pictures/`.
    save(fig, "model_quality_%s" % task, outdir, "Quality of the models · %s" % task)
    plt.close(fig)
    print("  %-14s дорожек %d   холст %.0f×%.0f мм   левый край подписей %.4f"
          % (task, n_lanes, width / MM, height / MM, label_left_edge(fig, ax)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=FIGURES)
    ap.add_argument("--tasks", nargs="+", default=["regression", "classification"],
                    choices=["regression", "classification"])
    args = ap.parse_args(argv)
    for task in args.tasks:
        draw(task, load(task), args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
