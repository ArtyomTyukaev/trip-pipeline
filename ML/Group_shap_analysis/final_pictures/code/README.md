# Код финальных графиков — это СНИМОК, а не рабочая копия

Здесь лежат скрипты, которыми нарисованы картинки из `../` и `../../SHAP_ANALYSY.md`,
скопированные **2026-08-31**.

**Правка файла в этой папке ничего не меняет.** Ни один из них отсюда не запускается: пути внутри
вычисляются от собственного расположения файла, а импорты уходят в другие папки проекта, которых
здесь нет. Снимок нужен, чтобы через год было видно, каким кодом получены именно эти картинки.

Чинить и перерисовывать — только в канонических местах:

| файл здесь | канонический путь | что рисует |
|---|---|---|
| `plot_group_shap_heatmap.py` | `analysis/group_shap/plotting/` | 12 хитмапов |
| `plot_group_shap_budget.py` | `analysis/group_shap/plotting/` | `group_shap_budget` |
| `plot_group_shap_budget_by_experiment.py` | `analysis/group_shap/plotting/` | разбивка по экспериментам |
| `plot_model_quality.py` | `analysis/group_shap/plotting/` | две картинки качества |
| `plot_group_shap_exp_vs_splicing.py` | `analysis/group_shap/plotting/` | скаттер, драйвер |
| `group_shap_shares.py` | `analysis/group_shap/plotting/` | загрузка `data/`, палитра, `save()` |
| `plot_exp_vs_splicing.py` | **`analysis/group_ablation/plotting/`** | скаттер, сама отрисовка |

Последний живёт в другом подпроекте: `draw()` общая с ablation-версией картинки, и это сделано
нарочно — цвета групп обязаны совпадать, а запуск одной и той же функции единственный способ это
гарантировать, а не декларировать.

## Чего здесь НЕТ

Снимок неполон намеренно. За пределами остались `analysis/experiment_labels.py` (имена генов и их
порядок), `analysis/group_definitions.py`, `analysis/feature_groups.py`,
`analysis/group_ablation/plotting/{common,group_shares}.py`,
`analysis/single_feature_ablation/plotting/plot_shap_matrix.py` и
`settings/custom_functions/shap_heatmaps{,_vertical}.py`. Полная копия означала бы продублировать
половину проекта, и эта половина немедленно начала бы расходиться с оригиналом.

## Как перерисовать всё из `../`

Из корня `TRIP_SCRIPTS`, окружение `TRIPml`:

```bash
G=analysis/group_shap
python $G/plotting/plot_group_shap_heatmap.py               --outdir $PWD/$G/figures
python $G/plotting/plot_group_shap_budget.py                --outdir $PWD/$G/figures
python $G/plotting/plot_group_shap_budget_by_experiment.py  --outdir $PWD/$G/figures --fixed-order
python $G/plotting/plot_model_quality.py                    --outdir $PWD/$G/figures
python $G/plotting/plot_group_shap_exp_vs_splicing.py       --outdir $PWD/$G/figures --regression-only
```

Затем скопировать нужные файлы в `../`. Данные читаются из `analysis/group_shap/data/`
(прогон 2026-08-30, 72 модели); переменная `GROUP_SHAP_DATA` позволяет указать другую папку.

## Два флага, добавленные ради этих картинок

Оба в `plot_exp_vs_splicing.py`, оба по умолчанию сохраняют прежнее поведение — ablation-версия
картинки после их добавления перерисовывается **побайтово идентично**, это проверено по `md5sum`.

- `band=False` — снимает полосу ±20 % вокруг диагонали и три строки в углу, которые её цитировали.
- `median=False` — снимает строку `median group at ×…`.
- `tasks=[...]`, `labels=False` — одна строка панелей вместо двух и вариант без подписей точек.
