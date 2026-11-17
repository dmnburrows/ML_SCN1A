# AlphaGenome Atlas：下载指定基因的 RNA 分数

核查日期：2026-11-03。本文依据 Atlas 论文的方法、Table S1、Data Availability，及当前官方 Atlas SDK。

脚本的 API 请求尚未运行：本次会话没有配置 AlphaGenome API key。代码只查现成的 Atlas 分数，不运行模型推理。

## 运行

需要 Python 3.11 或更新版本、AlphaGenome API key，以及相应服务访问权限。

```bash
python -m pip install -U alphagenome
python download_atlas_rna.py --metadata-only --out atlas_rna_metadata
```

脚本会隐藏输入 API key，不保存 key。先打开 `atlas_rna_metadata/available_tracks.csv`，核对实际可用的 biosample、assay 和 ontology_curie。一个 biosample 可能有多个实验 tracks。

例如，下载人 SCN1A P1a 内 100 bp 的所有单碱基替换，对 SCN1A 的 RNA 效应：

```bash
python download_atlas_rna.py --gene SCN1A --region chr2:166149101-166149200 --out SCN1A_P1a_example
```

这里的坐标是 GRCh38、1-based inclusive。脚本内部转为 SDK 的 0-based half-open。此 100 bp 仅为查询示例，不是 P1a 的完整范围，也不代表 P1a 定义已经由 Atlas 验证。

指定组织时，在 metadata 中找到准确 ontology_curie 后添加 `--ontology 实际ID`。不指定则保留该基因/区间所返回的所有 tracks；不要把它们平均成不明细胞类型的单一分数。

## 输出

- `available_tracks.csv`：该 scorer 可用的全部 track metadata。
- `rows.csv`：每个分数行的 variant、gene_name、gene_id 和输入背景 interval（如 API 返回）。
- `tracks.csv`：每一列对应的 track metadata。
- `raw_scores.csv.gz`：行是 gene/variant，列是 track；数值为原始 RNA 效应分数。
- `quantile_scores.csv.gz`：API 返回校准分数时才生成，不能按表达倍数解释。
- `scores.long.csv.gz`：已合并行 metadata、track metadata、raw_score、quantile_score 的长表。
- `manifest.json`：请求条件、软件版本、检索时间及分数含义。API key 不会写入该文件。

`--scorer RNA_SEQ_ACTIVE` 可另取 active-allele 分数，其含义是 REF/ALT 平均外显子信号的较大值，不是表达 log-fold change。

## 解释

RNA_SEQ 是 GeneMaskLFCScorer：

`raw_score = ln(mean_exon_ALT + 0.001) - ln(mean_exon_REF + 0.001)`。

正值表示 ALT 的平均预测外显子 RNA 信号较高，负值较低。它是基因及 track 特异的预测；不是 TPM、原始 reads、临床致病概率，也不能直接量化 P1a/P1b/P1c 各自的转录本丰度。

`query_interval` 选择变异位置；`gene_names` 选择评价哪一个基因。把查询区间设为 100 bp，不会把 RNA 汇总范围缩为该 100 bp。

Atlas 是预计算数据，输入约 1 Mb 背景及注释 mask 已由数据生成流程确定。查询不会重新计算不同上下文或补充未包含的基因。缺失数据不能解释为零效应。

字母图的参考碱基 contribution = 三种替换 raw_score 的负平均。Atlas 还提供乘 active-allele 信号的 rescaled 版本；该图值不能指数还原为 RNA 倍数。脚本导出的是独立变异评分及 metadata，不将其冒充网页字母图的数值。

## 来源

- Atlas 论文：Cheng et al., 2026, https://doi.org/10.64898/2026.09.16.26363192
- 官方完整报告：https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf
- 方法：报告 pp.18–19；contribution 及 rescaling：pp.31–32；Data Availability：p.46；Table S1：p.69。
- Atlas SDK：https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/atlas/atlas.py
- 基因汇总实现：https://github.com/google-deepmind/alphagenome_research/blob/main/src/alphagenome_research/model/variant_scoring/gene_mask.py
- 基因外显子 mask：https://github.com/google-deepmind/alphagenome_research/blob/main/src/alphagenome_research/model/variant_scoring/gene_mask_extractor.py
- 官方下载入口：https://alphagenome.google/downloads

论文 Data Availability 将原始 Atlas features 列为 API only。下载区的 AVI 综合分、splicing 或 AVI feature attribution 文件不能替代 gene-specific RNA_SEQ 分数。门户需要登录；本次未核实登录后的导出按钮。
