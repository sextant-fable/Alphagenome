# C. elegans AlphaGenome-like Data Format Report

## 1. 原论文/官方 AlphaGenome 数据形式

AlphaGenome 的核心任务可以概括为：

> 给定一段参考基因组 DNA 序列，预测这段序列上不同实验/组织/细胞状态对应的功能基因组信号轨道。

官方 research repo 说明模型实现是 JAX，并提供了从 TFRecord 读取训练数据的 loader。也就是说，原始训练数据的容器是压缩 TFRecord，但模型学习的本质数据结构是：

- 输入：DNA sequence，一般是 1 Mb 级别上下文，编码为 one-hot 矩阵。
- 标签：多个 functional genomics output tracks，例如 RNA-seq、DNase-seq、ATAC-seq、CAGE、ChIP-seq 等。
- 每个 track 有 metadata：output type、assay、biosample/tissue、ontology term、strand 等。
- 每个训练样本有 interval 信息：chromosome、start、end。
- 每类输出有 mask，表示哪些 track 在该样本/物种/上下文中有效。

我们当前项目只做其中的 RNA_SEQ 子任务，因此目标不是复刻全部 AlphaGenome 多模态输出，而是先把 C. elegans RNA-seq 数据整理成同一种“DNA sequence -> per-base RNA-seq tracks”的训练输入。

## 2. 原始数据

原始数据来自 Google Drive 文件夹和本地下载目录：

- `training_input_bigwig/`
- `Samples.xlsx`

原始信号文件：

- 20 个 `.bw` / bigWig 文件。
- 文件名为 `SRR7443583.bw` 到 `SRR7443602.bw`。
- 文件名 stem 与 Excel 表格中的 `ID` 列一一对应。

原始 metadata 内容：

- 物种：`Caenorhabditis elegans`
- 实验类型：`mRNA`，在本项目中映射为 AlphaGenome 的 `RNA_SEQ`
- 数据来源：`Warner/Murray 2019`
- 性别：`Hermaphrodite`
- 组织/细胞富集类型：
  - `Muscle (hlh-1)`：10 个样本
  - `Intestine (end1)`：8 个样本
  - `Pharynx (pha-4)`：2 个样本
- 发育阶段：
  - `Embryo T0`
  - `Embryo T1`
  - `Embryo T2`
  - `Embryo T3`
  - `Embryo T4`

BigWig 坐标系统：

- 参考基因组为 C. elegans `WBcel235`
- 染色体/contig：
  - `I`: 15072434
  - `II`: 15279421
  - `III`: 13783801
  - `IV`: 17493829
  - `V`: 20924180
  - `X`: 17718942
  - `MtDNA`: 13794

所有 20 个 bigWig 文件都通过 QC，染色体集合和长度一致。

## 3. 第一步转换：样本级 metadata 标准化

我们先把 Excel metadata 和 bigWig 文件名匹配起来，得到每个原始样本对应的标准化 track metadata：

输出文件：

- `alphagenome_custom/metadata/track_metadata.tsv`

主要字段如下：

- `organism`: `caenorhabditis_elegans`
- `output_type`: `RNA_SEQ`
- `assay`: `RNA-seq`
- `data_source`: `Warner/Murray 2019`
- `biosample_name`: tissue/cell type，例如 `Muscle (hlh-1)`
- `stage`: embryo time point，例如 `Embryo T2`
- `sex`: `Hermaphrodite`
- `sample_id`: SRR ID
- `file_path`: 对应的 bigWig 路径
- `strand`: `.`，表示 unstranded

这一步相当于把用户原始 Excel 表转成 AlphaGenome 风格的 track metadata。

## 4. 第二步转换：按 biological context 合并重复

AlphaGenome 的输出 track 更接近“一个 biological context 下的一个实验轨道”，而不是每个 SRR replicate 单独作为一个最终输出。

因此我们按以下字段分组：

- output type
- assay
- strand
- tissue/biosample
- developmental stage
- sex
- condition
- data source

20 个原始样本被合并为 11 个 RNA_SEQ tracks：

| Group | Tissue | Stage | Replicates |
| --- | --- | --- | --- |
| RNA_SEQ_001 | Intestine (end1) | Embryo T0 | 2 |
| RNA_SEQ_002 | Intestine (end1) | Embryo T1 | 2 |
| RNA_SEQ_003 | Intestine (end1) | Embryo T2 | 2 |
| RNA_SEQ_004 | Intestine (end1) | Embryo T3 | 1 |
| RNA_SEQ_005 | Intestine (end1) | Embryo T4 | 1 |
| RNA_SEQ_006 | Muscle (hlh-1) | Embryo T0 | 2 |
| RNA_SEQ_007 | Muscle (hlh-1) | Embryo T1 | 2 |
| RNA_SEQ_008 | Muscle (hlh-1) | Embryo T2 | 2 |
| RNA_SEQ_009 | Muscle (hlh-1) | Embryo T3 | 2 |
| RNA_SEQ_010 | Muscle (hlh-1) | Embryo T4 | 2 |
| RNA_SEQ_011 | Pharynx (pha-4) | Embryo T4 | 2 |

重复样本的处理方式：

- 对同一 group 内的 replicate bigWig 做逐碱基平均。
- 当前保留原 bigWig 的 existing signal scale，没有额外做 RPM/log normalization。
- 因为 Excel/metadata 中没有明确 normalization 单位，所以这一点后续需要结合原始数据说明进一步确认。

输出文件：

- 分组定义：`alphagenome_custom/metadata/track_groups.tsv`
- grouped bigWig：`alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_001.bw` 到 `RNA_SEQ_011.bw`
- grouped metadata：`alphagenome_custom/metadata/track_metadata_grouped.tsv`

最终模型的 RNA_SEQ 输出通道数就是 11。

## 5. 第三步转换：参考基因组与区间切分

为了从 DNA 序列预测信号，我们下载并准备了与 bigWig 坐标一致的参考基因组：

- FASTA：`alphagenome_custom/reference/genome.fa`
- FASTA index：`alphagenome_custom/reference/genome.fa.fai`
- GTF annotation：`alphagenome_custom/reference/annotation.gtf`

FASTA 染色体长度与 bigWig 染色体长度完全一致。

训练区间设计：

- window size：`1048576 bp`，即 `2^20`，约 1 Mb。
- stride：`524288 bp`，即 `2^19`，50% overlap。
- 采用 chromosome-level split，避免同一染色体相邻 overlapping windows 同时出现在 train 和 valid/test 中。

Split 结果：

| Split | Chromosomes | Number of intervals |
| --- | --- | --- |
| train | I, II, III, IV | 116 |
| valid | V | 39 |
| test | X | 33 |
| excluded | MtDNA | 0 |

输出文件：

- `alphagenome_custom/intervals/train.bed`
- `alphagenome_custom/intervals/valid.bed`
- `alphagenome_custom/intervals/test.bed`

## 6. 第四步转换：构建训练样本

对每个 genomic interval，我们做两件事：

1. 从参考 FASTA 取出该区间 DNA 序列，并做 one-hot 编码。
2. 从 11 个 grouped bigWig 中取出同一区间的 RNA-seq signal，并按 track 维度堆叠。

每个训练样本的逻辑结构：

```text
dna_sequence: [1048576, 4]
rna_seq: [1048576, 11]
rna_seq_mask: [1, 11]
rna_seq_strand: [1, 11]
interval_chromosome
interval_start
interval_end
```

字段解释：

- `dna_sequence`: DNA one-hot，A/C/G/T 四列；N 或非 ACGT 碱基为全 0。
- `rna_seq`: 每个位置、每个 RNA_SEQ track 的信号值。
- `rna_seq_mask`: 当前 11 条 RNA_SEQ track 都有效，所以全为 True。
- `rna_seq_strand`: 当前数据是 unstranded；在我们的 NPZ 中用内部编码保存，后续若转官方 proto/TFRecord 需要映射到官方 strand enum。
- `interval_*`: 记录该样本来自哪个基因组区间。

## 7. 当前最终保存格式

官方 AlphaGenome released training data 使用压缩 TFRecord 作为数据容器；我们当前本地没有 TensorFlow，并且后续计划用 PyTorch，所以使用框架无关的 NPZ 容器保存相同逻辑字段。

也就是说：

- 语义结构对齐 AlphaGenome：DNA one-hot + interval metadata + RNA_SEQ tracks + mask/strand。
- 文件容器暂时不是 TFRecord，而是 per-interval `.npz`。
- 后续如果必须接官方 TFRecord loader，可以从这些 NPZ 再转 TFRecord。
- 如果用 PyTorch，则可以直接读取 NPZ，不需要 TensorFlow。

最终数据目录：

- train：`alphagenome_custom/datasets/rna_seq_npz_train`
- valid：`alphagenome_custom/datasets/rna_seq_npz_valid`
- test：`alphagenome_custom/datasets/rna_seq_npz_test`

数据规模：

| Split | Examples | Disk size | RNA_SEQ dtype |
| --- | --- | --- | --- |
| train | 116 | 5.4G | float32 |
| valid | 39 | 1.8G | float32 |
| test | 33 | 1.5G | float32 |

RNA_SEQ 使用 `float32`，不是 `float16`，因为全量转换时发现部分 signal 超过 float16 的有限范围。为了避免标签溢出损坏，正式数据使用 float32。

## 8. PyTorch 读取时的格式

为了适配 PyTorch 的 Conv1d/sequence model，读取时会把数组转成 channel-first：

```text
dna_sequence: [4, 1048576]
rna_seq: [11, 1048576]
rna_seq_mask: [11, 1]
rna_seq_strand: [11]
```

batch 后变成：

```text
dna_sequence: [B, 4, 1048576]
rna_seq: [B, 11, 1048576]
rna_seq_mask: [B, 11, 1]
```

PyTorch 入口文件：

- `scripts/torch_rna_seq_dataset.py`
- `scripts/torch_smoke_train.py`

默认训练读取时对 target 做 `log1p`，原因是原始 bigWig signal 动态范围较大。原始 NPZ 文件中仍保存未 log-transform 的 float32 signal。

## 9. 质量控制结果

已经完成的 QC：

- 20 个原始 bigWig 全部可读。
- 所有 bigWig 的染色体集合和长度一致。
- 参考 FASTA 与 bigWig 染色体长度完全匹配。
- 11 个 grouped bigWig 全部可读。
- train/valid/test NPZ 全量读回成功。
- `rna_seq` 中 `NaN/Inf` 数量均为 0。

读回后的全局 signal 范围：

| Split | min | max | mean | non-finite count |
| --- | --- | --- | --- | --- |
| train | 0 | 506671 | 26.7011 | 0 |
| valid | 0 | 179308 | 22.5533 | 0 |
| test | 0 | 27789 | 14.6807 | 0 |

## 10. 可以向老师强调的结论

这一步工作的核心是把原始的 C. elegans RNA-seq bigWig 数据，转换成 AlphaGenome 风格的 sequence-to-track 训练输入。

具体来说：

1. 原始数据是 20 个样本级 bigWig signal tracks，加一个 Excel metadata 表。
2. 我们先把 Excel metadata 标准化成 AlphaGenome-like track metadata。
3. 再按 tissue/stage/assay/source 等 biological context 合并 replicate。
4. 20 个原始样本最终变成 11 个 RNA_SEQ 输出 tracks。
5. 用 WBcel235 参考基因组切 1 Mb 区间，提取 DNA one-hot。
6. 对每个 1 Mb 区间，从 11 个 grouped bigWig 中提取对应 RNA-seq signal。
7. 最终每个训练样本是 `DNA [1048576,4] -> RNA_SEQ [1048576,11]`。
8. 文件容器当前是 NPZ，方便 PyTorch；逻辑字段与 AlphaGenome 训练数据结构对齐，后续可以再转 TFRecord。
