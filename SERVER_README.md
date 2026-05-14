# EEHPC 集群使用指南

## 集群概览

| 项目 | 说明 |
|------|------|
| **集群名称** | eehpc |
| **登录节点** | login1 |
| **操作系统** | Rocky Linux 9.4 |
| **调度系统** | Slurm |
| **最大作业时长** | 7 天 (gpu3: 30天) |
| **用户组** | eehpc, team1 |

---

## 存储配置

| 存储位置 | 容量 | 配额 | 用途 |
|----------|------|------|------|
| `/home/$USER` | 204 TB (NFS) | 300 GB/用户 | 代码、小文件 |
| `/data/team1/$USER` | 146 TB (NFS) | 500 GB/team1用户 | **大型数据集、模型权重** |

> **本项目数据路径**: `/data/team1/zelinli6/Alphagenome/`

---

## GPU 节点配置

| 分区 | 节点 | GPU 型号 | 显存 | GPU数/节点 | CPU核/节点 |
|------|------|----------|------|-----------|-----------|
| **gpu1** (默认) | gpu1 | RTX 2080 Ti | 11 GB | 8 | 32 |
| **gpu1** | gpu2 | RTX 4080 SUPER | 16 GB | 8 | 32 |
| **gpu1** | gpu3 | RTX 4090 | 24 GB | 8 | 48 |
| **gpu2** | gpu11-12 | A100 80GB PCIe | 80 GB | 4 | 64 |
| **gpu3** (team1保留) | gpu13-16 | A100 40GB PCIe | 40 GB | 3 | 256 |

### A100 资源汇总

| 分区 | 节点数 | 总GPU数 | 单卡显存 | 可跨节点 |
|------|--------|---------|----------|----------|
| gpu2 | 2 | 8 张 | 80 GB | 否 (单节点最多2张) |
| gpu3 | 4 | 12 张 | 40 GB | 是 (跨4节点最多12张) |

---

## QoS 资源限制

| QoS | 分区 | 最大GPU/作业 | 最大CPU/作业 | 最大作业数/用户 | 最大提交数/用户 |
|-----|------|-------------|-------------|----------------|----------------|
| rtx4 | gpu1 | 4 | 32 | 4 | 6 |
| a100 | gpu2 | 2 | 16 | 2 | 3 |
| cpu1 | cpu1 | - | 128 | 4 | 6 |

---

## 重要规则

### 1. 禁止在登录节点运行程序
login1 是登录节点，没有任何计算能力。在上面跑程序会导致节点崩溃，影响所有人。提交作业前务必确认自己在计算节点上。

### 2. 资源用完及时释放
- `sbatch` 作业结束后自动释放资源
- `srun --pty bash` 获取的交互式 shell **不会自动释放**，必须手动 `exit` 或 `scancel <jobid>`
- `srun python xxx.py` 进程结束后自动释放

### 3. 计算节点无外网
计算节点无法访问互联网。模型、数据、Python 包必须在登录节点提前下载好。

---

## 提交作业方式

### 方式一：sbatch（推荐，生产环境首选）

创建脚本 `train.sh`：

```bash
#!/bin/bash
#SBATCH -p gpu2                    # 分区
#SBATCH --gres=gpu:2               # 申请2张GPU
#SBATCH --cpus-per-task=16         # CPU核数（必须多申请，默认只有2）
#SBATCH --job-name=my_train        # 作业名称
#SBATCH --output=logs/%j.out       # 标准输出日志
#SBATCH --error=logs/%j.err        # 错误日志

# 加载环境
module load cuda/12.1.0
source /data/team1/zelinli6/miniconda/bin/activate your_env

# 运行
python train.py
```

提交：
```bash
sbatch train.sh
```

### 方式二：srun 直接运行（调试/单次任务）

```bash
# 正确做法：一步到位
srun -p gpu2 --gres=gpu:2 --cpus-per-task=16 python train.py

# 运行脚本
srun -p gpu2 --gres=gpu:2 --cpus-per-task=16 bash scripts/my_script.sh
```

```bash
# 申请2张A100 80GB
srun -p gpu2 --gres=gpu:2 --cpus-per-task=16 python train.py

# 申请跨节点多卡 (gpu3, 最多12张40GB A100)
srun -p gpu3 --gres=gpu:3 --cpus-per-task=32 python distributed_train.py
```

### 方式三：交互式 Shell（仅调试用）

```bash
salloc -p gpu2 --gres=gpu:1 --cpus-per-task=8
# 进入分配节点后：
python debug.py
# 调试完毕必须退出！
exit
```

> **注意**: `srun --pty bash` 获取 shell 后不会自动释放资源，必须 `exit` 或 `scancel <jobid>`！

---

## 常用命令速查

```bash
# 查看分区和节点状态
sinfo
sinfo -Nel                  # 详细信息

# 查看作业队列
squeue                       # 所有作业
squeue -u $USER              # 我的作业

# 取消作业
scancel <JOBID>

# 查看GPU使用情况（在计算节点内）
nvidia-smi

# 查看CPU信息（在计算节点内）
lscpu
```

---

## --cpus-per-task 说明

**务必手动指定 `--cpus-per-task`**。Slurm 默认只给 2 个 CPU 进程，远不够驱动 GPU 计算。普通个人电脑都有 16 进程，建议：

| 场景 | 推荐 --cpus-per-task |
|------|---------------------|
| 单GPU轻量任务 | 8 |
| 单GPU训练 | 16 |
| 多GPU训练 | 16-32 |
| gpu3 多节点 | 32-64 |

---

## GPU 调用示例

### gpu1 (RTX 2080Ti / 4080 / 4090)

```bash
# RTX 4090（指定节点 gpu3）
srun -p gpu1 --gres=gpu:2 --cpus-per-task=16 -w gpu3 python train.py

# RTX 4080
srun -p gpu1 --gres=gpu:1 --cpus-per-task=8 python train.py
```

### gpu2 (A100 80GB, 2节点8张)

```bash
# 单节点双卡
sbatch -p gpu2 --gres=gpu:2 --cpus-per-task=16 train.sh

# 单卡
srun -p gpu2 --gres=gpu:1 --cpus-per-task=16 python evaluate.py
```

### gpu3 (A100 40GB, 4节点12张，team1专用，可跨节点)

```bash
# 跨节点多卡分布式训练
sbatch -p gpu3 --gres=gpu:3 --cpus-per-task=32 --nodes=4 train_dist.sh

# 单节点三卡
srun -p gpu3 --gres=gpu:3 --cpus-per-task=32 python train.py
```

---

## 环境配置

### CUDA 模块

```bash
module load cuda/12.1.0    # 推荐
module load cuda/13.0.2    # 最新（默认）
module unload cuda          # 卸载
module purge                # 清除所有模块
```

可用版本: `10.2.2`, `11.0.2`, `11.3.0`, `11.8.0`, `12.1.0`, `12.8.0`, `12.9.1`, `13.0.2`

### Conda 环境

```bash
source /data/team1/zelinli6/miniconda/bin/activate your_env_name
```

---

## 作业管理最佳实践

1. **记住 JobID**: 提交作业后记录返回的 JobID，方便追踪和取消
2. **不用就 scancel**: 及时取消不再需要的作业，释放资源给其他人
3. **多用 sbatch**: 生产环境用 `sbatch`，自动释放、有日志、可排队
4. **检查节点**: 跑程序前确认 `hostname` 输出是 `gpuXX` 而不是 `login1`
5. **日志输出**: sbatch 脚本中将日志输出到 `logs/` 目录便于回溯

---

## 常见错误

| 错误 | 原因 | 解决 |
|------|------|------|
| 在 login1 跑程序卡死 | 登录节点无GPU | 用 `srun/sbatch` 提交到计算节点 |
| GPU内存不足 (OOM) | batch size太大/模型太大 | 减小batch size或用A100 |
| 网络不通 | 计算节点无外网 | 在login1预下载数据/模型/包 |
| 训练很慢 | CPU核不够 | 增加 `--cpus-per-task` |
| 资源未释放 | srun --pty bash后未exit | `scancel <jobid>` |

---

## 本项目结构

```
Alphagenome/
├── scripts/                    # 训练/数据处理脚本
│   ├── build_rna_seq_npz_dataset.py
│   ├── torch_rna_seq_dataset.py
│   ├── torch_smoke_train.py
│   ├── prepare_reference_and_intervals.py
│   ├── average_bigwig_groups.py
│   ├── download_gdrive_bigwigs.py
│   ├── qc_bigwigs.py
│   ├── summarize_bigwig_signals.py
│   ├── read_rna_seq_npz_dataset.py
│   └── smoke_test_training_input.py
├── training_input_bigwig/       # BigWig训练数据 (20个SRR文件)
├── alphagenome_custom/          # 自定义配置/数据
│   ├── intervals/
│   ├── metadata/
│   ├── reference/
│   ├── datasets/
│   └── tracks/
├── shared/                      # 共享资源
├── remote_inventory/            # 远程数据清单
├── Samples.xlsx                 # 样本表
└── SERVER_README.md             # 本文件
```
