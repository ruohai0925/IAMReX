# Transolver_v2 —— 修改后的 Transolver + 过拟合排除实验

本目录是 `../Transolver/` 的第二版，目标有两个：**给出一个容量与数据量匹配的 Transolver**，
以及**用一套可复现的实验排除过拟合的可能性**。原 `../Transolver/` 不动，结果仍可对照。

> ⚠️ **这些脚本在本机未执行过**（按要求只交付、不运行）。首次使用请先跑自检：
> `python3 transolver_slim.py`，应打印三档预设的参数量并通过形状断言。

```
Transolver_v2/
├── transolver_slim.py   # 修改后的模型：自包含 Transolver，容量可调
├── train_v2.py          # 训练脚本：三分划分 / 早停 / 正则 / 学习曲线 / 随机标签 / 外部测试
└── README.md            # 本文件
```

---

## 1. 模型改了什么

### 1.1 去掉 physicsnemo 依赖

`transolver_slim.py` 是 `physicsnemo.models.transolver.Transolver`（standard 分支、
非结构化网格）的**自包含忠实重实现**，逐层对齐原版计算图：

```
preprocess(Mlp) → N×TransolverBlock → 末层输出 out_dim
Block: fx = Attn(LN(fx)) + fx ;  fx = MLP(LN(fx)) + fx ;  末层再接 LN+Linear
PhysicsAttention(IrregularMesh):
    x_mid, fx_mid = Linear_x(x), Linear_fx(x)                → (B,N,H,D)
    slice_weights = softmax(Linear_slice(x_mid) / T, dim=-1) → (B,N,H,S)
    slice_token   = Σ_N (slice_weights / Σ_N slice_weights)·fx_mid → (B,H,S,D)
    out_slice     = SDPA(qkv(slice_token))                   → (B,H,S,D)
    out           = Linear_out( Σ_S slice_weights · out_slice )
```

温度参数初始化 0.5、前向 clamp 到 `[0.5, 5]`，切片先归一再聚合，`dim_head = n_hidden // n_head`
——均与原版一致。`forward(fx=..., embedding=...)` 签名不变，**既有推理代码可直接替换**。

只依赖 `torch`。原版链路需要 `physicsnemo → warp / s3fs / jaxtyping / einops`，
在干净环境中经常装不上（本机实测依次缺 `warp`、`s3fs`、`jaxtyping`）。

### 1.2 容量降到与数据量匹配

这是核心改动。原配置在本任务上是 **778,641 参数**，而训练集只有约 2,000~6,000 个**支撑笼**
——注意笼内 27 个点由同一个 3 维相位刚性决定（`rel_hat = 相位 + 整数偏移`），
**不是独立样本**，所以真实样本数是笼数而非 54,000 行。参数量是独立样本数的两个数量级，
这是审稿人最容易攻击的地方。

| 预设 | n_hidden | n_layers | n_head | slice_num | mlp_ratio | 用途 |
| :-- | --: | --: | --: | --: | --: | :-- |
| `orig` | 128 | 4 | 4 | 32 | 4 | 复现原配置，作对照 |
| **`slim`** | **48** | **2** | **4** | **8** | **2** | **推荐默认** |
| `tiny` | 32 | 2 | 4 | 4 | 2 | 极限压缩，测精度下限 |

精确参数量由 `python3 transolver_slim.py` 打印。也可用 `--n-hidden/--n-layers/...` 逐项覆盖。

**若 `slim` 与 `orig` 精度相当，就同时得到两个结论**：模型不是靠大容量记忆；
以及推理时延可以从原版的 ~5.8 ms/10k点 显著下降——正好解决 `../Test/benchmark_report.md`
里 Transolver 相对 MLP 唯一的劣势。

## 2. 训练脚本加了什么

| 项 | 旧版 `train_transolver_w.py` | 本版 `train_v2.py` |
| :-- | :-- | :-- |
| 数据划分 | train/val 两分 | **train/val/test 三分**，test 只在最后评估一次 |
| 划分方式 | 按 case 随机 | 按 case 随机 **或按相位卦限留出** |
| 停止条件 | 跑满 epochs | **早停**（`--patience`），恢复最优 val 权重 |
| 正则 | 无 | **AdamW weight decay**（默认 1e-4）+ **dropout** |
| 标准化 | 用训练集统计量 | 同（保持），并显式记录到 `norm_stats.npz` |
| 截断笼 | 用 edge padding 补到 27 | **直接丢弃**，避免伪造样本 |
| 诊断 | 无 | `--data-frac` 学习曲线、`--shuffle-labels` 随机标签对照 |
| 外部检验 | 无 | `--extra-test` 在从未见过的数据集上评估 |
| 报告 | R²/L2 | 加 **参数量/训练笼数之比**、**val/train 与 test/train 的 L2 比值** |

### 关于 test 曲线的说明

脚本每轮也会记录 test loss，但**仅用于画图**。早停判据、最优权重选择**只看 val**，
test 从不参与任何决策。这样 `loss_curve.png` 上能同时看到三条曲线，而 test 仍是干净的。

### 两种划分方式的区别（很重要）

- `--split case`（默认）：按算例随机分。但各 case 只差一个随机亚格偏移，
  **本质上是 iid 划分**，val 与 train 同分布。"val 不涨"只能排除经典过拟合，
  **不能证明泛化**。
- `--split octant`：把相位立方体 `[-0.5,0.5]³` 按三分量符号分成 8 个卦限，
  整块留出一个不参与训练。这是**输入空间内的真外推**，才真正有说服力。

## 3. 排除过拟合的完整实验方案

五个实验，每个回答一个具体质疑。建议按顺序跑。

### 实验 1：三分划分 + 早停（基线）

```bash
python3 train_v2.py --data ../dataset5workspace_shape/batch_data_shape_sphere \
                    --preset slim --out out_slim
```

**回答**："有没有经典过拟合？"
**看**：`val/train 相对L2 之比`≈1、test 指标与 val 相当、`loss_curve.png` 上 val 未掉头向上。

### 实验 2：容量对照

```bash
for p in orig slim tiny; do
  python3 train_v2.py --data <同上> --preset $p --out out_$p
done
```

**回答**："778k 参数 vs 几千个样本，是不是靠记忆？"
**看**：若 `slim`（约 1/20 参数量）精度与 `orig` 相当，则容量不是精度来源。

### 实验 3：相位卦限外推

```bash
python3 train_v2.py --data <同上> --preset slim --split octant --holdout-octant 0 --out out_octant
```

**回答**："换个随机划分就说泛化，是不是太宽松？"
**看**：test（留出卦限）的 R² 会比 iid 划分低——**这个落差本身就是诚实的证据**，
说明指标对划分方式敏感，而你知道敏感在哪。

### 实验 4：学习曲线（数据量够不够）

```bash
for f in 0.02 0.05 0.1 0.25 0.5 1.0; do
  python3 train_v2.py --data <同上> --preset slim --data-frac $f --out out_lc_$f
done
```

**回答**："精度高是不是因为数据量太小？"——这是**唯一能正面回答该问题**的实验。
**看**：把各 `out_lc_*/metrics.json` 的 `metrics.val.relL2` 对 `n_train_cages` 作 log-log 图。
**曲线若已平台化，则证明再加数据也不会更好，当前数据量足够。**

### 实验 5：随机标签对照

```bash
python3 train_v2.py --data <同上> --preset slim --shuffle-labels --epochs 2000 --patience 0 --out out_shuffle
```

**回答**："参数量远大于样本量，凭什么说没记忆？"（Zhang et al. 2017 的标准做法）
**看**：模型若能把随机标签也拟合到低 train loss，说明**有足够容量记忆**；
而真实标签下 train/val gap 极小 ⇒ 它**选择了学结构而非记忆**。

### 附加：不变性验证（顺带产出论文里的鲁棒性表格）

```bash
python3 train_v2.py --data ../dataset5workspace_shape/batch_data_aniso_k1 \
  --extra-test ../dataset5workspace_shape/batch_data_aniso_k3 \
               ../dataset5workspace_shape/batch_data_shape_oblate \
               ../dataset5workspace_shape/batch_data_shape_triaxial \
               ../dataset5workspace_shape/batch_data_rot_45 \
               ../dataset5workspace_shape/batch_data_rot_90 \
  --preset slim --out out_inv
```

脚本会打印各外部集的 R² 与**极差**。极差极小即证明权重映射对形状/旋转/各向异性不变。

> ⚠️ 注意措辞：这些外部集**不是 out-of-distribution 数据**（论证见
> `../dataset5workspace_shape/README.md` §5.1），所以能支撑的是**不变性**，
> **不能**写成"验证了泛化能力"。

## 4. 用哪份数据

三分划分需要足够的 case 数。**`../dataset3workspace/` 目前磁盘上只剩 2 个 case，不够用**
（它被 `-c 2` 重跑覆盖过，也与 `../MLP/out` 训练时的 10 case 对不上）。可用的：

| 数据集 | case 数 | 笼数 | 说明 |
| :-- | --: | --: | :-- |
| `../dataset5workspace_shape/batch_data_*`（9 个） | 各 10 | 各 2,000 | **推荐**，可 6/2/2 三分 |
| `../dataset3workspace_nonuniform/batch_data` | 10 | 2,000 | 可用 |
| `../dataset4workspace/batch_data` | 20 瞬间 | 4,000 | 运动数据，可做时间外推 |
| `../dataset3workspace/batch_data` | **2** | 400 | **不够三分** |

`--data` 支持传多个目录合并（会自动给 case_id 加目录名前缀防撞名）。
但注意：合并同分布的数据集只增加采样密度，不扩充分布（见 dataset5 的 README §5.1）。

## 5. 产物

```
out_xxx/
├── model_best.pt      # {"state_dict":…, "cfg":…}，cfg 内含完整模型超参，加载时无需再传
├── norm_stats.npz     # xm, xs, ym, ys —— 推理时必须复用
├── metrics.json       # 参数量/笼数比、三集指标、外部集指标、早停轮次、完整 loss 历史
└── plots/loss_curve.png
```

加载已训模型：

```python
import torch
from transolver_slim import TransolverSlim
ck = torch.load("out_slim/model_best.pt", map_location="cpu")
model = TransolverSlim(**ck["cfg"]); model.load_state_dict(ck["state_dict"]); model.eval()
```

## 6. 依赖

`torch`（CPU 即可）、`numpy`、`pandas`、`pyarrow`；`matplotlib` 可选（缺失自动跳过绘图）。
**不需要** physicsnemo / warp / s3fs / jaxtyping。

见 `requirements.txt`：`pip install -r requirements.txt`。
