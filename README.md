# Tumor growth expert test code

国家重点研发计划青年科学家项目“肿瘤生长模型中的模型分析和数值方法”（2021YFA1001200）的选定程序。A组负责人周珍楠，B组负责人刘鎏。

本仓库提供源码、固定配置、必要输入和运行说明。完整数学技术报告和会前结果由项目组另行提供，不属于运行代码的依赖。新计算结果在本机生成，未将大批历史输出收入仓库。

## 测试范围

| 编号 | 主题 | 当前入口与范围 |
|---|---|---|
| T01 | 径向传播及平面行波 | [A组运行说明](A/README.md)，包含现场短演示 |
| T02 | 平面前沿的稳定与失稳 | A组基线及参数/波数补算 |
| T03 | 坏死核心与内外边界 | A组核心、初生及可选空间初生探测 |
| T04 | 有限指数与自由边界参考 | A组有限指数、加密和约化独立核验 |
| T05 | 随机输入的统计传播 | [B/T05](B/T05/README.md)，单随机输入凹形初态 |
| T06 | 双保真、RFPS及三保真 | [B/T06](B/T06/README.md)，圆形初态营养扰动算例 |
| T07 | 贝叶斯参数反演 | [参考代码及限制](B/T07/README.md)，尚未作为核验完成的反演入口 |
| T08 | 数据识别与预测 | [数据转换与参考Notebook](B/T08/README.md)，尚未核验完整训练及独立预测 |

各项实际采用的配置、结果来源和检查范围见[结果索引](docs/RESULT_INDEX.md)。文件清单见[MANIFEST.json](MANIFEST.json)，当前限制见[KNOWN_ISSUES.md](KNOWN_ISSUES.md)。

## 最快开始

在仓库根目录执行，第一条只需要Python标准库：

```text
python scripts/verify_materials.py
python A/run.py --list
```

A组建议使用Python 3.12独立环境，并安装已有依赖清单：

```text
python -m venv .venv
# Windows PowerShell: .venv/Scripts/Activate.ps1
# Linux: source .venv/bin/activate
python -m pip install -r A/requirements.lock.txt
python A/demo.py --prepare
python A/demo.py --run
python A/demo.py --show-latest
```

T01/T03的C++计算在Linux中执行；Windows使用WSL2的`Ubuntu-24.04`。会前安装GCC、CMake、GSL、Eigen3：

```bash
sudo apt update
sudo apt install build-essential cmake libgsl-dev libeigen3-dev
```

T05/T06使用MATLAB。每次从新MATLAB会话进入相应目录，只加载该项求解器，避免两个同名`run_experiment`或`experiment_config`混用。例如：

```matlab
cd('B/T05')                  % 从仓库根目录开始
run_experiment('short')      % T=0.01，计算SG与16点SC并生成图
```

T06的完整默认计算规模较大，快速试用请显式缩小样本数：

```matlab
cd('B/T06')                  % 从另一个新的MATLAB会话、仓库根目录开始
opts = struct('T',0.031,'n_train',24,'n_pool',24,'n_eval',8,'n_stats',12);
run_experiment('fresh',8,opts)
run_reference(8,opts)
```

上述T06缩小流程已在R2026a检查，m=8总计约7分钟、m=64约13分钟，包含绘图和文件写入；设备不同耗时会变化。`short`模式本身只缩短时间，不自动减少样本数。完整500/1200样本计算不能按此耗时推定，也不默认作为现场短演示。

## 输出及复测

A组输出到`A/baseline/code/runs`、各测试`results`或`A/supplement/runs`；演示操作记录在`A/demo_runs`。B组输出到各测试目录的`results`或Notebook工作目录。它们默认被Git忽略。

T05/T06原入口会写固定名称的结果，正式计算和短试用应在不同克隆目录运行，或先单独归档输出，防止覆盖。代码、配置或输入变更后需重新生成清单；源码校验成功不等同于专家测试已通过。

现场使用会前准备好的本地副本和环境，无需现场下载依赖。专家核验记录、意见及签署在实际测试后归档。

## 来源

保留各来源程序的作者声明；来源与整理改动见[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)和[PROVENANCE.json](PROVENANCE.json)。这些说明不为第三方代码重新指定许可证。
