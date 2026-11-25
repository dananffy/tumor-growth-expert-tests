# T08输入数据和参考Notebook

`data/`提供实际输入`2D_RS_Real_1_t.xlsx`与对应NPZ。坐标为11×11网格，时刻为0.25、0.375、0.5、0.625、0.75；`u_sol`维度为(y,x,t)，大小11×11×5，数值为0或1。

## 转换数据

从仓库根目录执行：

```text
python -m pip install numpy pandas openpyxl
python B/T08/excel2npz.py
```

默认输入为本目录`data/2D_RS_Real_1_t.xlsx`，新输出写入本目录`outputs/`。可用`--input`及`--output`指定其他文件；与现有NPZ逐元素一致的转换已检查。生成的自定义输出应存到独立工作目录。

## 参考训练程序

`reference/_v_predict.ipynb`及`reference/deeponet_inverse_tumor.ipynb`分别为PINN和DeepONet参考代码。Notebook可从仓库根目录、本目录或`reference/`作为内核工作目录打开；数据路径自动指向包内输入，不强制挂载Google Drive。依赖示例见`requirements.txt`，尚非已完成完整训练的锁定环境。

本轮仅整理可移植路径、清除历史输出，并修正DeepONet的pandas导入次序、重复覆盖CSV及数组标题格式。未改变损失模型、训练时段或时间原点，未重新训练网络。

## 当前范围

- PINN筛选到0.75，与报告中训练截止0.5的主配置不同。
- DeepONet注释描述时间平移，但数据坐标尚未实际平移，初值损失使用t=0；采用前需统一。
- 输入没有0.875及1.0两个验证时刻，当前验证掩码选不到数据。尚未交付与报告两时刻半径误差对应的完整预测记录。
- 本示例数据及网络输入维度不直接等同论文的61×61配置。

因此本目录可用于输入转换和实现审阅，暂不作为完整独立预测已复现的证据。选定配置的训练/推断、观测基准和预测输出可随后补入；已有完整记录能够对应时不必重新训练所有实验。

数据来自合作者随本项目提供的T08示例；关联研究为“Data-Driven Parameter Identification for Tumor Growth Models”，原论文入口：https://arxiv.org/abs/2511.15940 。公开数据的来源说明和许可由项目组随采用材料维护。
