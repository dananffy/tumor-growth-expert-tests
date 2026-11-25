# A组运行说明

本目录保留已采用的基线源码及11月10日补算程序，不依赖原项目文件夹中的结果、Python安装包或云盘绝对路径。Python科学计算依赖另按`requirements.lock.txt`安装。

所有以下命令从仓库根目录执行。可在目录外通过相应脚本的绝对路径启动。

## 基线程序

```text
python A/run.py --list
python A/run.py --case T01_radial_64
python A/run.py --case T01_planar_wave
python A/run.py --case T02_planar_33
python A/run.py --case T03_core_128_smoke
python A/run.py --case T03_core_256
python A/run.py --case T03_emergence_256
python A/run.py --case T04_autophagy_dx0.04
```

`--list`列出全部配置。T03烟雾算例只到t=0.02；其余完整大网格及指数扫描可能耗时较长。T02平面前沿是主要测试，径向谱程序作为补充。T04仅为说明中列明的无坏死有限指数特例。

T01/T03首次运行会在Linux用户缓存目录构建C++程序，随后新输出保存回本仓库的运行目录。支持`--build-only`和`--skip-build`。Windows默认使用WSL的`Ubuntu-24.04`，Python后处理仍在已安装科学依赖的Windows环境执行。可选`T03_onset_probe`在Windows由包装器整体交WSL执行，只需Linux标准库及C++工具链。

## T01现场演示

```text
python A/demo.py --prepare
python A/demo.py --run
python A/demo.py --show-latest
```

`--prepare`提前构建；`--run`使用该构建生成新的半径数据、误差CSV和图；`--show-latest`只显示最近操作记录。标准配置N=64、dt=0.04、T=0.2、32个控制点。它不替代完整加密组。

## 补算程序

```text
python A/supplement/run_t01_supplements.py
python A/supplement/t02_parameter_scan.py
```

第一条生成五组T01独立时间/空间补算并后处理；Windows使用WSL构建和计算，Linux直接计算。第二条生成17组T02参数、波数及分辨率试验。两者均保留新的运行目录；不是只重绘保存结果。原有N256/dt=0.01基线由基线入口另行计算。

T01的N64/dt=0.0025补算存在已披露的节点定位警告，不据此宣称严格空间二阶。完整研究结论和压力本构约定以对应数学技术报告为准。

## 直接使用单项脚本

仓库通过`.gitkeep`保留单项脚本所需的空结果目录；包装器在实际运行时也会创建输出目录。自定义T04短计算可使用：

```text
python A/baseline/code/t04/run_t04.py --mode basic --dx 0.08 --dt 0.002 --gammas 8 --end 0.01
```

自定义短计算仅用于运行检查，不作为完整收敛证据。
