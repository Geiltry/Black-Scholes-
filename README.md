# Black-Scholes 欧式期权定价器

一个面向金融工程入门学习与基础研究的小型 Python 项目。实现了欧式看涨/看跌期权的 Black-Scholes-Merton 定价模型，并提供希腊字母、隐含波动率反解、蒙特卡洛验证与看跌看涨平价检验。

## 功能特性

- 欧式看涨 / 看跌期权理论定价
- 希腊字母：Delta、Gamma、Vega、Theta、Rho
- 基于二分法的隐含波动率反解
- 几何布朗运动下的蒙特卡洛定价与置信区间
- 看跌看涨平价关系检验

## 项目结构
Black-Scholes/
├── black_scholes.py # 核心定价库
├── demo.py # 演示脚本
├── requirements.txt # 依赖
├── .gitignore
└── README.MD
