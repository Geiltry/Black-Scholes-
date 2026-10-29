"""演示：定价、希腊字母、隐含波动率、蒙特卡洛、平价检验。"""
from black_scholes import black_scholes, greeks, implied_vol, monte_carlo_price
import math

S, K, T, r, sigma = 100.0, 100.0, 1.0, 0.05, 0.20
print(f"参数: S={S}, K={K}, T={T}, r={r}, sigma={sigma}\n")

c = black_scholes(S, K, T, r, sigma, "call")
p = black_scholes(S, K, T, r, sigma, "put")
print(f"看涨期权 BS 价格: {c:.4f}")
print(f"看跌期权 BS 价格: {p:.4f}")

iv = implied_vol(c, S, K, T, r, "call")
print(f"由看涨价格反解隐含波动率: {iv:.6f}\n")

mc, se = monte_carlo_price(S, K, T, r, sigma, "call")
print(f"蒙特卡洛看涨价格: {mc:.4f} ± {1.96*se:.4f} (95% CI)\n")

print("希腊字母（看涨）:")
for k, v in greeks(S, K, T, r, sigma, "call").items():
    print(f"  {k:6s} = {v:10.4f}")

parity = c - p - (S - K * math.exp(-r * T))
print(f"\n看跌看涨平价检验 C-P-(S-Ke^-rT) = {parity:.2e}")