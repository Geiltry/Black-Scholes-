"""Black-Scholes-Merton 欧式期权定价核心模块。"""
import math
import numpy as np

SQRT2 = math.sqrt(2.0)
INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)


def norm_cdf(x: float) -> float:
    """标准正态分布累积分布函数。"""
    return 0.5 * (1.0 + math.erf(x / SQRT2))


def norm_pdf(x: float) -> float:
    """标准正态分布概率密度函数。"""
    return INV_SQRT_2PI * math.exp(-0.5 * x * x)


def _d1_d2(S, K, T, r, sigma):
    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    return d1, d1 - sigma * math.sqrt(T)


def black_scholes(S, K, T, r, sigma, option_type="call") -> float:
    """欧式期权理论价格。

    S: 标的现价, K: 行权价, T: 到期时间(年), r: 无风险利率, sigma: 波动率
    """
    if T <= 0:  # 已到期，返回内在价值
        return max(S - K, 0.0) if option_type == "call" else max(K - S, 0.0)

    d1, d2 = _d1_d2(S, K, T, r, sigma)
    disc_K = K * math.exp(-r * T)
    if option_type == "call":
        return S * norm_cdf(d1) - disc_K * norm_cdf(d2)
    return disc_K * norm_cdf(-d2) - S * norm_cdf(-d1)


def greeks(S, K, T, r, sigma, option_type="call") -> dict:
    """返回 Delta, Gamma, Vega, Theta, Rho。"""
    d1, d2 = _d1_d2(S, K, T, r, sigma)
    disc_K = K * math.exp(-r * T)
    pdf_d1 = norm_pdf(d1)
    sqrt_T = math.sqrt(T)

    gamma = pdf_d1 / (S * sigma * sqrt_T)
    vega = S * pdf_d1 * sqrt_T                      # 对 sigma(1.00) 的敏感度
    theta_common = -S * pdf_d1 * sigma / (2 * sqrt_T)

    if option_type == "call":
        delta = norm_cdf(d1)
        theta = theta_common - r * disc_K * norm_cdf(d2)
        rho = K * T * math.exp(-r * T) * norm_cdf(d2)
    else:
        delta = norm_cdf(d1) - 1.0
        theta = theta_common + r * disc_K * norm_cdf(-d2)
        rho = -K * T * math.exp(-r * T) * norm_cdf(-d2)

    return {"Delta": delta, "Gamma": gamma, "Vega": vega,
            "Theta": theta, "Rho": rho}


def implied_vol(price, S, K, T, r, option_type="call",
                lo=1e-6, hi=5.0, tol=1e-8, max_iter=200) -> float:
    """用二分法由期权价格反解隐含波动率。"""
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if black_scholes(S, K, T, r, mid, option_type) > price:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def monte_carlo_price(S, K, T, r, sigma, option_type="call",
                      n_paths=200_000, seed=42):
    """几何布朗运动下的蒙特卡洛定价，返回 (价格, 标准误)。"""
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal(n_paths)
    ST = S * np.exp((r - 0.5 * sigma ** 2) * T + sigma * math.sqrt(T) * Z)
    payoff = np.maximum(ST - K, 0.0) if option_type == "call" \
        else np.maximum(K - ST, 0.0)
    disc = math.exp(-r * T)
    return disc * payoff.mean(), disc * payoff.std(ddof=1) / math.sqrt(n_paths)