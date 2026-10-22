"""Section 5 of the paper: monthly rainfall in 38 Indonesian provinces.

Reproduces every number in Section 5 (Table 1 and the text) from the files in
data/. Writes results/example_provinces.json. Run:
    python scripts/example_provinces.py              (about one minute)
    python scripts/example_provinces.py --all-cycles (also enumerates all simple
                                                      cycles of D_geo)
"""
import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np

from pattern_tools import (best_packing, correlation_weights, cycles, digraph, effects, fit_gstar,
                           fit_star, h, hamiltonian_cycle, knn_pattern, load_example,
                           monthly_anomalies, multiplier)

ROOT = Path(__file__).resolve().parents[1]
EAST = ["Maluku", "North Maluku", "Papua", "West Papua", "Southwest Papua",
        "Highland Papua", "South Papua", "Central Papua"]


def mutual_pairs_matching(A):
    """Largest number of vertex-disjoint cycles of length 2."""
    N = A.shape[0]
    U = nx.Graph()
    U.add_edges_from((i, j) for i in range(N) for j in range(i + 1, N) if A[i, j] and A[j, i])
    return len(nx.max_weight_matching(U, maxcardinality=True))


def structure_geo(A, names, all_cycles):
    N = len(names)
    G = digraph(A)
    east = [names.index(x) for x in EAST]
    west = [i for i in range(N) if i not in east]
    mu = names.index("North Maluku")
    out = {
        "strongly_connected": nx.is_strongly_connected(G),
        "arcs_west_to_east": [(names[i], names[j]) for i in west for j in east if A[i, j]],
        "arcs_east_to_west": [(names[i], names[j]) for i in east for j in west if A[i, j]],
        "longest_cycle_inside_east": max(len(c) for c in cycles(G.subgraph(east))),
    }
    # A cycle that meets both groups can use only North Maluku from the east,
    # so the circumference is decided by the two integer programs below.
    out["west_plus_north_maluku_hamiltonian"] = hamiltonian_cycle(A, west + [mu]) is not None
    c30 = hamiltonian_cycle(A, west)
    out["west_hamiltonian"] = c30 is not None
    out["circumference"] = len(c30) if (c30 and not out["west_plus_north_maluku_hamiltonian"]) else None
    out["cycle_through_30_west"] = [names[i] for i in c30] if c30 else None
    out["disjoint_2_cycles"] = mutual_pairs_matching(A)
    if all_cycles:
        n, longest = 0, 0
        for c in nx.simple_cycles(G):
            n, longest = n + 1, max(longest, len(c))
        out["all_simple_cycles"] = {"count": n, "longest": longest}
    return out


def structure_cor(W):
    A = (W > 0).astype(int)
    N = A.shape[0]
    return {
        "symmetric": bool((A == A.T).all()),
        "mean_out_degree": float(A.sum(1).mean()),
        "hamiltonian": hamiltonian_cycle(A, range(N)) is not None,
        "disjoint_2_cycles": mutual_pairs_matching(A),
    }


def homogeneous(Z, W):
    N = W.shape[0]
    a, b = fit_star(Z, W)
    d, rho = 1 / (1 - a), b / (1 - a)
    t_dir, t_tot, t_ind = effects(multiplier(W, a, b))
    univ_min = d / (1 - rho**N)          # Theorem 1.3, rho > 0
    univ_max = d * (1 + (2 * (N // 2) / N) * rho**2 / (1 - rho**2))
    return dict(a=a, b=b, rho=rho, d=d, tau_dir=t_dir, tau_tot=t_tot, tau_ind=t_ind,
                univ_min_tau_dir=univ_min, univ_max_tau_dir=univ_max,
                univ_tau_ind=[t_tot - univ_max, t_tot - univ_min], one_plus_rho=1 + rho)


def pattern_max_geo(A, rho, L=12):
    """Max of tr M / d over Sigma_D (Corollary 1.2) by set packing over cycles of
    length <= L. A vertex on a longer cycle adds at most rho^(L+1)/(1-rho^(L+1))."""
    N = A.shape[0]
    cyc = cycles(digraph(A), length_bound=L)
    w = [h(len(c), rho) for c in cyc]
    lower, chosen = best_packing(N, cyc, w, "max")
    tail = rho ** (L + 1) / (1 - rho ** (L + 1))
    upper, _ = best_packing(N, cyc, w, "max", slack=np.full(N, tail))
    return dict(n_cycles=len(cyc), value=N + lower, bound=N + upper, tail_per_vertex=tail,
                chosen_lengths=sorted(len(cyc[k]) for k in chosen))


def heterogeneous_geo(Z, W, A, L=16):
    """Theorem 1.1 with location-wise GSTAR estimates and W in Sigma_{D_geo}."""
    N = A.shape[0]
    a, b = fit_gstar(Z, W)
    d, r = 1 / (1 - a), b / (1 - a)
    cyc = cycles(digraph(A), length_bound=L)
    phi = [np.prod(r[list(c)]) / (1 - np.prod(r[list(c)])) * d[list(c)].sum() for c in cyc]
    rbar = np.abs(r).max()
    g = rbar ** (L + 1)
    tail = d * g / (1 - g)   # |share of vertex i| in phi(C) for any cycle longer than L
    out = dict(max_abs_a_plus_abs_b=float((np.abs(a) + np.abs(b)).max()),
               r_range=[float(r.min()), float(r.max())], n_negative_r=int((r < 0).sum()),
               tau_dir=effects(multiplier(W, a, b))[0], mean_d=float(d.mean()),
               n_cycles=len(cyc), tail_max=float(tail.max()))
    # D_geo is strongly connected, so every nonempty family is admissible.
    for sense, sign in (("min", -1), ("max", 1)):
        v, _ = best_packing(N, cyc, phi, sense, nonempty=True)
        bd, _ = best_packing(N, cyc, phi, sense, slack=sign * tail)
        out[f"{sense}_tau_dir"] = (d.sum() + v) / N        # attained by a 0-1 matrix
        out[f"{sense}_tau_dir_bound"] = (d.sum() + bd) / N  # rigorous outer bound
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-cycles", action="store_true")
    args = ap.parse_args()

    Y, prov, months = load_example(ROOT / "data")
    names = list(prov["name"])
    N = len(names)
    Z = monthly_anomalies(Y, months)
    A_geo = knn_pattern(prov["lat"].to_numpy(), prov["lon"].to_numpy(), k=4)
    W_geo = A_geo / A_geo.sum(1, keepdims=True)
    W_cor = correlation_weights(Z)

    res = {"data": {"N": N, "T": len(months), "first": str(months[0]), "last": str(months[-1])}}
    res["structure_geo"] = structure_geo(A_geo, names, args.all_cycles)
    res["structure_cor"] = structure_cor(W_cor)

    geo, cor = homogeneous(Z, W_geo), homogeneous(Z, W_cor)
    k = res["structure_geo"]["circumference"]
    geo["pattern_min_tau_dir"] = geo["d"] * (N - k + k / (1 - geo["rho"] ** k)) / N
    pm = pattern_max_geo(A_geo, geo["rho"])
    geo["pattern_max_tau_dir"] = [geo["d"] * pm["value"] / N, geo["d"] * pm["bound"] / N]
    geo["pattern_max_detail"] = pm
    geo["pattern_tau_ind"] = [geo["tau_tot"] - geo["pattern_max_tau_dir"][1],
                              geo["tau_tot"] - geo["pattern_min_tau_dir"]]
    geo["rho_power_30"] = geo["rho"] ** 30
    # D_cor is Hamiltonian and has 19 = N/2 disjoint 2-cycles: pattern bounds = universal bounds.
    cor["pattern_tau_ind"] = list(cor["univ_tau_ind"])
    res["homogeneous"] = {"geo": geo, "cor": cor}
    res["ratios"] = {
        "tau_ind_cor_over_geo": cor["tau_ind"] / geo["tau_ind"],
        "geo_universal": geo["univ_tau_ind"][1] / geo["univ_tau_ind"][0],
        "geo_pattern": geo["pattern_tau_ind"][1] / geo["pattern_tau_ind"][0],
        "cor_universal": cor["univ_tau_ind"][1] / cor["univ_tau_ind"][0],
    }
    res["heterogeneous_geo"] = heterogeneous_geo(Z, W_geo, A_geo)
    a_k, b_k = fit_gstar(Z, W_cor)
    res["heterogeneous_cor"] = {"max_abs_a_plus_abs_b": float((np.abs(a_k) + np.abs(b_k)).max())}

    (ROOT / "results").mkdir(exist_ok=True)
    with open(ROOT / "results" / "example_provinces.json", "w") as f:
        json.dump(res, f, indent=2, default=float)

    # ---------------------------------------------------------------- summary
    sg, sk = res["structure_geo"], res["structure_cor"]
    print(f"Data: N = {N}, T = {len(months)} months ({months[0]} to {months[-1]})")
    print(f"D_geo: strongly connected = {sg['strongly_connected']}, circumference = {sg['circumference']}, "
          f"west+North Maluku Hamiltonian = {sg['west_plus_north_maluku_hamiltonian']}, "
          f"disjoint 2-cycles = {sg['disjoint_2_cycles']}")
    print(f"       arcs west->east {sg['arcs_west_to_east']}; east->west {sg['arcs_east_to_west']}")
    if "all_simple_cycles" in sg:
        print(f"       all simple cycles: {sg['all_simple_cycles']}")
    print(f"D_cor: symmetric = {sk['symmetric']}, Hamiltonian = {sk['hamiltonian']}, "
          f"disjoint 2-cycles = {sk['disjoint_2_cycles']}, mean out-degree = {sk['mean_out_degree']:.1f}")
    print("\nTable 1")
    print(f"{'W':5s} {'a':>7s} {'b':>7s} {'rho':>7s} {'t_dir':>7s} {'t_ind':>7s} "
          f"{'pattern range':>18s} {'all-W range':>18s} {'1+rho':>7s}")
    for nm, R in (("geo", geo), ("cor", cor)):
        p, u = R["pattern_tau_ind"], R["univ_tau_ind"]
        print(f"{nm:5s} {R['a']:7.4f} {R['b']:7.4f} {R['rho']:7.4f} {R['tau_dir']:7.3f} {R['tau_ind']:7.3f} "
              f"  [{p[0]:.4f}, {p[1]:.4f}]   [{u[0]:.4f}, {u[1]:.4f}] {R['one_plus_rho']:7.4f}")
    r = res["ratios"]
    print(f"\nEstimated spillover ratio cor/geo = {r['tau_ind_cor_over_geo']:.3f}")
    print(f"Range ratios with fixed coefficients: geo all W {r['geo_universal']:.4f}, "
          f"geo pattern {r['geo_pattern']:.4f}, cor {r['cor_universal']:.4f}")
    print(f"geo: rho^30 = {geo['rho_power_30']:.1e}; pattern min t_dir {geo['pattern_min_tau_dir']:.6f} "
          f"vs all-W min {geo['univ_min_tau_dir']:.6f}")
    print(f"geo: pattern max t_dir in [{geo['pattern_max_tau_dir'][0]:.8f}, {geo['pattern_max_tau_dir'][1]:.8f}] "
          f"({pm['n_cycles']} cycles of length <= 12) vs all-W max {geo['univ_max_tau_dir']:.4f}")
    hg = res["heterogeneous_geo"]
    print(f"\nGSTAR, W_geo: max|a_i|+|b_i| = {hg['max_abs_a_plus_abs_b']:.3f} "
          f"(W_cor: {res['heterogeneous_cor']['max_abs_a_plus_abs_b']:.3f}); "
          f"r_i in [{hg['r_range'][0]:.3f}, {hg['r_range'][1]:.3f}], {hg['n_negative_r']} negative")
    print(f"   estimated t_dir = {hg['tau_dir']:.4f}; mean d_i = {hg['mean_d']:.4f}")
    print(f"   Theorem 1.1 over Sigma_D_geo ({hg['n_cycles']} cycles of length <= 16):")
    print(f"   min t_dir in [{hg['min_tau_dir_bound']:.5f}, {hg['min_tau_dir']:.5f}], "
          f"max t_dir in [{hg['max_tau_dir']:.5f}, {hg['max_tau_dir_bound']:.5f}]")


if __name__ == "__main__":
    main()
