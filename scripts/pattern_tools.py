"""Helper functions for the 38-province example (Section 5 of the paper).

Notation follows the paper:
    M = (I - diag(a) - diag(b) W)^{-1},  d_i = 1/(1-a_i),  r_i = b_i/(1-a_i),
    r(C) = prod_{j in C} r_j,  phi(C) = r(C)/(1-r(C)) * sum_{i in C} d_i,
    h(l) = l rho^l / (1 - rho^l)   (homogeneous case).
"""
import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix

# ---------------------------------------------------------------- matrices

def multiplier(W, a, b):
    """Long-run multiplier M = (I - diag(a) - diag(b) W)^{-1}."""
    N = W.shape[0]
    a = np.broadcast_to(np.asarray(a, float), (N,))
    b = np.broadcast_to(np.asarray(b, float), (N,))
    return np.linalg.inv(np.eye(N) - np.diag(a) - np.diag(b) @ W)


def effects(M):
    """Average direct, total and indirect (spillover) effects (LeSage-Pace)."""
    N = M.shape[0]
    t_dir = np.trace(M) / N
    t_tot = M.sum() / N
    return t_dir, t_tot, t_tot - t_dir


def h(l, rho):
    return l * rho**l / (1 - rho**l)


# ---------------------------------------------------------------- digraphs

def digraph(A):
    G = nx.DiGraph()
    G.add_nodes_from(range(A.shape[0]))
    G.add_edges_from(zip(*np.nonzero(A)))
    return G


def cycles(G, length_bound=None):
    return [tuple(c) for c in nx.simple_cycles(G, length_bound=length_bound)]


# ---------------------------------------------------------------- integer programs

def cycle_incidence(N, cyc):
    lens = np.array([len(c) for c in cyc])
    rows = np.concatenate([np.array(c) for c in cyc])
    cols = np.repeat(np.arange(len(cyc)), lens)
    return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, len(cyc)))


def best_packing(N, cyc, w, sense="max", nonempty=False, slack=None):
    """Optimise sum of w over vertex-disjoint subfamilies of `cyc`.

    slack: optional length-N vector. Each vertex not covered by a chosen cycle
    may add slack[i] to the objective. With slack[i] equal to an upper bound on
    the share of vertex i in any cycle longer than those in `cyc`, the result
    is a rigorous bound on the optimum over all cycles.
    Returns (value, indices of chosen cycles).
    """
    s = -1.0 if sense == "max" else 1.0
    Mx = cycle_incidence(N, cyc)
    m = len(cyc)
    w = np.asarray(w, float)
    if slack is None:
        c, A, n = s * w, Mx, m
    else:
        c = np.r_[s * w, s * np.asarray(slack, float)]
        A = csr_matrix(np.c_[Mx.toarray(), np.eye(N)])
        n = m + N
    cons = [LinearConstraint(A, 0, 1)]
    if nonempty:
        cons.append(LinearConstraint(np.r_[np.ones(m), np.zeros(n - m)][None, :], 1, np.inf))
    res = milp(c, constraints=cons, integrality=np.ones(n), bounds=Bounds(0, 1))
    if res.status != 0:
        raise RuntimeError(res.message)
    return s * res.fun, np.nonzero(res.x[:m] > 0.5)[0]


def hamiltonian_cycle(A, V, time_limit=300):
    """Directed Hamiltonian cycle on the vertex subset V (MTZ formulation).

    Returns the cycle as a list of vertices, or None if the integer program is
    infeasible (proof that no Hamiltonian cycle exists).
    """
    V = list(V)
    B = A[np.ix_(V, V)]
    n = len(V)
    arcs = list(zip(*np.nonzero(B)))
    m = len(arcs)
    nv = m + n
    rows, lb, ub = [], [], []
    for i in range(n):
        for side in (0, 1):  # out-degree 1 and in-degree 1
            r = np.zeros(nv)
            r[[e for e, ab in enumerate(arcs) if ab[side] == i]] = 1
            rows.append(r); lb.append(1); ub.append(1)
    for e, (p, q) in enumerate(arcs):  # subtour elimination
        if p == 0 or q == 0:
            continue
        r = np.zeros(nv)
        r[m + p], r[m + q], r[e] = 1, -1, n
        rows.append(r); lb.append(-np.inf); ub.append(n - 1)
    res = milp(np.zeros(nv), constraints=LinearConstraint(np.array(rows), lb, ub),
               integrality=np.r_[np.ones(m), np.zeros(n)],
               bounds=Bounds(np.r_[np.zeros(m), 0, np.ones(n - 1)],
                             np.r_[np.ones(m), 0, np.full(n - 1, n - 1)]),
               options={"time_limit": time_limit})
    if res.status == 2:  # infeasible
        return None
    if res.status != 0:
        raise RuntimeError(res.message)
    succ = {p: q for (p, q), x in zip(arcs, res.x[:m]) if x > 0.5}
    c = [0]
    while succ[c[-1]] != 0:
        c.append(succ[c[-1]])
    assert len(c) == n and all(B[c[i], c[(i + 1) % n]] for i in range(n))
    return [V[i] for i in c]


# ---------------------------------------------------------------- data

def load_example(data_dir):
    """Rainfall (T x N, mm/month), province table, and year-month labels."""
    rain = pd.read_csv(data_dir / "province_rainfall_monthly.csv", float_precision="round_trip")
    prov = pd.read_csv(data_dir / "provinces.csv", float_precision="round_trip")
    ids = list(rain.columns[1:])
    prov = prov.set_index("region_id").loc[ids].reset_index()
    return rain[ids].to_numpy(), prov, pd.PeriodIndex(rain["year_month"], freq="M")


def monthly_anomalies(Y, months):
    """Standardise each series by calendar month: (x - mean) / sd, sd with ddof=1."""
    df = pd.DataFrame(Y)
    return df.groupby(np.asarray(months.month)).transform(lambda s: (s - s.mean()) / s.std()).to_numpy()


def haversine_km(lat, lon, R=6371.0):
    p, q = np.radians(lat), np.radians(lon)
    s = (np.sin((p[:, None] - p) / 2) ** 2
         + np.cos(p[:, None]) * np.cos(p) * np.sin((q[:, None] - q) / 2) ** 2)
    return 2 * R * np.arcsin(np.sqrt(s))


def knn_pattern(lat, lon, k=4):
    D = haversine_km(lat, lon)
    N = len(lat)
    A = np.zeros((N, N), int)
    for i in range(N):
        o = np.argsort(D[i], kind="stable")
        A[i, o[o != i][:k]] = 1
    return A


def correlation_weights(Z):
    C = np.corrcoef(Z.T)
    np.fill_diagonal(C, 0)
    W = np.clip(C, 0, None)
    return W / W.sum(1, keepdims=True)


def fit_star(Z, W):
    """Pooled OLS for Z(t) = a Z(t-1) + b W Z(t-1) + e(t), no intercept."""
    Y, L = Z[1:], Z[:-1]
    WL = L @ W.T
    a, b = np.linalg.lstsq(np.c_[L.ravel(), WL.ravel()], Y.ravel(), rcond=None)[0]
    return a, b


def fit_gstar(Z, W):
    """Location-wise OLS for the GSTAR(1;1) model, no intercept. Returns (a, b)."""
    Y, L = Z[1:], Z[:-1]
    WL = L @ W.T
    ab = np.array([np.linalg.lstsq(np.c_[L[:, i], WL[:, i]], Y[:, i], rcond=None)[0]
                   for i in range(Z.shape[1])])
    return ab[:, 0], ab[:, 1]
