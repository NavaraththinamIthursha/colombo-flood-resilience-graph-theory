"""
Colombo road-network robustness study (Graph Theory research assignment)
Reproduces every number and figure used in the report.
Inputs : colombo_nodes.csv, colombo_roads.csv (OpenStreetMap extract, 1,737 nodes)
Outputs: results/*.csv, figures/*.png, results/key_numbers.json
Run    : python colombo_analysis.py
"""
import json, random, time
import numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from scipy.stats import spearmanr

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 120})
KN = {}                                   # key numbers for the report

# ---------------------------------------------------------------- 1. graph
nodes = pd.read_csv("colombo_nodes.csv"); roads = pd.read_csv("colombo_roads.csv")
pos = {int(r.osmid): (r.x, r.y) for r in nodes.itertuples()}
KN["road_records"] = len(roads); KN["self_loop_records"] = int((roads.u == roads.v).sum())
G0 = nx.from_pandas_edgelist(roads, "u", "v", edge_attr=True, create_using=nx.Graph())   # "Part A" graph
KN["partA_edges"] = G0.number_of_edges()

G = nx.Graph()                            # simple undirected graph, shortest length kept, self-loops dropped
for r in roads.itertuples():
    if r.u == r.v: continue
    L = float(r.length)
    if G.has_edge(r.u, r.v): G[r.u][r.v]["length"] = min(G[r.u][r.v]["length"], L)
    else: G.add_edge(r.u, r.v, length=L, name=r.name)
n, m = G.number_of_nodes(), G.number_of_edges(); P = n*(n-1)/2
KN.update(nodes=n, edges=m, connected=nx.is_connected(G),
          total_length_km=sum(d["length"] for *_, d in G.edges(data=True))/1000,
          mean_degree=2*m/n, density=nx.density(G),
          deg_hist={int(k): int(v) for k, v in pd.Series(dict(G.degree())).value_counts().sort_index().items()},
          cyclomatic=m-n+1)

# ---------------------------------------------------------------- 2. bridges
bridges = list(nx.bridges(G)); B = len(bridges); bridge_set = {frozenset(b) for b in bridges}
imp = {}
for u, v in bridges:
    H = G.copy(); H.remove_edge(u, v)
    s = len(nx.node_connected_component(H, u)); imp[(u, v)] = min(s, n-s)
bi = pd.DataFrame([{"u": u, "v": v, "nodes_cut_off": imp[(u, v)]} for u, v in bridges])
w = np.array([imp[b]*(n-imp[b])/P for b in bridges])            # pair loss of each bridge
KN.update(bridges=B, stubs=int((bi.nodes_cut_off == 1).sum()), b2_4=int(bi.nodes_cut_off.between(2, 4).sum()),
          b5p=int((bi.nodes_cut_off >= 5).sum()), max_cut=int(bi.nodes_cut_off.max()), Phi=float(w.sum()),
          cut_hist={int(k): int(v) for k, v in bi.nodes_cut_off.value_counts().sort_index().items()})
H = G.copy(); H.remove_edges_from(bridges)
blocks = list(nx.connected_components(H)); comp = {x: i for i, c in enumerate(blocks) for x in c}
T = nx.Graph(); T.add_nodes_from(range(len(blocks)))
for u, v in bridges: T.add_edge(comp[u], comp[v], bridge=(u, v))
leaves = sum(1 for x in T if T.degree(x) == 1)
KN.update(blocks=len(blocks), leaves=leaves, lower_bound=int(np.ceil(leaves/2)),
          largest_block=max(len(c) for c in blocks))
bi.to_csv("results/bridge_consequences.csv", index=False)

# ---------------------------------------------------------------- 3. betweenness
ebc_u = nx.edge_betweenness_centrality(G, normalized=True)
ebc_w = nx.edge_betweenness_centrality(G, normalized=True, weight="length")
get = lambda d, e: d[e] if e in d else d[(e[1], e[0])]
edges = list(G.edges())
E = pd.DataFrame({"u": [e[0] for e in edges], "v": [e[1] for e in edges],
                  "name": [G[e[0]][e[1]].get("name") for e in edges],
                  "ebc_hop": [get(ebc_u, e) for e in edges], "ebc_len": [get(ebc_w, e) for e in edges],
                  "is_bridge": [frozenset(e) in bridge_set for e in edges]})
E.to_csv("results/edge_betweenness_both.csv", index=False)
k = int(.05*m)
top_u = set(E.sort_values("ebc_hop", ascending=False).head(k).index)
top_w = set(E.sort_values("ebc_len", ascending=False).head(k).index)
KN.update(spearman_hop_len=float(spearmanr(E.ebc_hop, E.ebc_len)[0]), top5_k=k,
          top5_overlap=len(top_u & top_w)/k, bridges_in_top5=int(E.loc[list(top_u | top_w)].is_bridge.sum()),
          bridge_mean_ebc_hop=float(E[E.is_bridge].ebc_hop.mean()), nonbridge_mean_ebc_hop=float(E[~E.is_bridge].ebc_hop.mean()))
top_names = (E.sort_values("ebc_hop", ascending=False).head(10)[["name", "ebc_hop"]].fillna("(unnamed)"))
top_names_w = (E.sort_values("ebc_len", ascending=False).head(10)[["name", "ebc_len"]].fillna("(unnamed)"))
top_names.to_csv("results/top10_hop.csv", index=False); top_names_w.to_csv("results/top10_len.csv", index=False)

# ---------------------------------------------------------------- 4. delay index (EDI)
rng = random.Random(42); ids = list(G.nodes)
pairs = [tuple(rng.sample(ids, 2)) for _ in range(1500)]
base, usage = {}, {}
for p in pairs:
    path = nx.dijkstra_path(G, p[0], p[1], weight="length")
    base[p] = nx.path_weight(G, path, "length")
    for a, b in zip(path[:-1], path[1:]): usage.setdefault(frozenset((a, b)), []).append(p)
mean_trip = float(np.mean(list(base.values()))); KN["mean_trip_m"] = mean_trip
edi = {}
for e, ps in usage.items():
    a, b = tuple(e); data = G[a][b]; G.remove_edge(a, b); extra = 0
    for p in ps:
        try: extra += nx.dijkstra_path_length(G, p[0], p[1], weight="length") - base[p]
        except nx.NetworkXNoPath: pass
    G.add_edge(a, b, **data); edi[e] = extra/len(pairs)
D = pd.DataFrame([{"u": tuple(e)[0], "v": tuple(e)[1], "edi_pct": 100*val/mean_trip} for e, val in edi.items()])
D["is_bridge"] = [frozenset((r.u, r.v)) in bridge_set for r in D.itertuples()]
D["ebc_hop"] = [get(ebc_u, (r.u, r.v)) for r in D.itertuples()]
D["ebc_len"] = [get(ebc_w, (r.u, r.v)) for r in D.itertuples()]
nbD = D[~D.is_bridge]
KN.update(edi_edges_used=len(D), edi_nonbridge=len(nbD),
          rho_edi_hop=float(spearmanr(nbD.ebc_hop, nbD.edi_pct)[0]), rho_edi_len=float(spearmanr(nbD.ebc_len, nbD.edi_pct)[0]),
          edi_max=float(nbD.edi_pct.max()), edi_median=float(nbD.edi_pct.median()), edi_p95=float(nbD.edi_pct.quantile(.95)))
D.to_csv("results/edge_delay_index.csv", index=False)

# ---------------------------------------------------------------- 5. disruption scenarios (Part A, 1 draw each + 100-draw means)
def metrics(g, name):
    cc = list(nx.connected_components(g)); big = max(len(c) for c in cc)
    return {"scenario": name, "edges": g.number_of_edges(), "components": len(cc),
            "largest_nodes": big, "largest_pct": 100*big/g.number_of_nodes()}
def reach(g, J): return 100*sum(nx.has_path(g, a, b) for a, b in J)/len(J)

random.seed(42)
J = [tuple(random.sample(ids, 2)) for _ in range(500)]
Gb = G.copy(); Gb.remove_edges_from(bridges)
tgt_hop = [(E.u[i], E.v[i]) for i in sorted(top_u)]; tgt_len = [(E.u[i], E.v[i]) for i in sorted(top_w)]
Gh = G.copy(); Gh.remove_edges_from(tgt_hop); Gl = G.copy(); Gl.remove_edges_from(tgt_len)
rows = []
for name, g in [("Baseline", G), ("All 398 bridges closed", Gb), ("Top-5% hop betweenness closed", Gh), ("Top-5% length betweenness closed", Gl)]:
    r = metrics(g, name); r["reach_pct"] = reach(g, J); rows.append(r)
rnd = []
for _ in range(100):
    g = G.copy(); g.remove_edges_from(random.sample(edges, k)); r = metrics(g, "x"); r["reach_pct"] = reach(g, J); rnd.append(r)
rd0 = pd.DataFrame(rnd)
rows.append({"scenario": f"{k} random roads closed (mean of 100)", "edges": m-k, "components": rd0.components.mean(),
             "largest_nodes": rd0.largest_nodes.mean(), "largest_pct": rd0.largest_pct.mean(), "reach_pct": rd0.reach_pct.mean()})
S = pd.DataFrame(rows); S.to_csv("results/disruption_scenarios.csv", index=False)
KN["rand116_comp_sd"] = float(rd0.components.std()); KN["rand116_reach_sd"] = float(rd0.reach_pct.std())

# ---------------------------------------------------------------- 6. augmentation ILP
lon = np.array([pos[i][0] for i in ids]); lat = np.array([pos[i][1] for i in ids])
X = (lon-lon.mean())*111320*np.cos(np.radians(6.925)); Y = (lat-lat.mean())*110574
tree = cKDTree(np.c_[X, Y]); R = 500; cand = []
for i, j in tree.query_pairs(R):
    a, b = ids[i], ids[j]
    if comp[a] == comp[b] or G.has_edge(a, b): continue
    d = float(np.hypot(X[i]-X[j], Y[i]-Y[j]))
    path = nx.shortest_path(T, comp[a], comp[b])
    cand.append((a, b, 1.2*d, [T[path[q]][path[q+1]]["bridge"] for q in range(len(path)-1)]))
Cn = len(cand); KN["candidates"] = Cn
bidx = {frozenset(b): i for i, b in enumerate(bridges)}
A = lil_matrix((B, Cn))
for j, (a, b, c, cv) in enumerate(cand):
    for e in cv: A[bidx[frozenset(e)], j] = 1
A = csr_matrix(A); cost = np.array([c[2] for c in cand])
def solve(obj, tl=120):
    t = time.time(); r = milp(obj, constraints=LinearConstraint(A, lb=np.ones(B)), integrality=np.ones(Cn), bounds=Bounds(0, 1), options={"time_limit": tl})
    return r, time.time()-t
r1, t1 = solve(np.ones(Cn)); sel_count = [j for j in range(Cn) if r1.x[j] > .5]
r2, t2 = solve(cost);        sel_length = [j for j in range(Cn) if r2.x[j] > .5]
KN.update(ilp_count_links=len(sel_count), ilp_count_km=float(cost[sel_count].sum()/1000), ilp_count_status=str(r1.message), ilp_count_time=t1,
          ilp_len_links=len(sel_length), ilp_len_km=float(cost[sel_length].sum()/1000), ilp_len_status=str(r2.message), ilp_len_time=t2,
          cand_mean_len=float(cost.mean()), ilp_count_mean_len=float(cost[sel_count].mean()), ilp_len_mean_len=float(cost[sel_length].mean()))

K = 200; Ad = A.tocsc(); covered = np.zeros(B, bool); hist = []
for q in range(K):
    gain = Ad.T.dot((~covered)*(w+1e-9)); j = int(np.argmax(gain))
    if gain[j] <= 0: break
    covered |= Ad[:, j].toarray().ravel() > 0; hist.append((q+1, int(covered.sum()), float(w[~covered].sum())))
hG = pd.DataFrame(hist, columns=["links", "covered", "loss_left"]); hG.to_csv("results/greedy_budget_curve.csv", index=False)
rg = np.random.default_rng(42); runs = []
for rep in range(30):
    order = rg.permutation(Cn); cov = np.zeros(B, bool); h = []
    for q in range(K):
        cov |= A[:, order[q]].toarray().ravel() > 0; h.append((cov.sum(), w[~cov].sum()))
    runs.append(h)
hR = np.array(runs).mean(0); w0 = w.sum()
for f in (.5, .8, .9, .95):
    KN[f"greedy_{int(f*100)}"] = int(hG[hG.loss_left <= (1-f)*w0].iloc[0].links)
for kk in (50, 100, 150, 179):
    KN[f"budget{kk}_greedy_bridges_left"] = int(B-hG.covered[kk-1]); KN[f"budget{kk}_random_bridges_left"] = float(B-hR[kk-1, 0])
    KN[f"budget{kk}_greedy_loss_left_pct"] = float(100*hG.loss_left[kk-1]/w0); KN[f"budget{kk}_random_loss_left_pct"] = float(100*hR[kk-1, 1]/w0)
# does a greedy set of 179 links also remove all bridges?  (no - report honestly)

# ---------------------------------------------------------------- 7. robustness of improved networks
def add_links(sel):
    g = G.copy()
    for j in sel: a, b, c, _ = cand[j]; g.add_edge(a, b, length=c)
    return g
G_opt = add_links(sel_count); G_len = add_links(sel_length)
rr = np.random.default_rng(7); G_rand = add_links(list(rr.choice(Cn, size=len(sel_count), replace=False)))
KN.update(bridges_left_opt=len(list(nx.bridges(G_opt))), bridges_left_len=len(list(nx.bridges(G_len))), bridges_left_rand=len(list(nx.bridges(G_rand))))
random.seed(42); rows = []
for name, g in [("Baseline", G), ("Random links (179)", G_rand), ("ILP min-count (179 links)", G_opt), ("ILP min-length (220 links)", G_len)]:
    Eg = list(g.edges())
    for frac in (0.05, 0.10, 0.15):
        for rep in range(100):
            h = g.copy(); h.remove_edges_from(random.sample(Eg, int(len(Eg)*frac))); r = metrics(h, name)
            rows.append((name, frac, r["components"], r["largest_pct"]))
rd = pd.DataFrame(rows, columns=["scenario", "frac", "components", "lcc"])
sm = rd.groupby(["scenario", "frac"]).agg(comp=("components", "mean"), comp_sd=("components", "std"), lcc=("lcc", "mean"), lcc_sd=("lcc", "std")).reset_index()
sm.to_csv("results/robustness_with_links.csv", index=False)
# targeted attack: remove baseline's top-5% hop-betweenness roads from each network
att = []
for name, g in [("Baseline", G), ("Random links (179)", G_rand), ("ILP min-count (179 links)", G_opt), ("ILP min-length (220 links)", G_len)]:
    h = g.copy(); h.remove_edges_from(tgt_hop); r = metrics(h, name); r["reach_pct"] = reach(h, J); att.append(r)
pd.DataFrame(att).to_csv("results/targeted_attack_with_links.csv", index=False)
# mean trip length change from added links (detour/benefit)
sub = pairs[:300]
def mean_trip_len(g): return float(np.mean([nx.dijkstra_path_length(g, a, b, weight="length") for a, b in sub]))
KN.update(trip_base=mean_trip_len(G), trip_opt=mean_trip_len(G_opt), trip_len=mean_trip_len(G_len))

pd.DataFrame(hist).to_csv("results/_hist.csv", index=False)
json.dump(KN, open("results/key_numbers.json", "w"), indent=1, default=str)

# ---------------------------------------------------------------- figures
def draw_base(ax): 
    for u, v in G.edges(): ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#cfcfcf", lw=.5, zorder=1)
def finish(ax, t): ax.set_aspect("equal"); ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_title(t)

fig, ax = plt.subplots(figsize=(8, 7.5)); draw_base(ax)
for r in bi.itertuples():
    s = r.nodes_cut_off; c = "#f4a261" if s == 1 else ("#e63946" if s <= 4 else "#6a040f")
    ax.plot([pos[r.u][0], pos[r.v][0]], [pos[r.u][1], pos[r.v][1]], color=c, lw=1.2 if s == 1 else 2.4, zorder=2)
ax.legend(handles=[Line2D([0], [0], color="#cfcfcf", lw=1, label=f"Non-bridge road ({m-B:,})"),
                   Line2D([0], [0], color="#f4a261", lw=1.5, label=f"Bridge, cuts off 1 node ({KN['stubs']})"),
                   Line2D([0], [0], color="#e63946", lw=2.4, label=f"Bridge, cuts off 2-4 nodes ({KN['b2_4']})"),
                   Line2D([0], [0], color="#6a040f", lw=2.4, label=f"Bridge, cuts off 5-10 nodes ({KN['b5p']})")], loc="lower left", fontsize=8)
finish(ax, "Colombo road network: bridges coloured by consequence"); plt.tight_layout(); plt.savefig("figures/fig1_bridge_map.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8)); vc = bi.nodes_cut_off.value_counts().sort_index()
ax[0].bar(vc.index.astype(str), vc.values, color="#e76f51"); ax[0].set_yscale("log")
ax[0].set_xlabel("Nodes cut off by closing the bridge"); ax[0].set_ylabel("Number of bridges (log)"); ax[0].set_title("(a) Consequence of bridge closure")
for x, y in zip(vc.index.astype(str), vc.values): ax[0].text(x, y*1.1, str(y), ha="center", fontsize=8)
nb_ = E[~E.is_bridge]; br_ = E[E.is_bridge]
ax[1].scatter(nb_.ebc_hop, nb_.ebc_len, s=5, alpha=.5, c="#264653", label="non-bridge"); ax[1].scatter(br_.ebc_hop, br_.ebc_len, s=10, c="#e63946", label="bridge")
ax[1].set_xlabel("Hop-count edge betweenness"); ax[1].set_ylabel("Length-weighted edge betweenness"); ax[1].set_title("(b) Hop vs length-weighted betweenness"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("figures/fig2_impact_betweenness.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for a, col, t in [(ax[0], "ebc_hop", "Hop-count betweenness"), (ax[1], "ebc_len", "Length-weighted betweenness")]:
    a.scatter(nbD[col], nbD.edi_pct, s=5, alpha=.5, c="#2a9d8f"); a.set_xlabel(t); a.set_ylabel("Edge Delay Index (% of mean trip)")
    a.set_title(f"Spearman rho = {spearmanr(nbD[col], nbD.edi_pct)[0]:.2f}")
plt.tight_layout(); plt.savefig("figures/fig3_edi.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 3, figsize=(12, 3.8)); lab = ["Baseline", "All bridges", "Top-5%\nhop EBC", "Top-5%\nlength EBC", "116 random\n(mean)"]
cols = ["#6c757d", "#6a040f", "#e63946", "#f4a261", "#2a9d8f"]
ax[0].bar(lab, S.components, color=cols); ax[0].set_ylabel("Connected components"); ax[0].set_title("(a) Fragmentation")
ax[1].bar(lab, S.reach_pct, color=cols); ax[1].set_ylabel("Reachable journeys (%)"); ax[1].set_ylim(0, 105); ax[1].set_title("(b) Journey reachability (500 O-D pairs)")
rs = pd.read_csv("results/robustness_with_links.csv"); d = rs[rs.scenario == "Baseline"]
ax[2].errorbar(d.frac*100, d.lcc, yerr=d.lcc_sd, marker="o", capsize=3, c="#264653"); ax[2].set_xlabel("Randomly closed roads (%)"); ax[2].set_ylabel("Largest component (% of nodes)"); ax[2].set_title("(c) Random closure, 100 runs")
for a in ax[:2]: a.tick_params(axis="x", labelsize=7)
plt.tight_layout(); plt.savefig("figures/fig4_scenarios.png", dpi=150); plt.close()

fig, ax = plt.subplots(figsize=(8, 7.5)); draw_base(ax)
for j in sel_length:
    a, b, _, _ = cand[j]; ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#1d7874", lw=1.8, zorder=3)
ax.legend(handles=[Line2D([0], [0], color="#cfcfcf", lw=1, label="Existing roads"), Line2D([0], [0], color="#1d7874", lw=2, label=f"Proposed links ({len(sel_length)}, min total length)")], loc="lower left", fontsize=8)
finish(ax, "Minimum-length 2-edge-connectivity augmentation"); plt.tight_layout(); plt.savefig("figures/fig5_augmentation_map.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(hG.links, B-hG.covered, c="#264653", label="Greedy (consequence-weighted)"); ax[0].plot(np.arange(1, K+1), B-hR[:, 0], c="#e76f51", ls="--", label="Random candidate links (mean of 30)")
ax[0].axvline(KN["lower_bound"], c="k", ls=":", lw=1); ax[0].set_xlabel("New links added"); ax[0].set_ylabel("Bridges remaining"); ax[0].set_title("(a) Bridges remaining"); ax[0].legend(fontsize=8)
ax[1].plot(hG.links, 100*hG.loss_left/w0, c="#264653"); ax[1].plot(np.arange(1, K+1), 100*hR[:, 1]/w0, c="#e76f51", ls="--")
ax[1].set_xlabel("New links added"); ax[1].set_ylabel("Residual single-bridge vulnerability (%)"); ax[1].set_title("(b) Consequence-weighted vulnerability")
plt.tight_layout(); plt.savefig("figures/fig6_budget_curve.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for name, c in [("Baseline", "#6c757d"), ("Random links (179)", "#e76f51"), ("ILP min-count (179 links)", "#1d7874")]:
    d = sm[sm.scenario == name]
    ax[0].errorbar(d.frac*100, d.lcc, yerr=d.lcc_sd, marker="o", capsize=3, label=name, c=c); ax[1].errorbar(d.frac*100, d.comp, yerr=d.comp_sd, marker="o", capsize=3, label=name, c=c)
ax[0].set_xlabel("Randomly closed roads (%)"); ax[0].set_ylabel("Largest component (% of nodes)"); ax[0].set_title("(a) Largest component"); ax[0].legend(fontsize=8)
ax[1].set_xlabel("Randomly closed roads (%)"); ax[1].set_ylabel("Connected components"); ax[1].set_title("(b) Fragmentation")
plt.tight_layout(); plt.savefig("figures/fig7_robustness.png", dpi=150); plt.close()
print("done")
