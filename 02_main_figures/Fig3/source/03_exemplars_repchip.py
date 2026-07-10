"""Pick exemplar programs by RCTD-weight correlation (OLIGO, ENDO) and layer specificity,
and select representative chip (good layer coverage + size). Also write rep-chip per-bin
tables for panels a,b,c. Subsample/grid as needed."""
import pyarrow.parquet as pq
import numpy as np
import pandas as pd
from collections import defaultdict

base = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"
work = "CORTEX_PROGRAM_ROOT/scripts/fig2/"
PROGS = [f"program_{i}" for i in range(1,61)]
CELLTYPES = ["AST","CHANDELIER","ENDO","ET","L2-L3 IT LINC00507","L3-L4 IT RORB",
 "L4-L5 IT RORB","L6 CAR3","L6 CT","L6 IT","L6B","LAMP5","MICRO","NDNF","NP",
 "OLIGO","OPC","PAX6","PVALB","SST","VIP","VLMC"]

meta = pd.read_parquet(work+"_meta_cache.parquet", columns=["bin","chip","x","y","region","majorDomain"]).set_index("bin")

# ---- correlation program vs celltype weight (streamed, sufficient-stats) ----
# we compute Pearson via running sums on pass bins only
pf_p = pq.ParquetFile(base+"spatial_bin50_program_score.parquet")
pf_r = pq.ParquetFile(base+"spatial_bin50_rctd_weights.parquet")
it_p = pf_p.iter_batches(batch_size=400000, columns=["bin"]+PROGS)
it_r = pf_r.iter_batches(batch_size=400000, columns=["bin"]+CELLTYPES+["rctd_pass_mask"])

n=0.0
sx=np.zeros(60); sy=np.zeros(22)
sxx=np.zeros(60); syy=np.zeros(22)
sxy=np.zeros((60,22))
for bp, br in zip(it_p, it_r):
    dp=bp.to_pandas().set_index("bin"); dr=br.to_pandas().set_index("bin")
    mask = dr["rctd_pass_mask"].to_numpy().astype(bool)
    X=dp[PROGS].to_numpy(dtype=np.float64)[mask]
    Y=dr[CELLTYPES].to_numpy(dtype=np.float64)[mask]
    n+=X.shape[0]
    sx+=X.sum(0); sy+=Y.sum(0)
    sxx+=(X*X).sum(0); syy+=(Y*Y).sum(0)
    sxy+=X.T@Y
mx=sx/n; my=sy/n
cov=sxy/n - np.outer(mx,my)
sdx=np.sqrt(sxx/n - mx**2); sdy=np.sqrt(syy/n - my**2)
corr=cov/np.outer(sdx,sdy)
C=pd.DataFrame(corr, index=PROGS, columns=CELLTYPES)
C.to_csv(work+"program_celltype_corr.tsv", sep="\t")
print("rctd pass bins:", int(n))
print("\nTop programs for OLIGO:\n", C["OLIGO"].sort_values(ascending=False).head(5).to_string())
print("\nTop programs for ENDO:\n", C["ENDO"].sort_values(ascending=False).head(5).to_string())
print("\nTop programs for VLMC:\n", C["VLMC"].sort_values(ascending=False).head(5).to_string())
print("\nProgram 37/26/56 corr to OLIGO/ENDO:\n", C.loc[["program_37","program_26","program_56"],["OLIGO","ENDO","VLMC","AST"]].to_string())

# layer-specific picks from earlier file
spec = pd.read_csv(work+"program_layer_specificity.tsv", sep="\t")
def top_for_layer(ly, k=1):
    return spec[spec.peak_layer==ly].sort_values("contrast",ascending=False).head(k)["program"].tolist()
print("\nL2:",top_for_layer("L2",3)," L4:",top_for_layer("L4",3)," L6:",top_for_layer("L6",3))

# inhibitory program: best corr to SST/PVALB/VIP/LAMP5
inh = C[["SST","PVALB","VIP","LAMP5","NDNF"]].max(1).sort_values(ascending=False)
print("\nTop inhibitory-assoc programs:\n", inh.head(5).to_string())
