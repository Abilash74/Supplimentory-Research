import networkx as nx
import pandas as pd
import numpy as np
import os
import random
import stellargraph as sg
from stellargraph.data import UnsupervisedSampler
from stellargraph.mapper import Attri2VecLinkGenerator, Attri2VecNodeGenerator
from stellargraph.layer import Attri2Vec, link_classification
from tensorflow import keras
from sklearn.metrics import roc_auc_score
import xgboost as xgb

# Set seed for reproducibility
seed = 42
np.random.seed(seed)
random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

data_dir = "data/"  # working directory
walk_length = 10
number_of_walks = 3
batch_size = 32
epochs = 100
layer_sizes = [16]
kfold = 5

# Read edge list
edgelist = pd.read_csv(os.path.join(data_dir, "pair_Sa"),
                       sep="\t", header=None, names=["source", "target", "label", "eo1", "eo2"])

# Read node list
file_nodes = pd.read_csv(os.path.join(data_dir, "pair_content_Sa"), sep="\t")
nodes = file_nodes.set_index('ID')
feats = nodes.columns[:-1]  # Exclude label column

# Create graph
G_all_nx = nx.from_pandas_edgelist(edgelist, edge_attr="label")
nx.set_node_attributes(G_all_nx, "label", "label")
G_all = sg.StellarGraph.from_networkx(G_all_nx, node_features=nodes[feats])

# Store probabilities for each method
synergy_probs = {"Hadamard": [], "L1-norm": [], "L2-norm": [], "Average": [], "Classification-based": []}
antagonistic_probs = {"Hadamard": [], "L1-norm": [], "L2-norm": [], "Average": [], "Classification-based": []}
lbl_p_synergistic, lbl_p_antagonistic = [], []

# Train Attri2Vec model
generator = Attri2VecNodeGenerator(G_all, batch_size)
attri2vec = Attri2Vec(layer_sizes=layer_sizes, generator=generator, bias=False, normalize=None)
x_inp, x_out = attri2vec.in_out_tensors()
embedding_model = keras.Model(inputs=x_inp, outputs=x_out)

for k in range(kfold):
    print(f"Processing Fold {k+1}/{kfold}...")

    # Read train and test interaction data
    edge_t = pd.read_csv(os.path.join(data_dir, f"{k+1}_t"), sep="\t",
                         header=None, names=["source", "target", "label", "eo1", "eo2"])
    edge_p = pd.read_csv(os.path.join(data_dir, f"{k+1}_p"), sep="\t",
                         header=None, names=["source", "target", "label", "eo1", "eo2"])

    # Convert labels from [1,2,3] → [0,1,2] for XGBoost compatibility
    edge_t["label"] = edge_t["label"] - 1
    edge_p["label"] = edge_p["label"] - 1

    # Generate embeddings
    node_ids = nodes.index
    node_gen = generator.flow(node_ids)
    emb = embedding_model.predict(node_gen, workers=4, verbose=1)

    # Function to extract edge features
    def extract_features(edge_data):
        feat_h, feat_L1, feat_L2, feat_av, feat_cmp_av = [], [], [], [], []
        for i in range(len(edge_data)):
            try:
                n1 = np.where(nodes.index == edge_data["source"][i])[0][0]
                n2 = np.where(nodes.index == edge_data["target"][i])[0][0]

                feat_h.append(np.ravel(emb[n1] * emb[n2]))  # Hadamard
                feat_L1.append(np.ravel(np.abs(emb[n1] - emb[n2])))  # L1-norm
                feat_L2.append(np.ravel((emb[n1] - emb[n2]) ** 2))  # L2-norm
                feat_av.append(np.ravel((emb[n1] + emb[n2]) / 2))  # Average
                feat_cmp_av.append((nodes.loc[edge_data["source"][i], feats] +
                                    nodes.loc[edge_data["target"][i], feats]) / 2)  # Classification-based
            except IndexError:
                print(f"Skipping invalid edge: {edge_data.iloc[i].to_dict()}")
        
        return np.array(feat_h), np.array(feat_L1), np.array(feat_L2), np.array(feat_av), np.array(feat_cmp_av)

    # Extract features for training and testing
    h_feat_t, L1_feat_t, L2_feat_t, av_feat_t, cmp_av_feat_t = extract_features(edge_t)
    h_feat_p, L1_feat_p, L2_feat_p, av_feat_p, cmp_av_feat_p = extract_features(edge_p)

    # Prepare labels
    lbl_p = edge_p["label"].values

    # **Separate Labels for Two Models**
    lbl_p_synergistic.extend((lbl_p == 1).astype(int))  # 1 for synergy, 0 for others
    lbl_p_antagonistic.extend((lbl_p == 2).astype(int))  # 1 for antagonistic, 0 for others

    # Function to train XGBoost for Binary Classification
    def train_xgb(train_X, train_y, test_X):
        xgb_model = xgb.XGBClassifier(use_label_encoder=False, eval_metric="logloss", random_state=seed)
        xgb_model.fit(train_X, train_y)
        pred_proba = xgb_model.predict_proba(test_X)[:, 1]  # Only positive class probability
        return pred_proba

    # **Train Separate Models for Each Feature Extraction Method**
    synergy_probs["Hadamard"].append(train_xgb(h_feat_t, (edge_t["label"] == 1).astype(int), h_feat_p))
    synergy_probs["L1-norm"].append(train_xgb(L1_feat_t, (edge_t["label"] == 1).astype(int), L1_feat_p))
    synergy_probs["L2-norm"].append(train_xgb(L2_feat_t, (edge_t["label"] == 1).astype(int), L2_feat_p))
    synergy_probs["Average"].append(train_xgb(av_feat_t, (edge_t["label"] == 1).astype(int), av_feat_p))
    synergy_probs["Classification-based"].append(train_xgb(cmp_av_feat_t, (edge_t["label"] == 1).astype(int), cmp_av_feat_p))

    antagonistic_probs["Hadamard"].append(train_xgb(h_feat_t, (edge_t["label"] == 2).astype(int), h_feat_p))
    antagonistic_probs["L1-norm"].append(train_xgb(L1_feat_t, (edge_t["label"] == 2).astype(int), L1_feat_p))
    antagonistic_probs["L2-norm"].append(train_xgb(L2_feat_t, (edge_t["label"] == 2).astype(int), L2_feat_p))
    antagonistic_probs["Average"].append(train_xgb(av_feat_t, (edge_t["label"] == 2).astype(int), av_feat_p))
    antagonistic_probs["Classification-based"].append(train_xgb(cmp_av_feat_t, (edge_t["label"] == 2).astype(int), cmp_av_feat_p))

# Compute AUC scores separately
print("\nAUC (Synergistic vs. Rest):")
for method in synergy_probs:
    print(f"{method}: {roc_auc_score(lbl_p_synergistic, np.hstack(synergy_probs[method])):.4f}")

print("\nAUC (Antagonistic vs. Rest):")
for method in antagonistic_probs:
    print(f"{method}: {roc_auc_score(lbl_p_antagonistic, np.hstack(antagonistic_probs[method])):.4f}")

# ============================================================
# SINGLE ROC PLOT (XGBoost – Average Operator)
# Synergy (red) + Antagonism (blue)
# High Quality, Helvetica, MD-style
# ============================================================

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import numpy as np

# --------------------- FONT & STYLE (Helvetica) ---------------------
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["figure.dpi"] = 300

# --------------------- DATA ---------------------
true_synergy = np.array(lbl_p_synergistic)
pred_synergy = np.hstack(synergy_probs["Average"])

true_ant = np.array(lbl_p_antagonistic)
pred_ant = np.hstack(antagonistic_probs["Average"])

# Compute ROC
fpr_s, tpr_s, _ = roc_curve(true_synergy, pred_synergy)
auc_s = auc(fpr_s, tpr_s)

fpr_a, tpr_a, _ = roc_curve(true_ant, pred_ant)
auc_a = auc(fpr_a, tpr_a)

# --------------------- PLOT ---------------------
fig, ax = plt.subplots(figsize=(6, 6))

# ROC – Synergy (red)
ax.plot(
    fpr_s, tpr_s,
    color='red', lw=2.5,
    solid_capstyle='round', solid_joinstyle='round'
)

# ROC – Antagonism (blue)
ax.plot(
    fpr_a, tpr_a,
    color='blue', lw=2.5,
    solid_capstyle='round', solid_joinstyle='round'
)

# Diagonal reference
ax.plot([0, 1], [0, 1], '--', color='gray', lw=1)

# Axes
ax.set_xlabel("1 - Specificity", fontsize=14, fontweight="bold")
ax.set_ylabel("Sensitivity", fontsize=14, fontweight="bold")
ax.set_title("ROC Curve – XGBoost (Average Operator)", fontsize=16, fontweight="bold")

# MD-style grid
ax.grid(alpha=0.25, linestyle='--', linewidth=0.6)

# Tick styling
ax.tick_params(axis='both', labelsize=12, width=1.2)

# Border thickness
for spine in ax.spines.values():
    spine.set_linewidth(1.2)

# Legend (Helvetica)
ax.legend(fontsize=12, frameon=False, loc="lower right")

plt.tight_layout()
plt.show()

# --------------------- SAVE ---------------------
fig.savefig("XGB_ROC_Average_SinglePlot.png", dpi=600, bbox_inches="tight")
# fig.savefig("XGB_ROC_Average_SinglePlot.pdf", dpi=600, bbox_inches="tight")
