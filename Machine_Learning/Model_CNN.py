# Cross-validation for
# 'Prediction of antibacterial interaction between essential oils via graph embedding approach'

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
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import auc
import matplotlib.pyplot as plt

# ======================================================================
# PARAMETERS
# ======================================================================

data_dir = "data/"
walk_length = 10
number_of_walks = 3
batch_size = 32
epochs = 100
layer_sizes = [16]
kfold = 5

# ======================================================================
# READ DATA
# ======================================================================

edgelist = pd.read_csv(os.path.join(data_dir, "pair_Sa"),
                       sep="\t", header=None,
                       names=["source", "target", "label", "eo1", "eo2"])

file_nodes = pd.read_csv(os.path.join(data_dir, "pair_content_Sa"), sep="\t")
nodes = file_nodes.set_index("ID")
feats = nodes.columns[:-1]

G_all_nx = nx.from_pandas_edgelist(edgelist, edge_attr="label")
nx.set_node_attributes(G_all_nx, "label", "label")
G_all = sg.StellarGraph.from_networkx(G_all_nx, node_features=nodes[feats])

# ======================================================================
# STORAGE
# ======================================================================

lbl1_p = []
lbl2_p = []
av_prob1 = []

# We will store the CNN model from last fold for ROC
clf_avg_model = None

# ======================================================================
# CNN TRAINING FUNCTION
# ======================================================================

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, Flatten, Dense
from tensorflow.keras.optimizers import Adam

def build_cnn(input_dim):
    model = Sequential([
        Conv1D(filters=16, kernel_size=3, activation='relu', input_shape=(input_dim, 1)),
        Flatten(),
        Dense(8, activation='relu'),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(0.001), loss='binary_crossentropy')
    return model

# ======================================================================
# K-FOLD LOOP
# ======================================================================

for k in range(kfold):

    print(f"Processing Fold {k+1}/{kfold}...")

    edge_gt = pd.read_csv(os.path.join(data_dir, f"{k+1}_t12"),
                          sep="\t", header=None,
                          names=["source", "target", "label", "eo1", "eo2"])

    edge_t = pd.read_csv(os.path.join(data_dir, f"{k+1}_t"),
                          sep="\t", header=None,
                          names=["source", "target", "label", "eo1", "eo2"])

    edge_p = pd.read_csv(os.path.join(data_dir, f"{k+1}_p"),
                          sep="\t", header=None,
                          names=["source", "target", "label", "eo1", "eo2"])

    file_node_t = pd.read_csv(os.path.join(data_dir, f"{k+1}_tn12"), sep="\t")
    node_t = file_node_t.set_index("ID")

    # Build training graph
    G_t_nx = nx.from_pandas_edgelist(edge_gt, edge_attr="label")
    nx.set_node_attributes(G_t_nx, "label", "label")
    G_t = sg.StellarGraph.from_networkx(G_t_nx, node_features=node_t[feats])

    sampler = UnsupervisedSampler(G_t, nodes=list(G_t.nodes()),
                                  length=walk_length, number_of_walks=number_of_walks)

    generator = Attri2VecLinkGenerator(G_t, batch_size)
    attri2vec = Attri2Vec(layer_sizes=layer_sizes, generator=generator,
                          bias=False, normalize=None)

    x_inp, x_out = attri2vec.in_out_tensors()
    prediction = link_classification(output_dim=1, output_act="sigmoid")(x_out)

    model = keras.Model(inputs=x_inp, outputs=prediction)
    model.compile(optimizer=keras.optimizers.Adam(1e-2),
                  loss="binary_crossentropy")

    model.fit(generator.flow(sampler),
              epochs=epochs, verbose=0, shuffle=True)

    # Embedding model
    embedding_model = keras.Model(inputs=x_inp[0], outputs=x_out[0])
    node_ids = nodes.index
    node_gen = Attri2VecNodeGenerator(G_all, batch_size).flow(node_ids)
    emb = embedding_model.predict(node_gen, verbose=0)

    # ==================================================================
    # FEATURE EXTRACTION
    # ==================================================================

    av_feat_t = []
    av_feat_p = []

    for i in range(len(edge_t)):
        s = edge_t.iloc[i]["source"]
        t = edge_t.iloc[i]["target"]
        av_feat_t.append(((nodes.loc[s, feats] + nodes.loc[t, feats]) / 2).values.astype(float))

    for i in range(len(edge_p)):
        s = edge_p.iloc[i]["source"]
        t = edge_p.iloc[i]["target"]
        av_feat_p.append(((nodes.loc[s, feats] + nodes.loc[t, feats]) / 2).values.astype(float))

        if edge_p["label"][i] == 1:
            lbl1_p.append(1); lbl2_p.append(0)
        elif edge_p["label"][i] == 2:
            lbl1_p.append(0); lbl2_p.append(1)
        else:
            lbl1_p.append(0); lbl2_p.append(0)

    av_feat_t = np.array(av_feat_t)
    av_feat_p = np.array(av_feat_p)

    # ==================================================================
    # CNN MODEL (TRAIN + PREDICT FOR TEST)
    # ==================================================================

    y_train = (edge_t["label"] == 1).astype(int).values
    X_train = av_feat_t.reshape(len(av_feat_t), av_feat_t.shape[1], 1)
    X_test = av_feat_p.reshape(len(av_feat_p), av_feat_p.shape[1], 1)

    cnn = build_cnn(av_feat_t.shape[1])
    cnn.fit(X_train, y_train, epochs=10, batch_size=32, verbose=0)

    # Save model for final ROC block
    clf_avg_model = cnn

    # Test predictions
    av_prob1.extend(cnn.predict(X_test).flatten())

# ======================================================================
# PRINT AUC
# ======================================================================

auc_avg = roc_auc_score(lbl1_p, av_prob1)
print("\nAUC (CNN – Average Operator):", auc_avg)
# ===============================================================
# ===============================================================
# COMBINED ROC PLOT – CNN (Average Operator)
# Synergistic & Antagonistic on ONE plot
# ===============================================================

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# High-quality Helvetica settings
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["figure.dpi"] = 300
plt.rcParams["axes.linewidth"] = 1.3

# =========================
# Synergistic vs Rest
# =========================
true_synergy = np.array(lbl_p_synergistic)
pred_synergy = np.hstack(synergy_probs["Average"])

fpr_s, tpr_s, _ = roc_curve(true_synergy, pred_synergy)
auc_s = auc(fpr_s, tpr_s)

# =========================
# Antagonistic vs Rest
# =========================
true_antag = np.array(lbl_p_antagonistic)
pred_antag = np.hstack(antagonistic_probs["Average"])

fpr_a, tpr_a, _ = roc_curve(true_antag, pred_antag)
auc_a = auc(fpr_a, tpr_a)

# =========================
# PLOT BOTH
# =========================
fig, ax = plt.subplots(figsize=(6.4, 6.4))

# Synergy curve (red)
ax.plot(
    fpr_s, tpr_s,
    color="red", lw=2.5,
    label=f"Synergy"
)

# Antagonism curve (blue)
ax.plot(
    fpr_a, tpr_a,
    color="blue", lw=2.5,
    label=f"Antagonism"
)

# Diagonal reference
ax.plot([0, 1], [0, 1], '--', color="gray", lw=1)

# Labels & Title
ax.set_xlabel("1 - Specificity", fontsize=14, fontweight="bold")
ax.set_ylabel("Sensitivity", fontsize=14, fontweight="bold")
ax.set_title("ROC Curve – CNN (Average Operator)", fontsize=16, fontweight="bold")

# Grid + styling
ax.grid(alpha=0.25, linestyle='--', linewidth=0.6)
ax.legend(fontsize=12, frameon=False)

# Thicker spine borders
for spine in ax.spines.values():
    spine.set_linewidth(1.2)

fig.tight_layout()

# SAVE HIGH QUALITY FIG
fig.savefig("CNN_ROC_Combined_Average.png", dpi=600, bbox_inches="tight")

plt.show()
