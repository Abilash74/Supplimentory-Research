# Cross-validation for
# 'Prediction of antibacterial interaction between essential oils via graph embedding approach'
#

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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

import scipy.stats as st

data_dir = "data/"   # working directory
walk_length = 10
number_of_walks = 3
batch_size = 32
epochs = 100
layer_sizes = [16]
kfold = 5

### Read edge list 
edgelist = pd.read_csv( os.path.join(data_dir, "pair_Sa"), 
  sep="\t", header=None, names=["source", "target", "label", "eo1", "eo2"])

### Read node list (The last column is the label)
file_nodes = pd.read_csv( os.path.join(data_dir, "pair_content_Sa"), sep="\t" )
nodes = file_nodes.set_index('ID')
feats = nodes.columns[0:(len(nodes.columns)-1)]

G_all_nx = nx.from_pandas_edgelist(edgelist, edge_attr="label")
nx.set_node_attributes(G_all_nx, "label", "label")
G_all = sg.StellarGraph.from_networkx(G_all_nx, node_features=nodes[feats])

# print(G_all.info())


lbl0 = []; lbl1_p = []; lbl2_p = [];
h_pred1 = [];  L1_pred1 = []; L2_pred1 = []; av_pred1 = []; 
h_prob1 = [];  L1_prob1 = []; L2_prob1 = []; av_prob1 = []; 
h_prob2 = [];  L1_prob2 = []; L2_prob2 = []; av_prob2 = []; 
cmp_av_prob1 = []; cmp_av_pred1 = []; cmp_av_prob2 = []; 

for k in range( kfold ):
  ### Read splitted (train/predict) interaction data
  edge_gt = pd.read_csv( os.path.join(data_dir, str(k+1) + "_t12"),
    sep="\t", header=None, names=["source", "target", "label", "eo1", "eo2"] )
  edge_t = pd.read_csv( os.path.join(data_dir,  str(k+1) + "_t"),
    sep="\t", header=None, names=["source", "target", "label", "eo1", "eo2"] )
  edge_p = pd.read_csv( os.path.join(data_dir,  str(k+1) + "_p"),
    sep="\t", header=None, names=["source", "target", "label", "eo1", "eo2"] )
  file_node_t = pd.read_csv( os.path.join(data_dir, str(k+1) + "_tn12"), sep="\t")
  node_t = file_node_t.set_index('ID')
  G_t_nx = nx.from_pandas_edgelist(edge_gt, edge_attr="label")
  nx.set_node_attributes(G_t_nx, "label", "label")
  G_t_node_features = node_t[feats]
  G_t = sg.StellarGraph.from_networkx(G_t_nx, node_features=G_t_node_features)
  unsupervised_samples = UnsupervisedSampler(
    G_t, nodes=list(G_t.nodes()), length=walk_length, number_of_walks=number_of_walks
  )
  generator = Attri2VecLinkGenerator(G_t, batch_size)
  attri2vec = Attri2Vec(layer_sizes=layer_sizes, generator=generator, bias=False, normalize=None)
  ### Build the model and expose input and output sockets of attri2vec, for node pair inputs:
  x_inp, x_out = attri2vec.in_out_tensors()
  prediction = link_classification(
    output_dim=1, output_act="sigmoid", edge_embedding_method="ip"
  )(x_out)
  model = keras.Model(inputs=x_inp, outputs=prediction)
  model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-2),
    loss=keras.losses.binary_crossentropy,
    metrics=[keras.metrics.binary_accuracy],
  )
  history = model.fit(
    generator.flow(unsupervised_samples),
    epochs=epochs,
    verbose=2,
    use_multiprocessing=False,
    workers=1,
    shuffle=True,
  )
  x_inp_src = x_inp[0]
  x_out_src = x_out[0]
  embedding_model = keras.Model(inputs=x_inp_src, outputs=x_out_src)
  ### Node embedding
  node_ids = nodes.index
  node_gen = Attri2VecNodeGenerator(G_all, batch_size).flow(node_ids)
  emb = embedding_model.predict(node_gen, workers=4, verbose=1)
  ### Selection of binary operator
  h_feat_t = [];  h_feat_p = [];  L1_feat_t = []; L1_feat_p = []; 
  L2_feat_t = []; L2_feat_p = []; av_feat_t = []; av_feat_p = []; 
  cmp_av_feat_t = []; cmp_av_feat_p = []; 
  ### Edge fearure
  for i in range(len(edge_t)):
    n1 = np.where( nodes.index == edge_t["source"][i] )[0].tolist()
    n2 = np.where( nodes.index == edge_t["target"][i] )[0].tolist()
    h_feat_t.append( np.ravel( emb[n1] * emb[n2] ) )
    L1_feat_t.append( np.ravel(np.abs(emb[n1] - emb[n2]) ))
    L2_feat_t.append( np.ravel((emb[n1] - emb[n2]) ** 2 ))
    av_feat_t.append( np.ravel((emb[n1] + emb[n2]) / 2 ))
    cmp_av_feat_t.append( (nodes[feats].loc[edge_t["source"][i],:] + nodes[feats].loc[edge_t["target"][i],:]) / 2)
  for i in range(len(edge_p)):
    n1 = np.where( nodes.index == edge_p["source"][i] )[0].tolist()
    n2 = np.where( nodes.index == edge_p["target"][i] )[0].tolist()
    h_feat_p.append( np.ravel( emb[n1] * emb[n2] ) )
    L1_feat_p.append( np.ravel(np.abs(emb[n1] - emb[n2]) ))
    L2_feat_p.append( np.ravel((emb[n1] - emb[n2]) ** 2 ))
    av_feat_p.append( np.ravel((emb[n1] + emb[n2]) / 2 ))
    cmp_av_feat_p.append( (nodes[feats].loc[edge_p["source"][i],:] + nodes[feats].loc[edge_p["target"][i],:]) / 2)
    if edge_p["label"][i] == 1:
      lbl1_p.append( 1 )
      lbl2_p.append( 0 )
    elif edge_p["label"][i] == 2:
      lbl1_p.append( 0 )
      lbl2_p.append( 1 )
    else:
      lbl1_p.append( 0 )
      lbl2_p.append( 0 )
  lbl0.extend( edge_p["label"] )
  ### Hadamard
  clf_edge_pred_from_feat = LogisticRegression(verbose=0, solver="lbfgs", class_weight='balanced', max_iter=500, multi_class='multinomial')
  clf_edge_pred_from_feat.fit(h_feat_t, edge_t["label"])
  edge_pred_from_feat = clf_edge_pred_from_feat.predict_proba(h_feat_p)
  positive_class_index = 0
  negative_class_index = 1
  if clf_edge_pred_from_feat.classes_[1] == 1:
    positive_class_index = 1
  elif clf_edge_pred_from_feat.classes_[2] == 1:
    positive_class_index = 2
  if clf_edge_pred_from_feat.classes_[0] == 2:
    negative_class_index = 0
  elif clf_edge_pred_from_feat.classes_[2] == 2:
    negative_class_index = 2
  h_prob1.extend( edge_pred_from_feat[:, positive_class_index] )
  h_pred1.extend( clf_edge_pred_from_feat.predict(h_feat_p) )
  h_prob2.extend( edge_pred_from_feat[:, negative_class_index] )
  ### L1-norm
  clf_edge_pred_from_feat = LogisticRegression(verbose=0, solver="lbfgs", class_weight='balanced', max_iter=500, multi_class='multinomial')
  clf_edge_pred_from_feat.fit(L1_feat_t, edge_t["label"])
  edge_pred_from_feat = clf_edge_pred_from_feat.predict_proba(L1_feat_p)
  positive_class_index = 0
  negative_class_index = 1
  if clf_edge_pred_from_feat.classes_[1] == 1:
    positive_class_index = 1
  elif clf_edge_pred_from_feat.classes_[2] == 1:
    positive_class_index = 2
  if clf_edge_pred_from_feat.classes_[0] == 2:
    negative_class_index = 0 
  elif clf_edge_pred_from_feat.classes_[2] == 2:
    negative_class_index = 2
  L1_prob1.extend( edge_pred_from_feat[:, positive_class_index] )
  L1_pred1.extend( clf_edge_pred_from_feat.predict(L1_feat_p) )
  L1_prob2.extend( edge_pred_from_feat[:, negative_class_index] )
  ### L2-norm
  clf_edge_pred_from_feat = LogisticRegression(verbose=0, solver="lbfgs", class_weight='balanced', max_iter=500, multi_class='multinomial')
  clf_edge_pred_from_feat.fit(L2_feat_t, edge_t["label"])
  edge_pred_from_feat = clf_edge_pred_from_feat.predict_proba(L2_feat_p)
  positive_class_index = 0
  negative_class_index = 1
  if clf_edge_pred_from_feat.classes_[1] == 1:
    positive_class_index = 1
  elif clf_edge_pred_from_feat.classes_[2] == 1:
    positive_class_index = 2
  if clf_edge_pred_from_feat.classes_[0] == 2:
    negative_class_index = 0
  elif clf_edge_pred_from_feat.classes_[2] == 2:
    negative_class_index = 2
  L2_prob1.extend( edge_pred_from_feat[:, positive_class_index] )
  L2_pred1.extend( clf_edge_pred_from_feat.predict(L2_feat_p) )
  L2_prob2.extend( edge_pred_from_feat[:, negative_class_index] )
  ### average
  clf_edge_pred_from_feat = LogisticRegression(verbose=0, solver="lbfgs", class_weight='balanced', max_iter=500, multi_class='multinomial')
  clf_edge_pred_from_feat.fit(av_feat_t, edge_t["label"])
  edge_pred_from_feat = clf_edge_pred_from_feat.predict_proba(av_feat_p)
  positive_class_index = 0
  negative_class_index = 1
  if clf_edge_pred_from_feat.classes_[1] == 1:
    positive_class_index = 1
  elif clf_edge_pred_from_feat.classes_[2] == 1:
    positive_class_index = 2
  if clf_edge_pred_from_feat.classes_[0] == 2:
    negative_class_index = 0
  elif clf_edge_pred_from_feat.classes_[2] == 2:
    negative_class_index = 2
  av_prob1.extend( edge_pred_from_feat[:, positive_class_index] )
  av_pred1.extend( clf_edge_pred_from_feat.predict(av_feat_p) )
  av_prob2.extend( edge_pred_from_feat[:, negative_class_index] )
  ### Traditional classification-based (average operator)
  clf_edge_pred_from_feat = LogisticRegression(verbose=0, solver="lbfgs", class_weight='balanced', max_iter=500, multi_class='multinomial')
  clf_edge_pred_from_feat.fit(cmp_av_feat_t, edge_t["label"])
  edge_pred_from_feat = clf_edge_pred_from_feat.predict_proba(cmp_av_feat_p)
  positive_class_index  = 0
  negative_class_index = 1
  if clf_edge_pred_from_feat.classes_[1] == 1:
    positive_class_index = 1
  elif clf_edge_pred_from_feat.classes_[2] == 1:
    positive_class_index = 2
  if clf_edge_pred_from_feat.classes_[0] == 2:
    negative_class_index = 0
  elif clf_edge_pred_from_feat.classes_[2] == 2:
    negative_class_index = 2
  cmp_av_prob1.extend( edge_pred_from_feat[:, positive_class_index] )
  cmp_av_pred1.extend( clf_edge_pred_from_feat.predict(cmp_av_feat_p) )
  cmp_av_prob2.extend( edge_pred_from_feat[:, negative_class_index] )

print("AUC (synergistic-versus-rest):\n")
print("Hadamard\t"   + str(roc_auc_score(lbl1_p, h_prob1)))
print("L1-norm\t"  + str(roc_auc_score(lbl1_p, L1_prob1)))
print("L2-norm\t"  + str(roc_auc_score(lbl1_p, L2_prob1)))
print("average\t" + str(roc_auc_score(lbl1_p, av_prob1)))
print("classication-based\t" + str(roc_auc_score(lbl1_p, cmp_av_prob1)))

print("AUC (antagonistic-versus-rest):\n")
print("Hadamard\t"   + str(roc_auc_score(lbl2_p, h_prob2)))
print("L1-norm\t"  + str(roc_auc_score(lbl2_p, L1_prob2)))
print("L2-norm\t"  + str(roc_auc_score(lbl2_p, L2_prob2)))
print("average\t" + str(roc_auc_score(lbl2_p, av_prob2)))
print("classication-based\t" + str(roc_auc_score(lbl2_p, cmp_av_prob2)))


# ===============================================================
# COMBINED ROC PLOT – MAR Logistic Regression (Average Operator)
# Synergistic & Antagonistic ROC Curves in ONE Figure
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
# Synergy (label 1 vs rest)
# =========================
true_synergy = np.array(lbl1_p)           # true synergy labels (0/1)
pred_synergy = np.array(av_prob1)         # MAR probability for synergy

fpr_s, tpr_s, _ = roc_curve(true_synergy, pred_synergy)
auc_s = auc(fpr_s, tpr_s)


# =========================
# Antagonism (label 2 vs rest)
# =========================
true_antag = np.array(lbl2_p)             # true antagonism labels (0/1)
pred_antag = np.array(av_prob2)           # MAR probability for antagonism

fpr_a, tpr_a, _ = roc_curve(true_antag, pred_antag)
auc_a = auc(fpr_a, tpr_a)


# =========================
# PLOT BOTH
# =========================
fig, ax = plt.subplots(figsize=(6.4, 6.4))

# ROC curve – Synergy (red)
ax.plot(
    fpr_s, tpr_s,
    color="red", lw=2.5,
    label=f"Synergy"
)

# ROC curve – Antagonism (blue)
ax.plot(
    fpr_a, tpr_a,
    color="blue", lw=2.5,
    label=f"Antagonism"
)

# Reference diagonal
ax.plot([0, 1], [0, 1], '--', color="gray", lw=1)

# Axis labels + title
ax.set_xlabel("1 - Specificity", fontsize=14, fontweight="bold")
ax.set_ylabel("Sensitivity", fontsize=14, fontweight="bold")
ax.set_title("ROC Curve – MLR Model (Average Operator)", fontsize=16, fontweight="bold")

# Grid + spines
ax.grid(alpha=0.25, linestyle='--', linewidth=0.6)
for spine in ax.spines.values():
    spine.set_linewidth(1.2)

# Legend
ax.legend(frameon=False, fontsize=12)

fig.tight_layout()

# Save the figure
fig.savefig("MLR_ROC_Combined_Average.png", dpi=600, bbox_inches="tight")

plt.show()
