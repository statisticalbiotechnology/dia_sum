#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 11:28:15 2022

@author: ptruong
"""


def get_clusters(ctrl_lt_file = "proteins.1vs2.tsv",
                 ctrl_st_file = "proteins.1vs3.tsv",
                 lt_st_file = "proteins.2vs3.tsv"):
    cols = ['q_value', 'log2_fold_change', "upregulated"]
    
    ctrl_lt = parse_triqler(ctrl_lt_file).set_index("protein")
    ctrl_st = parse_triqler(ctrl_st_file).set_index("protein")
    lt_st = parse_triqler(lt_st_file).set_index("protein")
    
    ctrl_lt["upregulated"] = ctrl_lt.log2_fold_change > 0
    ctrl_st["upregulated"] = ctrl_st.log2_fold_change > 0 
    lt_st["upregulated"] = lt_st.log2_fold_change > 0
    
    
    ctrl_lt = ctrl_lt[cols].rename({"q_value":"ctrl-lt:q_value", "upregulated":"ctrl-lt:upregulated", "log2_fold_change":"ctrl-lt:log2_fold_change"}, axis = 1)
    ctrl_st = ctrl_st[cols].rename({"q_value":"ctrl-st:q_value", "upregulated":"ctrl-st:upregulated", "log2_fold_change":"ctrl-st:log2_fold_change"}, axis = 1)
    lt_st = lt_st[cols].rename({"q_value":"lt-st:q_value", "upregulated":"lt-st:upregulated", "log2_fold_change":"lt-st:log2_fold_change"}, axis = 1)
    
    df = pd.concat([ctrl_lt, ctrl_st, lt_st], axis = 1)
    
    # Clusters 
    C1 = ((df["ctrl-lt:upregulated"] == True) & (df["ctrl-st:upregulated"] == True) & (df["lt-st:upregulated"] == False))
    C2 = ((df["ctrl-lt:upregulated"] == False) & (df["ctrl-st:upregulated"] == False) & (df["lt-st:upregulated"] == True))
    C3 = (df["lt-st:upregulated"] == False)      
    C4 = (df["lt-st:upregulated"] == True)
    C5 = (df["ctrl-st:upregulated"] == True)
    C6 = (df["ctrl-st:upregulated"] == False)
    return C1, C2, C3, C4, C5, C6

ctrl_lt = parse_triqler("proteins.1vs2.tsv").set_index("protein")
ctrl_st = parse_triqler("proteins.1vs3.tsv").set_index("protein")
lt_st = parse_triqler("proteins.2vs3.tsv").set_index("protein")


fdr_threshold = 0.05
# Check ctrl-lt C1,C2
ctrl_lt[C1 & (ctrl_lt.q_value < fdr_threshold)]
ctrl_lt[C2 & (ctrl_lt.q_value < fdr_threshold)]

# Check lt-st C3, c4
lt_st[C1 & (lt_st.q_value < fdr_threshold)]
lt_st[C2 & (lt_st.q_value < fdr_threshold)]

# Check ctrl-st C5, C6
ctrl_lt[C1 & (ctrl_lt.q_value < fdr_threshold)]
ctrl_lt[C2 & (ctrl_lt.q_value < fdr_threshold)]





