#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


import pandas as pd
from collections import defaultdict

def updateMeta(config):
    _sanity_checks(config)
    if config['metasheet'].endswith(".csv"):
        metadata = pd.read_csv(config['metasheet'], sep=',')
    if config['metasheet'].endswith(".tsv"):
        metadata = pd.read_csv(config['metasheet'], sep='\t')
    if config['metasheet'].endswith(".txt"):
        metadata = pd.read_csv(config['metasheet'], sep=' ')

    config["comparisons"] = ["%s" %(row["comp_id"]) for index, row in metadata.iterrows()]
    config["comps"] = _get_comp_info(metadata)
#     config["metacols"] = [c for c in metadata.columns if c.lower()[:4] != 'comp']
#     config["file_info"] = { sampleName : config["samples"][sampleName] for sampleName in metadata.index }
    config["ordered_sample_list"] = [ sample for sample in config["samples"]]
    return config


def _sanity_checks(config):
    #metasheet pre-parser: converts dos2unix, catches invalid chars
    _invalid_map = {'\r':'\n', '(':'.', ')':'.', ' ':'_', '/':'.', '$':''}
    _meta_f = open(config['metasheet'])
    _meta = _meta_f.read()
    _meta_f.close()

    _tmp = _meta.replace('\r\n','\n')
    #check other invalids
    for k in _invalid_map.keys():
        if k in _tmp:
            _tmp = _tmp.replace(k, _invalid_map[k])

    #did the contents change?--rewrite the metafile
    if _meta != _tmp:
        #print('converting')
        _meta_f = open(config['metasheet'], 'w')
        _meta_f.write(_tmp)
        _meta_f.close()


def _get_comp_info(metadata):
#     comps_info = defaultdict(dict)
    comps_info = {}
    for index, row in metadata.iterrows():
        comp = row["comp_id"]
        normal = row["normal"]
        tumor  = row["tumor"]

        comps_info[comp] = {}
        comps_info[comp]['normal'] = normal
        comps_info[comp]['tumor'] = tumor
        if "rna" in row:
            rna  = row["rna"]
            comps_info[comp]['rna'] = rna
    return comps_info

