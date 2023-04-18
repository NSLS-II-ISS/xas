from tiled.client import from_uri, from_profile
from tiled.client.node import Node
from tiled.queries import Key
from tiled.client.cache import Cache

import pandas as pd

from collections import UserDict, namedtuple

from xas.xdash_math import LarchCalculator, calc_mus


_LABEL_DICT = {'mu': 'mu',
               'normalized': 'mu norm',
               'flattened': 'mu flat'}

class DataManager:
    def __init__(self, source_node: Node):
        self.source_node = source_node
        self.processed_data = {}
        
        self.data_types = ("xas", "xes", "rixs")


    def _process_data(self, uid, channel, params=None):
        if params is None:
            params = dict()

        data = self.source_node[uid].read()
        metadata = self.source_node[uid].metadata
        # data_type = metadata["type"]
        data_type = "xas"  # set to xas for now until md is labeled in tiled
        
        if data_type.lower() in self.data_types:
            if data_type.lower() == "xas":
                calc_mus(data)
                
                energy = data["energy"]
                mu_in = data[channel]

                norm_result = pd.DataFrame()
                norm_result["energy"] = energy
                norm_result["norm"], norm_parameters = LarchCalculator.normalize(energy, 
                                                                                  mu_in, 
                                                                                  flatten_output=False, 
                                                                                  return_norm_parameters=True,
                                                                                  **params)
                
                norm_result["flat"] = LarchCalculator.normalize(energy, mu_in, **norm_parameters)
                return {"data": norm_result, "parameters": norm_parameters}
        else:
            raise ValueError("unsupported data type")

    def get_metadata(self, uid):
        return self.source_node[uid].metadata

    def get_data(self, uid, channel, kind='mu', processing_parameters=None):
        energy = self.get_raw_data(uid)["energy"]
        if kind == "mu":
            mu = self.get_raw_data(uid)[channel]
        elif kind == "normalized":
            mu = self.get_processed_data(uid, channel, processing_parameters=processing_parameters)["norm"]
        elif kind == "flattened":
            mu = self.get_processed_data(uid, channel, processing_parameters=processing_parameters)["flat"]
        return energy, mu

    def get_plotting_data(self, uid, channel, kind='mu', processing_parameters=None):
        energy, mu = self.get_data(uid, channel, kind=kind, processing_parameters=processing_parameters)
        scan_id = self.get_metadata(uid)["scan_id"]
        label = f"{channel} {kind} {scan_id}" if kind != 'mu' else f"{channel} {scan_id}"
        return energy, mu, label

    def get_raw_data(self, uid):
        # bad, fix later
        df = self.source_node[uid].read()
        calc_mus(df)
        return df
    
    def get_processed_data(self, uid, channel, processing_parameters=None):
        if uid in self.processed_data.keys():
            # check if channel has not yet been calculated or if the processing parameters are new
            if channel not in self.processed_data[uid].keys() or processing_parameters != self.get_processing_parameters(uid, channel):
                self.processed_data[uid][channel] = self._process_data(uid, channel, params=processing_parameters)
        else:
            self.processed_data[uid] = dict()
            self.processed_data[uid][channel] = self._process_data(uid, channel, params=processing_parameters)
        return self.processed_data[uid][channel]["data"]
    
    def get_processing_parameters(self, uid, channel):
        if uid in self.processed_data.keys():
            if channel not in self.processed_data[uid].keys():
                self.processed_data[uid][channel] = self._process_data(uid, channel)
        else:
            self.processed_data[uid] = dict()
            self.processed_data[uid][channel] = self._process_data(uid, channel)
        params = self.processed_data[uid][channel]["parameters"]
        return params


def get_iss_sandbox():
    username=input("BNL Username: ")
    # cache maximum of 2GB in RAM
    try:
        tiled_catalog = from_profile('nsls2', verify=False, username=username, cache=Cache.in_memory(2e9))
    except:
        tiled_catalog = from_uri("https://localhost:8008", verify=False, username=username, cache=Cache.in_memory(2e9))
    # tiled_catalog = from_uri("https://localhost:8008", verify=False, username=username, cache=Cache.in_memory(2e9))
    iss_sandbox_node = tiled_catalog["iss"]["sandbox"]
    return iss_sandbox_node


def filter_node_by_metadata_key(node, key, value):
    key = str(key)
    value = str(value)
    return node.search(Key(key) == value)


def filter_node_for_proposal(node, year, cycle, proposal, cutoff=50):
    year = str(year)
    cycle = str(cycle)
    proposal = str(proposal)
    # try:
    r = node.search(Key('year') == year).search(Key('cycle') == cycle).search(Key('proposal') == proposal)
    cutoff_scan_id = r[r.keys()[len(r) - cutoff]].metadata['scan_id']
    return r.search(Key('scan_id') >= cutoff_scan_id)
    # except:
    # return node.search(Key('year') == year).search(Key('cycle') == cycle).search(Key('PROPOSAL') == proposal)

def min_value_for_metadata_key(node: Node, key):
    """Loop over entries in node and find the minimum value of metadata[key]"""
    return min([client_obj.metadata[key] for client_obj in node.values_indexer if key in client_obj.metadata.keys()], default=None)


def sort_nodes_by_metadata_key(list_of_nodes: list[Node], sort_key, node_labels:list=None, reverse_order=False):
    """Returns nodes sorted by the minimum value for the `sort_key` within each node.
    If `node_labels` is provided then the labels are sorted in the same order as the sorted nodes
    and also returned."""    
    if node_labels is None:
        # weird lambda returning tuple is a hacky way to push `None`s to the end of sorted list
        return sorted(list_of_nodes, key=lambda node: (min_value_for_metadata_key(node, sort_key) is None, min_value_for_metadata_key(node, sort_key)), 
                      reverse=reverse_order)
    else:
        sorted_nodes_and_labels = sorted(
            zip(list_of_nodes, node_labels), key=lambda node: (min_value_for_metadata_key(node[0], sort_key) is None, min_value_for_metadata_key(node[0], sort_key)),
            reverse=reverse_order,
        )
        sorted_nodes, sorted_labels = zip(*sorted_nodes_and_labels)
        return list(sorted_nodes), list(sorted_labels)


def group_node_by_metadata_key(node, key, return_values=False):
    """Return a list of sub-nodes each of which correspond to a unique value for the specified key.
    Also effectively removes any entries which do not have values for the specified metadata key."""
    unique_values = sorted(list(set(v.metadata[key] for v in node.values_indexer if key in v.metadata.keys())))
    grouped_nodes = [node.search(Key(key) == uv) for uv in unique_values]
    if return_values:
        return grouped_nodes, unique_values
    else:
        return grouped_nodes


def get_df_for_uid(node, uid):
    singleton_node = node.search(Key("uid") == uid)
    dfc = singleton_node[singleton_node.keys()[0]]
    return dfc.read()


def get_unique_values_for_key(node, key):
    values = []
    for uid in node:
        if key in node[uid].metadata.keys():
            value = node[uid].metadata[key]
            values.append(value)
            next_node = node.search(Key(key) != value)
            if len(next_node) > 0:
                next_values = get_unique_values_for_key(next_node, key)
                values.extend(next_values)
        return values








