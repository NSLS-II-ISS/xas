from tiled.client import from_uri
from tiled.client.node import Node
from tiled.queries import Key
from tiled.client.cache import Cache

from collections import UserDict

# Outlining TiledReader object
# may be absolete? https://blueskyproject.io/tiled/tutorials/caching.html
#
# def PROCESS_MY_DATASET(data):
#     data_type = 'xas'

#     if data_type == 'xas':
#         default_kwargs = {.....}
#         processed_data = larch_normalize_this_somehow(data, **default_kwargs) # ideally at some point AIML people will take care of kwargs and whatnot


#     return processed_data,

# def function_for_xas_processing():
#     pass


# mapping_of_processing = {'xas' : function_for_xas_processing}


# class StoredTiledEntry(UserDict):
#     def __init__(self, client_obj):
#         super().__init__()
#         self.client_obj = client_obj
#         # self.source_type = type(client_obj)
#         # self.data = client_obj.read()
#         # self.metadata = client_obj.metadata

#     # def read(self):
#     #     """implemented to closer mimic tiled client object funcionality"""
#     #     return self.data
    
#     def __repr__(self):
#         return f"Stored data and metadata from {self.source_type}"

#     def __getitem__(self, kind):
#         if kind in self.data.keys():
#             return self.data[kind]
#         else:
#             if kind == 'binned':
#                 data = self.client_obj.read()
#                 metadata = self.client_obj.metadata

#             elif kind == 'processed':
#                 data, metadata = PROCESS_MY_DATASET(self.data['binned'])
#                 self.data[kind] = {'data': data, 'metadata': metadata}
#             else:
#                 raise KeyError('NOT IMPLEMENTED')
#                 # see if you have it

#     # def process(self):


# class TiledReader(UserDict):
#     """Custom dictionary for accessing data from a source tiled client node.
#     When data is pulled from the node then it will be saved in the `TiledReader.data` 
#     attribute and pulled from there in the future."""

#     def __init__(self, source_node: Node):
#         super().__init__()
#         # self.data = {}
#         self.source_node = source_node

#     def __setitem__(self, *args):
#         raise RuntimeError("Cannot set uid values using TiledReader")
    
#     def __getitem__(self, uid: str):
#         if uid in self.data.keys():
#             return self.data[uid]
#         elif uid in self.source_node.keys():
#             self.data[uid] = StoredTiledEntry(self.source_node[uid])
#             return self.data[uid]
#         else:
#             raise KeyError("Could not find uid in locally stored data or source node")

        


def get_iss_sandbox():
    username=input("BNL Username: ")
    # cache maximum of 2GB in RAM
    tiled_catalog = from_uri("https://localhost:8008", verify=False, username=username, cache=Cache.in_memory(2e9))
    iss_sandbox_node = tiled_catalog["iss"]["sandbox"]
    return iss_sandbox_node


def filter_node_by_metadata_key(node, key, value):
    key = str(key)
    value = str(value)
    return node.search(Key(key) == value)


def filter_node_for_proposal(node, year, cycle, proposal):
    year = str(year)
    cycle = str(cycle)
    proposal = str(proposal)
    return node.search(Key('year') == year).search(Key('cycle') == cycle).search(Key('PROPOSAL') == proposal)


def min_value_for_metadata_key(node: Node, key):
    """Loop over metadata for entries in node and find the minimum value of metadata[key]"""
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
    """Return a list of sub-nodes each of which correspond to a unique value for the specified key"""
    unique_values = sorted(list(set(v.metadata[key] for v in node.values_indexer)))
    grouped_nodes = [node.search(Key(key) == uv) for uv in unique_values]
    if return_values:
        return grouped_nodes, unique_values
    else:
        return grouped_nodes


def get_df_for_uid(node, uid):
    singleton_node = node.search(Key("uid") == uid)
    dfc = singleton_node[singleton_node.keys()[0]]
    return dfc.read()


class StoredTiledEntry:
    def __init__(self, client_obj):
        self.source_type = type(client_obj)
        self.data = client_obj.read()
        self.metadata = client_obj.metadata
        # self.processed_data = {}

    def read(self):
        """implemented to closer mimic tiled client object funcionality"""
        return self.data
    
    def __repr__(self):
        return f"Stored data and metadata from {self.source_type}"


class TiledReader(UserDict):
    """Custom dictionary for accessing data from a source tiled client node.
    When data is pulled from the node it will also be saved in the `TiledReader.data` 
    attribute and pulled from there in the future."""

    def __init__(self, source_node: Node):
        super().__init__()
        self.source_node = source_node

    def __setitem__(self, *args):
        raise RuntimeError("Cannot set uid values using TiledReader")
    
    def __getitem__(self, uid: str):
        if uid in self.data.keys():
            return self.data[uid]
        elif uid in self.source_node.keys():
            self.data[uid] = StoredTiledEntry(self.source_node[uid])
            return self.data[uid]
        else:
            raise KeyError("Could not find uid in locally stored data or source node")
