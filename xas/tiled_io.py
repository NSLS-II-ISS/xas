from tiled.client import from_uri
from tiled.client.node import Node
from tiled.queries import Key

from collections import UserDict


class StoredTiledEntry:
    def __init__(self, client_obj):
        self.source_type = type(client_obj)
        self.data = client_obj.read()
        self.metadata = client_obj.metadata

    def read(self):
        """implemented to closer mimic tiled client object funcionality"""
        return self.data
    
    def __repr__(self):
        return f"Stored data and metadata from {self.source_type}"


class TiledReader(UserDict):
    """Custom dictionary for accessing data from a source tiled client node.
    When data is pulled from the node then it will be saved in the `TiledReader.data` 
    attribute and pulled from there in the future."""

    def __init__(self, source_node: Node):
        self.data = {}
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

        


def get_iss_sandbox():
    username=input("BNL Username: ")
    tiled_catalog = from_uri("https://localhost:8008", verify=False, username=username)
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


def sort_nodes_by_metadata_key(list_of_nodes: list, sort_key):
    return list_of_nodes.sort(key=(lambda node: node[sort_key]))


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


class TiledViewer:
    def __init__(self, tiled_client):
        self.tiled_client = tiled_client
