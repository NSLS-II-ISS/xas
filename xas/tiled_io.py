from tiled.client import from_uri
from tiled.queries import Key
# TILED_CATALOG = from_uri("https://localhost:8008", verify=False, username=input("BNL Username: \n"))
# ISS_SANDBOX = TILED_CATALOG["iss"]["sandbox"]


def get_iss_sandbox():
    username=input("BNL Username: ")
    tiled_catalog = from_uri("https://localhost:8008", verify=False, username=username)
    iss_sandbox_node = tiled_catalog["iss"]["sandbox"]
    return iss_sandbox_node


def filter_node_for_proposal(node, year, cycle, proposal):
    year = str(year)
    cycle = str(cycle)
    proposal = str(proposal)
    return node.search(Key('year') == year).search(Key('cycle') == cycle).search(Key('PROPOSAL') == proposal)


def sort_node_by_metadata_key(node, key, return_values=False):
    """Return a list of sub-nodes each of which correspond to a unique value for the specified key"""
    unique_values = sorted(list(set(v.metadata[key] for v in node.values_indexer)))
    sorted_nodes = [node.search(Key(key) == uv) for uv in unique_values]
    if return_values:
        return sorted_nodes, unique_values
    else:
        return sorted_nodes


def get_df_for_uid(node, uid):
    singleton_node = node.search(Key("uid") == uid)
    dfc = singleton_node[singleton_node.keys()[0]]
    return dfc.read()


class TiledViewer:
    def __init__(self, tiled_client):
        self.tiled_client = tiled_client
