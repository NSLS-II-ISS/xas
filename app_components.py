import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

import plotly.graph_objects as go

from xas.tiled_io import group_node_by_metadata_key

def build_scangroup_interactable(scangroup_node, group_label):
    select_all = html.Div([
        dbc.Checkbox(id={"type": "select_all", "group": group_label}, style={"display": "inline-block"}),
        html.Div("select all", style={"display": "inline-block", "padding": "3px"}),
    ])
    
    scan_labels = [html.Div([
            dbc.Checkbox(id={"type": "scan_check", "uid": k, "group": group_label}, style={"display": "inline-block"}),
            html.Div(v.metadata["scan_id"],
                     style={"display": "inline-block", "padding": "3px"},), 
            html.Br(),
        ]) 
        for k, v in scangroup_node.items_indexer
    ]
    return [select_all] + scan_labels
    # return scan_labels


def build_nested_accordion(base_node, sort_keys: list[str], _node_label=""):
    current_key = sort_keys[0]
    next_nodes, next_labels = group_node_by_metadata_key(base_node, current_key, return_values=True)
    next_level_keys = sort_keys[1:]

    # reached final level of sorting
    if len(next_level_keys) == 0:
        accordion_items = [
            dbc.AccordionItem(
                build_scangroup_interactable(sg_node, group_label=(_node_label+sg_label)),
                title=sg_label
            )
            for sg_node, sg_label in zip(next_nodes, next_labels)
        ]

    # recursively build next level of structure
    else:
        accordion_items = [
            dbc.AccordionItem(
                build_nested_accordion(sub_node, next_level_keys, _node_label=(_node_label+sub_label)),
                title=sub_label
            )
            for sub_node, sub_label in zip(next_nodes, next_labels)
        ]

    return dbc.Accordion(accordion_items, 
                         start_collapsed=True, 
                         always_open=True,)


def build_proposal_accordion(proposal_node, sort_key):
    if sort_key == "sample_name":
        proposal_accordion = build_nested_accordion(proposal_node, ("sample_name", "monochromator_scan_uid"))
    elif sort_key == "monochromator_scan_uid":
        proposal_accordion = build_nested_accordion(proposal_node, ("monochromator_scan_uid", "sample_name"))
    return html.Div(proposal_accordion, style={"max-height": "700px", "overflow-y": "scroll"})


visualization_tab = dbc.Tab([
        dcc.Store(id="previous_plot_data"),
        dbc.Row(dcc.Graph(figure=go.Figure(layout={"height": 800}), id="spectrum_plot")),
    ],
    label="Visualization",
)