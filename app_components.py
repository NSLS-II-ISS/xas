import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

import plotly.graph_objects as go

from xas.tiled_io import sort_node_by_metadata_key

def build_scangroup_interactable(scangroup_node):
    select_all = html.Div(dbc.Button("select all", color="secondary", id={"type": "select_all_btn"}),
                          style={"padding-bottom": "10px"})
    scan_labels = [html.Div([
            html.Div(v.metadata["scan_id"],
                     style={"display": "inline-block", "padding": "3px"},), 
            # dbc.Checkbox(id={"type": "scan_check", "group": group_label}, style={"display": "inline-block"}),
            dbc.Checkbox(id={"type": "scan_check", "uid": k}, style={"display": "inline-block"}),
            html.Br(),
        ]) 
        for k, v in scangroup_node.items_indexer
    ]
    return [select_all] + scan_labels
    # return scan_labels


def build_nested_accordions(base_node, *sort_keys):
    current_key = sort_keys[0]
    next_nodes, next_labels = sort_node_by_metadata_key(base_node, current_key, return_values=True)
    next_level_keys = sort_keys[1:]
    if len(next_level_keys) == 0:
        accordion_items = [
            dbc.AccordionItem(
                build_scangroup_interactable(sg_node),
                title=sg_label
            )
            for sg_node, sg_label in zip(next_nodes, next_labels)
        ]
    else:
        accordion_items = [
            dbc.AccordionItem(
                build_nested_accordions(sub_node, *next_level_keys),
                title=sub_label
            )
            for sub_node, sub_label in zip(next_nodes, next_labels)
        ]
    return dbc.Accordion(accordion_items, start_collapsed=True, always_open=True)



def build_proposal_accordion(proposal_node, sort_key):
    if sort_key == "sample_name":
        proposal_accordion = build_nested_accordions(proposal_node, "sample_name", "monochromator_scan_uid")
    elif sort_key == "monochromator_scan_uid":
        proposal_accordion = build_nested_accordions(proposal_node, "monochromator_scan_uid", "sample_name")
    return proposal_accordion


visualization_tab = dbc.Tab(
    [
    # dbc.Row(
    #     dbc.Col(
    #         html.Table([
    #             html.Thead(html.Tr([html.Th(" "), html.Th("Scan"), html.Th("mut"), html.Th("muf"), html.Th("mur")])),
    #         ], style={"width": "50%"}, id="scan_table"),
    #     ), 
    # ),
    dbc.Row(dcc.Graph(figure=go.Figure(layout={"height": 800}), id="spectrum_plot")),
    ],
    label="Visualization",
)