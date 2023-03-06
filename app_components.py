import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

import plotly.graph_objects as go

from xas.tiled_io import sort_node_by_metadata_key

def build_scangroup_interactable(scangroup_node):
    scan_labels = [html.Div([
            html.Div(f"{v.metadata['scan_id']}",
                     style={"display": "inline-block", "padding": "5px"},), 
            dbc.Button("+",
                        id={"type": "plus_btn", "uid": k},
                        color="success", 
                        size="sm"),
                        # style={"background-color": "white"}),
            html.Br(),
        ]) 
        for k, v in scangroup_node.items_indexer
    ]
    return scan_labels


def build_sample_accordion(sample_node):
    scangroup_nodes = sort_node_by_metadata_key(sample_node, "monochromator_scan_uid")
    sample_accordion_items = [
        dbc.AccordionItem(
            build_scangroup_interactable(sg_node),
            title=f"scan group {i+1}"
        )
        for i, sg_node in enumerate(scangroup_nodes)
    ]
    return dbc.Accordion(sample_accordion_items, start_collapsed=True, always_open=True)


def build_scan_accordion(scan_node):
    sample_nodes, sample_labels = sort_node_by_metadata_key(scan_node, "sample_name", return_values=True)
    scan_accordion_items = [
        dbc.AccordionItem(
            build_scangroup_interactable(smp_node),
            title=smp_label
        )
        for smp_node, smp_label in zip(sample_nodes, sample_labels)
    ]
    return dbc.Accordion(scan_accordion_items, start_collapsed=True, always_open=True)


def build_proposal_accordion(proposal_node, sort_key):
    sub_nodes, sub_labels = sort_node_by_metadata_key(proposal_node, sort_key, return_values=True)
    if sort_key == "sample_name":
        build_sub_accordion = build_sample_accordion
    elif sort_key == "monochromator_scan_uid":
        build_sub_accordion = build_scan_accordion
    else:
        raise ValueError("Unsupported sorting key")
    proposal_accordion_items = [
        dbc.AccordionItem(
            build_sub_accordion(sub_node),
            title=sub_name
        )
        for sub_node, sub_name in zip(sub_nodes, sub_labels)
    ]
    return dbc.Accordion(proposal_accordion_items, start_collapsed=True, always_open=True)


visualization_tab = dbc.Tab(
    [
    dbc.Row(
        dbc.Col(
            html.Table([
                html.Thead(html.Tr([html.Th(" "), html.Th("Scan"), html.Th("mut"), html.Th("muf"), html.Th("mur")])),
            ], style={"width": "50%"}, id="scan_table"),
        ), 
    ),
    dbc.Row(dcc.Graph(figure=go.Figure(layout={"height": 800}), id="spectrum_plot")),
    ],
    label="Visualization",
)