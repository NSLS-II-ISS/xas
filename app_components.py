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


def build_nested_accordion(base_node, groupby_keys: list[str], _node_label=""):
    current_key = groupby_keys[0]
    next_nodes, next_labels = group_node_by_metadata_key(base_node, current_key, return_values=True)
    next_level_keys = groupby_keys[1:]

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

    return dbc.Accordion(accordion_items, start_collapsed=True, always_open=True,)


def build_proposal_accordion(proposal_node, groupby_keys):
    proposal_accordion = build_nested_accordion(proposal_node, groupby_keys)
    return html.Div(proposal_accordion, style={"max-height": "700px", "overflow-y": "scroll"})


# class FilterInput:
#     def __init__(self, filters=None) -> None:
#         if filters is not None:
#             self.filters = filters

#     def build_new_input




def build_filter_input(filter_index):
    key_input = dbc.Input(id={"type": "filter_key_input", "index": filter_index},
                          placeholder="metadata key")
    value_input = dbc.Input(id={"type": "filter_value_input", "index": filter_index},
                            placeholder="value")
    delete_button = dbc.Button("X", 
                               id={"type": "filter_delete_btn", "index": filter_index},
                               color="light",)
    
    key_value_inputgroup = dbc.InputGroup([
        key_input,
        dbc.InputGroupText(":"),
        value_input,
    ])
    
    return dbc.Row([
        dbc.Col(key_value_inputgroup),
        dbc.Col(delete_button, width=1),
    ],)
    # ], id={"type": "filter_input_row", "index": filter_index})


visualization_tab = dbc.Tab([
        dcc.Store(id="previous_plot_data"),
        dbc.Row(dcc.Graph(figure=go.Figure(layout={"height": 800}), id="spectrum_plot")),
        dbc.Row([
            dbc.Col(
                dbc.Button("plot", id="plot_btn", style={"width": "100%"}),
                width=4,
            ),
            dbc.Col(
                dbc.Button("clear figure", id="clear_btn", style={"width": "100%"}),
                width=4
            )
        ], justify="center")
    ],
    label="Visualization",
)

normalization_scheme_panel = dbc.Card([
    html.Div("XAS Normalization Parameters", className="mb-3"),
    html.Div([
        dbc.InputGroup([
            dbc.InputGroupText(["E", html.Sub("0")]),
            dbc.Input(id="xas_e0_input", type="number"),
            dbc.InputGroupText("[eV]"),
        ]),
        html.Div("Pre-edge range"),
        dbc.InputGroup([
            dbc.Input(id="xas_pre_edge_start_input", type="number"),
            dbc.InputGroupText("⮕"),
            dbc.Input(id="xas_pre_edge_stop_input", type="number"),
            dbc.InputGroupText("[eV]"),
        ]),
        html.Div("Post-edge range"),
        dbc.InputGroup([
            dbc.Input(id="xas_post_edge_start_input", type="number"),
            dbc.InputGroupText("⮕"),
            dbc.Input(id="xas_post_edge_stop_input", type="number"),
            dbc.InputGroupText("[eV]"),
        ], class_name="mb-2"),
        dbc.InputGroup([
            dbc.InputGroupText("Polynom order"),
            dbc.Input(id="xas_polynom_order_input", type="number"),
        ]),
    ], style={"padding-bottom": "8px"}),
    html.Div([
        dbc.RadioItems(
            options=[
                {"label": "mu", "value": "mu"},
                {"label": "normalized", "value": "normalized"},
                {"label": "flattened", "value": "flattened"},
            ],
            value="mu",
            id="xas_normalization_radioitems",
        ),
    ]),
],
body=True,
id="norm_scheme_panel")

