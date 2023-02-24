import dash
from dash import html, dcc, Input, Output, State, ALL
import dash_bootstrap_components as dbc

import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from tiled.queries import Key
from xas.tiled_io import ISS_SANDBOX, filter_node_for_proposal, sort_node_by_metadata_key, get_df_for_uid
from xas.analysis import check_scan

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "new ISS app"


def build_scangroup_interactable(scangroup_node):
    scan_labels_and_checklists = [html.Div([
        html.Div(f"{v.metadata['scan_id']}", id={"type": "scan_id_disp", "uid": k}),
        dcc.Checklist([
                {"label": html.Div("mut"), "value": "mut"},
                {"label": html.Div("muf"), "value": "muf"},
                {"label": html.Div("mur"), "value": "mur"},
            ],
            id={"type": "scan_checklist", "uid": k},
            ),
        html.Br()
    ])
        for k, v in scangroup_node.items_indexer]
    return scan_labels_and_checklists


def build_sample_accordion(sample_node):
    scangroup_nodes = sort_node_by_metadata_key(sample_node, "monochromator_scan_uid")
    sample_accordion_items = [
        dbc.AccordionItem(
            build_scangroup_interactable(sg_node),
            title=f"scan group {i+1}"
        )
        for i, sg_node in enumerate(scangroup_nodes)
    ]
    return dbc.Accordion(sample_accordion_items, start_collapsed=True)


def build_scan_accordion(scan_node):
    sample_nodes, sample_labels = sort_node_by_metadata_key(scan_node, "sample_name", return_values=True)
    scan_accordion_items = [
        dbc.AccordionItem(
            build_scangroup_interactable(smp_node),
            title=smp_label
        )
        for smp_node, smp_label in zip(sample_nodes, sample_labels)
    ]
    return dbc.Accordion(scan_accordion_items, start_collapsed=True)


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
    return dbc.Accordion(proposal_accordion_items)


app.layout = dbc.Container([
    html.H1("XDash",
            style={
                "textAlign": "center",
                "font-size": "400%",
                }),
    dbc.Row([
        dbc.Col([
            dbc.Row(dbc.Col(html.Div("Search by proposal"))),
            dbc.Row([
                dbc.Col(dbc.Input(id="year_input", placeholder="year")),
                dbc.Col(dbc.Input(id="cycle_input", placeholder="cycle")),
                dbc.Col(dbc.Input(id="proposal_input", placeholder="proposal")),
                dbc.Col(dbc.Button("search", id="search_btn", n_clicks=0))
            ], style={"padding-bottom": "10px"}),
            dbc.Row([
                    dbc.Col(dcc.Dropdown([
                        {"label": "sample", "value": "sample_name"},
                        {"label": "scan plan", "value": "monochromator_scan_uid"},
                    ], placeholder="Sort by...", id="sort_dropdown"
                    )),
                    dbc.Col(dbc.Button("sort", id="sort_btn")),
            ], style={"padding-bottom": "10px"}),
            dbc.Row(dbc.Col(id="accordion_loc")),
        ], width=3),
        dbc.Col([
            dbc.Row(dcc.Graph(id="spectrum_plot")),
            dbc.Row(dbc.Col(
                dbc.Button("plot selected spectra", id="plot_btn", class_name="d-grid gap-2 col-3 mx-auto"),
            )),
        ]),
    ]),
], fluid=True)


@app.callback(
    Output("accordion_loc", "children"),
    Input("search_btn", "n_clicks"),
    Input("sort_btn", "n_clicks"),
    State("sort_dropdown", "value"),
    State("year_input", "value"),
    State("cycle_input", "value"),
    State("proposal_input", "value"),
)
def show_proposal_accordion(n_search_clicks, n_sort_clicks, dropdown_choice, year, cycle, proposal):
    if n_search_clicks == 0:
        return
    if dropdown_choice is None:
        dropdown_choice = "sample_name" 
    return build_proposal_accordion(filter_node_for_proposal(ISS_SANDBOX, year, cycle, proposal), sort_key=dropdown_choice)


@app.callback(
    Output("spectrum_plot", "figure"),
    Input("plot_btn", "n_clicks"),
    State({"type": "scan_checklist", "uid": ALL}, "value"),
    State({"type": "scan_checklist", "uid": ALL}, "id"),
)
def plot_selected_spectra(click, checklist_values, checklist_id_dicts):
    fig = go.Figure(layout={"height": 800})
    for cv, cd in zip(checklist_values, checklist_id_dicts):
        if cv is None:
            continue
        uid = cd["uid"]
        df = ISS_SANDBOX[uid].read()
        df["mut"] = -np.log(df["it"] / df["i0"])
        df["mur"] = -np.log(df["ir"] / df["it"])
        df["muf"] = df["iff"] / df["i0"]
        for val in cv:
            fig.add_scatter(x=df["energy"], y=df[val])
    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
