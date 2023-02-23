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


def build_proposal_accordion(proposal_node):
    sample_nodes, sample_names = sort_node_by_metadata_key(proposal_node, "sample_name", return_values=True)
    proposal_accordion_items = [
        dbc.AccordionItem(
            build_sample_accordion(smp_node),
            title=smp_name
        )
        for smp_node, smp_name in zip(sample_nodes, sample_names)
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
                dbc.Col(dbc.Button("search", id="search_btn"))
            ]),
            html.Br(),
            dbc.Row(dbc.Col(dbc.Button("check scans", color="success", id="check_btn"))),
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
    State("year_input", "value"),
    State("cycle_input", "value"),
    State("proposal_input", "value"),
)
def search_click(click, year, cycle, proposal):
    return build_proposal_accordion(filter_node_for_proposal(ISS_SANDBOX, year, cycle, proposal))


@app.callback(
    Output("spectrum_plot", "figure"),
    Input("plot_btn", "n_clicks"),
    State({"type": "scan_checklist", "uid": ALL}, "value"),
    State({"type": "scan_checklist", "uid": ALL}, "id"),
)
def plot_selected_spectra(click, checklist_values, checklist_id_dicts):
    fig = go.Figure()
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

# not working
@app.callback(
    Output({"type": "scan_checklist", "uid": ALL}, "children"),
    Input("check_btn", "n_clicks"),
    State({"type": "scan_checklist", "uid": ALL}, "id")
)
def run_check_scans(click, checklist_id_dicts):
    print("pp")
    return [[
        {"label": html.Div(",", style={"color": "red"}), "value": "mut"},
        {"label": html.Div("muf", style={"color": "red"}), "value": "muf"},
        {"label": html.Div("mur", style={"color": "red"}), "value": "mur"},
    ]] * len(checklist_id_dicts)


if __name__ == "__main__":
    app.run_server(debug=True)
