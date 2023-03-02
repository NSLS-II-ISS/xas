import dash
from dash import html, dcc, Input, Output, State, ALL
import dash_bootstrap_components as dbc

import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from tiled.queries import Key
from xas.tiled_io import get_iss_sandbox, filter_node_for_proposal
from xas.analysis import check_scan

from app_components import build_proposal_accordion

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "new ISS app"


app.layout = dbc.Container([
    dcc.Store(id="bnl_username"),
    dbc.Modal(
        dbc.ModalHeader(""),
        dbc.ModalBody(dbc.Input(id="username_input")),
    ),  
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
            dbc.Row(
                dbc.Col(dbc.Spinner(html.Div(id="accordion_loc"), color="primary"))
            ),
        ], width=3),
        dbc.Col([
            dbc.Row(dcc.Graph(figure=go.Figure(layout={"height": 800}), id="spectrum_plot")),
            dbc.Row(dbc.Col(
                dbc.Button("plot selected spectra", id="plot_btn", class_name="d-grid gap-2 col-3 mx-auto"),
            )),
        ], width=9),
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
    Input({"type": "plus_btn", "uid": ALL}, "n_clicks"),
    State("spectrum_plot", "figure"),
    prevent_initial_call=True,
)
def add_spectrum_to_plot(click, current_fig):
    uid = dash.ctx.triggered_id["uid"]
    print(dash.ctx.triggered_id)
    fig = go.Figure(current_fig)

    df = ISS_SANDBOX[uid].read()
    df["mut"] = -np.log(df["it"] / df["i0"])
    df["mur"] = -np.log(df["ir"] / df["it"])
    df["muf"] = df["iff"] / df["i0"]
    
    for mu in ("mut", "muf", "mur"):
        fig.add_scatter(x=df["energy"], y=df[mu])
    
    return fig


# @app.callback(
#     Output("spectrum_plot", "figure"),
#     Input("plot_btn", "n_clicks"),
#     State({"type": "scan_checklist", "uid": ALL}, "value"),
#     State({"type": "scan_checklist", "uid": ALL}, "id"),
# )
# def plot_selected_spectra(click, checklist_values, checklist_id_dicts):
#     fig = go.Figure(layout={"height": 800})
#     for cv, cd in zip(checklist_values, checklist_id_dicts):
#         if cv is None:
#             continue
#         uid = cd["uid"]
        # df = ISS_SANDBOX[uid].read()
        # df["mut"] = -np.log(df["it"] / df["i0"])
        # df["mur"] = -np.log(df["ir"] / df["it"])
        # df["muf"] = df["iff"] / df["i0"]
#         for val in cv:
#             fig.add_scatter(x=df["energy"], y=df[val])
#     return fig


if __name__ == "__main__":
    ISS_SANDBOX = get_iss_sandbox()
    app.run_server(debug=True)
