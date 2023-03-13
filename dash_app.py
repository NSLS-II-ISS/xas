import dash
from dash import html, dcc, Input, Output, State, ALL, MATCH
import dash_bootstrap_components as dbc

import numpy as np
import plotly.graph_objects as go

from xas.tiled_io import get_iss_sandbox, filter_node_for_proposal
from xas.analysis import check_scan

from app_components import build_proposal_accordion, visualization_tab

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "new ISS app"


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
            dbc.Row([
                dbc.Col(dbc.Spinner(html.Div(id="accordion_loc"), color="primary")),
                dbc.Col([
                    dbc.Row(
                        dbc.Card([
                            dbc.Checklist(
                                options=[
                                    {"label": "mut", "value": "mut"},
                                    {"label": "muf", "value": "muf"},
                                    {"label": "mur", "value": "mur"},
                                ],
                                id="channel_checklist",
                            ),
                            dbc.Button("see more", color="link", size="sm", id="more_channels_btn"),
                            ],
                            body=True
                        ),
                    style={"padding-bottom": "10px"}),
                    dbc.Row(dbc.Button("plot", id="plot_btn"), style={"padding-bottom": "10px"}),
                    dbc.Row(dbc.Button("clear figure", id="clear_btn"))
                ]),
            ]),
        ], width=3),
        dbc.Col([
            dbc.Tabs([
                visualization_tab,
            ]),
        ], width=9),
    ]),
    dbc.Row(html.Div("test text"))
# ], fluid=True)
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


def calc_mus(df):
    df["mut"] = -np.log(df["it"] / df["i0"])
    df["mur"] = -np.log(df["ir"] / df["it"])
    df["muf"] = df["iff"] / df["i0"]


@app.callback(
    Output("spectrum_plot", "figure"),
    Input("plot_btn", "n_clicks"),
    Input("clear_btn", "n_clicks"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "value"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "id"),
    State("spectrum_plot", "figure"),
    State("channel_checklist", "value"),
    prevent_initial_call=True,
)
def update_plot(plot_click, clear_click, selected_scans, selected_scan_id_dicts, current_fig, selected_channels):
    fig = go.Figure(current_fig)
    
    if dash.ctx.triggered_id == "clear_btn":
        fig.data = ()

    if dash.ctx.triggered_id == "plot_btn":
        if selected_channels is not None:
            for scan_selected, id_dict in zip(selected_scans, selected_scan_id_dicts):
                if scan_selected is True:
                    uid = id_dict["uid"]
                    scan_id = ISS_SANDBOX[uid].metadata["scan_id"]
                    df = ISS_SANDBOX[uid].read()
                    calc_mus(df)
                    for ch in selected_channels:

                        # check spectrum isn't already plotted
                        if f"{scan_id} {ch}" not in [trace.name for trace in fig.data]:
                            fig.add_scatter(x=df["energy"], y=df[ch], name=f"{scan_id} {ch}")
    return fig
        

@app.callback(
    Output({"type": "scan_check", "uid": ALL, "group": MATCH}, "value"),
    Input({"type": "select_all", "group": MATCH}, "value"),
    prevent_initial_call=True,
)
def select_all_scans_in_group(select_all_chk):
    if select_all_chk is True:
        return tuple(True for _ in range(len(dash.ctx.outputs_list)))    
    else:
        return tuple(False for _ in range(len(dash.ctx.outputs_list)))    


if __name__ == "__main__":
    ISS_SANDBOX = get_iss_sandbox()
    app.run_server(debug=True)
