import dash
from dash import html, dcc, Input, Output, State, ALL, MATCH
import dash_bootstrap_components as dbc

import numpy as np
import plotly.graph_objects as go

from tiled.queries import Key
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
            dbc.Row(
                dbc.Col(dbc.Spinner(html.Div(id="accordion_loc"), color="primary"))
            ),
        ], width=3),
        dbc.Col([
            dbc.Tabs([
                visualization_tab,
            ]),
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
    Output("scan_table", "children"),
    Input({"type": "plus_btn", "uid": ALL}, "n_clicks"),
    State("scan_table", "children"),
    State({"type": "table_row", "uid": ALL}, "id"),
    prevent_initial_call=True,
)
def add_scan_to_table(click, table_rows, row_id_dicts):
    uid = dash.ctx.triggered_id["uid"]
    scan_id = ISS_SANDBOX[uid].metadata["scan_id"]
    new_row = html.Tr([html.Td(scan_id),
                     html.Td(html.Div(dbc.Checkbox(value=False, id={"type": "mut_check", "uid": uid}))),
                     html.Td(html.Div(dbc.Checkbox(value=False, id={"type": "muf_check", "uid": uid}))),
                     html.Td(html.Div(dbc.Checkbox(value=False, id={"type": "mur_check", "uid": uid}))),
                     ], id={"type": "table_row", "uid": uid})
    if uid not in [d["uid"] for d in row_id_dicts]:
        table_rows.append(new_row)
    return table_rows


def remove_trace_from_fig(name_to_remove, fig: go.Figure):
    current_traces = fig.data
    fig.data = tuple(trace for trace in current_traces if trace.name != name_to_remove)


@app.callback(
    Output("spectrum_plot", "figure"),
    Input({"type": "mut_check", "uid": ALL}, "value"),
    Input({"type": "muf_check", "uid": ALL}, "value"),
    Input({"type": "mur_check", "uid": ALL}, "value"),
    State("spectrum_plot", "figure"),
    prevent_initial_call=True,   
)
def add_spectrum_to_plot(mut_chk, muf_chk, mur_chk,
                         current_fig):
    print(dash.ctx.triggered_id)
    uid = dash.ctx.triggered_id["uid"]
    scan_id = ISS_SANDBOX[uid].metadata["scan_id"]
    fig = go.Figure(current_fig)
    
    if (dash.ctx.triggered_id["type"] == "mut_check"):
        if dash.ctx.triggered[0]["value"] is True:
            df = ISS_SANDBOX[uid].read()
            df["mut"] = -np.log(df["it"] / df["i0"])
            fig.add_scatter(x=df["energy"], y=df["mut"], name=f"{scan_id} mut")
        if dash.ctx.triggered[0]["value"] is False:
            remove_trace_from_fig(f"{scan_id} mut", fig)

    if (dash.ctx.triggered_id["type"] == "muf_check"):
        if  dash.ctx.triggered[0]["value"] is True:
            df = ISS_SANDBOX[uid].read()
            df["muf"] = df["iff"] / df["i0"]
            fig.add_scatter(x=df["energy"], y=df["muf"], name=f"{scan_id} muf")
        if dash.ctx.triggered[0]["value"] is False:
            remove_trace_from_fig(f"{scan_id} muf", fig)

    if (dash.ctx.triggered_id["type"] == "mur_check"):
        if  dash.ctx.triggered[0]["value"] is True:
            df = ISS_SANDBOX[uid].read()
            df["mur"] = -np.log(df["ir"] / df["it"])
            fig.add_scatter(x=df["energy"], y=df["mur"], name=f"{scan_id} mur")
        if dash.ctx.triggered[0]["value"] is False:
            remove_trace_from_fig(f"{scan_id} mur", fig)
    
    return fig

# @app.callback(
#     Output("spectrum_plot", "figure"),
#     Input({"type": "plus_btn", "uid": ALL}, "n_clicks"),
#     State("spectrum_plot", "figure"),
#     prevent_initial_call=True,
# )
# def add_spectrum_to_plot(click, current_fig):
#     uid = dash.ctx.triggered_id["uid"]
#     fig = go.Figure(current_fig)

#     scan_id = ISS_SANDBOX[uid].metadata["scan_id"]
#     df = ISS_SANDBOX[uid].read()
#     df["mut"] = -np.log(df["it"] / df["i0"])
#     df["mur"] = -np.log(df["ir"] / df["it"])
#     df["muf"] = df["iff"] / df["i0"]
    
#     for mu in ("mut", "muf", "mur"):
#         fig.add_scatter(x=df["energy"], y=df[mu], name=f"{scan_id} {mu}")
    
#     return fig


if __name__ == "__main__":
    ISS_SANDBOX = get_iss_sandbox()
    app.run_server(debug=True)
