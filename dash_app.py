import dash
from dash import html, dcc, Input, Output, State, ALL, MATCH
import dash_bootstrap_components as dbc

import numpy as np
import larch
from larch.xafs import pre_edge
import plotly.graph_objects as go
import itertools

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
                dbc.Col(
                    dcc.Dropdown(
                    options = [
                        {"label": "sample", "value": "sample_name"},
                        {"label": "scan", "value": "monochromator_scan_uid"},
                    ], 
                    value = [
                        "sample_name",
                        "monochromator_scan_uid",
                    ],
                    placeholder="Group by...", 
                    id="groupby_dropdown",
                    multi=True
                )),
                dbc.Col(dbc.Button("apply", id="apply_btn")),
            ]),
            dbc.Row([
                dbc.Col(
                    dbc.Button("add filter", color="link", size="sm"),
                    width=2
                ),
                html.Div(id="filters_loc")
            ], align="start",
            ),
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
                            dbc.Button("see more", 
                                        color="link",
                                        size="sm",
                                        n_clicks=0,
                                        id="change_visible_channels_btn"),
                            ],
                            body=True
                        ),
                    style={"padding-bottom": "10px"}),
                    dbc.Row(dbc.Button("plot", id="plot_btn"), style={"padding-bottom": "10px"}),
                    dbc.Row(dbc.Button("clear figure", id="clear_btn"), style={"padding-bottom": "10px"}),
                    dbc.Row([
                        dcc.Store(id="xas_normalization_scheme"),
                        dbc.Card([
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
                            ]),
                            html.Div([
                                dbc.Checklist(
                                    options= {
                                        "label": "Normalized", "value": "Normalized",
                                        "label": "Flattened", "value": "Flattened",
                                    },
                                    id="xas_normalization_checklist"
                                )
                            ])

                        ],
                        body=True,
                        id="norm_scheme_panel")
                    ])
                ], style={"max-height": "700px", "overflow-y": "auto"}),
            ]),
        ], width=4),
        dbc.Col([
            dbc.Tabs([
                visualization_tab,
            ]),
        ], width=8),
    ],
    style={"max-height": "800px", "overflow-y": "visible"}),
    dbc.Row(html.Div("test text"))
], fluid=True)


@app.callback(
    Output("accordion_loc", "children"),
    Input("search_btn", "n_clicks"),
    Input("apply_btn", "n_clicks"),
    State("groupby_dropdown", "value"),
    State("year_input", "value"),
    State("cycle_input", "value"),
    State("proposal_input", "value"),
)
def show_proposal_accordion(n_search_clicks, n_apply_clicks, dropdown_choice, year, cycle, proposal):
    if n_search_clicks == 0:
        return
    if not dropdown_choice:  # check if empty or None
        dropdown_choice = ("sample_name", "monochromator_scan_uid", )
    print(dropdown_choice)
    return build_proposal_accordion(filter_node_for_proposal(ISS_SANDBOX, year, cycle, proposal), groupby_keys=dropdown_choice)


def calc_mus(df):
    """Shorthand function to calculate mut, muf, and mur 'on the fly'"""
    df["mut"] = -np.log(df["it"] / df["i0"])
    df["mur"] = -np.log(df["ir"] / df["it"])
    df["muf"] = df["iff"] / df["i0"]


@app.callback(
    Output("xas_normalization_scheme", "data"),
    Input("xas_e0_input", "value"),
    Input("xas_pre_edge_start_input", "value"),
    Input("xas_pre_edge_stop_input", "value"),
    Input("xas_post_edge_start_input", "value"),
    Input("xas_post_edge_start_input", "value"),
    Input("xas_polynom_order_input", "value"),
)
def update_normalization_scheme(
    e0_input,
    pre_edge_start_input, 
    pre_edge_stop_input,
    post_edge_start_input,
    post_edge_stop_input,
    post_edge_polynom_order_input,
    ):
    """Returns dict of `larch.xafs.pre_edge` keyword-argument pairs
    to be stored as json in a `dcc.Store` object."""
    print(
        dict(
            # step and nvict could be implemented as inputs later
            e0=e0_input,
            step=None,
            pre1=pre_edge_start_input,
            pre2=pre_edge_stop_input,
            norm1=post_edge_start_input,
            norm2=post_edge_stop_input,
            nnorm=post_edge_polynom_order_input,
            nvict=None,
        )
    )
    return dict(
        # step and nvict could be implemented as inputs later
        e0=e0_input,
        step=None,
        pre1=pre_edge_start_input,
        pre2=pre_edge_stop_input,
        norm1=post_edge_start_input,
        norm2=post_edge_stop_input,
        nnorm=post_edge_polynom_order_input,
        nvict=None,
    )




@app.callback(
    Output("spectrum_plot", "figure"),
    Output("previous_plot_data", "data"),
    Input("plot_btn", "n_clicks"),
    Input("clear_btn", "n_clicks"),
    # Input("norm_view_toggle", "value"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "value"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "id"),
    State("spectrum_plot", "figure"),
    State("previous_plot_data", "data"),
    State("channel_checklist", "value"),

    State("xas_normalization_checklist", "value"),
    
    
    prevent_initial_call=True,
)
# def update_plot(plot_click, clear_click, normalized_view, selected_scans, selected_scan_id_dicts, current_fig, previous_data, selected_channels):
def update_plot(plot_click, clear_click, selected_scans, selected_scan_id_dicts, current_fig, previous_data, selected_channels, normalization_selections):
    fig = go.Figure(current_fig)
    updated_previous_data = fig.data
    
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
                        # if "Normalization" in normalization_selections:

                        # check spectrum isn't already plotted
                        if f"{scan_id} {ch}" not in [trace.name for trace in fig.data]:
                            fig.add_scatter(x=df["energy"], y=df[ch], name=f"{scan_id} {ch}")

    if dash.ctx.triggered_id == "norm_view_toggle":
        if normalized_view is True:
            current_data = fig.data
            fig.data = ()
            for trace in current_data:
                raw_larch_group = larch.Group(energy=np.array(trace.x), mu=np.array(trace.y))
                norm_larch_group = larch.Group()
                pre_edge(raw_larch_group, group=norm_larch_group)
                fig.add_scatter(x=raw_larch_group.energy, y=norm_larch_group.flat, name=f"{trace.name} norm")
        else:
            # if previous_data is None:
            #     raise dash.exceptions.PreventUpdate
            # else:
            fig.data = ()
            for trace_data in previous_data:
                fig.add_scatter(**trace_data)

    return fig, updated_previous_data
        

@app.callback(
    Output("channel_checklist", "options"),
    Output("change_visible_channels_btn", "children"),
    Input("change_visible_channels_btn", "n_clicks"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "value"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "id"),
    State("change_visible_channels_btn", "children"),
    prevent_initial_call=True,
)
def change_visible_channels(n_channel_clicks, selected_scans, scan_id_dicts, current_btn_text):
    default_options = [
        {"label": "mut", "value": "mut"},
        {"label": "muf", "value": "muf"},
        {"label": "mur", "value": "mur"},
    ]

    if current_btn_text == "see more":
        if any(selected for selected in selected_scans):
            selected_uids = [id_dict["uid"] for id_dict in itertools.compress(scan_id_dicts, selected_scans)]
            selected_scan_df_cols = [set(ISS_SANDBOX[uid].read().keys()) for uid in selected_uids]

            # flatten into set of all unique column names
            other_channels = set.union(*selected_scan_df_cols)
            
            channel_btn_text = "see less"
            new_options = [{"label": ch, "value": ch} for ch in sorted(other_channels)]

        else:
            new_options = []

        channel_options = default_options + new_options
    
    else:
        channel_options = default_options
        channel_btn_text = "see more"
    
    return channel_options, channel_btn_text


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
