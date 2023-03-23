import dash
from dash import html, dcc, Input, Output, State, ALL, MATCH
import dash_bootstrap_components as dbc

import plotly.graph_objects as go
from itertools import compress  # basically numpy bool array casting using python iterables

from xas import tiled_io
from xas.tiled_io import filter_node_by_metadata_key, filter_node_for_proposal
from xas.analysis import check_scan

from app_components import build_proposal_accordion, build_filter_input, visualization_tab, normalization_scheme_panel
from app_math import calc_mus, LarchCalculator

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
                dbc.Col(
                    dbc.Button("search", id="search_btn", n_clicks=0, style={"width": "100%"}),
                    width=2,
                    # style={"text-align": "right"},
                    ),
            ]),
            dbc.Row([
                html.Div(id="filters_loc"),
                dbc.Col(
                    dbc.Button("add filter",
                               id="add_filter_btn", 
                               color="link", 
                               size="sm"),
                    width=2,
                ),
            ], align="start",
            ),
            dbc.Row([
                dbc.Col([
                    dbc.Label("Group by"),
                    dcc.Dropdown(
                    options = [
                        {"label": "sample", "value": "sample_name"},
                        {"label": "scan", "value": "monochromator_scan_uid"},
                    ], 
                    value = [
                        "sample_name",
                        "monochromator_scan_uid",
                    ],
                    id="groupby_dropdown",
                    multi=True
                    ),
                ]),
                dbc.Col([
                    dbc.Label("Sort by"),
                    dbc.Input(id="sort_input"),
                ]),
                dbc.Col(
                    dbc.Button("apply", id="apply_btn"),
                    # align="end",
                    width=2,
                ),
            ], align="end", class_name="mb-3"),
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
                    # dbc.Row(dbc.Button("plot", id="plot_btn"), style={"padding-bottom": "10px"}),
                    # dbc.Row(dbc.Button("clear figure", id="clear_btn"), style={"padding-bottom": "10px"}),
                    dbc.Row([
                        dcc.Store(id="xas_normalization_scheme"),
                        normalization_scheme_panel,
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
    State({"type": "filter_key_input", "index": ALL}, "value"),
    State({"type": "filter_value_input", "index": ALL}, "value"),
)
def show_proposal_accordion(n_search_clicks, n_apply_clicks, dropdown_choice, year, cycle, proposal, other_filter_keys, other_filter_values):
    proposal_node = filter_node_for_proposal(ISS_SANDBOX, year, cycle, proposal)
    
    if other_filter_keys and other_filter_values:
        for key, value in zip(other_filter_keys, other_filter_values):
            if key and value:
                proposal_node = filter_node_by_metadata_key(proposal_node, key, value)

    if n_search_clicks == 0:
        return
    if not dropdown_choice:  # check if empty or None
        dropdown_choice = ("sample_name", "monochromator_scan_uid", )

    return build_proposal_accordion(proposal_node, groupby_keys=dropdown_choice)


@app.callback(
    Output("filters_loc", "children"),
    Output({"type": "filter_delete_btn", "index": ALL}, "id"),
    Input("add_filter_btn", "n_clicks"),
    Input({"type": "filter_delete_btn", "index": ALL}, "n_clicks"),
    State({"type": "filter_delete_btn", "index": ALL}, "id"),
    State("filters_loc", "children"),
    prevent_initial_callback=True,
)
def update_filters(add_filter_click, delete_filter_click, current_filter_id_dicts, current_filters):

    updated_id_dicts = current_filter_id_dicts
    updated_filters = current_filters
    
    if dash.ctx.triggered_id == "add_filter_btn":
        if current_filters is None:
            new_filter = build_filter_input(filter_index=0)
            updated_filters = [new_filter]
        else: 
            new_filter = build_filter_input(filter_index=len(current_filters))
            updated_filters.append(new_filter)
    
    # TODO fix index updating
    if isinstance(dash.ctx.triggered_id, dict):
        if dash.ctx.triggered_id["type"] == "filter_delete_btn":
            delete_filter_index = dash.ctx.triggered_id["index"]

            updated_filters.pop(delete_filter_index)
            updated_id_dicts.pop(delete_filter_index)

            for new_index, id_dict in enumerate(updated_id_dicts):
                print(new_index)
                id_dict.update({"index": new_index})
                
    print(updated_id_dicts)
    return updated_filters, updated_id_dicts


@app.callback(
    Output("xas_normalization_scheme", "data"),
    Input("xas_e0_input", "value"),
    Input("xas_pre_edge_start_input", "value"),
    Input("xas_pre_edge_stop_input", "value"),
    Input("xas_post_edge_start_input", "value"),
    Input("xas_post_edge_stop_input", "value"),
    Input("xas_polynom_order_input", "value"),
    prevent_initial_callback=True,
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
    to be stored as json in a `dcc.Store` object"""
    larch_pre_edge_kwargs = dict(
        # step and nvict could be implemented as inputs later
        e0=e0_input,
        step=None,
        pre1=pre_edge_start_input,
        pre2=pre_edge_stop_input,
        norm1=post_edge_start_input,
        norm2=post_edge_stop_input,
        nnorm=post_edge_polynom_order_input,
        nvict=0,  # for some reason this is the only pre_edge keyword that doesn't default to None
    )
    return larch_pre_edge_kwargs


# TODO implement plot undo button using stored previous data
@app.callback(
    Output("spectrum_plot", "figure"),
    Output("previous_plot_data", "data"),
    Input("plot_btn", "n_clicks"),
    Input("clear_btn", "n_clicks"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "value"),
    State({"type": "scan_check", "uid": ALL, "group": ALL}, "id"),
    State("spectrum_plot", "figure"),
    State("previous_plot_data", "data"),
    State("channel_checklist", "value"),

    State("xas_normalization_scheme", "data"),
    State("xas_normalization_radioitems", "value"),
    
    prevent_initial_call=True,
)
def update_plot(
    plot_click,
    clear_click,
    selected_scans,
    selected_scan_id_dicts,
    current_fig, 
    previous_data,
    selected_channels,
    larch_normalization_kwargs,
    xas_normalization_selection,
    ):
    fig = go.Figure(current_fig)
    updated_previous_data = fig.data
    
    if dash.ctx.triggered_id == "clear_btn":
        fig.data = ()

    if dash.ctx.triggered_id == "plot_btn":
        if selected_channels is not None:
            for id_dict in compress(selected_scan_id_dicts, selected_scans):
                uid = id_dict["uid"]
                scan_id = SANDBOX_READER[uid].metadata["scan_id"]
                df = SANDBOX_READER[uid].read()
                calc_mus(df)

                for ch in selected_channels:
                    
                    mu_label = f"{scan_id} {ch}"
                    if xas_normalization_selection == "mu":
                        mu_plot = df[ch]
                    elif xas_normalization_selection == "normalized":
                        mu_plot = LarchCalculator.normalize(df["energy"], df[ch], flatten_output=False, **larch_normalization_kwargs)
                        mu_label += " norm"
                    elif xas_normalization_selection == "flattened":
                        mu_plot = LarchCalculator.normalize(df["energy"], df[ch], flatten_output=True, **larch_normalization_kwargs)
                        mu_label += " flat"

                    # check spectrum isn't already plotted
                    if mu_label not in [trace.name for trace in fig.data]:
                        fig.add_scatter(x=df["energy"], y=mu_plot, name=mu_label)

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
            selected_uids = [id_dict["uid"] for id_dict in compress(scan_id_dicts, selected_scans)]
            selected_scan_df_cols = [set(SANDBOX_READER[uid].read().keys()) for uid in selected_uids]

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
        return tuple(True for _ in dash.ctx.outputs_list)
    else:
        return tuple(False for _ in dash.ctx.outputs_list)


if __name__ == "__main__":
    ISS_SANDBOX = tiled_io.get_iss_sandbox()
    SANDBOX_READER = tiled_io.TiledReader(ISS_SANDBOX)
    app.run_server(debug=True)
