import pandas as pd
import plotly.express as px  # (version 4.7.0 or higher)
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output  # pip install dash (version 2.0.0 or higher)
import json
import os

app = Dash(__name__, requests_pathname_prefix='/siamang/', assets_folder='/data/www/wimpouw/web/siamang/assets/')

#server = app.server
# get the directory path where this script is located
dir_path = os.path.dirname(os.path.abspath(__file__))

# set the current working directory to this directory
os.chdir(dir_path)
#print(dir_path)
# -- Import and clean data (importing csv into pandas)
df = pd.read_csv("assets/data.csv")

# ------------------------------------------------------------------------------
# App layout
app.layout = html.Div([

    html.H1("Multimodal Vocalizations Siamang", style={'text-align': 'center'}),
    html.H2("Pouw, Kehy, Gamba, Ravignani", style={'text-align': 'center'}),
    dcc.Dropdown(id="mode_rep",
                 options=[
                     {"label": "envelope", "value": "envelope"}],#,{"label": "F0", "value": "F0"}
                 multi=False,
                 value="envelope",
                 style={'width': "40%"}
                 ),

    html.Div(id='output_container', children=[], style = {'color': 'white'}),
    html.Br(),
    dcc.Graph(id='MY_XY_Map', figure={},style={'width': '40%', 'textAlign': 'center','display': 'inline-block'}),
    html.Video(controls=True, id='videoplayer', src='', style={'width': '30%', 'textAlign': 'center', 'display': 'inline-block', 'vertical-align': 'top'}, autoPlay=True),
    html.Video(controls=True, id='videoplayer2', src='', style={'width': '20%', 'textAlign': 'center', 'display': 'inline-block', 'vertical-align': 'top'}, autoPlay=True)
])

# ------------------------------------------------------------------------------
# Connect the Plotly graphs with Dash Components
@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='MY_XY_Map', component_property='figure'),
     Output(component_id ='videoplayer', component_property='src'),
    Output(component_id ='videoplayer2', component_property='src')],
    [Input(component_id='mode_rep', component_property='value'), Input('MY_XY_Map', 'clickData')]
)

def update_graph(option_slctd, clickData):
    container = "Click on any point to inspect the multimodal event associated with it.  \
                \n We only show one camera angle for this visualization, and the motion tracking performance for that camera angle. \
                \n For some data points there will be no motion tracking video because the individual was not visible or the tracking was poor. Further note that motion tracking was applied on bounding boxed portion of the original camera (see methods).  \
                \n You subselected events for: {}".format(option_slctd)
    dff = df.copy()
    #select variables
    Y = dff['peak_amp']
    if option_slctd=="F0":
        Y=dff['peak_F0']
    if option_slctd=="envelope":
        Y=dff['peak_amp']
    # Plotly Express
    fig = px.scatter(
        data_frame=dff,
        x=dff['peak_acc'],
        y= Y,
        color = dff['individual'],
        opacity=0.75,
        #color='Pct of Colonies Impacted',
        hover_data=['vidnames'],
        color_continuous_scale=px.colors.sequential.YlOrRd,
        template='plotly_dark', trendline="ols", trendline_scope="overall",
        labels={
            "peak_acc": "Body-scaled magnitude of the global peak in thorax acceleration<br>(peak acc)",
            "peak_amp": "The nearest maximum amplitude envelope peak<br>(nearest peak envelope)",
            "peak_F0": "The nearest maximum F0 peak<br>(nearest peak F0)",
            "individual": "Individual"
        })
    fig.update_traces(marker_size=15)

    #get the video clicked on
    check = str(clickData)
    converted_to_legal_json = check.replace("'", '"')
    test = eval(converted_to_legal_json)
    src = ''
    src2 = ''
    if test is not None:
        test = eval(str(converted_to_legal_json))
        test = eval(str(test['points']))
        test = eval(str(test[0]))
        vidname = 'gopro_orange_opp_'+test['customdata'][0]
        src = 'assets/' + vidname + '.mp4'
        #also add in the second tracked video
        src2 = 'snip_' + vidname + '.mp4'
        src2 = 'assets/' + src2
        src2 = src2.replace(".mp4", "DLC_resnet50_SiamangSimpleV2Jan31shuffle1_500000_filtered_labeled-converted.mp4")

    return container, fig, src, src2

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    app.run_server(debug=True)
