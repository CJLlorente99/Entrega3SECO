from abc import abstractproperty
import argparse
import math
from numpy import sign
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import math

df = pd.read_csv("test1ESC1-MOTOR3POS", skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
fig = go.Figure()
fig.add_trace(go.Line(x=df['time'], y=(df['data']), name='Posición 1ms'))

df = pd.read_csv("test1ESC5-MOTOR3POS", skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
fig.add_trace(go.Line(x=df['time'], y=df['data'], name='Posición 5ms'))

df = pd.read_csv("test1ESC10-MOTOR3POS", skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
fig.add_trace(go.Line(x=df['time'], y=df['data'], name='Posición 10ms'))

fig.update_layout(
    xaxis_title="Tiempo [s]",
    yaxis_title="Posición [V]",
    title={
    'text': "Comparacion con diferentes periodos",
    'y':0.9,
    'x':0.5,
    'xanchor': 'center',
    'yanchor': 'middle'}
)

fig.write_image("comparation1msvs5ms.pdf")
