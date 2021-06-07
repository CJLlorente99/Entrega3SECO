from abc import abstractproperty
import argparse
import math
from numpy import sign
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import glob
import os
import math
import re

tolerancia = 0.02*math.pi
trmax = 0.3
tsmax = 0.5
Mpmax = 1.15*math.pi
Mpmin = 1.08*math.pi

parser = argparse.ArgumentParser(description='From data to plots. IMPORTANT, files inside inputFolder should have the following format -> WHATEVERX-MOTOR3Y. Accepted values for X are ESC,PAR,RAM or SIN. Accepted values for Y are POS,REF,U and USAT')
parser.add_argument('inputFolder', type=str, help='Path of the folder with csv texts')

args = parser.parse_args()

if not os.path.exists(args.inputFolder + '/PAR'):
    try:
        os.makedirs(args.inputFolder + '/PAR')
    except:
        exit()

if not os.path.exists(args.inputFolder + '/SIN'):
    try:
        os.makedirs(args.inputFolder + '/SIN')
    except:
        exit()

if not os.path.exists(args.inputFolder + '/RAM'):
    try:
        os.makedirs(args.inputFolder + '/RAM')
    except:
        exit()

if not os.path.exists(args.inputFolder + '/ESC'):
    try:
        os.makedirs(args.inputFolder + '/ESC')
    except:
        exit()

inputFiles = glob.glob(args.inputFolder + '/*')
escFiles = []
ramFiles = []
parFiles = []
sinFiles = []

for fileName in inputFiles:
    # filter files

    if re.search('MOTOR',fileName):
        if re.search('ESC', fileName):
            escFiles.append(fileName)
        elif re.search('PAR',fileName):
            parFiles.append(fileName)
        elif re.search('RAM',fileName):
            ramFiles.append(fileName)
        elif re.search('SIN',fileName):
            sinFiles.append(fileName)


esc = {'ESC':escFiles, 'PAR':parFiles, 'RAM':ramFiles, 'SIN':sinFiles}

for signal in esc:
    for file in esc[signal]:

        if re.search('REF',file):
            continue

        if re.search('USAT',file):
            continue

        # deal with position and ref
        if re.search('POS',file):
            dest = args.inputFolder + '/' + signal + '/' + 'respuestaReal' + signal

            df = pd.read_csv(file, skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
            fig = go.Figure()
            fig.add_trace(go.Line(x=df['time'], y=(df['data']), name='Posición'))
            
            df = pd.read_csv(file[:-3]+'REF', skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
            fig.add_trace(go.Line(x=df['time'], y=df['data'], name='Referencia'))

            if signal == 'ESC':
                fig.add_hline(y=Mpmax, line_color='green', line_width=1, line_dash='dash')
                fig.add_hline(y=Mpmin, line_color='green', line_width=1, line_dash='dash')
                fig.add_hline(y=math.pi+tolerancia, line_color='red', line_width=1, line_dash='dash')
                fig.add_hline(y=math.pi-tolerancia, line_color='red', line_width=1, line_dash='dash')
                fig.add_vline(x=trmax, line_width=1, line_dash="dash", line_color="black")
                fig.add_vline(x=tsmax, line_width=1, line_dash="dash", line_color="black")

            fig.write_image(dest + ".pdf")

        if re.search('U',file):
            dest = args.inputFolder + '/' + signal + '/' + 'voltajeReal' + signal

            df = pd.read_csv(file, skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
            fig = go.Figure()
            fig.add_trace(go.Line(x=df['time'], y=(df['data']), name='U'))

            df = pd.read_csv(file+'SAT', skipinitialspace=True, skip_blank_lines=True, delim_whitespace=True,header=None, names=['time','data'])
            fig.add_trace(go.Line(x=df['time'], y=df['data'], name='USAT'))

            fig.update_layout(
                xaxis_title="Tiempo [s]",
                yaxis_title="Voltaje [V]",
                title={
                'text': "voltajeReal" + signal,
                'y':0.9,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'middle'}
            )

            fig.write_image(dest + ".pdf")
