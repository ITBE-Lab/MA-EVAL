from MA import *
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
import random
import gc
import os
import math
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.layouts import row, column, layout
from bokeh.palettes import d3
from bokeh.models import LinearAxis, Range1d, LinearColorMapper, FixedTicker, BasicTicker, Grid, ColorBar
from bokeh.models.formatters import FuncTickFormatter
from bokeh.layouts import gridplot
from bokeh.io import save
import numpy as np
import colorsys
from statistics import mean, median, stdev
from sys import stderr
from check_accuracy import *
from scipy import stats
from sklearn import linear_model
from subprocess import call
from simulate_pacbio import ReadlengthProvider
import json
import glob
from Bio import SeqIO
from read_simulation import *
from pacBioReads import *
from bam_bai import *
from bokeh_util import *
import glob

def get_ambiguity_distribution(reference, min_len=10, max_len=20, num_queries = 100000):
    def get_all_queries(l):
        if l <= 0:
            yield ""
        else:
            for q in get_all_queries(l-1):
                yield q + "A"
                yield q + "C"
                yield q + "G"
                yield q + "T"

    def get_random_queries(l, amount):
        for _ in range(amount):
            q = ""
            for _ in range(l):
                char = random.randint(1,4)
                if char == 1:
                    q += "A"
                elif char == 2:
                    q += "C"
                elif char == 3:
                    q += "G"
                elif char == 4:
                    q += "T"
            yield q

    def get_random_pos_queries(l, amount, pack):
        if 4**l <= amount: # in this case there are not enough samples...
            print("computing all samples")
            for q in get_all_queries(l):
                yield q
        else: # compute the random unique samples
            print("computing samples from random location")
            already_seen = set()
            already_seen.add("")
            max_pos = pack.unpacked_size() - l - 1
            tries = 0
            for index in range(amount):
                sequence = ""
                while sequence in already_seen:
                    tries += 1
                    pos = random.randint(0, max_pos)
                    while pack.is_bridging(pos, l):
                        pos = random.randint(0, max_pos)
                    sequence = str(pack.extract_from_to(pos, pos+l))
                assert(sequence != "")
                already_seen.add(sequence)
                if (index+1) % 10000 == 0:
                    print(100*(index+1)/amount, "% num_tries:", tries/10000)
                    tries = 0
                yield sequence

    fm_index = FMIndex()
    fm_index.load(reference)

    pack = Pack()
    pack.load(reference)

    r1max = 10
    r2size = 15

    indices1 = []
    indices2 = []
    for i in range(r1max):
        indices1.append(i)
    for i in range(r2size):
        indices2.append(2**i+r1max)
    data1 = []
    data2 = []
    for l in range(min_len, max_len):
        print(l, " of [", min_len, ", ", max_len, ")")
        num_queries_actual = 0
        data1.append( [] )
        for _ in range(r1max):
            data1[-1].append(0.0)
        data2.append( [] )
        for _ in range(r2size):
            data2[-1].append(0.0)
        for q in get_random_pos_queries(l, num_queries, pack):
            ambiguity = fm_index.get_ambiguity(NucSeq(q))
            #in this case an ambiguity of 0 might occur
            if l <= 4**num_queries and ambiguity == 0:
                continue
            #if not fm_index.test_sa_interval(NucSeq(q), pack):
            #    print(q)
            #    print("found error")
            #    exit()
            if ambiguity == 0:
                print(ambiguity)
                print(q)
            if ambiguity < r1max:
                data1[-1][ambiguity] += 1.0
            elif int(math.log2(ambiguity - r1max+1)) < r2size:
                data2[-1][int(math.log2(ambiguity - r1max+1))] += 1.0
            num_queries_actual += 1
        for index, number in enumerate(data1[-1]):
            data1[-1][index] = number / num_queries_actual
        for index, number in enumerate(data2[-1]):
            data2[-1][index] = number / num_queries_actual


    print(indices1)
    print(indices2)
    print("[copy here]")
    print( (min_len, max_len, data1, data2) )#the tuple to be copied
    print("[copy here end]")

    color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=0,
                    high=1
                )

    plot = figure(title="ambiguity on human genome",
            x_range=(0,r1max), y_range=(min_len, max_len),
            x_axis_label='ambiguity', y_axis_label='sequence length',
            plot_width=700, plot_height=500,
            min_border_bottom=10, min_border_top=10,
            min_border_left=10, min_border_right=15,
            tools=["save"]
        )
    plot.image(image=[data1], color_mapper=color_mapper,
            dh=[max_len - min_len], dw=[r1max], x=[0], y=[min_len])

    plot2 = figure(x_range=(r1max,2**r2size+r1max), y_range=(min_len, max_len),
            min_border_bottom=10, min_border_top=10,
            min_border_left=20, min_border_right=15,
            plot_width=500, plot_height=500,tools=[],
            x_axis_type="log"
        )
    plot2.image(image=[data2], color_mapper=color_mapper,
            dh=[max_len - min_len], dw=[2**r2size+r1max], x=[r1max], y=[min_len])

    font = "Helvetica"
    font_size = '15pt'
    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
    color_bar.major_label_text_font=font
    color_bar.major_label_text_font_size=font_size
    plot.add_layout(color_bar, 'left')

    plot.legend.label_text_font=font
    plot.legend.label_text_font_size=font_size
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.axis.axis_label_text_font_size=font_size
    plot.axis.major_label_text_font_size=font_size

    plot2.legend.label_text_font=font
    plot2.legend.label_text_font_size=font_size
    plot2.yaxis.visible = False
    plot2.axis.axis_label_text_font=font
    plot2.axis.major_label_text_font=font
    plot2.axis.axis_label_text_font_size=font_size
    plot2.axis.major_label_text_font_size=font_size

    show(gridplot( [[plot, plot2]] ))

get_ambiguity_distribution("/MAdata/genome/zebrafish")