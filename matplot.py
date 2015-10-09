from PySide import QtCore, QtGui
import numpy as np
import math
import os
import csv
from matplotlibwidget import maplot, scatter
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from Bio import SeqIO
import hdf5_gen


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)



items3prime = 'MazF none,MazF PNK,MqsR none,MqsR PNK,Log none,Log PNK,Stat none,Stat PNK,stat_exo- none,stat_exo- PNK'#'MazF PNK-,MazF PNK+,MqsR PNK-,MqsR PNK+,Log PNK-,Log PNK+,Stat PNK-,Stat PNK+,stat_exo- PNK-,stat_exo- PNK+'
items5prime = 'MazF none,MazF PNK,MazF TAP,MqsR none,MqsR PNK,MqsR TAP,Log none,Log PNK,Log TAP,Stat none,Stat PNK,Stat TAP,stat_exo- none,stat_exo- PNK,stat_exo- TAP'

# Setting up the whole interface
class Plotter(QtGui.QWidget):
    def setupUi(self):

        self.plotButton = QtGui.QPushButton('Plot', self)
        QtCore.QObject.connect(self.plotButton, QtCore.SIGNAL('clicked()'), self.PlotFunc)

        self.widget1 = maplot(self)
        self.toolbar1 = NavigationToolbar(self.widget1.canvas, self)
        self.widget1.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))

        self.widget2 = scatter(self)
        self.widget2.setObjectName(_fromUtf8("Scatter"))
        self.toolbar2 = NavigationToolbar(self.widget2.canvas, self)
        self.widget2.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))

        self.proc1combo = QtGui.QComboBox(self)
        self.proc1combo.addItems(items5prime.split(','))

        self.proc2combo = QtGui.QComboBox(self)
        self.proc2combo.addItems(items5prime.split(','))

        self.primcombo = QtGui.QComboBox(self)
        self.primcombo.addItems('5prime,3prime'.split(','))
        self.primcombo.currentIndexChanged[str].connect(self.on_combo_prime_change)

        self.subcombo = QtGui.QComboBox(self)
        self.subcombo.addItems('16S,23S'.split(','))

        self.valcombo = QtGui.QComboBox(self)
        self.valcombo.addItems('Relative,Absolute values'.split(','))

        self.text = QtGui.QTextEdit(self)
        self.text.setMaximumWidth(400)
        self.clearButton = QtGui.QPushButton('clear', self)
        QtCore.QObject.connect(self.clearButton, QtCore.SIGNAL('clicked()'), self.text.clear)


        pbox = QtGui.QVBoxLayout()
        pbox.setSpacing(1)
        pbox.addWidget(self.toolbar1)
        pbox.addWidget(self.widget1)
        pbox.addWidget(self.toolbar2)
        pbox.addWidget(self.widget2)

        tbbox = QtGui.QHBoxLayout()
        tbbox.addWidget(self.plotButton)
        tbbox.addWidget(self.clearButton)

        cbbox = QtGui.QVBoxLayout()
        cbbox.addWidget(self.proc1combo)
        cbbox.addWidget(self.proc2combo)
        cbbox.addWidget(self.primcombo)
        cbbox.addWidget(self.subcombo)
        cbbox.addWidget(self.valcombo)

        bbox = QtGui.QVBoxLayout()
        bbox.addLayout(cbbox)
        bbox.addWidget(self.text)
        bbox.addLayout(tbbox)

        abox = QtGui.QHBoxLayout()
        abox.addLayout(bbox)
        abox.addLayout(pbox)

        self.setLayout(abox)

    data_dic = {}

    #=====================================================================================================
    # Function that retrieves data from HDF according to the selected parameters and draws requested plots
    def PlotFunc(self):

        global data_dic

        subunit = ['y_pos_16S', 'y_pos_23S']
        subunit_MA = ['MA_y_pos_16S', 'MA_y_pos_23S']
        subunitc = ['colour_16S', 'colour_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        plotname = ['16S Scatterplot', '23S Scatterplot']

        #================================================================================
        # Name conversion to compensate for different labels in the cboxes and in the hdf
        name_conv = {'MazF': 'MazF2h',
                    'MqsR': 'MqsR2h',
                    'Log': 'MG1655log',
                    'Stat': 'MG1655stats',
                    'stat_exo-': 'delta3',
                    'none': 'PNK-',
                    'PNK': 'PNK+'}

        proc1 = self.proc1combo.currentText()
        proc2 = self.proc2combo.currentText()
        prim_select = self.primcombo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')
        prim_select = str(prim_select)

        #===============================================================
        # Retrieving data from HDF. Changing colors for both processings
        if prim_select == '5prime':
            data_dic = hdf5_gen.get_hdf_data([prim_select, prim_select], [name_conv[proc1[0]], name_conv[proc2[0]]],[proc1[1], proc2[1]],['_input_', 'some'])
        else:
            data_dic = hdf5_gen.get_hdf_data([prim_select, prim_select], [name_conv[proc1[0]], name_conv[proc2[0]]],[name_conv[proc1[1]], name_conv[proc2[1]]],['_input_', 'some'])
        colors = []
        for i in range(len(data_dic["data1"][subunitc[sub_select]])):
            c = data_dic["data1"][subunitc[sub_select]][i]
            if c == 'b':
                colors.append('black')
            else:
                colors.append(c)

        #===============================================================
        # Drawing the MA plot
        MA_X_16S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1'][subunit_MA[sub_select]], data_dic['data2'][subunit_MA[sub_select]])]
        MA_Y_16S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1'][subunit_MA[sub_select]], data_dic['data2'][subunit_MA[sub_select]])]
        self.widget1.canvas.ax.clear()
        self.widget1.canvas.ax.set_ylabel('M', fontsize=20)
        self.widget1.canvas.ax.set_xlabel('A', fontsize=20)
        series1 = self.widget1.canvas.ax.scatter(MA_X_16S, MA_Y_16S, alpha=0.5, c=colors, linewidths=( 0, 0, 0), picker=True, label='Datapoints')
        self.widget1.canvas.ax.set_ylim(min(MA_Y_16S)-2,max(MA_Y_16S) + 2)
        mazf = self.widget1.canvas.ax.scatter(0,min(MA_Y_16S)-200, alpha=0.5, c='red', marker = 'o', label = ' _ACA')
        mqsr = self.widget1.canvas.ax.scatter(0,min(MA_Y_16S)-200, alpha=0.5, c='cyan', marker = 'o', label = 'G_CB')
        self.widget1.canvas.ax.legend(handles=[series1,mazf,mqsr],loc='best', scatterpoints = 1)
        self.widget1.canvas.draw()

        #===============================================================
        # Drawing the erlative and absolute graphs. Value_select 0 - relative, 1 - absolute
        if value_select == 0:
            self.widget2.canvas.ax.clear()
            hundred_1 = 0
            hundred_2 = 0
            #===============================================================
            # Selection of locations for 100% for 5-prim and 3-prim

            if prim_select == '5prime':
                if sub_select == 0:
                    area=range(-6,7)
                else:
                    area=range(0,14)
                for i in area:
                    index = data_dic['data1'][nucl_data[sub_select]].index(float(i))
                    hundred_1 += data_dic['data1'][subunit[sub_select]][index]
                    hundred_2 += data_dic['data2'][subunit[sub_select]][index]
            elif prim_select == '3prime':
                if sub_select == 0:
                    area=range(1541,1550)
                else:
                    area=range(2901,2906)
                for i in area:
                    index = data_dic['data1'][nucl_data[sub_select]].index(float(i))
                    hundred_1 += data_dic['data1'][subunit[sub_select]][index]
                    hundred_2 += data_dic['data2'][subunit[sub_select]][index]

            heights_1 = [(x / hundred_1 * 100) for x in data_dic['data1'][subunit[sub_select]]]
            heights_2 = [(x / hundred_2 * 100) for x in data_dic['data2'][subunit[sub_select]]]
            self.widget2.canvas.ax.set_ylabel('Relative percentage of reads', fontsize='large')
            self.widget2.canvas.ax.set_xlabel('Nucleotide Position', fontsize='large')
            self.widget2.canvas.ax.set_title(plotname[sub_select], fontsize='large')
            series1 = self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], heights_1, alpha=0.5, facecolors=data_dic['data1'][subunitc[sub_select]],
                                                     picker=True, marker = data_dic['data1']['symbol'], label=' '.join(proc1), linewidth='1')
            series2 = self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], heights_2, alpha=0.5, facecolors=data_dic['data2'][subunitc[sub_select]],
                                                     picker=True, marker = data_dic['data2']['symbol'], label=' '.join(proc2), linewidth='1')
            self.widget2.canvas.ax.set_ylim(-20,max(heights_1 + heights_2) + 10)
            mazf = self.widget2.canvas.ax.scatter(0,-1000, alpha=0.5, c='red', marker = 'o', label = ' _ACA')
            mqsr = self.widget2.canvas.ax.scatter(0,-1000, alpha=0.5, c='cyan', marker = 'o', label = 'G_CB')
            self.widget2.canvas.ax.legend(handles=[series1, series2, mazf, mqsr],loc='best', scatterpoints = 1)

            self.widget2.canvas.draw()

        elif value_select == 1:
            self.widget2.canvas.ax.clear()
            self.widget2.canvas.ax.set_ylabel('Read Counts', fontsize='large')
            self.widget2.canvas.ax.set_xlabel('Nucleotide Position', fontsize='large')
            self.widget2.canvas.ax.set_title(plotname[sub_select], fontsize='large')
            series1 = self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], data_dic['data1'][subunit[sub_select]], alpha=0.5, c=data_dic['data1'][subunitc[sub_select]],
                                                     picker=True, marker = data_dic['data1']['symbol'], label=' '.join(proc1))
            self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], [-1 * data for data in data_dic['data1'][subunit_neg[sub_select ]]], alpha=0.5, c=data_dic['data1'][subunitc[sub_select]],
                                                     picker=True, marker = data_dic['data1']['symbol'])
            series2 = self.widget2.canvas.ax.scatter(data_dic['data2'][nucl_data[sub_select]], data_dic['data2'][subunit[sub_select]], alpha=0.5, c=data_dic['data2'][subunitc[sub_select]],
                                                     picker=True, marker = data_dic['data2']['symbol'], label=' '.join(proc2))
            self.widget2.canvas.ax.scatter(data_dic['data2'][nucl_data[sub_select]], [-1 * data for data in data_dic['data2'][subunit_neg[sub_select ]]], alpha=0.5, c=data_dic['data2'][subunitc[sub_select]],
                                                     picker=True, marker = data_dic['data2']['symbol'])
            self.widget2.canvas.fig.tight_layout()
            max_height = max(data_dic['data1'][subunit[sub_select]] + data_dic['data2'][subunit[sub_select]])
            self.widget2.canvas.ax.set_ylim(-0.2*max_height,max_height)
            mazf = self.widget2.canvas.ax.scatter(0,-0.3*max_height, alpha=0.5, c='red', marker = 'o', label = ' _ACA')
            mqsr = self.widget2.canvas.ax.scatter(0,-0.3*max_height, alpha=0.5, c='cyan', marker = 'o', label = 'G_CB')
            self.widget2.canvas.ax.legend(handles=[series1, series2, mazf, mqsr],loc='best', scatterpoints = 1)
            self.widget2.canvas.draw()

    #==========================================================================
    # Change labels in cboxes requried for 3-prime and 5-prime
    def on_combo_prime_change(self, index):
        text = self.primcombo.currentText()
        if text == '5prime':
            self.proc2combo.clear()
            self.proc2combo.addItems(items5prime.split(','))
            self.proc1combo.clear()
            self.proc1combo.addItems(items5prime.split(','))
        else:
            self.proc2combo.clear()
            self.proc2combo.addItems(items3prime.split(','))
            self.proc1combo.clear()
            self.proc1combo.addItems(items3prime.split(','))


    #======================================================
    # Return information when datapoint is clicked
    def onpick(self, event):

        global data_dic

        proc1 = self.proc1combo.currentText()
        proc2 = self.proc2combo.currentText()
        prim_select = self.primcombo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')
        prim_select = str(prim_select)

        subunit = ['y_pos_16S', 'y_pos_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        fastas = ['data/16S_new.fasta', 'data/23S_new.fasta']

        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.
        for fasta in SeqIO.parse(fastas[sub_select], "fasta"):
            fasta = fasta.seq
            pass

        # Relative ja absolute joonestamine. Value_select 0 - relative, 1 - absolute
        if value_select == 1:
            nucleotide_pos = np.take(data_dic['data1'][nucl_data[sub_select]], ind)
            sample_1_pos_read_count = np.take(data_dic['data1'][subunit[sub_select]], ind)
            sample_2_pos_read_count = np.take(data_dic['data2'][subunit[sub_select]], ind)
            sample_1_neg_read_count = np.take(data_dic['data1'][subunit_neg[sub_select]], ind)
            sample_2_neg_read_count = np.take(data_dic['data2'][subunit_neg[sub_select]], ind)

            for array_ind in range(len(nucleotide_pos)):
                position = int(nucleotide_pos[array_ind]) + 114
                if prim_select == '5prime':
                    sequence = str(fasta[position-3:position]+ '_' + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                else:
                    sequence = str(fasta[position-3:position] + '<b>' + fasta[position] + '</b>' + '_' + fasta[position+1:position+4])
                sequence = ("").join([x if x != 'T' else 'U'for x in sequence ])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,} <br>{5} '+' strand: {1:,}<br>{6} '+' strand {2:,}\
                                <br>{5} '-' strand: {3:,}<br>{6} '-' strand {4:,}<br>=======================<br>"\
                                .format(int(nucleotide_pos[array_ind]),
                                        sample_1_pos_read_count[array_ind],
                                        sample_2_pos_read_count[array_ind],
                                        sample_1_neg_read_count[array_ind],
                                        sample_2_neg_read_count[array_ind],
                                        ' '.join(proc1),
                                        ' '.join(proc2),
                                        sequence)

                self.text.insertHtml(message)
                self.text.moveCursor(QtGui.QTextCursor.Start)

        elif value_select == 0:

            hundred_1 = 0
            hundred_2 = 0
            #===============================================================
            # Selection of locations for 100% for 5-prim and 3-prim

            if prim_select == '5prime':
                if sub_select == 0:
                    area=range(-6,7)
                else:
                    area=range(0,14)
                for i in area:
                    index = data_dic['data1'][nucl_data[sub_select]].index(float(i))
                    hundred_1 += data_dic['data1'][subunit[sub_select]][index]
                    hundred_2 += data_dic['data2'][subunit[sub_select]][index]
            elif prim_select == '3prime':
                if sub_select == 0:
                    area=range(1541,1550)
                else:
                    area=range(2901,2906)
                for i in area:
                    index = data_dic['data1'][nucl_data[sub_select]].index(float(i))
                    hundred_1 += data_dic['data1'][subunit[sub_select]][index]
                    hundred_2 += data_dic['data2'][subunit[sub_select]][index]
            heights_1 = [(x / hundred_1 * 100) for x in data_dic['data1'][subunit[sub_select]]]
            heights_2 = [(x / hundred_2 * 100) for x in data_dic['data2'][subunit[sub_select]]]
            ind = event.ind


            nucleotide_pos = np.take(data_dic['data1'][nucl_data[sub_select]], ind)
            sample_1_pos_read_count = np.take(heights_1, ind)
            sample_2_pos_read_count = np.take(heights_2, ind)
            sample_1_pos_abs_count = np.take(data_dic['data1'][subunit[sub_select]], ind)
            sample_2_pos_abs_count = np.take(data_dic['data2'][subunit[sub_select]], ind)

            for array_ind in range(len(nucleotide_pos)):
                position = int(nucleotide_pos[array_ind]) + 114
                if prim_select == '5prime':
                    print
                    if fasta[position:position+3] in ['ACA','GCU']:
                        threeprimeside = '<i>' + fasta[position+1:position+4] +'<\i>'
                    else:
                        threeprimeside = fasta[position+1:position+4]

                    sequence = str(fasta[position-3:position]+ '_' + '<b>' + fasta[position] + '</b>' + threeprimeside)
                else:
                    if fasta[position+1:position+4] in ['ACA','GCU']:
                        threeprimeside = '<i>' + fasta[position+1:position+4] +'<\i>'
                    else:
                        threeprimeside = fasta[position+1:position+4]
                    sequence = str(fasta[position-4:position-1]+ '_' + '<b>' + fasta[position-1] + '</b>' + '_' + fasta[position:position+3])
                #===============================================
                # This awkward piece of code effectively removes the position 0 from the sequence,
                # as position 0 does not really exist in DNA.
                reported_position = int(nucleotide_pos[array_ind])
                if reported_position <= 0:
                    reported_position -= 1
                sequence = ("").join([x if x != 'T' else 'U'for x in sequence ])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,}<br>Percentages:<br>{3}: {1}<br>{4}:\
                            {2}<br>Absolute values:<br>{3}: {5:,}<br>{4}: {6:,}<br>=======================<br>"\
                            .format(reported_position,
                                    round(sample_1_pos_read_count[array_ind],4),
                                    round(sample_2_pos_read_count[array_ind], 4),
                                    ' '.join(proc1), ' '.join(proc2),
                                    int(sample_1_pos_abs_count[array_ind]),
                                    int(sample_2_pos_abs_count[array_ind]),
                                    sequence)

                self.text.insertHtml(message)
                self.text.moveCursor(QtGui.QTextCursor.Start)

# Enables saving data from the plotter program. Not developed for the moment as there is no need for it.
    def save(self):

        proc1 = self.proc1combo.currentText()
        proc2 = self.proc2combo.currentText()
        prim_select = self.primcombo.currentText()
        sub_select = self.subcombo.currentIndex()
        #value_select = self.valcombo.currentIndex() #Necessary if you want to be able to save data in both percentage and absolute value
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')
        prim_select = str(prim_select)

        subunit = ['y_pos_16S', 'y_pos_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']

        global data_dic

        path, _ = QtGui.QFileDialog.getSaveFileName(self, 'Save data', os.getcwd(), selectedFilter='*.csv')

        index = ['data1', 'data2']
        with open(path, 'w') as outfile:
            write = csv.writer(outfile, dialect='excel')
            col1 = data_dic['data1'][nucl_data[sub_select]]
            col2 = data_dic['data1'][subunit[sub_select]]
            col3 = data_dic['data1'][subunit_neg[sub_select]]
            col4 = data_dic['data2'][nucl_data[sub_select]]
            col5 = data_dic['data2'][subunit[sub_select]]
            col6 = data_dic['data2'][subunit_neg[sub_select]]
            data = zip(col1, col2, col3, col4, col5, col6)
            for line in data:
                write.writerow(line)



