from PyQt4 import QtCore, QtGui

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

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1366, 1024)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))

        self.plotButton = QtGui.QPushButton(self.centralwidget)
        self.plotButton.setGeometry(QtCore.QRect(40, 600, 121, 27))
        self.plotButton.setObjectName(_fromUtf8("plotButton"))

        self.widget1 = maplot(self.centralwidget)
        self.widget1.setGeometry(QtCore.QRect(400, 20, 800, 400))
        self.widget1.setObjectName(_fromUtf8("MA plot"))

        self.widget2 = scatter(self.centralwidget)
        self.widget2.setGeometry(QtCore.QRect(400, 420, 800, 400))
        self.widget2.setObjectName(_fromUtf8("Scatter"))

        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1366, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))

        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.combo1 = QtGui.QComboBox(self.centralwidget)
        self.combo1.move(40, 100)
        self.combo1.addItems('MazF2h none,MazF2h PNK,MazF2h TAP,MqsR2h none,MqsR2h PNK,MqsR2h TAP'.split(','))

        self.combo2 = QtGui.QComboBox(self.centralwidget)
        self.combo2.move(40, 150)
        self.combo2.addItems('MazF2h none,MazF2h PNK,MazF2h TAP,MqsR2h none,MqsR2h PNK,MqsR2h TAP'.split(','))

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "RNA-Tools", None))
        self.plotButton.setText(_translate("MainWindow", "Plot", None))

from matplotlibwidget import maplot, scatter
