import sys
import matplotlib
matplotlib.use("Qt4Agg")
from matplot import *
from PySide import QtGui
import matplotlib.pyplot as plt

plt.style.use('ggplot')


class GUIForm(QtGui.QMainWindow):

    def __init__(self, parent=None):

        super(GUIForm, self).__init__()

        self.initUI()


    def initUI(self):

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setCentralWidget(self.ui)

        exitAction = QtGui.QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

        compareAction = QtGui.QAction('Compare/Plot', self)
        compareAction.triggered.connect(self.set)

        saveAction = QtGui.QAction('&Save plot data', self)
        saveAction.triggered.connect(self.save_data)

        self.menuBar = self.menuBar()
        fileMenu = self.menuBar.addMenu('&File')
        fileMenu.addAction(exitAction)
        fileMenu.addAction(saveAction)
        editMenu = self.menuBar.addMenu('&Edit')
        editMenu.addAction(compareAction)

        self.show()

    def set(self):
        status = self.ui.stackedLayout.currentIndex()

        if status == 0:
            self.ui.stackedLayout.setCurrentIndex(1)
        else:
            self.ui.stackedLayout.setCurrentIndex(0)

    def save_data(self):
        status = self.ui.stackedLayout.currentIndex()
        if status == 0:
            self.ui.plotter.save()
        else:
            self.ui.compare.save()

class Ui_MainWindow(QtGui.QWidget):

    def setupUi(self, MainWindow):
        self.plotter = Plotter()
        self.plotter.setupUi()
        #self.compare = Compare()
        #self.compare.setupUi()
        self.stackedLayout = QtGui.QStackedLayout()
        self.stackedLayout.addWidget(self.plotter)
        #self.stackedLayout.addWidget(self.compare)
        self.stackedLayout.setCurrentIndex(0)

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addLayout(self.stackedLayout)

        self.setLayout(mainLayout)





if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    sys.exit(app.exec_())