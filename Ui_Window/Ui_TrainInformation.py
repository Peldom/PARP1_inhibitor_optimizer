# -*- coding: utf-8 -*-


'''
author: bin

训练信息子窗体Ui设计类

'''

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_TrainInformation(QtWidgets.QWidget):
    def __init__(self, inf):
        super().__init__()
        self.setMinimumSize(400, 400)

        verticalLayout = QtWidgets.QVBoxLayout(self)
        scrollArea = QtWidgets.QScrollArea()
        scrollArea.setWidgetResizable(True)
        

        formLayoutWidget = QtWidgets.QWidget()
        formLayout = QtWidgets.QFormLayout(formLayoutWidget)

        for (key, value) in inf:
            formLayout.addRow(
                QtWidgets.QLabel(key, self),
                QtWidgets.QLabel(str(value), self))


        scrollArea.setWidget(formLayoutWidget)
        verticalLayout.addWidget(scrollArea)
