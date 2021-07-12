# -*- coding: utf-8 -*-


'''
author: bin

训练指标子窗体Ui设计类

'''

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Logo(QtWidgets.QDockWidget):
    def __init__(self, Form):
        super().__init__(Form)
        self.setFeatures(QtWidgets.QDockWidget.AllDockWidgetFeatures)
        self.setWindowTitle("首页")

        verticalLayoutWidget = QtWidgets.QWidget()
        verticalLayoutWidget.setMaximumSize(700, 700)
        verticalLayout = QtWidgets.QVBoxLayout(verticalLayoutWidget)
        verticalLayout.setContentsMargins(10, 10, 10, 10)

        lable = QtWidgets.QLabel(self)
        lable.setScaledContents(True)
        # lable.resize(500, 500)
        lable.setPixmap(QtGui.QPixmap("./img/logo.jpg"))

        verticalLayout.addWidget(lable)

        self.setWidget(verticalLayoutWidget)

