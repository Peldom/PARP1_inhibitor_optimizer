# -*- coding: utf-8 -*-


'''
author: bin

训练指标子窗体类

'''

import sys
import json
import re
import os

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Particular(QtWidgets.QDockWidget):
    def __init__(self, Form, setting):
        super().__init__(Form)
        self.setFeatures(QtWidgets.QDockWidget.AllDockWidgetFeatures)
        self.setWindowTitle("模型训练指标")

        self.settingPath = Form.settingPath

        verticalLayoutWidget = QtWidgets.QWidget()
        verticalLayout = QtWidgets.QVBoxLayout(verticalLayoutWidget)
        verticalLayout.setContentsMargins(10, 10, 10, 10)

        
        self.comboBox = QtWidgets.QComboBox(verticalLayoutWidget)
        self.comboBox.activated[int].connect(self.helpLoad)
        # self.comboBox.setObjectName("particularComboBox")
        verticalLayout.addWidget(self.comboBox)

        
        self.lineEdit = QtWidgets.QLineEdit("", verticalLayoutWidget)
        action = QtWidgets.QAction(QtGui.QIcon("img/add.ico"), None, verticalLayoutWidget)
        action.setObjectName("addBtn")
        action.triggered.connect(self.on_addBtn_clicked)
        self.lineEdit.addAction(action, QtWidgets.QLineEdit.TrailingPosition)

        verticalLayout.addWidget(self.lineEdit)


        self.lable = QtWidgets.QLabel("说明", verticalLayoutWidget)
        verticalLayout.addWidget(self.lable)


        verticalLayout.addStretch(1)


        lable = QtWidgets.QLabel("指标名：训练步数", verticalLayoutWidget)
        verticalLayout.addWidget(lable)


        self.metrics, self.steps = [], []
        self.formLayoutWidget = QtWidgets.QWidget()
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        verticalLayout.addWidget(self.formLayoutWidget)


        verticalLayout.addStretch(1)


        horizontalLayoutWidget = QtWidgets.QWidget()
        horizontalLayout = QtWidgets.QHBoxLayout(horizontalLayoutWidget)

        self.removeBtn = QtWidgets.QPushButton("删除指标", horizontalLayoutWidget)
        self.removeBtn.setObjectName("removeBtn")
        horizontalLayout.addWidget(self.removeBtn)


        self.resetBtn = QtWidgets.QPushButton("重置指标", horizontalLayoutWidget)
        self.resetBtn.setObjectName("resetBtn")
        horizontalLayout.addWidget(self.resetBtn)

        verticalLayout.addWidget(horizontalLayoutWidget)


        self.setWidget(verticalLayoutWidget)


        self.key, self.help, self.data = [], [], []
        for key, value in setting["metrics"].items():
            self.comboBox.addItem(value["name"] if value["name"] else key)
            self.key.append(key)
            self.help.append(value["help"])
            self.data.append(value["value"])
            if value["recommend"]:
                self.addParticular(key, value["value"])

        QtCore.QMetaObject.connectSlotsByName(self)
        

    def helpLoad(self, index):
        self.lable.setText("帮助：" + self.help[index])
        self.lineEdit.setText(str(self.data[index]))


    def addParticular(self, key, value):
        self.metrics.append(key)
        self.steps.append(value)
        self.formLayout.addRow(
            QtWidgets.QLabel(key, self.formLayoutWidget),
            QtWidgets.QLabel(str(value), self.formLayoutWidget))


    def getParticular(self):
        return {
            "metrics" : self.metrics,
            "steps" : self.steps
        }


    def particularReset(self, setting):
        while self.formLayout.rowCount() > 0:
            self.on_removeBtn_clicked()
        for key, value in setting["metrics"].items():
            if value["recommend"] :
                self.addParticular(key, value["value"])


    @QtCore.pyqtSlot()
    def on_resetBtn_clicked(self):
        with open(self.settingPath,'r') as setFile:
            setting = json.load(setFile)
        self.particularReset(setting)
        _ = QtWidgets.QMessageBox.information(self,'提示', "训练指标重置完毕")


    @QtCore.pyqtSlot()
    def on_addBtn_clicked(self):
        try:
            self.addParticular(
                self.key[
                    self.comboBox.currentIndex()],
                int(self.lineEdit.text().strip())
            )
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))

    
    @QtCore.pyqtSlot()
    def on_removeBtn_clicked(self):
        index = self.formLayout.rowCount()
        if index > 0:
            self.formLayout.removeRow(index - 1)
            self.metrics.pop()
            self.steps.pop()


    


