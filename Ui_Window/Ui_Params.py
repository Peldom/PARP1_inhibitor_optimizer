# -*- coding: utf-8 -*-


'''
author: bin

训练参数子窗体类

'''

import sys
import json
import re
import os

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Params(QtWidgets.QDockWidget):
    def __init__(self, Form, setting):
        super().__init__(Form)
        self.setFeatures(QtWidgets.QDockWidget.AllDockWidgetFeatures)

        self.setWindowTitle("模型训练参数")

        self.settingPath = Form.settingPath

        verticalLayoutWidget = QtWidgets.QWidget()
        verticalLayoutWidget.setMinimumSize(350, 400)
        verticalLayout = QtWidgets.QVBoxLayout(verticalLayoutWidget)

        # 参数区

        scrollArea = QtWidgets.QScrollArea()
        scrollArea.setWidgetResizable(True)
        scrollArea.setWidget(self.initScrollAreaUi(setting))
        verticalLayout.addWidget(scrollArea)

        # 底部按钮

        horizontalLayoutWidget = QtWidgets.QWidget()
        horizontalLayout = QtWidgets.QHBoxLayout(horizontalLayoutWidget)

        self.checkBtn = QtWidgets.QPushButton("检查参数", horizontalLayoutWidget)
        self.checkBtn.setObjectName("paramsCheckBtn")
        horizontalLayout.addWidget(self.checkBtn)

        self.resetBtn = QtWidgets.QPushButton("重置参数", horizontalLayoutWidget)
        self.resetBtn.setObjectName("paramsResetBtn")
        horizontalLayout.addWidget(self.resetBtn)

        verticalLayout.addWidget(horizontalLayoutWidget)
        self.setWidget(verticalLayoutWidget)
        QtCore.QMetaObject.connectSlotsByName(self)


    def initScrollAreaUi(self, setting):
        formLayoutWidget = QtWidgets.QWidget()
        formLayout = QtWidgets.QFormLayout(formLayoutWidget)


        self.label, self.Field = {}, {}

        def formLAddRow(key, value):
            self.label[key] = QtWidgets.QLabel(
                value["name"] if value["name"] else key, 
                formLayoutWidget)
            self.Field[key] = QtWidgets.QLineEdit(str(value["value"]), formLayoutWidget)
            self.lineEditAddHelp(self.Field[key], value["help"])
            formLayout.addRow(self.label[key], self.Field[key])

        formLAddRow("name", setting["name"])
        for key, value in setting["params"].items():
            formLAddRow(key, value)

        # 数据集

        self.label["data_set"] = QtWidgets.QLabel(setting["data_set"]["name"], formLayoutWidget)

        verticalLayout = QtWidgets.QVBoxLayout()

        self.choseDataSetBtn = QtWidgets.QPushButton("选择数据集", formLayoutWidget)
        self.choseDataSetBtn.setObjectName("choseDataSetBtn")
        verticalLayout.addWidget(self.choseDataSetBtn)

        self.Field["data_set"] = QtWidgets.QLineEdit(str(setting["data_set"]["value"]), formLayoutWidget)
        self.lineEditAddHelp(self.Field["data_set"], setting["data_set"]["help"])
        verticalLayout.addWidget(self.Field["data_set"])

        formLayout.addRow(self.label["data_set"], verticalLayout)

        return formLayoutWidget


    def lineEditAddHelp(self, lineEdit, value):
        action = QtWidgets.QAction(QtGui.QIcon("img/question.ico"), value, lineEdit)
        action.triggered.connect(self.on_iconShowHelp_clicked)
        action.setStatusTip(value)
        lineEdit.addAction(action, QtWidgets.QLineEdit.TrailingPosition)



    @QtCore.pyqtSlot()
    def on_iconShowHelp_clicked(self):
        _ = QtWidgets.QMessageBox.information(self,'提示信息', self.sender().text())
        

    @QtCore.pyqtSlot()
    def on_choseDataSetBtn_clicked(self):
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择数据集", 
            "./", ("dataset (*.csv *.smi)"))
        if not fileName:
            return
        self.Field["data_set"].setText(fileName)


    def TypeTrans(self, key, pType, settingWidget, std):
        def checkNull(x):
            if x == "None":
                return None
            else:
                raise Exception("此处可空")

        def checkPath(x):
            if os.path.exists(x):
                return x
            else:
                raise Exception("路径错误或不存在")
        
        switch = {
            "str": lambda x: x,
            "int": lambda x: int(x),
            "float": lambda x: float(x),
            "null": checkNull,
            "path": checkPath,
            "int null": lambda x: None if x == "None" else int(x),
            "list-int": lambda x: list(map(int, re.split(r",\s",x.strip("[]"))))
        }

        try:
            if  ("check" in std) and re.match(std["check"], 
                settingWidget.Field[key].text().strip()) == None :
                raise Exception("请详细参看说明")
            return switch[pType](settingWidget.Field[key].text().strip())
        except Exception as e:
            raise Exception(settingWidget.label[key].text() + " 格式错误 " + str(e))



    def getTranParams(self, stdSetting):
        tSetting = {}
        tSetting["name"] = self.TypeTrans("name", stdSetting["name"]["type"], 
            self, stdSetting["name"])
        tSetting["params"] = {}
        for key, value in stdSetting["params"].items(): 
            tSetting["params"][key] = self.TypeTrans(key, value["type"], self, value)
        tSetting["data_set"] = self.TypeTrans("data_set", stdSetting["data_set"]["type"], 
            self, stdSetting["data_set"])
        return tSetting
    

    @QtCore.pyqtSlot()    
    def on_paramsCheckBtn_clicked(self):
        try:
            with open(self.settingPath,'r') as setFile:
                stdSetting = json.load(setFile)
            _ = self.getTranParams(stdSetting)
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            _ = QtWidgets.QMessageBox.information(self,'提示', "参数正确")


    def paramsReset(self, setting):
        self.Field["name"].setText(str(setting["name"]["value"]))
        for key, value in setting["params"].items():
            self.Field[key].setText(str(value["value"]))
        self.Field["data_set"].setText(str(setting["data_set"]["value"]))


    @QtCore.pyqtSlot()
    def on_paramsResetBtn_clicked(self):
        with open(self.settingPath,'r') as setFile:
            setting = json.load(setFile)
        self.paramsReset(setting)
        _ = QtWidgets.QMessageBox.information(self,'提示', "训练参数重置完毕")



