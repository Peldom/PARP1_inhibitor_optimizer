# -*- coding: utf-8 -*-


'''
author: bin

设置页面

'''

import time
import sys
import json
import re
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Window.Ui_Params import Ui_Params
from Ui_Window.Ui_Particular import Ui_Particular


class Ui_Setting(QtWidgets.QDockWidget):
    def __init__(self, Form):
        super().__init__(Form)
        self.setFeatures(QtWidgets.QDockWidget.AllDockWidgetFeatures)
        self.setWindowTitle("模型设置")

        self.settingPath = "./Ui_Window/stdSetting.json"
        with open(self.settingPath,'rb') as setFile:
            setting = json.load(setFile)

        self.trainParams = Ui_Params(self, setting)
        self.trainParticular = Ui_Particular(self, setting)


        formLayoutWidget = QtWidgets.QWidget()
        formLayout = QtWidgets.QFormLayout(formLayoutWidget)
        formLayout.setContentsMargins(10, 10, 10, 10)


        horizontalLayout = QtWidgets.QHBoxLayout()

        self.checkBtn = QtWidgets.QPushButton("模型训练参数", formLayoutWidget)
        self.checkBtn.setObjectName("paramsBtn")
        horizontalLayout.addWidget(self.checkBtn)

        self.resetBtn = QtWidgets.QPushButton("模型训练指标", formLayoutWidget)
        self.resetBtn.setObjectName("particularBtn")
        horizontalLayout.addWidget(self.resetBtn)

        formLayout.addRow(horizontalLayout)


        horizontalLayout = QtWidgets.QHBoxLayout()

        self.resetBtn = QtWidgets.QPushButton("重置设置", formLayoutWidget)
        self.resetBtn.setObjectName("settingResetBtn")
        horizontalLayout.addWidget(self.resetBtn)

        # formLayout.addRow(horizontalLayout)


        # horizontalLayout = QtWidgets.QHBoxLayout()

        self.loadBtn = QtWidgets.QPushButton("加载设置", formLayoutWidget)
        self.loadBtn.setObjectName("settingLoadBtn")
        horizontalLayout.addWidget(self.loadBtn)

        self.saveBtn = QtWidgets.QPushButton("保存设置", formLayoutWidget)
        self.saveBtn.setObjectName("settingSaveBtn")
        horizontalLayout.addWidget(self.saveBtn)

        formLayout.addRow(horizontalLayout)


        self.setWidget(formLayoutWidget)

        QtCore.QMetaObject.connectSlotsByName(self)


    def show(self):
        super().show()
        self.trainParams.show()
        self.trainParticular.show()


    def close(self):
        super().close()
        self.trainParams.close()
        self.trainParticular.close()


    @QtCore.pyqtSlot()
    def on_paramsBtn_clicked(self):
        self.trainParams.show()


    @QtCore.pyqtSlot()
    def on_particularBtn_clicked(self):
        self.trainParticular.show()

    
    @QtCore.pyqtSlot()
    def on_settingResetBtn_clicked(self):
        with open(self.settingPath,'rb') as setFile:
            setting = json.load(setFile)
        self.trainParams.paramsReset(setting)
        self.trainParticular.particularReset(setting)
        _ = QtWidgets.QMessageBox.information(self,'提示', "训练参数，训练指标重置完毕")


    @QtCore.pyqtSlot()
    def on_settingLoadBtn_clicked(self):
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择设置", 
            "./setting/", ("setting (*.json)"))
        if not fileName:
            return
        # else:
        try:
            # 参数
            with open(fileName,'rb') as setFile:
                setting = json.load(setFile)
            self.trainParams.Field["name"].setText(str(setting["name"]))
            for key, value in setting["params"].items():
                self.trainParams.Field[key].setText(str(value))
            self.trainParams.Field["data_set"].setText(str(setting["data_set"]))

            # 指标
            while self.trainParticular.formLayout.rowCount() > 0:
                self.trainParticular.on_particularRemoveBtn_clicked()
            for (key, value) in zip(
                setting["metrics"]["metrics"], 
                setting["metrics"]["steps"]):
                self.trainParticular.addParticular(key, value)
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            _ = QtWidgets.QMessageBox.information(self,'提示', "加载成功")


    @QtCore.pyqtSlot()
    def on_settingSaveBtn_clicked(self):
        self.getSetting()
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "保存设置",
            './setting/setting' + time.strftime("%Y%m%d%H%M%S", time.localtime()) + ".json",
            ("setting (*.json)"))
        if not fileName:
            return
        with open(fileName,'w+') as setFile:
            json.dump(self.getSetting(), setFile, indent=4, ensure_ascii=False)
        _ = QtWidgets.QMessageBox.information(self,'提示', "保存成功")


    def getSetting(self):
        try:
            with open(self.settingPath,'rb') as setFile:
                stdSetting = json.load(setFile)
            setting = self.trainParams.getTranParams(stdSetting)
            setting["metrics"] = self.trainParticular.getParticular()
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            return setting

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = Ui_Setting(None)
    w.show()
    sys.exit(app.exec_())