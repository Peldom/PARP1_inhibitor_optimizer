# -*- coding: utf-8 -*-


'''
author: bin

页面基础逻辑类

'''

import json
import time
import sys
import os
import re

from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Main import Ui_Main

class mainWindow(Ui_Main):
    def __init__(self):
        self.settingPath = "./Ui_Window/stdSetting.json"

        super().__init__(self.settingPath)

        



    # 通用函数


    @QtCore.pyqtSlot()
    def on_iconShowHelp_clicked(self):
        _ = QtWidgets.QMessageBox.information(self,'提示信息', self.sender().text())


    def getSetting(self):
        try:
            with open(self.settingPath,'r') as setFile:
                stdSetting = json.load(setFile)
            tSetting = self._getTranParams(stdSetting)
            tSetting["metrics"] = self._getParticular(stdSetting)
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            self.uerSetting = tSetting


    def _TypeTrans(self, key, pType, settingWidget, std):
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
            


    # 训练参数设置页面


    @QtCore.pyqtSlot()
    def on_choseDataSetBtn_clicked(self):
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择数据集", 
            "./", ("dataset (*.csv *.smi)"))
        if not fileName:
            return
        self.trainParams.Field["data_set"].setText(fileName)


    def _getTranParams(self, stdSetting):
        tSetting = {}
        tSetting["name"] = self._TypeTrans("name", stdSetting["name"]["type"], 
            self.trainParams, stdSetting["name"])
        tSetting["params"] = {}
        for key, value in stdSetting["params"].items(): 
            tSetting["params"][key] = self._TypeTrans(key, value["type"], self.trainParams, value)
        tSetting["data_set"] = self._TypeTrans("data_set", stdSetting["data_set"]["type"], 
            self.trainParams, stdSetting["data_set"])
        return tSetting
    

    @QtCore.pyqtSlot()    
    def on_paramsCheckBtn_clicked(self):
        try:
            with open(self.settingPath,'r') as setFile:
                stdSetting = json.load(setFile)
            _ = self._getTranParams(stdSetting)
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            _ = QtWidgets.QMessageBox.information(self,'提示', "参数正确")


    def _paramsReset(self, setting):
        self.trainParams.Field["name"].setText(str(setting["name"]["value"]))
        for key, value in setting["params"].items():
            self.trainParams.Field[key].setText(str(value["value"]))
        self.trainParams.Field["data_set"].setText(str(setting["data_set"]["value"]))


    @QtCore.pyqtSlot()
    def on_paramsResetBtn_clicked(self):
        with open(self.settingPath,'r') as setFile:
            setting = json.load(setFile)
        self._paramsReset(setting)
        _ = QtWidgets.QMessageBox.information(self,'提示', "训练参数重置完毕")


    # 训练指标设置页面


    def _getParticular(self, stdSetting):
        tSetting = {}
        for key, value in self.trainParticular.label.items(): 
            if value.isChecked() :
                tSetting[key] = self._TypeTrans(key, "int", self.trainParticular, stdSetting["metrics"][key])
        return tSetting
        

    @QtCore.pyqtSlot()
    def on_particularCheckBtn_clicked(self):
        try:
            with open(self.settingPath,'r') as setFile:
                stdSetting = json.load(setFile)
            _ = self._getParticular(stdSetting)
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            _ = QtWidgets.QMessageBox.information(self,'提示', "训练指标正确")


    def _particularReset(self, setting):
        for key, value in setting["metrics"].items():
            if value["recommend"] :
                self.trainParticular.label[key].setChecked(True)
                self.trainParticular.Field[key].setText(
                    str(value["value"]) if value["value"] else "")
            else:
                self.trainParticular.label[key].setChecked(False)
                self.trainParticular.Field[key].setText(None)


    @QtCore.pyqtSlot()
    def on_particularResetBtn_clicked(self):
        with open(self.settingPath,'r') as setFile:
            setting = json.load(setFile)
        self._particularReset(setting)
        _ = QtWidgets.QMessageBox.information(self,'提示', "训练指标重置完毕")


    # 工具栏

    @QtCore.pyqtSlot()
    def on_settingTool_clicked(self):
        self.trainSetting.show()
        self.trainParams.show()
        self.trainParticular.show()
    
    
    # 主页面


    @QtCore.pyqtSlot()
    def on_settingCheckBtn_clicked(self):
        try:
            self.getSetting()
        except Exception as err:
            _ = QtWidgets.QMessageBox.warning(self,'警告',str(err))
        else:
            _ = QtWidgets.QMessageBox.information(self,'提示', "训练参数，训练指标正确")

    
    @QtCore.pyqtSlot()
    def on_settingResetBtn_clicked(self):
        with open(self.settingPath,'r') as setFile:
            setting = json.load(setFile)
        self._paramsReset(setting)
        self._particularReset(setting)
        _ = QtWidgets.QMessageBox.information(self,'提示', "训练参数，训练指标重置完毕")


    @QtCore.pyqtSlot()
    def on_settingLoadBtn_clicked(self):
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择设置", 
            "./setting/", ("setting (*.json)"))
        if not fileName:
            return
        try:
            # 参数
            with open(fileName,'r') as setFile:
                self.uerSetting = json.load(setFile)
            self.trainParams.Field["name"].setText(str(self.uerSetting["name"]))
            for key, value in self.uerSetting["params"].items():
                self.trainParams.Field[key].setText(str(value))
            self.trainParams.Field["data_set"].setText(str(self.uerSetting["data_set"]))

            # 指标
            for key, value in self.trainParticular.label.items():
                if key in self.uerSetting["metrics"]:
                    value.setChecked(True)
                    self.trainParticular.Field[key].setText(
                        str(self.uerSetting["metrics"][key]))
                else:
                    value.setChecked(False)
                    self.trainParticular.Field[key].setText(None)
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
            json.dump(self.uerSetting, setFile, indent=4, ensure_ascii=False)
        _ = QtWidgets.QMessageBox.information(self,'提示', "保存成功")






if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = mainWindow()
    w.show()

    sys.exit(app.exec_())