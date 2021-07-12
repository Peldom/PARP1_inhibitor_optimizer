# -*- coding: utf-8 -*-


'''
author: bin

模型训练线程类

'''


from PyQt5 import QtCore, QtGui, QtWidgets

from ORGANIC.model.test import ORGANIC


class Td_Organic(QtCore.QThread):
    signal_inf = QtCore.pyqtSignal(list)
    signal_progress = QtCore.pyqtSignal(str, int, int)
    signal_finished = QtCore.pyqtSignal(str, str)
    signal_draw = QtCore.pyqtSignal(str, dict)
    signal_message = QtCore.pyqtSignal(str, str)


    def __init__(self, name, setting):
        super().__init__()
        self.setting = setting
        self.setObjectName(name)

        

    def run(self):
        ''' 此函数为ORGANIC的调用函数

            注意：ORGANIC需按样例修改，增加信号机制

            此处使用需添加 self.repeater(model) 保证线程间信号传递
        '''

        for i in self.range(self.setting["params"]["TRAIN_LOOP"]):
            model = ORGANIC(self.setting["name"], self.setting["params"])

            self.repeater(model) # 此处不可省略
            
            model.load_training_set(self.setting["data_set"])
            model.set_training_program(
                    self.setting["metrics"]["metrics"], 
                    self.setting["metrics"]["steps"])

            model.load_metrics()
            model.train()
        
        

    def repeater(self, model):
        """ 信号中继站
        
            model：ORGANIC的一个实例
        """
        model.signal_inf[list].connect(
            lambda inf: self.signal_inf.emit(
                inf +
                [("TOTAL_BATCH", sum(self.setting["metrics"]["steps"]))] +
                [(x, y) for x, y in self.setting["params"].items()] +
                [(x, y) for x, y in zip(
                    self.setting["metrics"]["metrics"], 
                    self.setting["metrics"]["steps"])]))
        model.signal_progress[str, int, int].connect(
            lambda name, step, loop: self.signal_progress.emit(
                name, step, loop))
        model.signal_draw[str, dict].connect(
            lambda name, tmp: self.signal_draw.emit(name, tmp))
        model.signal_message[str, str].connect(
            lambda name, message: self.signal_message.emit(name, message))
        
        


    def range(self, loop):
        step = 0
        while step < loop:
            self.signal_progress.emit("global", step, loop)
            yield step
            step += 1
        self.signal_progress.emit("global", step, loop)
        self.signal_finished.emit(self.objectName(), "模型训练完成")


    
    def test(self):
        print("here")