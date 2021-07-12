

import numpy as np
import random
import time

from PyQt5 import QtCore



class ORGANIC(QtCore.QObject):
    """Main class, where every interaction between the user
    and the backend is performed.
    """

    signal_inf = QtCore.pyqtSignal(list)
    signal_progress = QtCore.pyqtSignal(str, int, int)
    signal_draw = QtCore.pyqtSignal(str, dict)
    signal_message = QtCore.pyqtSignal(str, str)

    def range(self, name, iterable):
        """ 进度条控制

            name:对应进度条名

            iterable:整数或可迭代对象
        """
        if type(iterable) == int:
            loop = iterable
            iterable = range(iterable)
        else:
            loop = len(iterable)
        step = 0
        for item in iterable:
            self.signal_progress.emit(name, step, loop)
            yield item
            step += 1
        self.signal_progress.emit(name, step, loop)
    

    def __init__(self, name, params={}):
        super().__init__()

        # Set parameters 此处建议如下简化，可不修改
        self.PREFIX = name
        self.PRETRAIN_GEN_EPOCHS = params['PRETRAIN_GEN_EPOCHS']
        self.PRETRAIN_DIS_EPOCHS = params['PRETRAIN_DIS_EPOCHS']
        self.GEN_ITERATIONS = params['GEN_ITERATIONS']
        self.GEN_BATCH_SIZE = params['GEN_BATCH_SIZE']
        self.SEED = params['SEED']
        random.seed(self.SEED)
        np.random.seed(self.SEED)
        self.DIS_BATCH_SIZE = params['DIS_BATCH_SIZE']
        self.DIS_EPOCHS = params['DIS_EPOCHS']
        self.EPOCH_SAVES = params['EPOCH_SAVES']
        self.GEN_EMB_DIM = params['GEN_EMB_DIM']
        self.GEN_HIDDEN_DIM = params['GEN_HIDDEN_DIM']
        self.SAMPLE_NUM = params['SAMPLE_NUM']
        self.BIG_SAMPLE_NUM = params['BIG_SAMPLE_NUM'] if params['BIG_SAMPLE_NUM'] else self.SAMPLE_NUM * 5
        self.LAMBDA = params['LAMBDA']
        self.MAX_LENGTH = params['MAX_LENGTH']
        self.DIS_EMB_DIM = params['DIS_EMB_DIM']
        self.DIS_FILTER_SIZES = params['DIS_FILTER_SIZES']
        self.DIS_NUM_FILTERS = params['DIS_FILTER_SIZES']
        self.DIS_DROPOUT = params['DIS_DROPOUT']
        self.DIS_L2REG = params['DIS_L2REG']
        self.START_TOKEN = params['START_TOKEN']
        self.CHK_PATH = params['CHK_PATH']


    def load_training_set(self, file):

        # 模拟运行，需删除
        self.train_samples = "123456"
        self.POSITIVE_NUM = self.NUM_EMB = random.randint(0, 100)
        to_use = list(map(str,np.random.normal(size=100)))

        # 进程通信 information 
        self.signal_inf.emit([
            ('模型名称' , self.PREFIX),
            ('训练数据集大小' , len(self.train_samples)),
            ('最大数据长度' , self.MAX_LENGTH),
            ('平均数据长度' , np.mean([len(s) for s in to_use])),
            ('Num valid data points' , self.POSITIVE_NUM),
            ('Size of alphabet' , self.NUM_EMB)])


    # 此处无改动
    def set_training_program(self, metrics=None, steps=None):
        """Sets a program of metrics and epochs
        for training the model and generating molecules.

        """

        # 模拟运行，需删除
        self.TOTAL_BATCH = sum(steps)
        


    # 此处无改动
    def load_metrics(self):
        """Loads the metrics."""

        pass

    
    def pretrain(self):
        """Pretrains generator and discriminator."""

        # 模拟运行，需删除
        pre_gen_loss = np.random.normal(size=self.PRETRAIN_GEN_EPOCHS)
        pre_dis_loss = np.random.normal(size=self.PRETRAIN_DIS_EPOCHS)
        pre_dis_accuracy = np.random.lognormal(size=self.PRETRAIN_DIS_EPOCHS)


        for epoch in self.range("pre_gen", self.PRETRAIN_GEN_EPOCHS):
            
            # 模拟运行，需删除
            time.sleep(0.01)
            loss = pre_gen_loss[epoch]

            # 进程通信 
            self.signal_draw.emit("pre_gen",{
                "pre_gen_loss":loss})



        for i in self.range("pre_dis", self.PRETRAIN_DIS_EPOCHS):
            
            # 模拟运行，需删除
            time.sleep(0.02)
            loss = pre_dis_loss[i]
            accuracy = pre_dis_accuracy[i]

            # 进程通信 
            self.signal_draw.emit("pre_dis",{
                "pre_dis_loss":loss,
                "pre_dis_accuracy":accuracy})




    def train(self, ckpt_dir='checkpoints/'):
        """Trains the model. If necessary, also includes pretraining."""

        self.pretrain()

        # 模拟运行，需删除
        path = "./try"

        # 进程通信
        self.signal_message.emit("pre", "预训练保存于：{}".format(path))


    
        for nbatch in self.range("batch", self.TOTAL_BATCH):

            for it in self.range("gen", self.GEN_ITERATIONS):
                
                # 模拟运行，需删除
                time.sleep(0.1)
                rewards = np.random.normal(size=self.SAMPLE_NUM)
                
                # 进程通信 
                self.signal_draw.emit("gen",{
                    "gen_mean":np.mean(rewards),
                    "gen_std":np.std(rewards),
                    "gen_min":np.min(rewards),
                    "gen_max":np.max(rewards)})


            # 模拟运行，需删除
            dis_loss = np.random.normal(size=self.DIS_EPOCHS)
            dis_accuracy = np.random.lognormal(size=self.DIS_EPOCHS)

            # generate for discriminator
            for i in self.range("dis", self.DIS_EPOCHS):

                # 模拟运行，需删除
                time.sleep(0.1)
                loss = dis_loss[i]
                accuracy = dis_accuracy[i]

                # 进程通信 
                self.signal_draw.emit("dis",{   
                    "dis_loss":loss,
                    "dis_accuracy":accuracy})

         
                    
        # 进程通信 save models
        self.signal_message.emit("train", "模型保存于：{}".format(path))







