import tensorflow as tf
print(tf.__version__)
a=tf.constant([2,3,5],dtype=tf.float32)
b=tf.nn.softmax(a)
print(b)