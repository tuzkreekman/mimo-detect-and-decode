from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import os

def generate_data(n = 10000, l = 0.1):
  '''
  Fill in your implementation here.
  '''
  X = np.random.normal(size=(n,6))
  scale = l
  N1 = X[:,:5] + X[:,1:] + np.random.exponential(scale=scale,size=(n,5)) 
  N2 = np.ndarray((n,6))
  N2[:,:5] = N1 
  N2[:,1:] = N2[:,1:] + N1 
  N2 = N2 + np.random.exponential(scale=scale,size=(n,6)) 
  Y = N2[:,:5] + N2[:,1:] + np.random.exponential(scale=scale,size=(n,5)) 
  return X,Y


def generate_batch(data, batch_size = 100):
  x, y = data
  indices = np.random.choice(len(x), batch_size, replace=False)
  return x[indices], y[indices]


class model:
  def __init__(self, num_layer = 1, num_neuron = 100):
    self.num_layer = num_layer
    self.num_neuron = num_neuron

  def layer(self, x, num_neuron, i, activation):
    input_shape = x.get_shape()
    W = tf.get_variable("weights_%d" %i, 
      [input_shape[-1], num_neuron], initializer = tf.truncated_normal_initializer(stddev = 0.1))
    b = tf.get_variable("bias_%d" %i, 
      [1, num_neuron], initializer = tf.truncated_normal_initializer(stddev = 0.1))
    y = tf.matmul(x, W) + b
    return activation(y)
  
  def construct_graph(self):
    tf.reset_default_graph()
    self.x = tf.placeholder(tf.float32, [None, 5])
    h = self.x
    '''
    Fill in your implementation here.
    '''
    for i in range(self.num_layer):
        h = self.layer(h,self.num_neuron,i,tf.nn.relu)
    #output layer
    h = self.layer(h,6,self.num_layer,lambda x:x)
    y = h

    self.prediction = y
    self.y_ = tf.placeholder(tf.float32, [None, 6])
    self.loss = tf.reduce_mean(tf.square(y - self.y_))
    self.train_step = tf.train.AdamOptimizer(1e-4).minimize(self.loss)
    
    return tf.Session(graph=tf.get_default_graph())
    
  def train(self, sess, data):
    tf.global_variables_initializer().run(session=sess)
    losses = []
    batch_size = 1000
    for _ in range(1000):
      batch_xs, batch_ys = generate_batch(data, batch_size)
      loss, _ = sess.run([self.loss, self.train_step], feed_dict={self.x: batch_xs, self.y_: batch_ys})
      losses.append(loss)
    plt.plot(losses)
    plt.title('training loss')
    plt.show()

  def test(self, sess, t):
    return sess.run(self.prediction, feed_dict = {self.x: t})

if __name__ == '__main__':

  '''
  Fill in your implementation here.
  '''
  X_train, Y_train = generate_data(10000,.1)
  X_test , Y_test  = generate_data(1000 ,.1)
  layers = 2
  neurons = 150
  m = model(layers, neurons)
  sess = m.construct_graph()
  m.train(sess, (Y_train, X_train))
  X_pred = m.test(sess, Y_test)
  print("Neural network rmse:",np.linalg.norm(X_pred-X_test)/1000)
  cov = np.cov(np.concatenate([X_train,Y_train],axis=1).T)
  cov_y = cov[6:,6:]
  cov_yx = cov[6:,:6]
  h = np.dot(np.linalg.inv(cov_y),cov_yx)
  X_lin = np.apply_along_axis(lambda x:np.dot(h.T,x),1,Y_test)
  print("Linear estimate rmse:",np.linalg.norm(X_lin-X_test)/1000)



  

