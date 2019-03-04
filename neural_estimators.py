# EE 278 Autumn 2017
# Balaji Prabhakar, Zi Yin, Pin Pin Tea-mangkornpan

from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import os

# Feel free to edit any parts of the code. 

def generate_data_AWGN(n = 1000, P = 2, N = 1):
  '''
  Fill in your implementation here.
  '''
  x = np.random.normal(0,np.sqrt(P),(n,1))
  z = np.random.normal(0,np.sqrt(N),(n,1))
  y = x+z
  return y, x

def generate_data_laplace(n = 1000, lambda_x = 2, lambda_n = 1):
  '''
  Fill in your implementation here.
  '''
  x = np.random.laplace(0,1/lambda_x,(n,1))
  z = np.random.normal(0,1/lambda_n,(n,1))
  y = x+z
  return y, x

def generate_data_unif_cube(n=1000, P = 2, N = 1):
  '''
  Fill in your implementation here.
  '''
  x = np.random.uniform(-P,P,(n,1))
  z = np.random.normal(-N,N,(n,1))
  y = x**3+z
  return y, x

def generate_batch(data, batch_size = 100):
  x, y = data
  indices = np.random.choice(len(x), batch_size, replace=False)
  return x[indices], y[indices]


class model:
  def __init__(self, num_layer = 1, num_neuron = 100):
    self.num_layer = num_layer
    self.num_neuron = num_neuron

  def layer(self, x, num_neuron, layer_id, activation):
    input_shape = x.get_shape()
    W = tf.get_variable("weights_%d" %layer_id, 
      [input_shape[-1], num_neuron], initializer = tf.truncated_normal_initializer(stddev = .1))
    b = tf.get_variable("bias_%d" %layer_id, 
      [1, num_neuron], initializer = tf.truncated_normal_initializer(stddev = .1))
    y = tf.matmul(x, W) + b
    return activation(y)

  def construct_graph(self):
    tf.reset_default_graph()
    self.x = tf.placeholder(tf.float32, [None, 1])
    h = self.x
    
	#hidden layers
    for i in range(self.num_layer):
        h = self.layer(h,self.num_neuron,i,tf.nn.relu)
    #output layer
    h = self.layer(h,1,self.num_layer,lambda x:x)
    y = h

    self.prediction = y
    self.y_ = tf.placeholder(tf.float32, [None, 1])
    self.loss = tf.reduce_mean(tf.square(y - self.y_))
    self.train_step = tf.train.GradientDescentOptimizer(1e-3).minimize(self.loss)
    return tf.Session(graph=tf.get_default_graph())
    
  def train(self, sess, data):
    tf.global_variables_initializer().run(session=sess)
    losses = []
    batch_size = 1000
    for _ in range(1000):
      batch_xs, batch_ys = generate_batch(data, batch_size)
      loss, _ = sess.run([self.loss, self.train_step], feed_dict={self.x: batch_xs, self.y_: batch_ys})
      losses.append(loss)
    # plt.plot(losses)
    # plt.title('training loss')
    # plt.show()

  def test(self, sess, t):
    return sess.run(self.prediction, feed_dict = {self.x: t})

  def draw(self, y, xhat, x, lin_est, name):
    fig = figure(1)
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, color = 'red', label='X')
    ax1.scatter(xhat, y, color='blue', label='Xhat')
    if lin_est is not None:
	    ax1.plot(lin_est, y, color='green', label='Linear estimate')
    ax1.grid(True)
    ax1.set_xlabel('Xhat and X')
    ax1.set_ylabel('Y')
    ax1.set_title('{} - {} layers each with {} neurons'.format(name, self.num_layer, self.num_neuron))
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    show()

  def draw_mse(self, x, y):
    fig = figure(1)
    ax1 = fig.add_subplot(111)
    ax1.plot(x,y)
    ax1.set_ylabel('MSE')
    ax1.set_xlabel('Number of layers')
    show()

def get_mse(predictions, actual):
	return np.mean((predictions-actual)**2)

if __name__ == '__main__':
	  #part (a)
	  P=10
	  N=1
	  Y_train, X_train = generate_data_AWGN(10000,P,N)
	  Y_test,   X_test = generate_data_AWGN(1000, P,N)
	  layers = 1
	  neurons = 100
	  m1 = model(layers,neurons)
	  sess = m1.construct_graph()
	  m1.train(sess, (Y_train, X_train))
	  X_pred = m1.test(sess, Y_test)
	  X_lmmse_est = P/(P+N)*Y_test
	  m1.draw(Y_test, X_pred, X_test, X_lmmse_est, 'AWGN')
	  print("[Gaussian] Estimator using NN: MSE=",get_mse(X_pred,X_test))
	  print("[Gaussian] Estimator using Linear Estimator: MSE=",get_mse(X_lmmse_est,X_test))
	  print("[Gaussian] NN had", layers, "hidden layers and", neurons, "neurons per hidden layer.")
	  mse = {}
	  for layers in range(1,6):
	      m1 = model(layers,neurons)
	      sess = m1.construct_graph()
	      m1.train(sess, (Y_train, X_train))
	      X_pred = m1.test(sess, Y_test)
	      mse[layers] = get_mse(X_pred,X_test)
	  plt.plot(list(mse.keys()),[mse[m] for m in mse],'.-',label='Neural network MSE',markersize=14)
	  plt.plot(list(mse.keys()),[get_mse(X_lmmse_est,X_test) for m in mse],label='Linear MSE')
	  plt.ylim([.8,1])
	  plt.legend()
	  plt.ylabel('MSE')
	  plt.xlabel('Hidden layers')
	  plt.title('MSE of Estimator for Gaussian Data')
	  plt.show()
	  
	  
	  #part (b)
	  lambda_x =  .5
	  lambda_z = 1.0
	  b_x = 1/lambda_x
	  b_z = 1/lambda_z
	  Y_train, X_train = generate_data_laplace(10000,lambda_x,lambda_z)
	  Y_test,   X_test = generate_data_laplace(1000, lambda_x,lambda_z)
	  P = 2*(b_x)**2
	  N = 2*(b_z)**2	 
	  layers = 4
	  neurons = 100
	  m2 = model(layers,neurons)
	  sess = m2.construct_graph()
	  m2.train(sess, (Y_train, X_train))
	  X_pred = m2.test(sess, Y_test)
	  X_lmmse_est = P/(P+N)*Y_test
	  m2.draw(Y_test, X_pred, X_test, X_lmmse_est, 'Laplace')
	  print("[Laplace] Estimator using NN: MSE=",get_mse(X_pred,X_test))
	  print("[Laplace] Estimator using Linear Estimator: MSE=",get_mse(X_lmmse_est,X_test))
	  print("[Laplace] NN had", layers, "hidden layers and", neurons, "neurons per hidden layer.")
	  mse = {}
	  for layers in range(1,6):
	      m2 = model(layers,neurons)
	      sess = m2.construct_graph()
	      m2.train(sess, (Y_train, X_train))
	      X_pred = m2.test(sess, Y_test)
	      mse[layers] = get_mse(X_pred,X_test)
	  plt.plot(list(mse.keys()),[mse[m] for m in mse],'.-',label='Neural network MSE',markersize=14)
	  plt.plot(list(mse.keys()),[get_mse(X_lmmse_est,X_test) for m in mse],label='Linear MSE')
	  plt.ylim([.8,1])
	  plt.legend()
	  plt.ylabel('MSE')
	  plt.xlabel('Hidden layers')
	  plt.title('MSE of Estimator for Laplace Data')
	  plt.show()

	  # part (c)
	  P=4
	  N=.5
	  Y_train, X_train = generate_data_unif_cube(10000,P,N)
	  Y_test,   X_test = generate_data_unif_cube(1000, P,N)
	  layers = 5
	  neurons = 200
	  m3 = model(layers,neurons)
	  sess = m3.construct_graph()
	  m3.train(sess, (Y_train, X_train))
	  X_pred = m3.test(sess, Y_test)
	  m3.draw(Y_test, X_pred, X_test, None, 'Uniform')
	  print("[Uniform] Estimator using NN: MSE=",get_mse(X_pred,X_test))
	  print("[Uniform] NN had", layers, "hidden layers and", neurons, "neurons per hidden layer.")
	  mse = {}
	  for layers in range(1,15):
	      m3 = model(layers,neurons)
	      sess = m3.construct_graph()
	      m3.train(sess, (Y_train, X_train))
	      X_pred = m3.test(sess, Y_test)
	      mse[layers] = get_mse(X_pred,X_test)
	  plt.plot(list(mse.keys()),[mse[m] for m in mse],'.-',label='Neural network MSE',markersize=14)
	  plt.ylim([.8,1])
	  plt.legend()
	  plt.ylabel('MSE')
	  plt.xlabel('Hidden layers')
	  plt.title('MSE of Estimator for Cubic Uniform Data')
	  plt.show()
  
  
  

