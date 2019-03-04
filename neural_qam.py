from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import os

def generate_data(n = 10000, d=2):
  # onehot is the one-hot representation of H. 
  # For example, if H = 2, then onehot = [0,0,1,0]. 
  '''
  Fill in your implementation here.
  '''
  s = [[1,1],[-1,1],[-1,-1],[1,-1]]
  P = np.random.randint(0,4,size=(n,),dtype=np.int32)
  y = [np.array(s[i])*d/2 for i in P]
  y = np.array(y) + np.random.normal(scale=1,size=(n,2)) 
  onehot = np.zeros((n,4))
  for (i,j) in enumerate(P):
    onehot[i,j] = 1
  return y, onehot


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
    W = tf.get_variable("weights_%d"%i, 
      [input_shape[-1], num_neuron], initializer = tf.truncated_normal_initializer(stddev = .1))
    b = tf.get_variable("bias_%d" %i, 
      [1, num_neuron], initializer = tf.truncated_normal_initializer(stddev = .1))
    y = tf.matmul(x, W) + b
    return activation(y)
  
  def construct_graph(self):
    tf.reset_default_graph()
    self.x = tf.placeholder(tf.float32, [None, 2])
    h = self.x
    '''
    Fill in your implementation here.
    '''
    #hidden layers
    for i in range(self.num_layer):
        h = self.layer(h,self.num_neuron,i,tf.nn.relu)
    #output layer
    h = self.layer(h,4,self.num_layer,tf.nn.softmax)
    y = h
    
    self.prediction = y
    self.y_ = tf.placeholder(tf.float32, [None, 4])
    self.loss = tf.reduce_mean(-tf.reduce_sum(self.y_ * tf.log(y), reduction_indices=[1]))
    self.train_step = tf.train.GradientDescentOptimizer(0.5).minimize(self.loss)
    
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
  d = 3
  Y_train, P_train = generate_data(10000,d)
  Y_test , P_test  = generate_data(1000 ,d)
  layers = 3
  neurons = 100
  m = model(layers, neurons)
  sess = m.construct_graph()
  m.train(sess, (Y_train, P_train))
  P_pred = m.test(sess, Y_test)
  S_test = np.argmax(P_test,axis=1)
  S_pred = np.argmax(P_pred,axis=1)
  accuracy = np.mean(S_test==S_pred)
  print("Sent s0:",np.sum(P_test[:,0]==1),"; predicted:",np.sum(S_pred==0))
  print("Sent s1:",np.sum(P_test[:,1]==1),"; predicted:",np.sum(S_pred==1))
  print("Sent s2:",np.sum(P_test[:,2]==1),"; predicted:",np.sum(S_pred==2))
  print("Sent s3:",np.sum(P_test[:,3]==1),"; predicted:",np.sum(S_pred==3))
  print("Layers %d, neurons %d: %.2f" %(layers,neurons,accuracy))
  Y_test , P_test  = generate_data(1000 ,d)
  S_test = np.argmax(P_test,axis=1)
  S_pred = np.apply_along_axis(lambda x: 2*(x[1]<0) + (x[0]>0 and x[1]<0) + (x[0]<0 and x[1]>0),1,Y_test)
  accuracy = np.mean(S_test==S_pred)
  print("Sent s0:",np.sum(P_test[:,0]==1),"; predicted:",np.sum(S_pred==0))
  print("Sent s1:",np.sum(P_test[:,1]==1),"; predicted:",np.sum(S_pred==1))
  print("Sent s2:",np.sum(P_test[:,2]==1),"; predicted:",np.sum(S_pred==2))
  print("Sent s3:",np.sum(P_test[:,3]==1),"; predicted:",np.sum(S_pred==3))
  print("Maximum likelihood: %.2f" %(accuracy))
  
  
  

  

