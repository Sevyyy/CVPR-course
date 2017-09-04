#import packages
import numpy as np
import pandas as pd
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import struct

#function for load data
def load_image(filename):
	file = open(filename, 'rb')
	buf = file.read()
	head = struct.unpack_from('>IIII', buf, 0)
	n = head[1]
	row = head[2]
	col = head[3]
	offset = struct.calcsize('>IIII')
	bit = n*row*col
	bits = '>' + str(bit) + 'B'
	image = struct.unpack_from(bits, buf, offset)
	file.close()
	image = np.reshape(image, [n, row*col])
	return image

def load_label(filename):
	file = open(filename, 'rb')
	buf = file.read()
	head = struct.unpack_from('>II', buf, 0)
	n = head[1]
	offset = struct.calcsize('>II')
	bits = '>' + str(n) + 'B'
	label = struct.unpack_from(bits, buf, offset)
	file.close()
	label = np.reshape(label, [n,1])
	return label

def load_mnist_data():
	train_image_file = "data/train-images.idx3-ubyte"
	train_label_file = "data/train-labels.idx1-ubyte"
	test_image_file = "data/t10k-images.idx3-ubyte"
	test_label_file = "data/t10k-labels.idx1-ubyte"

	#load the train image
	train_image = load_image(train_image_file)
	train_label = load_label(train_label_file)
	test_image = load_image(test_image_file)
	test_label = load_label(test_label_file)

	return train_image, train_label, test_image, test_label


def main():
	#load my data
	data = np.loadtxt("./data/mynum", delimiter=",")
	print data.shape
	#load mnist data
	train_image, train_label, test_image, test_label = load_mnist_data()
	#construct model
	adb = AdaBoostClassifier(learning_rate = 0.6)
	adb = adb.fit(train_image, train_label)
	predictions = adb.predict(test_image)
	#count the correct answer
	count = 0
	for i in range(predictions.size):
		if predictions[i] == test_label[i]:
			count += 1
	print count
	print count*1.0/10000
	#test my number data
	pred = adb.predict(data)
	print pred

if __name__ == "__main__":
	main()

