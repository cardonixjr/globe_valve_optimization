import matplotlib.pyplot as plt
import numpy as np

def myPlot(function, a,b, marks = []):
		
	fig = plt.figure()

	x_values = np.arange(a,b+(b-a)/25,(b-a)/25)
	y_values = [function(x) for x in x_values]

	#calcula o minimo e maximo da funcao
	minimo = maximo = function(x_values[0])
	for x in x_values:
		if function(x) < minimo: minimo = function(x)
		if function(x) > maximo: maximo = function(x)
	
	plt.xlim(a,b)
	plt.ylim(minimo, maximo)
	plt.plot(x_values, y_values)
	if marks != []:	
		for mark in marks:
			plt.scatter(mark[0], mark[1], s=16, c='red')
		plt.xlabel("x")
	plt.ylabel("f(x)")
	
	plt.draw()
	plt.pause(0.5)
	input()
	plt.close('all')

