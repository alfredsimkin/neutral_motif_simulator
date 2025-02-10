'''
A hybrid between original and version2. Like version2, uses a continuous
colormap. Like version 1, prints only points with data, and uses methods that
I can somewhat understand.
'''
print 'hello'

import random
import custom
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy
from matplotlib import colors
import cPickle
import bisect
import math

point_radius=5
distance=3 #defines the distance over which to look for additional points
data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/paper_stats/'
real_mirs={line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in open('/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best')}
stat_list=['sites', 'gain', 'loss']
sim_count_dict=cPickle.load(open(data_folder+'full_sim_counts', 'rb'))
real_count_dict=cPickle.load(open(data_folder+'real_counts', 'rb'))
real_mir_dict, sim_mir_dict={},{}
eightmer_list=custom.count_in_base('AAAAAAAA', 4, 'ACGTz')

def calculate_distance(first_point, second_point):
	x_distance=abs(first_point[0]-second_point[0])
	y_distance=abs(first_point[1]-second_point[1])
	return (x_distance**2+y_distance**2)**0.5

def get_good_points(first_point, distance):
	good_points=[]
	for x in range(-distance, distance+1):
		for y in range(-distance, distance+1):
			new_point=(first_point[0]+x, first_point[1]+y)
			dist=calculate_distance(first_point, new_point)
			if dist<=distance:
				good_points.append(new_point)
	return good_points

def calculate_density(first_point, points, distance):
	'''
	returns the number of points within distance units of first_point given a dictionary of points
	'''
	count=0
	good_points=get_good_points(first_point, distance)
	for point in good_points:
		if point in points:
			count+=points[point]
	return count

def make_point_dict(x_values, y_values):
	'''
	returns the number of points present for each coordinate pair in a graph
	'''
	point_dict={}
	for value_number, x_value in enumerate(x_values):
		point=(x_value, y_values[value_number])
		if point not in point_dict:
			point_dict[point]=0
		point_dict[point]+=1
	return point_dict

def expand_point(point_list, point_radius):
	import copy
	'''
	expands points to be visible because individual points are too small to see
	points are square for convenience (radius is not euclidean)
	'''
	new_list=copy.deepcopy(point_list)
	for y_number, y in enumerate(point_list):
		for x_number, x in enumerate(y):
			if x!=0:
				for new_y in range(y_number-point_radius, y_number+point_radius+1):
					for new_x in range(x_number-point_radius, x_number+point_radius+1):
						if new_list[new_y][new_x]==0:
							new_list[new_y][new_x]=x
	return new_list
			
lims=[[0,3000],[0,3000]]
#lims=[[3609,4019],[1307,1917]]
for stat_number in range(1,3):
	stat=stat_list[stat_number]
	rank_dict={}
	sim_count_list, real_count_list, sim_mir_list, real_mir_list=[],[],[],[]
	for eightmer in real_count_dict:
		sim_count_mean, sim_count_var=custom.meanvar(sim_count_dict[eightmer][stat_number])
		sim_count_mean=int(sim_count_mean)
		sim_count_list.append(sim_count_mean)
		real_count_list.append(real_count_dict[eightmer][stat_number])
		if eightmer in real_mirs:
			sim_mir_list.append(sim_count_mean)
			real_mir_list.append(real_count_dict[eightmer][stat_number])
	custom.scatter_it(sim_count_list, real_count_list, 'sim_'+stat_list[stat_number], 'real_'+stat_list[stat_number], [[0,3000],[0,3000]])
	point_dict=make_point_dict(sim_count_list, real_count_list)
	xmax=max(sim_count_list)
	ymax=max(real_count_list)
	density_list=[[0.0]*(xmax+1+point_radius) for point in range(ymax+1+point_radius)]
	print 'getting densities'
	output_file=open(stat_list[stat_number]+'_individual_points', 'w')
	for point in point_dict:
		x_point, y_point=point
		count=calculate_density(point, point_dict, distance)
		output_file.write('\t'.join(map(str, [x_point, y_point, count]))+'\n')
		density_list[y_point][x_point]=math.log(count, 2)+1 #format looks weird but plt.pcolormesh plots y coords as first coordinate

	plt.xlim(lims[0])
	plt.ylim(lims[1])
	ax = plt.gca()
	ax.xaxis.set_tick_params(width=2)
	ax.yaxis.set_tick_params(width=2)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(direction='out')

	print 'got densities'
	density_list=expand_point(density_list, point_radius)
	cmap=plt.cm.gist_earth_r
	print 'step 1'
	density_array=numpy.array(density_list)
	#cPickle.dump(density_array, open('density_array', 'wb'))
	#heatmap = plt.pcolormesh(density_array, cmap=cmap) #too slow with pdf output
	heatmap = plt.pcolormesh(density_array, cmap=cmap, rasterized=True, vmax=9)
	print 'step 2'
	plt.colorbar(heatmap)
	print 'step 3'
	plt.savefig(stat+'.pdf', dpi=1200, bbox_inches='tight')
	print 'step 4'
	plt.clf()
