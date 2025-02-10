'''
takes output from analyze_single_removals_real_vs_fake.py and graphs observed
in black and expected in pale blue
'''
import cPickle
import sys
import matplotlib as mpl
mpl.use('Agg') #comment this in if using a nongraphical interface, comment out to enable graphical interface
import matplotlib.pyplot as plt

data_folder=sys.argv[1]

sim_file=open(data_folder+'fake_distributions', 'rb')
real_file=open(data_folder+'real_distributions', 'rb')

def add_to_plot(score_list, the_color):
	x_list=[float(item_number)/len(score_list) for item_number, item in enumerate(score_list)]
#	x_list=[float(item_number) for item_number, item in enumerate(score_list)]	
	plt.plot(x_list, score_list, color=the_color)
	plt.plot(x_list, [0 for number in range(len(score_list))], color=[1.0, 0, 0])
	plt.margins(0,0)
	return x_list, score_list

def add_points(point_dict, x_list, y_list, eightmer, stat, rep):
	if eightmer not in point_dict:
		point_dict[eightmer]={}
	if stat not in point_dict[eightmer]:
		point_dict[eightmer][stat]={}
	point_dict[eightmer][stat][rep]=[x_list, y_list]
	return point_dict

def gather_stats(input_file):
	try:
		stat_list=cPickle.load(input_file)
		eightmer, stat, rep, gene_scores=stat_list
		scores=[score[0] for score in sorted(gene_scores)]
		return eightmer, stat, rep, scores
	except EOFError:
		return False

def get_new_point(point_one, point_two, x_value):
#	print point_one, point_two, x_value
	distance=x_value-point_one[0]
	fraction=distance/(point_two[0]-point_one[0])
	new_y=point_one[1]+fraction*(point_two[1]-point_one[1])
	return new_y

def calculate_new_points(real_xlist, real_ylist, fake_xlist, fake_ylist, nested_fakes):
	for point_number, real_x in enumerate(real_xlist):
		for fake_number, fake_x in enumerate(fake_xlist):
			if fake_x==real_x:
				nested_fakes[point_number].append(fake_ylist[fake_number])
				break
			elif fake_x>real_x:
				nested_fakes[point_number].append(get_new_point([fake_xlist[fake_number-1], fake_ylist[fake_number-1]], [fake_x, fake_ylist[fake_number]], real_x))
				break
	return nested_fakes

real_dict={}
output=1
while output:
	output=gather_stats(real_file)
	if output:
		eightmer, stat, rep, scores=output
		real_dict[(eightmer, stat)]=scores

output=1
#gain_limits=[[-0.02, 1.02],[-3,1]]
#loss_limits=[[-0.02, 1.02],[-12,12]]
gain_limits=[-3,8]
loss_limits=[-12,12]
point_dict={}

while output:
	output=gather_stats(sim_file)
	if output:
		eightmer, stat, rep, scores=output
		print 'graphing', eightmer, stat, rep
		x_list, y_list=add_to_plot(scores, [0.5, 0.5, 1.0])
		point_dict=add_points(point_dict, x_list, y_list, eightmer, stat, rep)
		if rep==99:
			x_list, y_list=add_to_plot(real_dict[(eightmer, stat)], [0.0, 0.0, 0.0])
			point_dict=add_points(point_dict, x_list, y_list, eightmer, stat, 'real')
			plt.title(eightmer)
			if stat=='gains':
				#plt.xlim(gain_limits[0])
				#plt.ylim(gain_limits[1])
				plt.ylim(gain_limits)
			if stat=='losses':
				#plt.xlim(loss_limits[0])
				#plt.ylim(loss_limits[1])
				plt.ylim(loss_limits)
			ax = plt.gca()
			ax.xaxis.set_tick_params(width=2)
			ax.yaxis.set_tick_params(width=2)
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('bottom')
			ax.tick_params(direction='out')
			plt.savefig(data_folder+eightmer+'_'+stat+'_for_figure.pdf')
			plt.clf()
cPickle.dump(point_dict, open(data_folder+'point_dict', 'wb'))

#point_dict=cPickle.load(open(data_folder+'point_dict', 'rb'))
output_file=open(data_folder+'summary_stats', 'w')
for eightmer in point_dict:
	print eightmer
	for stat in point_dict[eightmer]:
		real_xlist, real_ylist=point_dict[eightmer][stat]['real']
#		print 'real', real_xlist, real_ylist
		nested_fakes=[[] for number in real_xlist]
		for rep in point_dict[eightmer][stat]:
			if rep!='real':
				fake_xlist, fake_ylist=point_dict[eightmer][stat][rep]
				print 'rep is', rep, len(fake_xlist)#, fake_xlist, fake_ylist
				nested_fakes=calculate_new_points(real_xlist, real_ylist, fake_xlist, fake_ylist, nested_fakes)
#				print 'new fakes are', nested_fakes
		smaller_count, bigger_count=0,0
		for point_number, real_y in enumerate(real_ylist):
#			print real_y, nested_fakes[point_number]
			if real_y<min(nested_fakes[point_number]):
				smaller_count+=1
			if real_y>max(nested_fakes[point_number]):
				bigger_count+=1
		output_file.write('\t'.join(map(str, [eightmer, stat, smaller_count, bigger_count, len(real_xlist)]))+'\n')
