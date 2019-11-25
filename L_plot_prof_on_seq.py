import time
import tempfile
import os
import pkgutil
import collections

import numpy as np
import pandas as pd
import peakutils
import statsmodels.api as sm


import matplotlib.pyplot as plt
from matplotlib import rc,gridspec,rcParams
from matplotlib.widgets import Slider, Button, CheckButtons, RadioButtons

from scipy.optimize import minimize 
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
from scipy.ndimage.filters import gaussian_filter as gfilt


from PIL import Image, ImageDraw, ImageFont



##Functions to plot profiles on sequence.
#######################

def plot_prof_on_seq(data=None,data_file=None,DNAseq,prof_names=None,\
			graphshow=False,pngfileout=None,title='',prof_columns="Intensity",seq_column="Site",rescale=None,zero_at=0,ylab=None,colors=None,\
			colorb={'A':'#0b0','T':'#0b0','G':'#00b','C':'#00b'},colorf={'A':'#fafafa','T':'#fafafa','G':'#fafafa','C':'#fafafa'},\
			plot_options={'linewidth':1.0,'markersize':8.0,'figsize':(12,3),'fontsize':None,'legendloc':'upper right'}):
	"""
	Function makes a plot of some values along the DNA sequence.
	data_file - csv file with a dataframe.
	prof_columns - list of column names to plot.
	rescale - None - no normalization,
		'every' - every profile separtely devide by its max value, (Y/Ymax),
		'every_min_max' - (Y-Ymin)/(Ymax-Ymin) for every profile,
		'together' - divide all the profiles by highest among max values,
		'together_min_max' - (Y-Ymin)/(Ymax-Ymin), where min and max are global for all profiles,
		'fit' - a linear fit between profiles without intersect will be made,
			and then they will be rescaled by the highest max value (as in 'together' mode).
		'fit_min_max' - a normal linear fit between profiles (with intersect) will be made
			and then they will be rescaled to [0,1] segement from global min to global max (same as together_min_max). 
	colors - list of colors for plottin several datasets.
	colorb - colors for DNA bases
	"""
	if plot_options['fontsize']:
		rcParams.update({'font.size': plot_options['fontsize']})

	#Let's make an image of a sequence
	temp = tempfile.TemporaryFile()
	im=Image.new('RGB',(24*len(DNAseq),40))
	aspect=(24.0*len(DNAseq))/40.0
	draw = ImageDraw.Draw(im)
	
	#Open font from package in a tricky way, independent of package installation mode
	#temp2 = tempfile.TemporaryFile()
	#fontfile = pkgutil.get_data('hydroid', 'pkgdata/cnrb.otf')
	#temp2.write(fontfile)
	#temp2.seek(0)
	font = ImageFont.truetype('cnrb.otf', 40)
	####

	for i,l in zip(range(len(DNAseq)),DNAseq):
		draw.rectangle(((24*i,0),(24*(i+1),40)),fill=colorb[l],outline=colorb[l])
		draw.text((24*i+1, 0), l,fill=colorf[l],font=font)
	im.save(temp, "PNG")
	# im.save('temp.png', "PNG")
	temp.seek(0)

	prof_columns=[prof_columns] if isinstance(prof_columns,basestring) else list(prof_columns)
	prof_names=([prof_names] if isinstance(prof_names,basestring) else list(prof_names)) if prof_names else None
	colors=list(colors) if colors else None

  if(data==None):
	  data=pd.read_csv(data_file,comment='#')
	# print data
	yvalues=data[prof_columns].values

	if rescale=='together':
		yvalues=yvalues/np.nanmax(yvalues)
	elif rescale=='together_min_max':
		yvalues=yvalues-np.nanmin(yvalues)
		yvalues=yvalues/np.nanmax(yvalues)
	elif rescale=='every':
		yvalues=yvalues/np.nanmax(yvalues,axis=0)
	elif rescale=='every_min_max':
		yvalues=yvalues-np.nanmin(yvalues,axis=0)
		yvalues=yvalues/np.nanmax(yvalues,axis=0)
	elif rescale=='fit':
		coefs=[]
		for c,i in zip(prof_columns,range(len(prof_columns))):
			fitres = sm.OLS(yvalues[:,i], yvalues[:,0],missing='drop').fit()
			coefs.append(fitres.params[0])
		yvalues=yvalues/np.array(coefs)
		yvalues=yvalues/np.nanmax(yvalues)
	elif rescale=='fit_min_max':
		coefs=[]
		const=[]
		for c,i in zip(prof_columns,range(len(prof_columns))):
			fitres = sm.OLS(yvalues[:,i],sm.add_constant(yvalues[:,0], prepend=False) ,missing='drop').fit()
			coefs.append(fitres.params[0])
			const.append(fitres.params[1])
		yvalues=(yvalues-np.array(const))/np.array(coefs)
		yvalues=yvalues-np.nanmin(yvalues)
		yvalues=yvalues/np.nanmax(yvalues)

	seqimgheight=12/3/aspect*np.nanmax(yvalues)*1.05
	vlinestart=seqimgheight/(np.nanmax(yvalues)*1.05+seqimgheight)

	fig=plt.figure(figsize=plot_options['figsize'])
	ax1 = plt.subplot()

	ind=np.array(map(lambda x:int(x[0:-1]),data[seq_column].values))-zero_at

	ax1.yaxis.grid(linewidth=0.5)
	#Let's draw minor grid ticks
	ticklabs=[]
	for ls in range(0-zero_at,len(DNAseq)-zero_at):
		if ls%5==0:
			ax1.axvline(ls,linewidth=1 if ls%10==0 else 0.5,color='black',ymin=vlinestart)
		if ls%10==0:
			ticklabs.append(ls)

	for c,i in zip(prof_columns,range(len(prof_columns))):
		ax1.plot(ind, yvalues[:,i], '.b-',color= colors[i] if colors else None,label=prof_names[i] if prof_names else c,linewidth=plot_options['linewidth'],markersize=plot_options['markersize'])
	
	ax1.set_xticks(ticklabs)
	ax1.set_ylabel((ylab if ylab else 'Value')+(", Normalized (%s)"%rescale if rescale else ""))
	ax1.set_title(title)
	ax1.set_ylim((-seqimgheight,np.nanmax(yvalues)*1.05))
	ax1.set_xlim((0.5-zero_at,len(DNAseq)+0.5-zero_at))
	ax1.set_yticks(ax1.get_yticks()[np.array(ax1.get_yticks())>=0])
	# legend=ax1.legend(loc="upper right",shadow=False)
	legend=ax1.legend(loc=plot_options['legendloc'],shadow=False)
	# legend=ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	legend.get_frame().set_facecolor('#FFFFFF')
	legend.get_frame().set_alpha(1)

	seqimg = Image.open(temp)
	# im = OffsetImage(seqimg, zoom=0.26, resample=True)
	# ab = AnnotationBbox(im, (1,1), xycoords='data',\
	# xybox=(0.5,0),box_alignment=(0., 0.),frameon=False)
	# ax1.add_artist(ab)
	# ax1.imshow(seqimg,interpolation='nearest', aspect=aspect,extent=(0.5,len(DNAseq)+0.5,0.0,0.2))
	# ax1.imshow(seqimg,interpolation='nearest', aspect=2000)
	ax1.imshow(seqimg,interpolation='bilinear',aspect='auto',extent=(0.5-zero_at,len(DNAseq)+0.5-zero_at,-seqimgheight,0.0))
	


	fig.tight_layout()

	if(pngfileout):
		plt.savefig(os.path.splitext(pngfileout)[0]+'.png',dpi=(300))
	if(graphshow):
		plt.show()
	plt.close(fig)
	temp.close()
	temp2.close()
