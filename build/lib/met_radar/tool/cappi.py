# Module: Radar imager: contoured altitude plan position indicator(CAPPI) 
# Author: bo_fan@qq.com(Thanks for any bug report)
# Version: 20230915 v1.0

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib import cm,colormaps
from matplotlib.colors import ListedColormap,BoundaryNorm
import cmaps

class cappi(object):
    '''
    radarplot=cappi(xx,yy,cc,field,savepath)

    1.intro
    -----------------------------------------------------------
        plot radar variables, including ppi, ppi in map and so 
        on

    2.inputs
    -----------------------------------------------------------
        xx: gate_x after meshgrid <array, nx*ny> 
        yy: gate_y after meshgrid <array, nx*ny>
        cc: field variable <array, nx*ny>
        field: field name as 'REF', 'KDP', 'ZDR', 'VEL', 'RHO'
        savepath: optional, path for saving output figure 
                  such as: savepath=outdir+str(swp+1)+'_swp='+
                           str(sweep[swp])+'deg_'+field+'.png'

    3.outputs
    -----------------------------------------------------------
        plt.show or savefig in savepath
    '''
    
    def __init__(self,xx,yy,cc,field,savename=None):
        
        self.field_map={'REF':'reflectivity',
                        'dBZ':'reflectivity',
                        'KDP':'specific_differential_phase',
                        'ZDR':'differential_reflectivity',
                        'VEL':'velocity',
                        'CC':'cross_correlation_ratio',
                        'RHO':'cross_correlation_ratio',
                        'rhohv':'cross_correlation_ratio',}
        self.field_lim={'REF':[-10,70],
                        'dBZ':[-10,70],
                        'KDP':[-2,5],
                        'ZDR':[-1,8],
                        'VEL':[-30,30],
                        'CC':[0,1],
                        'RHO':[0,1],
                        'rhohv':[0,1],}
        self.field_levs={'REF':_perfect_levs(np.linspace(-10,70,17)),
                         'dBZ':_perfect_levs(np.linspace(-10,70,17)),
                         'KDP':_perfect_levs(np.linspace(-2,5,15)),
                         'ZDR':_perfect_levs(np.linspace(-1,8,19)),
                         'VEL':_perfect_levs(np.linspace(-30,30,13)),
                         'CC':[0,.1,.3,.5,.6,.7,.8,.85,.9,.92,.94,.95,.96,.97,.98,.99,1],
                         'RHO':[0,.1,.3,.5,.6,.7,.8,.85,.9,.92,.94,.95,.96,.97,.98,.99,1],
                         'rhohv':[0,.1,.3,.5,.6,.7,.8,.85,.9,.92,.94,.95,.96,.97,.98,.99,1],}
        self.field_cmap={'REF':radarCmap(self,'dBZ'),
                         'dBZ':radarCmap(self,'dBZ'),
                         'KDP':radarCmap(self,'KDP'),
                         'ZDR':radarCmap(self,'ZDR'),
                         'VEL':radarCmap(self,'VEL'),
                         'CC':radarCmap(self,'rhohv'),
                         'RHO':radarCmap(self,'rhohv'),
                         'rhohv':radarCmap(self,'rhohv'),}
        self.field_locs={'REF':self.field_levs['dBZ'],
                         'dBZ':self.field_levs['dBZ'],
                         'KDP':self.field_levs['KDP'],
                         'ZDR':self.field_levs['ZDR'],
                         'VEL':self.field_levs['VEL'],
                         'CC':np.linspace(0,1,len(self.field_levs['rhohv'])),
                         'RHO':np.linspace(0,1,len(self.field_levs['rhohv'])),
                         'rhohv':np.linspace(0,1,len(self.field_levs['rhohv'])),}
        self.xx,self.yy,self.cc,self.field,self.savename=xx,yy,cc,field,savename
        

    def draw(self):

        xx,yy,cc,field,savename=self.xx,self.yy,self.cc,self.field,self.savename
        plt.close()
        fig,ax1=plt.figure(figsize=(4, 3.5)),plt.subplot(111)
        im=plt.pcolormesh(xx,yy,cc,cmap=self.field_cmap[field][0],
                          norm=self.field_cmap[field][1]) 
                        #   vmin=field_lim[field][0],vmax=field_lim[field][1])
        # im=plt.contourf(xx,yy,cc,cmap=field_cmap[field][0],
        #                 norm=field_cmap[field][1],levels=field_levs[field]) 
        # 获取ax1的边界框
        bbox = ax1.get_position()

        # 计算相对比例
        width = bbox.width      # axes的宽度
        height = bbox.height    # axes的高度

        # 使用axes宽度的比例
        pad_ratio = 0.005        # 间距占axes宽度的2%
        cbar_width_ratio = 0.02  # 颜色条宽度占axes宽度的3%

        cax = fig.add_axes([
            bbox.x1 + pad_ratio * width,  # 右边界 + 2%宽度作为间距
            bbox.y0,                      # 底部对齐
            cbar_width_ratio * width,     # 颜色条宽度为axes宽度的3%
            bbox.height                   # 与axes同高
        ])
        
        cb=ColorbarBase(cax,cmap=self.field_cmap[field][0],
                        norm=self.field_cmap[field][1],ticks=self.field_levs[field])
        if savename is not None:
            plt.savefig(savename, dpi=600)
        return  ax1,cax

def radarCmap(self,field):

    # match different radar field variables to different colormaps
    if self.field_map[field]==self.field_map['dBZ']:
        c1=[255,255,255,256]
        c2=[0,236,236,256]
        c3=[1,160,246,256]
        c4=[0,0,246,256]
        c5=[0,255,0,256]
        c6=[0,200,0,256]
        c7=[0,144,0,256]
        c8=[255,255,0,256]
        c9=[231,192,0,256]
        c10=[255,144,0,256]
        c11=[255,0,0,256]
        c12=[214,0,0,256]
        c13=[192,0,0,256]
        c14=[221,94,253,256]
        c15=[191,2,238,256]
        c16=[85,1,106,256]
        c=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16]
        cmap = ListedColormap(np.array(c)/256)
        norm= BoundaryNorm(self.field_levs[field],cmap.N)
    elif self.field_map[field]==self.field_map['rhohv']:
        # cmap = cm.jet
        try:
            cmap=cmaps.Carbone42
        except:
            cmap=colormaps['Carbone42']
        colorslist=cmap(np.linspace(0,1,len(self.field_levs[field])-1)) #cmap(np.exp(np.linspace(0,2,len(field_levs[field])-1))/np.exp(2))
        cmap = ListedColormap(colorslist)
        norm = BoundaryNorm(self.field_levs[field],cmap.N)
        # cmap,norm= from_levels_and_colors(field_lev[field],colorslist)
    else:
        try:
            cmap=cmaps.Carbone42
        except:
            cmap=colormaps['Carbone42']
        colorslist=cmap(np.linspace(0,1,len(self.field_levs[field])-1)) 
        cmap = ListedColormap(colorslist)
        norm = BoundaryNorm(self.field_levs[field],cmap.N)
    return cmap,norm

def _perfect_levs(levs):
    # transfer levs(list or array) to a new list
    # to draw beatiful ticklabels, like ' 0 0.1 1' instead of '0.0 0.1 1.0' with ugly '.0' 
    levs_list=list(levs)
    for i in range(len(levs)):
        if levs[i]%1==0: # is integer
            levs_list[i]=int(levs[i])
    return levs_list

