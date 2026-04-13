# Module: 2D freely pcolor  
# Author: bo_fan@qq.com(Thanks for any bug report)
# Version: 20251217 v2.0

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib import cm,colormaps
from matplotlib.colors import Normalize
import cmaps

def todraw(xx,yy,cc,cmap='jet',lim=None,savepath=None,ax=None):
    '''
    f=todraw(ax,xx,yy,cc,figsize=[4,3.5],cmap='jet',lim=None,savepath=None)

    1.intro
    -----------------------------------------------------------
        plot radar variables, including ppi, ppi in map and so 
        on

    2.inputs
    -----------------------------------------------------------
        xx: gate_x after meshgrid <array, nx*ny> 
        yy: gate_y after meshgrid <array, nx*ny>
        cc: field variable <array, nx*ny>
        cmap: colormap, default None
        lim: color limits, default None
        savepath: optional, path for saving output figure 
                  such as: savepath=outdir+str(swp+1)+'_swp='+
                           str(sweep[swp])+'deg_'+field+'.png'

    3.outputs
    -----------------------------------------------------------
        plt.show or savefig in savepath
    '''
    cmap=plt.get_cmap(cmap)
    plt.close()
    if ax is None:
        fig,ax=plt.figure(figsize=[4,3.5]),plt.subplot(111)
    else:
        plt.sca(ax)
        fig=plt.gcf()
    if lim is None:
        vmin,vmax=cc.min(),cc.max()
    else:
        vmin,vmax=lim[0],lim[1]
    im=plt.pcolormesh(xx,yy,cc,cmap=cmap,vmin=vmin,vmax=vmax)
    # 获取ax1的边界框
    bbox = ax.get_position()

    # 计算相对比例
    width = bbox.width      # axes的宽度
    height = bbox.height    # axes的高度

    # 使用axes宽度的比例
    pad_ratio = 0.01        # 间距占axes宽度的2%
    cbar_width_ratio = 0.05  # 颜色条宽度占axes宽度的3%

    cax = fig.add_axes([
        bbox.x1 + pad_ratio * width,  # 右边界 + 2%宽度作为间距
        bbox.y0,                      # 底部对齐
        cbar_width_ratio * width,     # 颜色条宽度为axes宽度的3%
        bbox.height                   # 与axes同高
    ])
    cb=ColorbarBase(cax,cmap=cmap,norm=Normalize(vmin=im.norm.vmin,vmax=im.norm.vmax))
    if savepath is not None:
        plt.savefig(savepath, dpi=600)
    return ax,cax

