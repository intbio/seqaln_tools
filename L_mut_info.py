#!/usr/bin/env python
"""
L_mut_info.py - library caclulate mutual information of two sequences
via R.
Sequences should be simple strigns of capital letters.
Only 20 amino acids + gaps (-) are allowed.

"""
__author__="Alexey Shaytan"

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com

def mi(a,b):
    """
    a - sequnece as string
    b - sequence as string,
    returns mi in bits.
    """

    aa="ACDEFGHIKLMNPQRSTVWY-"
    ro.r('aalev=c("%s")'%'","'.join(aa))
    ro.r('comblev=apply(expand.grid(aalev,aalev),1,function(x) paste0(x[1],x[2]))')

    ro.r('x1=c("%s")'%'","'.join(a))
    ro.r('x2=c("%s")'%'","'.join(b))
    ro.r('library(entropy)')
    ro.r('countsx1=table(factor(x1,levels=aalev))')
    ro.r('countsx2=table(factor(x2,levels=aalev))')
    ro.r('entx1=entropy(countsx1,unit="log2")')
    ro.r('entx2=entropy(countsx2,unit="log2")')
    ro.r('p=paste0(x1,x2)')
    ro.r('pcounts=table(factor(p,levels=comblev))')
    ro.r('entcomb=entropy(pcounts,unit="log2")')
    mi=ro.r('mi=entx1+entx2-entcomb')[0]
    return mi


if __name__ == '__main__':

    a="KKKKKKKKAKKK"
    b="DDSDDDDDDDDD"
    print mi(a,b)
