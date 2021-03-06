#!/usr/bin/env python3
#
# hdf5 tools for Illumina
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# January 2019

import MultiScaleData as MSD
import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.colors as _colors
import yaml as _yaml

def OpenCached(filename,cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds

def plot(ds,n_layer=None,log=False,**options):
    _plt.gca().set_aspect(1)

    for i,layer in reversed(list(enumerate(ds[:n_layer]))):
        n = ds[i].shape[0]
        buff = ds._attrs['layers'][i]['buffer']
        N = n/2 - buff
        ind = _np.arange(-N-1,N+1)+0.5
        X,Y = _np.meshgrid(ind,ind)

        if 'vmin' not in options:
            options['vmin'] = min(list(map(_np.min,ds)))
        if 'vmax' not in options:
            options['vmax'] = max(list(map(_np.max,ds)))

        if log:
            options['norm'] = _colors.LogNorm(
                vmin=options['vmin'],
                vmax=options['vmax']
            )

        psize = ds.pixel_size(i)/1000.
        _plt.pcolor(
            X*psize,
            Y*psize,
            layer[buff:n-buff,buff:n-buff][::-1],
            **options
        )

def from_domain(params,data=None):
    if isinstance(params,str):
        with open(params) as f:
            params = _yaml.safe_load(f)
    attrs = { k:v for k,v in list(params.items()) \
        if k not in ['extents','observers'] }
    attrs['obs_lat'] = [ d['latitude'] for d in params['observers'] ]
    attrs['obs_lon'] = [ d['longitude'] for d in params['observers'] ]
    attrs['obs_x'] = [ d['x'] for d in params['observers'] ]
    attrs['obs_y'] = [ d['y'] for d in params['observers'] ]
    attrs['layers'] = [ {k:v for k,v in list(d.items()) if k != 'layer'} \
        for d in params['extents'] ]
    return MSD.MultiScaleData(attrs,data)
