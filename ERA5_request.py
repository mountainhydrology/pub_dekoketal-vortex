#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 10:42:12 2019

@author: kok00042
"""

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'variable':[
            'temperature','u_component_of_wind','v_component_of_wind',
            'vertical_velocity'
        ],
        'pressure_level':[
            '200','225','250',
            '300','350','400',
            '450','500','550'
        ],
        'year':'1994',
        'month':[
            '06','07','08'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'area'          : [40, 60, 50, 70], # North, West, South, East.
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    },
    '1994_box1.nc')