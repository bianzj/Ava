import cdsapi
import numpy as np
import os
from subprocess import Popen

##############################################
##### download polar orbiting satellite data
#############################################
# k = 0
# target_dir = 'I:/vj1_raw'
# order_number = '501495362'
# while k<100:
#     p = Popen('wget -e robots=off -m -np -R .html,.tmp -nH -c --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/orders/'+order_number+'/" --header "Authorization: Bearer 61A02C42-F341-11EA-B719-FDECBABD0E0F" -P '+target_dir)
#     stdout, stderr = p.communicate()
#     k = k+1
#
# print('success!')


############################################
##### down ERA5 reanalysis data
############################################
#!/usr/bin/env python

c = cdsapi.Client()
kmonth_ = [1,2,3,4,5,6,7,8,9,10,11,12]
# kday_ = np.linspace(13,29,17)
kday_ = np.asarray([31,28,31,30,31,30,31,31,30,31,30,31])
for k1 in range(len(kmonth_)):
    kmonth = kmonth_[k1]
    for k2 in range(kday_[k1]):
        # kday = kday_[k2]
        kday = k2+1

        # '2m_temperature', 'skin_temperature', 'surface_solar_radiation_downward_clear_sky',
        # 'surface_solar_radiation_downwards', 'surface_thermal_radiation_downward_clear_sky',
        # 'surface_thermal_radiation_downwards',
        # 'uv_visible_albedo_for_diffuse_radiation', 'uv_visible_albedo_for_direct_radiation',

        outfile = 'F:/ERA5_2019_EVAP/ERA5.single-level.2019%02d'%kmonth+'%02d'%kday+'.nc'
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', 'snow_density',
                     'soil_temperature_level_1', 'surface_latent_heat_flux', 'surface_pressure',
                      'surface_sensible_heat_flux', 'temperature_of_snow_layer', 'volumetric_soil_water_layer_1',
                ],
                'year': '2019',
                'month': '%02d'%kmonth,
                'day': '%02d'%kday,
                'time': [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ],
                'format': 'netcdf',
            },
            outfile)