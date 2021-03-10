import cdsapi
import numpy as np

c = cdsapi.Client()

kmonth_ = [1]
# kday_ = np.linspace(13,29,17)
kday_ = np.asarray([14])
for k1 in range(3):
    kmonth = kmonth_[k1]
    for k2 in range(1):
        # kday = kday_[k2]
        kday = 16

        outfile = 'G:\ERA5\ERA5.single-level.2018%02d'%kmonth+'%02d'%kday+'.nc'
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': 'total_column_water_vapour',
                'year': '2018',
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