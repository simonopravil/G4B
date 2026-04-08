import ee


def Reducer1(Sentinel):
    # Define reducers
    median        = ee.Reducer.median().unweighted()
    sd          = ee.Reducer.stdDev().unweighted()
    percentiles = ee.Reducer.percentile([10, 25, 50, 75, 90]).unweighted()
    allMetrics  = median.combine(sd, sharedInputs=True).combine(percentiles, sharedInputs=True)
    stm = Sentinel.reduce(allMetrics)
    return stm
def Reducer2(Sentinel):
    # Define reducers
    median        = ee.Reducer.median().unweighted()
    sd          = ee.Reducer.stdDev().unweighted()
    allMetrics  = median.combine(sd, sharedInputs=True)
    stm = Sentinel.reduce(allMetrics)
    return stm

def prepareSTMs(Sentinel):
    Spring = (Sentinel
        .filter(ee.Filter.calendarRange(3,5,'month'))
        .select(['NDVI','NDMI','EVI', 'SWIR_ratio', 'MSAVI'],
        ['Spring_NDVI','Spring_NDMI','Spring_EVI', 'Spring_SWIR_ratio', 'Spring_MSAVI'])
    )
    Summer = (Sentinel
        .filter(ee.Filter.calendarRange(6,8,'month'))
        .select(['NDVI','NDMI','EVI', 'SWIR_ratio', 'MSAVI'],
        ['Summer_NDVI','Summer_NDMI','Summer_EVI', 'Summer_SWIR_ratio', 'Summer_MSAVI'])
    )
    Autumn = (Sentinel
        .filter(ee.Filter.calendarRange(9,10,'month'))
        .select(['NDVI','NDMI','EVI', 'SWIR_ratio', 'MSAVI'],
        ['Autumn_NDVI','Autumn_NDMI','Autumn_EVI', 'Autumn_SWIR_ratio', 'Autumn_MSAVI'])
    )
    SummerAutumn = (Sentinel
        .filter(ee.Filter.calendarRange(6,10,'month'))
        .select(['NDVI','NDMI','EVI', 'SWIR_ratio', 'MSAVI'],
        ['SummerAutumn_NDVI','SummerAutumn_NDMI','SummerAutumn_EVI', 'SummerAutumn_SWIR_ratio', 'SummerAutumn_MSAVI'])
    )
    Spring_Summer = (Sentinel
        .filter(ee.Filter.calendarRange(3,8,'month'))
        .select(['NDVI','NDMI','EVI', 'SWIR_ratio', 'MSAVI'],
        ['Spring_Summer_NDVI','Spring_Summer_NDMI','Spring_Summer_EVI', 'Spring_Summer_SWIR_ratio', 'Spring_Summer_MSAVI'])
    )

    stm = Reducer1(Sentinel)
    stms = ee.Image([stm, Reducer2(Spring), Reducer2(Summer), Reducer2(Autumn), Reducer2(SummerAutumn), Reducer2(Spring_Summer)]).multiply(100).toInt()
    return stms

def addAuxiliary(stms):
  # Load land cover image
  landcover = ee.Image('projects/g4bproject/assets/LC_con_RF').rename('landcover')
  # Define the landcover classes we're interested in
  classes = {
        1: 'urban',
        2: 'cropland',
        3: 'forest',
        4: 'shrub',
        5: 'grass',
        6: 'bare',
        7: 'water',
        8: 'wet',
  }
  # Create masks for each land cover class
  chm_data_dict = {
    'urban': landcover.eq(1).multiply(ee.Image.constant(1)),  # Assuming 1 is urban class
    'cropland': landcover.eq(2).multiply(ee.Image.constant(1)),  # Assuming 2 is cropland
    'grassland': landcover.eq(5).multiply(ee.Image.constant(1)),  # Assuming 3 is grassland
    'forest': landcover.eq(3).multiply(ee.Image.constant(1))  # Assuming 4 is forest
  }
  # Process each land cover class
  landcover_proximity_bands = []
  for name, chm_data in chm_data_dict.items():
    # Calculate proximity metrics for each class
    chm010m_mean = chm_data.reduceNeighborhood(
      reducer=ee.Reducer.sum(),
      kernel=ee.Kernel.circle(radius=1),
      #inputWeight0"mask"  #Removed extra quote
      inputWeight="mask"
    ).rename(f'{name}_prox10m')

    chm090m_mean = chm_data.reduceNeighborhood(
      reducer=ee.Reducer.sum(),
      kernel=ee.Kernel.circle(radius=9),
      #inputWeight0"mask" #Removed extra quote
      inputWeight="mask"
    ).rename(f'{name}_prox90m')

    chm360m_mean = chm_data.reduceNeighborhood(
      reducer=ee.Reducer.sum(),
      kernel=ee.Kernel.circle(radius=36),
      #inputWeight0"mask" #Removed extra quote
      inputWeight="mask"
    ).rename(f'{name}_prox360m')
    chm500m_mean = chm_data.reduceNeighborhood(
      reducer=ee.Reducer.sum(),
      kernel=ee.Kernel.circle(radius=50),
      #inputWeight0"mask" #Removed extra quote
      inputWeight="mask"
    ).rename(f'{name}_prox500m')

    # Add to bands list
    landcover_proximity_bands.extend([chm010m_mean, chm090m_mean, chm360m_mean,chm500m_mean])

  #terrain
  image = ee.Image()
  # terrain
  # Mosaic the collection to get a single global image
  alos = ee.ImageCollection("JAXA/ALOS/AW3D30/V4_1")
  elevation = alos.select('DSM').mosaic().rename('elevation')
  slope = ee.Terrain.slope(elevation).rename('slope')
  aspect = ee.Terrain.aspect(elevation).rename('aspect')
  # TPI
  focal_mean = elevation.focalMean(5, 'square')
  tpi = elevation.subtract(focal_mean).rename('tpi').multiply(100)

  solar_dni= ee.Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/global_solar_atlas/dni_LTAy_AvgDailyTotals').rename('solar').multiply(100) #longterm average of direct normal irradiation
  #climate taken from chelsa check  https://chelsa-climate.org/, described in technical specs: https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf
  clim_cmi=ee.Image("projects/g4bproject/assets/climate/cmi_mean").rename('clim_cmi')# Mean monthly climate moisture index
  clim_gf5=ee.Image("projects/g4bproject/assets/climate/gdgfgd5").rename('clim_gf5')# First growing degree day above 5°C
  clim_gd5=ee.Image("projects/g4bproject/assets/climate/gdd5").rename('clim_gd5')# Growing degree days heat sum above 5°C
  clim_b01=ee.Image("projects/g4bproject/assets/climate/bio1").rename('clim_b01')# mean annual air temperature
  clim_b12=ee.Image("projects/g4bproject/assets/climate/bio12").rename('clim_b12')# bio12 annual precipitation amount
  #soils
  soil_prd=ee.Image("projects/g4bproject/assets/soils/soil_cac").rename('soil_prd')#soil calcium
  soil_dth=ee.Image("projects/g4bproject/assets/soils/soil_dth").rename('soil_dth')#soil depth

  #nightlight, intactness, canopy height
  nightlight = (ee.ImageCollection('NOAA/VIIRS/DNB/ANNUAL_V22')
              .filter(ee.Filter.date('2022-01-01', '2023-01-01'))
              .select('median')
              .first()
              .rename('nightlight')
              .multiply(100)
          )
  chm = ee.Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1").rename('chm')
  chm_sd = ee.Image("users/nlang/ETH_GlobalCanopyHeightSD_2020_10m_v1").rename('chm_sd').multiply(100)
  #we switching this to the 10m resolution outcome
  chm_m = (ee.ImageCollection("projects/ee-chm-eu-2019/assets/planet_chm_2019") #meta CHM https://gee-community-catalog.org/projects/meta_trees/?h=canopy#dataset-citation
              .mosaic()
              .rename('chm_m'))


  smallTrees = chm_m.gt(10).And(chm_m.lt(20)).multiply(ee.Image.constant(1))
  bigTrees = chm_m.gt(20).multiply(ee.Image.constant(1))
  #proximity to forest

  smTno_009m = smallTrees.reduceNeighborhood(reducer=ee.Reducer.sum(),kernel=ee.Kernel.circle(3),inputWeight="mask").rename('smTno_009m').multiply(100);
  smTno_045m = smallTrees.reduceNeighborhood(reducer=ee.Reducer.sum(),kernel=ee.Kernel.circle(15),inputWeight="mask").rename('smTno_045m').multiply(100);
  smTno_090m = smallTrees.reduceNeighborhood(reducer=ee.Reducer.sum(),kernel=ee.Kernel.circle(30),inputWeight="mask").rename('SmTno_090m').multiply(100);
  bgTno_009m = bigTrees.reduceNeighborhood(reducer=ee.Reducer.sum(),kernel=ee.Kernel.circle(3),inputWeight="mask").rename('bgTno_009m').multiply(100);
  bgTno_045m = bigTrees.reduceNeighborhood(reducer=ee.Reducer.sum(),kernel=ee.Kernel.circle(15),inputWeight="mask").rename('bgTno_045m').multiply(100);
  bgTno_090m = bigTrees.reduceNeighborhood(reducer=ee.Reducer.sum(),kernel=ee.Kernel.circle(30),inputWeight="mask").rename('bgTno_090m').multiply(100);

  #chm009m_mean = chm.reduceNeighborhood(reducer=ee.Reducer.mean(),kernel=ee.Kernel.circle(radius=1)).rename('prox10m').multiply(100)
  #chm090m_mean = chm.reduceNeighborhood(reducer=ee.Reducer.mean(),kernel=ee.Kernel.circle(radius=9)).rename('prox90m').multiply(100)
  #chm360m_mean = chm.reduceNeighborhood(reducer=ee.Reducer.mean(),kernel=ee.Kernel.circle(radius=36)).rename('prox360m').multiply(100)
  #chm009m_stdD = chm.reduceNeighborhood(reducer=ee.Reducer.stdDev(),kernel=ee.Kernel.circle(radius=1)).rename('prox10m_std').multiply(100)
  #chm090m_stdD = chm.reduceNeighborhood(reducer=ee.Reducer.stdDev(),kernel=ee.Kernel.circle(radius=9)).rename('prox90m_std').multiply(100)
  #chm360m_stdD = chm.reduceNeighborhood(reducer=ee.Reducer.stdDev(),kernel=ee.Kernel.circle(radius=36)).rename('prox360m_std').multiply(100)
  return stms.addBands([elevation, slope, aspect, tpi,
                        solar_dni,landcover,clim_cmi, clim_gf5, clim_gd5, clim_b01, clim_b12,
                        #soil_prd, soil_dth,
                        nightlight, chm_sd,chm_m,
                        smTno_009m, smTno_045m, smTno_090m,bgTno_009m, bgTno_045m, bgTno_090m,
                        #chm009m_mean,chm090m_mean,chm360m_mean,chm009m_stdD,chm090m_stdD,chm360m_stdD,
                        landcover_proximity_bands
                        ]).toInt()

def addAuxiliary2(stms, mowing_ds):
    # 1. Convert the Collection into one Image with many bands
    # If the collection has 5 images with 6 bands each,
    # this creates 1 image with 30 bands.
    mowing_img = mowing_ds.toBands()

    # 2. Now apply Image methods (unmask and cast to Integer)
    mowing_median = mowing_ds.median()
    # 3. Add these many bands to your existing stms Image
    return stms.addBands([mowing_img, mowing_median])