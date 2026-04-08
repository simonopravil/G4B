import ee

def getLandcoverBands():
  landcover = ee.Image('projects/g4bproject/assets/LC_con_RF').rename('landcover')

  masks = {
    'urban': landcover.eq(1).multiply(ee.Image.constant(1)),
    'cropland': landcover.eq(2).multiply(ee.Image.constant(1)),
    'grassland': landcover.eq(5).multiply(ee.Image.constant(1)),
    'forest': landcover.eq(3).multiply(ee.Image.constant(1))
  }

  radii = [
    ('prox10m', 1),
    ('prox90m', 9),
    ('prox360m', 36),
    ('prox500m', 50)
  ]

  proximity_bands = []
  for name, mask in masks.items():
    for suffix, radius in radii:
      band = mask.reduceNeighborhood(
        reducer=ee.Reducer.sum(),
        kernel=ee.Kernel.circle(radius=radius),
        inputWeight="mask"
      ).rename(f'{name}_{suffix}')
      proximity_bands.append(band)

  return landcover, proximity_bands


def getTerrainBands():
  alos = ee.ImageCollection("JAXA/ALOS/AW3D30/V4_1")
  elevation = alos.select('DSM').mosaic().rename('elevation')
  slope = ee.Terrain.slope(elevation).rename('slope')
  aspect = ee.Terrain.aspect(elevation).rename('aspect')
  focal_mean = elevation.focalMean(5, 'square')
  tpi = elevation.subtract(focal_mean).rename('tpi').multiply(100)

  return [elevation, slope, aspect, tpi]


def getSolarBand():
  return ee.Image(
    'projects/earthengine-legacy/assets/projects/sat-io/open-datasets/global_solar_atlas/dni_LTAy_AvgDailyTotals'
  ).rename('solar').multiply(100)


def getClimateBands():
  clim_cmi = ee.Image("projects/g4bproject/assets/climate/cmi_mean").rename('clim_cmi')
  clim_gf5 = ee.Image("projects/g4bproject/assets/climate/gdgfgd5").rename('clim_gf5')
  clim_gd5 = ee.Image("projects/g4bproject/assets/climate/gdd5").rename('clim_gd5')
  clim_b01 = ee.Image("projects/g4bproject/assets/climate/bio1").rename('clim_b01')
  clim_b12 = ee.Image("projects/g4bproject/assets/climate/bio12").rename('clim_b12')

  return [clim_cmi, clim_gf5, clim_gd5, clim_b01, clim_b12]


def getSoilBands():
  soil_prd = ee.Image("projects/g4bproject/assets/soils/soil_cac").rename('soil_prd')
  soil_dth = ee.Image("projects/g4bproject/assets/soils/soil_dth").rename('soil_dth')

  return [soil_prd, soil_dth]


def getNightlightBand():
  return (
    ee.ImageCollection('NOAA/VIIRS/DNB/ANNUAL_V22')
      .filter(ee.Filter.date('2022-01-01', '2023-01-01'))
      .select('median')
      .first()
      .rename('nightlight')
      .multiply(100)
  )


def getCanopyAndTreeBands():
  chm_sd = ee.Image("users/nlang/ETH_GlobalCanopyHeightSD_2020_10m_v1") \
    .rename('chm_sd') \
    .multiply(100)

  chm_m = (
    ee.ImageCollection("projects/ee-chm-eu-2019/assets/planet_chm_2019")
      .mosaic()
      .rename('chm_m')
  )

  smallTrees = chm_m.gt(10).And(chm_m.lt(20)).multiply(ee.Image.constant(1))
  bigTrees = chm_m.gt(20).multiply(ee.Image.constant(1))

  smTno_009m = smallTrees.reduceNeighborhood(
    reducer=ee.Reducer.sum(),
    kernel=ee.Kernel.circle(3),
    inputWeight="mask"
  ).rename('smTno_009m').multiply(100)

  smTno_045m = smallTrees.reduceNeighborhood(
    reducer=ee.Reducer.sum(),
    kernel=ee.Kernel.circle(15),
    inputWeight="mask"
  ).rename('smTno_045m').multiply(100)

  smTno_090m = smallTrees.reduceNeighborhood(
    reducer=ee.Reducer.sum(),
    kernel=ee.Kernel.circle(30),
    inputWeight="mask"
  ).rename('SmTno_090m').multiply(100)

  bgTno_009m = bigTrees.reduceNeighborhood(
    reducer=ee.Reducer.sum(),
    kernel=ee.Kernel.circle(3),
    inputWeight="mask"
  ).rename('bgTno_009m').multiply(100)

  bgTno_045m = bigTrees.reduceNeighborhood(
    reducer=ee.Reducer.sum(),
    kernel=ee.Kernel.circle(15),
    inputWeight="mask"
  ).rename('bgTno_045m').multiply(100)

  bgTno_090m = bigTrees.reduceNeighborhood(
    reducer=ee.Reducer.sum(),
    kernel=ee.Kernel.circle(30),
    inputWeight="mask"
  ).rename('bgTno_090m').multiply(100)

  tree_proximity_bands = [
    smTno_009m, smTno_045m, smTno_090m,
    bgTno_009m, bgTno_045m, bgTno_090m
  ]

  return chm_sd, chm_m, tree_proximity_bands


def addAuxiliary(stms):
  landcover, landcover_proximity_bands = getLandcoverBands()
  terrain_bands = getTerrainBands()
  solar_dni = getSolarBand()
  climate_bands = getClimateBands()
  soil_bands = getSoilBands()
  nightlight = getNightlightBand()
  chm_sd, chm_m, tree_proximity_bands = getCanopyAndTreeBands()

  return stms.addBands(
    terrain_bands +
    [solar_dni, landcover] +
    climate_bands + soil_bands +
    [nightlight, chm_sd, chm_m] +
    tree_proximity_bands +
    landcover_proximity_bands
  ).toInt()



def addManagement(stms):
  mowing = ee.ImageCollection(
    'projects/ee-simonopravil/assets/DisertationProject/Mowing/Mowing_Intensity_Alps_Carpathians'
  )

  def mask_no_grass(img):
    return img.updateMask(img.select('b7').neq(0)).updateMask(img.neq(999))

  mowing_ds = (
    mowing
      .map(mask_no_grass)
      .select(['b1', 'b6'], ['DOY1', 'n_mowing'])
      .filter(ee.Filter.stringStartsWith('system:index', 'Carpathians'))
  )

  mowing_img = mowing_ds.toBands()
  mowing_median = mowing_ds.median()

  return stms.addBands([mowing_img, mowing_median])