import ee

originalNames = ['B2','B3','B4', 'B5', 'B6', 'B7','B8','B11','B12']
commonNames = ['blue', 'green', 'red', 'R1', 'R2', 'R3','nir','swir1', 'swir2']

def unique_values(collection, field):
    values = ee.Dictionary(collection.reduceColumns(ee.Reducer.frequencyHistogram(), [field]).get('histogram')).keys()
    return values

def daily_mosaics(imgs):

  def simplifyDate(img):
    d = ee.Date(img.get('system:time_start'))
    day = d.get('day')
    m = d.get('month')
    y = d.get('year')
    simpleDate = ee.Date.fromYMD(y,m,day)
    return img.set('simpleTime',simpleDate.millis())

  imgs = imgs.map(simplifyDate)
  days = unique_values(imgs,'simpleTime')

  def do_mosaic(d):
    d = ee.Number.parse(d)
    d = ee.Date(d)
    t = imgs.filterDate(d,d.advance(1,'day'))
    f = ee.Image(t.first())
    t = t.mosaic()
    t = t.set('system:time_start',d.millis())
    t = t.copyProperties(f)
    return t

  imgs = days.map(do_mosaic)

  return ee.ImageCollection.fromImages(imgs)

def addVariables(img):
  date = img.set({'DATE': ee.Date(img.get('system:time_start')).format('YYYY-MM-dd')})
  ndvi = img.normalizedDifference(['nir', 'red']).rename('NDVI')
  ndmi = img.normalizedDifference(['red', 'nir']).rename('NDMI')
  evi = img.expression('2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))',{'nir': img.select('nir'), 'red': img.select('red'), 'blue': img.select('blue')}).rename('EVI')
  swir_ratio = img.expression('swir2/swir1',{'swir2': img.select('swir2'), 'swir1': img.select('swir1')}).rename('SWIR_ratio')
  msavi = img.expression('(2 * nir + 1 - sqrt(pow((2 * nir + 1), 2) - 8 * (nir - red)) ) / 2', {'nir': img.select('nir'), 'red': img.select('red')}).rename('MSAVI')
  nbr1 = img.expression('(nir - swir2) / (nir + swir2)',{'nir': img.select('nir'), 'swir2': img.select('swir2')}).rename('NBR1')
  nbr2 = img.expression('(swir1 - swir2) / (swir1 + swir2)',{'swir2': img.select('swir2'), 'swir1': img.select('swir1')}).rename('SWIR_ratio').rename('NBR2')
  # We can add multiple more indiceis
  return img.addBands([ndvi,ndmi,evi, swir_ratio, msavi,nbr1,nbr2])

def getSentinel(aoi, start, end, cloudTresh = 80, clearTresh = 0.60):
    """
    Retrieves a Sentinel-2 image collection for a specified Area of Interest (AOI)
    within a date range, applies cloud filtering, and adds additional bands for
    further analysis.

    Parameters:
    aoi (ee.Geometry or ee.Feature): The Area of Interest (AOI) for which the
                                     Sentinel-2 images will be retrieved.
    start (str): The start year for filtering the images. Should be in 'YYYY' format.
    end (str): The end year for filtering the images. Should be in 'YYYY' format.
    cloudTresh (float): The threshold for filtering images based on the percentage
                        of cloud cover. Images with cloud coverage below this
                        percentage are selected.
    clearTresh (float): The threshold for filtering based on the cloud score. Used to create a mask.

    Returns:
    ee.ImageCollection: A processed and filtered collection of Sentinel-2 images
                        for the specified AOI and date range, with daily mosaics.
    """

    csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')

    Sentinel = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
                .filterBounds(aoi)
                .filterDate(f'{start}-03-01', f'{end}-11-01')
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudTresh))
                .linkCollection(csPlus, ['cs_cdf'])
                .map(lambda img:
                     img.addBands(img.select('cs_cdf').gte(clearTresh).rename('mask')))
                .select(originalNames, commonNames)
                .map(addVariables)
    )
    # Convert the filtered and processed Sentinel collection to daily mosaics.
    Sentinel = daily_mosaics(Sentinel)
    # Return the final processed Sentinel-2 image collection.
    return Sentinel


