import numpy
import matplotlib.pyplot as plt
import matplotlib.image
import gdal
import osr
import ogr
import psycopg2
import os
import sys
import conn_param

sys.setrecursionlimit(5000)

# for a pixel p=(x,y) return a list of tuples of the 8 neighbors
def neighborhood(p):
    return [(p[0]-1, p[1]-1),\
            (p[0]-1, p[1]+0),\
            (p[0]-1, p[1]+1),\
            (p[0]+0, p[1]+1),\
            (p[0]+1, p[1]+1),\
            (p[0]+1, p[1]+0),\
            (p[0]+1, p[1]-1),\
            (p[0]+0, p[1]-1)]
			

# check if pixel p=(x,y) is in image (boud check)
def in_image(p,w,h):
    return p[0]>=0 and p[1]>=0 and p[0]<w and p[1]<h

def detect_down(ref_p, height, threshold, min, target):
	new_pt=[]
	nbh = neighborhood(ref_p)
	#new_pt.empty()
	for n in nbh:
		if in_image(n,w,h) and target[n]==0 and height[n] < height[ref_p] and height[n]>=min:
			target[n] = 1
			new_pt.append(n)
	for n in nbh:
		if in_image(n,w,h) and target[n]==0 and height[n] >= height[ref_p] and height[n] < height[ref_p]+threshold and height[n]>=min:
			nbh2 = neighborhood(n)
			for n2 in nbh2:
				if in_image(n2,w,h) and target[n]==0 and height[n2] <= height[n]:
					target[n] = 1
					new_pt.append(n)
	return new_pt
	
def detect_up(ref_p, height, threshold, max, target):
	new_pt=[]
	nbh = neighborhood(ref_p)
	#new_pt.empty()
	for n in nbh:
		if in_image(n,w,h) and target[n]==0 and height[n] >= height[ref_p] and height[n]<=max:
			target[n] = 1
			new_pt.append(n)
	for n in nbh:
		if in_image(n,w,h) and target[n]==0 and height[n] < height[ref_p] and  height[n] < height[ref_p]-threshold and height[n]<=max:
			nbh2 = neighborhood(n)
			for n2 in nbh2:
				if in_image(n2,w,h) and target[n]==0 and height[n2] >= height[n]:
					target[n] = 1
					new_pt.append(n)
	return new_pt

def ski_slope_down(ref_p, height, threshold, min, target):
	np = detect_down(ref_p, height, threshold, min, target)
	for p in np:
		ski_slope_down(p, height, threshold, min, target)
		
def ski_slope_up(ref_p, height, threshold, max, target):
	np = detect_up(ref_p, height, threshold, max, target)
	for p in np:
		ski_slope_up(p, height, threshold, max, target)

'''
test = numpy.array([(100,99,100,102,102), (102,102,98, 103, 103), (100,95,96,97,100), (94,95,96,97,98), (100,100,95,100,100)])

w,h=test.shape
testnew = numpy.zeros((w, h)).astype(numpy.bool)

sp=(0,2)
testnew[sp]=1

ski_slope_down(sp, test,2,90, testnew)

print test
print testnew
'''
# load complete raster
img = gdal.Open('C:\ds_test_data\ign_mnt25_alpes.tif')
band1 = img.GetRasterBand(1)

rastinit = img.GetGeoTransform()

#x,y geographic reference matrix
imgx=numpy.zeros((1,img.RasterXSize)).astype(numpy.float)
imgy=numpy.zeros((img.RasterYSize,1)).astype(numpy.float)

for i in range(0,imgx.shape[1]):
	imgx[0,i]=rastinit[0]+(i*rastinit[1])


for i in range(0,imgy.shape[0]):
	imgy[i,0]=rastinit[3]+(i*rastinit[5])
	
#create final shapefile
driver = ogr.GetDriverByName('ESRI Shapefile')
if os.path.exists('ds_gravitaires.shp'):
	driver.DeleteDataSource('ds_gravitaires.shp')
skiarea = driver.CreateDataSource('ds_gravitaires.shp')
srs = osr.SpatialReference()
srs.ImportFromEPSG(2154)
ds = skiarea.CreateLayer('ds_gravitaires', srs)
ds.CreateField(ogr.FieldDefn('ind', ogr.OFTString))

#Connect to DB
myconn = psycopg2.connect("host="+conn_param.host+" dbname="+conn_param.dbname+" user="+conn_param.user+" password="+conn_param.password)

#load resorts extent
resort=myconn.cursor()
resort.execute("""
with a as (
	select indicatif_station, st_buffer(st_envelope(st_union(the_geom)),2000, 'endcap=square') geom
	from stations.geo_rm_sta_alpes_ind
	group by indicatif_station
	)

select indicatif_station, st_xmin(geom) xmin, st_ymin(geom), st_xmax(geom) xmax, st_ymax(geom) ymax
from a
where indicatif_station in ('3811')
order by indicatif_station;""")

for sta in resort:
	print sta[0]
	
	for i in range(0,imgx.shape[1]):
		if sta[1]-imgx[0,i]>0 and sta[1]-imgx[0,i+1]<0:
			mincol = i
		if sta[3]-imgx[0,i]>0 and sta[1]-imgx[0,i+1]<0:
			maxcol = i
	for i in range(0,imgy.shape[0]):
		if imgy[i,0]-sta[4]>0 and imgy[i+1,0]-sta[4]<0:
			minrow = i
		if imgy[i,0]-sta[2]>0 and imgy[i+1,0]-sta[2]<0:
			maxrow = i

	height = band1.ReadAsArray(mincol, minrow, maxcol-mincol, maxrow-minrow)

	# get width and heigth of image
	w,h = height.shape
	
	print "raster extracted"
	print w, h
	raw_input()

	# load starting points
	cur=myconn.cursor()
	query = """
	select gid, indicatif_station,
	st_x(st_startpoint(the_geom)) as sx, st_y(st_startpoint(the_geom)) as sy,
	st_x(st_endpoint(the_geom)) as ex, st_y(st_endpoint(the_geom)) as ey 
	from stations.geo_rm_sta_alpes_ind
	where indicatif_station = %s
	order by gid
	limit 1;"""
	ind = sta[0]
	#print ind
	cur.execute(query, (ind,))

	#sort top and bottom points
	pth = []
	ptb = []

	for rm in cur:
		for i in range(0,imgx.shape[1]):
			if rm[2]-imgx[0,i]>0 and rm[2]-imgx[0,i+1]<0:
				scol = i-mincol
			if rm[4]-imgx[0,i]>0 and rm[4]-imgx[0,i+1]<0:
				ecol = i-mincol
		for i in range(0,imgy.shape[0]):
			if imgy[i,0]-rm[3]>0 and imgy[i+1,0]-rm[3]<0:
				srow = i -minrow
			if imgy[i,0]-rm[5]>0 and imgy[i+1,0]-rm[5]<0:
				erow = i-minrow
		
		sp = (srow, scol)
		ep = (erow, ecol)

		if in_image(sp, w, h) and in_image(ep, w, h):
			if height[sp] > height[ep]:
				pth.append(sp)
				ptb.append(ep)
			else:
				pth.append(ep)
				ptb.append(sp)
		else:
			print "decoupage tout pourri"
	
	#get min and max of ski area
	hb = []
	for p in ptb:
		hb.append(height[p])
	minh = min(hb)

	hh = []
	for p in pth:
		hh.append(height[p])
	maxh = max(hh)
	
	print "points sorted"

	threshold = 0 # allow for sligth descent from bottom or ascent from top

	# initialise result image where pixels in domain will be tagged to True
	fromtop=numpy.zeros((w, h)).astype(numpy.int)
	frombot=numpy.zeros((w, h)).astype(numpy.int)

	#ski area computation
	count=0
	
	for p in pth:
		ski_slope_down(p, height, threshold, minh, fromtop)
	
	print "down done"
	raw_input()
	
	for p in ptb:
		ski_slope_up(p, height, threshold, maxh, frombot)
	
	ski_area=frombot*fromtop
	
	#create raster from ski_area array
	driver = gdal.GetDriverByName("GTiff")
	skirast = driver.Create('test.tif', maxcol-mincol, maxrow-minrow, 1, gdal.GDT_Byte)	
	skirast.SetGeoTransform((imgx[(0,mincol)], rastinit[1], 0, imgy[(minrow,0)], 0, rastinit[5]))
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(2154)
	skirast.SetProjection(srs.ExportToWkt())
	skirast.GetRasterBand(1).WriteArray(ski_area)
	
	print "raster saved"
	
	#create tmp shapefile for poligonyzation
	driver = ogr.GetDriverByName('ESRI Shapefile')
	if os.path.exists('tmp.shp'):
		driver.DeleteDataSource('tmp.shp')
	tmpshp = driver.CreateDataSource('tmp.shp')
	tmplayer = tmpshp.CreateLayer('tmp', srs)
	tmplayer.CreateField(ogr.FieldDefn("val", ogr.OFTInteger))
	
	gdal.Polygonize(skirast.GetRasterBand(1), None, tmplayer, 0, [], None)
	skirast = None
	
	tmplayer.SetAttributeFilter("val = 1")

	defn = ds.GetLayerDefn()
	newds = ogr.Feature(defn)
	newds.SetField('ind', sta[0])
	i=0
	if tmplayer.GetFeatureCount()>1:
		end=tmplayer.GetFeatureCount()
		newgeom = "MULTIPOLYGON ("
		for feat in tmplayer:
			i=i+1
			geom = feat.GetGeometryRef()
			geomwkt=geom.ExportToWkt()
			if i==1:
				newgeom = newgeom+"("+geomwkt[9:len(geomwkt)-1]+")"
			else:
				newgeom = newgeom+", ("+geomwkt[9:len(geomwkt)-1]+")"
		newgeom = newgeom+")"
		#print newgeom
		geom = ogr.CreateGeometryFromWkt(newgeom)
	else:
		for feat in tmplayer:
			geom = feat.GetGeometryRef()
		
	newds.SetGeometry(geom)
	ds.CreateFeature(newds)
	
	tmpshp = None

# display result
plt.imshow(ski_area*height)
for i in range(0,len(pth)):
	plt.plot(pth[i][1], pth[i][0], 'ro', c='b')
	plt.plot(ptb[i][1], ptb[i][0], 'ro', c='r')
	plt.plot([pth[i][1], ptb[i][1]], [pth[i][0], ptb[i][0]],'-')
plt.show()