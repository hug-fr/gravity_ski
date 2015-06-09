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

#sys.setrecursionlimit(5000)

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
	new_pt = []
	if height[ref_p] > min:
	
	#First tests to see if D8 matrix is inside raster
		if ref_p[0] - 1 >= 0:
			rmin = ref_p[0] - 1
		else:
			rmin = 0
		if ref_p[0] + 2 <= w:
			rmax = ref_p[0] + 2
		else:
			rmax = w
		if ref_p[1] - 1 >= 0:
			cmin = ref_p[1] - 1
		else:
			cmin = 0
		if ref_p[1] + 2 <= h:
			cmax = ref_p[1] + 2
		else:
			cmax = h
		
		#D8 comparison level 1	
		ref_mat = height[rmin:rmax,cmin:cmax]
		tmptarget = target[rmin:rmax,cmin:cmax]
		tmpsize = tmptarget.shape
		targetzero = numpy.zeros((tmpsize[0], tmpsize[1]), numpy.int8)
		targetzero = targetzero + 1 * numpy.logical_and(numpy.logical_and((ref_mat - height[ref_p] < 0),tmptarget==0), ref_mat >= min)
		# print ref_p, height[ref_p]
		# print rmin, rmax, cmin, cmax
		# print tmptarget, ref_mat
		# print tmpsize
		# raw_input()
		
		#Fill target et retrieve new points
		for n in numpy.argwhere(targetzero==1):
			p = (n[0] + rmin, n[1] + cmin)
			target[p] = 1
			new_pt.append(p)
		
		tmptarget = tmptarget + targetzero
		#First test to see if level 2 D8 matrix is inside raster
		for n in numpy.argwhere(tmptarget==0):
			p = (n[0] + rmin, n[1] + cmin)
			if p[0] - 1 >= 0:
				rmin2 = p[0] - 1
			else:
				rmin2 = 0
			if p[0] + 2 <= w:
				rmax2 = p[0] + 2
			else:
				rmax2 = w
			if p[1] - 1 >= 0:
				cmin2 = p[1] - 1
			else:
				cmin2 = 0
			if p[1] + 2 <= h:
				cmax2 = p[1] + 2
			else:
				cmax2 = h
			#D8 comparison level 2
			if height[p] <= height[ref_p] + threshold:
				ref_mat2 = height[rmin2:rmax2,cmin2:cmax2]
				tmptarget2 = target[rmin2:rmax2,cmin2:cmax2]
				tmpsize = tmptarget2.shape
				targetzero = numpy.zeros((tmpsize[0], tmpsize[1]), numpy.int8)
				targetzero = targetzero + 1 * numpy.logical_and(numpy.logical_and((ref_mat2 - height[p] < 0), tmptarget2 == 0), ref_mat2 >= min)
				
				#Fill target and retrieve new points
				if 1 in targetzero:
					target[p] = 1
					for o in numpy.argwhere(targetzero==1):
						q = (rmin2 + o[0], cmin2 + o[1])
						target[q] = 1
						new_pt.append(q)
	return new_pt
	
def detect_up(ref_p, height, threshold, max, target):
	new_pt = []
	if height[ref_p] < max:
	
	#First tests to see if D8 matrix is inside raster
		if ref_p[0] - 1 >= 0:
			rmin = ref_p[0] - 1
		else:
			rmin = 0
		if ref_p[0] + 2 <= w:
			rmax = ref_p[0] + 2
		else:
			rmax = w
		if ref_p[1] - 1 >= 0:
			cmin = ref_p[1] - 1
		else:
			cmin = 0
		if ref_p[1] + 2 <= h:
			cmax = ref_p[1] + 2
		else:
			cmax = h
		
		#D8 comparison level 1	
		ref_mat = height[rmin:rmax,cmin:cmax]
		tmptarget = target[rmin:rmax,cmin:cmax]
		tmpsize = tmptarget.shape
		targetzero = numpy.zeros((tmpsize[0], tmpsize[1]), numpy.int8)
		targetzero = targetzero + 1 * numpy.logical_and(numpy.logical_and((ref_mat - height[ref_p] > 0),tmptarget==0), ref_mat <= max)
		# print ref_p, height[ref_p]
		# print rmin, rmax, cmin, cmax
		# print tmptarget, ref_mat
		# print tmpsize
		# raw_input()
		
		#Fill target et retrieve new points
		for n in numpy.argwhere(targetzero==1):
			p = (n[0] + rmin, n[1] + cmin)
			target[p] = 1
			new_pt.append(p)
		
		tmptarget = tmptarget + targetzero
		#First test to see if level 2 D8 matrix is inside raster
		for n in numpy.argwhere(tmptarget==0):
			p = (n[0] + rmin, n[1] + cmin)
			if p[0] - 1 >= 0:
				rmin2 = p[0] - 1
			else:
				rmin2 = 0
			if p[0] + 2 <= w:
				rmax2 = p[0] + 2
			else:
				rmax2 = w
			if p[1] - 1 >= 0:
				cmin2 = p[1] - 1
			else:
				cmin2 = 0
			if p[1] + 2 <= h:
				cmax2 = p[1] + 2
			else:
				cmax2 = h
			#D8 comparison level 2
			if height[p] >= height[ref_p] - threshold:
				ref_mat2 = height[rmin2:rmax2,cmin2:cmax2]
				tmptarget2 = target[rmin2:rmax2,cmin2:cmax2]
				tmpsize = tmptarget2.shape
				targetzero = numpy.zeros((tmpsize[0], tmpsize[1]), numpy.int8)
				targetzero = targetzero + 1 * numpy.logical_and(numpy.logical_and((ref_mat2 - height[p] > 0), tmptarget2 == 0), ref_mat2 <= max)
				
				#Fill target and retrieve new points
				if 1 in targetzero:
					target[p] = 1
					for o in numpy.argwhere(targetzero==1):
						q = (rmin2 + o[0], cmin2 + o[1])
						target[q] = 1
						new_pt.append(q)
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
	
# #create final shapefile
# driver = ogr.GetDriverByName('ESRI Shapefile')
# if os.path.exists('ds_gravitaires.shp'):
	# driver.DeleteDataSource('ds_gravitaires.shp')
# skiarea = driver.CreateDataSource('ds_gravitaires.shp')
# srs = osr.SpatialReference()
# srs.ImportFromEPSG(2154)
# ds = skiarea.CreateLayer('ds_gravitaires', srs)
# ds.CreateField(ogr.FieldDefn('ind', ogr.OFTString))

#Connect to DB
myconn = psycopg2.connect("host="+conn_param.host+" dbname="+conn_param.dbname+" user="+conn_param.user+" password="+conn_param.password)
resort=myconn.cursor()

#create table to store rm geom
resort.execute("""
drop table if exists stations.geo_dsa_building;
create table stations.geo_dsa_building (
	ind varchar(10),
	id_rm varchar(255),
	frombot_geom geometry,
	fromtop_geom geometry)
""")
myconn.commit()

#load resorts extent
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
	#raw_input()

	# load starting points
	cur=myconn.cursor()
	query = """
	select gid, indicatif_station,
	st_x(st_startpoint(the_geom)) as sx, st_y(st_startpoint(the_geom)) as sy,
	st_x(st_endpoint(the_geom)) as ex, st_y(st_endpoint(the_geom)) as ey 
	from stations.geo_rm_sta_alpes_ind
	where indicatif_station = %s
	order by gid;"""
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
		
		sp = (rm[0], srow, scol)
		ep = (rm[0], erow, ecol)
		
		#print sp[1:3],height[sp[1:3]], ep[1:3],height[ep[1:3]]

		if in_image(sp[1:3], w, h) and in_image(ep[1:3], w, h):
			if height[sp[1:3]] > height[ep[1:3]]:
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
		hb.append(height[p[1:3]])
	minh = min(hb)

	hh = []
	for p in pth:
		hh.append(height[p[1:3]])
	maxh = max(hh)
	
	print "points sorted"

	threshold = 5 # allow for sligth descent from bottom or ascent from top

	#ski area computation
	count=0
	
	for i in range(0,len(pth)):
		id_rm = pth[i][0]
		p = (pth[i][1],pth[i][2])
		# initialise result image where pixels fromtop will be tagged to True
		fromtop=numpy.zeros((w, h)).astype(numpy.int)		
		ski_slope_down(p, height, threshold, minh, fromtop)
		
		#create tmp raster
		driver = gdal.GetDriverByName("GTiff")
		skirast = driver.Create('test.tif', maxcol-mincol, maxrow-minrow, 1, gdal.GDT_Byte)	
		skirast.SetGeoTransform((imgx[(0,mincol)], rastinit[1], 0, imgy[(minrow,0)], 0, rastinit[5]))
		srs = osr.SpatialReference()
		srs.ImportFromEPSG(2154)
		skirast.SetProjection(srs.ExportToWkt())
		skirast.GetRasterBand(1).WriteArray(fromtop)
		
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
		for feat in tmplayer:
			geom = feat.GetGeometryRef()
			fromtop_wkt=geom.ExportToWkt()
			
		tmpshp = None
			
		# initialise result image where pixels in frombot will be tagged to True
		frombot=numpy.zeros((w, h)).astype(numpy.int)
		p = (ptb[i][1],pth[i][2])
		ski_slope_up(p, height, threshold, maxh, frombot)
		
		#create raster from ski_area array
		driver = gdal.GetDriverByName("GTiff")
		skirast = driver.Create('test.tif', maxcol-mincol, maxrow-minrow, 1, gdal.GDT_Byte)	
		skirast.SetGeoTransform((imgx[(0,mincol)], rastinit[1], 0, imgy[(minrow,0)], 0, rastinit[5]))
		srs = osr.SpatialReference()
		srs.ImportFromEPSG(2154)
		skirast.SetProjection(srs.ExportToWkt())
		skirast.GetRasterBand(1).WriteArray(frombot)
		
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
		#Get geom WKT
		tmplayer.SetAttributeFilter("val = 1")
		for feat in tmplayer:
			geom = feat.GetGeometryRef()
			frombot_wkt=geom.ExportToWkt()
		
		query = """
		insert into stations.geo_dsa_building
		select %s::varchar, %s::varchar, ST_GeomFromText(%s,2154), ST_GeomFromText(%s,2154)
		"""
		cur.execute(query,(sta[0],id_rm,frombot_wkt,fromtop_wkt,))
		myconn.commit()
	
		tmpshp = None
		
		print id_rm, " done"
# # display result
# plt.imshow(ski_area*height)
# for i in range(0,len(pth)):
	# plt.plot(pth[i][1], pth[i][0], 'ro', c='b')
	# plt.plot(ptb[i][1], ptb[i][0], 'ro', c='r')
	# plt.plot([pth[i][1], ptb[i][1]], [pth[i][0], ptb[i][0]],'-')
# plt.show()