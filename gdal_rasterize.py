# -*- coding: utf-8 -*-
'''
credits to http://gis.stackexchange.com/users/2856/luke
in http://gis.stackexchange.com/questions/16837/how-can-i-turn-a-shapefile-into-a-mask-and-calculate-the-mean

Usage:
       gdal_rasterize.py [-p prototype_raster] [-u]
       [-at] [-burn value] | [-a attribute]
       [-where expression] | [-sql select_statement]
       [-of format] [-a_srs srs_def] [-co \"NAME=VALUE\"]*
       [-a_nodata value] [-init value]
       [-te xmin ymin xmax ymax] [-tr xres yres] [-ts width height]
       [-ot {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/
            CInt16/CInt32/CFloat32/CFloat64}]
       <src_datasource> <dst_filename>
Where:
    -p          Path to a prototype raster
    -u          Update the existing destination raster
    -at         Enables the ALL_TOUCHED rasterization option so that
                all pixels touched by lines or polygons will be
                updated not just those one the line render path, or
                whose center point is within the polygon. Defaults to
                disabled for normal rendering rules.
    -burn       A fixed value to burn into a band for all objects.
    -a          Identifies an attribute field on the features to be
                used for a burn in value.
    -where      An optional SQL WHERE style query expression to be
                applied to select features to burn in from the input
                layer.
    -sql        An SQL statement to be evaluated against the
                datasource to produce a virtual layer of features to
                be burned in.
    -of         The output format. Default is GeoTIFF (GTIFF). Use the
                short format name.
    -a_nodata   Assign a specified nodata value to output bands.
    -init       Pre-initialize the output image with this value.
                However, it is not marked as the nodata value in the output
                file.
    -a_srs      Override the projection for the output file. If not specified,
                the projection of the input vector file will be used if
                available. If incompatible projections between input and output
                files, no attempt will be made to reproject features. The
                srs_def may be any of the usual GDAL/OGR forms, complete WKT,
                PROJ.4, EPSG:n or a file containing the WKT.
    -co         Passes a creation option to the output format driver. Multiple
                -co options may be listed. See format specific documentation
                for legal creation options for each format.
    -te         set georeferenced extents. The values must be expressed in
                georeferenced units. If not specified, the extent of the output
                file will be the extent of the vector layers.
    -tr         set target resolution. The values must be expressed in
                georeferenced units. Both must be positive values.
    -tap        Align the coordinates of the extent of the output file to the
                values of the -tr, such that the aligned extent includes the
                minimum extent.
    -ts         Set output file size in pixels and lines. Note that -ts cannot
                be used with -tr
    -ot         For the output bands to be of the indicated data type.
                Defaults to Float64
    -w          Intermediate working format for rasters.
                Defaults to 'MEM'.
'''
#---Imports
from osgeo import gdal, gdalconst, ogr
import sys, os,  math, tempfile

gdal.UseExceptions()
ogr.UseExceptions()

def create_raster(xsize, ysize, driver='MEM',tmpfile='', gt=None, srs_wkt=None, nodata=None,  init=None, datatype=gdal.GDT_Byte):
    # Create a memory raster to rasterize into.
    out_ds = gdal.GetDriverByName(driver).Create(tmpfile, xsize, ysize, 1 ,datatype)
    if init is not None:out_ds.GetRasterBand(1).Fill(init)
    if nodata is not None:out_ds.GetRasterBand(1).SetNoDataValue(nodata)
    if gt:out_ds.SetGeoTransform(gt)
    if srs_wkt:out_ds.SetProjection(srs_wkt)
    return out_ds

def rasterize(src,dst,update=False,
              prototype=None,
              nodata=-9999,init=-9999,
              te=None,tr=None, ts=None,
              sql=None, where=None,
              out_format='GTIFF', out_type='Byte',
              out_srs=None, co=[], working_format='MEM',
              **kwargs):

    #Some examples
    #dstds=rasterize(src ,dst,ts=[500, 500], out_type='Int16',  burn_values=[37])
    #dstds=rasterize(src ,'',prototype=prot, nodata=-999, init=-999, ts=[500, 500], out_type='Int16', options=["ATTRIBUTE=test"],  out_format='MEM')
    #rasterize(src ,dst, tr=[0.01, 0.01], nodata=-999, init=-999, out_type='Int16', options=["ATTRIBUTE=test"])
    #rasterize(src ,dst, tr=[0.01, 0.01], options=["ATTRIBUTE=test"],  out_format='ECW',  co=['TARGET=10'])
    #rasterize(src ,dst, where='test=1', update=True, burn_values = [10])
    #rasterize(src ,dst, tr=[0.01, 0.01], te=[143.5, -43.5, 148.5, -39.5], out_type='Int16', options=["ATTRIBUTE=test"])

    #Temporary files
    if working_format=='MEM':tmpfd,tmpfile=[None,'']
    else:tmpfd,tmpfile=tempfile.mkstemp()

    #Open the vector layer to rasterize from
    src_ds = ogr.Open(src)
    if not src_ds: raise RuntimeError,'\'%s\' does not exist in the file system.' % src
    if sql:src_lyr=ds.ExecuteSQL(sql)
    else:
        src_lyr=src_ds.GetLayer()
        if where:src_lyr.SetAttributeFilter(where)
    xmin,xmax,ymin,ymax=src_lyr.GetExtent()
    src_ext=xmin,ymin,xmax,ymax
    try:
        src_wkt=src_lyr.GetSpatialRef().ExportToWkt()
        if not out_srs:out_srs=src_wkt
    except:
        out_srs='LOCAL_CS["arbitrary"]' # From http://svn.osgeo.org/gdal/trunk/autotest/alg/rasterize.py
    datatype=gdal.GetDataTypeByName(out_type)

    #Get a GDAL Dataset to rasterize into
    dst_driver=None
    if update: #Can we update an existing raster
        dstds = gdal.Open(dst, gdalconst.GA_Update)
    else:
        gdal.ErrorReset()
        if prototype:
            protds=gdal.Open(prototype)
            dst_driver=protds.GetDriver()
            dstds = create_raster(protds.RasterXSize, protds.RasterYSize, working_format,
                                  tmpfile,protds.GetGeoTransform(), protds.GetProjection(),
                                  datatype=datatype, nodata=nodata, init=init)
            del protds
        elif te and tr:
            out_gt=(te[0],tr[0],0,te[3],0, -tr[1])
            out_cols=int(math.ceil((te[2]-te[0])/tr[0]))
            out_rows=int(math.ceil((te[3]-te[1])/tr[1]))
            dstds = create_raster(out_cols, out_rows, working_format, tmpfile,
                                  out_gt, out_srs, datatype=datatype, nodata=nodata,
                                  init=init)
        elif te and ts:
            xres=(te[2]-te[0])/float(ts[0])
            yres=(te[3]-te[1])/float(ts[1])
            out_gt=(te[0],xres,0,te[3],0, -yres)
            dstds = create_raster(ts[0], ts[1], working_format, tmpfile,
                                  out_gt, out_srs, datatype=datatype, nodata=nodata,
                                  init=init)
        elif ts:
            xres=(src_ext[2]-src_ext[0])/float(ts[0])
            yres=(src_ext[3]-src_ext[1])/float(ts[1])
            out_gt=(src_ext[0],xres,0,src_ext[3],0, -yres)
            dstds = create_raster(ts[0], ts[1], working_format, tmpfile,
                                  out_gt, out_srs, datatype=datatype, nodata=nodata,
                                  init=init)
        elif tr:
            out_cols=int(math.ceil((src_ext[2]-src_ext[0])/float(tr[0])))
            out_rows=int(math.ceil((src_ext[3]-src_ext[1])/float(tr[1])))
            out_gt=(src_ext[0],tr[0],0,src_ext[3],0,-tr[1])
            dstds = create_raster(out_cols, out_rows, working_format, tmpfile,
                                  out_gt, out_srs, datatype=datatype, nodata=nodata,
                                  init=init)
        else:
            raise RuntimeError, '%s does not exist and neither a prototype raster nor appropriate options were specified!'%dst

        if dst_driver  is None:
            try:dst_driver = gdal.GetDriverByName(out_format)
            except:raise RuntimeError,'Format driver %s not found, pick a supported driver.' % out_format

    err = gdal.RasterizeLayer(dstds, [1], src_lyr, **kwargs)
    if err != 0:
        raise RuntimeError,  "error rasterizing layer: %s" % err
    if not update:
        dstds=dst_driver.CreateCopy(dst,dstds, 1, co)
        try:
            os.close(tmpfd)
            os.unlink(tmpfile)
        except:pass

    dstds.FlushCache()
    return dstds

if __name__ == '__main__':
    try:
        kwargs={}
        args=map(str,  sys.argv[1:])
        src,dst=args[-2:];del args[-2:]
        while args: #This is used instead of getopt or optparse
                    #to support gdal non standard long option name
                    #format "-option value" instead of "--option=value"
            arg=args.pop(0)
            if arg == '-p':
                kwargs['prototype']=args.pop(0)
            elif arg == '-u':
                kwargs['update']=True
            elif arg == '-sql':
                kwargs['sql']=args.pop(0)
            elif arg == '-w':
                kwargs['where']=args.pop(0)
            elif arg == '-at':
                kwargs['options'].append('ALL_TOUCHED=TRUE')
            elif arg == '-burn':
                kwargs['burn_values']=[int(args.pop(0))]
            elif arg == '-a':
                if kwargs.has_key(options):kwargs['options'].append()
                else:kwargs['options']=['ATTRIBUTE=%s'%args.pop(0)]
            elif arg == '-a_srs':
                kwargs['out_srs']=args.pop(0)
            elif arg == '-co':
                if kwargs.has_key(co):kwargs['co'].append(args.pop(0))
                else:kwargs['co']=[args.pop(0)]
            elif arg == '-a_nodata':
                kwargs['nodata']=args.pop(0)
            elif arg == '-init':
                kwargs['init']=args.pop(0)
            elif arg == '-te':
                kwargs['te']=map(float, args[:4]);del args[:4]
            elif arg == '-tr':
                kwargs['tr']=map(int, args[:2]);del args[:2]
            elif arg == '-ts':
                kwargs['ts']=map(int, args[:2]);del args[:2]
            elif arg == '-of':
                kwargs['out_format']=args.pop(0)
            elif arg == '-ot':
                kwargs['out_type']=args.pop(0)
            elif arg == '-w':
                kwargs['working_format']=args.pop(0)
            else:raise Exception
    except:
        print >> sys.stderr,__doc__
        sys.exit(1)

    dstds=rasterize(str(src) ,str(dst),  **kwargs)
