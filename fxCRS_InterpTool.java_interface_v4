package com.lundin.prowess.tool.fxCRS_Interp;

import org.javaseis.properties.DataDomain;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Arrays;
import java.util.List; 
import java.util.ArrayList;
import java.util.Collections;
import java.lang.Math;

import com.lgc.gpr.util.ParameterSet;
import java.lang.Math;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import com.lundin.tools.ParamsMethods.*;
import edu.mines.jtk.la.DMatrix;
import edu.mines.jtk.util.ArrayMath;
import com.lundin.tools.Matlab;

import com.lgc.prodesk.geometry.Geometry;
import com.lgc.gpr.util.Format;
import com.lgc.gpr.util.ParameterSet;
import com.lgc.prodesk.hdr.HdrEntry;
import com.lgc.prodesk.hdr.CatalogedHdrEntry;
import com.lgc.prodesk.hdr.HdrCatalog;
import com.lgc.prowess.hdr.HdrCatalogImpl;
import com.lgc.prodesk.seisdata.DataContext;
import com.lgc.prowess.exec.JobContext;
import com.lgc.prowess.hdr.Hdr;
import com.lgc.prowess.seisdata.SeisData;
import com.lgc.prowess.tool.InitPhaseException;
import com.lgc.prowess.tool.ToolContext;
import com.lgc.prowess.tool.SimpleTool;

import java.io.IOException;
import com.lgc.prowess.tool.ExecPhaseException;
import com.lgc.gpr.util.ParameterNotFoundException;
import com.lgc.prowess.seisdata.SeisUtil;
import org.javaseis.grid.BinGrid;
import org.javaseis.grid.IGridAccessor;
import org.javaseis.util.SeisException;
import com.lgc.prodesk.hdr.DimensionHdr;
import com.lgc.prowess.javaseis.smart.SmartFramework;
import com.lgc.prowess.seisdata.SeismicDataset;
import com.lgc.prowess.tool.DatasetParms;
import com.lgc.prowess.seisdata.SeisUtil;
import com.lgc.prowess.seisdata.DataContextModifier;
import com.lgc.prowess.seisdata.DataContextImpl;
import org.javaseis.properties.AxisLabel;
import com.lgc.prowess.util.StringArrays;

public class fxCRS_InterpTool extends SimpleTool {


  // Use a log instead of stdout and stderr
  private static final Logger LOG = Logger.getLogger("com.lundin.prowess.tool.sagc");
  private float _minu;
  private float _maxu;
  private float _maxstr;
  private float _maxspatstr;
  private float _maxtrace;
  private float _min_offset;
  private float _max_offset;
  private float _step_offset;
  
  private int _output_mode;
  private int _ref;
  
  private int _apt_iline;
  private int _apt_xline;
  private int _apt_offset_min;
  private int _apt_offset_max;
  private double[] _apt_offset_contp;
  private double[] _apt_offset_val_contp;
  
  
  // define the geometry instance
  private Geometry _geometry; 
  private double _ref_point_x;
  private double _ref_point_y;
  private int    _ref_point_inline;
  private int    _ref_point_xline;
  private double _azimuth;
  private double _inline_step;
  private double _xline_step;
  private double _time_step;
 
  private CatalogedHdrEntry xlex;
  private CatalogedHdrEntry ilex;
  private CatalogedHdrEntry _dimensionHdrHypercube;
  private CatalogedHdrEntry _dimensionHdrVolume;
  private CatalogedHdrEntry _dimensionHdrFrame;
  private CatalogedHdrEntry _dimensionHdrTrace;
  
  private CatalogedHdrEntry _sidewayDimensionHdrHypercube;
  private CatalogedHdrEntry _sidewayDimensionHdrVolume;
  private CatalogedHdrEntry _sidewayDimensionHdrFrame;
  private CatalogedHdrEntry _sidewayDimensionHdrTrace;
  
  // define the new header entry used to flag the CRS interpolated trace 
  private CatalogedHdrEntry _catalogedHdrEntry_1; 
  private CatalogedHdrEntry _catalogedHdrEntry_2; 
   
   
  // define the necessary trace header catalog instances
  private CatalogedHdrEntry _QC_Header;
    
    
  // define the parameters instances
  private SeismicDataset _veldataset;
  private SeismicDataset _crsParadataset_A1;
  private SeismicDataset _crsParadataset_A2;
  private SeismicDataset _crsParadataset_B11;
  private SeismicDataset _crsParadataset_B12;
  private SeismicDataset _crsParadataset_B22;
  
  // define the smart framework instance for velocity/parameter volumesmode
  private SmartFramework _vel_fw;
  
  // define the range of parameter file
  private int[] _inline_Vel = new int[2];
  private int[] _xline_Vel  = new int[2];  
  private int[] _sample_Vel = new int[2]; 
  
  // define super-gather inline/xline
  private int _IL_SG;
  private int _XL_SG;
  
  
  // define if display the QC log
  int qc_log_flag = 0;
  
  
  @Override
  public DataContext init(JobContext jobContext, ToolContext toolContext,
      DataContext dataContext) throws InitPhaseException {
    
    //----------------------------------------------------------------------------------------------------------
    // create new header to flag CRS interpolation generated traces
    //----------------------------------------------------------------------------------------------------------
    	
    // create a new data context in order to save the new header
    DataContextImpl newDataContext 	= 	new DataContextImpl(dataContext);
    
    // create the new header catalog 
    HdrCatalogImpl  newHdrCatalog   	= 	newDataContext.getHdrCatalogImpl();
    
    //create the new header (entry)
    String hdrName_1 			= 	"CRS_flag";
    HdrEntry hdrEntry_1 		= 	HdrEntry.makeEntry(hdrName_1, "Flag of CRS interpolated trace", 1, HdrEntry.INT_FORMAT);
    _catalogedHdrEntry_1 		= 	newHdrCatalog.addEntry(hdrEntry_1);
    
    String hdrName_2 			= 	"CRS_type";
    HdrEntry hdrEntry_2 		= 	HdrEntry.makeEntry(hdrName_2, "Flag of CRS output data type", 1, HdrEntry.INT_FORMAT);
    _catalogedHdrEntry_2 		= 	newHdrCatalog.addEntry(hdrEntry_2);	

    // Get the user parameters for this tool
    ParameterSet menuParms = toolContext.getMenuParms();
    
    _apt_iline       		= menuParms.getInt("APT_IL_PARM", 5);
    _apt_xline       		= menuParms.getInt("APT_XL_PARM", 5);
    //_apt_offset_min  		= menuParms.getInt("APT_OFFSET_MIN_PARM", 1);
    //_apt_offset_max  		= menuParms.getInt("APT_OFFSET_MAX_PARM", 5);  
    _output_mode 		= menuParms.getInt("OUTPUT_MODE_PARM", 1);
    
    _min_offset			= menuParms.getFloat("OFFSET_MIN_PARM", -3475.0f);
    _step_offset		= menuParms.getFloat("OFFSET_STEP_PARM", 50.0f);
    _max_offset			= menuParms.getFloat("OFFSET_MAX_PARM", 3475.0f);
    
    _apt_offset_contp		= StringArrays.stringToDoubleArray(menuParms.getString("_APT_OFFSET_CONTP_PARM", "1,70,71,140"));
    _apt_offset_val_contp	= StringArrays.stringToDoubleArray(menuParms.getString("_APT_OFFSET_CONTP_MAX_PARM", "5,5,5,5"));


    //HdrCatalog hdrCatalog = dataContext.getHdrCatalog();
    
    
   if (!newHdrCatalog.entryExists("smhxline")) {
      	throw new InitPhaseException("The input data is not supergathers!");
    }
   
    // convenience class for accessing header entries
    xlex 		= newHdrCatalog.getEntry("smhxline");
    ilex 		= newHdrCatalog.getEntry("smhiline");

    // get the geometry setting from datacontext
    _geometry 		= dataContext.getGeometry();
    _ref_point_x     	= _geometry.xILine1Start;
    _ref_point_y	= _geometry.yILine1Start;
    _ref_point_inline   = _geometry.minILine;
    _ref_point_xline    = _geometry.minXLine;
    _inline_step	= _geometry.dCdpXLine;
    _xline_step		= _geometry.dCdpILine;
    _azimuth		= _geometry.azimuth;
    _time_step          = dataContext.getSampleInterval();
    
    // convert the ms to s for the sample ratio (ex. 4ms -> 0.004s)
    if(_time_step  >=  1){
    	_time_step = _time_step * 0.001;
    }
    
    // convert the azimuth angle to the angle range (0-180degree) 
    if (_azimuth  <  0){
        _azimuth = 360 - Math.abs(_azimuth);
    }
    
    if(_azimuth  >=  180){
    	_azimuth = _azimuth - 180;
    }
    

    LOG.info(" -------- QC of geometry:  azimuth ------------_");
    LOG.info(" QC of geometry:  _ref_point_x ="       + _ref_point_x);
    LOG.info(" QC of geometry:  _ref_point_y ="       + _ref_point_y);
    LOG.info(" QC of geometry:  _ref_point_inline ="  + _ref_point_inline);
    LOG.info(" QC of geometry:  _ref_point_xline ="   + _ref_point_xline);    
    LOG.info(" QC of geometry:  _azimuth ="           + _azimuth);
    LOG.info(" QC of geometry:  _inline_step ="       + _inline_step);
    LOG.info(" QC of geometry:  _xline_step ="        + _xline_step);
    LOG.info(" QC of geometry:  _time_step ="         + _time_step);

    // check the input data dimensions
    String[] filePaths_vel         = DatasetParms.getFilePaths(menuParms, fxCRS_InterpProc.VEL_DATASET,jobContext.getAreaLinePath());
    String[] filePaths_crspara_A1  = DatasetParms.getFilePaths(menuParms, fxCRS_InterpProc.CRS_PARAMETER_A1,jobContext.getAreaLinePath());
    String[] filePaths_crspara_A2  = DatasetParms.getFilePaths(menuParms, fxCRS_InterpProc.CRS_PARAMETER_A2,jobContext.getAreaLinePath());
    String[] filePaths_crspara_B11 = DatasetParms.getFilePaths(menuParms, fxCRS_InterpProc.CRS_PARAMETER_B11,jobContext.getAreaLinePath());
    String[] filePaths_crspara_B12 = DatasetParms.getFilePaths(menuParms, fxCRS_InterpProc.CRS_PARAMETER_B12,jobContext.getAreaLinePath());
    String[] filePaths_crspara_B22 = DatasetParms.getFilePaths(menuParms, fxCRS_InterpProc.CRS_PARAMETER_B22,jobContext.getAreaLinePath());
    
    String filePath_vel 	 = (filePaths_vel.length > 0) ? filePaths_vel[0] : null;
    String filePath_crspara_A1   = (filePaths_crspara_A1.length  > 0) ? filePaths_crspara_A1[0]  : null;
    String filePath_crspara_A2   = (filePaths_crspara_A2.length  > 0) ? filePaths_crspara_A2[0]  : null;    
    String filePath_crspara_B11  = (filePaths_crspara_B11.length > 0) ? filePaths_crspara_B11[0] : null;
    String filePath_crspara_B12  = (filePaths_crspara_B12.length > 0) ? filePaths_crspara_B12[0] : null;        
    String filePath_crspara_B22  = (filePaths_crspara_B22.length > 0) ? filePaths_crspara_B22[0] : null;    
    
    
    if (filePath_vel == null || filePath_crspara_A1 == null || filePath_crspara_A2 == null || filePath_crspara_B11 == null || filePath_crspara_B12 == null || filePath_crspara_B22 == null) 
    throw new InitPhaseException("No dataset was selected");

    try {
    
      // load the velocity/parameters datasets based on the provided data paths 
      _veldataset         = SeisUtil.openDataset(filePath_vel);
      _crsParadataset_A1  = SeisUtil.openDataset(filePath_crspara_A1);
      _crsParadataset_A2  = SeisUtil.openDataset(filePath_crspara_A2);
      _crsParadataset_B11 = SeisUtil.openDataset(filePath_crspara_B11);
      _crsParadataset_B12 = SeisUtil.openDataset(filePath_crspara_B12);
      _crsParadataset_B22 = SeisUtil.openDataset(filePath_crspara_B22);
      
      // check the volume names of the velocity/parameters
      LOG.info("Dataset Descriptive Name="  + _veldataset.getDescName());
      LOG.info("CRS Para Descriptive Name=" + _crsParadataset_A1.getDescName());
      LOG.info("CRS Para Descriptive Name=" + _crsParadataset_A2.getDescName());
      LOG.info("CRS Para Descriptive Name=" + _crsParadataset_B11.getDescName());
      LOG.info("CRS Para Descriptive Name=" + _crsParadataset_B12.getDescName());      
      LOG.info("CRS Para Descriptive Name=" + _crsParadataset_B22.getDescName());      
      
      // Get the volume range of the velocity/parameters
      _vel_fw = _veldataset.getDataContext().getSmartFramework();  
      
      _inline_Vel[0]  = (int)_vel_fw.getMinFrame();  
      _inline_Vel[1]  = (int)_vel_fw.getMaxFrame();
       _xline_Vel[0]  = (int)_vel_fw.getMinTrace();  
       _xline_Vel[1]  = (int)_vel_fw.getMaxTrace();
      _sample_Vel[0]  = (int)_vel_fw.getMinSample(); 
      _sample_Vel[1]  = (int)_vel_fw.getMaxSample();
      
      LOG.info("The velocity volume's information: dimension= " + _vel_fw.getDataDimensions()); 
      LOG.info("The velocity volume's information: min-frame= " + _vel_fw.getMinFrame()); 
      LOG.info("The velocity volume's information: max-frame= " + _vel_fw.getMaxFrame()); 
      LOG.info("The velocity volume's information: min-trace= " + _vel_fw.getMinTrace()); 
      LOG.info("The velocity volume's information: max-trace= " + _vel_fw.getMaxTrace()); 
      
    } catch (IOException e) {
      throw new InitPhaseException(e);
    }
    
    //SmartFramework smartFramework1 = _veldataset.getDataContext().getSmartFramework(); 

    /* TODO: Is this check necessary? 
    SmartFramework smartFramework = _veldataset.getDataContext().getSmartFramework();
    if (smartFramework == null)
      throw new InitPhaseException("Dataset '" + filePath_vel
                                   + "' does not have a well-formed JavaSeis framework");
    */
    
    DimensionHdr[] inputDimensionHdrs = dataContext.getSeisGrid().getHdrs();

    // ***** CODE ADDED HERE *****
    _dimensionHdrTrace = newHdrCatalog.getEntry(inputDimensionHdrs[1].getHdrEntry());
    if (_dimensionHdrTrace == null)  // This should never happen.
      throw new InitPhaseException("Header entry '" + inputDimensionHdrs[1].getHdrEntry()
                                   + "' is missing from the input data");
                                   
                                   
    LOG.info("------Input prestack data dimensionHdrTrace key-----------" + _dimensionHdrTrace );  
    
    
    _dimensionHdrFrame = newHdrCatalog.getEntry(inputDimensionHdrs[2].getHdrEntry());
    if (_dimensionHdrFrame == null)  // This should never happen.
      throw new InitPhaseException("Header entry '" + inputDimensionHdrs[2].getHdrEntry()
                                   + "' is missing from the input data");
                                   
    LOG.info("------Input prestack data dimensionHdrFrame key-----------" + _dimensionHdrFrame );  

    // Note that we have to check to see if the dataset is more than 3-dimensional.
    if (inputDimensionHdrs.length > 3) {
      _dimensionHdrVolume = newHdrCatalog.getEntry(inputDimensionHdrs[3].getHdrEntry());
      if (_dimensionHdrVolume == null)
        throw new InitPhaseException("Header entry '" + inputDimensionHdrs[3].getHdrEntry()
                                     + "' is missing from the input data");
    } else {
      // This is a 3-dimensional dataset.  No problem.
      _dimensionHdrVolume = null;
    }
    
    LOG.info("------Input prestack data dimensionHdrVolume key-----------" + _dimensionHdrVolume ); 

    // Note that we have to check to see if the dataset is more than 4-dimensional.
    if (inputDimensionHdrs.length > 4) {
      _dimensionHdrHypercube = newHdrCatalog.getEntry(inputDimensionHdrs[4].getHdrEntry());
      if (_dimensionHdrHypercube == null)
        throw new InitPhaseException("Header entry '" + inputDimensionHdrs[4].getHdrEntry()
                                     + "' is missing from the input data");
    } else {
      // This is a 3-dimensional or 4-dimensional dataset.  No problem.
      _dimensionHdrHypercube = null;
    }
    
    LOG.info("------Input prestack data dimensionHdrHypercube key-----------" + _dimensionHdrHypercube ); 

    // ***** CODE ADDED HERE *****
    //==============================================================================================
    
    
    // No change to the input data context.        
 
   //return dataContext;
   return newDataContext;
  } // context


  @Override
  public SeisData exec(SeisData data) {

    // Get the volume/frame for this ensemble.  We will assume that the first trace is
    // representative.
    
    int[] hdr0 = data.getHdr(0);
    //int[][] hdr = data.getHdrBuffer();
    //int[] hdr0  = hdr[0];
    
    
    long hypercubeKey;
    if (_dimensionHdrHypercube == null) {
      hypercubeKey = SmartFramework.IMPLIED_LOGICAL_VALUE;
    } else {
      hypercubeKey = _dimensionHdrHypercube.getLongVal(hdr0);
    }

    long volumeKey;
    if (_dimensionHdrVolume == null) {
      volumeKey = SmartFramework.IMPLIED_LOGICAL_VALUE;
    } else {
      volumeKey = _dimensionHdrVolume.getLongVal(hdr0);
    }

    long frameKey;
    if (_dimensionHdrFrame == null) {
      frameKey = SmartFramework.IMPLIED_LOGICAL_VALUE;
    } else {
      frameKey = _dimensionHdrFrame.getLongVal(hdr0);
    }  
    
    long traceKey;
    if (_dimensionHdrTrace == null) {
      traceKey = SmartFramework.IMPLIED_LOGICAL_VALUE;
    } else {
      traceKey = _dimensionHdrTrace.getLongVal(hdr0);
    }    

    LOG.info("Input prestack data hypercubeKey key " + hypercubeKey);     
    LOG.info("Input prestack data volume key ="  + volumeKey); 
    LOG.info("Input prestack data frameKey key=" + frameKey); 
    LOG.info("Input prestack data traceKey key=" + traceKey); 
  
    
    // define variables to hold input CRS parameters 
    SeisData sidewaysVelData     = null;
    SeisData sidewaysCRSData_A1  = null;
    SeisData sidewaysCRSData_A2  = null;
    SeisData sidewaysCRSData_B11 = null;
    SeisData sidewaysCRSData_B12 = null;
    SeisData sidewaysCRSData_B22 = null;
    
    // get the indices of inline and xline for the velocity and CRS parameters, 
    // be aware of the main data are prestack data (4D) and velocity & CRS parameters are poststack volumes(3D) 
    long hypercubeKey_poststak 	= 0;
    long volumeKey_poststak 	= 0;
    long frameKey_poststak  	= volumeKey;   //inline index
    
    // read the sideway datas and output the frame objects 
    try {
      // only read selected frame, the whole inline
      sidewaysVelData     = _veldataset.readFrame(hypercubeKey_poststak, volumeKey_poststak , frameKey_poststak); 
      sidewaysCRSData_A1  = _crsParadataset_A1.readFrame(hypercubeKey_poststak,  volumeKey_poststak , frameKey_poststak); 
      sidewaysCRSData_A2  = _crsParadataset_A2.readFrame(hypercubeKey_poststak,  volumeKey_poststak , frameKey_poststak); 
      sidewaysCRSData_B11 = _crsParadataset_B11.readFrame(hypercubeKey_poststak, volumeKey_poststak , frameKey_poststak); 
      sidewaysCRSData_B12 = _crsParadataset_B12.readFrame(hypercubeKey_poststak, volumeKey_poststak , frameKey_poststak); 
      sidewaysCRSData_B22 = _crsParadataset_B22.readFrame(hypercubeKey_poststak, volumeKey_poststak , frameKey_poststak);
      
    } catch (IOException e) {
      LOG.warning("Error reading hypercube=" + hypercubeKey + " volume=" + volumeKey
                  + " frame=" + frameKey
                  + " in the sideways dataset");
                  
      sidewaysVelData = null;  // For clarity.
    }


    if (sidewaysVelData == null || sidewaysCRSData_A1 == null || sidewaysCRSData_A2 == null || sidewaysCRSData_B11 == null || sidewaysCRSData_B12 == null || sidewaysCRSData_B22 == null ) {
      // Usually that just means that this data is not present.
      
      LOG.info("================ Job debugging flag, one of parameters are not presented ================================");
      
      if (LOG.isLoggable(Level.FINE)) {
        LOG.fine("Unable to read hypercube=" + hypercubeKey + " volume=" + volumeKey
                 + " frame=" + frameKey
                 + " in the sideways dataset");
      }

    } else {
      // We have both data and sidewaydaa, proceed to process now.
      
      // ***** CODE ADDED HERE *****
      
	//----------------------------------------------------------------------------------------------------------
      	//  ***** get the prestack CRS gathers   *****
	//----------------------------------------------------------------------------------------------------------
	
	// get the prestack seismic data 
      	
      	float[][] data2d 	= data.getTraces();	  // get all the non-zero traces 
      	//float[][] data2d 	= data.getTraceBuffer();  // get all the traces which including the zero traces  	
      	

    	// get the seismic data dimensions (be aware that the data in promax is transposed)
        int ntr_gth 		= data2d.length;    // get the number of rows
        int ntimes_gth		= data2d[0].length; // get the number of coloumn
        
        // get the trace header of seismic data
        int[][] header_gth	= data.getHdrs();        // get all the non-zero traces' header
        //int[][] header_gth	= data.getHdrBuffer();   // get all the traces' header which including the zero traces  
        
        //get the SG-inline an SG-xline of current cmp
    	_IL_SG = data.getHdrCatalog().getEntry(HdrEntry.SG_ILINE).getIntVal(header_gth[0]);
    	_XL_SG = data.getHdrCatalog().getEntry(HdrEntry.SG_XLINE).getIntVal(header_gth[0]);

	// define the maximal offset planes ( ## these hardcoded paramters need to be set from the module interface, which need to be updated later ## ) 
	int offset_plane_max 	  = (int) ((_max_offset - _min_offset)/_step_offset + 1);
	
	// define the used offset-binning plan ( ## these hardcoded paramters need to be set from the module interface, which need to be updated later ##)
	double[] offset_BinCenter = Matlab.linSpace(_min_offset, _max_offset, offset_plane_max);
	
	// define the CRS interpolation used maximal offset-aperture(half aperture).  ( ## these hardcoded paramters need to be set from the module interface, which need to be updated later ## ) 
	// int offset_aperture_max	  = 5;
	
	// define the CRS interpolation used maximal offset-aperture(half aperture).  
        int [] offset_aperture_max    = new int [offset_plane_max]; 
        offset_aperture_max 	      = Matlab.array_double2int(Matlab.interpolateLinear(_apt_offset_val_contp, Matlab.array_double2int(_apt_offset_contp), Matlab.array_double2int(Matlab.linSpace(1, offset_plane_max, offset_plane_max)))); 
        

        //QC the loaded parameters
        if (qc_log_flag == 1){
      		LOG.info("------------------- Basic information of input seismic data ----------------------------------" );
      		LOG.info("++++++++++++++++++++ QC: seismic data dimension (samples * traces): "  + ntimes_gth + " * " + ntr_gth + "++++++++++++++++++++++++");
		LOG.info("++++++++++++++++++++ QC:  _min_offset: "  		+ _min_offset + "++++++++++++++++++++++++");
		LOG.info("++++++++++++++++++++ QC:  _max_offset : " 		+ _max_offset + "++++++++++++++++++++++++");
		LOG.info("++++++++++++++++++++ QC:  _step_offset: " 		+ _step_offset + "++++++++++++++++++++++++");
		LOG.info("++++++++++++++++++++ QC:  offset_plane_max  : " 	+ offset_plane_max + "++++++++++++++++++++");
		for (int i = 0; i< offset_plane_max; i++){
		
			LOG.info("----- offset_aperture_max: QC sample: "  + i + " value: " + offset_aperture_max[i] + "------" );
		
		}
	}

	
	//----------------------------------------------------------------------------------------------------------
      	//   ***** get the reference stacking velocity and CRS parameters (2d array) for current CRS supergather   *****
      	//----------------------------------------------------------------------------------------------------------
      	
      	// get the velocity and parameters arrays for the whole frame
      	float[][] vels2d 	= sidewaysVelData.getTraces();
      	float[][] CRS_para_A1 	= sidewaysCRSData_A1.getTraces();
      	float[][] CRS_para_A2 	= sidewaysCRSData_A2.getTraces();
      	float[][] CRS_para_B11 	= sidewaysCRSData_B11.getTraces();
      	float[][] CRS_para_B12 	= sidewaysCRSData_B12.getTraces();
      	float[][] CRS_para_B22 	= sidewaysCRSData_B22.getTraces();
      	
      	// get the dimension of parameter files (be aware that the data in promax is transposed)
        int ncmps_para 		= vels2d.length;    
        int ntimes_para 	= vels2d[0].length; 
      	
	// define the velocity and parameters arrays for current CMP only
	double[] velRef_array 		= new double[ntimes_gth];
	double[] A1_array		= new double[ntimes_gth];
	double[] A2_array		= new double[ntimes_gth];
	double[] B11_array		= new double[ntimes_gth];
	double[] B12_array		= new double[ntimes_gth];
	double[] B22_array		= new double[ntimes_gth];
	
	//get the index of Xline from the velocity cube
	int[] _xline_array 			= Matlab.linSpaceInt(_xline_Vel[0],_xline_Vel[1]);
	ArrayList<Integer> ind_xline_list 	= Matlab.find(_xline_array,_XL_SG);
	int ind_xline 				= ind_xline_list.get(0);
	
	// get the parameter based on the indecies		
	velRef_array  		= Matlab.getRow(Matlab.subset(vels2d,      ind_xline, ind_xline, 0, ntimes_gth-1),0);
	A1_array  		= Matlab.getRow(Matlab.subset(CRS_para_A1, ind_xline, ind_xline, 0, ntimes_gth-1),0);
	A2_array  		= Matlab.getRow(Matlab.subset(CRS_para_A2, ind_xline, ind_xline, 0, ntimes_gth-1),0);
	B11_array  		= Matlab.getRow(Matlab.subset(CRS_para_B11,ind_xline, ind_xline, 0, ntimes_gth-1),0);
	B12_array  		= Matlab.getRow(Matlab.subset(CRS_para_B12,ind_xline, ind_xline, 0, ntimes_gth-1),0);
	B22_array  		= Matlab.getRow(Matlab.subset(CRS_para_B22,ind_xline, ind_xline, 0, ntimes_gth-1),0);
	
	
	//QC the loaded parameters
	if (qc_log_flag == 1){
		
		LOG.info("vel_array (QC samples: 100,300,500):" + velRef_array[100] + " : " + velRef_array[300] + " : " + velRef_array[500]);
		LOG.info("A1_array (QC samples: 100,300,500):"  + A1_array[100]     + " : " + A1_array[300]     + " : " + A1_array[500]);
		LOG.info("A2_array (QC samples: 100,300,500):"  + A2_array[100]     + " : " + A2_array[300]     + " : " + A2_array[500]);
		LOG.info("B11_array (QC samples: 100,300,500):" + B11_array[100]    + " : " + B11_array[300]    + " : " + B11_array[500]);
		LOG.info("B12_array (QC samples: 100,300,500):" + B12_array[100]    + " : " + B12_array[300]    + " : " + B12_array[500]);
		LOG.info("B22_array (QC samples: 100,300,500):" + B22_array[100]    + " : " + B22_array[300]    + " : " + B22_array[500]);
		LOG.info("B22_array (QC samples: 100,300,500):" + B22_array[100]    + " : " + B22_array[300]    + " : " + B22_array[500]);
 	}

	//----------------------------------------------------------------------------------------------------------
	//   ***** extract the CRS needed trace header from the input prestack data  *****
	//---------------------------------------------------------------------------------------------------------- 
        
        // add the header entry 
    	_QC_Header = data.getHdrCatalog().getEntry(HdrEntry.OFB_NO);
	
	// extract trace headers from the data 
	//int[] hdr_test 			= new int[ntr_gth];
	int[] 	 _Inline		= new int[ntr_gth];
	int[] 	 _Xline			= new int[ntr_gth];
	int[] 	 _trc_type		= new int[ntr_gth];
	double[] _offsets 		= new double[ntr_gth]; 
	int[]    _offsets_bin		= new int[ntr_gth];
	double[] _SourceX		= new double[ntr_gth];
	double[] _SourceY		= new double[ntr_gth];
	double[] _GroupX		= new double[ntr_gth];
	double[] _GroupY		= new double[ntr_gth];
	double[] _MidPointX		= new double[ntr_gth];
	double[] _MidPointY		= new double[ntr_gth];	
	double[] _MidPointX_local	= new double[ntr_gth];
	double[] _MidPointY_local	= new double[ntr_gth];
	double[] _CDPX			= new double[ntr_gth];
	double[] _CDPY			= new double[ntr_gth];
	double[] _CDPX_local		= new double[ntr_gth];
	double[] _CDPY_local		= new double[ntr_gth];
	double[] _CDP2MidPoint_dist	= new double[ntr_gth];   // trace header to save the distance between actual midpoint and bin center

     	for (int i = 0; i < ntr_gth; i++){ 
     	
     		// get the trace-header for every trace
		//hdr_test[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.OFB_NO).getIntVal(header_gth[i]);
		_Inline[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.ILINE_NO).getIntVal(header_gth[i]);
		_Xline[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.XLINE_NO).getIntVal(header_gth[i]);	
		_trc_type[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.TRC_TYPE).getIntVal(header_gth[i]);
		_offsets[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.OFFSET).getDoubleVal(header_gth[i]);
		_offsets_bin[i] 		=  data.getHdrCatalog().getEntry(HdrEntry.DISKITER).getIntVal(header_gth[i]);  // The header DISKITER is used for temporarily store the offset number
		_SourceX[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.SOU_XD).getDoubleVal(header_gth[i]);
		_SourceY[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.SOU_YD).getDoubleVal(header_gth[i]);
		_GroupX[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.REC_XD).getDoubleVal(header_gth[i]);
		_GroupY[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.REC_YD).getDoubleVal(header_gth[i]);
		_MidPointX[i] 			=  0.5 * ( _SourceX[i] + _GroupX[i]);
		_MidPointY[i] 			=  0.5 * ( _SourceY[i] + _GroupY[i]);
		_CDPX[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.CDP_XD).getDoubleVal(header_gth[i]);
		_CDPY[i] 			=  data.getHdrCatalog().getEntry(HdrEntry.CDP_YD).getDoubleVal(header_gth[i]);
		
		// derive the local coordinates by transformation
		double[][] _cdp_local 		=  GeometryTransform(_CDPX[i], _CDPY[i], _ref_point_x, _ref_point_y, _ref_point_inline, _ref_point_xline, _azimuth, _inline_step, _xline_step);
		double[][] _MidPoint_local 	=  GeometryTransform(_MidPointX[i], _MidPointY[i], _ref_point_x, _ref_point_y, _ref_point_inline, _ref_point_xline, _azimuth, _inline_step, _xline_step);
		_CDPX_local[i] 			=  _cdp_local[0][0];
		_CDPY_local[i] 			=  _cdp_local[1][0]; 
		_MidPointX_local[i] 		=  _MidPoint_local[0][0];
		_MidPointY_local[i] 		=  _MidPoint_local[1][0]; 	
		
		//derive the distance between center of bin and mid-point location
		_CDP2MidPoint_dist[i]		=  Math.sqrt(Math.pow(_CDPX_local[i]-_MidPointX_local[i],2) + Math.pow(_CDPY_local[i]-_MidPointY_local[i],2));			
      	}
      	
	
	//----------------------------------------------------------------------------------------------------------
	//  ***** derive the CRS calculation needed parameters  *****
	//----------------------------------------------------------------------------------------------------------
	
	// get the local coordinates of reference-cmp
	int []_Rcmp_trcs_indices 	=  	findIndex_FromArrays(_Inline, _IL_SG, _Xline, _XL_SG );
	int []_trcs_indices_uniq	=  	Matlab.unique(_Rcmp_trcs_indices);
	double _Rcmp_X  		=  	_CDPX_local[_trcs_indices_uniq[0]]; // refernce point is the center of the grid of the reference CMP
	double _Rcmp_Y  		=  	_CDPY_local[_trcs_indices_uniq[0]];
	
	// get the velocity of reference-cmp
	double[] _Rcmp_vel 		= 	velRef_array;

	// get the spatial displament and half-offset
	double[] _delta_X		=	new double[ntr_gth];
	double[] _delta_Y		=	new double[ntr_gth];
	double[] _offset_h		=	new double[ntr_gth];
	
	for (int i = 0; i < ntr_gth; i++) {
		_delta_X[i] 	= _Rcmp_X - _MidPointX_local[i];
		_delta_Y[i] 	= _Rcmp_Y - _MidPointY_local[i];
		_offset_h[i]	= _offsets[i]/2;
	}	 	
	
	// QC the parameters
	if (qc_log_flag == 1){
	
		LOG.info("+++++++++++++ Current SG_inline is: " + _IL_SG + " and SG_xline is: " +  _XL_SG + "+++++++++++++" );
		LOG.info("+++++++++++++ Current reference CMP_X is: " + _Rcmp_X  + " and reference CMP_Y is: " +  _Rcmp_Y + "+++++++++++++" );
		LOG.info("+++++++++++++ Current reference velocity (sample 100) is: " + _Rcmp_vel[100]  + " and sample 500 is: " + _Rcmp_vel[500]);
		LOG.info("+++++++++++++ Current QC coordinate displament: _delta_X[1500] is: "  + _delta_X[1500-1]   + " and _delta_Y[1500] is: "  +  _delta_Y[1500-1]  + "+++++++++++++" );
		LOG.info("+++++++++++++ Current QC coordinate displament: _offset_h[1500] is: " + _offset_h[1500-1]  + " and _offset_h[2000] is: " +  _offset_h[2000] + "+++++++++++++" );
	}
	
	
	//----------------------------------------------------------------------------------------------------------
	//  ***** derive the CRS travel time *****
	//----------------------------------------------------------------------------------------------------------
	
	// define arrays to save the CRS gathers
	double[][] gather_org = new double[ntr_gth][ntimes_gth];  
	double[][] gather_crs = new double[ntr_gth][ntimes_gth];  
	double[][] gather_nmo = new double[ntr_gth][ntimes_gth];  
	
	//define arrays to save the calculated travel time
	double[][] _t_crs_2d  = new double[ntr_gth][ntimes_gth];   
	
	//set the vector t0
	 double[] _t0_array = Matlab.linSpace(0, (ntimes_gth-1)*_time_step, ntimes_gth);
	 
	 //QC the input parameters' dimension
	 if (qc_log_flag == 1){
	 
		 LOG.info("===================== QC input parameters' dimension ==============================================" );
		 LOG.info("+++++++++++++ Current QC time range: is: " + _t0_array [0]  + " : " + _t0_array[ntimes_gth-1] + "+++++++++++++" );
		 LOG.info("+++++++++++++ _t0 array: " 	    + _t0_array.length 	+  "+++++++++++++" );
		 LOG.info("+++++++++++++ _delta_X array: "  + _delta_X.length   +  "+++++++++++++" );
		 LOG.info("+++++++++++++ _delta_Y array: "  + _delta_Y.length   +  "+++++++++++++" );
		 LOG.info("+++++++++++++ A1_array array: "  +  A1_array.length  +  "+++++++++++++" );
		 LOG.info("+++++++++++++ A2_array array: "  +  A2_array.length  +  "+++++++++++++" );
		 LOG.info("+++++++++++++ B11_array array: " +  B11_array.length +  "+++++++++++++" );
		 LOG.info("+++++++++++++ B12_array array: " +  B12_array.length +  "+++++++++++++" );
		 LOG.info("+++++++++++++ B22_array array: " +  B22_array.length +  "+++++++++++++" );
		 LOG.info("+++++++++++++ _offset_h array: " +  _offset_h.length +  "+++++++++++++" );
		 LOG.info("+++++++++++++ velRef    array: " +  velRef_array.length +  "+++++++++++++" );
	 }
	
	//---------------------------------------------------------------------------------------------------------- 
	//  *****  Loop through every trace to calculate CRS travel time and derive moveout corrected traces   ***** 
	//----------------------------------------------------------------------------------------------------------
	
	for (int i = 0; i < ntr_gth; i++) {
	
		// get the single trace
		double[] trace	=	Matlab.getRow(Matlab.matrix_float2double(data2d),i);
		
		// get the half-offset and displacement for the conisdered trace
		double _h 	=	_offset_h[i];	
		double _dx 	=	_delta_X[i];	 
		double _dy 	=	_delta_Y[i];
	
		// calculate 3D ZO CRS travel time 
		double[] _t_crs = crs_ZO_3D_travelTime(_t0_array, A1_array, A2_array, B11_array, B12_array, B22_array, velRef_array, _dx, _dy, _h);
		
		// QC the calculated CRS travel-time
		if (qc_log_flag == 1){ 
		
			// define the QC traces
			int qc_trace  = 1000;
			int qc_sample = 200; 
		
			LOG.info("===================== QC CRS travel time ==============================================" );
			//LOG.info("+++++++++++++ QC at Trace: " + qc_trace + " Sample " + String.valueOf(qc_sample) +  ": " + _t_crs[qc_sample-1]  +  "+++++++++++++" );
			LOG.info("+++++++++++++ QC at Trace: " + qc_trace + " Sample: " + qc_sample +  "+++++++++++++" );
			LOG.info("+++++++++++++ Inline/Xline:" + _Inline[qc_trace -1]   + "/" + _Xline[qc_trace -1] +  "+++++++++++++" );
			LOG.info("+++++++++++++ half-offset: " + _h  + "+++++++++++++" );
			LOG.info("+++++++++++++ velocity: "    + velRef_array[qc_sample-1]  + "+++++++++++++" );
			LOG.info("+++++++++++++ dx: 	     " + _dx + "+++++++++++++" );
			LOG.info("+++++++++++++ dy:          " + _dy + "+++++++++++++" );
			LOG.info("+++++++++++++ A1:          " + A1_array[qc_sample-1]  + "+++++++++++++" );
			LOG.info("+++++++++++++ A2:          " + A2_array[qc_sample-1]  + "+++++++++++++" );
			LOG.info("+++++++++++++ B11:         " + B11_array[qc_sample-1] + "+++++++++++++" );
			LOG.info("+++++++++++++ B12:         " + B12_array[qc_sample-1] + "+++++++++++++" );
			LOG.info("+++++++++++++ B22:         " + B22_array[qc_sample-1] + "+++++++++++++" );
			LOG.info("+++++++++++++ t0:          " + _t0_array[qc_sample-1] + "+++++++++++++" );
			LOG.info("+++++++++++++ tCRS:        " + _t_crs[qc_sample-1]    + "+++++++++++++" );
		}
		
		// derive the CRS movout corrected trace with Linear interpolation	
		gather_crs[i] = Matlab.FastLinearInterp(_t0_array,trace,_t_crs,(double)0);
		
		//save the travel time to travel-time array
		_t_crs_2d[i]  = _t_crs;		
	}
	
	
	
	//---------------------------------------------------------------------------------------------------------- 
	//  *****  Apply CRS based interpolation  ***** 
	//----------------------------------------------------------------------------------------------------------
	
	// define output seismic array
	double[][] traces_org      = new double [offset_plane_max][ntimes_gth]; 
	double[][] traces_intp 	   = new double [offset_plane_max][ntimes_gth];
	double[][] traces_org_rmo  = new double [offset_plane_max][ntimes_gth]; 
	double[][] traces_intp_rmo = new double [offset_plane_max][ntimes_gth];
		
	
	// define the traceheader array
	int[]    _iline_org       = new int [offset_plane_max]; 
	int[]    _xline_org       = new int [offset_plane_max];
	double[] _SourceX_org     = new double [offset_plane_max]; 
	double[] _SourceY_org     = new double [offset_plane_max];
	double[] _GroupX_org      = new double [offset_plane_max];
	double[] _GroupY_org      = new double [offset_plane_max];
	double[] _MidPointX_org   = new double [offset_plane_max];
	double[] _MidPointY_org   = new double [offset_plane_max];
	double[] _CDPX_org  	  = new double [offset_plane_max];
	double[] _CDPY_org        = new double [offset_plane_max];
	double[] _Offset_org      = new double [offset_plane_max];
	int[]    _OffsetBin_org   = new int [offset_plane_max];	
	int[]    _TraceType_org   = new int [offset_plane_max];	
	
	int[]    _iline_intp      = new int [offset_plane_max]; 
	int[]    _xline_intp      = new int [offset_plane_max];
	double[] _SourceX_intp    = new double [offset_plane_max]; 
	double[] _SourceY_intp    = new double [offset_plane_max];
	double[] _GroupX_intp     = new double [offset_plane_max];
	double[] _GroupY_intp     = new double [offset_plane_max];
	double[] _MidPointX_intp  = new double [offset_plane_max];
	double[] _MidPointY_intp  = new double [offset_plane_max];
	double[] _CDPX_intp  	  = new double [offset_plane_max];
	double[] _CDPY_intp       = new double [offset_plane_max];
	double[] _Offset_intp     = new double [offset_plane_max];
	int[]    _OffsetBin_intp  = new int [offset_plane_max];	
	int[]    _TraceType_intp  = new int [offset_plane_max];	
	

	// interpolation
	for (int ioff = 0; ioff < offset_plane_max; ioff++) {
	
		//get actual finite offset number (offset number starting from 1)
		int ioff_plane      	= ioff + 1; 
	
		//derive the original trace(s) at considered finite offset
		int[] trc_flag 	    	= findIndex_FromArrays(_Inline, _IL_SG, _Xline, _XL_SG, _offsets_bin, ioff_plane);
		
		// set the inline/xline and offset number for the output traceheader
		_iline_intp[ioff]   	= _IL_SG;
		_xline_intp[ioff]   	= _XL_SG;
		_OffsetBin_intp[ioff]   = ioff_plane;
		_iline_org[ioff]   	= _IL_SG;
		_xline_org[ioff]   	= _XL_SG;
		_OffsetBin_org[ioff]    = ioff_plane;
		
		// QC the found existed traces of original considered cmp gather
		//LOG.info("+++++++++++++ QC before CRS interpolation----- finite-Offset number: " + ioff + ", --- found the total existed traces: " + trc_flag.length + "+++++++++++++" );
		
		// Case-1: No traces are presented in current offsetplane, then apply interpolation for the missing trace.
		if (trc_flag.length == 0){

			// initiate the variables
			int sum_trcs_flag 	= 0;
			int offset_interp_aprt 	= 1;
			
			// derive the trace indices based on defined CRS aperture
			while (sum_trcs_flag == 0) {
			
				//define the offset range for summation
				int offset_apert_start = ioff_plane  - offset_interp_aprt;
		               	int offset_apert_end   = ioff_plane  + offset_interp_aprt;
		               
		               	//adjust the offset range near the edge
		               	if(offset_apert_start <1){
		                	offset_apert_start = 1;
		                }
		                if(offset_apert_end > offset_plane_max){
		                    	offset_apert_end = offset_plane_max;
		                }
		                	
		                //select traces within the aperture
				ArrayList<Integer> trc_found_flag_arraylist     = Matlab.findBetween(Matlab.array_int2double(_offsets_bin), (double)offset_apert_start, (double) offset_apert_end);
				
				if (trc_found_flag_arraylist.size()==0){
					sum_trcs_flag = 0;
					offset_interp_aprt++;
				}else{
					int[] trc_found_flag 			= Matlab.arraylist2array_int(trc_found_flag_arraylist);
					sum_trcs_flag 				= trc_found_flag.length;
					
					// apply the summation to derive the interpolated trace 
					traces_intp[ioff] 			= Matlab.sum(Matlab.getMultipleRows(gather_crs,trc_found_flag), 0);
					traces_intp[ioff] 			= ArrayMath.div(traces_intp[ioff],(double) trc_found_flag.length);
					
					// output zero trace for current location of the original data
					traces_org[ioff] 			= ArrayMath.mul(traces_intp[ioff],(double)0);
			
					//check if the half-aperture is larger than the pre-defined maximal allowed traces
					if ( offset_interp_aprt > offset_aperture_max[ioff]){
					
						// output zero trace for current location of the original data
						traces_intp[ioff] 			= Matlab.linSpace((double) 0, (double) 0, ntimes_gth);
						traces_org[ioff]  			= ArrayMath.mul(traces_intp[ioff],(double)0);
						break;
					}
				} 	
			}
			

			// --- derive the trace header for current interpolated trace ---

			//get the CDP/CMP coordinates
			int[] trc_flag_tmp		=   findIndex_FromArrays(_Inline, _IL_SG, _Xline, _XL_SG);   //find the trace indices of current CDP location.
			_CDPX_intp[ioff] 		=   Matlab.select(_CDPX,trc_flag_tmp[0],trc_flag_tmp[0])[0]; //get bin_center coordinates from 1st trace (all found traces shall have identical coord)
			_CDPY_intp[ioff] 		=   Matlab.select(_CDPY,trc_flag_tmp[0],trc_flag_tmp[0])[0];
			_MidPointX_intp[ioff] 		=   _CDPX_intp[ioff];
			_MidPointY_intp[ioff] 		=   _CDPY_intp[ioff];			
			
			//get the norminal offset value	
			_Offset_intp[ioff]		=   offset_BinCenter[ioff]; 
			
			//derive the shot/receiver coordinates based on survey azimuth and norminal offset
			double[] shot_receiver_coords 	=   shot_receiver_coord_generate(_CDPX_intp[ioff], _CDPY_intp[ioff], _azimuth, _Offset_intp[ioff]);	
			_SourceX_intp[ioff] 		=   shot_receiver_coords[0];
			_SourceY_intp[ioff] 		=   shot_receiver_coords[1];
			_GroupX_intp[ioff] 		=   shot_receiver_coords[2];
			_GroupY_intp[ioff] 		=   shot_receiver_coords[3];
			
			//set the trace type for the interpolated trace
			_TraceType_intp[ioff]		=   9;
			
			// --- derive the trace header for current original trace ---
			_SourceX_org[ioff] 		=   _SourceX_intp[ioff] ;
			_SourceY_org[ioff] 		=   _SourceY_intp[ioff];
			_GroupX_org[ioff] 		=   _GroupX_intp[ioff] ;
			_GroupY_org[ioff] 		=   _GroupY_intp[ioff];
			_MidPointX_org[ioff] 		=   _MidPointX_intp[ioff];
			_MidPointY_org[ioff] 		=   _MidPointY_intp[ioff];
			_CDPX_org[ioff] 		=   _CDPX_intp[ioff];
			_CDPY_org[ioff] 		=   _CDPY_intp[ioff];
			_Offset_org[ioff]		=   _Offset_intp[ioff];
			_TraceType_org[ioff]		=   0;
			
		}else{ // Case-2: One or multiple traces found in current offsetplane, keep the nearest trace or apply CRS interpolation with small aperture to derive the output.
			
			
			//find the trace which is mostly closed to the bin-center
			ArrayList<Double> _CDP2MidPoint_dist_ArrayList 		=  Matlab.select(_CDP2MidPoint_dist,trc_flag);
			double[] _CDP2MidPoint_dist_traces			=  Matlab.arraylist2array_double(_CDP2MidPoint_dist_ArrayList);
			int trc_nearest_flag					=  trc_flag[Matlab.min_index(_CDP2MidPoint_dist_traces)]; 
			
			// get the seismic trace in current location for original data
			traces_org[ioff]					=   Matlab.getRow(gather_crs, trc_nearest_flag);
			
			// treatment of existing trace, option(1): keep the original trace.
			if (_output_mode == 1 || _output_mode == 3 ) { 
				
				// get the seismic trace in current location
				traces_intp[ioff]				=   traces_org[ioff];

				//get the trace header for current trace
				_SourceX_org[ioff] 				=   Matlab.select(_SourceX,trc_nearest_flag,trc_nearest_flag)[0];
				_SourceY_org[ioff] 				=   Matlab.select(_SourceY,trc_nearest_flag,trc_nearest_flag)[0];
				_GroupX_org[ioff] 				=   Matlab.select(_GroupX,trc_nearest_flag,trc_nearest_flag)[0];
				_GroupY_org[ioff] 				=   Matlab.select(_GroupY,trc_nearest_flag,trc_nearest_flag)[0];
				_MidPointX_org[ioff] 				=   Matlab.select(_MidPointX,trc_nearest_flag,trc_nearest_flag)[0];
				_MidPointY_org[ioff] 				=   Matlab.select(_MidPointY,trc_nearest_flag,trc_nearest_flag)[0];
				_CDPX_org[ioff] 				=   Matlab.select(_CDPX,trc_nearest_flag,trc_nearest_flag)[0];
				_CDPY_org[ioff] 				=   Matlab.select(_CDPY,trc_nearest_flag,trc_nearest_flag)[0];
				_Offset_org[ioff]				=   Matlab.select(_offsets,trc_nearest_flag,trc_nearest_flag)[0];
				_TraceType_org[ioff]				=   Matlab.select(_trc_type,trc_nearest_flag,trc_nearest_flag)[0];
				
				// --- derive the trace header for interpolated trace ---
				_SourceX_intp[ioff] 				= _SourceX_org[ioff];
				_SourceY_intp[ioff] 				= _SourceY_org[ioff];
				_GroupX_intp[ioff]  				= _GroupX_org[ioff];
				_GroupY_intp[ioff]  				= _GroupY_org[ioff];
				_MidPointX_intp[ioff]  				= _MidPointX_org[ioff];
				_MidPointY_intp[ioff]  				= _MidPointY_org[ioff];
				_CDPX_intp[ioff]       				= _CDPX_org[ioff];
				_CDPY_intp[ioff]       				= _CDPY_org[ioff];
				_Offset_intp[ioff]     				= _Offset_org[ioff];
				_TraceType_intp[ioff]  				= _TraceType_org[ioff];
				

			} else{ // treatment of existing trace,  option(2): apply CRS interpolation to derive the output new trace
			
				//define the offset range for summation. (only use neaby offset planes for interpolation)
				int offset_apert_start = ioff_plane  - offset_aperture_max[ioff];
                        	int offset_apert_end   = ioff_plane  + offset_aperture_max[ioff];
                        	
                        	//adjust the offset range near the edge
                        	if(offset_apert_start <1){
                            		offset_apert_start = 1;
                        	}
                        	if(offset_apert_end > offset_plane_max){
                            		offset_apert_end = offset_plane_max;
                        	}
                        	
                                //select traces within the aperture
				ArrayList<Integer> trc_found_flag_arraylist     =   Matlab.findBetween(Matlab.array_int2double(_offsets_bin), (double)offset_apert_start, (double) offset_apert_end);
				int[] trc_found_flag 			        =   Matlab.arraylist2array_int(trc_found_flag_arraylist);

				// apply the arithmatic averege to derive the interpolated trace 
				double[][] crs_selected_traces_tmp 		=   Matlab.getMultipleRows(gather_crs,trc_found_flag);
				traces_intp[ioff] 				=   ArrayMath.div(Matlab.sum(crs_selected_traces_tmp, 0),(double) trc_found_flag.length);
				
				//--- get the trace header for current trace based on survey azimuth and norminal offset ---
				//get the cdp and midpoint coordinates
				_CDPX_intp[ioff] 				=   Matlab.select(_CDPX,trc_nearest_flag,trc_nearest_flag)[0];
				_CDPY_intp[ioff] 				=   Matlab.select(_CDPY,trc_nearest_flag,trc_nearest_flag)[0];
				_MidPointX_intp[ioff] 				=   _CDPX_intp[ioff];
				_MidPointY_intp[ioff] 				=   _CDPY_intp[ioff];
				
				//get the norminal offset value	DISKITER
				_Offset_intp[ioff]				=   offset_BinCenter[ioff]; 
				
				//derive the shot/receiver coordinates based on survey azimuth and norminal offset
				double[] shot_receiver_coords 			=   shot_receiver_coord_generate(_CDPX_intp[ioff], _CDPY_intp[ioff], _azimuth, _Offset_intp[ioff]);	
				_SourceX_intp[ioff] 				=   shot_receiver_coords[0];
				_SourceY_intp[ioff] 				=   shot_receiver_coords[1];
				_GroupX_intp[ioff] 				=   shot_receiver_coords[2];
				_GroupY_intp[ioff] 				=   shot_receiver_coords[3];
				
				//set the trace type for the existing trace
				_TraceType_intp[ioff]				=   9;
				
				// --- derive the trace header for current original trace ---
				_SourceX_org[ioff] 				=   Matlab.select(_SourceX,trc_nearest_flag,trc_nearest_flag)[0];
				_SourceY_org[ioff] 				=   Matlab.select(_SourceY,trc_nearest_flag,trc_nearest_flag)[0];
				_GroupX_org[ioff] 				=   Matlab.select(_GroupX,trc_nearest_flag,trc_nearest_flag)[0];
				_GroupY_org[ioff] 				=   Matlab.select(_GroupY,trc_nearest_flag,trc_nearest_flag)[0];
				_MidPointX_org[ioff] 				=   Matlab.select(_MidPointX,trc_nearest_flag,trc_nearest_flag)[0];
				_MidPointY_org[ioff] 				=   Matlab.select(_MidPointY,trc_nearest_flag,trc_nearest_flag)[0];
				_CDPX_org[ioff] 				=   Matlab.select(_CDPX,trc_nearest_flag,trc_nearest_flag)[0];
				_CDPY_org[ioff] 				=   Matlab.select(_CDPY,trc_nearest_flag,trc_nearest_flag)[0];
				_Offset_org[ioff]				=   Matlab.select(_offsets,trc_nearest_flag,trc_nearest_flag)[0];
				_TraceType_org[ioff]				=   1;
				
			}
		}
		
		//--- derive the moveout-reversed CMP gathers for the output ---
		
		// derive the reference travel time for current finite-offset plane. (N.B. CRS moveout and NMO moveout are identical at the reference cmp location)
		double[] _t_ref_org  = nmo_travelTime(_t0_array, velRef_array, _Offset_org[ioff]/2.0);	
		double[] _t_ref_intp = nmo_travelTime(_t0_array, velRef_array, _Offset_intp[ioff]/2.0);	
		
		// derive the moveout-reversed CMP gathers.
		traces_org_rmo[ioff]  = Matlab.FastLinearInterp(_t_ref_org,traces_org[ioff], _t0_array,(double)0); 
		traces_intp_rmo[ioff] = Matlab.FastLinearInterp(_t_ref_intp,traces_intp[ioff],_t0_array,(double)0); 

		//qc the trace header of derived CMP offset-plane
		if (qc_log_flag == 1){ 
			LOG.info("------------- QC of offset-plane: "   + ioff_plane + "-----------------------------------------------------" );
			LOG.info("+++++++++++++ SourceX    : " + _SourceX_intp[ioff]    + "+++++++++++++" );
			LOG.info("+++++++++++++ SourceY    : " + _SourceY_intp[ioff]    + "+++++++++++++" );
			LOG.info("+++++++++++++ GroupX     : " + _GroupX_intp[ioff]     + "+++++++++++++" );
			LOG.info("+++++++++++++ GroupY     : " + _GroupY_intp[ioff]     + "+++++++++++++" );
			LOG.info("+++++++++++++ MidPointX  : " + _MidPointX_intp[ioff]  + "+++++++++++++" );
			LOG.info("+++++++++++++ MidPointY  : " + _MidPointY_intp[ioff]  + "+++++++++++++" );
			LOG.info("+++++++++++++ CDPX       : " + _CDPX_intp[ioff]       + "+++++++++++++" );
			LOG.info("+++++++++++++ CDPY       : " + _CDPY_intp[ioff]       + "+++++++++++++" );
			LOG.info("+++++++++++++ Offset     : " + _Offset_intp[ioff]     + "+++++++++++++" );
		}
	}
	
	
	
	//----------------------------------------------------------------------------------------------------------
	//   ***** generate the output of cmp with original traces + interpolated missing traces  *****
	//----------------------------------------------------------------------------------------------------------
	

	//derive the original trace(s) at considered reference cmp location
	int[] trc_exist_flag	= findIndex_FromArrays(_Inline, _IL_SG, _Xline, _XL_SG);
	int[] trc_intp_flag	= findIndex_FromArrays(_TraceType_intp, 9);
	int trc_output_length	= trc_exist_flag.length + trc_intp_flag.length;
	
	
	//qc the trace header of derived CMP offset-plane
	if (qc_log_flag == 1){ 
		LOG.info("================  Statistics of applied CRS interpolation ================ " );
		LOG.info("+++++++++++++ Number of original exisiting offset traces:   : " + trc_exist_flag.length   + "+++++++++++++" );
		LOG.info("+++++++++++++ Number of CRS interpolated offset traces:     : " + trc_intp_flag.length    + "+++++++++++++" );
	}	
	
	//derive the original trace(s) at considered reference cmp location
	double[][] traces_org_plus_interp  = new double [trc_output_length][ntimes_gth]; 
	
	// derive the corresponding trace header array
	int[]    _iline_org_pI       = new int [trc_output_length]; 
	int[]    _xline_org_pI       = new int [trc_output_length];
	double[] _SourceX_org_pI     = new double [trc_output_length]; 
	double[] _SourceY_org_pI     = new double [trc_output_length];
	double[] _GroupX_org_pI      = new double [trc_output_length];
	double[] _GroupY_org_pI      = new double [trc_output_length];
	double[] _MidPointX_org_pI   = new double [trc_output_length];
	double[] _MidPointY_org_pI   = new double [trc_output_length];
	double[] _CDPX_org_pI  	     = new double [trc_output_length];
	double[] _CDPY_org_pI        = new double [trc_output_length];
	double[] _Offset_org_pI      = new double [trc_output_length];
	int[]    _OffsetBin_org_pI   = new int [trc_output_length];	
	int[]    _TraceType_org_pI   = new int [trc_output_length];	
	
	
	if (_output_mode == 3){ 


		//generate the output and trace header
		for (int i=0;i<trc_output_length;i++) {
	
	
			if(i<trc_exist_flag.length){
			
				// get the original existing traces
				traces_org_plus_interp[i] =  Matlab.getRow(gather_crs, trc_exist_flag[i]);
			
				// get corresponding trace header
				_iline_org_pI[i]	= _IL_SG;
				_xline_org_pI[i]       	= _XL_SG;
				_SourceX_org_pI[i]   	= _SourceX[trc_exist_flag[i]]; 
				_SourceY_org_pI[i]   	= _SourceY[trc_exist_flag[i]];
				_GroupX_org_pI[i]      	= _GroupX[trc_exist_flag[i]];
				_GroupY_org_pI[i]      	= _GroupY[trc_exist_flag[i]];
				_MidPointX_org_pI[i]   	= _MidPointX[trc_exist_flag[i]];
				_MidPointY_org_pI[i]   	= _MidPointY[trc_exist_flag[i]];
				_CDPX_org_pI[i]  	= _CDPX[trc_exist_flag[i]];
				_CDPY_org_pI[i]        	= _CDPY[trc_exist_flag[i]];
			
				_Offset_org_pI[i]      	= _offsets[trc_exist_flag[i]];
				_OffsetBin_org_pI[i]   	= _offsets_bin[trc_exist_flag[i]];	
				_TraceType_org_pI[i]   	= _trc_type[trc_exist_flag[i]];


			}else {
				// get the crs interpolated missing traces
				traces_org_plus_interp[i] =  Matlab.getRow(traces_intp, trc_intp_flag[i-trc_exist_flag.length]);
			
				// get corresponding trace header
				_iline_org_pI[i]	= _IL_SG;
				_xline_org_pI[i]       	= _XL_SG;
				_SourceX_org_pI[i]   	= _SourceX_intp[trc_intp_flag[i-trc_exist_flag.length]]; 
				_SourceY_org_pI[i]   	= _SourceY_intp[trc_intp_flag[i-trc_exist_flag.length]];
				_GroupX_org_pI[i]      	= _GroupX_intp[trc_intp_flag[i-trc_exist_flag.length]];
				_GroupY_org_pI[i]      	= _GroupY_intp[trc_intp_flag[i-trc_exist_flag.length]];
			
				_MidPointX_org_pI[i]   	= _MidPointX_intp[trc_intp_flag[i-trc_exist_flag.length]];
				_MidPointY_org_pI[i]   	= _MidPointY_intp[trc_intp_flag[i-trc_exist_flag.length]];
				_CDPX_org_pI[i]  	= _CDPX_intp[trc_intp_flag[i-trc_exist_flag.length]];
				_CDPY_org_pI[i]        	= _CDPY_intp[trc_intp_flag[i-trc_exist_flag.length]];
			
				_Offset_org_pI[i]      	= _Offset_intp[trc_intp_flag[i-trc_exist_flag.length]];
				_OffsetBin_org_pI[i]   	= _OffsetBin_intp[trc_intp_flag[i-trc_exist_flag.length]];	
				_TraceType_org_pI[i]   	= _TraceType_intp[trc_intp_flag[i-trc_exist_flag.length]];
		
			}
	
		}
	}
	//----------------------------------------------------------------------------------------------------------
	//   ***** set the trace header for the interpolated data  *****
	//----------------------------------------------------------------------------------------------------------

	if ( _output_mode == 1 || _output_mode == 2 ) {
	
		// define output planes
		int output_flag = 0; 
	
	     	for (int ind = 0; ind < ntr_gth; ind++){ 
	     	
	     		// set the internal index
	     		int i = ind;
	     	
			// find the index 
			if( ind < offset_plane_max){
				i = ind;				
				output_flag = 1;
				
			} else if ( ind >= offset_plane_max && ind < 2*offset_plane_max){
				
				i = ind    -  1*offset_plane_max;
				output_flag = 2;

			} else if ( ind >=    2* offset_plane_max && ind < 3*offset_plane_max){
				
				i = ind    -  2*offset_plane_max;
				output_flag = 3;
					
			} else if ( ind >=    3* offset_plane_max && ind < 4*offset_plane_max){
				
				i = ind    -  3*offset_plane_max;
				output_flag = 4;
			
			} else if ( ind >=    4* offset_plane_max && ind < 5*offset_plane_max){
				
				i = ind    -  4*offset_plane_max;
				output_flag = 5;

			}else{
				i = ind;
			}
	     		
	     		if (ind < 5*offset_plane_max){
	     		
	     			if ( output_flag == 2 || output_flag == 4 ){

					data.getHdrCatalog().getEntry(HdrEntry.ILINE_NO).setIntVal(header_gth[ind],_iline_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.XLINE_NO).setIntVal(header_gth[ind],_xline_intp[i]);	
					data.getHdrCatalog().getEntry(HdrEntry.OFFSET).setDoubleVal(header_gth[ind],_Offset_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.AOFFSET).setDoubleVal(header_gth[ind],Math.abs(_Offset_intp[i]));
					data.getHdrCatalog().getEntry(HdrEntry.DISKITER).setIntVal(header_gth[ind],_OffsetBin_intp[i]);
					//data.getHdrCatalog().getEntry(HdrEntry.OFB_NO).setIntVal(header_gth[ind],_OffsetBin_intp[i]);   
					//should not redefine this header, otherwise it will change the original frame.
			
					data.getHdrCatalog().getEntry(HdrEntry.SOU_XD).setDoubleVal(header_gth[ind],_SourceX_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.SOU_YD).setDoubleVal(header_gth[ind],_SourceY_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_XD).setDoubleVal(header_gth[ind],_GroupX_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_YD).setDoubleVal(header_gth[ind],_GroupY_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_XD).setDoubleVal(header_gth[ind],_CDPX_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_YD).setDoubleVal(header_gth[ind],_CDPY_intp[i]);
			
					data.getHdrCatalog().getEntry(HdrEntry.SOU_X).setFloatVal(header_gth[ind],(float) _SourceX_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.SOU_Y).setFloatVal(header_gth[ind],(float) _SourceY_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_X).setFloatVal(header_gth[ind],(float) _GroupX_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_Y).setFloatVal(header_gth[ind],(float) _GroupY_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_X).setFloatVal(header_gth[ind],(float) _CDPX_intp[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_Y).setFloatVal(header_gth[ind],(float) _CDPY_intp[i]);
					
					data.getHdrCatalog().getEntry(HdrEntry.CHAN).setFloatVal(header_gth[ind],ind+1);
					data.getHdrCatalog().getEntry(HdrEntry.TRC_TYPE).setIntVal(header_gth[ind],_TraceType_intp[i]);
					if (_TraceType_intp[i] == 9){
						_catalogedHdrEntry_1.setIntVal(header_gth[ind], 1);
					}else{
						_catalogedHdrEntry_1.setIntVal(header_gth[ind], 0);
					}
					_catalogedHdrEntry_2.setIntVal(header_gth[ind], output_flag);
				}
				
				if ( output_flag == 1 || output_flag == 3 || output_flag == 5){
		     			
					data.getHdrCatalog().getEntry(HdrEntry.ILINE_NO).setIntVal(header_gth[ind],_iline_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.XLINE_NO).setIntVal(header_gth[ind],_xline_org[i]);	
					data.getHdrCatalog().getEntry(HdrEntry.OFFSET).setDoubleVal(header_gth[ind],_Offset_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.AOFFSET).setDoubleVal(header_gth[ind],Math.abs(_Offset_org[i]));
					data.getHdrCatalog().getEntry(HdrEntry.DISKITER).setIntVal(header_gth[ind],_OffsetBin_org[i]);
					//data.getHdrCatalog().getEntry(HdrEntry.OFB_NO).setIntVal(header_gth[ind],_OffsetBin_intp[i]);   
					//should not redefine this header, otherwise it will change the original frame.
			
					data.getHdrCatalog().getEntry(HdrEntry.SOU_XD).setDoubleVal(header_gth[ind],_SourceX_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.SOU_YD).setDoubleVal(header_gth[ind],_SourceY_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_XD).setDoubleVal(header_gth[ind],_GroupX_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_YD).setDoubleVal(header_gth[ind],_GroupY_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_XD).setDoubleVal(header_gth[ind],_CDPX_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_YD).setDoubleVal(header_gth[ind],_CDPY_org[i]);
			
					data.getHdrCatalog().getEntry(HdrEntry.SOU_X).setFloatVal(header_gth[ind],(float) _SourceX_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.SOU_Y).setFloatVal(header_gth[ind],(float) _SourceY_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_X).setFloatVal(header_gth[ind],(float) _GroupX_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.REC_Y).setFloatVal(header_gth[ind],(float) _GroupY_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_X).setFloatVal(header_gth[ind],(float) _CDPX_org[i]);
					data.getHdrCatalog().getEntry(HdrEntry.CDP_Y).setFloatVal(header_gth[ind],(float) _CDPY_org[i]);
					
					data.getHdrCatalog().getEntry(HdrEntry.CHAN).setFloatVal(header_gth[ind],ind+1);
					data.getHdrCatalog().getEntry(HdrEntry.TRC_TYPE).setIntVal(header_gth[ind],_TraceType_org[i]);
					_catalogedHdrEntry_1.setIntVal(header_gth[ind], 0);
					_catalogedHdrEntry_2.setIntVal(header_gth[ind], output_flag);
					
					
				}

			}else {
				_catalogedHdrEntry_2.setIntVal(header_gth[ind], 0);
			}

	      	}
	      	
      	}
      	
      	
      	if ( _output_mode == 3 ) {
      	
      		// define output planes
		int output_flag = 6; 
										
		for (int ind = 0; ind < ntr_gth; ind++){ 
	     	
	     		// set the internal index
	     		int i = ind;
		
	      		if (ind < trc_output_length){
			     			
				data.getHdrCatalog().getEntry(HdrEntry.ILINE_NO).setIntVal(header_gth[ind],_iline_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.XLINE_NO).setIntVal(header_gth[ind],_xline_org_pI[i]);	
				data.getHdrCatalog().getEntry(HdrEntry.OFFSET).setDoubleVal(header_gth[ind],_Offset_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.AOFFSET).setDoubleVal(header_gth[ind],Math.abs(_Offset_org_pI[i]));
				data.getHdrCatalog().getEntry(HdrEntry.DISKITER).setIntVal(header_gth[ind],_OffsetBin_org_pI[i]);
			
				data.getHdrCatalog().getEntry(HdrEntry.SOU_XD).setDoubleVal(header_gth[ind],_SourceX_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.SOU_YD).setDoubleVal(header_gth[ind],_SourceY_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.REC_XD).setDoubleVal(header_gth[ind],_GroupX_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.REC_YD).setDoubleVal(header_gth[ind],_GroupY_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.CDP_XD).setDoubleVal(header_gth[ind],_CDPX_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.CDP_YD).setDoubleVal(header_gth[ind],_CDPY_org_pI[i]);
			
				data.getHdrCatalog().getEntry(HdrEntry.SOU_X).setFloatVal(header_gth[ind],(float) _SourceX_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.SOU_Y).setFloatVal(header_gth[ind],(float) _SourceY_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.REC_X).setFloatVal(header_gth[ind],(float) _GroupX_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.REC_Y).setFloatVal(header_gth[ind],(float) _GroupY_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.CDP_X).setFloatVal(header_gth[ind],(float) _CDPX_org_pI[i]);
				data.getHdrCatalog().getEntry(HdrEntry.CDP_Y).setFloatVal(header_gth[ind],(float) _CDPY_org_pI[i]);
				
				data.getHdrCatalog().getEntry(HdrEntry.CHAN).setFloatVal(header_gth[ind],ind+1);
				data.getHdrCatalog().getEntry(HdrEntry.TRC_TYPE).setIntVal(header_gth[ind],_TraceType_org_pI[i]);
				
				if (_TraceType_org_pI[i] == 9){
					_catalogedHdrEntry_1.setIntVal(header_gth[ind], 1);
				}else{
					_catalogedHdrEntry_1.setIntVal(header_gth[ind], 0);
				}
				_catalogedHdrEntry_2.setIntVal(header_gth[ind], output_flag);				
				
							

	      		}else{
	      		
				_catalogedHdrEntry_2.setIntVal(header_gth[ind], 0);
	      		}
	      }
      	
      	
      	}
      	
	//----------------------------------------------------------------------------------------------------------	
	// ***** QC DISPLAY THE OTHER CHECKED VARIABLES IN LOG *****
	//----------------------------------------------------------------------------------------------------------
	
	if (qc_log_flag == 1){ 
	
	      	LOG.info("================ START DISPLAY THE CHECKED VARIABLES IN LOG FILE ================================");
	      	LOG.info("Found the dimension of input CRS sup-gather ensemble (traces x times):" + ntr_gth + "*" + ntimes_gth + " in the job");
	      	LOG.info("Found the dimension of input parameter files (traces x times):" + ncmps_para + "*" + ntimes_para + " in the job");
	      	LOG.info("Found the dimension of trace header of prestack data (traces x times):" + header_gth.length + "*" + header_gth[0].length + " in the job");
	      	//LOG.info("Found the querying trace header:" +  _QC_Header + " with value range:" + _QC_Header.getFloatVal(header_gth[0]) + " : " + _QC_Header.getFloatVal(header_gth[ntr_gth-1])  + " in the job ");
	      	//LOG.info("Found the querying trace header:" +  _QC_Header + " with value range:" + hdr_test[0] + " : " + hdr_test[ntr_gth-1]  + " in the job ");
	      	LOG.info("Found the querying trace header OFFSET:  with value range:"     + _offsets[0] + " : " + _offsets[ntr_gth-1]  + " in the job ");
	      	LOG.info("Found the querying trace header OFFSET_BIN:  with value range:" + _offsets_bin[0] + " : " + _offsets_bin[ntr_gth-1]  + " in the job ");
	      	LOG.info("Found the querying trace header Source-X:  with value range:"   + _SourceX[0] + " : " + _SourceX[ntr_gth-1]  + " in the job ");
	      	LOG.info("Found the querying trace header Source-Y:  with value range:"   + _SourceY[0] + " : " + _SourceY[ntr_gth-1]  + " in the job ");
	      	LOG.info("Found the querying trace header Receiver-X:  with value range:" + _GroupX[0]  + " : " + _GroupX[ntr_gth-1]   + " in the job ");
	      	LOG.info("Found the querying trace header Receiver-Y:  with value range:" + _GroupY[0]  + " : " + _GroupY[ntr_gth-1]   + " in the job ");
	      	LOG.info("Found the querying trace header CDP-X:  with value range:"      + _CDPX[0]    + " : " + _CDPX[ntr_gth-1]   + " in the job ");
	      	LOG.info("Found the querying trace header CDP-Y:  with value range:"      + _CDPY[0]    + " : " + _CDPY[ntr_gth-1]   + " in the job ");
	      	LOG.info("Found the querying trace header inline:  with value range:"     + _Inline[0]  + " : " + _Inline[ntr_gth-1]   + " in the job ");
	      	LOG.info("Found the querying trace header xline:  with value range:"      + _Xline[0]   + " : " + _Xline[ntr_gth-1]   + " in the job ");
		
		LOG.info("================ QC Seismic Trace Header ================================");
		int qc_trace_n = 1500-1; int qc_sample_n = 400;
		LOG.info("Found the querying trace header Trcace index Number:  with value range:" + qc_trace_n			);
		LOG.info("QC the trace header nth:" + qc_trace_n + " from the total trace header: " + 	_MidPointX.length	);
		
		LOG.info("Found the querying trace header inline:      with value :"      + _Inline[qc_trace_n]  	   	);
	      	LOG.info("Found the querying trace header xline:       with value :"      + _Xline[qc_trace_n]   	   	);
	      	LOG.info("Found the querying trace header OFFSET:      with value :"      + _offsets[qc_trace_n] 	   	);
	      	LOG.info("Found the querying trace header Source-X:    with value :"      + _SourceX[qc_trace_n] 	   	);
	      	LOG.info("Found the querying trace header Source-Y:    with value :"      + _SourceY[qc_trace_n] 	   	);
	      	LOG.info("Found the querying trace header Receiver-X:  with value :"      + _GroupX[qc_trace_n]      		);
	      	LOG.info("Found the querying trace header Receiver-Y:  with value :"      + _GroupY[qc_trace_n]  	   	);
	      	LOG.info("Found the querying trace header CMP-X:       with value :"      + _MidPointX[qc_trace_n]   		);
	      	LOG.info("Found the querying trace header CMP-Y:       with value :"      + _MidPointY[qc_trace_n]   		);      	
	      	LOG.info("Found the querying trace header CDP-X:       with value :"      + _CDPX[qc_trace_n]    	   	);
	      	LOG.info("Found the querying trace header CDP-Y:       with value :"      + _CDPY[qc_trace_n]    	   	);
	      	LOG.info("Found the querying trace header CMP-X_local: with value :"      + _MidPointX_local[qc_trace_n]   	);
	      	LOG.info("Found the querying trace header CMP-Y_local: with value :"      + _MidPointY_local[qc_trace_n]   	);      	
	      	LOG.info("Found the querying trace header CDP-X_local: with value :"      + _CDPX_local[qc_trace_n]   	 	);
	      	LOG.info("Found the querying trace header CDP-Y_local: with value :"      + _CDPY_local[qc_trace_n]    		);
	      	
		LOG.info("================ QC Velocity/parameter volume range ================================");
	      	LOG.info("The velocity volume's information: dimension= " + _vel_fw.getDataDimensions()); 
	      	LOG.info("The velocity volume's information: min-frame= " + _vel_fw.getMinFrame()); 
	      	LOG.info("The velocity volume's information: max-frame= " + _vel_fw.getMaxFrame()); 
	      	LOG.info("The velocity volume's information: min-trace= " + _vel_fw.getMinTrace()); 
	      	LOG.info("The velocity volume's information: max-trace= " + _vel_fw.getMaxTrace()); 
	      	LOG.info("The velocity volume's information: min-sample=" + _vel_fw.getMinSample()); 
	      	LOG.info("The velocity volume's information: max-sample=" + _vel_fw.getMaxSample()); 
	      	
	      	//LOG.info("The velocity volume QC information: Inline/Xline/sample: " + _Inline[qc_trace_n] + "/" + _Xline[qc_trace_n] + "/" + qc_sample_n + " = " + vel_array[qc_trace_n][qc_sample_n]); 
	      	LOG.info("================ COMPLETE DISPLAY THE CHECKED VARIABLES IN LOG FILE =============================");
	
	}
	
	
     //----------------------------------------------------------------------------------------------------------
     // Output by replacing the data to calculated semblance
     //----------------------------------------------------------------------------------------------------------     
     for (int iR = 0; iR<ntr_gth; iR++){ 

		for (int iC = 0; iC<ntimes_gth; iC++){ 
				
                // choose the output (1: the regularized + interpolated CMP gather with/without Moveout Correction 2: the interpolated CMP gather with/without Moveout Correction,and 3: the moveout corrected CRS gather.) 

				if(_output_mode == 4){
					
					data2d[iR][iC] = (float)gather_crs[iR][iC];
				
				}else if(_output_mode == 3){
					
					if( iR < trc_output_length){
						data2d[iR][iC] = (float)traces_org_plus_interp[iR][iC];
					}else{
						data2d[iR][iC] = 0;
					}
				
				}else{

					if( iR < offset_plane_max){
					
						data2d[iR][iC] = (float)traces_org[iR][iC];
				
					} else if ( iR >= 1* offset_plane_max && iR < 2*offset_plane_max){
				
						data2d[iR][iC] = (float)traces_intp[iR-1*offset_plane_max][iC];

					} else if ( iR >= 2* offset_plane_max && iR < 3*offset_plane_max){
				
						data2d[iR][iC] = (float)traces_org_rmo[iR-2*offset_plane_max][iC];
					
					} else if ( iR >= 3* offset_plane_max && iR < 4*offset_plane_max){
				
						data2d[iR][iC] = (float)traces_intp_rmo[iR-3*offset_plane_max][iC];
						
					} else if ( iR >= 4* offset_plane_max && iR < 5*offset_plane_max){
				
						data2d[iR][iC] = (float)_t_crs_2d[iR-4*offset_plane_max][iC];
					}
					
					else{
						data2d[iR][iC] = 0;
					}
				}
		}
			
     }	 

	// free the sideway data from memory
      	sidewaysVelData.free();
      	sidewaysCRSData_A1.free();
      	sidewaysCRSData_A2.free();
      	sidewaysCRSData_B11.free();
      	sidewaysCRSData_B12.free();
      	sidewaysCRSData_B22.free();
   
    }

    return data;       
  }
  
  
  /**
	(1) Method for coordination transformation of acquisition coordinates to local (processing) coordinates
  */

  private double[][] GeometryTransform(double x, double y, double ref_point_x0, double ref_point_y0, int ref_point_inline, int ref_point_xline, double rotation_angle, double inline_step, double xline_step){
  
  
        // convert the degree to radiant angle
        double theta = rotation_angle * Math.PI / 180.0 ;
        
        // get the input coordinates array
        double[][] coord_input 	= {{x - ref_point_x0 , 0},{y - ref_point_y0, 0}};
        
        // derive the transformation matrix
	double[][] Trans_mat 	= { {Math.cos(theta), Math.sin(theta)}, {-Math.sin(theta), Math.cos(theta)}};
	
	//transorm the original acquisition coordination system to processing coordinate system
	DMatrix Trans_mat_DM 	= new DMatrix(Trans_mat);  
	DMatrix coord_input_DM 	= new DMatrix(coord_input );
	DMatrix coord_output_DM	= Trans_mat_DM.times(coord_input_DM);
	double[][] coord_output = coord_output_DM.get();
	  
	//derive the inline and xline number in the processing coordinate system (Only use for the QC purpose)
	int inline_qc   = (int)Math.round(coord_output[1][0]/inline_step) + ref_point_inline;
	int xline_qc	= (int)Math.round(coord_output[0][0]/xline_step)  + ref_point_xline;
	
	
	// qc the result if needed, default without QC.
	int qc_flag = 0; 
	if (qc_flag == 1){	
	
		LOG.info("Debugging Coord_input :"       + "{"  + coord_input[0][0]  + " and " + coord_input[0][1]    + "}" + "{"  + coord_input[1][0]   + " and " + coord_input[1][1]  + "}");
		LOG.info("Debugging Trans_max   :"       + "{"  + Trans_mat[0][0]    + " and " + Trans_mat[0][1]      + "}" + "{"  + Trans_mat[1][0]     + " and " + Trans_mat[1][1]    + "}");
		LOG.info("Debugging coord_output:"       + "{"  + coord_output[0][0] + " and " + coord_output[0][1]   + "}" + "{"  + coord_output[1][0]  + " and " + coord_output[1][1] + "}");
		LOG.info("Debugging inline_qc :"         +  inline_qc);
		LOG.info("Debugging xline_qc :"          +  xline_qc);
	
		LOG.info("Debugging Coord_input dimension:"      + coord_input.length + " X " + coord_input[0].length);
		LOG.info("Debugging Trans_max dimension:"        + Trans_mat.length +   " X " + Trans_mat[0].length);
	}
	
	//output the results
	return coord_output;
  }
  
  /**
	(2) Method for selecting with multiple conditions 
  */

  private int [] findIndex_FromArrays(int[] aData, int aValue, int[] bData, int bValue ) {

	ArrayList<Integer> IndFound = new ArrayList<Integer>();
	
	
	int i;
	int len_a = aData.length;
	int len_b = bData.length;
	
	if (len_a != len_b){
		LOG.info("The input arrays are not in the same length !"); 
	}
	
	for (i=0; i<len_a; i++) {
		if (aData[i]== aValue && bData[i]== bValue ) { 
			IndFound.add(i);
		}
	}
	
	// output the index array
	int [] IndFoundArray = new int [IndFound.size()];
	for (i=0; i<IndFound.size(); i++) {
	IndFoundArray[i] = IndFound.get(i);
	}
	   

	return IndFoundArray;	   
	   
  }
  
    private int [] findIndex_FromArrays(int[] aData, int aValue, int[] bData, int bValue, int[] cData, int cValue ) {

	ArrayList<Integer> IndFound = new ArrayList<Integer>();
	
	
	int i;
	int len_a = aData.length;
	int len_b = bData.length;
	int len_c = cData.length;
	
	if (len_a != len_b || len_a != len_c ){
		LOG.info("The input arrays are not in the same length !"); 
	}
	   
	for (i=0; i<len_a; i++) {
		if (aData[i]== aValue && bData[i]== bValue && cData[i]== cValue) { 
			IndFound.add(i);
		}
	}
	
	// output the index array
	int [] IndFoundArray = new int [IndFound.size()];
	
	for (i=0; i<IndFound.size(); i++) {
		IndFoundArray[i] = IndFound.get(i);
	}
	   

	return IndFoundArray;	   
	   
  }
  
  
      private int [] findIndex_FromArrays(int[] aData, int aValue) {

	ArrayList<Integer> IndFound = new ArrayList<Integer>();
	
	int i;
	int len_a = aData.length;
	   
	for (i=0; i<len_a; i++) {
		if (aData[i]== aValue) { 
			IndFound.add(i);
		}
	}
	
	// output the index array
	int [] IndFoundArray = new int [IndFound.size()];
	
	for (i=0; i<IndFound.size(); i++) {
		IndFoundArray[i] = IndFound.get(i);
	}
	   
	return IndFoundArray;	   
	   
  }
  

  /**
	(3) Method for 3D ZO CRS travel-time calculation
  */

  private double [] crs_ZO_3D_travelTime(double[] t0, double[] A1, double[] A2, double[] B11, double[] B12, double[] B22, double[] vel, double dx, double dy, double h) {

	// get the length of time samples
	int n_samples = t0.length;
	
	// create output travel time array
	double[] t_crs = new double[n_samples];
	
	// derive the 3D ZO CRS travel time
	for(int i = 0; i < n_samples; i++){
	
		// option-1: travel time calculation.
		// double t_crs2 = Math.pow((t0[i] + dx * A1[i] + dy * A2[i]),2) + B11[i] * Math.pow(dx,2) + B22[i] * Math.pow(dy,2) + 2 * B12[i] * dx * dy + Math.pow((2 * h / vel[i]),2);
		
		// option-2: travel time calculation with neglected the cross-term. (the approximated mixed partial derivative is always very noisy). 
		double t_crs2 = Math.pow((t0[i] + dx * A1[i] + dy * A2[i]),2) + B11[i] * Math.pow(dx,2) + B22[i] * Math.pow(dy,2) + 0 * (2 * B12[i] * dx * dy) + Math.pow((2 * h / vel[i]),2);
		
		if (t_crs2 > 0){	
			t_crs[i]= Math.sqrt(t_crs2);
		}
		else{
			t_crs[i]= Double.NaN;
		}
	} 
	
	return t_crs; 

  }

  /**
	(4) Method for normal moveout travel-time calculation
  */

  private double [] nmo_travelTime(double[] t0, double[] vel, double h) {

	// get the length of time samples
	int n_samples = t0.length;
	
	// create output travel time array
	double[] t_nmo = new double[n_samples];
	
	// derive the NMO travel time
	for(int i = 0; i < n_samples; i++){

		// normal moveout travel time calculation with 2nd order approximation
		t_nmo[i]= Math.sqrt(Math.pow(t0[i],2) + Math.pow((2 * h / vel[i]),2));

	}
	return t_nmo; 
  } 

  /**
	(5) Method for travel-time regularization by replacing the missing travel
  */

  private double [] NanVal_replace(double [] input_array, double [] reference_array){

	// get the length of refernce travel-time array
	int n_samples = input_array.length;
	
	// set the nan value count
	int nanVal_count = 0;
	
	// check if two arrays have identical length
	if(n_samples != reference_array.length){
		System.out.println("Error: input arrays have different length!! ");
	}
	else{
		for(int i = 0; i < n_samples; i++){
			if(Double.isNaN(input_array[i])){
				input_array[i] = reference_array[i];
				
				nanVal_count++;
			}
		}
	}
	return input_array; 
  }

  /**
	(6) derive the nominal shot & receiver coordinates based on mid-point, azimuth, and offset information.
  */
  
  private double [] shot_receiver_coord_generate(double CMP_X, double CMP_Y, double azimuth, double offset){
  	
  	// check if the azimuth is positive value
  	assert azimuth < 0 : "Input azimuth shall be positive value !! ";
  
  	// generate the output array
  	double[] shot_receiver_coord = new double [4];
  	double SOU_X, SOU_Y, REC_X, REC_Y, alpha, alpha_rad;
  	
  	if (azimuth >= 0 && azimuth < 90.0) { 
  		
  		//derive the application angle
  		alpha = 90.0 - azimuth; 
  		
	  	//convert angle to radian
		alpha_rad = Math.toRadians(alpha); 
	  	
	  	// derive the shot's coordinates
	  	SOU_X = CMP_X + offset/2.0 * Math.cos(alpha_rad);
	  	SOU_Y = CMP_Y + offset/2.0 * Math.sin(alpha_rad);
	  	
	  	// derive the receiver's coordinates 
	  	REC_X = CMP_X - offset/2.0 * Math.cos(alpha_rad);
	  	REC_Y = CMP_Y - offset/2.0 * Math.sin(alpha_rad);
  	}
  	else if (azimuth >= 90.0 && azimuth < 180.0) {
  	
  	
  		//derive the application angle
  		alpha = azimuth - 90.0; 
  		
	  	//convert angle to radian
		alpha_rad = Math.toRadians(alpha); 
	  	
	  	// derive the shot's coordinates
	  	SOU_X = CMP_X + offset/2.0 * Math.cos(alpha_rad);
	  	SOU_Y = CMP_Y - offset/2.0 * Math.sin(alpha_rad);
	  	
	  	// derive the receiver's coordinates 
	  	REC_X = CMP_X - offset/2.0 * Math.cos(alpha_rad);
	  	REC_Y = CMP_Y + offset/2.0 * Math.sin(alpha_rad);


  	}
  	else if (azimuth >= 180.0 && azimuth < 270.0) {
  	
  		 //derive the application angle
  		alpha = 270.0 - azimuth; 
  		
	  	//convert angle to radian
		alpha_rad = Math.toRadians(alpha); 
	  	
	  	// derive the shot's coordinates
	  	SOU_X = CMP_X - offset/2.0 * Math.cos(alpha_rad);
	  	SOU_Y = CMP_Y - offset/2.0 * Math.sin(alpha_rad);
	  	
	  	// derive the receiver's coordinates 
	  	REC_X = CMP_X + offset/2.0 * Math.cos(alpha_rad);
	  	REC_Y = CMP_Y + offset/2.0 * Math.sin(alpha_rad);

  	}
  	else {
  	
  		//derive the application angle
  		alpha = azimuth -270.0; 
  		
	  	//convert angle to radian
		alpha_rad = Math.toRadians(alpha); 
	  	
	  	// derive the shot's coordinates
	  	SOU_X = CMP_X - offset/2.0 * Math.cos(alpha_rad);
	  	SOU_Y = CMP_Y + offset/2.0 * Math.sin(alpha_rad);
	  	
	  	// derive the receiver's coordinates 
	  	REC_X = CMP_X + offset/2.0 * Math.cos(alpha_rad);
	  	REC_Y = CMP_Y - offset/2.0 * Math.sin(alpha_rad);
  	
  	}
  	
  	// output the derived coordinates
  	shot_receiver_coord[0] =  SOU_X;
  	shot_receiver_coord[1] =  SOU_Y;
  	shot_receiver_coord[2] =  REC_X;
  	shot_receiver_coord[3] =  REC_Y;
  	
  	return shot_receiver_coord; 
  }
  
  
   /**
	(7) Method for selecting with multiple conditional ranges 
  */
  
      private int [] findIndex_FromArrays_range(int[] aData, int aValueMin, int aValueMax, int[] bData, int bValueMin, int bValueMax  ) {

	ArrayList<Integer> IndFound = new ArrayList<Integer>();
	
	
	int i;
	int len_a = aData.length;
	int len_b = bData.length;
	
	// check if the input arrays have identical length
	assert len_a != len_b : " The input arrays are not in the same length ! ";   
	   
	for (i=0; i<len_a; i++) {
		if (aData[i]>= aValueMin && aData[i]<= aValueMax && bData[i]>= bValueMin && bData[i]<= bValueMax) { 
			IndFound.add(i);
		}
	}
	
	// output the index array
	int [] IndFoundArray = new int [IndFound.size()];
	
	for (i=0; i<IndFound.size(); i++) {
		IndFoundArray[i] = IndFound.get(i);
	}

	return IndFoundArray;	   
  }


  // Any action which must be taken for cleanup of resources should be
  // put in this method.
  @Override
  public void completeNormally() {
  }
  // Any action which must be taken if the tool is told to abort
  // execution should be put in this method.
  @Override
  public void abort() {
  }

} // tool
