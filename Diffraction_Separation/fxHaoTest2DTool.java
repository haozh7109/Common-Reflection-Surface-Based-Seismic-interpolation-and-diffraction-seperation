
package com.lundin.prowess.tool.fxHaoTest2D;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import java.util.List; 
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.lang.Math;

import com.lgc.gpr.util.ParameterSet;
import com.lgc.prodesk.hdr.CatalogedHdrEntry;
import com.lgc.prodesk.hdr.DimensionHdr;
import com.lgc.prodesk.hdr.HdrCatalog;
import com.lgc.prodesk.seisdata.DataContext;
import com.lgc.prowess.exec.JobContext;
import com.lgc.prowess.javaseis.smart.SmartFramework;
import com.lgc.prowess.seisdata.SeisData;
import com.lgc.prowess.seisdata.SeisUtil;
import com.lgc.prowess.seisdata.SeismicDataset;
import com.lgc.prowess.tool.DatasetParms;
import com.lgc.prowess.tool.InitPhaseException;
import com.lgc.prowess.tool.SimpleTool;
import com.lgc.prowess.tool.ToolContext;

import com.lundin.tools.Matlab;
import edu.mines.jtk.la.DMatrix;
import edu.mines.jtk.util.ArrayMath; 

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import com.lgc.prowess.util.StringArrays;

/**
 * This tool will have to be ThreadTool if you want more complicated behavior, such as
 * interleaving the data that is read sideways with the data that is coming down the flow.
 * <p>fxHaoTest2D
 * This tool reads data sideways and subtracts it from the data that is coming
 * down the flow.
 */

public class fxHaoTest2DTool extends SimpleTool {

  private static final Logger LOG = Logger.getLogger("com.lundin.prowess.tool.fxHaoTest2D");
  private SeismicDataset _dataset;
  private CatalogedHdrEntry _dimensionHdrHypercube;
  private CatalogedHdrEntry _dimensionHdrVolume;
  private CatalogedHdrEntry _dimensionHdrFrame;
  private CatalogedHdrEntry _dimensionHdrTrace;
  private CatalogedHdrEntry _sidewaysDimensionHdrTrace;
  private HashMap<Integer,float[]> _traceMap = new HashMap<Integer,float[]>();

  private int _THREADS;
  private double _TIME_SAMPLE_PARM;
  private double _SPACE_SAMPLE_PARM;	
  private double _ALPHA_START; 
  private double _ALPHA_STEP; 
  private double _ALPHA_END; 
  private double[] _ALPHA; 
  private double _THETA_START; 
  private double _THETA_STEP; 
  private double _THETA_END; 
  private double[] _THETA;
  private int _CDS_APERTURE_TIME_SAMPLES;
  private double[] _TIME_CONTROL_POINT;
  private double[] _APER_CONTROL_POINTS_STR;
  private double _WATER_VELOCITY;
  private double _CDS_CALC_START;

  @Override
  public DataContext init(JobContext jobContext, ToolContext toolContext,
                          DataContext dataContext) throws InitPhaseException {

    // Get menu parameters	
    ParameterSet menuParms = toolContext.getMenuParms();

    _THREADS 		= menuParms.getInt("THREADS_PARM", 15);
    _TIME_SAMPLE_PARM	= menuParms.getFloat("TIME_SAMPLE_PARM", 0.002f);
    _SPACE_SAMPLE_PARM	= menuParms.getFloat("SPACE_SAMPLE_PARM", 6.25f);	

    _ALPHA_START 	= menuParms.getFloat("ALPHA_START_PARM", -20.0f);
    _ALPHA_STEP		= menuParms.getFloat("ALPHA_STEP_PARM", 0.25f);
    _ALPHA_END 		= menuParms.getFloat("ALPHA_END_PARM", 20.0f);	
    _ALPHA 		= new double[]{_ALPHA_START, _ALPHA_STEP, _ALPHA_END};

    _THETA_START 	= menuParms.getFloat("THETA_START_PARM", 0.0f);
    _THETA_STEP		= menuParms.getFloat("THETA_STEP_PARM", 0.0f);
    _THETA_END 		= menuParms.getFloat("THETA_END_PARM", 0.0f);	
    _THETA 		= new double[]{_THETA_START, _THETA_STEP, _THETA_END};

    _CDS_APERTURE_TIME_SAMPLES	= menuParms.getInt("CDS_APERTURE_TIME_SAMPLES_PARM",3);	

    String TIME_CONTROL_POINTS_STR  	= menuParms.getString("TIME_CONTROL_POINTS_PARM", "0,0.3,2.5");
    _TIME_CONTROL_POINT			= StringArrays.stringToDoubleArray(TIME_CONTROL_POINTS_STR);

    String APER_CONTROL_POINTS_STR	= menuParms.getString("APER_CONTROL_POINTS_PARM", "0,400,4000");
    _APER_CONTROL_POINTS_STR		= StringArrays.stringToDoubleArray(APER_CONTROL_POINTS_STR);


    _WATER_VELOCITY	= menuParms.getFloat("WATER_VELOCITY_PARM", 1500.0f);
    _CDS_CALC_START	= menuParms.getFloat("CDS_CALC_START_PARM", 0.5f);


    //ParameterSet menuParms = toolContext.getMenuParms();
    //String[] filePaths = DatasetParms.getFilePaths(menuParms, jobContext.getAreaLinePath());
    //String[] filePaths = DatasetParms.getFilePaths(menuParms, fxHaoTest2DProc.DATASET,jobContext.getAreaLinePath());

    String[] filePaths = DatasetParms.getFilePaths(menuParms, fxHaoTest2DProc.DATASET,jobContext.getAreaLinePath());

    String filePath = (filePaths.length > 0) ? filePaths[0] : null;
    if (filePath == null) throw new InitPhaseException("No dataset was selected");

    try {
      _dataset = SeisUtil.openDataset(filePath);
      LOG.info("Dataset Descriptive Name=" + _dataset.getDescName());
    } catch (IOException e) {
      throw new InitPhaseException(e);
    }

    /* TODO: Is this check necessary?
    SmartFramework smartFramework = _dataset.getDataContext().getSmartFramework();
    if (smartFramework == null)
      throw new InitPhaseException("Dataset '" + filePath
                                   + "' does not have a well-formed JavaSeis framework");
    */

    HdrCatalog hdrCatalog = dataContext.getHdrCatalog();

    DimensionHdr[] inputDimensionHdrs = dataContext.getSeisGrid().getHdrs();

    // ***** CODE ADDED HERE *****
    _dimensionHdrTrace = hdrCatalog.getEntry(inputDimensionHdrs[1].getHdrEntry());
    if (_dimensionHdrTrace == null)  // This should never happen.
      throw new InitPhaseException("Header entry '" + inputDimensionHdrs[1].getHdrEntry()
                                   + "' is missing from the input data");

    _dimensionHdrFrame = hdrCatalog.getEntry(inputDimensionHdrs[2].getHdrEntry());
    if (_dimensionHdrFrame == null)  // This should never happen.
      throw new InitPhaseException("Header entry '" + inputDimensionHdrs[2].getHdrEntry()
                                   + "' is missing from the input data");

    // Note that we have to check to see if the dataset is more than 3-dimensional.
    if (inputDimensionHdrs.length > 3) {
      _dimensionHdrVolume = hdrCatalog.getEntry(inputDimensionHdrs[3].getHdrEntry());
      if (_dimensionHdrVolume == null)
        throw new InitPhaseException("Header entry '" + inputDimensionHdrs[3].getHdrEntry()
                                     + "' is missing from the input data");
    } else {
      // This is a 3-dimensional dataset.  No problem.
      _dimensionHdrVolume = null;
    }

    // Note that we have to check to see if the dataset is more than 4-dimensional.
    if (inputDimensionHdrs.length > 4) {
      _dimensionHdrHypercube = hdrCatalog.getEntry(inputDimensionHdrs[4].getHdrEntry());
      if (_dimensionHdrHypercube == null)
        throw new InitPhaseException("Header entry '" + inputDimensionHdrs[4].getHdrEntry()
                                     + "' is missing from the input data");
    } else {
      // This is a 3-dimensional or 4-dimensional dataset.  No problem.
      _dimensionHdrHypercube = null;
    }

    // ***** CODE ADDED HERE *****
    try {
      HdrCatalog sidewaysHdrCatalog = _dataset.getDataContext().getHdrCatalog();
      _sidewaysDimensionHdrTrace = sidewaysHdrCatalog.getEntry(inputDimensionHdrs[1].getHdrEntry());
      if (_sidewaysDimensionHdrTrace == null)
        throw new InitPhaseException("Header entry '" + inputDimensionHdrs[1].getHdrEntry()
                                     + "' is missing from the 'sideways' data");
    } catch (IOException e) {
      e.printStackTrace();
      throw new InitPhaseException(e);
    }

    // No change to the input data context.
    return dataContext;
  }


  @Override
  public SeisData exec(SeisData data) {

    // Get the volume/frame for this ensemble.  We will assume that the first trace is
    // representative.
    int[] hdr0 = data.getHdr(0);

    long volumeKey;
    if (_dimensionHdrVolume == null) {
      volumeKey = SmartFramework.IMPLIED_LOGICAL_VALUE;
    } else {
      volumeKey = _dimensionHdrVolume.getLongVal(hdr0);
    }

    long hypercubeKey;
    if (_dimensionHdrHypercube == null) {
      hypercubeKey = SmartFramework.IMPLIED_LOGICAL_VALUE;
    } else {
      hypercubeKey = _dimensionHdrHypercube.getLongVal(hdr0);
    }

    long frameKey = _dimensionHdrFrame.getLongVal(hdr0);

    SeisData sidewaysData = null;
    try {
      sidewaysData = _dataset.readFrame(hypercubeKey, volumeKey, frameKey);
    } catch (IOException e) {
      LOG.warning("Error reading hypercube=" + hypercubeKey + " volume=" + volumeKey
                  + " frame=" + frameKey
                  + " in the sideways dataset");
      sidewaysData = null;  // For clarity.
    }

    if (sidewaysData == null) {
      // Usually that just means that this data is not present.
      if (LOG.isLoggable(Level.FINE)) {
        LOG.fine("Unable to read hypercube=" + hypercubeKey + " volume=" + volumeKey
                 + " frame=" + frameKey
                 + " in the sideways dataset");
      }

    } else {
      // We have both data and sidewaysdata, proceed to process now.
      
      // ***** CODE ADDED HERE *****
	

      	// Setup for multithreading
	int _threads = _THREADS;
	int cores = _THREADS; //Runtime.getRuntime().availableProcessors();
	int numberOfThreads = _threads < 1 ? cores : _threads;
	LOG.info("Node has "+cores+" cores and will use " + numberOfThreads + " threads");
	
	//ExecutorService service = Executors.newFixedThreadPool(cores);
	final BlockingQueue<Runnable> queue = new ArrayBlockingQueue<>(500);
	ExecutorService service = new ThreadPoolExecutor(numberOfThreads, numberOfThreads, 0L, TimeUnit.MILLISECONDS, queue, new ThreadPoolExecutor.CallerRunsPolicy());


      // ------ get the stacked 2d data and velocity file  ------
      float[][] data2d 	= data.getTraces();
      float[][] vels2d 	= sidewaysData.getTraces();

      // ------ get the size of the 2d frame  ------
      int data_Ntsamps 	= data.getTraceLen();
      int data_Ntrs   	= data.countTraces();

      // ------ define the temporay variables which link to the input predefined parameters (to be set in the promax window) ------
      final double dt	= _TIME_SAMPLE_PARM;   // sample ratio in time (ms)
      final double dx	= _SPACE_SAMPLE_PARM;  // sample ratio in space (m)


      // ------ get the actual size of the 2d seismc and data, assumes that the 2d Frame is triple length in the input  ------
      int data_Ntsamps_seis = (int) (data_Ntsamps / 3); //get the seismic data length
	
      float[][] data2d_seis =  new float [data_Ntrs][data_Ntsamps_seis];
      for (int iR = 0; iR<data_Ntrs; iR++){ 
	for (int iC = 0; iC<data_Ntsamps_seis; iC++){ 
      		data2d_seis[iR][iC] = (float)data2d[iR][iC];	
	}
      }


      
      double [] data_midpoints = new double [data_Ntrs]; // define the relative spatial coordinates of input data.
      for (int n = 0;n<data_Ntrs;n++) {
      data_midpoints[n] = n * dx;
      }
	
      //define the spatial coordinates of output data. (CDS will be calculated on such grids) 	
      double [] cds_midpoints = data_midpoints; 

      // === or select part of data for CDS test ==
      //int test_cds_locations    = 30;
      //double [] cds_midpoints 	= new double [test_cds_locations];
      //for (int iloc = 0; iloc < test_cds_locations; iloc++){
      //		cds_midpoints[iloc] = data_midpoints[iloc+(data_Ntrs/2)];     // pick locations for the CDS calculation test.	
      //} 

      // define the CDS aperture in mid-point coordinate.
      int[] contPoint_t = new int [_TIME_CONTROL_POINT.length]; 
      for (int m = 0;m<_TIME_CONTROL_POINT.length;m++) {
      	contPoint_t[m] = (int)(_TIME_CONTROL_POINT[m]/_TIME_SAMPLE_PARM);
      }
      
      double[] contPoint_aptx 	= _APER_CONTROL_POINTS_STR;

      
      // define the time-sample series
      int [] data_samples = new int [data_Ntsamps_seis]; 
      for (int m = 0;m<data_Ntsamps_seis;m++) {
      	data_samples[m] = m;
      }
      
      // derive the time-variant CDS spatial aperture by linear interpolation 	
      double [] apert_x = new double [data_Ntsamps_seis]; 
      apert_x 		=  Matlab.interpolateLinear(contPoint_aptx, contPoint_t, data_samples); 
      
      // define the CDS temporal aperture (samples, total temporal aperture size is 2*apert_t+1) .
      final int apert_t = _CDS_APERTURE_TIME_SAMPLES; 

      // ------ define the CDS parameter searching space: emergencing angle(alpha) and structure dips (theta) ------
      // define the searching parameter space for alpha:
      double alpha_min 	= _ALPHA[0];  // the min alpha angle
      double alpha_step	= _ALPHA[1];  // the alpha increment
      double alpha_max  = _ALPHA[2];  // the max alpha angle
      
      int alpha_k	= (int)((alpha_max - alpha_min)/alpha_step) + 1;
      double [] alpha   = new double [alpha_k];
 
      for (int k = 0;k<alpha_k;k++) {    
      alpha[k] = alpha_min + k * alpha_step;	
      }

      // define the searching parameter space for theta:
      double theta_min 	= _THETA[0]; 	// the min dip angle
      double theta_step	= _THETA[1];	// the dip increment   
      double theta_max  = _THETA[2]; 	// the max dip angle
     
      int theta_k	= (int)((theta_max - theta_min)/theta_step) + 1;
      double [] theta   = new double [theta_k];
 
      for (int k = 0;k<theta_k;k++) {    
      theta[k] = theta_min + k * theta_step;	
      }

      // ------ define the misc. parameters for CDS calculation ------
      final double vel_shallow 	= _WATER_VELOCITY; // define the shallow layer's velocity for CDS calculation
      final int tCal_Samp_start	= (int) (_CDS_CALC_START / _TIME_SAMPLE_PARM);  // define the begining time sample for CDS calculation


      // ------  CDS search based on semblance calculation   ------ 

      // define the variables to save the output from calculation
      final double [][] SemblancePara  	= new double[data_Ntrs][data_Ntsamps_seis];
      final double [][] A_para 	  	= new double[data_Ntrs][data_Ntsamps_seis];
      final double [][] C_para 	  	= new double[data_Ntrs][data_Ntsamps_seis];

      // define the variables to save all the output from calculation	
      final double [][] CDS_para 	= new double[data_Ntrs][data_Ntsamps];	

      // semblance based common diffraction surface (CDS) search, the multi-threads parallel executing.
	int total_cmps = cds_midpoints.length;


	for (int ICMP = 0;ICMP<total_cmps;ICMP++){

		//define the fixed variables		
		final int icmp 				= ICMP;
		final float[][] seis2d		     	= data2d_seis; 
		final float[][] seis_vel2d		= vels2d;
		final double [] seis_midpoints		= data_midpoints;
		final double [] cds_location		= {cds_midpoints[icmp]};

		final double [] seis_apert_x		= apert_x;
		final int seis_apert_t			= apert_t;
		final double [] seis_alpha   		= alpha;
		final double [] seis_theta   		= theta;

		//find the cmp location of current CDS 
		ArrayList<Integer> ind_cmp_list = Matlab.find(data_midpoints,cds_midpoints[icmp]);
	       	final int ind_cmp = ind_cmp_list.get(0);
 


		//prarallel running		
		service.execute(new Runnable() {
                	@Override
                	public void run(){ 
				double [][] CDSCal_tmp = midpoint_CDS_search(seis2d, seis_midpoints, dt, cds_location, seis_apert_x, seis_apert_t, seis_alpha, seis_theta, seis_vel2d, vel_shallow, tCal_Samp_start);
					for (int it = 0, Nt = CDSCal_tmp[0].length; it<Nt; it++){				
						CDS_para[ind_cmp][it] = CDSCal_tmp[ind_cmp][it];
					} 
      			}


		});
	}


      	try {
		service.shutdown();
		service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	} catch (InterruptedException e) {
		e.printStackTrace();
	}


      // Output by replacing the data to calculated semblance
      for (int iR = 0; iR<data_Ntrs; iR++){ 
		for (int iC = 0; iC<data_Ntsamps; iC++){ 
			data2d[iR][iC] = (float)CDS_para[iR][iC];	
		}
      }
	

      // ***** DISPLAY THE CHECKED VARIABLES IN LOG *****
      //LOG.info("================ START DISPLAY THE CHECKED VARIABLES IN LOG FILE ================================");
      //LOG.info("Found the dimension of input ensemble:" + data_Ntsamps_seis + "*" + data_Ntrs + " in the test job");
      //LOG.info("================ COMPLETE DISPLAY THE CHECKED VARIABLES IN LOG FILE =============================");


      // // Be sure that we don't keep a reference to any traces.
      // _traceMap.clear();

      // The reference to the sideways data is about to be lost - we have
      // to free it first.
      sidewaysData.free();
    }

    return data;

	   
  }




/**
===================================
Method for CDS parameter searching
===================================
*/

private double[][] midpoint_CDS_search(float[][] data2d, double[] data_midpoints, double dt, double[] cds_midpoints, double[] apert_x, int apert_t, double[] alpha, double theta[], float[][] vels2d, double vel_shallow, int tCal_Samp_start ) {
         

	//get dimensions of data
	int size_x 	= data2d.length;
	int size_t 	= data2d[0].length;
	int size_x_cds 	= cds_midpoints.length;
        
	//define variables to save outputs from the calculation
	//outputs at cds locations (if CDS is applied on decimated grid) 
	double [][] semblance_para_tmp 	= new double [size_x_cds][size_t];
	double [][] A_para_tmp 		= new double [size_x_cds][size_t];
	double [][] C_para_tmp 		= new double [size_x_cds][size_t];

	//outputs at data locations (if CDS is applied on data grid) 
	double [][] semblance_para 	= new double [size_x][size_t];
	double [][] A_para 		= new double [size_x][size_t];
	double [][] C_para 		= new double [size_x][size_t];
	

	// ------ start the CDS parameters searching, loop-1 through spatial coordinates  ------ 
	for (int ix = 0; ix<size_x_cds; ix++) {    
			
			//get the x-coordinate of searching location 
			double x0 = cds_midpoints[ix];
	      	
			//get the index of x0 in the data 
			ArrayList<Integer> ind_x0_list = Matlab.find(data_midpoints,x0);
	       		int ind_x0 = ind_x0_list.get(0);
		
		// ------ loop-2 through temporal coordinates  ------ 
		for (int it = 0; it<size_t; it++) { 
		//for (int it = 100; it<size_t;  it=it+100) {
			if (it >= tCal_Samp_start) {

				double t0 = (it+1) * dt;
				// ***** DISPLAY THE LIVE PROGRESS LOG ***** (MUST BE REMOVED IN PRODUCTION!!)	
				//LOG.info("++++++++++++++ The CDS is calculating at CMP location: " + ind_x0 + " CMP coordinate: " + " at CMP coordinate: " +  cds_midpoints[ix] + " with time sample: " + it + " in the job +++++++++++++++++");
				if ((it%10) == 0) {				
					LOG.info("++++++++++++++ The CDS parameter searching is calculating at CMP location: " + ind_x0 + " with time sample: " + it + " in the job +++++++++++++++++");
				}

				//extract the traces within the midpoint aperture.
				ArrayList<Integer> ind_x_apt_AList = Matlab.findBetween(data_midpoints,x0 - apert_x[it],x0 + apert_x[it]);
				int ind_bound_min     = (int)ind_x_apt_AList.get(0);
				int ind_bound_max     =  ind_bound_min + ind_x_apt_AList.size() - 1;
				double[][] data2d_apt = Matlab.subset(data2d, ind_bound_min, ind_bound_max, 0, size_t-1);


				//get the x-coordinates of selected data
				double [] midpoints_apt = Matlab.select(data_midpoints,ind_bound_min,ind_bound_max );
			
				//LOG.info("----------- The midpoints_apt has defined aperture range: " + midpoints_apt[0] + "to" + midpoints_apt[midpoints_apt.length-1] + " in the job -----------------");

				//define the local semblance matrix (dimension: N_alpha x N_theta) 
				double[][] semblance_local 	= new double [alpha.length][theta.length];
				double[][] A_local 		= new double [alpha.length][theta.length];
				double[][] C_local 		= new double [alpha.length][theta.length];

				// ------ loop-3 through A parameters  -----
				for (int ia = 0, size_a = alpha.length; ia<size_a; ia++){
			
			
					//convert angle to radian
					double alpha_rad = Math.toRadians(alpha[ia]); 
			

					//derive A parameter			
					double A = (2*Math.sin(alpha_rad))/vel_shallow;
				
					// ------ loop-4 through B parameters  -----
					for (int ib = 0, size_b = theta.length; ib<size_b; ib++){

						//convert angle to rad
						double theta_rad = Math.toRadians(theta[ib]);
					
						//derive C parameter
						//double C =  (4./ (Math.pow(vels2d[ix][it] * Math.cos(theta_rad),2))) *  (1-Math.pow((Math.sin(alpha_rad)),2));
						double C =  (4./ (Math.pow(vels2d[ind_x0][it] * Math.cos(theta_rad),2))) *  (1-Math.pow((Math.sin(alpha_rad)),2));					
						// calculate semblances based on current CDS parameters (A and C)
						semblance_local[ia][ib] = CDS_semblance(data2d_apt, A, C, midpoints_apt, x0, t0, apert_t, dt);	
						A_local[ia][ib] = A;
						C_local[ia][ib] = C;
					}
				}
		
				// derive the calculated best parameters
				double[] max_ind 	= Matlab.max2d(semblance_local);
				int max_indrow		= (int)max_ind[1]; 		
				int max_indCol		= (int)max_ind[2];
		
				// assign derived values to the tmp outputs
				semblance_para_tmp[ix][it]	=  semblance_local[max_indrow][max_indCol];
			        A_para_tmp[ix][it]		=  A_local[max_indrow][max_indCol]; 
				C_para_tmp[ix][it]		=  C_local[max_indrow][max_indCol];

			}

			else {
				// assign derived values to the tmp outputs
				semblance_para_tmp[ix][it]	=  0;
		        	A_para_tmp[ix][it]		=  0; 
				C_para_tmp[ix][it]		=  0;
			}

			// assign derived values to the final outputs
			semblance_para[ind_x0][it]	=  semblance_para_tmp[ix][it];
			A_para[ind_x0][it]		=  A_para_tmp[ix][it]; 
			C_para[ind_x0][it]		=  C_para_tmp[ix][it];			
					
		}

	}

	//combine 3 parameters to one frame and return as the output of CDS parameter search. 
	double [][] CDS_para 	= new double [size_x][size_t*3]; // define the output

						
	for (int ix = 0; ix<size_x; ix++){
		for (int it = 0, Nt = size_t*3; it<Nt; it++){	
			if (it < size_t) {	
		
				CDS_para[ix][it] = semblance_para[ix][it];
			}
			else if(it >= size_t && it < 2*size_t ){

				CDS_para[ix][it] = A_para[ix][it-size_t];

			}
			else{
				CDS_para[ix][it] = C_para[ix][it-(size_t*2)];
			}
			
		}
	} 	

	return CDS_para;

 }




/**
Method for semblance parameter calculation. 
*/

private double CDS_semblance(double[][] data2d_apt, double A, double C, double[] midpoints_apt, double x0, double t0, int apert_t, double dt) {

	//define the semblance variable
	double semblance_val = 0; 

	//get dimensions of data
	int size_x 	= data2d_apt.length;
	int size_t 	= data2d_apt[0].length;

	//define the 2d matrix to save the extracted seismic amplitude
  	double[][] amp  = new double[apert_t*2+1][size_x];

	//define the trace count
	int count_n = 0;

	// ------ loop-1 through traces inside aperture  -----
	for (int ix = 0; ix<size_x; ix++) {
		
		//get the x-coordinate
		double x  = midpoints_apt[ix]; // the length of midpoints_apt must equal to size_x 
		
		//get the constant values for inner loops
		double dm   	= x - x0;
		double a_dm 	= A * dm;
		double c_dm_2 	= C * Math.pow(dm,2);

		// ------ loop-2 through time samples inside aperture  -----
		for (int it = -apert_t, t_upper_bound = apert_t+1; it<t_upper_bound; it++) {

			//get the reference value t
			double t0_new = t0 + it*dt;
			
			
			//get the cds travel time
			double t_cds  = Math.sqrt(Math.pow((t0_new + a_dm),2) + c_dm_2);

			//get the amplitude by linear interpolation
			double smp_loc  = t_cds/dt;
			int smp_loc_up  = (int)Math.ceil(smp_loc);   //upper boundary of sample location
			int smp_loc_low = (int)Math.floor(smp_loc);  //lower boundary of sample location

			//linear interpolation 
			if (smp_loc_low >= 0 && smp_loc_up < size_t){

				//define the weight
				double weight 		=  smp_loc - (double)smp_loc_low;
				amp[it + apert_t][ix] 	=  data2d_apt[ix][smp_loc_up] * weight + data2d_apt[ix][smp_loc_low] * (1-weight);
			}

		}
		//update trace count
		count_n++;
	}


	// ------  Apply hanning window and calculate the semblance ------------
	double [][] weights_hannWindow2d = HannWindow2d(amp.length,amp[0].length);

	//multiplication of extracted amplitude with 2D hanning window
	double[][] amp_taped = ArrayMath.mul(amp,weights_hannWindow2d);	

	//derive the semblance value
	semblance_val  = Semblance_calculation(amp_taped,count_n);

	//output the calculated semblance value
	return semblance_val;

}

/**
Method for 1D hanning window calculation (output a 1D array)
*/

private double [] HannWindow(int length_n){

	//define variables for the weights 
	double [] wind_weights = new double [length_n];
	
	//calculate the window function values
	for (int iN = 0; iN < length_n ; iN++){

		wind_weights[iN] = 0.5 * ( 1 - Math.cos(2 * Math.PI * ((double)iN /(length_n-1)))); //be aware of the integter to double conversion. 

	}
	
	//output weights
	return wind_weights;
}


/**
Method for 2D hanning window calculation (output a 2D array)
*/
private double [][] HannWindow2d(int length_n, int length_m){

	//define variables for the weights 
	double [][] wind2d_weights = new double [length_n][length_m];
	
	//calculate the 2D window function values
	for (int iN = 0; iN < length_n ; iN++){

		for (int iM = 0; iM < length_m ; iM++){

			wind2d_weights[iN][iM] = 0.5 * ( 1 - Math.cos(2 * Math.PI * ((double)iN /(length_n-1)))); //be aware of the integter to double conversion.  
		}

	}
	
	//output weights
	return wind2d_weights;
}


/**
Method for semblance value calculation
*/
private double Semblance_calculation(double[][] AmpArray, int countN){	

	// get the sizes of the 2d amplitude array
	double NR = AmpArray.length;
	double NC = AmpArray[0].length;

	// calculate the numerator of the semblance fraction
	double semb_numerator = 0;

	for (int iR = 0; iR < NR ; iR++){
		double row_sum = 0;
		for (int iC = 0; iC < NC ; iC++){
			row_sum = row_sum + AmpArray[iR][iC];
		}
		semb_numerator = semb_numerator + Math.pow(row_sum,2);  
	}

	// calculate the denorminator of the semblance fraction
	double semb_denorminator = 0;
	for (int iR = 0; iR < NR ; iR++){
		for (int iC = 0; iC < NC ; iC++){
			semb_denorminator = semb_denorminator + Math.pow(AmpArray[iR][iC],2);
		}
	}

	// calculate the semblance value
	double Semblance_Val = (semb_numerator  / semb_denorminator) / countN; //apply the normalization fraction: 1/countN, as the CDS apterture is increasing with depth)

	return Semblance_Val;
}

  private void close() {
    try {
      _dataset.close();
    } catch (IOException e) {
      // Silently eat the exception.
    }
  }


  @Override
  public void completeNormally() {
    this.close();
  }

  @Override
  public void abort() {
    this.close();
  }

}
