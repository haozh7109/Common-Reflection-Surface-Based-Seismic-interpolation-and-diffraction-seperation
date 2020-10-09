
package com.lundin.prowess.tool.fxHaoTest2D;

import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

import com.lgc.prodesk.flowbuilder.pui.ParmProperty;
import com.lgc.prodesk.flowbuilder.pui.pw.PopChooseParm;
import com.lgc.prodesk.flowbuilder.pui.pw.ChooseParm;
import com.lgc.prodesk.flowbuilder.pui.pw.FloatParm;
import com.lgc.prodesk.flowbuilder.pui.pw.IntParm;
import com.lgc.prodesk.flowbuilder.pui.pw.PWParm;
import com.lgc.prodesk.flowbuilder.pui.pw.PWProc;
import com.lgc.prodesk.flowbuilder.pui.pw.LabelParm;
import com.lgc.prodesk.seisdata.DataContext;
import com.lgc.prowess.proj.DbDataset;
import com.lgc.prowess.tool.DatasetParms;
import com.lgc.prodesk.flowbuilder.pui.TypeDesc;
import com.lgc.prodesk.flowbuilder.pui.pw.StringParm;


public class fxHaoTest2DProc extends PWProc  {


	//Parameter: dataset
	static final String DATASET = "DATASET";
	private DatasetParms _datasetParms = new DatasetParms(this, DATASET);


	//thread number
	private static PWParm _warning = new LabelParm("If you get OutOfMemoryError, reduce threads");
	private PWParm _threads = new IntParm(this, "THREADS_PARM", "Number of threads", "Below 1 means all cores. If you get OutOfMemory and cannot allocate more, check the log to see the number of cores on the node and put this below that.", 15);

	//sample ratio in time and space
	private PWParm _time_smp_rate  = new FloatParm(this, "TIME_SAMPLE_PARM", "Time sample ratio (second)", "Time sample ratio", 0.002f);
	private PWParm _space_smp_rate = new FloatParm(this, "SPACE_SAMPLE_PARM", "Space sample ratio (meter)", "Space sample ratio", 6.25f);


	/* CDS parameter: Alpha, Grouped togheter */
	private PWParm _alpha_title = new LabelParm("Set the CDS parameter Alpha (degree))", false, LabelParm.BELOW);

	private PWParm _alpha_start = 
		new FloatParm(this, "ALPHA_START_PARM",
			"alpha min",
			"",
			-20.0f);

	private PWParm _alpha_spacer_1 = new LabelParm();

	private PWParm _alpha_step = 
		new FloatParm(this, "ALPHA_STEP_PARM",
			"alpha step",
			"",
			0.25f);

	private PWParm _alpha_spacer_2 = new LabelParm();

	private PWParm _alpha_end = 
		new FloatParm(this, "ALPHA_END_PARM",
			"alpha max",
			"",
			20.0f);
	/* end */


	/* CDS parameter: Theta, Grouped togheter */
	private PWParm _theta_title = new LabelParm("Set the CDS parameter Theta (degree))", false, LabelParm.BELOW);

	private PWParm _theta_start = 
		new FloatParm(this, "THETA_START_PARM",
			"alpha min",
			"",
			0.0f);

	private PWParm _theta_spacer_1 = new LabelParm();

	private PWParm _theta_step = 
		new FloatParm(this, "THETA_STEP_PARM",
			"theta step",
			"",
			0.0f);

	private PWParm _theta_spacer_2 = new LabelParm();

	private PWParm _theta_end = 
		new FloatParm(this, "THETA_END_PARM",
			"theta max",
			"",
			0.0f);
	/* end */


	//CDS aperture in time (samples)
	private PWParm _aperture_time = new IntParm(this, "CDS_APERTURE_TIME_SAMPLES_PARM", "CDS aperture in time (sample)", "CDS aperture in time samples", 3);
	

	//CDS aperture in space (distance)
	private StringParm _time_control_points = new StringParm(this, "TIME_CONTROL_POINTS_PARM", "CDS aperture in space: time control points (second)", "Enter time control points in second","0,0.3,2.5" );
	private StringParm _aper_control_points = new StringParm(this, "APER_CONTROL_POINTS_PARM", "CDS aperture in space: distance control points (meter)", "Enter aperture control points in meter","0,400,4000" );



	//additional CDS parameters (Water velocity and start time of CDS calculation) 
	private PWParm _water_velocity  = new FloatParm(this, "WATER_VELOCITY_PARM", "Water velocity (meter/second)", "water velocity (meter/second)", 1500.0f);
	private PWParm _cds_cal_start   = new FloatParm(this, "CDS_CALC_START_PARM", "Start time of CDS calucluation (second)", "start time of CDS calucluation", 0.5f);




	private PWParm[] _contents;
	public fxHaoTest2DProc() {

		addParmGroup(_alpha_start, _alpha_spacer_1, _alpha_step, _alpha_spacer_2, _alpha_end);
		addParmGroup(_theta_start, _theta_spacer_1, _theta_step, _theta_spacer_2, _theta_end);

		List<PWParm> items = new ArrayList<>();
		items.addAll(Arrays.asList(_datasetParms.getMenuParms()));
		items.add(_threads);
		items.add(_time_smp_rate);
		items.add(_space_smp_rate);
		items.add(_aperture_time);
		items.add(_time_control_points);
		items.add(_aper_control_points);
		items.add(_water_velocity);
		items.add(_cds_cal_start);
		items.add(new LabelParm());
		items.add(new LabelParm("Alpha"));
		items.add(_alpha_title);
		items.addAll(Arrays.asList(_alpha_start, _alpha_spacer_1, _alpha_step, _alpha_spacer_2, _alpha_end));
		items.add(new LabelParm());
		items.add(new LabelParm("Theta"));
		items.add(_theta_title);
		items.addAll(Arrays.asList(_theta_start, _theta_spacer_1, _theta_step, _theta_spacer_2, _theta_end));
		_contents = items.toArray(new PWParm[0]);
	}

  
	@Override
	public PWParm[] getMenuParms() {
		return _contents;
	}


	@Override
	public void runProcRules(DataContext dataContextIn){
		_datasetParms.runProcRules();
	}

}
