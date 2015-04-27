/****************************************************************************

Author/Copyright:      	Patrik Steinhoff

Creation Date:	       	End of 2012
	
Description:	       	creates .vtk files from tracks, detected particles
						and sorts them into stacks
	
Routines contained:     PostToParaview_c,

****************************************************************************/
#include "ptv_is_2_unstructured_vtk.h"
/*Define some global variables:*/
int	seq_first, seq_last;	       			/* 1. and last img of seq */
double  X_lay[2], Zmin_lay[2], Zmax_lay[2];	/* illu. layer data */
trackparameters tpar;						/* tracking parameters */
FILE  *fpp;
char   	seq_name[4][128];      				/* sequence names */
char	path_to_WF[256], pathBuffer[256];
char	velo_par[256], sequence_par[256],track_par[256],criteria_par[256],point_offset[256];

int main(int argc,char **argv)
    //int argc;			/* Number of command-line arguments. */
    //char **argv;		/* Values of command-line arguments. */
{
	char Buffer[80];
	int count, size;
	int start;
	int Exportflag[4];
	int i;
	FILE *fpp;
	sprintf(path_to_WF,"D://EBC0043_DFG_3D-PTV_pst/WorkingFolder_EBC/");//adjust your default path to your working-folder here.
	strcpy(velo_par, "parameters/velo_trans.par");
	strcpy(sequence_par, "parameters/sequence.par");
	strcpy(track_par, "parameters/track.par");
	strcpy(criteria_par, "parameters/criteria.par");
	strcpy(point_offset, "parameters/3d_point_offset.par");
	 printf("==============================================\n"); 
	 printf("\\_\\_\\_\\_\\_\\_\\_\\_      \\_	Particle\n"); 
	 printf("\\_   \\_  \\_    \\_    \\_		Tracking\n"); 
	 printf("\\_\\_\\_   \\_     \\_  \\_		Velocimetry\n"); 
	 printf("\\_       \\_      \\_\\_  	\n"); 
	 printf("\\_       \\_       \\_  \n"); 
	 printf("==============================================\n"); 

	 printf("==============================================\n");
	 printf("\\_\\_\\_ \\_\\_  \\_\\_\\_	Energy efficient\n");
	 printf("\\_     \\_  \\ \\_		Buildings and indoor	\n");
	 printf("\\_\\_\\_ \\_\\_\\ \\_		Climate,	\n");
	 printf("\\_     \\_  \\ \\_		RWTH Aachen	\n");
	 printf("\\_\\_\\_ \\_\\_\\ \\_\\_\\_	University	\n");
	 printf("==============================================\n");
	 printf("===  Export ptv-data to .vtk data format:  ===\n");
	 printf("==============================================\n");
	 printf("Enter -h for additional information.\n");
	 printf("Please select your option(s):\n");
	 gets(Buffer);
	 size = strlen(Buffer);
	 for(i=0;i<4;i++)Exportflag[i]=0;
	 start = 0;
	 while (start !=1)
	 {
		 count = 0;
		 if(size ==0) start=1;
		 while(size>count) 
		 {
			switch( Buffer[count] ) {
				case 'a':
				fprintf( stderr, "Enabling export of 3d-points.\n" );
				//flag for 3d points
				Exportflag[0]=1;
				start = 1;
				break;

				case 'b':
				fprintf( stderr, "Enabling export of tracked points with Ux,Uy,Uz.\n" );
				//flag for tracked points
				Exportflag[1]=1;
				start = 1;
				break;

				case 'c':
				fprintf( stderr, "Enabling export of trajectories.\n" );
				//flag for trajects
				Exportflag[2]=1;
				start = 1;
				break;

				case 'd':
				fprintf( stderr, "Enabling export of stacks.\n");
				//flag for stacks
				Exportflag[3]=1;
				start = 1;
				break;

				case 'h':
				default:
				printf( "\nUsage: -[abcdh]\n\n"
						"Options:\n"
						"\n**************************************************************\n"
						"\nEBC implemented a export function from ptv.is files into vtk data format,\n which works well with the postprocessing tool PARAVIEW.\n\n"
						"Ensure that there is enough memory (RAM) available.\n\n"
						"   -a"
						"        Enable export of 3d coordinates from all detected points\n (rt_is.number).\n\n"
						"   -b"
						"        Enable export of tracked points with Ux,Uy,Uz.\n\n"
						"   -c"
						"        Enable export of trajectories.\n"
						"        This contains to types:/\n"
						"        (i)Export of trajects between two consecutive images.\n"
						"        (ii)Export of trajects in to folders according to their length\n\n"	
						"   -d"
						"        Enable export of stacks of trajectories.\n\n"
						"   -h"
						"        Prints this usage informaton.\n\n"
					  );
					
			}
			count+=1;
		}
		 sprintf(Buffer, "");
		 if(start == 0)
		 {
			printf("Please select your option(s):\n");
			gets(Buffer);
			size = strlen(Buffer);
		 }
		 if (start ==1)
		 {
			 //ToDo: Read config file in same directory as exe-file and ask for path
			 printf("Please press \"Return\" to start export with configuration from ptc2vtk.config.\n");
			 printf("Or enter path to your working folder and use manual settings:\n");
			 printf("Format:\n");
			 printf("Volume://folder1/folder2/../YourWorkingFolder/\n\n");
			 gets(pathBuffer); //not used yet.
			 if(strlen(pathBuffer)==0)
			 {
					strcat(pathBuffer, "ptv2vtk.config");
					fpp = fopen(pathBuffer,"r");
					if(fpp == 0)
					{
						printf("Could not find file: ptv2vtk.config\n");
						start = 0;
					}
					else
					{			
						for(i=0; i<4; i++) 
						{
							fscanf(fpp, "%d\n", &Exportflag[i]);
						}
						fscanf (fpp,"%s\n", &path_to_WF);
						fclose (fpp);
					}
			 }
		 }
		if (start ==1)
		 {
			 PostToParaview_c(argc,argv, path_to_WF, Exportflag);
		 }
	}
    //PostToParaview_c(argc,argv, path_to_WF, Exportflag);
	printf("Press any key to quit application.");
	gets(Buffer);
    return 0;			/* Needed only to prevent compiler warning. */
}

int PostToParaview_c( int argc, const char** argv, char base_path[256],int flag4export[4])
{
  int i,anz, anz1, anz2,anz3, m, j, k,l, q, r, tmp_ID, nOfP, nVerts,nPoints, time, timeID, lineOfP, tmpVeloPos, minSelect;
  int x_offset, y_offset, z_offset;
  int dummy;
  FILE *fp1, *fp2,*fp22, *fp88, *fp66, *fp77, *fp99,*fp55, *fp111;
  char val[256],val1[256], val2[256], val3[256],val4[256];
  char path[256], selectionPath[256], tmp_path[512],tmp_file_name[256];
  vector *line1, *line2, *line3, *veloPos1;
  double color, ymin=0, ymax=0, cubes;
  double mx, my, mz, dx, dy, dz, du, ux, uy,uz, dl, rotz, rotx;
  double Ux,Uy,Uz,Um;
  int connector, nStack, stacksize;
  double MeanU, tmpValue;
  float interframeTime;
  double U_meanMax;
  int stackOnly = 0;

  readseqtrackcrit(base_path);
  U_meanMax = tpar.Vmax;

  strcpy(tmp_path, base_path);
  strcat(tmp_path, velo_par);
	 fp77=fopen(tmp_path,"r");
	 fscanf(fp77, "%f\n", &interframeTime);
	 fscanf(fp77, "%d\n", &minSelect);
	 fscanf(fp77, "%d\n", &stacksize); 
	 //fscanf(fp77, "%d\n", &stackOnly); 
	 fclose(fp77);

  strcpy(tmp_path, base_path);
  strcat(tmp_path, point_offset);
	 fp22= fopen(tmp_path,"r");
	 fscanf (fp22, "%d\n", &x_offset);
	 fscanf (fp22, "%d\n", &y_offset);
	 fscanf (fp22, "%d\n", &z_offset);
	 fclose(fp22);
  //************************************************************************************//
  //create vtk frame file in xml-style:
   /* open file for line elements */
  strcpy(tmp_path, base_path);
  strcpy(tmp_file_name, "PostProcessing/");
  strcat(tmp_path, tmp_file_name);
  CreateDirectoryA(tmp_path,NULL);

  strcpy(tmp_path, base_path);
  strcpy(tmp_file_name, "PostProcessing/Frame.vtu");
  strcat(tmp_path, tmp_file_name);
  fp88 = fopen (tmp_path, "w");
  fprintf(fp88, "<?xml version=\"1.0\"?> \n");
  fprintf(fp88, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">   \n");
  fprintf(fp88, "<UnstructuredGrid>   \n");
  fprintf(fp88, "<Piece NumberOfPoints=\"8\" NumberOfCells=\"1\" >   \n");
  fprintf(fp88, " <PointData Scalars=\"scalars\">  \n");
  fprintf(fp88, " <DataArray type=\"Float32\" Name=\"outVect\" NumberOfComponents=\"3\" format=\"ascii\">  \n");
  fprintf(fp88, " -1 -1 -1 \n -1 -1 0 \n -1 0 -1 \n -1 0 0 \n 0 -1 -1 \n 0 -1 0 \n 0 0 -1 \n 0 0 0 \n");
  fprintf(fp88, "  </DataArray> \n");
  fprintf(fp88, "  </PointData> \n");
  fprintf(fp88, "  <CellData Scalars=\"scalars\"> \n");
  fprintf(fp88, "  <DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\"> \n");
  fprintf(fp88, "  0 \n");
  fprintf(fp88, "  </DataArray> \n");
  fprintf(fp88, "  <DataArray type=\"Float32\" Name=\"cVects\" NumberOfComponents=\"3\" format=\"ascii\"> \n");
  fprintf(fp88, "  -1 -1 -1 \n");
  fprintf(fp88, "   </DataArray> \n");
  fprintf(fp88, "  </CellData> \n");
  fprintf(fp88, "  <Points> \n");
  fprintf(fp88, "  <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"> \n");
  //Here the point data of the hexa has to be inserted!
  //fprintf(fp88, "   \n");
   /* create boundaries from object volume */
 // volumedimension(&X_lay[1], &X_lay[0], &ymax, &ymin, &Zmax_lay[1], &Zmin_lay[0]);
  ymin=y_offset;
  ymax=y_offset + (sqrt(X_lay[0]*X_lay[0])+sqrt(X_lay[1]*X_lay[1]))/2;
  //
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[0] + x_offset , ymin, Zmin_lay[0] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[0] + x_offset , ymax, Zmin_lay[0] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[1] + x_offset , ymin, Zmin_lay[0] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[1] + x_offset , ymax, Zmin_lay[0] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[0] + x_offset , ymin, Zmax_lay[1] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[0] + x_offset , ymax, Zmax_lay[1] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[1] + x_offset , ymin, Zmax_lay[1] + z_offset);
  fprintf(fp88, "    %5.2f %5.2f %5.2f \n", X_lay[1] + x_offset , ymax, Zmax_lay[1] + z_offset);
  //
  fprintf(fp88, " </DataArray>  \n");
  fprintf(fp88, "  </Points> \n");
  fprintf(fp88, "  <Cells> \n");
  fprintf(fp88, "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> \n");
  fprintf(fp88, "  0 4 2 6 1 5 3 7 \n");
  fprintf(fp88, "   </DataArray> \n");
  fprintf(fp88, "  <!-- 0 9 3 12 1 10 4 13 --> \n");
  fprintf(fp88, "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> \n");
  fprintf(fp88, "  8 \n");
  fprintf(fp88, "   </DataArray> \n");
  fprintf(fp88, "  <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> \n");
  fprintf(fp88, "  11 \n");
  fprintf(fp88, "  </DataArray> \n");
  fprintf(fp88, "   </Cells>\n");
  fprintf(fp88, "  </Piece> \n");
  fprintf(fp88, "  </UnstructuredGrid> \n");
  fprintf(fp88, "  </VTKFile> \n");

  		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, "parameters/velo_trans.par");
		strcat(tmp_path, tmp_file_name);
	fp77=fopen(tmp_path,"r");
	   fscanf(fp77, "%f\n", &interframeTime);//time between two frames in ms to calc velocity 
	   fscanf(fp77, "%d\n", &minSelect);
	fclose(fp77);

    /* read trackfile from ptv and create vectorfield-data */
 if(flag4export[1] == 1)
 {
  for (i=seq_first; i<seq_last;i++)
    {
      if (i < 10)             sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)       sprintf (val, "res/ptv_is.%2d",  i);
	  else if (i <1000)       sprintf (val, "res/ptv_is.%3d",  i);
      else       sprintf (val, "res/ptv_is.%4d",  i);
  
   color = ((double)(i-seq_first))/((double)(seq_last-1-seq_first));

      //printf("Create VTK-Data, read file: %s\n", val);   
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");
      
      fscanf (fp1,"%d\n", &anz1);    
      line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
      for (j=0;j<anz1;j++)
		  {
			fscanf (fp1, "%d\n", &line1[j].p);
			fscanf (fp1, "%d\n", &line1[j].n);
			fscanf (fp1, "%lf\n", &line1[j].x1);
			fscanf (fp1, "%lf\n", &line1[j].y1);
			fscanf (fp1, "%lf\n", &line1[j].z1);
		  }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)             sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)       sprintf (val, "res/ptv_is.%2d",  i+1);
	  else if (i+1 <1000)       sprintf (val, "res/ptv_is.%3d",  i+1);
      else       sprintf (val, "res/ptv_is.%4d",  i+1);
      
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");     
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
      
      for (j=0;j<anz2;j++)
		  {
			fscanf (fp1, "%d\n", &line2[j].p); 
			fscanf (fp1, "%d\n", &line2[j].n);
			fscanf (fp1, "%lf\n", &line2[j].x1);
			fscanf (fp1, "%lf\n", &line2[j].y1);
			fscanf (fp1, "%lf\n", &line2[j].z1);
		  }
      fclose (fp1);
      
	//calc number of lines:
	nOfP=0;
    for(j=0;j<anz1;j++) 
	{ 	
		m = line1[j].n;
		if (m >= 0 && m <= anz2 && line2[m].p>0)
		{
			//count connected particles
			nOfP = nOfP++;
		}
	 }
     // fprintf(fp2, "# time step %d\n", i);

	//added if-condition to ensure usage of valid data
/*	if(nOfP>1)
	{*/

	  if (i < 10)		sprintf (val2, "PostProcessing/prd_Umean_%1d.vtk", i);
      else if (i < 100)	sprintf (val2, "PostProcessing/prd_Umean_%2d.vtk", i);
	  else if (i< 1000)	sprintf (val2, "PostProcessing/prd_Umean_%3d.vtk", i);
      else				sprintf (val2, "PostProcessing/prd_Umean_%4d.vtk", i);
	  	strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val2);
		strcat(tmp_path, tmp_file_name);
	  fp66 = fopen(tmp_path,"w");

	  //fprintf(fp88, "timeStep #%d, Xm, Ym, Zm, UmeanPlus, UmeanMinus, dx, dy, dz, dl, dt \n", i); //write time step to csv
	  //write data to vtk file in simple legacy format
	  fprintf(fp66, "# vtk DataFile Version 2.0 \n");
	  fprintf(fp66, "Xm, Ym, Zm, Uvecs, timestep %d \n",i);
	  fprintf(fp66, "ASCII \n \n");
	  
	  fprintf(fp66, "DATASET POLYDATA POINTS %d float \n", nOfP );

		nOfP=0;
      for(j=0;j<anz1;j++) 
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < anz2 && line2[m].p>0) 
		{	
	

			/* create cubes of last time step */	 
			if(i ==seq_last-1) {	  }
	  
			/* create last cubes of ending trajectories  */
	  
			if(line2[m].n <= 0) {	  }

			//calc middle and translation
			mx=(line1[j].x1+line2[m].x1)/2 + x_offset;
			my=(line1[j].y1+line2[m].y1)/2 + y_offset; 
			mz=(line1[j].z1+line2[m].z1)/2 + z_offset; 
			//be careful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
			dx=line1[j].x1-line2[m].x1;
			dy=line1[j].y1-line2[m].y1;
			dz=line1[j].z1-line2[m].z1;
			du=sqrt(dx*dx+dy*dy);
			dl=sqrt(dx*dx+dy*dy+dz*dz);
			MeanU=dl/interframeTime;

			//fprintf(fp88, "%d, %7.3f, %7.3f, %7.3f, %7.6f, %7.6f, %7.3f, %7.3f, %7.3f, %7.3f,   %7.3f \n",i, mx, my, mz, MeanU, MeanU, dx, dy, dz, dl, interframeTime);
		
			fprintf(fp66, "%7.3f %7.3f %7.3f \n", mx, my, mz);
			//fprintf(fp66, "%3.6f, %3.6f, %3.5f, %3.5f, %3.5f, %3.5f, %5.3f \n", MeanU, MeanU, dx, dy, dz, dl, interframeTime);

			rotz=0;
			if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
			rotx=0;
			if(du == 0.0) {rotx=M_PI/2;}

			if(du != 0.0) 
			{
				if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
				if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
			}
		}	
      }
  
  
	  fclose(fp66);
      strcpy(val, "");
      free(line1); free(line2);
/*	 } //end if nOfP>1
	else
	{
	  printf("Skipped time-step! Too few points in data.\n");  
	}*/
   }  /* end of first sequence loop */
  
  
 // fclose(fp2);
  fclose(fp88);
  
  printf("Inter fame time: %f ms\n", interframeTime);
  for (i=seq_first; i<seq_last;i++)
    {
      if (i < 10)           sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)     sprintf (val, "res/ptv_is.%2d", i);
	  else if (i < 1000)	sprintf (val, "res/ptv_is.%3d", i);
      else					sprintf (val, "res/ptv_is.%4d", i);
  
      color = ((double)(i-seq_first))/((double)(seq_last-1-seq_first));

      printf("Create VTK-Data, read file: %s\n", val);   
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");
      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) malloc((anz1+1)*sizeof (vector));
      for (j=0;j<anz1;j++) 
		  {
			fscanf (fp1, "%d\n", &line1[j].p);
			fscanf (fp1, "%d\n", &line1[j].n);
			fscanf (fp1, "%lf\n", &line1[j].x1);
			fscanf (fp1, "%lf\n", &line1[j].y1);
			fscanf (fp1, "%lf\n", &line1[j].z1);
		  }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)         sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)   sprintf (val, "res/ptv_is.%2d", i+1);
	  else if (i+1 < 1000)  sprintf (val, "res/ptv_is.%3d", i+1);
      else					sprintf (val, "res/ptv_is.%4d", i+1);
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");          
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) malloc((anz2+1)*sizeof (vector));
      
      for (j=0;j<anz2;j++) 
		  {
				fscanf (fp1, "%d\n", &line2[j].p);
				fscanf (fp1, "%d\n", &line2[j].n);
				fscanf (fp1, "%lf\n", &line2[j].x1);
				fscanf (fp1, "%lf\n", &line2[j].y1);
				fscanf (fp1, "%lf\n", &line2[j].z1);
		  }
      fclose (fp1);
      
    //  fprintf(fp2, "# time step %d\n", i);

	  if (i < 10)           sprintf (val2, "PostProcessing/prd_Umean_%1d.vtk", i);
      else if (i < 100)     sprintf (val2, "PostProcessing/prd_Umean_%2d.vtk", i);
	  else if (i < 1000)    sprintf (val2, "PostProcessing/prd_Umean_%3d.vtk", i);
      else					sprintf (val2, "PostProcessing/prd_Umean_%4d.vtk", i);

// add information of vertices
	  	strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
	   fp66 = fopen(tmp_path,"a");
	  
	  nOfP=0;
	 
      for(j=0;j<anz1;j++) 
	  { 	
			m = line1[j].n;
			if (m >= 0 && m < anz2 && line2[m].p>0)
			{
				//count connected particles
				nOfP = nOfP++;
			}
	  }
	   nVerts=0;
	  nVerts=2*nOfP;
	 fprintf(fp66, "\n VERTICES %d %d \n",nOfP, nVerts);
	fclose(fp66);  
	
	fp66 = fopen(tmp_path,"a");
	  
	  nOfP=0;
      for(j=0;j<anz1;j++) 
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < anz2 && line2[m].p>0)
			{
			  //count connected particles
				fprintf(fp66,"1 %d \n", nOfP);
				nOfP = nOfP++;
			}

	  }

	fclose(fp66);

//append velocity-data part 1
	  fp66 = fopen(tmp_path,"a");

	  nOfP=0;
      for(j=0;j<anz1;j++) 
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < anz2 && line2[m].p>0)
			{
			  //count connected particles
			  nOfP = nOfP++;
			}
	  }

	 fprintf(fp66, "POINT_DATA %d \n", nOfP);
	 fprintf(fp66, "SCALARS MeanU double 1 \n");
	 fprintf(fp66, "LOOKUP_TABLE default\n");
	
 
	nOfP=0;

      for(j=0;j<anz1;j++)
	  { 	
			m = line1[j].n;
			if (m >= 0 && m < anz2 && line2[m].p>0) 
			{		 /* create cubes of last time step */	 
					 if(i ==seq_last-1) {   }
	  				 /* create last cubes of ending trajectories  */
					 if(line2[m].n <= 0) {	 }

					 mx=(line1[j].x1+line2[m].x1)/2;
					 my=(line1[j].y1+line2[m].y1)/2; 
					 mz=(line1[j].z1+line2[m].z1)/2; 
					 //be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
					 dx=line1[j].x1-line2[m].x1;
					 dy=line1[j].y1-line2[m].y1;
					 dz=line1[j].z1-line2[m].z1;
					 du=sqrt(dx*dx+dy*dy);
					 dl=sqrt(dx*dx+dy*dy+dz*dz);
					 MeanU=dl/interframeTime;
					 //print mean U with check of max Velo
					 if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))fprintf(fp66, "-1.00000\n");
					 else fprintf(fp66, "%3.5f \n", MeanU);

					rotz=0;
					if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
					rotx=0;
					if(du == 0.0) {rotx=M_PI/2;}

					if(du != 0.0) 
						{
						  if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
						  if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
						  if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
						  if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
						  if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
						  if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
						  if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
						  if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
						}
	    	  	   
			}	
      }
  
      //fprintf(fp2, "# end of time step %d\n\n", i); 
	  fclose(fp66);
//append velocity data 2 (vecs)
	   fp66 = fopen(tmp_path,"a");

	  nOfP=0;
      for(j=0;j<anz1;j++) 
		  { 	
			m = line1[j].n;
			if (m >= 0 && m < anz2)
				{
				  //count connected particles
				  nOfP = nOfP++;
				}
		  }
	 fprintf(fp66,"\n VECTORS Umean float \n");

	nOfP=0;
	for(j=0;j<anz1;j++)
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < anz2 && line2[m].p>0) 
		{	
	 

			 /* create cubes of last time step */	 
			 if(i ==seq_last-1) { }
	  		 /* create last cubes of ending trajectories  */
	  		 if(line2[m].n <= 0) { }

			mx=(line1[j].x1+line2[m].x1)/2;
			my=(line1[j].y1+line2[m].y1)/2; 
			mz=(line1[j].z1+line2[m].z1)/2; 
			dx=line1[j].x1-line2[m].x1;
			dy=line1[j].y1-line2[m].y1;
			dz=line1[j].z1-line2[m].z1;
			//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
			du=sqrt(dx*dx+dy*dy);
			dl=sqrt(dx*dx+dy*dy+dz*dz);
			MeanU=dl/interframeTime;
			ux= -dx/interframeTime;
			uy= -dy/interframeTime;
			uz= -dz/interframeTime;
			//print velos
			fprintf(fp66, "%3.5f %3.5f %3.5f \n", ux, uy,uz);

			rotz=0;
			if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
			rotx=0;
			if(du == 0.0) {rotx=M_PI/2;}

			if(du != 0.0) 
				{
				  if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				  if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				  if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				  if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				  if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				  if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
				  if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				  if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
				}
	    	  	   
		}	
	}
  
fclose(fp66);
//append particle-data part 1
	  fp66 = fopen(tmp_path,"a");

	  nOfP=0;
      for(j=0;j<anz1;j++) 
		  { 	
			m = line1[j].n;
			if (m >= 0 && m < anz2 && line2[m].p>0)
			{
			  //count connected particles
			  nOfP = nOfP++;
			}
		  }

	 //fprintf(fp66, "\n POINT_DATA %d \n", nOfP);
	 fprintf(fp66, "SCALARS NextParticleID double 1 \n");
	 fprintf(fp66, "LOOKUP_TABLE default\n");
	
 
		nOfP=0;

      for(j=0;j<anz1;j++)
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < anz2 && line2[m].p>0) 
		{	
			fprintf(fp66, "%d \n", m);
			/* create cubes of last time step */	 
			if(i ==seq_last-1) {   }
			/* create last cubes of ending trajectories  */
			if(line2[m].n <= 0) {	 }

			mx=(line1[j].x1+line2[m].x1)/2;
			my=(line1[j].y1+line2[m].y1)/2; 
			mz=(line1[j].z1+line2[m].z1)/2; 
			//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
			dx=line1[j].x1-line2[m].x1;
			dy=line1[j].y1-line2[m].y1;
			dz=line1[j].z1-line2[m].z1;
			du=sqrt(dx*dx+dy*dy);
			dl=sqrt(dx*dx+dy*dy+dz*dz);
			MeanU=dl/interframeTime;
		
			rotz=0;
			if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
			rotx=0;
			if(du == 0.0) {rotx=M_PI/2;}

			if(du != 0.0) 
				{
				  if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				  if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				  if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				  if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				  if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				  if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
				  if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				  if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
				}
		} 	  	   
	}	

fclose(fp66);

// This part does not work fine in paraview,
//because the number of previous particles does not
//necessarily match with the actual number nOfP

//Code to append particle-data part 2
	  fp66 = fopen(tmp_path,"a");

/*	  nOfP=0;
      for(j=0;j<anz1;j++) { 	
	m = line1[j].p;
	if (m >= 0)
	{
	  //count connected particles
	  nOfP = nOfP++;
	}
	  }
*/
		fprintf(fp66, "SCALARS PreviousParticleID double 1 \n");
	 fprintf(fp66, "LOOKUP_TABLE default\n");
	
 
		nOfP=0;

      for(j=0;j<anz1;j++)
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < anz2 && line2[m].p>0) 
		{	
			 fprintf(fp66, "%d \n", line1[j].p);
			 // create cubes of last time step	 
			 if(i ==seq_last-1) {   }
	  		 // create last cubes of ending trajectories
			 if(line2[m].n <= 0) {	 }

			 mx=(line1[j].x1+line2[m].x1)/2;
			 my=(line1[j].y1+line2[m].y1)/2; 
			 mz=(line1[j].z1+line2[m].z1)/2; 
			 //be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
			 dx=line1[j].x1-line2[m].x1;
			 dy=line1[j].y1-line2[m].y1;
			 dz=line1[j].z1-line2[m].z1;
			 du=sqrt(dx*dx+dy*dy);
			 dl=sqrt(dx*dx+dy*dy+dz*dz);
			 MeanU=dl/interframeTime;
		
			rotz=0;
			if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
			rotx=0;
			if(du == 0.0) {rotx=M_PI/2;}

			if(du != 0.0) 
			{
			  if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
			  if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
			  if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
			  if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
			  if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
			  if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
			  if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
			  if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
			}
	    	  	   
		}
      }

	  //Start here, if you like to add more specific pointdata
  
	  fclose(fp66);

      strcpy(val, "");
      free(line1); free(line2);
    }  /* end of sequence loop */
}  

  //fprintf(fp2, "# trajectories finished\n");
  
  //fclose(fp2);
  //fclose(fp88);

 

  //Tcl_Eval(interp, ".text delete 2");
  //Tcl_Eval(interp, ".text insert 2 \"Points & Velocities written to vtk-Files\"");
  //Tcl_Eval(interp, "update idletasks");  

  //Tcl_Eval(interp, ".text delete 3");
  //Tcl_Eval(interp, ".text insert 3\"Generate Trajectories\"");
  //Tcl_Eval(interp, "update idletasks"); 

  printf("Points & Velocities written to vtk-Files.\n");
  printf("Generate Trajectories.\n");
//**********************************************************************//
  //Trial to export time resolved trajectories into paraview data format // P.Steinhoff, 20.04.2012
if(flag4export[2]==1)
{
  //read interframe-time
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, "parameters/velo_trans.par");
		strcat(tmp_path, tmp_file_name);
   fp77=fopen(tmp_path,"r");
	  fscanf(fp77, "%f\n", &interframeTime);
	  fscanf(fp77, "%d\n", &minSelect);
   fclose(fp77);

		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, "PostProcessing/Trajects/");
		strcat(tmp_path, tmp_file_name);
    CreateDirectoryA(tmp_path,NULL);
    /* read trackfile from ptv and create vectorfield-data */
 
  for (i=seq_first; i<seq_last;i++)
    {
      if (i < 10)           sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)     sprintf (val, "res/ptv_is.%2d", i);
      else if (i < 1000)    sprintf (val, "res/ptv_is.%3d", i);
	  else					sprintf (val, "res/ptv_is.%4d", i);
  
      //printf("Create VTK-Data, read file: %s\n", val); 
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");
      
      color = ((double)(i-seq_first))/((double)(seq_last-1-seq_first));

      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
      for (j=0;j<anz1;j++) 
	  {
		fscanf (fp1, "%d\n", &line1[j].p);
		fscanf (fp1, "%d\n", &line1[j].n);
		fscanf (fp1, "%lf\n", &line1[j].x1);
		fscanf (fp1, "%lf\n", &line1[j].y1);
		fscanf (fp1, "%lf\n", &line1[j].z1);
      }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)         sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)   sprintf (val, "res/ptv_is.%2d", i+1);
      else if (i+1 < 1000)  sprintf (val, "res/ptv_is.%3d", i+1);
	  else					sprintf (val, "res/ptv_is.%4d", i+1);
      
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");     
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
      
      for (j=0;j<anz2;j++) 
	  {
		fscanf (fp1, "%d\n", &line2[j].p); 
		fscanf (fp1, "%d\n", &line2[j].n);
		fscanf (fp1, "%lf\n", &line2[j].x1);
		fscanf (fp1, "%lf\n", &line2[j].y1);
		fscanf (fp1, "%lf\n", &line2[j].z1);
      }
      fclose (fp1);
      
     // fprintf(fp2, "# time step %d\n", i);

	  if (i < 10)			sprintf (val3, "PostProcessing/Trajects/Traject_%1d.vtk",i);
      else if (i < 100)     sprintf (val3, "PostProcessing/Trajects/Traject_%2d.vtk",i);
	  else if (i < 1000)    sprintf (val3, "PostProcessing/Trajects/Traject_%3d.vtk",i);
      else					sprintf (val3, "PostProcessing/Trajects/Traject_%4d.vtk",i);

		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val3);
		strcat(tmp_path, tmp_file_name);
	  fp99 = fopen(tmp_path,"w");
	  if(!fp99)printf("Could not create file %s.\n",val3);
	  //fprintf(fp88, "timeStep #%d, Xm, Ym, Zm, UmeanPlus, UmeanMinus, dx, dy, dz, dl, dt \n", i); //write time step to csv
	  //write data to vtk file in simple legacy format
	  fprintf(fp99, "# vtk DataFile Version 2.1 \n");
	  fprintf(fp99, "Trajects from timestep %d with color option \n",i);
	  fprintf(fp99, "ASCII \n \n");
	  
	 //calc number of lines:
	  nOfP=0;
      for(j=0;j<anz1;j++) 
	  { 	
		m = line1[j].n;
		if (m >= 0 && m < Ms)
		{
		  //count connected particles
		  nOfP = nOfP++;
		}
	  }

 fprintf(fp99, "DATASET POLYDATA POINTS %d float \n", 2*nOfP );

 veloPos1 = (vector *) malloc((anz1+1)*sizeof(vector)); //Changed allocated memory to 500*xyz due to problems in allocating enough memory and writing data
 if(!veloPos1)printf("Error malloc veloPos1.\n");
		
 nOfP=0;
 for(j=0;j<anz1;j++)
 { 	
	m = line1[j].n;
	if (m >= 0 && m < anz2) 
	{	
			/* create cubes of last time step */	 
			if(i ==seq_last-1) {	  }
	  
			/* create last cubes of ending trajectories  */
			if(line2[m].n <= 0) {	  }

			fprintf(fp99, "%5.3f %5.3f %5.3f \n", line1[j].x1 + x_offset, line1[j].y1 + y_offset, line1[j].z1 + z_offset);
			fprintf(fp99, "%5.3f %5.3f %5.3f \n", line2[m].x1 + x_offset, line2[m].y1 + y_offset, line2[m].z1 + z_offset);
			//What happens here is not necessary for trajects...
			mx=(line1[j].x1+line2[m].x1)/2;
			my=(line1[j].y1+line2[m].y1)/2; 
			mz=(line1[j].z1+line2[m].z1)/2; 
			//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
			dx=line1[j].x1-line2[m].x1;
			dy=line1[j].y1-line2[m].y1;
			dz=line1[j].z1-line2[m].z1;
			du=sqrt(dx*dx+dy*dy);
			dl=sqrt(dx*dx+dy*dy+dz*dz);
			MeanU=dl/interframeTime;
			veloPos1[j].Ux=(-dx/interframeTime);
			veloPos1[j].Uy=(-dy/interframeTime);
			veloPos1[j].Uz=(-dz/interframeTime);
			//print mean U with check of max Velo
			if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))veloPos1[j].Um=-1;
		    //if (MeanU > 2*tpar.Vmax)veloPos1[j].Um=-1;
		    else veloPos1[j].Um=MeanU;
				

			rotz=0;
			if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
			rotx=0;
			if(du == 0.0) {rotx=M_PI/2;}

			if(du != 0.0) 
			 {
				 if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				 if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				 if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
				 if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
				 if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				 if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
				 if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
				 if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
			 }
	}	
 }
	   fclose(fp99);
  
// add information of trajects for t and t+1
	fp55 = fopen(tmp_path,"a");
	nOfP=0; 
    for(j=0;j<anz1;j++) 
	{ 	
		m = line1[j].n;
		if (m >= 0 && m < Ms)
		{
		  //count connected particles
			nOfP = nOfP++;
		}
	}
	nVerts=0;
	nVerts=3*nOfP;
	fprintf(fp55, "\n LINES %d %d \n",nOfP, nVerts);
	fclose(fp55);  
	
	fp55 = fopen(tmp_path,"a"); 
	nOfP=0;
    for(j=0;j<anz1;j++) 
	{ 	
		m = line1[j].n;
		if (m >= 0 && m < Ms)
		{
		  //count connected particles
			fprintf(fp55,"2 %d %d \n", nOfP, nOfP+1);
			nOfP = nOfP+2;
		}

	}

	fclose(fp55);

// add/define time-step specific color for line
	   fp55 = fopen(tmp_path,"a");
	  
	  nOfP=0;
	 
      for(j=0;j<anz1;j++) { 	
	m = line1[j].n;
	if (m >= 0 && m < Ms)
	{
	  //count connected particles
		nOfP = nOfP++;
	}
	
	  }
	  
	  fprintf(fp55, "POINT_DATA %d \n", 2*nOfP);
	  
	 //fprintf(fp55, "\n LINES %d %d \n",nOfP, nVerts);
	fclose(fp55);  
	/*
	fp55 = fopen(val3,"a");
	  fprintf(fp55, "SCALARS colorID double 1 \n");
	  fprintf(fp55, "LOOKUP_TABLE default\n");
	  nOfP=0;
      for(j=0;j<anz1;j++) { 	
	m = line1[j].n;
	if (m >= 0)
	{
		color = ((double)(i-seq_first))/((double)(seq_last-1-seq_first));
		fprintf(fp55,"%d \n", color);
		fprintf(fp55,"%d \n", color);
		nOfP = nOfP++;
	}

	  }

	fclose(fp55);
	*/
	//append mean veloData to "lines"
	
	fp55 = fopen(tmp_path,"a");
	fprintf(fp55, "SCALARS Umean double 1 \n");
	fprintf(fp55, "LOOKUP_TABLE default\n");

      for(j=0;j<anz1;j++) 
		{ 	
			 m = line1[j].n;
			 if (m >= 0 && m < Ms) 
			 {
				fprintf(fp55,"%2.5f \n", veloPos1[j].Um);
				fprintf(fp55,"%2.5f \n", veloPos1[j].Um);
			 }
		}
    fclose(fp55);

	//append mean x-veloData to "lines"
	
	fp55 = fopen(tmp_path,"a");
	fprintf(fp55, "SCALARS Umean_X double 1 \n");
	fprintf(fp55, "LOOKUP_TABLE default\n");
    for(j=0;j<anz1;j++) 
	{ 	
		 m = line1[j].n;
		 if (m >= 0 && m < Ms)
		 {
			 fprintf(fp55,"%2.5f \n", veloPos1[j].Ux);
			 fprintf(fp55,"%2.5f \n", veloPos1[j].Ux);
		 }
    }
    fclose(fp55);

	fp55 = fopen(tmp_path,"a");
	fprintf(fp55, "SCALARS Umean_Y double 1 \n");
	fprintf(fp55, "LOOKUP_TABLE default\n");
    for(j=0;j<anz1;j++) 
	{ 	
		m = line1[j].n;
		if (m >= 0 && m < Ms)
		{
			fprintf(fp55,"%2.5f \n", veloPos1[j].Uy);
			fprintf(fp55,"%2.5f \n", veloPos1[j].Uy);
		}
    }
    fclose(fp55);

	fp55 = fopen(tmp_path,"a");
	fprintf(fp55, "SCALARS Umean_Z double 1 \n");
	fprintf(fp55, "LOOKUP_TABLE default\n");
    for(j=0;j<anz1;j++) 
	{ 	
		m = line1[j].n;
		if (m >= 0 && m < Ms)
		{
			fprintf(fp55,"%2.5f \n", veloPos1[j].Uz);
			fprintf(fp55,"%2.5f \n", veloPos1[j].Uz);
		}
	}
    fclose(fp55);
	 
    strcpy(val, "");
    free(line1); 
	free(line2);
	free(veloPos1);
    }  /* end of sequence loop */

  //end of trial

//*******************************************************************************//
 //Trial to export trajectories sorted with/by ID into paraview data format // P.Steinhoff, 26.04.2012

  //read interframe-time 
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, "parameters/velo_trans.par");
		strcat(tmp_path, tmp_file_name);
   fp77=fopen(tmp_path,"r");
	   fscanf(fp77, "%f\n", &interframeTime);
	   fscanf(fp77, "%d\n", &minSelect);
   fclose(fp77);

    /* read trackfile from ptv and create traject-data */
  for (i=seq_first; i<seq_last;i++)
    {
	  q=0;
      if (i < 10)           sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)     sprintf (val, "res/ptv_is.%2d", i);
	  else if (i < 1000)    sprintf (val, "res/ptv_is.%3d", i);
      else					sprintf (val, "res/ptv_is.%4d", i);

		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      printf("Create VTK-Data for trajectories, read file: %s\n", val);     
      fp1 = fopen (tmp_path, "r");
      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
      for (j=0;j<anz1;j++) 
	  {
		fscanf (fp1, "%d\n", &line1[j].p);
		fscanf (fp1, "%d\n", &line1[j].n);
		fscanf (fp1, "%lf\n", &line1[j].x1);
		fscanf (fp1, "%lf\n", &line1[j].y1);
		fscanf (fp1, "%lf\n", &line1[j].z1);
      }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)         sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)   sprintf (val, "res/ptv_is.%2d", i+1);
      else if (i+1 < 1000)  sprintf (val, "res/ptv_is.%3d", i+1);
	  else					sprintf (val, "res/ptv_is.%4d", i+1);
      	strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val);
		strcat(tmp_path, tmp_file_name);
      fp1 = fopen (tmp_path, "r");     
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
      
      for (j=0;j<anz2;j++) 
		  {
			fscanf (fp1, "%d\n", &line2[j].p); 
			fscanf (fp1, "%d\n", &line2[j].n);
			fscanf (fp1, "%lf\n", &line2[j].x1);
			fscanf (fp1, "%lf\n", &line2[j].y1);
			fscanf (fp1, "%lf\n", &line2[j].z1);
		  }
      fclose (fp1);
      
//start for-loop to write data in file

	for(l=0;l<anz1;l++) 
	{
		m = line1[l].n;

		if(m >= 0 && m <= anz2 && line1[l].p < 0) //to check if a new traject beginns
		{
		  tmpVeloPos=1;	  	
		  //count points of traject
		  nPoints=1;
		  q=q+1; //every traject becomes an ID
		  lineOfP=line1[l].n;
	  
		  //for(timeID=i+1;timeID <= seq_last && lineOfP >= 0 && lineOfP < Ms;timeID++)
		  for(timeID=i+1;timeID<=seq_last && lineOfP>=0 && lineOfP<= Ms;timeID++)
		  {
	 			time=timeID;

				if (time < 10)          sprintf (val, "res/ptv_is.%1d", time);
				else if (time < 100)    sprintf (val, "res/ptv_is.%2d", time);
				else if (time < 1000)   sprintf (val, "res/ptv_is.%3d", time);
				else					sprintf (val, "res/ptv_is.%4d", time);
		
				strcpy(tmp_path, base_path);
				strcpy(tmp_file_name, val);
				strcat(tmp_path, tmp_file_name);
				fp1 = fopen (tmp_path, "r");
				fscanf (fp1,"%d\n", &anz3);
      
				line3 = (vector *) malloc ((anz3+1)*sizeof (vector));
				for (j=0;j<anz3;j++)
					{
						fscanf (fp1, "%d\n", &line3[j].p);
						fscanf (fp1, "%d\n", &line3[j].n);
						fscanf (fp1, "%lf\n", &line3[j].x1);
						fscanf (fp1, "%lf\n", &line3[j].y1);
						fscanf (fp1, "%lf\n", &line3[j].z1);
					}

				strcpy(val, "");     
				fclose (fp1);

				if(line3[lineOfP].n <0 && line3[lineOfP].p>=0)//handle last Point of traject
				{
					// fprintf(fp111, "%5.3f %5.3f %5.3f \n", line3[lineOfP].x1, line3[lineOfP].y1, line3[lineOfP].z1);
					//q=q+1; 
				}

				if(lineOfP >= 0 && lineOfP < Ms) 
				{
					lineOfP=line3[lineOfP].n;		//update line of next particle-position
					nPoints=nPoints+1;				//for each valid link set counter for length of traject +1
					//add check if line3[lineOfP].n points at a real link--> line4[line3[lineOfP].n].p==lineOfP
				}
		  
			  free(line3);//clear vector
		   } //end of timeID-loop

	  if(nPoints>=minSelect)
	  {
 
		//ensure trajects to be sorted by length 
		//sprintf (selectionPath, "PostProcessing/Trajects/Selection%3d/",nPoints); 
		//CreateDirectoryA(selectionPath,NULL);
 
	   if(nPoints<10)
	   {
			//ensure trajects to be sorted by length 
			sprintf (selectionPath, "PostProcessing/Trajects/Selection%1d/",nPoints); 
			CreateDirectoryA(selectionPath,NULL);
			if(i<10)
			{
							if (q < 10)				sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%1d_length_%1d_ID_%1d.vtk",nPoints,i,nPoints,q);
							else if (q < 100)		sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%1d_length_%1d_ID_%2d.vtk",nPoints,i,nPoints,q);
							else if (q < 1000)		sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%1d_length_%1d_ID_%3d.vtk",nPoints,i,nPoints,q);
							else					sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%1d_length_%1d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
			else if(i<100)
			{
						if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%2d_length_%1d_ID_%1d.vtk",nPoints,i,nPoints,q);
						else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%2d_length_%1d_ID_%2d.vtk",nPoints,i,nPoints,q);
						else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%2d_length_%1d_ID_%3d.vtk",nPoints,i,nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%2d_length_%1d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
			else if(i<1000)
			{
						if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%3d_length_%1d_ID_%1d.vtk",nPoints,i,nPoints,q);
						else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%3d_length_%1d_ID_%2d.vtk",nPoints,i,nPoints,q);
						else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%3d_length_%1d_ID_%3d.vtk",nPoints,i,nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%3d_length_%1d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
			else	
			{
						if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%4d_length_%1d_ID_%1d.vtk",nPoints,i,nPoints,q);
						else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%4d_length_%1d_ID_%2d.vtk",nPoints,i,nPoints,q);
						else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%4d_length_%1d_ID_%3d.vtk",nPoints,i,nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Selection%1d/Traject_time_%4d_length_%1d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
	   }
	   else //more than 99 points are not supported yet (P.Steinhoff)
	   {
			//ensure trajects to be sorted by length 
			sprintf (selectionPath, "PostProcessing/Trajects/Selection%2d/",nPoints); 
			CreateDirectoryA(selectionPath,NULL);
			if(i<10)
			{
							if (q < 10)				sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%1d_length_%2d_ID_%1d.vtk",nPoints,i,nPoints,q);
							else if (q < 100)		sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%1d_length_%2d_ID_%2d.vtk",nPoints,i,nPoints,q);
							else if (q < 1000)		sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%1d_length_%2d_ID_%3d.vtk",nPoints,i,nPoints,q);
							else					sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%1d_length_%2d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
			else if(i<100)
			{
						if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%2d_length_%2d_ID_%1d.vtk",nPoints,i,nPoints,q);
						else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%2d_length_%2d_ID_%2d.vtk",nPoints,i,nPoints,q);
						else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%2d_length_%2d_ID_%3d.vtk",nPoints,i,nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%2d_length_%2d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
			else if(i<1000)
			{
						if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%3d_length_%2d_ID_%1d.vtk",nPoints,i,nPoints,q);
						else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%3d_length_%2d_ID_%2d.vtk",nPoints,i,nPoints,q);
						else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%3d_length_%2d_ID_%3d.vtk",nPoints,i,nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%3d_length_%2d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
			else
			{
						if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%4d_length_%2d_ID_%1d.vtk",nPoints,i,nPoints,q);
						else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%4d_length_%2d_ID_%2d.vtk",nPoints,i,nPoints,q);
						else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%4d_length_%2d_ID_%3d.vtk",nPoints,i,nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Selection%2d/Traject_time_%4d_length_%2d_ID_%4d.vtk",nPoints,i,nPoints,q);
			}
	   }
	  
	}//end if nPoints>minSelect
	else //write to alternative directory
	{
		if(i<10)
		{
						if (q < 10)				sprintf (val4, "PostProcessing/Trajects/Traject_time_%1d/Traject_length_%5d_ID_%1d.vtk",i,10000+nPoints,q);
						else if (q < 100)		sprintf (val4, "PostProcessing/Trajects/Traject_time_%1d/Traject_length_%5d_ID_%2d.vtk",i,10000+nPoints,q);
						else if (q < 1000)		sprintf (val4, "PostProcessing/Trajects/Traject_time_%1d/Traject_length_%5d_ID_%3d.vtk",i,10000+nPoints,q);
						else					sprintf (val4, "PostProcessing/Trajects/Traject_time_%1d/Traject_length_%5d_ID_%4d.vtk",i,10000+nPoints,q);
		}
		else if(i<100)
		{
					if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Traject_time_%2d/Traject_length_%5d_ID_%1d.vtk",i,10000+nPoints,q);
					else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Traject_time_%2d/Traject_length_%5d_ID_%2d.vtk",i,10000+nPoints,q);
					else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Traject_time_%2d/Traject_length_%5d_ID_%3d.vtk",i,10000+nPoints,q);
					else					sprintf (val4, "PostProcessing/Trajects/Traject_time_%2d/Traject_length_%5d_ID_%4d.vtk",i,10000+nPoints,q);
		}
		else if(i<1000)
		{
					if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Traject_time_%3d/Traject_length_%5d_ID_%1d.vtk",i,10000+nPoints,q);
					else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Traject_time_%3d/Traject_length_%5d_ID_%2d.vtk",i,10000+nPoints,q);
					else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Traject_time_%3d/Traject_length_%5d_ID_%3d.vtk",i,10000+nPoints,q);
					else					sprintf (val4, "PostProcessing/Trajects/Traject_time_%3d/Traject_length_%5d_ID_%4d.vtk",i,10000+nPoints,q);
		}
		else
		{
					if (q < 10)             sprintf (val4, "PostProcessing/Trajects/Traject_time_%4d/Traject_length_%5d_ID_%1d.vtk",i,10000+nPoints,q);
					else if (q < 100)       sprintf (val4, "PostProcessing/Trajects/Traject_time_%4d/Traject_length_%5d_ID_%2d.vtk",i,10000+nPoints,q);
					else if (q < 1000)      sprintf (val4, "PostProcessing/Trajects/Traject_time_%4d/Traject_length_%5d_ID_%3d.vtk",i,10000+nPoints,q);
					else					sprintf (val4, "PostProcessing/Trajects/Traject_time_%4d/Traject_length_%5d_ID_%4d.vtk",i,10000+nPoints,q);
		}
	  }//end else of directory declaration


	  //Create time directories:
		if(i<10)		sprintf (path, "PostProcessing/Trajects/Traject_time_%1d/",i);
		else if(i<100)	sprintf (path, "PostProcessing/Trajects/Traject_time_%2d/",i);
		else if(i<1000) sprintf (path, "PostProcessing/Trajects/Traject_time_%3d/",i);
		else			sprintf (path, "PostProcessing/Trajects/Traject_time_%4d/",i);
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, path);
		strcat(tmp_path, tmp_file_name);
		CreateDirectoryA(tmp_path,NULL); 
	
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, val4);
		strcat(tmp_path, tmp_file_name);
	  fp111= fopen(tmp_path,"w");
	  if(!fp111)printf("Error creating file: %s\n",tmp_path);
	  //write data to vtk file in simple legacy format
	  fprintf(fp111, "# vtk DataFile Version 2.1 \n");
	  fprintf(fp111, "Trajects from timestep %d  \n",i);
	  fprintf(fp111, "ASCII \n \n");
	  fclose(fp111);

	  //Header for number of points
	  fp111= fopen(tmp_path,"a");
	  fprintf(fp111, "DATASET POLYDATA POINTS %d float\n", nPoints);
	  //fclose(fp111);

	  //Starting-point of traject
	  //fp111= fopen(val4,"a");
	  fprintf(fp111, "%5.3f %5.3f %5.3f \n", line1[l].x1 + x_offset, line1[l].y1 + y_offset, line1[l].z1 + z_offset);
	  //fclose(fp111);

	  //Create a Matrix to calculate Velocities
	  //Write first line (inital position)
	  if (nOfP == 0) nOfP=1;
	  veloPos1 = (vector *) malloc ((nPoints +2)*sizeof (vector)); //Changed allocated memory to 500*xyz due to problems in allocating enough memory and writing data
	  //veloPos1 = (vector *) calloc ( nPoints , sizeof (vector)+1);
	  veloPos1[0].x1=line1[l].x1;
	  veloPos1[0].y1=line1[l].y1;
	  veloPos1[0].z1=line1[l].z1;

	  lineOfP=line1[l].n;

	  //fp111= fopen(val4,"a");
	 
		for(timeID=i+1;timeID<=seq_last && lineOfP >= 0 && lineOfP <= Ms;timeID++)
			{ 
				//get_next_Particle_Position_function of (lineOfP,timeID);
				time=timeID;

				if (time < 10)          sprintf (val, "res/ptv_is.%1d", time);
				else if (time < 100)    sprintf (val, "res/ptv_is.%2d", time);
				else if (time < 1000)   sprintf (val, "res/ptv_is.%3d", time);
				else					sprintf (val, "res/ptv_is.%4d", time);
				anz3=0;
				strcpy(tmp_path, base_path);
				strcpy(tmp_file_name, val);
				strcat(tmp_path, tmp_file_name);
				fp1 = fopen (tmp_path, "r");
				fscanf (fp1,"%d\n", &anz3);
      
				line3 = (vector *) malloc ((anz3+1)* sizeof (vector));
				for (j=0;j<anz3;j++) 
					{
						fscanf (fp1, "%d\n", &line3[j].p);
						fscanf (fp1, "%d\n", &line3[j].n);
						fscanf (fp1, "%lf\n", &line3[j].x1);
						fscanf (fp1, "%lf\n", &line3[j].y1);
						fscanf (fp1, "%lf\n", &line3[j].z1);
					}

				strcpy(val, "");     
				fclose (fp1);

				//append point data of one single traject
				fprintf(fp111, "%5.3f %5.3f %5.3f \n", line3[lineOfP].x1 + x_offset, line3[lineOfP].y1 + y_offset, line3[lineOfP].z1 + z_offset);
		
				//add position-data to vector/matrix
				veloPos1[tmpVeloPos].x1=line3[lineOfP].x1;
				veloPos1[tmpVeloPos].y1=line3[lineOfP].y1;
				veloPos1[tmpVeloPos].z1=line3[lineOfP].z1; 
	
				//handle last Point of traject
				 if(line3[lineOfP].n <0 && line3[lineOfP].p>=0)
					 {
						// fprintf(fp111, "%5.3f %5.3f %5.3f \n", line3[lineOfP].x1, line3[lineOfP].y1, line3[lineOfP].z1);
						//q=q+1; 
					 }
				 //update line of next particle-position
				 lineOfP=line3[lineOfP].n;
				 tmpVeloPos=tmpVeloPos++;
				 free(line3);
			}

		//fclose(fp111);

		  //add lines
		//fp111= fopen(val4,"a");
		  fprintf(fp111, "LINES %d %d\n", nPoints-1, 3*(nPoints-1));
	  
		  for(nOfP=1; nOfP < nPoints; nOfP++)
		  {
			  fprintf(fp111,"2 %d %d \n", nOfP-1, nOfP);
		  }

		//fclose(fp111);

		  //add particle/path data
		//fp111= fopen(val4,"a");
		  fprintf(fp111, "POINT_DATA %d \n", nPoints);

		  fprintf(fp111, "SCALARS Umean_xyz double 1\n");
		  fprintf(fp111, "LOOKUP_TABLE default\n");

		  for(nOfP=1; nOfP <=nPoints; nOfP++)
		  {
			  if(nOfP<nPoints)
			  {
				  dx=veloPos1[nOfP].x1-veloPos1[nOfP-1].x1;
				  dy=veloPos1[nOfP].y1-veloPos1[nOfP-1].y1;
				  dz=veloPos1[nOfP].z1-veloPos1[nOfP-1].z1;
				  MeanU = sqrt(dx*dx+dy*dy+dz*dz)/interframeTime;
				  if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))fprintf(fp111, "-1.00000\n");
				  else fprintf(fp111,"%2.5f \n", MeanU);
			  }
			  else
			  {
				  fprintf(fp111,"-0.00100 \n");
			  }
		  }
		//fclose(fp111);
	 

		  // add mean Ux for lines here...
		//fp111= fopen(val4,"a");
		  // Only one declaration of point-data-size and type
		  //fprintf(fp111, "\n POINT_DATA %d \n", nPoints);
		  fprintf(fp111, "SCALARS Umeanx double 1\n");
		  fprintf(fp111, "LOOKUP_TABLE default\n");

		  for(nOfP=1; nOfP <=nPoints; nOfP++)
		  {
			  if(nOfP<nPoints)
			  {
					fprintf(fp111,"%2.5f \n", (veloPos1[nOfP].x1-veloPos1[nOfP-1].x1)/interframeTime);
			  }
			  else
			  {
					fprintf(fp111,"-2.00000 \n");
			  }
		  }
		//fclose(fp111);

		  // add mean Uy for lines here...
		//fp111= fopen(val4,"a");
		  // Only one declaration of point-data-size and type
		  //fprintf(fp111, "\n POINT_DATA %d \n", nPoints);
		  fprintf(fp111, "SCALARS Umeany double 1\n");
		  fprintf(fp111, "LOOKUP_TABLE default\n");

		  for(nOfP=1; nOfP <=nPoints; nOfP++)
		  {
			  if(nOfP<nPoints)
			  {
					fprintf(fp111,"%2.5f \n", (veloPos1[nOfP].y1-veloPos1[nOfP-1].y1)/interframeTime);
			  }
			  else
			  {
					fprintf(fp111,"-2.00000 \n");
			  }
		  }
		//fclose(fp111);

		  // add mean Uz for lines here...
		//fp111= fopen(val4,"a");
		  // Only one declaration of point-data-size and type
		  //fprintf(fp111, "\n POINT_DATA %d \n", nPoints);
		  fprintf(fp111, "SCALARS Umeanz double 1\n");
		  fprintf(fp111, "LOOKUP_TABLE default\n");

		  for(nOfP=1; nOfP <=nPoints; nOfP++)
		  {
			  if(nOfP<nPoints)
			  {
					fprintf(fp111,"%2.5f \n", (veloPos1[nOfP].z1-veloPos1[nOfP-1].z1)/interframeTime);
			  }
			  else
			  {
					fprintf(fp111,"-2.00000 \n");
			  }
		  }
		//fclose(fp111);

		  // add mean Uxyz for lines here...
		//fp111= fopen(val4,"a");
		  // Only one declaration of point-data-size and type
		  //fprintf(fp111, "\n POINT_DATA %d \n", nPoints);
		  fprintf(fp111, "SCALARS Particle_Path double 1\n");
		  fprintf(fp111, "LOOKUP_TABLE default\n");

		  for(nOfP=1; nOfP <= nPoints; nOfP++)
		  {
			  tmpValue=((double)(nOfP))/((double)(nPoints));
			  fprintf(fp111,"%2.5f \n", tmpValue);
		  }

		fclose(fp111);
	  
		  // add new scalar_data here: *****
		  //
		  //*****
		  //
		free(veloPos1);

	   }
	 }	  //end of for l-loop 
	free(line1);
	free(line2);
	}	  //end of for i-loop 
   //}
//End of ID-export
}
/*----------------------------------------------------------------*/
   //Export stacks for analysis of patterns
   //ToDo: Read stacksize
if(flag4export[3]==1)
{
    nStack=0;
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, "parameters/velo_trans.par");
		strcat(tmp_path, tmp_file_name);
	fp77=fopen(tmp_path,"r");
	fscanf(fp77, "%f\n", &interframeTime);
	fscanf(fp77, "%d\n", &minSelect);
	fscanf(fp77, "%d\n", &stacksize); 
	//stacksize = 100; //Transfered to PTV GUI
	fclose(fp77);
		strcpy(tmp_path, base_path);
		strcpy(tmp_file_name, "PostProcessing/Stacks/");
		strcat(tmp_path, tmp_file_name);
    CreateDirectoryA(tmp_path,NULL);
   printf("Writing stacks for pattern analysis.\n");
   if(stacksize>=(seq_last-seq_first))stacksize = (seq_last-seq_first);
   printf("Used stacksize: %d\n", stacksize);
   for(i=seq_first; i<seq_last-1; i+=stacksize)
   {
	   nOfP=0;//reset line counter
	   	//Create file to write out stack data of connections
		if (nStack+1 < 10)			sprintf (val3, "PostProcessing/Stacks/Stack_%1d.vtk",nStack+1);
		else if (nStack+1 < 100)	sprintf (val3, "PostProcessing/Stacks/Stack_%2d.vtk",nStack+1);
		else if (nStack+1 < 1000)	sprintf (val3, "PostProcessing/Stacks/Stack_%3d.vtk",nStack+1);
		else						sprintf (val3, "PostProcessing/Stacks/Stack_%4d.vtk",nStack+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val3);
			strcat(tmp_path, tmp_file_name);
		fp99 = fopen(tmp_path,"w");
	    if(!fp99)printf("Could not create file %s.\n",tmp_path);
		//Write header of data file
		fprintf(fp99, "# vtk DataFile Version 2.1 \n");
		fprintf(fp99, "Trajects from stack %d with color option \n",nStack+1);
		fprintf(fp99, "ASCII \n \n");

	   for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
	   {
		//read j ptv_is-file
		if		(j < 10)		sprintf (val, "res/ptv_is.%1d", j);
		else if	(j < 100)		sprintf (val, "res/ptv_is.%2d", j);
		else if	(j < 1000)		sprintf (val, "res/ptv_is.%3d", j);
		else					sprintf (val, "res/ptv_is.%4d", j);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz1);
		line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
		for (k=0;k<anz1;k++) {
		fscanf (fp1, "%d\n",  &line1[k].p);
		fscanf (fp1, "%d\n",  &line1[k].n);
		fscanf (fp1, "%lf\n", &line1[k].x1);
		fscanf (fp1, "%lf\n", &line1[k].y1);
		fscanf (fp1, "%lf\n", &line1[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//read j+1
		if		(j+1 < 10)		sprintf (val, "res/ptv_is.%1d", j+1);
		else if	(j+1 < 100)		sprintf (val, "res/ptv_is.%2d", j+1);
		else if	(j+1 < 1000)	sprintf (val, "res/ptv_is.%3d", j+1);
		else					sprintf (val, "res/ptv_is.%4d", j+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
		for (k=0;k<anz2;k++) {
		fscanf (fp1, "%d\n",  &line2[k].p);
		fscanf (fp1, "%d\n",  &line2[k].n);
		fscanf (fp1, "%lf\n", &line2[k].x1);
		fscanf (fp1, "%lf\n", &line2[k].y1);
		fscanf (fp1, "%lf\n", &line2[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//start operating part
	    // calc/sum up number of lines for stack:
		for(k=0;k<anz1;k++)
		{ 	
		m = line1[k].n;
		if (m >= 0 && m <= anz2){nOfP = nOfP++;}//count connected particles
		}
		//end operating part
		free(line1);
		free(line2);
	   }
	fprintf(fp99, "DATASET POLYDATA POINTS %d float \n", 2*nOfP );
	//veloPos1 = (vector *) malloc ( 500*(nOfP+1) * (sizeof (vector)+1)); 
	//if(!veloPos1) printf("Error malloc veloPos1.\n");
	   for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
	   {
		//read j ptv_is-file
		if		(j < 10)		sprintf (val, "res/ptv_is.%1d", j);
		else if	(j < 100)		sprintf (val, "res/ptv_is.%2d", j);
		else if	(j < 1000)		sprintf (val, "res/ptv_is.%3d", j);
		else					sprintf (val, "res/ptv_is.%4d", j);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz1);
		line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
		for (k=0;k<anz1;k++) {
		fscanf (fp1, "%d\n",  &line1[k].p);
		fscanf (fp1, "%d\n",  &line1[k].n);
		fscanf (fp1, "%lf\n", &line1[k].x1);
		fscanf (fp1, "%lf\n", &line1[k].y1);
		fscanf (fp1, "%lf\n", &line1[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//read j+1
		if		(j+1 < 10)		sprintf (val, "res/ptv_is.%1d", j+1);
		else if	(j+1 < 100)		sprintf (val, "res/ptv_is.%2d", j+1);
		else if	(j+1 < 1000)	sprintf (val, "res/ptv_is.%3d", j+1);
		else					sprintf (val, "res/ptv_is.%4d", j+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
		for (k=0;k<anz2;k++) {
		fscanf (fp1, "%d\n",  &line2[k].p);
		fscanf (fp1, "%d\n",  &line2[k].n);
		fscanf (fp1, "%lf\n", &line2[k].x1);
		fscanf (fp1, "%lf\n", &line2[k].y1);
		fscanf (fp1, "%lf\n", &line2[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//start operating part
	    // write points of lines for stack:
		for(k=0;k<anz1;k++)
		{ 	
		 m = line1[k].n;
		 if (m >= 0 && m <= anz2)
		 {
			if( line1[k].x1 < -100000 || line1[k].y1 < -100000 || line1[k].z1 < -100000||
				line2[m].x1 < -100000 || line2[m].y1 < -100000 || line2[m].z1 < -100000)
			{
				line1[k].x1=0.0;	line1[k].y1=0.0;	line1[k].z1=0.0;
				line2[m].x1=0.0;	line2[m].y1=0.0;	line2[m].z1=0.0;
				fprintf(fp99, "%5.3f %5.3f %5.3f \n", line1[k].x1 + x_offset, line1[k].y1 + y_offset, line1[k].z1 + z_offset);
				fprintf(fp99, "%5.3f %5.3f %5.3f \n", line2[m].x1 + x_offset, line2[m].y1 + y_offset, line2[m].z1 + z_offset); 
			}
			else if(line1[k].x1 > 100000 || line1[k].y1 > 100000 || line1[k].z1 > 100000|| 
					line2[m].x1 > 100000 || line2[m].y1 > 100000 || line2[m].z1 > 100000)
			{
				line1[k].x1=0.0;	line1[k].y1=0.0;	line1[k].z1=0.0;
				line2[m].x1=0.0;	line2[m].y1=0.0;	line2[m].z1=0.0;
				fprintf(fp99, "%5.3f %5.3f %5.3f \n", line1[k].x1 + x_offset, line1[k].y1 + y_offset, line1[k].z1 + z_offset);
				fprintf(fp99, "%5.3f %5.3f %5.3f \n", line2[m].x1 + x_offset, line2[m].y1 + y_offset, line2[m].z1 + z_offset); 
			}
			else
			{
				 fprintf(fp99, "%5.3f %5.3f %5.3f \n", line1[k].x1 + x_offset, line1[k].y1 + y_offset, line1[k].z1 + z_offset);
				 fprintf(fp99, "%5.3f %5.3f %5.3f \n", line2[m].x1 + x_offset, line2[m].y1 + y_offset, line2[m].z1 + z_offset); 
			}
		 }
		}
		//end operating part
		free(line1);
		free(line2);
	   }
	 nVerts=0;
	 nVerts=3*nOfP;
	 fprintf(fp99, "\n LINES %d %d \n",nOfP, nVerts);
	 connector = 0;
	 for(k=0;k<nOfP;k++) //establish connections from Points
	 { 	
		fprintf(fp99,"2 %d %d \n", connector, connector+1);
		connector = connector+2;
	 }
	 fprintf(fp99, "POINT_DATA %d \n", 2*nOfP);
	//append mean veloData to "lines"
	fprintf(fp99, "SCALARS Umean double 1 \n");
	fprintf(fp99, "LOOKUP_TABLE default\n");
    for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
	   {
		//read j ptv_is-file
		if		(j < 10)		sprintf (val, "res/ptv_is.%1d", j);
		else if	(j < 100)		sprintf (val, "res/ptv_is.%2d", j);
		else if	(j < 1000)		sprintf (val, "res/ptv_is.%3d", j);
		else					sprintf (val, "res/ptv_is.%4d", j);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz1);
		line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
		for (k=0;k<anz1;k++) {
		fscanf (fp1, "%d\n",  &line1[k].p);
		fscanf (fp1, "%d\n",  &line1[k].n);
		fscanf (fp1, "%lf\n", &line1[k].x1);
		fscanf (fp1, "%lf\n", &line1[k].y1);
		fscanf (fp1, "%lf\n", &line1[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//read j+1
		if		(j+1 < 10)		sprintf (val, "res/ptv_is.%1d", j+1);
		else if	(j+1 < 100)		sprintf (val, "res/ptv_is.%2d", j+1);
		else if	(j+1 < 1000)	sprintf (val, "res/ptv_is.%3d", j+1);
		else					sprintf (val, "res/ptv_is.%4d", j+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
		for (k=0;k<anz2;k++) {
		fscanf (fp1, "%d\n",  &line2[k].p);
		fscanf (fp1, "%d\n",  &line2[k].n);
		fscanf (fp1, "%lf\n", &line2[k].x1);
		fscanf (fp1, "%lf\n", &line2[k].y1);
		fscanf (fp1, "%lf\n", &line2[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//start operating part
	    // write points of lines for stack:
		for(k=0;k<anz1;k++)
		{ 	
		 m = line1[k].n;
		 if (m >= 0 && m <= anz2)
		 {   		//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
				dx=line1[k].x1-line2[m].x1;
				dy=line1[k].y1-line2[m].y1;
				dz=line1[k].z1-line2[m].z1;
				du=sqrt(dx*dx+dy*dy);
				dl=sqrt(dx*dx+dy*dy+dz*dz);
				MeanU=dl/interframeTime;
				Ux=(-dx/interframeTime);
				Uy=(-dy/interframeTime);
				Uz=(-dz/interframeTime);
			    //print mean U with check of max Velo
		        if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))Um=-1;
		        else Um=MeanU;
		 fprintf(fp99,"%2.5f \n", Um);
		 fprintf(fp99,"%2.5f \n", Um);
		 }
		}
		//end operating part
		free(line1);
		free(line2);
	   }
	//append mean x-veloData to "lines"
	fprintf(fp99, "SCALARS Umean_x double 1 \n");
	fprintf(fp99, "LOOKUP_TABLE default\n");
    for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
	   {
		//read j ptv_is-file
		if		(j < 10)		sprintf (val, "res/ptv_is.%1d", j);
		else if	(j < 100)		sprintf (val, "res/ptv_is.%2d", j);
		else if	(j < 1000)		sprintf (val, "res/ptv_is.%3d", j);
		else					sprintf (val, "res/ptv_is.%4d", j);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");    
		fscanf (fp1,"%d\n", &anz1);
		line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
		for (k=0;k<anz1;k++) {
		fscanf (fp1, "%d\n",  &line1[k].p);
		fscanf (fp1, "%d\n",  &line1[k].n);
		fscanf (fp1, "%lf\n", &line1[k].x1);
		fscanf (fp1, "%lf\n", &line1[k].y1);
		fscanf (fp1, "%lf\n", &line1[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//read j+1
		if		(j+1 < 10)		sprintf (val, "res/ptv_is.%1d", j+1);
		else if	(j+1 < 100)		sprintf (val, "res/ptv_is.%2d", j+1);
		else if	(j+1 < 1000)	sprintf (val, "res/ptv_is.%3d", j+1);
		else					sprintf (val, "res/ptv_is.%4d", j+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
		for (k=0;k<anz2;k++) {
		fscanf (fp1, "%d\n",  &line2[k].p);
		fscanf (fp1, "%d\n",  &line2[k].n);
		fscanf (fp1, "%lf\n", &line2[k].x1);
		fscanf (fp1, "%lf\n", &line2[k].y1);
		fscanf (fp1, "%lf\n", &line2[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//start operating part
	    // write points of lines for stack:
		for(k=0;k<anz1;k++)
		{ 	
		 m = line1[k].n;
		 if (m >= 0 && m <= anz2)
		 {   		//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
				dx=line1[k].x1-line2[m].x1;
				dy=line1[k].y1-line2[m].y1;
				dz=line1[k].z1-line2[m].z1;
				du=sqrt(dx*dx+dy*dy);
				dl=sqrt(dx*dx+dy*dy+dz*dz);
				MeanU=dl/interframeTime;
				Ux=(-dx/interframeTime);
				Uy=(-dy/interframeTime);
				Uz=(-dz/interframeTime);
			    //print mean U with check of max Velo
				if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))Um=-1;
		        else Um=MeanU;
		 fprintf(fp99,"%2.5f \n", Ux);
		 fprintf(fp99,"%2.5f \n", Ux);
		 }
		}
		//end operating part
		free(line1);
		free(line2);
	   }
	//append mean y-veloData to "lines"
	fprintf(fp99, "SCALARS Umean_y double 1 \n");
	fprintf(fp99, "LOOKUP_TABLE default\n");
    for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
	   {
		//read j ptv_is-file
		if		(j < 10)		sprintf (val, "res/ptv_is.%1d", j);
		else if	(j < 100)		sprintf (val, "res/ptv_is.%2d", j);
		else if	(j < 1000)		sprintf (val, "res/ptv_is.%3d", j);
		else					sprintf (val, "res/ptv_is.%4d", j);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz1);
		line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
		for (k=0;k<anz1;k++) {
		fscanf (fp1, "%d\n",  &line1[k].p);
		fscanf (fp1, "%d\n",  &line1[k].n);
		fscanf (fp1, "%lf\n", &line1[k].x1);
		fscanf (fp1, "%lf\n", &line1[k].y1);
		fscanf (fp1, "%lf\n", &line1[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//read j+1
		if		(j+1 < 10)		sprintf (val, "res/ptv_is.%1d", j+1);
		else if	(j+1 < 100)		sprintf (val, "res/ptv_is.%2d", j+1);
		else if	(j+1 < 1000)	sprintf (val, "res/ptv_is.%3d", j+1);
		else					sprintf (val, "res/ptv_is.%4d", j+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
		for (k=0;k<anz2;k++) {
		fscanf (fp1, "%d\n",  &line2[k].p);
		fscanf (fp1, "%d\n",  &line2[k].n);
		fscanf (fp1, "%lf\n", &line2[k].x1);
		fscanf (fp1, "%lf\n", &line2[k].y1);
		fscanf (fp1, "%lf\n", &line2[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//start operating part
	    // write points of lines for stack:
		for(k=0;k<anz1;k++)
		{ 	
		 m = line1[k].n;
		 if (m >= 0 && m <= anz2)
		 {   		//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
				dx=line1[k].x1-line2[m].x1;
				dy=line1[k].y1-line2[m].y1;
				dz=line1[k].z1-line2[m].z1;
				du=sqrt(dx*dx+dy*dy);
				dl=sqrt(dx*dx+dy*dy+dz*dz);
				MeanU=dl/interframeTime;
				Ux=(-dx/interframeTime);
				Uy=(-dy/interframeTime);
				Uz=(-dz/interframeTime);
			    //print mean U with check of max Velo
				if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))Um=-1;
		        else Um=MeanU;
		 fprintf(fp99,"%2.5f \n", Uy);
		 fprintf(fp99,"%2.5f \n", Uy);
		 }
		}
		//end operating part
		free(line1);
		free(line2);
	   }
	//append mean z-veloData to "lines"
	fprintf(fp99, "SCALARS Umean_z double 1 \n");
	fprintf(fp99, "LOOKUP_TABLE default\n");
    for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
	   {
		//read j ptv_is-file
		if		(j < 10)		sprintf (val, "res/ptv_is.%1d", j);
		else if	(j < 100)		sprintf (val, "res/ptv_is.%2d", j);
		else if	(j < 1000)		sprintf (val, "res/ptv_is.%3d", j);
		else					sprintf (val, "res/ptv_is.%4d", j);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");    
		fscanf (fp1,"%d\n", &anz1);
		line1 = (vector *) malloc ((anz1+1)*sizeof (vector));
		for (k=0;k<anz1;k++) {
		fscanf (fp1, "%d\n",  &line1[k].p);
		fscanf (fp1, "%d\n",  &line1[k].n);
		fscanf (fp1, "%lf\n", &line1[k].x1);
		fscanf (fp1, "%lf\n", &line1[k].y1);
		fscanf (fp1, "%lf\n", &line1[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//read j+1
		if		(j+1 < 10)		sprintf (val, "res/ptv_is.%1d", j+1);
		else if	(j+1 < 100)		sprintf (val, "res/ptv_is.%2d", j+1);
		else if	(j+1 < 1000)	sprintf (val, "res/ptv_is.%3d", j+1);
		else					sprintf (val, "res/ptv_is.%4d", j+1);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector *) malloc ((anz2+1)*sizeof (vector));
		for (k=0;k<anz2;k++) {
		fscanf (fp1, "%d\n",  &line2[k].p);
		fscanf (fp1, "%d\n",  &line2[k].n);
		fscanf (fp1, "%lf\n", &line2[k].x1);
		fscanf (fp1, "%lf\n", &line2[k].y1);
		fscanf (fp1, "%lf\n", &line2[k].z1);
        }
        strcpy(val, "");     
        fclose (fp1);
		//start operating part
	    // write points of lines for stack:
		for(k=0;k<anz1;k++)
		{ 	
		 m = line1[k].n;
		 if (m >= 0 && m <= anz2)
		 {   		//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
				dx=line1[k].x1-line2[m].x1;
				dy=line1[k].y1-line2[m].y1;
				dz=line1[k].z1-line2[m].z1;
				du=sqrt(dx*dx+dy*dy);
				dl=sqrt(dx*dx+dy*dy+dz*dz);
				MeanU=dl/interframeTime;
				Ux=(-dx/interframeTime);
				Uy=(-dy/interframeTime);
				Uz=(-dz/interframeTime);
			    //print mean U with check of max Velo
				if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))Um=-1;
		        else Um=MeanU;
		 fprintf(fp99,"%2.5f \n", Uz);
		 fprintf(fp99,"%2.5f \n", Uz);
		 }
		}
		//end operating part
		free(line1);
		free(line2);
	   }


	   fclose(fp99);
	   nStack++;//counter for stacks
   }
 }

/*----------------------------------------------------------------*/


   /**Added export to paraview for correspondending points**/
   if(flag4export[0]==1)//option "a"
   {
	 for (i=seq_first; i<seq_last;i++)
	 {
	   if		(i+1 < 10)		sprintf (val, "res/rt_is.%1d", i+1);
	   else if	(i+1 < 100)		sprintf (val, "res/rt_is.%2d",  i+1);
	   else if	(i+1 < 1000)	sprintf (val, "res/rt_is.%3d",  i+1);
	   else						sprintf (val, "res/rt_is.%4d",  i+1);
      
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
		fp1 = fopen (tmp_path, "r");     
	   fscanf (fp1,"%d\n", &anz);
	   line1 = (vector *) malloc ((anz+1)*sizeof (vector));
      
	   for (j=0;j<anz;j++) //Read data to vector
	   {
		fscanf (fp1, "%d\n", &dummy); 
		fscanf (fp1, "%lf\n", &line1[j].x1);
		fscanf (fp1, "%lf\n", &line1[j].y1);
		fscanf (fp1, "%lf\n", &line1[j].z1);
		fscanf (fp1, "%d\n", &dummy); 
		fscanf (fp1, "%d\n", &dummy); 
		fscanf (fp1, "%d\n", &dummy);
		fscanf (fp1, "%d\n", &dummy); 
	   }
	   fclose (fp1);

	   //Print data to file
		if		(i < 10)	sprintf (val2, "PostProcessing/prd_3dPoints_%1d.vtk", i);
		else if (i < 100)   sprintf (val2, "PostProcessing/prd_3dPoints_%2d.vtk", i);
		else if (i < 1000)  sprintf (val2, "PostProcessing/prd_3dPoints_%3d.vtk", i);
		else				sprintf (val2, "PostProcessing/prd_3dPoints_%4d.vtk", i);
		printf("Create VTK-Data of 3d-points, create file: %s\n", val2);
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val2);
			strcat(tmp_path, tmp_file_name);
		fp66 = fopen(tmp_path,"w");
		//write data to vtk file in simple legacy format
		fprintf(fp66, "# vtk DataFile Version 2.0 \n");
		fprintf(fp66, "Xm, Ym, Zm, Uvecs, timestep %d \n",i);
		fprintf(fp66, "ASCII \n \n");
		fprintf(fp66, "DATASET POLYDATA POINTS %d float \n", anz+2 );
		fprintf(fp66, "00000.000 00000.000 00000.000 \n");
		fprintf(fp66, "00001.000 00001.000 00001.000 \n");
		for(j=0;j<anz;j++)
		{
			fprintf(fp66, "%5.3f %5.3f %5.3f \n", line1[j].x1+x_offset, line1[j].y1+y_offset, line1[j].z1+z_offset);
		}
		fprintf(fp66, "\n VERTICES %d %d \n", anz+2, 2*(anz+2));
		fprintf(fp66, "1 0 \n");
		fprintf(fp66, "1 0 \n");
		for(j=0;j<anz;j++)
		{
			fprintf(fp66, "1 %d \n", j);
		}
		fprintf(fp66, "POINT_DATA %d \n", anz+2);
		fprintf(fp66, "SCALARS pnr double 1 \n");
		fprintf(fp66, "LOOKUP_TABLE default\n");
		fprintf(fp66, "0 \n0 \n");
		for(j=0;j<anz;j++)
		{
			fprintf(fp66, "%d \n", j);
		}

		fclose (fp66);
		free(line1);

	 }
 }
 /*End of point-export*/
 printf("Export finished.\n");
  sprintf(val, "...done");
  //Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  //Tcl_Eval(interp, ".text delete 4");
  //Tcl_Eval(interp, ".text insert 4 $tbuf");
  
  //return TCL_OK;
  return 1;
}

void readseqtrackcrit (char base_path[256])
{
  int i_img;
  char tmp_path[256];
  /* reads pixfiles and try to track particles in imagespace
     over the sequence */
  strcpy(tmp_path, base_path);
  strcat(tmp_path, sequence_par);
  fpp = fopen(tmp_path,"r");
  for (i_img=0; i_img<4; i_img++) {
    fscanf (fpp, "%s\n", seq_name[i_img]);
  }
  /* name of sequence */
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);

  strcpy(tmp_path, base_path);
  strcat(tmp_path, track_par);
  fpp = fopen(tmp_path,"r");
  fscanf (fpp, "%lf\n", &tpar.dvxmin);
  fscanf (fpp, "%lf\n", &tpar.dvxmax);
  fscanf (fpp, "%lf\n", &tpar.dvymin);
  fscanf (fpp, "%lf\n", &tpar.dvymax);
  fscanf (fpp, "%lf\n", &tpar.dvzmin);
  fscanf (fpp, "%lf\n", &tpar.dvzmax);
  fscanf (fpp, "%lf\n", &tpar.dangle);
  fscanf (fpp, "%lf\n", &tpar.dacc);
  fscanf (fpp,"%d\n", &tpar.add);
  fscanf (fpp,"%i\n", &tpar.nframe); //max nParticle
  fscanf (fpp,"%i\n", &tpar.deltaT); //delta t
  fscanf (fpp,"%f\n", &tpar.Vmax); //max. velo
  fscanf (fpp,"%f\n", &tpar.link_tol); // linking tol.
  fscanf (fpp,"%f\n", &tpar.jump_tol); // jumping tol.
  fscanf (fpp,"%i\n", &tpar.flag[0]); //gluing
  fscanf (fpp,"%i\n", &tpar.flag[1]); //track
  fscanf (fpp,"%i\n", &tpar.flag[2]); //iterate
  fscanf (fpp,"%i\n", &tpar.flag[3]); //remove Kints
  fscanf (fpp,"%i\n", &tpar.flag[4]); //write binary
  fscanf (fpp,"%i\n", &tpar.flag[5]); //write Traj
  fscanf (fpp,"%i\n", &tpar.flag[6]); //output
  fscanf (fpp,"%i\n", &tpar.flag[7]); //sameCamConf
  //flag for printing searchquader & points
  fscanf (fpp,"%i\n", &tpar.showAdds);
  /*
    fscanf (fpp,"%d\n", &tpar.dsumg);
    fscanf (fpp,"%d\n", &tpar.dn);
    fscanf (fpp,"%d\n", &tpar.dnx);
    fscanf (fpp,"%d\n", &tpar.dny);
  */
  fclose (fpp);
  tpar.glue=tpar.flag[0];
  tpar.track=tpar.flag[1];
  tpar.iter=tpar.flag[2];
  tpar.remKints=tpar.flag[3];
  tpar.writeBin=tpar.flag[4];
  tpar.writeTraj=tpar.flag[5];
  tpar.out=tpar.flag[6];
  tpar.sameCamConf=tpar.flag[7];

  /* read illuminated layer data */
  strcpy(tmp_path, base_path);
  strcat(tmp_path, criteria_par);
  fpp = fopen(tmp_path,"r");
  fscanf (fpp, "%lf\n", &X_lay[0]);
  fscanf (fpp, "%lf\n", &Zmin_lay[0]);
  fscanf (fpp, "%lf\n", &Zmax_lay[0]);
  fscanf (fpp, "%lf\n", &X_lay[1]);
  fscanf (fpp, "%lf\n", &Zmin_lay[1]);
  fscanf (fpp, "%lf\n", &Zmax_lay[1]);
  fclose (fpp);
}