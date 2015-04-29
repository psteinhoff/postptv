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


int main(int argc,char **argv)
    //int argc;			/* Number of command-line arguments. */
    //char **argv;		/* Values of command-line arguments. */
{
	char Buffer[80];
	int count, size;
	int start;
	int Exportflag[5];
	int i;
	FILE *fpp;
	sprintf(path_to_WF,"");//if useful: adjust your default path to your working-folder here.
	/*strcpy(velo_par,		"parameters/velo_trans.par");
	strcpy(sequence_par,	"parameters/sequence.par");
	strcpy(track_par,		"parameters/track.par");
	strcpy(criteria_par,	"parameters/criteria.par");
	strcpy(point_offset,	"parameters/3d_point_offset.par");
	 */
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
	 for(i=0;i<5;i++)Exportflag[i]=0;
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

				case 'e':
				fprintf( stderr, "Enabling export of trajectories sorted by ID.\n" );
				//flag for trajects
				Exportflag[4]=1;
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
						"        (i)Export of trajects between two consecutive images.\n"
						"   -d"
						"        Enable export of stacks of trajectories.\n\n"
						"   -e"
						"        Enable export of trajectories.\n"
						"        (ii)Export of trajects in to folders according to their length\n\n"
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
			clearLine(18,70);
			PostToParaview_c (argc,argv, path_to_WF, Exportflag);
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
	int dummy, count;
	FILE *fp1, *fp2,*fp22, *fp88, *fp66, *fp77, *fp99,*fp55, *fp111;
	char val[256],val1[256], val2[256], val3[256],val4[256];
	char resultFolder[256], path[256], selectionPath[256], tmp_path[512],tmp_file_name[256],tmpPath[512];
	vector *line1, *line2, *line3, *veloPos1;
	double color, ymin=0, ymax=0, cubes;
	double mx, my, mz, dx, dy, dz, du, ux, uy,uz, dl, rotz, rotx;
	double Ux,Uy,Uz,Um;
	int connector, nStack, stacksize;
	double MeanU, tmpValue;
	float interframeTime;
	double U_meanMax;
	int stackOnly = 0;
	vector **pointDataArray, *timeStep;


	sprintf(resultFolder, "res/");

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
	fclose(fp88);

	strcpy(tmp_path, base_path);
	strcpy(tmp_file_name, "parameters/velo_trans.par");
	strcat(tmp_path, tmp_file_name);

	fp77=fopen(tmp_path,"r");
		fscanf(fp77, "%f\n", &interframeTime);	//time between two frames in ms to calc velocity 
		fscanf(fp77, "%d\n", &minSelect);		//minimal length of trajectory for processing
	fclose(fp77);

	//allocate Memory and read all ptv_is-files into pointDataArray
	pointDataArray =(vector **) malloc( (seq_last-seq_first+1)*sizeof(timeStep));
	if (pointDataArray)
	{	
		for(i=0;i < (seq_last-seq_first+1); i++)
		{
			compose_name_plus_nr (resultFolder, "ptv_is.",seq_first + i, val);  

			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);

			fp1 = fopen (tmp_path, "r");
			if (! fp1)  showErrorMessageOnFileOpening(tmp_path);
			else
			{
				fscanf (fp1,"%d\n", &anz1);    
				pointDataArray[i] = (vector *)malloc((anz1) *sizeof(*pointDataArray[i]));
				pointDataArray[i][0].type=anz1;
				for (j=0;j<anz1;j++)
					{
						fscanf (fp1, "%d\n",	&pointDataArray[i][j].p);
						fscanf (fp1, "%d\n",	&pointDataArray[i][j].n);
						fscanf (fp1, "%lf\n",	&pointDataArray[i][j].x1);
						fscanf (fp1, "%lf\n",	&pointDataArray[i][j].y1);
						fscanf (fp1, "%lf\n",	&pointDataArray[i][j].z1);
					}
				strcpy(val, "");     
				fclose (fp1);
			}
			progressBar(18,60,(int)((i*100)/(seq_last-seq_first)));	
		}	
	}
	printf("Succesfully loaded all ptv.is-files.\n");
	clearLine(18,70);
	//process your data here:
	if(flag4export[1]==1)
	{
		//**********// create prd_Umean-files //**********//
		printf("Write points & velocities to vtk-Files.\n");
		for(i=0;i < (seq_last-seq_first); i++)
			{	
				//calc number of lines:
				nOfP = countConnectedParticles (pointDataArray[i][0].type,		pointDataArray[i], 
												pointDataArray[i+1][0].type,	pointDataArray[i+1]);
				compose_name_plus_nr_str("PostProcessing/prd_Umean_",".vtk",seq_first+i,val2);
				fp66 = fopen(val2,"w");
				if (! fp66)  showErrorMessageOnFileOpening(val2);
				//write data to vtk file in simple legacy format
				else
				{
					fprintf(fp66, "# vtk DataFile Version 2.0 \n");
					fprintf(fp66, "Xm, Ym, Zm, Uvecs, timestep %d \n",seq_first+i);
					fprintf(fp66, "ASCII \n \n");  
					fprintf(fp66, "DATASET POLYDATA POINTS %d float \n", nOfP );
					//write pont data
					for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type && pointDataArray[i+1][m].p>0) 
						{	
							/* create cubes of last time step */	 
							if(i ==seq_last-1) {	  }
	  
							/* create last cubes of ending trajectories  */
							if(pointDataArray[i+1][m].n <= 0) {	  }

							//calc middle and translation
							mx=(pointDataArray[i][j].x1+pointDataArray[i+1][m].x1)/2 + x_offset;
							my=(pointDataArray[i][j].y1+pointDataArray[i+1][m].y1)/2 + y_offset; 
							mz=(pointDataArray[i][j].z1+pointDataArray[i+1][m].z1)/2 + z_offset; 
							//be careful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
							fprintf(fp66, "%7.3f %7.3f %7.3f \n", mx, my, mz);
	
						}	
					}

					//append infomation about vertices
					nVerts=2*nOfP;
					fprintf(fp66, "\n VERTICES %d %d \n",nOfP, nVerts);
					count=0;
					for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type && pointDataArray[i+1][m].p>0) 
							{
								//count connected particles
								fprintf(fp66,"1 %d \n", count);
								count = count++;
							}
					}

					//append velocity-data part 1
					fprintf(fp66, "POINT_DATA %d \n", nOfP);
					fprintf(fp66, "SCALARS MeanU double 1 \n");
					fprintf(fp66, "LOOKUP_TABLE default\n");
					for(j=0;j<pointDataArray[i][0].type;j++)
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type && pointDataArray[i+1][m].p>0) 
						{		 /* create cubes of last time step */	 
								if(i ==seq_last-1) {   }
	  							/* create last cubes of ending trajectories  */
								if(pointDataArray[i+1][m].n <= 0) {	 }

								mx = (pointDataArray[i][j].x1 + pointDataArray[i+1][m].x1) /2;
								my = (pointDataArray[i][j].y1 + pointDataArray[i+1][m].y1) /2; 
								mz = (pointDataArray[i][j].z1 + pointDataArray[i+1][m].z1) /2; 
								//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...
								dx = pointDataArray[i][j].x1 - pointDataArray[i+1][m].x1;
								dy = pointDataArray[i][j].y1 - pointDataArray[i+1][m].y1;
								dz = pointDataArray[i][j].z1 - pointDataArray[i+1][m].z1;
								du = sqrt(dx*dx+dy*dy);
								dl = sqrt(dx*dx+dy*dy+dz*dz);
								MeanU = dl/interframeTime;
								//print mean U with check of max Velo
								if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))fprintf(fp66, "-1.00000\n"); //write artificial value if value is out of specified  range
								else fprintf(fp66, "%3.5f \n", MeanU);	   
						}	
					}

					//append velocity data 2 (vecs)
					fprintf(fp66,"\n VECTORS Umean float \n");
					for(j=0;j<pointDataArray[i][0].type;j++)
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type && pointDataArray[i+1][m].p>0) 
						{	
							/* create cubes of last time step */	 
							if(i ==seq_last-1) { }
							/* create last cubes of ending trajectories  */
							if(pointDataArray[i+1][m].n <= 0) { }

							dx = pointDataArray[i][j].x1 - pointDataArray[i+1][m].x1;
							dy = pointDataArray[i][j].y1 - pointDataArray[i+1][m].y1;
							dz = pointDataArray[i][j].z1 - pointDataArray[i+1][m].z1;
							//be careuful!!! moving backwards! 1-2, 2-3, 3-4... n-(n+1)....means moving backwards...

							ux= -dx/interframeTime;
							uy= -dy/interframeTime;
							uz= -dz/interframeTime;
							//print velos
							fprintf(fp66, "%3.5f %3.5f %3.5f \n", ux, uy,uz); 
						}	
					}

					//append particle-data part 1
					fprintf(fp66, "SCALARS NextParticleID double 1 \n");
					fprintf(fp66, "LOOKUP_TABLE default\n");
					for(j=0;j<pointDataArray[i][0].type;j++)
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type && pointDataArray[i+1][m].p>0) 
						{	
							fprintf(fp66, "%d \n", m);
						}	
					}

					//append particle-data part 2
					fprintf(fp66, "SCALARS PreviousParticleID double 1 \n");
					fprintf(fp66, "LOOKUP_TABLE default\n");
					for(j=0;j<pointDataArray[i][0].type;j++)
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type && pointDataArray[i+1][m].p>0) 
						{	
							fprintf(fp66, "%d \n", pointDataArray[i][j].p);   	
						}	
					}

					//perfom next action with data


				fclose(fp66);

				}
				progressBar(18,60,(int)((i*100)/(seq_last-seq_first)));	
			}
		printf("Points & Velocities written to vtk-Files.\n");
		clearLine(18,70);
	}
//**********************************************************************//
	//**********// create prd_Umean-files //**********//
	if(flag4export[2] ==1)
	{
		printf("Write Trajectories.\n");
		CreateDirectoryA("PostProcessing/Trajects/",NULL);
		for(i=0;i < (seq_last-seq_first); i++)
			{
				compose_name_plus_nr_str("PostProcessing/Trajects/Traject_",".vtk",seq_first + i,val3);
				fp99 = fopen(val3,"w");
					if (! fp99)  showErrorMessageOnFileOpening(val3);
				//write data to vtk file in simple legacy format
				fprintf(fp99, "# vtk DataFile Version 2.1 \n");
				fprintf(fp99, "Trajects from timestep %d with color option \n",i);
				fprintf(fp99, "ASCII \n \n");
			
				//nOfP = countConnectedParticles (pointDataArray[i][0].type,		pointDataArray[i], 
				//									pointDataArray[i+1][0].type,	pointDataArray[i+1]);

				nOfP = countValidatedParticleConnections(pointDataArray[i], pointDataArray[i+1][0].type);
				veloPos1 = (vector *) malloc((pointDataArray[i][0].type+1)*sizeof(vector)); //Changed allocated memory to 500*xyz due to problems in allocating enough memory and writing data
				if(!veloPos1)printf("Error malloc veloPos1.\n");

				fprintf(fp99, "DATASET POLYDATA POINTS %d float \n", 2*nOfP );
				for(j=0;j<pointDataArray[i][0].type;j++)
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type) 
						{	
								/* create cubes of last time step */	 
								if(i ==seq_last-1) {	  }
	  
								/* create last cubes of ending trajectories  */
								if(pointDataArray[i+1][m].n <= 0) {	  }

								fprintf(fp99, "%5.3f %5.3f %5.3f \n",
									pointDataArray[i][j].x1 + x_offset, 
									pointDataArray[i][j].y1 + y_offset, 
									pointDataArray[i][j].z1 + z_offset);
								fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
									pointDataArray[i+1][m].x1 + x_offset, 
									pointDataArray[i+1][m].y1 + y_offset, 
									pointDataArray[i+1][m].z1 + z_offset);

								dx = pointDataArray[i][j].x1 - pointDataArray[i+1][m].x1;
								dy = pointDataArray[i][j].y1 - pointDataArray[i+1][m].y1;
								dz = pointDataArray[i][j].z1 - pointDataArray[i+1][m].z1;

								dl = sqrt(dx*dx+dy*dy+dz*dz);
								MeanU=dl/interframeTime;
								veloPos1[j].Ux=(-dx/interframeTime);
								veloPos1[j].Uy=(-dy/interframeTime);
								veloPos1[j].Uz=(-dz/interframeTime);
							
								if (MeanU > 2*sqrt(3*(tpar.Vmax)*(tpar.Vmax)))veloPos1[j].Um=-1;
								else veloPos1[j].Um=MeanU;

						}	
					}
				nVerts=3*nOfP;
				fprintf(fp99, "\n LINES %d %d \n",nOfP, nVerts);
				count=0;
				for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type)
						{
							fprintf(fp99,"2 %d %d \n", count, count+1);
							count = count+2;
						}
					}

				fprintf(fp99, "POINT_DATA %d \n", 2*nOfP);

				//append mean veloData to "lines"
				fprintf(fp99, "SCALARS Umean double 1 \n");
				fprintf(fp99, "LOOKUP_TABLE default\n");
				for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type) 
						{
							fprintf(fp99,"%2.5f \n", veloPos1[j].Um);
							fprintf(fp99,"%2.5f \n", veloPos1[j].Um);
						}
					}

				//append mean x-veloData to "lines"
				fprintf(fp99, "SCALARS Umean_X double 1 \n");
				fprintf(fp99, "LOOKUP_TABLE default\n");
				for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type) 
						{
							fprintf(fp99,"%2.5f \n", veloPos1[j].Ux);
							fprintf(fp99,"%2.5f \n", veloPos1[j].Ux);
						}
					}

				//append mean y-veloData to "lines"
				fprintf(fp99, "SCALARS Umean_Y double 1 \n");
				fprintf(fp99, "LOOKUP_TABLE default\n");
				for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type) 
						{
							fprintf(fp99,"%2.5f \n", veloPos1[j].Uy);
							fprintf(fp99,"%2.5f \n", veloPos1[j].Uy);
						}
					}

				//append mean z-veloData to "lines"
				fprintf(fp99, "SCALARS Umean_Z double 1 \n");
				fprintf(fp99, "LOOKUP_TABLE default\n");
				for(j=0;j<pointDataArray[i][0].type;j++) 
					{ 	
						m = pointDataArray[i][j].n;
						if (m >= 0 && m < pointDataArray[i+1][0].type) 
						{
							fprintf(fp99,"%2.5f \n", veloPos1[j].Uz);
							fprintf(fp99,"%2.5f \n", veloPos1[j].Uz);
						}
					}
				fclose(fp99);
				free (veloPos1);
				progressBar(18,60,(int)((i*100)/(seq_last-seq_first)));	
			}
		printf("Trajectories written to vtk-Files.\n");
		clearLine(18,70);
	}
//**********************************************************************//
	//**********// Export trajectories sorted with/by ID //**********//	
	if(flag4export[4]==1)
	{
		printf("Export trajectories sorted by ID.\n");
		/*idea: look for starting tracks in each time step and trace them until end*/
		q = 0;
		veloPos1 = (vector *) malloc ((seq_last - seq_first +1)*sizeof (vector));
			if(!veloPos1)printf("Error malloc veloPos1.\n");
		for (i=0; i< (seq_last-seq_first);i++)
		{
			//printf("Create VTK-Data for trajectories, read file: %s\n", val);   
			for(l=0;l<pointDataArray[i][0].type;l++)
			{
				m=pointDataArray[i][l].n;
				if(m >= 0 && m <= pointDataArray[i+1][0].type && pointDataArray[i][l].p < 0) //to check if a new traject beginns
				{
					tmpVeloPos	=	1;
					nPoints		=	1;
					q			=	q+1; //every trajectory gets an individual ID
					lineOfP		=	pointDataArray[i][l].n;

					for(timeID = i+1; timeID<= (seq_last-seq_first) && lineOfP >=0 && lineOfP <= pointDataArray[timeID][0].type;timeID++)
					{
						time = timeID;
					
						if(pointDataArray[time][lineOfP].n <0 && pointDataArray[time][lineOfP].p>=0)//handle last Point of traject
						{	}
						if(lineOfP >= 0 && lineOfP <= pointDataArray[time][0].type) 
						{
							lineOfP=pointDataArray[time][lineOfP].n;		//update line of next particle-position
							nPoints=nPoints+1;				//for each valid link set counter for length of traject +1
						}
					}
				}
				strcpy(tmpPath,"");
				strcpy(path,"");
				compose_name_plus_nr_str("PostProcessing/Trajects/Traject_time_","/",seq_first+i,path);
				if(nPoints>=minSelect)
				{
					compose_name_plus_nr_str("PostProcessing/Trajects/Selection","/",nPoints,path);
					CreateDirectoryA(path,NULL);

					sprintf(tmpPath,"%sTraject_time_",path);						//PostProcessing/Trajects/Selection#/Traject_time_
					compose_name_plus_nr_str(tmpPath,"_length_",seq_first+i,path);	//PostProcessing/Trajects/Selection#/Traject_time_#_length_
					compose_name_plus_nr_str(path,"_ID_",nPoints,tmpPath);			//PostProcessing/Trajects/Selection#/Traject_time_#_length_#_ID_
					compose_name_plus_nr_str(tmpPath,".vtk",q,val4);				//PostProcessing/Trajects/Selection#/Traject_time_#_length_#_ID_#.vtk
					  
				}//end if nPoints>minSelect
				else //write to alternative directory
				{
					//Create time directories:
					compose_name_plus_nr_str("PostProcessing/Trajects/Traject_time_","/",seq_first+i,path); //PostProcessing/Trajects/Traject_time_#/
					CreateDirectoryA(path,NULL);

					sprintf(tmpPath,"%sTraject_length_",path);				//PostProcessing/Trajects/Traject_time_#/Traject_length_
					compose_name_plus_nr_str(tmpPath,"_ID_",nPoints,path);	//PostProcessing/Trajects/Traject_time_#/Traject_length_#_ID_
					compose_name_plus_nr_str(path,".vtk",q,val4);			//PostProcessing/Trajects/Traject_time_#/Traject_length_#_ID_#.vtk

				}

				fp111= fopen(val4,"w");
				if(!fp111)printf("Error creating file: %s\n",val4);
				//write data to vtk file in simple legacy format
				fprintf(fp111, "# vtk DataFile Version 2.1 \n");
				fprintf(fp111, "Trajects from timestep %d  \n",seq_first+i);
				fprintf(fp111, "ASCII \n \n");
				fprintf(fp111, "DATASET POLYDATA POINTS %d float\n", nPoints);
				fprintf(fp111, "%5.3f %5.3f %5.3f \n", 
						pointDataArray[i][l].x1 + x_offset, 
						pointDataArray[i][l].y1 + y_offset, 
						pointDataArray[i][l].z1 + z_offset);
				//Create a Matrix to calculate Velocities
				if (nOfP == 0) nOfP=1;
				//reset veloPos1 with zeros
				for(j=0;j<(seq_last - seq_first);j++)
				{
					veloPos1[j].n = -1;
					veloPos1[j].p = -2;
					veloPos1[j].type = 0;
					veloPos1[j].Um = 0;
					veloPos1[j].Ux = 0;
					veloPos1[j].Uy = 0;
					veloPos1[j].Uz = 0;
					veloPos1[j].x1 = 0;
					veloPos1[j].y1 = 0;
					veloPos1[j].z1 = 0;			
				}
				veloPos1[0].x1=pointDataArray[i][l].x1;
				veloPos1[0].y1=pointDataArray[i][l].y1;
				veloPos1[0].z1=pointDataArray[i][l].z1;
				lineOfP=pointDataArray[i][l].n;
				for(timeID=i+1;timeID<=(seq_last-seq_first) && lineOfP >= 0 && lineOfP <= Ms;timeID++)
				{
					time = timeID;
					//append point data of one single traject
					fprintf(fp111, "%5.3f %5.3f %5.3f \n", 
							pointDataArray[time][lineOfP].x1 + x_offset, 
							pointDataArray[time][lineOfP].y1 + y_offset, 
							pointDataArray[time][lineOfP].z1 + z_offset);
					//add trajectory position-data to array
					veloPos1[tmpVeloPos].x1=pointDataArray[time][lineOfP].x1;
					veloPos1[tmpVeloPos].y1=pointDataArray[time][lineOfP].y1;
					veloPos1[tmpVeloPos].z1=pointDataArray[time][lineOfP].z1;
					//handle last Point of traject
					if(pointDataArray[time][lineOfP].n <0 && pointDataArray[time][lineOfP].p>=0)
					{ }
					//update line of next particle-position
					lineOfP=pointDataArray[time][lineOfP].n;
					tmpVeloPos = tmpVeloPos++;
				}
				tmpVeloPos	=	1;//reset counter

				//add lines
				fprintf(fp111, "LINES %d %d\n", nPoints-1, 3*(nPoints-1));
				for(nOfP=1; nOfP < nPoints; nOfP++)
				{
					fprintf(fp111,"2 %d %d \n", nOfP-1, nOfP);
				}

				//add particle/path data
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

				// add mean Ux
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

				// add mean Uy
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

				// add mean Uz
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

				// add mean Uxyz
				fprintf(fp111, "SCALARS Particle_Path double 1\n");
				fprintf(fp111, "LOOKUP_TABLE default\n");

				for(nOfP=1; nOfP <= nPoints; nOfP++)
				{
					tmpValue=((double)(nOfP))/((double)(nPoints));
					fprintf(fp111,"%2.5f \n", tmpValue);
				}

				//..add additional scalar values here..
			
				fclose(fp111);
			}
			
			progressBar(18,60,(int)((i*100)/(seq_last-seq_first)));			
		}
		free (veloPos1);
	 	printf("Trajectories exported sorted by ID to vtk-Files.\n");
		clearLine(18,70);
	}
//**********************************************************************//
//**********// Export stacks for analysis of patterns //**********//	
	if(flag4export[3]==1)
	{
		printf("Writing stacks for pattern analysis.\n");
		CreateDirectoryA("PostProcessing/Stacks/",NULL);	
		nStack = 0;
		if (stacksize > seq_last) stacksize=seq_last; //limit stacksize to values that make sense
		for(i=0; i < (seq_last - seq_first); i+=stacksize)
		{
			nOfP = 0;
			if(i+stacksize+seq_first > seq_last) stacksize = seq_last - seq_first - i; //adapt size of last stack 
			compose_name_plus_nr_str("PostProcessing/Stacks/Stack_",".vtk",nStack+1,val3);
			fp99 = fopen(val3,"w");
			if(!fp99)printf("Could not create file %s.\n",val3);
			//Write header of data file
			fprintf(fp99, "# vtk DataFile Version 2.1 \n");
			fprintf(fp99, "Trajects from stack %d with color option \n",nStack+1);
			fprintf(fp99, "ASCII \n \n");

			for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
			{
				//j -> line1 // j+1 -> line2
				// calc/sum up number of lines for stack:
				for(k=0;k<pointDataArray[j][0].type;k++)
				{ 	
					m = pointDataArray[j][k].n;
					if (m >= 0 && m <= pointDataArray[j+1][0].type){nOfP = nOfP++;}//count connected particles
				}
			}
			fprintf(fp99, "DATASET POLYDATA POINTS %d float \n", 2*nOfP );

			for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
			{
				//j -> line1 // j+1 -> line2
				// write points of lines for stack:
				for(k=0;k<pointDataArray[j][0].type;k++)
				{ 	
					m = pointDataArray[j][k].n;
					if (m >= 0 && m <= pointDataArray[j+1][0].type)
					{
						if( pointDataArray[j][k].x1 < -100000 || pointDataArray[j][k].y1 < -100000 || pointDataArray[j][k].z1 < -100000||
						pointDataArray[j+1][m].x1 < -100000 || pointDataArray[j+1][m].y1 < -100000 || pointDataArray[j+1][m].z1 < -100000)
						{
							pointDataArray[j][k].x1=0.0;	pointDataArray[j][k].y1=0.0;	pointDataArray[j][k].z1=0.0;
							pointDataArray[j+1][m].x1=0.0;	pointDataArray[j+1][m].y1=0.0;	pointDataArray[j+1][m].z1=0.0;
							fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
								pointDataArray[j][k].x1 + x_offset, 
								pointDataArray[j][k].y1 + y_offset, 
								pointDataArray[j][k].z1 + z_offset);
							fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
								pointDataArray[j+1][m].x1 + x_offset, 
								pointDataArray[j+1][m].y1 + y_offset, 
								pointDataArray[j+1][m].z1 + z_offset); 
						}
						else if(pointDataArray[j][k].x1 > 100000 || pointDataArray[j][k].y1 > 100000 || pointDataArray[j][k].z1 > 100000|| 
							pointDataArray[j+1][m].x1 > 100000 || pointDataArray[j+1][m].y1 > 100000 || pointDataArray[j+1][m].z1 > 100000)
						{
							pointDataArray[j][k].x1=0.0;	pointDataArray[j][k].y1=0.0;	pointDataArray[j][k].z1=0.0;
							pointDataArray[j+1][m].x1=0.0;	pointDataArray[j+1][m].y1=0.0;	pointDataArray[j+1][m].z1=0.0;
							fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
								pointDataArray[j][k].x1 + x_offset, 
								pointDataArray[j][k].y1 + y_offset, 
								pointDataArray[j][k].z1 + z_offset);
							fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
								pointDataArray[j+1][m].x1 + x_offset, 
								pointDataArray[j+1][m].y1 + y_offset, 
								pointDataArray[j+1][m].z1 + z_offset); 
						}
						else
						{
							 fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
								 pointDataArray[j][k].x1 + x_offset, 
								 pointDataArray[j][k].y1 + y_offset, 
								 pointDataArray[j][k].z1 + z_offset);
							 fprintf(fp99, "%5.3f %5.3f %5.3f \n", 
								 pointDataArray[j+1][m].x1 + x_offset, 
								 pointDataArray[j+1][m].y1 + y_offset, 
								 pointDataArray[j+1][m].z1 + z_offset); 
						}
					}
				}
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
				//j -> line1 // j+1 -> line2
				// calc/sum up number of lines for stack:
				for(k=0;k<pointDataArray[j][0].type;k++)
				{ 	
					m = pointDataArray[j][k].n;
					if (m >= 0 && m <= pointDataArray[j+1][0].type)
					{
						dx=pointDataArray[j][k].x1-pointDataArray[j+1][m].x1;
						dy=pointDataArray[j][k].y1-pointDataArray[j+1][m].y1;
						dz=pointDataArray[j][k].z1-pointDataArray[j+1][m].z1;
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
			}

			//append mean x-veloData to "lines"
			fprintf(fp99, "SCALARS Umean_x double 1 \n");
			fprintf(fp99, "LOOKUP_TABLE default\n");
			for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
			{
				//j -> line1 // j+1 -> line2
				// calc/sum up number of lines for stack:
				for(k=0;k<pointDataArray[j][0].type;k++)
				{ 	
					m = pointDataArray[j][k].n;
					if (m >= 0 && m <= pointDataArray[j+1][0].type)
					{
						dx=pointDataArray[j][k].x1-pointDataArray[j+1][m].x1;
						Ux=(-dx/interframeTime);
						fprintf(fp99,"%2.5f \n", Ux);
						fprintf(fp99,"%2.5f \n", Ux);
					}
				}
			}

			//append mean y-veloData to "lines"
			fprintf(fp99, "SCALARS Umean_y double 1 \n");
			fprintf(fp99, "LOOKUP_TABLE default\n");
			for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
			{
				//j -> line1 // j+1 -> line2
				// calc/sum up number of lines for stack:
				for(k=0;k<pointDataArray[j][0].type;k++)
				{ 	
					m = pointDataArray[j][k].n;
					if (m >= 0 && m <= pointDataArray[j+1][0].type)
					{
						dx=pointDataArray[j][k].y1-pointDataArray[j+1][m].y1;
						Uy=(-dy/interframeTime);
						fprintf(fp99,"%2.5f \n", Uy);
						fprintf(fp99,"%2.5f \n", Uy);
					}
				}
			}

			//append mean z-veloData to "lines"
			fprintf(fp99, "SCALARS Umean_z double 1 \n");
			fprintf(fp99, "LOOKUP_TABLE default\n");
			for(j=i; j<i+stacksize-2; j++)//inner loop to run through stack data only
			{
				//j -> line1 // j+1 -> line2
				// calc/sum up number of lines for stack:
				for(k=0;k<pointDataArray[j][0].type;k++)
				{ 	
					m = pointDataArray[j][k].n;
					if (m >= 0 && m <= pointDataArray[j+1][0].type)
					{
						dz=pointDataArray[j][k].z1-pointDataArray[j+1][m].z1;
						Uz=(-dz/interframeTime);
						fprintf(fp99,"%2.5f \n", Uz);
						fprintf(fp99,"%2.5f \n", Uz);
					}
				}
			}
			fclose(fp99);
			nStack++;//counter for stacks
			
			progressBar(18,60,(int)((i*100)/(seq_last-seq_first)));			
		}
		printf("Finished writing stacks for pattern analysis.\n");
		clearLine(18,70);
	}
//**********************************************************************//
//**********// do some  scientific stuff //**********//	
	

	//free allocated memory
	if (pointDataArray)
	{
		for(i=0;i < (seq_last-seq_first+1); i++)
		{
			free(pointDataArray[i]);
		}	
		free (pointDataArray);
	}

//**********************************************************************//
   /**Added export to paraview for correspondending points**/
	if(flag4export[0]==1)//option "a"
	{
		for (i=seq_first; i<seq_last;i++)
		{
			compose_name_plus_nr (resultFolder, "rt_is.",i+1, val);
      
			strcpy(tmp_path, base_path);
			strcpy(tmp_file_name, val);
			strcat(tmp_path, tmp_file_name);
			fp1 = fopen (tmp_path, "r");   
			if (! fp1)  showErrorMessageOnFileOpening(tmp_path);
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

			compose_name_plus_nr_str("PostProcessing/prd_3dPoints_",".vtk",i,val2);

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
	return 1;
}