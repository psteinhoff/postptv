#include <ptv_is_2_unstructured_vtk.h>


/* tools / makros */
int round (double x)
/* round input double x to integer */
{
 if (x>=0)
   return((int)(x+.5));
 else
   return((int)(x-.5));
}

FILE *fopen_r (char *filename)
{
	//char	filename[256];
/*	tries to open a file;
	gives a message, if it cannot open it
	and waits until it has been created 	 */
  FILE	*fpr;
  int  	count;

  fpr = fopen (filename, "r");
  if ( ! fpr)
    {
      printf ("could not open %s, please create this file\n", filename);

      /* wait until file can be opened */
      while ( ! fpr)	fpr = fopen (filename, "r");

      /* wait until file really created */
      for (count=0; count<100000; count++);
    }

  return (fpr);
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

/*************************************************************************/
void compose_name_plus_nr (char *basename, char *str, int nr, char *filename)
{
  char	nr_ch[256];

  if		(nr < 10)		sprintf (nr_ch, "%1d", nr);
  else if	(nr < 100)		sprintf (nr_ch, "%2d", nr);
  else if   (nr < 1000)		sprintf (nr_ch, "%3d", nr);
  else if   (nr < 10000)	sprintf (nr_ch, "%4d", nr);
  else						sprintf (nr_ch, "%5d", nr);

  sprintf (filename, "%s%s%s", basename, str, nr_ch);
}

void compose_name_plus_nr_str (char *basename, char *str, int nr, char *filename)
{
  char	nr_ch[256];

  if		(nr < 10)		sprintf (nr_ch, "%1d", nr);
  else if	(nr < 100)		sprintf (nr_ch, "%2d", nr);
  else if   (nr < 1000)		sprintf (nr_ch, "%3d", nr);
  else if   (nr < 10000)	sprintf (nr_ch, "%4d", nr);
  else						sprintf (nr_ch, "%5d", nr);

  sprintf (filename, "%s%s%s", basename, nr_ch, str);
}

void showErrorMessageOnFileOpening (char *filename)
{
	printf("Can't open ascii file %s!\n",filename);
	getch();
}

int countConnectedParticles (int size1, vector *line1, int size2, vector *line2)
{
	int j,nOfP,m;
	nOfP=0;
	for(j=0;j<size1;j++) 
	{ 	
		m = line1[j].n;
		if (m >= 0 && m <= size2 && line2[m].p>0)
		{
			//count connected particles
			nOfP = nOfP++;
		}
	 }
	return nOfP;
}

int countValidatedParticleConnections(vector *line, int upperLimit)
{
	int i;
	int j,nOfP,m;
	nOfP=0;

	for(j=0;j<line[0].type;j++) 
	{ 	
		m = line[j].n;
		if (m >= 0 && m <= upperLimit)
		{
			//count connected particles
			nOfP = nOfP++;
		}
	 }
	return nOfP;
}

/***********************************************************************/
void printSimpleText(int x, int y, char const *text)
{
    COORD	target = {x, y};
    DWORD	written;
 
    WriteConsoleOutputCharacterA(GetStdHandle(STD_OUTPUT_HANDLE), text,
                                            strlen(text),
                                            target, &written);
}

void printSimpleLine(int endOfline,int line,int textOffset, char const *newText)
{
	int	i;
	
	for(i = 0; i < endOfline; ++i)
		{     
			printSimpleText(i,line, " ");
		}
	printSimpleText(textOffset,line, newText);

}

void printHeader()
{	
	int lineEnd,offset;

	lineEnd = 50;
	offset = 0;

	//system("cls");
	printSimpleLine(lineEnd,1,offset, "====================================================");
	printSimpleLine(lineEnd,2,offset, "\\_\\_\\_\\_  \\_  \\_\\_\\_\\_\\_  \\_   Rheinisch-");
	printSimpleLine(lineEnd,3,offset, "\\_   \\_\\_ \\__ \\_  \\_  \\_  \\_   Westfaelische");
	printSimpleLine(lineEnd,4,offset, "\\_\\_\\_  \\_\\_\\_\\_  \\_  \\_\\_\\_   Technische");
	printSimpleLine(lineEnd,5,offset, "\\_  \\_   \\__ \\__  \\_  \\_  \\_   Hochschule");
	printSimpleLine(lineEnd,6,offset, "\\_   \\_   \\_  \\_  \\_  \\_  \\_   Aachen");
	printSimpleLine(lineEnd,7,offset, "====================================================");
	/*printf("\\_\\_\\_\\_  \\_  \\_\\_\\_\\_\\_  \\_   Rheinisch-	\n"); 
	printf("\\_   \\_\\_ \\__ \\_  \\_  \\_  \\_   Westfaelische\n"); 
	printf("\\_\\_\\_  \\_\\_\\_\\_  \\_  \\_\\_\\_   Technische	\n"); 
	printf("\\_  \\_   \\__ \\__  \\_  \\_  \\_   Hochschule	\n"); 
	printf("\\_   \\_   \\_  \\_  \\_  \\_  \\_   Aachen		\n"); 
	printf("============================================\n"); */
	//#7
	printSimpleLine(lineEnd,8,offset, "\\_\\_\\_ \\_\\_  \\_\\_\\_            Energy efficient");
	printSimpleLine(lineEnd,9,offset, "\\_     \\_  \\ \\_                Buildings and indoor");
	printSimpleLine(lineEnd,10,offset, "\\_\\_\\_ \\_\\_\\ \\_                Climate");
	printSimpleLine(lineEnd,11,offset, "\\_     \\_  \\ \\_");
	printSimpleLine(lineEnd,12,offset, "\\_\\_\\_ \\_\\_\\ \\_\\_\\_");
	printSimpleLine(lineEnd,13,offset, "====================================================");
	/*printf("\\_\\_\\_ \\_\\_  \\_\\_\\_	Energy efficient\n");
	printf("\\_     \\_  \\ \\_		Buildings and indoor	\n");
	printf("\\_\\_\\_ \\_\\_\\ \\_		Climate	\n");
	printf("\\_     \\_  \\ \\_			\n");
	printf("\\_\\_\\_ \\_\\_\\ \\_\\_\\_		\n");
	printf("============================================\n");*/
	//#13
	printSimpleLine(lineEnd,14,offset, "\\_\\_\\_\\_\\_\\_\\_\\_      \\_       Particle");
	printSimpleLine(lineEnd,15,offset, "\\_   \\_  \\_    \\_    \\_        Tracking");
	printSimpleLine(lineEnd,16,offset, "\\_\\_\\_   \\_     \\_  \\_         Velocimetry");
	printSimpleLine(lineEnd,17,offset, "\\_       \\_      \\_\\_");
	printSimpleLine(lineEnd,18,offset, "\\_       \\_       \\_");
	printSimpleLine(lineEnd,19,offset, "====================================================");
	/*printf("\\_\\_\\_\\_\\_\\_\\_\\_      \\_	Particle\n"); 
	printf("\\_   \\_  \\_    \\_    \\_		Tracking\n"); 
	printf("\\_\\_\\_   \\_     \\_  \\_		Velocimetry\n"); 
	printf("\\_       \\_      \\_\\_  	\n"); 
	printf("\\_       \\_       \\_  \n"); 
	printf("============================================\n"); */
	//#19

}

void progressBar(int line, int width, int progressInPercent)
{
	char	textToPrint[128];
	char	tmp[128];
	int		i, progress;

	if (width < 10) width = 10;
	if (width > 60) width = 60;
	sprintf(textToPrint,"%3d %%",progressInPercent); 
	
	//strcpy(tmp,""); 
	memset(tmp,0,sizeof(tmp));
	progress = (int)(progressInPercent*width/100);
	tmp[0]='[';
	for(i=1;i<60;i++)
		{
			if (i< width)tmp[i]='.';
			else tmp[i]=' ';
		}
	for(i=1;i<progress+1;i++)
		{
			tmp[i]='='; 
		}
	tmp[width]=']';
	strcat(textToPrint,tmp);
	printSimpleText(1,line,textToPrint);
}

void clearLine(int line, int width)
{
	printSimpleLine(width,line,0,"");
}